"""Variables du tracé : modifications du tracé de référence (le DXF).

Plutôt que de redessiner l'échappement de toutes pièces, chaque variable
déplace ou réoriente une partie du tracé lu dans le DXF ; à valeurs nulles,
on retrouve exactement le dessin. C'est ce que l'optimisation fait varier.

Variables globales :
  - entraxe roue-ancre, entraxe ancre-balancier (mm, valeurs absolues) ;
  - rayon de l'ellipse (mm, distance de son centre à l'axe du balancier) ;
  - largeur de l'entrée de fourchette (mm) : flancs et cornes s'écartent
    symétriquement ;
  - course de l'ancre (°) : les butées tournent autour de l'axe de l'ancre.

Variables de chaque palette (entrée, sortie), reconstruite avec ses congés :
  - position (°) : rotation de la palette autour de l'axe de l'ancre ; elle
    répartit repos et chute ;
  - inclinaison du plan d'impulsion (°) : rotation autour de l'arête de
    repos ; elle règle la levée ;
  - longueur du plan d'impulsion (µm) : l'arête de sortie avance ou recule
    le long du plan ;
  - tirage (°) : rotation de la face de repos autour de l'arête de repos.
"""

import copy
from dataclasses import asdict, dataclass, fields

import numpy as np

from .dxf_import import Trace
from .geometrie import Corps, Element, echantillonner_arc, rotation

DEG = np.pi / 180.0


@dataclass
class VariablesTrace:
    entraxe_roue_ancre: float = np.nan        # mm ; nan : celui du dessin
    entraxe_ancre_balancier: float = np.nan   # mm
    rayon_ellipse: float = np.nan             # mm
    largeur_fourchette: float = np.nan        # mm
    course_ancre: float = np.nan              # °
    position_entree: float = 0.0              # °
    impulsion_entree: float = 0.0             # °
    longueur_entree: float = 0.0              # µm
    tirage_entree: float = 0.0                # °
    position_sortie: float = 0.0
    impulsion_sortie: float = 0.0
    longueur_sortie: float = 0.0
    tirage_sortie: float = 0.0

    def en_dict(self):
        return asdict(self)


#: libellés et unités, pour l'interface et l'optimisation
LIBELLES = {
    "entraxe_roue_ancre": ("Entraxe roue-ancre", "mm"),
    "entraxe_ancre_balancier": ("Entraxe ancre-balancier", "mm"),
    "rayon_ellipse": ("Rayon de l'ellipse", "mm"),
    "largeur_fourchette": ("Largeur d'entrée de fourchette", "mm"),
    "course_ancre": ("Course de l'ancre", "°"),
    "position_entree": ("Position palette d'entrée", "°"),
    "impulsion_entree": ("Inclinaison plan d'impulsion d'entrée", "°"),
    "longueur_entree": ("Longueur plan d'impulsion d'entrée", "µm"),
    "tirage_entree": ("Tirage d'entrée", "°"),
    "position_sortie": ("Position palette de sortie", "°"),
    "impulsion_sortie": ("Inclinaison plan d'impulsion de sortie", "°"),
    "longueur_sortie": ("Longueur plan d'impulsion de sortie", "µm"),
    "tirage_sortie": ("Tirage de sortie", "°"),
}


def _intersection(p, u, q, v):
    """Intersection des droites p + s.u et q + t.v."""
    M = np.array([[u[0], -v[0]], [u[1], -v[1]]])
    s, _ = np.linalg.solve(M, q - p)
    return p + s * u


def _unit(v):
    return v / np.hypot(*v)


def _conge(C, a, b, r):
    """Congé de rayon r au coin C entre les demi-droites de directions
    unitaires a et b (partant de C). Rend (T_a, T_b, points de l'arc de T_a à
    T_b, centre)."""
    cos_t = np.clip(a @ b, -1.0, 1.0)
    theta = np.arccos(cos_t)
    d = r / np.tan(theta / 2)
    Ta, Tb = C + d * a, C + d * b
    bis = _unit(a + b)
    centre = C + (r / np.sin(theta / 2)) * bis
    a0 = np.arctan2(*(Ta - centre)[::-1])
    a1 = np.arctan2(*(Tb - centre)[::-1])
    span = np.mod(a1 - a0, 2 * np.pi)
    if span <= np.pi:
        pts = echantillonner_arc(centre, r, a0, a0 + span)
    else:
        pts = echantillonner_arc(centre, r, a1, a1 + (2 * np.pi - span))[::-1]
    return Ta, Tb, pts, centre


def mesurer(trace):
    """Valeurs de référence des variables globales, lues sur le tracé."""
    t = trace
    c_ell = next(e.centre for e in t.balancier.elements if e.est_arc and e.rayon > 0.1)
    flancs = [e for e in t.ancre.elements if e.etiquette.startswith("flanc d'entrée")]
    bag = t.ancre.element("baguette").points
    u = _unit(bag[-1] - bag[0])
    n = np.array([-u[1], u[0]])
    largeur = abs((flancs[0].points.mean(0) - flancs[1].points.mean(0)) @ n)
    from .tables import course_butees
    a = course_butees(t)
    return {"entraxe_roue_ancre": float(np.hypot(*(t.roue.axe - t.ancre.axe))),
            "entraxe_ancre_balancier": float(np.hypot(*(t.balancier.axe - t.ancre.axe))),
            "rayon_ellipse": float(np.hypot(*(c_ell - t.balancier.axe))),
            "largeur_fourchette": float(largeur),
            "course_ancre": float((max(a) - min(a)) / DEG)}


def _translater(corps, delta, axe_aussi=True):
    c = copy.deepcopy(corps)
    for e in c.elements:
        e.points = e.points + delta
        if e.centre is not None:
            e.centre = e.centre + delta
    if axe_aussi:
        c.axe = c.axe + delta
    c._reconstruire()
    return c


def _palette(ancre, cote, d_pos, d_imp, d_long, d_tir):
    """Reconstruit la palette `cote` de l'ancre avec les modifications."""
    el = {lab: ancre.element(f"{lab} {cote}") for lab in
          ("face de repos", "plan d'impulsion", "dos de palette", "arête de repos", "arête de sortie")}
    R, I, D = el["face de repos"].points, el["plan d'impulsion"].points, el["dos de palette"].points
    uR, uI, uD = _unit(R[-1] - R[0]), _unit(I[-1] - I[0]), _unit(D[-1] - D[0])
    L = _intersection(R[0], uR, I[0], uI)              # arête de repos (vive)
    X = _intersection(I[0], uI, D[0], uD)              # arête de sortie (vive)
    base_R = max((R[0], R[-1]), key=lambda p: np.hypot(*(p - L)))
    base_D = max((D[0], D[-1]), key=lambda p: np.hypot(*(p - X)))
    dir_R, dir_D = _unit(base_R - L), _unit(base_D - X)
    r_rep, r_sor = el["arête de repos"].rayon, el["arête de sortie"].rayon
    # bras de ce côté : droite qui porte les bases
    bras = min((e for e in ancre.elements if e.etiquette == "bras"),
               key=lambda e: np.hypot(*(e.points[-1] - base_R)) + np.hypot(*(e.points[0] - base_R)))
    A0 = ancre.axe
    bout = max((bras.points[0], bras.points[-1]), key=lambda p: np.hypot(*(p - A0)))
    u_bras = _unit(bout - A0)

    rot = lambda v, a: rotation(v, np.zeros(2), a)
    # 1. position : toute la palette tourne autour de l'axe de l'ancre
    L, X = rotation(L, A0, d_pos * DEG), rotation(X, A0, d_pos * DEG)
    dir_R, dir_D = rot(dir_R, d_pos * DEG), rot(dir_D, d_pos * DEG)
    # 2. plan d'impulsion : rotation autour de l'arête de repos, puis longueur
    X = L + rot(X - L, d_imp * DEG)
    X = X + (d_long * 1e-3) * _unit(X - L)
    # 3. tirage : la face de repos tourne autour de l'arête de repos
    dir_R = rot(dir_R, d_tir * DEG)
    # bases sur le bras
    bR = _intersection(L, dir_R, A0, u_bras)
    bD = _intersection(X, dir_D, A0, u_bras)
    # congés
    T1, T2, arc_R, c_R = _conge(L, dir_R, _unit(X - L), r_rep)
    T3, T4, arc_S, c_S = _conge(X, _unit(L - X), dir_D, r_sor)
    nouveaux = [
        Element(np.array([bR, T1]), f"face de repos {cote}"),
        Element(arc_R, f"arête de repos {cote}", centre=c_R, rayon=r_rep),
        Element(np.array([T2, T3]), f"plan d'impulsion {cote}"),
        Element(arc_S, f"arête de sortie {cote}", centre=c_S, rayon=r_sor),
        Element(np.array([T4, bD]), f"dos de palette {cote}"),
    ]
    # le bras va jusqu'à la base la plus éloignée de l'axe
    fin = max((bR, bD), key=lambda p: np.hypot(*(p - A0)))
    nouveau_bras = Element(np.array([A0, fin]), "bras")
    reste = [e for e in ancre.elements
             if not e.etiquette.endswith(" " + cote) and e is not bras]
    return reste + [nouveau_bras] + nouveaux


def construire(ref: Trace, var: VariablesTrace) -> Trace:
    """Nouveau tracé, le tracé de référence modifié par les variables."""
    m = mesurer(ref)
    A0 = ref.ancre.axe

    def cible(nom):
        v = getattr(var, nom)
        return m[nom] if v is None or not np.isfinite(v) else v

    # roue : glisse le long de la ligne des centres
    u_r = _unit(ref.roue.axe - A0)
    roue = _translater(ref.roue, (cible("entraxe_roue_ancre") - m["entraxe_roue_ancre"]) * u_r)
    # balancier : idem ; l'ellipse s'éloigne ou se rapproche de son axe
    u_b = _unit(ref.balancier.axe - A0)
    balancier = _translater(ref.balancier,
                            (cible("entraxe_ancre_balancier") - m["entraxe_ancre_balancier"]) * u_b)
    c_ell = next(e.centre for e in balancier.elements if e.est_arc and e.rayon > 0.1)
    u_e = _unit(c_ell - balancier.axe)
    balancier = _translater(balancier, (cible("rayon_ellipse") - m["rayon_ellipse"]) * u_e,
                            axe_aussi=False)
    # ancre : palettes, puis fourchette
    ancre = copy.deepcopy(ref.ancre)
    for cote in ("entrée", "sortie"):
        ancre.elements = _palette(ancre, cote, getattr(var, f"position_{_s(cote)}"),
                                  getattr(var, f"impulsion_{_s(cote)}"),
                                  getattr(var, f"longueur_{_s(cote)}"),
                                  getattr(var, f"tirage_{_s(cote)}"))
        ancre._reconstruire()
    bag = ancre.element("baguette").points
    u = _unit(bag[-1] - bag[0])
    n = np.array([-u[1], u[0]])
    dl = 0.5 * (cible("largeur_fourchette") - m["largeur_fourchette"])
    for e in ancre.elements:
        if e.etiquette.startswith(("flanc d'entrée", "corne")):
            cote_e = np.sign((e.points.mean(0) - A0) @ n)
            e.points = e.points + cote_e * dl * n
    ancre._reconstruire()
    # butées : tournent autour de l'axe de l'ancre, symétriquement
    d_course = 0.5 * (cible("course_ancre") - m["course_ancre"]) * DEG
    butees = []
    for b in ref.butees:
        c = copy.deepcopy(b)
        # tourner la butée dans le sens qui éloigne la fourchette d'elle
        s = -1.0 if (c.axe - A0)[0] * (bag[-1] - A0)[1] - (c.axe - A0)[1] * (bag[-1] - A0)[0] > 0 else 1.0
        for e in c.elements:
            e.points = rotation(e.points, A0, s * d_course)
            e.centre = rotation(e.centre, A0, s * d_course)
        c.axe = rotation(c.axe, A0, s * d_course)
        c._reconstruire()
        butees.append(c)
    return Trace(roue, ancre, balancier, butees, ref.nombre_dents, ref.rayon_tete,
                 ref.sens_roue, ref.transformation)


def _s(cote):
    return "entree" if cote == "entrée" else "sortie"


def noms_variables():
    return [f.name for f in fields(VariablesTrace)]
