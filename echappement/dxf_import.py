"""Lecture du tracé d'échappement à ancre suisse depuis un DXF.

Le DXF (export de CAO, un seul calque, sans blocs) ne dit pas quelle pièce est
quoi : les pièces sont reconnues par leur GÉOMÉTRIE, en regroupant les
éléments qui se touchent bout à bout.

  - la roue d'échappement : le plus gros contour ;
  - l'ancre : le contour qui part de son axe (bras, palettes, baguette),
    complété par les traits isolés de la fourchette ;
  - l'ellipse du balancier : le contour fermé fait d'arcs de grand rayon ;
  - cinq cercles : les axes de la roue, de l'ancre et du balancier, et les
    deux butées (même rayon).

Le sens de rotation de la roue et la palette d'entrée se déduisent de la forme
des palettes : sur chacune, la dent atteint d'abord l'arête de repos, puis
soulève la palette en glissant sur le plan d'impulsion vers l'arête de sortie.
L'arête de repos est donc celle des deux extrémités du plan d'impulsion qui est
la plus éloignée du centre de la roue. La palette d'entrée est celle que la
dent rencontre en premier.

Le tracé est ramené à la VUE STANDARD : balancier en haut, roue en bas,
palette d'entrée à gauche. Si le DXF n'y est pas déjà, il est retourné haut-bas
et/ou gauche-droite ; les deux ensemble font une rotation de 180°. La
transformation appliquée est rendue dans Trace.transformation.

Les éléments dessinés deux fois sont écartés à la lecture.
"""

import collections
from dataclasses import dataclass

import ezdxf
import numpy as np

from .geometrie import Corps, Element, echantillonner_arc

#: deux extrémités plus proches que cela sont jointes (mm)
TOL_JONCTION = 1e-5


@dataclass(eq=False)
class Trace:
    """Les pièces de l'échappement, dans la pose dessinée, en vue standard."""
    roue: Corps
    ancre: Corps
    balancier: Corps
    butees: list
    nombre_dents: int
    rayon_tete: float
    sens_roue: int               # sens de rotation de la roue (+1 trigo, -1 horaire)
    transformation: str = "aucune"

    @property
    def pas(self):
        return 2 * np.pi / self.nombre_dents


def _lire_entites(chemin):
    """Éléments du DXF, tels que dessinés, doublons écartés : (segments et
    arcs sous forme brute, cercles)."""
    vus, bruts, cercles = set(), [], []
    for e in ezdxf.readfile(chemin).modelspace():
        t = e.dxftype()
        if t == "LINE":
            a, b = tuple(np.array(e.dxf.start)[:2]), tuple(np.array(e.dxf.end)[:2])
            cle = ("L",) + tuple(np.round(sorted([a, b]), 6).ravel())
            brut = ("L", np.array([a, b]))
        elif t == "ARC":
            c = np.array(e.dxf.center)[:2]
            cle = ("A",) + tuple(np.round([*c, e.dxf.radius, e.dxf.start_angle,
                                           e.dxf.end_angle], 6))
            brut = ("A", c, e.dxf.radius, np.radians(e.dxf.start_angle),
                    np.radians(e.dxf.end_angle))
        elif t == "CIRCLE":
            cercles.append((np.array(e.dxf.center)[:2], float(e.dxf.radius)))
            continue
        else:
            continue
        if cle not in vus:
            vus.add(cle)
            bruts.append(brut)
    return bruts, cercles


def _elements(bruts, sx, sy):
    """Éléments échantillonnés, symétrie (sx, sy) appliquée (±1 par axe)."""
    s = np.array([sx, sy], dtype=float)
    out = []
    for b in bruts:
        if b[0] == "L":
            out.append(Element(points=b[1] * s))
        else:
            _, c, r, a0, a1 = b
            out.append(Element(points=echantillonner_arc(c, r, a0, a1) * s,
                               centre=c * s, rayon=r))
    return out


def _groupes(elements):
    """Regroupe les éléments qui se touchent bout à bout."""
    n = len(elements)
    parent = list(range(n))

    def racine(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    bouts = np.array([p for e in elements for p in (e.points[0], e.points[-1])])
    proprio = np.repeat(np.arange(n), 2)
    ordre = np.lexsort((bouts[:, 1], bouts[:, 0]))
    for k in range(len(ordre)):
        for m in range(k + 1, len(ordre)):
            a, b = bouts[ordre[k]], bouts[ordre[m]]
            if b[0] - a[0] > TOL_JONCTION:
                break
            if abs(b[1] - a[1]) <= TOL_JONCTION:
                ra, rb = racine(proprio[ordre[k]]), racine(proprio[ordre[m]])
                if ra != rb:
                    parent[ra] = rb
    g = collections.defaultdict(list)
    for i in range(n):
        g[racine(i)].append(elements[i])
    return sorted(g.values(), key=len, reverse=True)


def _cross(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def _etiqueter_palettes(ancre, centre_roue, axe_balancier):
    """Étiquette les faces des deux palettes, côté par côté (« gauche »,
    « droite »), et rend leurs arêtes {côté: (arête de repos, arête de sortie)}.

    Le plan d'impulsion est la face la plus proche du centre de la roue ; son
    extrémité la plus éloignée de ce centre est l'arête de repos ; la face de
    repos est l'autre face qui part de cette arête. Les bras, qui partent de
    l'axe de l'ancre, sont exclus."""
    r = lambda p: float(np.hypot(*(np.asarray(p) - centre_roue)))
    aretes = {}
    for cote, signe in (("gauche", -1), ("droite", +1)):
        lignes = [e for e in ancre.elements if not e.est_arc and not e.etiquette
                  and np.sign(e.points.mean(0)[0] - ancre.axe[0]) == signe
                  and min(np.hypot(*(e.points - ancre.axe).T)) > 1e-4
                  and e.longueur > 0.15
                  # côté roue, pas côté fourchette
                  and np.hypot(*(e.points.mean(0) - centre_roue)) <
                  np.hypot(*(e.points.mean(0) - axe_balancier))]
        impulsion = min(lignes, key=lambda e: np.mean([r(p) for p in e.points]))
        a, b = impulsion.points[0], impulsion.points[-1]
        coin_repos, coin_sortie = (a, b) if r(a) > r(b) else (b, a)
        repos = min((e for e in lignes if e is not impulsion),
                    key=lambda e: min(np.hypot(*(e.points - coin_repos).T)))
        impulsion.etiquette = f"plan d'impulsion {cote}"
        repos.etiquette = f"face de repos {cote}"
        for e in ancre.elements:
            if e.est_arc and not e.etiquette:
                if np.hypot(*(e.centre - coin_repos)) < 0.02:
                    e.etiquette = f"arête de repos {cote}"
                elif np.hypot(*(e.centre - coin_sortie)) < 0.02:
                    e.etiquette = f"arête de sortie {cote}"
        for e in lignes:
            if not e.etiquette:
                e.etiquette = f"dos de palette {cote}"
        aretes[cote] = (coin_repos, coin_sortie)
    return aretes


def _etiqueter_fourchette(traits, axe_ancre):
    """Flancs d'entrée (les deux traits parallèles à la baguette les plus
    proches de son axe), cornes (les deux suivants) et fond de fourchette."""
    baguette = next(e for e in traits if min(np.hypot(*(e.points - axe_ancre).T)) < 1e-4)
    baguette.etiquette = "baguette"
    u = baguette.points[-1] - baguette.points[0]
    u /= np.hypot(*u)
    if np.hypot(*(baguette.points[0] - axe_ancre)) > 1e-4:
        u = -u
    paralleles = [e for e in traits if e is not baguette
                  and abs(_cross(u, (e.points[-1] - e.points[0]) / e.longueur)) < 0.05]
    for e in traits:
        if e is not baguette and e not in paralleles:
            e.etiquette = "fond de fourchette"
    cote = lambda e: _cross(u, e.points.mean(0) - axe_ancre)
    paralleles.sort(key=lambda e: abs(cote(e)))
    for e in paralleles[:2]:
        e.etiquette = "flanc d'entrée " + ("droit" if cote(e) < 0 else "gauche")
    for e in paralleles[2:]:
        e.etiquette = "corne " + ("droite" if cote(e) < 0 else "gauche")


def _construire(bruts, cercles, sx, sy):
    """Reconnaît les pièces dans le dessin pris avec la symétrie (sx, sy).
    Rend (Trace, palette d'entrée « gauche »/« droite »)."""
    s = np.array([sx, sy], dtype=float)
    cercles = [(c * s, r) for c, r in cercles]
    groupes = _groupes(_elements(bruts, sx, sy))

    roue_el = groupes[0]
    pts = np.concatenate([e.points for e in roue_el])
    centre_roue = min((c for c, _ in cercles), key=lambda c: np.hypot(*(c - pts.mean(0))))

    rayons = collections.Counter(round(r, 6) for _, r in cercles)
    r_butee = next(r for r, n in rayons.items() if n == 2)
    butees_c = [(c, r) for c, r in cercles if round(r, 6) == r_butee]
    axes = [c for c, r in cercles if round(r, 6) != r_butee
            and np.hypot(*(c - centre_roue)) > 1e-6]

    def part_de(g, c):
        return any(min(np.hypot(*(e.points[[0, -1]] - c).T)) < 1e-4 for e in g)
    ancre_el = next(g for g in groupes[1:] if len(g) > 4)
    axe_ancre = next(c for c in axes if part_de(ancre_el, c))
    axe_balancier = next(c for c in axes if c is not axe_ancre)

    # ellipse : contour fait d'arcs de grand rayon ; le rayon du balancier
    # (trait de l'axe à l'ellipse) n'est pas une surface de contact
    ellipse_el = next(g for g in groupes if any(e.est_arc and e.rayon > 0.1 for e in g))
    for e in ellipse_el:
        e.etiquette = "ellipse"
    restants = [g for g in groupes if g is not roue_el and g is not ancre_el
                and g is not ellipse_el]
    traits = [e for g in restants for e in g if not part_de([e], axe_balancier)]
    baguette = [e for e in ancre_el if not e.est_arc
                and min(np.hypot(*(e.points - axe_ancre).T)) < 1e-4
                and np.hypot(*(e.points - axe_ancre).T).max() > 0.5
                and np.hypot(*(e.points - axe_balancier).T).min()
                < np.hypot(*(axe_ancre - axe_balancier))]
    _etiqueter_fourchette(traits + baguette[:1], axe_ancre)

    for e in roue_el:
        r = np.hypot(*(e.points - centre_roue).T)
        if r.min() > 2.3 and not e.est_arc and 0.1 < e.longueur < 0.14:
            e.etiquette = "plan de la dent"
        elif r.min() > 2.3:
            e.etiquette = "tête de dent"
        else:
            e.etiquette = "dent"
    roue = Corps("roue", centre_roue, roue_el)
    ancre = Corps("ancre", axe_ancre, ancre_el + traits)
    aretes = _etiqueter_palettes(ancre, centre_roue, axe_balancier)
    for e in ancre.elements:
        if not e.etiquette:
            e.etiquette = "bras" if min(np.hypot(*(e.points - axe_ancre).T)) < 1e-4 else "ancre"
    balancier = Corps("balancier", axe_balancier, ellipse_el)
    butees = []
    for c, r in sorted(butees_c, key=lambda cr: cr[0][0]):
        nom = "butée " + ("gauche" if c[0] < axe_ancre[0] else "droite")
        el = Element(points=echantillonner_arc(c, r, 0.0, 2 * np.pi), etiquette=nom,
                     centre=c, rayon=r)
        butees.append(Corps(nom, c, [el], fixe=True))

    # Sens de la roue : de l'arête de repos vers l'arête de sortie.
    # Palette d'entrée : la dent la rencontre avant l'autre.
    repos_g, sortie_g = aretes["gauche"]
    sens = int(np.sign(_cross(repos_g - centre_roue, sortie_g - centre_roue)))
    ang = {c: np.arctan2(*(aretes[c][0] - centre_roue)[::-1]) for c in aretes}
    ecart = (ang["droite"] - ang["gauche"] + np.pi) % (2 * np.pi) - np.pi
    entree = "gauche" if np.sign(ecart) == sens else "droite"
    sortie = "droite" if entree == "gauche" else "gauche"
    prefixes = ("plan d'impulsion", "face de repos", "arête de", "dos de palette")
    for e in ancre.elements:
        for cote, nom in ((entree, "entrée"), (sortie, "sortie")):
            if e.etiquette.startswith(prefixes) and e.etiquette.endswith(" " + cote):
                e.etiquette = e.etiquette[: -len(cote)] + nom

    r_pts = np.hypot(*(pts - centre_roue).T)
    rayon_tete = float(r_pts.max())
    pointes = pts[r_pts > rayon_tete - 0.01]
    a = np.sort(np.arctan2(pointes[:, 1] - centre_roue[1], pointes[:, 0] - centre_roue[0]))
    nombre = int(np.sum(np.diff(np.concatenate([a, [a[0] + 2 * np.pi]])) > np.radians(3)))
    return Trace(roue, ancre, balancier, butees, nombre, rayon_tete, sens), entree


def lire_trace(chemin):
    """Lit le DXF et rend le Trace en vue standard (balancier en haut, roue en
    bas, palette d'entrée à gauche), en retournant le dessin si besoin."""
    bruts, cercles = _lire_entites(chemin)
    sx, sy = 1, 1
    trace, entree = _construire(bruts, cercles, sx, sy)
    if trace.roue.axe[1] > trace.ancre.axe[1]:          # roue en haut : retourner
        sy = -1
        trace, entree = _construire(bruts, cercles, sx, sy)
    if entree == "droite":
        sx = -1
        trace, entree = _construire(bruts, cercles, sx, sy)
    trace.transformation = {(1, 1): "aucune", (1, -1): "symétrie haut-bas",
                            (-1, 1): "symétrie gauche-droite",
                            (-1, -1): "rotation de 180°"}[(sx, sy)]
    return trace
