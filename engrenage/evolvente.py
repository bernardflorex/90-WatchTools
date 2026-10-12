"""Génération de profils d'engrenage en développante de cercle.

Port de Fortran/Egc/Evolve.F90 (K. Trümpy). Lit un fichier .Cev et calcule,
pour chaque mobile, un demi-pas de denture (flanc + pieds + têtes) sous forme
d'une suite de NI points, origine (0,0) au centre du mobile.

Usage : python evolvente.py [fichier.json|fichier.Cev] [--ni 150] [--plot] [--png sortie.png] [--dents 3]
(--plot et --png dessinent les dents en prise et demandent matplotlib)
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Roue:
    mo: float  # module
    zz: int  # nombre de dents
    nw: float  # angle de pression (°)
    qe: float  # rapport creux/plein au diamètre primitif
    qg: float  # longueur du raccord de pied (en multiples du module)
    ma: float  # hauteur de départ de la développante sous le primitif (x module)
    me: float  # hauteur de fin de la développante au-dessus du primitif (x module)
    # données de pression de Hertz (facultatives, conservées mais pas encore utilisées)
    e_mod: float = None
    poisson: float = None
    epaisseur: float = None
    narcs: int = 3  # arcs de cercle sur la développante pour l'export .Eno


# Paramètres d'optimisation (onglet Optimisation, voir optim.py et recherche.py) : critères
# "couple" (variation, à minimiser), "jeu" (jeu maximal, à minimiser), "rendement" (perte, à
# minimiser), "hertz" (pression maximale, à minimiser). `libres` : indices des arcs libres par roue
# (None = tous). `deplacement_um` : écart de profil exploré autour du départ ; `jeu_max_admis` : 0
# = jeu du profil de départ à l'entraxe maximal ; `r_tete_min` : rayon minimal de l'arc de tête (mm,
# 0 = pas de limite) ; `dissymetrie` : 1 = réparer les interférences par le demi-profil arrière.
OPTIMISATION_DEFAUT = {"principal": "couple", "secondaire": "hertz", "ponderation": 0.5,
                       "equilibre": 0,
                       "deplacement_um": 10.0, "jeu_min_admis": 0.1, "jeu_max_admis": 0.0,
                       "r_tete_min": 0.0, "population": 24, "evaluations": 12000, "departs": 6,
                       "graine": 1, "dissymetrie": 1, "libres": {"m": None, "e": None}}


@dataclass
class Projet:
    roues: list  # [menante, menée]
    mu: float = 0.2  # coefficient de frottement de la simulation
    mt: float = 1.0  # couple sur la menante (mN.m)
    da: float = 0.0  # variation de l'entraxe (µm)
    da_neg: float = -30.0  # plage de variation d'entraxe étudiée (µm), pour la suite du programme
    da_pos: float = 30.0
    optimisation: dict = field(default_factory=lambda: json.loads(json.dumps(OPTIMISATION_DEFAUT)))


ROLES = ("menante", "menee")  # clés des roues dans le fichier .json
HERTZ = ("e_mod", "poisson", "epaisseur")
FORMAT_JSON = "engcre"


def lire_cev(chemin):
    """Lit un ancien fichier .Cev (import seulement) : ligne 1 en-têtes, lignes 2 et 3 :
    MO ZZ NW QE QG MA ME. Les lignes 5 et 6 (E-module, Poisson, épaisseur) sont facultatives."""
    lignes = Path(chemin).read_text(encoding="latin-1").splitlines()
    if len(lignes) < 3:
        raise ValueError("fichier .Cev incomplet (3 lignes au minimum)")
    roues = []
    for k, ligne in enumerate(lignes[1:3]):
        v = ligne.split()
        roue = Roue(
            mo=float(v[0]),
            zz=int(v[1]),
            nw=float(v[2]),
            qe=float(v[3]),
            qg=float(v[4]),
            ma=float(v[5]),
            me=float(v[6]),
        )
        try:
            h = [float(x) for x in lignes[4 + k].split()[:3]]
            if len(h) == 3:
                roue.e_mod, roue.poisson, roue.epaisseur = h
        except (IndexError, ValueError):
            pass
        roues.append(roue)
    return Projet(roues)


def ecrire_json(chemin, projet):
    """Écrit le projet en JSON structuré : {"format", "version", "roues": {role: {...}}}."""
    roues = {}
    for role, roue in zip(ROLES, projet.roues):
        d = {c: getattr(roue, c) for c in ("mo", "zz", "nw", "qe", "qg", "ma", "me", "narcs")}
        hertz = {c: getattr(roue, c) for c in HERTZ if getattr(roue, c) is not None}
        if hertz:
            d["hertz"] = hertz
        roues[role] = d
    doc = {"format": FORMAT_JSON, "version": 1, "roues": roues, "simulation": {"mu": projet.mu, "mt": projet.mt, "da": projet.da,
                                                 "da_neg": projet.da_neg, "da_pos": projet.da_pos},
           "optimisation": projet.optimisation}
    Path(chemin).write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def lire_json(chemin):
    doc = json.loads(Path(chemin).read_text(encoding="utf-8"))
    if not isinstance(doc, dict) or doc.get("format") != FORMAT_JSON:
        raise ValueError(f"pas un fichier {FORMAT_JSON} (clé \"format\" absente ou différente)")
    roues = []
    for role in ROLES:
        try:
            d = doc["roues"][role]
            roues.append(
                Roue(
                    mo=float(d["mo"]),
                    zz=int(d["zz"]),
                    nw=float(d["nw"]),
                    qe=float(d["qe"]),
                    qg=float(d["qg"]),
                    ma=float(d["ma"]),
                    me=float(d["me"]),
                    narcs=int(d.get("narcs", 3)),
                    **{c: float(v) for c, v in d.get("hertz", {}).items() if c in HERTZ},
                )
            )
        except (KeyError, TypeError) as e:
            raise ValueError(f"roue « {role} » incomplète : {e!r}") from e
    sim = doc.get("simulation", {})
    opti = {**json.loads(json.dumps(OPTIMISATION_DEFAUT)), **doc.get("optimisation", {})}
    return Projet(roues, mu=float(sim.get("mu", 0.2)), mt=float(sim.get("mt", 1.0)), optimisation=opti,
                  da=float(sim.get("da", 0.0)), da_neg=float(sim.get("da_neg", -30.0)),
                  da_pos=float(sim.get("da_pos", 30.0)))


def lire_projet(chemin):
    """.json (format actuel) ou .Cev (ancien format, lecture seule)."""
    return lire_json(chemin) if Path(chemin).suffix.lower() == ".json" else lire_cev(chemin)


def evolve(roue, ni=150):
    """Profil de `roue` : liste de `ni` points (x, y), du pied à la tête.

    Ordre des points (comme EVOLVE.F90, avec nk = ni // 10) :
      - arc de pied de cercle (tangent à l'axe du creux)
      - raccord rectiligne de pied de longueur qg*mo (2 points, si qg > 0)
      - développante (ne = ni - 2*nk points, resserrés autour du primitif)
      - arc de cercle de tête, se terminant sur l'axe de la dent (x = 0)

    `zz < 0` : denture intérieure (couronne). La dent pointe vers le centre, donc le pied est au
    *grand* rayon (sur l'axe du creux) et la tête au *petit* (sur l'axe de la dent) ; `ma` reste
    la hauteur de développante du côté du pied et `me` du côté de la tête. Tout se ramène à un
    signe `sg` : le rayon de départ de la développante est `rt - sg*ma*mo` et celui d'arrivée
    `rt + sg*me*mo`, et la développante est prise sur sa branche miroir (`sg*x`), parce que pour
    une couronne l'angle par rapport à l'axe de la dent doit *décroître* quand le rayon
    décroît, alors que pour une denture extérieure il décroît quand le rayon croît.
    """
    mo, qg = roue.mo, roue.qg
    sg = 1 if roue.zz > 0 else -1  # +1 denture extérieure, -1 couronne intérieure
    zz = abs(roue.zz)
    pi = math.pi
    nk = ni // 10
    ne = ni - 2 * nk
    ev = [[0.0, 0.0] for _ in range(ni)]

    # Paramètres de la développante (rayons toujours positifs, le sens vient de sg)
    rt = mo * zz / 2.0  # rayon primitif
    rb = rt * math.cos(pi * roue.nw / 180.0)  # rayon de base
    wa = math.sqrt((max(rt - sg * roue.ma * mo, rb) / rb) ** 2 - 1.0)  # début (côté pied)
    wt = math.sqrt((rt / rb) ** 2 - 1.0)  # intersection avec le primitif
    we = math.sqrt((max(rt + sg * roue.me * mo, rb) / rb) ** 2 - 1.0)  # fin (côté tête)

    # Développante ; n0 = indice du point sur le primitif (pied / tête)
    n0 = int((ne - 1) * (wt - wa) / (we - wa) + 1)
    for i in range(ne):
        if i <= n0:
            w = wa + i * (wt - wa) / n0
        else:
            w = wt + (i - n0) * (we - wt) / (ne - n0 - 1)
        ev[nk + i] = [
            sg * rb * (math.sin(w) - w * math.cos(w)),
            rb * (math.cos(w) + w * math.sin(w)),
        ]

    # Rotation de be = angle du point primitif + demi-épaisseur de dent
    be = math.atan(ev[nk + n0][0] / ev[nk + n0][1]) + pi / (zz * (1.0 + roue.qe))
    cb, sb = math.cos(be), math.sin(be)
    for i in range(nk, nk + ne):
        x, y = ev[i]
        ev[i] = [x * cb - y * sb, x * sb + y * cb]

    # Arc de cercle de tête (rayon rr, centre (0, ry)), tangent à la fin de la développante et
    # fermé sur l'axe de la dent. L'axe coupe le cercle en deux points, (0, ry+rr) et (0, ry-rr) :
    # la tête est le plus éloigné du centre pour une denture extérieure (w -> 0), le plus proche
    # pour une couronne (w -> pi).
    wd = pi / 2.0 - sg * we + be
    x, y = ev[ni - nk - 1]
    rr = -x / math.sin(wd)  # forme exacte et valable quel que soit le quadrant de wd
    ry = y + x / math.tan(wd)
    wf = 0.0 if sg > 0 else pi
    for i in range(1, nk + 1):
        w = wf + (nk - i) * (wd - wf) / nk
        ev[ni - nk + i - 1] = [-rr * math.sin(w), ry + rr * math.cos(w)]

    # Raccord rectiligne de pied : 2 points à qg*mo/2 et qg*mo du début de la développante
    wb = be - sg * wa  # opposé de l'angle de la tangente au départ de la développante
    if qg > 0:
        x0, y0 = ev[nk]
        for i in (1, 2):
            ev[nk + i - 3] = [
                x0 + sg * math.sin(wb) * qg * mo / i,
                y0 - sg * math.cos(wb) * qg * mo / i,
            ]
    ng = nk - 1 if qg > 0 else nk + 1  # dernier point (1-based) de l'arc de pied

    # Arc de cercle de pied : centre sur l'axe du creux (le profil y est donc perpendiculaire,
    # ce qui le raccorde sans cassure à la dent voisine), tangent au raccord. Le cercle coupe
    # l'axe du creux de part et d'autre du centre ; le pied est du côté du centre de la roue
    # pour une denture extérieure, du côté opposé pour une couronne (w_axe).
    ta = math.tan(pi / 2.0 - pi / zz)
    xg, yg = ev[ng - 1]
    mx = xg - (yg + xg * ta) / (ta + math.tan(wb))
    my = -mx * ta
    rr = math.hypot(mx - xg, my - yg)
    w_axe = -sg * pi / 2.0 + pi / zz
    wd = wb - w_axe
    for i in range(1, ng):
        w = wb - (ng - i) * wd / (ng - 1)
        ev[i - 1] = [mx + rr * math.cos(w), my + rr * math.sin(w)]

    return [(-x, y) for x, y in ev]


def ecrire_egr(chemin, roues, profils):
    """Écrit les profils au format .Egr (sans bloc Hertz), 2 blocs de NI points."""
    with open(chemin, "w", encoding="latin-1", newline="\n") as f:
        for roue, pts in zip(roues, profils):
            f.write(
                f" Developpante {roue.nw:6.2f} {roue.qe:7.2f}{roue.qg:7.2f}"
                f"{roue.ma:7.2f}{roue.me:7.2f}\n"
            )
            f.write(f"{roue.mo:10.6f}{roue.zz:6d}{1.0:10.6f}{0.0:10.6f}\n{len(pts):6d}\n")
            for x, y in pts:
                f.write(f"{x:12.7f}{y:12.7f}\n")


def _rot(pt, a):
    c, s = math.cos(a), math.sin(a)
    return (pt[0] * c - pt[1] * s, pt[0] * s + pt[1] * c)


def contour_dents(pts, zz, indices, angle0, centre, pts_dos=None):
    """Contour de dents complètes d'une roue, en coordonnées du dessin.

    `pts` est le demi-pas de evolve() : flanc droit (x>0) du pied à la tête,
    dent d'axe +y. La dent j est tournée de angle0 + j*2π/|zz| (sens trigonométrique)
    autour de `centre` -- la *valeur absolue* : une dent est parcourue de l'axe du creux d'un
    côté à celui de l'autre, donc les dents ne se raccordent bout à bout que dans ce sens-là,
    quel que soit le signe de `zz` (couronne intérieure). Le contour parcourt chaque dent de droite à gauche (flanc
    droit, puis son miroir) ; les pieds de dents voisines se rejoignent sur l'axe
    du creux. `pts_dos` : demi-profil arrière d'une dent dissymétrique ; le flanc droit
    (x>0, le flanc arrière) est alors dessiné avec lui, le flanc gauche (menant) avec `pts`.
    """
    gauche = [(-x, y) for x, y in reversed(pts)]
    dent = list(pts if pts_dos is None else pts_dos) + gauche[1:]
    contour = []
    for j in indices:
        a = angle0 + j * 2.0 * math.pi / abs(zz)
        contour += [(centre[0] + x, centre[1] + y) for x, y in (_rot(q, a) for q in dent)]
    return contour


def entraxe_signe(roues, da=0.0):
    """Entraxe signé rt1 + rt2 + `da` (mm) ; négatif avec une couronne intérieure (zz < 0), les
    deux centres étant alors du même côté du point primitif.

    `da` s'ajoute algébriquement, et **c'est un écartement au sens du jeu, pas de la distance
    entre axes** : `da > 0` desserre toujours la denture. Sur une denture extérieure les deux
    coïncident (plus d'entraxe = plus de jeu) ; sur une couronne non, puisque éloigner le pignon
    du centre de la couronne l'enfonce plus profond dans les dents et *réduit* le jeu -- avec
    `da > 0` l'entraxe géométrique |rt1 + rt2| diminue donc. C'est la convention dont dépend tout
    le reste : l'optimisation prend -dA comme le cas serré (contrainte `jeu_min`, anti-blocage) et
    +dA comme le cas desserré (`jeu_max`)."""
    return sum(r.mo * r.zz / 2.0 for r in roues) + da


def rayons_fonctionnement(roues, da=0.0):
    """Rayons primitifs de fonctionnement (proportionnels aux nombres de dents) pour l'entraxe
    `entraxe_signe(roues, da)` ; négatif pour une couronne intérieure."""
    r1, r2 = (r.mo * r.zz / 2.0 for r in roues)
    cy = entraxe_signe(roues, da) / (r1 + r2)
    return r1 * cy, r2 * cy


def _placement(roues, n_dents, angle1, angle2, da):
    """[(nombre de dents, indices des dents, angle de la dent 0, centre)] des deux roues,
    et [(centre, rayon primitif nominal)] pour le tracé des cercles primitifs.

    Les centres sont écartés du véritable entraxe rt1 + rt2 + da, avec des rayons de
    fonctionnement (`rayons_fonctionnement`) proportionnels aux nombres de dents : leur somme
    vaut l'entraxe par construction, donc ils sont toujours tangents (nécessaire pour la
    cinématique de l'engrènement). Les cercles primitifs *dessinés* ont eux le rayon nominal
    m*z/2, fixe : à entraxe nominal (da=0) ils sont tangents comme les précédents, mais un
    serrage (da<0) les fait se chevaucher et un écartement (da>0) les sépare, ce que les
    rayons de fonctionnement ne peuvent jamais montrer.
    """
    r1, r2 = rayons_fonctionnement(roues, da)
    r1_nom, r2_nom = (r.mo * r.zz / 2.0 for r in roues)
    # Les dents de chaque roue doivent pointer vers le point primitif. La roue 1 est centrée en
    # (0, -r1) : le point primitif est donc au-dessus d'elle si r1 > 0 (denture extérieure) et
    # en dessous si r1 < 0 (couronne menante), auquel cas son repère est retourné d'un demi-tour
    # -- le même `delta1` que `contact.Engrenement`, sans quoi on dessine les dents du côté opposé
    # et la couronne n'engrène avec rien. Symétriquement pour la roue 2, centrée en (0, r2).
    d1 = math.pi if roues[0].zz < 0 else 0.0
    d2 = math.pi if roues[1].zz > 0 else 0.0
    if angle2 is None:
        # décalée d'un demi-pas, pour qu'une dent de la roue 1 soit en face d'un creux de la roue 2
        angle2 = d2 - math.pi / abs(roues[1].zz)
    idx = range(-(n_dents // 2), n_dents - n_dents // 2)
    c1, c2 = (0.0, -r1), (0.0, r2)
    return ([(roues[0].zz, idx, angle1 + d1, c1), (roues[1].zz, idx, angle2, c2)],
            [(c1, r1_nom), (c2, r2_nom)])


def points_dents(pts, zz, indices, angle0, centre):
    """Points `pts` du demi-pas reportés sur les deux flancs de chaque dent (comme contour_dents)."""
    res = []
    for j in indices:
        a = angle0 + j * 2.0 * math.pi / abs(zz)
        for x, y in pts:
            for px in (x, -x):
                qx, qy = _rot((px, y), a)
                res.append((centre[0] + qx, centre[1] + qy))
    return res


def traits_dents(pts, normales, zz, indices, angle0, centre, demi, pts_dos=None, normales_dos=None):
    """Petits traits (segments) de longueur 2*`demi` en travers du profil à chacun de ses points,
    reportés sur les deux flancs de chaque dent (comme points_dents) ; le flanc droit (arrière)
    utilise `pts_dos` et `normales_dos` s'ils sont donnés (dent dissymétrique)."""
    segs = []
    flancs = ((1.0, pts if pts_dos is None else pts_dos, normales if normales_dos is None else normales_dos),
              (-1.0, pts, normales))
    for j in indices:
        a = angle0 + j * 2.0 * math.pi / abs(zz)
        for sg, pp, nn in flancs:
            for (x, y), (nx, ny) in zip(pp, nn):
                p1 = _rot((sg * (x - demi * nx), y - demi * ny), a)
                p2 = _rot((sg * (x + demi * nx), y + demi * ny), a)
                segs.append(((centre[0] + p1[0], centre[1] + p1[1]), (centre[0] + p2[0], centre[1] + p2[1])))
    return segs


def geometrie_engrenage(roues, profils, n_dents=3, angle1=0.0, angle2=None, da=0.0, profils_dos=None):
    """Deux roues, comme Drawkpo.F90 : roue 1 en bas (axe de dent vers le haut), roue 2 en
    haut tournée de 180°. Le point primitif est en (0,0) et l'entraxe vaut rt1 + rt2 + `da`.
    `angle1`, `angle2` : angle de l'axe de la dent 0 de chaque roue. Par défaut, la
    roue menée est décalée d'un demi-pas : une dent de la roue 1 se trouve en face d'un
    creux de la roue 2 (angle 180° - π/Z2). `profils_dos` : demi-profils arrière (dents dissymétriques).
    """
    placement, cercles = _placement(roues, n_dents, angle1, angle2, da)
    dos = profils_dos or (None, None)
    return [contour_dents(pr, *pl, pts_dos=d) for pr, pl, d in zip(profils, placement, dos)], cercles


PAS_POURCENT = 5.0  # graduation de l'échelle des courbes superposées (%)
MODULES_PAR_PAS = 0.6  # hauteur d'un pas de graduation, en modules (Fortran : 5 % = 20 mm / 300)
COURBES = (  # (clé, nom, couleur) des trois courbes superposées, comme W, C, N du Fortran
    ("rapport", "W  rapport de transmission", "tab:purple"),
    ("couple", "C  rapport de couple", "tab:brown"),
    ("rendement", "N  rendement", "tab:green"),
)


def y_pourcent(pourcent, mo):
    """Ordonnée (mm) d'une valeur en % sur l'échelle superposée : 100 % au point primitif."""
    return (pourcent - 100.0) / PAS_POURCENT * MODULES_PAR_PAS * mo


def points_courbes(superposition, mo):
    """Points (x, y) des courbes superposées, pour dimensionner la fenêtre."""
    if not superposition:
        return []
    xs, series, _ = superposition
    pts = [(0.0, y_pourcent(100.0 + k * PAS_POURCENT, mo)) for k in (-3, 3)]
    for valeurs in series.values():
        pts += [(x, y_pourcent(v, mo)) for x, v in zip(xs, valeurs) if v == v]
    return pts


def limites_engrenage(roues, profils, n_dents, positions, extra=(), da=0.0, profils_dos=None):
    """Fenêtre ((xmin, xmax), (ymin, ymax)) contenant les roues dans toutes les `positions`
    (et les points `extra`)."""
    xs, ys = [p[0] for p in extra], [p[1] for p in extra]
    for pos in positions:
        for cont in geometrie_engrenage(roues, profils, n_dents, pos.alpha1, pos.alpha2, da, profils_dos)[0]:
            xs += [p[0] for p in cont]
            ys += [p[1] for p in cont]
    mx, my = 0.05 * (max(xs) - min(xs)), 0.05 * (max(ys) - min(ys))
    return (min(xs) - mx, max(xs) + mx), (min(ys) - my, max(ys) + my)


def dessiner_engrenage(ax, roues, profils, n_dents=3, position=None, trajet=None, limites=None,
                       superposition=None, da=0.0, profils_arcs=None, points_tangence=None, paires=None,
                       marque_hertz=None, profils_dos=None):
    """Dessine sur `ax` `n_dents` dents de chaque roue, zoom sur la zone d'engrènement.

    `position` (contact.Position) : angles des roues et points de contact ; sinon les roues
    sont presque en prise, sans contact. `trajet` : points de la ligne d'engrènement.
    `limites` : ((xmin, xmax), (ymin, ymax)), sinon ajustées aux roues.
    `superposition` : (abscisses, {clé: valeurs en %}, indice courant) ; les trois courbes
    (COURBES) sont tracées au-dessus des roues à l'abscisse du point d'engrènement, avec
    100 % au point primitif (Draweng.F90). `da` : variation d'entraxe (mm). `profils_arcs` :
    profils approchés par des arcs de cercle (points échantillonnés), tracés en pointillés.
    `points_tangence` : [(points, normales[, points arrière, normales arrière])] par roue, dans le repère du profil ; les points de
    tangence des arcs sont marqués d'un petit trait en travers du profil, sur toutes les dents.
    `profils_dos` : demi-profils arrière des dents dissymétriques (None : dents symétriques).
    `marque_hertz` : ((x, y) sur la menante, (x, y) sur la menée), points du flanc menant où la
    pression de Hertz est maximale (repère de la roue), marqués d'une étoile sur chaque dent. `paires` : paire de dents utilisée à chaque
    position de la simulation, pour interrompre la ligne d'engrènement et les courbes quand on
    saute d'une paire à la suivante.
    """
    a1, a2 = (position.alpha1, position.alpha2) if position else (0.0, None)
    contours, cercles = geometrie_engrenage(roues, profils, n_dents, a1, a2, da, profils_dos)
    couleurs = ("tab:blue", "tab:red")
    noms = ("menante", "menée")
    for k, (cont, roue) in enumerate(zip(contours, roues)):
        ax.plot(*zip(*cont), color=couleurs[k], lw=1.4, label=f"{noms[k]}  Z={roue.zz}")
    if profils_arcs:
        for k, cont in enumerate(geometrie_engrenage(roues, profils_arcs, n_dents, a1, a2, da)[0]):
            ax.plot(*zip(*cont), color="black", lw=0.9, ls=(0, (4, 2)),
                    label="arcs de cercle (export)" if k == 0 else None, zorder=4)
    if points_tangence is not None:
        from matplotlib.collections import LineCollection

        placement = _placement(roues, n_dents, a1, a2, da)[0]
        segs = [sg for pn, pl in zip(points_tangence, placement)
                for sg in traits_dents(pn[0], pn[1], *pl, 0.05 * roues[0].mo, *pn[2:])]
        ax.add_collection(LineCollection(segs, colors="black", linewidths=1.3, zorder=7,
                                         label="points de tangence des arcs"))
    if marque_hertz is not None:
        placement = _placement(roues, n_dents, a1, a2, da)[0]
        etoiles = []
        for (x, y), (zz, idx, angle0, centre) in zip(marque_hertz, placement):
            for j in idx:
                px, py = _rot((x, y), angle0 + j * 2.0 * math.pi / abs(zz))
                etoiles.append((centre[0] + px, centre[1] + py))
        ax.plot(*zip(*etoiles), "*", color="red", mec="black", mew=0.6, ms=13, ls="",
                label="pression de Hertz max", zorder=8)
    # cercles primitifs (arcs autour du point primitif)
    for (c, r), col in zip(cercles, couleurs):
        pas = math.pi * roues[0].mo / r
        a = [math.pi / 2 * (1 if c[1] < 0 else -1) + t * pas * n_dents / 2 * 1.4
             for t in (i / 100 - 0.5 for i in range(101))]
        ax.plot([c[0] + r * math.cos(t) for t in a], [c[1] + r * math.sin(t) for t in a],
                "--", color=col, lw=0.6, alpha=0.6)
    if limites is None:
        limites = limites_engrenage(roues, profils, n_dents, [position or _SansContact()],
                                    points_courbes(superposition, roues[0].mo), da, profils_dos)
    ax.set_xlim(*limites[0])
    ax.set_ylim(*limites[1])
    if trajet is not None:
        _tracer_par_paires(ax, [q[0] for q in trajet], [q[1] for q in trajet], paires,
                           color="tab:orange", lw=2, alpha=0.8, label="ligne d'engrènement")
    if superposition:
        _dessiner_courbes(ax, roues[0].mo, superposition, paires)
    if position:
        ax.plot(*position.contact, "o", color="black", ms=7, label="point d'engrènement", zorder=5)
        if position.contact2:
            ax.plot(*position.contact2, "s", color="white", mec="black", ms=7,
                    label="2e paire en prise", zorder=5)
        titre = "engrènement d'une paire de dents"
    else:
        # la roue menée est décalée d'un demi-pas, depuis 180° si elle est extérieure et depuis 0°
        # si c'est une couronne intérieure (son centre est alors du même côté que celui de la roue 1)
        titre = (f"{n_dents} dents par roue, roue menée décalée de "
                 f"{'180°' if roues[1].zz > 0 else '0°'} - 1/2 pas")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")
    entraxe = f", ΔA = {da * 1000:+g} µm" if da else ""
    angle = f", NW={roues[0].nw}°" if roues[0].nw is not None else ""
    ax.set_title(f"{titre}\n(m={roues[0].mo}{angle}{entraxe})")


def _tracer_par_paires(ax, xs, ys, paires, **style):
    """Trace xs, ys ; quand la paire de dents change, le saut à la paire voisine est relié par un
    trait pointillé fin (comme Draweng.F90)."""
    if paires is None:
        ax.plot(xs, ys, **style)
        return
    label = style.pop("label", None)
    debut, premier = 0, True
    for i in range(1, len(xs) + 1):
        if i < len(xs) and paires[i] == paires[debut]:
            continue
        ax.plot(xs[debut:i], ys[debut:i], label=label if premier else None, **style)
        premier = False
        if i < len(xs):
            ax.plot(xs[i - 1:i + 1], ys[i - 1:i + 1], color=style.get("color"), lw=0.8, ls=":",
                    alpha=0.7, zorder=style.get("zorder"))
        debut = i


def _dessiner_courbes(ax, mo, superposition, paires=None):
    """Axe gradué (100 % au point primitif, sur la ligne des centres) et courbes W, C, N."""
    xs, series, i = superposition
    tic = 0.15 * mo
    ax.plot([0, 0], [y_pourcent(85.0, mo), y_pourcent(115.0, mo)], color="black", lw=0.8, zorder=3)
    for k in range(-3, 4):
        y = y_pourcent(100.0 + k * PAS_POURCENT, mo)
        ax.plot([-tic, tic], [y, y], color="black", lw=0.8, zorder=3)
        ax.text(1.6 * tic, y, f"{100 + k * PAS_POURCENT:g} %", fontsize=7, va="center",
                fontweight="bold" if k == 0 else "normal", zorder=3)
    for cle, nom, couleur in COURBES:
        ys = [y_pourcent(v, mo) for v in series[cle]]
        _tracer_par_paires(ax, xs, ys, paires, color=couleur, lw=1.6, label=nom, zorder=4)
        if ys[i] == ys[i]:
            ax.plot(xs[i], ys[i], "o", color=couleur, ms=5, zorder=6)


class SansContact:
    """Position sans contact : roues presque en prise (pour limites_engrenage)."""

    alpha1, alpha2 = 0.0, None


_SansContact = SansContact


def tracer(roues, profils, n_dents=3, png=None, afficher=True):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 7))
    dessiner_engrenage(ax, roues, profils, n_dents)
    if png:
        fig.savefig(png, dpi=150, bbox_inches="tight")
    if afficher:
        plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("fichier", nargs="?", default="test.Cev", help=".json ou ancien .Cev")
    p.add_argument("--ni", type=int, default=150, help="points par profil")
    p.add_argument("--plot", action="store_true", help="affiche les dents en prise")
    p.add_argument("--png", help="enregistre le dessin dans ce fichier")
    p.add_argument("--dents", type=int, default=3, help="dents dessinées par roue")
    a = p.parse_args()

    roues = lire_projet(a.fichier).roues
    profils = [evolve(r, a.ni) for r in roues]

    sortie = Path("Result") / (Path(a.fichier).stem + ".Egr")
    sortie.parent.mkdir(exist_ok=True)
    ecrire_egr(sortie, roues, profils)
    print(f"{len(profils)} profils de {a.ni} points -> {sortie}")
    if a.plot or a.png:
        tracer(roues, profils, a.dents, a.png, a.plot)


if __name__ == "__main__":
    sys.exit(main())
