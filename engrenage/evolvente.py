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
from dataclasses import dataclass
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


@dataclass
class Projet:
    roues: list  # [menante, menée]
    mu: float = 0.2  # coefficient de frottement de la simulation
    mt: float = 1.0  # couple sur la menante (mN.m)
    da: float = 0.0  # variation de l'entraxe (µm)


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
    doc = {"format": FORMAT_JSON, "version": 1, "roues": roues, "simulation": {"mu": projet.mu, "mt": projet.mt, "da": projet.da}}
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
    return Projet(roues, mu=float(sim.get("mu", 0.2)), mt=float(sim.get("mt", 1.0)),
                  da=float(sim.get("da", 0.0)))


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
    """
    mo, zz, qg = roue.mo, roue.zz, roue.qg
    pi = math.pi
    nk = ni // 10
    ne = ni - 2 * nk
    ev = [[0.0, 0.0] for _ in range(ni)]

    # Paramètres de la développante
    rt = mo * zz / 2.0  # rayon primitif
    rb = rt * math.cos(pi * roue.nw / 180.0)  # rayon de base
    wa = math.sqrt((max(rt - roue.ma * mo, rb) / rb) ** 2 - 1.0)  # début
    wt = math.sqrt((rt / rb) ** 2 - 1.0)  # intersection avec le primitif
    we = math.sqrt(((rt + roue.me * mo) / rb) ** 2 - 1.0)  # fin

    # Développante ; n0 = indice du point sur le primitif (pied / tête)
    n0 = int((ne - 1) * (wt - wa) / (we - wa) + 1)
    for i in range(ne):
        if i <= n0:
            w = wa + i * (wt - wa) / n0
        else:
            w = wt + (i - n0) * (we - wt) / (ne - n0 - 1)
        ev[nk + i] = [
            rb * (math.sin(w) - w * math.cos(w)),
            rb * (math.cos(w) + w * math.sin(w)),
        ]

    # Rotation de be = angle du point primitif + demi-épaisseur de dent
    be = math.atan(ev[nk + n0][0] / ev[nk + n0][1]) + pi / (zz * (1.0 + roue.qe))
    cb, sb = math.cos(be), math.sin(be)
    for i in range(nk, nk + ne):
        x, y = ev[i]
        ev[i] = [x * cb - y * sb, x * sb + y * cb]

    # Arc de cercle de tête (rayon rr, centre (0, ry))
    wd = pi / 2.0 - we + be
    ta = math.tan(wd)
    x, y = ev[ni - nk - 1]
    rr = -x * math.sqrt(1.0 + 1.0 / (ta * ta))
    ry = y + x / ta
    for i in range(1, nk + 1):
        w = (nk - i) * wd / nk
        ev[ni - nk + i - 1] = [-rr * math.sin(w), ry + rr * math.cos(w)]

    # Raccord rectiligne de pied : 2 points à qg*mo/2 et qg*mo du début de la développante
    wb = be - wa
    if qg > 0:
        x0, y0 = ev[nk]
        for i in (1, 2):
            ev[nk + i - 3] = [
                x0 + math.sin(wb) * qg * mo / i,
                y0 - math.cos(wb) * qg * mo / i,
            ]
    ng = nk - 1 if qg > 0 else nk + 1  # dernier point (1-based) de l'arc de pied

    # Arc de cercle de pied, tangent au raccord et à l'axe du creux
    ta = math.tan(pi / 2.0 - pi / zz)
    xg, yg = ev[ng - 1]
    mx = xg - (yg + xg * ta) / (ta + math.tan(wb))
    my = -mx * ta
    rr = math.hypot(mx - xg, my - yg)
    wd = wb + pi * (0.5 - 1.0 / zz)
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


def contour_dents(pts, zz, indices, angle0, centre):
    """Contour de dents complètes d'une roue, en coordonnées du dessin.

    `pts` est le demi-pas de evolve() : flanc droit (x>0) du pied à la tête,
    dent d'axe +y. La dent j est tournée de angle0 + j*2π/zz (sens trigonométrique)
    autour de `centre`. Le contour parcourt chaque dent de droite à gauche (flanc
    droit, puis son miroir) ; les pieds de dents voisines se rejoignent sur l'axe
    du creux.
    """
    gauche = [(-x, y) for x, y in reversed(pts)]
    dent = pts + gauche[1:]
    contour = []
    for j in indices:
        a = angle0 + j * 2.0 * math.pi / zz
        contour += [(centre[0] + x, centre[1] + y) for x, y in (_rot(q, a) for q in dent)]
    return contour


def rayons_fonctionnement(roues, da=0.0):
    """Rayons des cercles primitifs de fonctionnement pour un entraxe rt1 + rt2 + da (mm)."""
    r1, r2 = (r.mo * r.zz / 2.0 for r in roues)
    cy = (r1 + r2 + da) / (r1 + r2)
    return r1 * cy, r2 * cy


def geometrie_engrenage(roues, profils, n_dents=3, angle1=0.0, angle2=None, da=0.0):
    """Deux roues, comme Drawkpo.F90 : roue 1 en bas (axe de dent vers le haut), roue 2 en
    haut tournée de 180°. Le point primitif est en (0,0) et l'entraxe vaut rt1 + rt2 + `da`.
    `angle1`, `angle2` : angle de l'axe de la dent 0 de chaque roue. Par défaut, la
    roue menée est décalée d'un demi-pas : une dent de la roue 1 se trouve en face d'un
    creux de la roue 2 (angle 180° - π/Z2).
    """
    r1, r2 = rayons_fonctionnement(roues, da)
    z1, z2 = roues[0].zz, roues[1].zz
    if angle2 is None:
        angle2 = math.pi - math.pi / z2
    idx = range(-(n_dents // 2), n_dents - n_dents // 2)
    c1, c2 = (0.0, -r1), (0.0, r2)
    return [
        contour_dents(profils[0], z1, idx, angle1, c1),
        contour_dents(profils[1], z2, idx, angle2, c2),
    ], [(c1, r1), (c2, r2)]


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


def limites_engrenage(roues, profils, n_dents, positions, extra=(), da=0.0):
    """Fenêtre ((xmin, xmax), (ymin, ymax)) contenant les roues dans toutes les `positions`
    (et les points `extra`)."""
    xs, ys = [p[0] for p in extra], [p[1] for p in extra]
    for pos in positions:
        for cont in geometrie_engrenage(roues, profils, n_dents, pos.alpha1, pos.alpha2, da)[0]:
            xs += [p[0] for p in cont]
            ys += [p[1] for p in cont]
    mx, my = 0.05 * (max(xs) - min(xs)), 0.05 * (max(ys) - min(ys))
    return (min(xs) - mx, max(xs) + mx), (min(ys) - my, max(ys) + my)


def dessiner_engrenage(ax, roues, profils, n_dents=3, position=None, trajet=None, limites=None,
                       superposition=None, da=0.0, profils_arcs=None):
    """Dessine sur `ax` `n_dents` dents de chaque roue, zoom sur la zone d'engrènement.

    `position` (contact.Position) : angles des roues et points de contact ; sinon les roues
    sont presque en prise, sans contact. `trajet` : points de la ligne d'engrènement.
    `limites` : ((xmin, xmax), (ymin, ymax)), sinon ajustées aux roues.
    `superposition` : (abscisses, {clé: valeurs en %}, indice courant) ; les trois courbes
    (COURBES) sont tracées au-dessus des roues à l'abscisse du point d'engrènement, avec
    100 % au point primitif (Draweng.F90). `da` : variation d'entraxe (mm). `profils_arcs` :
    profils approchés par des arcs de cercle (points échantillonnés), tracés en pointillés.
    """
    a1, a2 = (position.alpha1, position.alpha2) if position else (0.0, None)
    contours, cercles = geometrie_engrenage(roues, profils, n_dents, a1, a2, da)
    couleurs = ("tab:blue", "tab:red")
    noms = ("menante", "menée")
    for k, (cont, roue) in enumerate(zip(contours, roues)):
        ax.plot(*zip(*cont), color=couleurs[k], lw=1.4, label=f"{noms[k]}  Z={roue.zz}")
    if profils_arcs:
        for k, cont in enumerate(geometrie_engrenage(roues, profils_arcs, n_dents, a1, a2, da)[0]):
            ax.plot(*zip(*cont), color="black", lw=0.9, ls=(0, (4, 2)),
                    label="arcs de cercle (export)" if k == 0 else None, zorder=4)
    # cercles primitifs (arcs autour du point primitif)
    for (c, r), col in zip(cercles, couleurs):
        pas = math.pi * roues[0].mo / r
        a = [math.pi / 2 * (1 if c[1] < 0 else -1) + t * pas * n_dents / 2 * 1.4
             for t in (i / 100 - 0.5 for i in range(101))]
        ax.plot([c[0] + r * math.cos(t) for t in a], [c[1] + r * math.sin(t) for t in a],
                "--", color=col, lw=0.6, alpha=0.6)
    if limites is None:
        limites = limites_engrenage(roues, profils, n_dents, [position or _SansContact()],
                                    points_courbes(superposition, roues[0].mo), da)
    ax.set_xlim(*limites[0])
    ax.set_ylim(*limites[1])
    if trajet is not None:
        ax.plot(*zip(*trajet), color="tab:orange", lw=2, alpha=0.8, label="ligne d'engrènement")
    if superposition:
        _dessiner_courbes(ax, roues[0].mo, superposition)
    if position:
        ax.plot(*position.contact, "o", color="black", ms=7, label="point d'engrènement", zorder=5)
        if position.contact2:
            ax.plot(*position.contact2, "s", color="white", mec="black", ms=7,
                    label="2e paire en prise", zorder=5)
        titre = "engrènement d'une paire de dents"
    else:
        titre = f"{n_dents} dents par roue, roue menée décalée de 180° - 1/2 pas"
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")
    entraxe = f", ΔA = {da * 1000:+g} µm" if da else ""
    ax.set_title(f"{titre}\n(m={roues[0].mo}, NW={roues[0].nw}°{entraxe})")


def _dessiner_courbes(ax, mo, superposition):
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
        ax.plot(xs, ys, color=couleur, lw=1.6, label=nom, zorder=4)
        if ys[i] == ys[i]:
            ax.plot(xs[i], ys[i], "o", color=couleur, ms=5, zorder=6)


class _SansContact:
    alpha1, alpha2 = 0.0, None


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
