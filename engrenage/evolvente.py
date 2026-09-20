"""Génération de profils d'engrenage en développante de cercle.

Port de Fortran/Egc/Evolve.F90 (K. Trümpy). Lit un fichier .Cev et calcule,
pour chaque mobile, un demi-pas de denture (flanc + pieds + têtes) sous forme
d'une suite de NI points, origine (0,0) au centre du mobile.

Usage : python evolvente.py [fichier.Cev] [--ni 150] [--plot]
"""

import argparse
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


def lire_cev(chemin):
    """Lit les deux roues (menante, menée) d'un fichier .Cev.

    Ligne 1 : en-têtes. Lignes 2 et 3 : MO ZZ NW QE QG MA ME.
    Les lignes suivantes (Hertz) ne sont pas nécessaires ici.
    """
    lignes = Path(chemin).read_text(encoding="latin-1").splitlines()
    roues = []
    for ligne in lignes[1:3]:
        v = ligne.split()
        roues.append(
            Roue(
                mo=float(v[0]),
                zz=int(v[1]),
                nw=float(v[2]),
                qe=float(v[3]),
                qg=float(v[4]),
                ma=float(v[5]),
                me=float(v[6]),
            )
        )
    return roues


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


def tracer(roues, profils):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, roue, pts in zip(axes, roues, profils):
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ".-", ms=3)
        ax.plot(0, 0, "k+")
        ax.set_aspect("equal")
        ax.set_title(f"m={roue.mo}  Z={roue.zz}  NW={roue.nw}°")
        ax.grid(True, alpha=0.3)
    plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("fichier", nargs="?", default="test.Cev")
    p.add_argument("--ni", type=int, default=150, help="points par profil")
    p.add_argument("--plot", action="store_true")
    a = p.parse_args()

    roues = lire_cev(a.fichier)
    profils = [evolve(r, a.ni) for r in roues]

    sortie = Path("Result") / (Path(a.fichier).stem + ".Egr")
    sortie.parent.mkdir(exist_ok=True)
    ecrire_egr(sortie, roues, profils)
    print(f"{len(profils)} profils de {a.ni} points -> {sortie}")
    if a.plot:
        tracer(roues, profils)


if __name__ == "__main__":
    sys.exit(main())
