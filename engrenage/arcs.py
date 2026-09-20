"""Approximation d'un profil par des arcs de cercle tangents, et export .Eno.

Repris de Krebax.F90 et Lagrmo.F90 (K. Trümpy) : le demi-pas du profil (points de
evolvente.evolve, `ZP` du Fortran) est remplacé par une suite d'arcs de cercle qui se
raccordent tangentiellement. Les arcs sont définis par leurs points de transition ; les
rayons et les centres s'en déduisent (`krbgs`) à partir de l'angle de la tangente au
départ (-pi/Z, dans le repère du profil). Les points de transition sont optimisés
(multiplicateur de Lagrange) pour minimiser la somme des carrés des distances du profil aux
arcs, avec une tangente horizontale (angle 0) à la fin, sur l'axe de la dent.
"""

import math
from dataclasses import dataclass

import numpy as np

NARCS_MAX = 9  # arcs sur la développante (MK-3 dans Engcre.F90, MK=12)


@dataclass
class ProfilArcs:
    points: np.ndarray  # (n,2) points de transition, dans le repère du profil (KR du Fortran)
    centres: np.ndarray  # (n-1,2) centres des arcs
    rayons: np.ndarray  # (n-1,) rayons signés
    ecart_max: float  # distance maximale d'un point du profil aux arcs (mm)
    indice_ecart: int  # indice du profil où elle est atteinte
    erreur_angle: float  # écart de la tangente finale par rapport à l'axe de la dent (rad)

    @property
    def n_points(self):
        return len(self.points)


def krbgs(a1, kr):
    """Arcs tangents passant par les points `kr`, pour l'angle a1 de la tangente au départ.

    Retourne (centres, rayons signés, angle de la tangente à la fin) : KRBGS du Fortran.
    """
    co, si = math.cos(a1), math.sin(a1)
    n = len(kr) - 1
    pm, ra = np.zeros((n, 2)), np.zeros(n)
    for i in range(n):
        dx, dy = kr[i + 1, 0] - kr[i, 0], kr[i + 1, 1] - kr[i, 1]
        rr = 0.5 * (dx * dx + dy * dy) / (dx * si - dy * co)
        pm[i] = kr[i, 0] + rr * si, kr[i, 1] - rr * co
        vx, vy = kr[i + 1, 0] - pm[i, 0], kr[i + 1, 1] - pm[i, 1]
        ra[i] = math.copysign(rr, vx * dy - vy * dx)
        co, si = vy / abs(rr), -vx / abs(rr)
    return pm, ra, math.atan2(-si, co)


def sqab(xy, kr, pm, ra):
    """Somme des carrés des distances des points `xy` à leur arc, distance maximale, indice.

    Les points sont pris dans l'ordre : on passe à l'arc suivant quand le point sort de
    l'arc courant (SQAB du Fortran).
    """
    a = xy[None, :, :] - pm[:, None, :]  # (arcs, points, 2)
    c1 = kr[:-1, None, :] - pm[:, None, :]
    c2 = kr[1:, None, :] - pm[:, None, :]

    def angle(c):
        return np.arctan2(a[..., 0] * c[..., 1] - a[..., 1] * c[..., 0],
                          a[..., 0] * c[..., 0] + a[..., 1] * c[..., 1])

    dehors = angle(c1) * angle(c2) > 0
    ecart = np.abs(np.abs(ra)[:, None] - np.hypot(a[..., 0], a[..., 1]))
    n_arcs = len(pm)
    arc, somme, maxi, imax = 0, 0.0, 0.0, 0
    for i in range(len(xy)):
        if dehors[arc, i]:
            arc = min(arc + 1, n_arcs - 1)
        d = ecart[arc, i]
        somme += d * d
        if d > maxi:
            maxi, imax = d, i
    return somme, maxi, imax


def lagrmo(aw, xy, kp):
    """Optimise les points de transition `kp` (les deux extrémités restent fixes) : LAGRMO."""
    nk = len(kp)
    dt = 1e-4  # variation relative des coordonnées de transition
    rx = nk - 2  # relaxation
    mn = 2 * nk - 3
    kr = kp.copy()
    dt *= math.hypot(kr[-1, 0] - kr[0, 0], kr[-1, 1] - kr[0, 1])
    meilleur = (math.inf, 0, 0.0)  # (ecart max, indice, erreur d'angle)
    for _ in range(100):
        pm, ra, a0 = krbgs(aw[0], kr)
        f0 = sqab(xy, kr, pm, ra)[0]
        kk = kr.copy()
        mt, bs = np.zeros((mn, mn)), np.zeros(mn)
        for n in range(1, mn):  # inconnues : x, y des points de transition intérieurs
            j, i = (n - 1) % 2, 1 + (n - 1) // 2
            f, a = {}, {}
            for s in (-1, 1):
                kk[i, j] = kr[i, j] + s * dt
                pm_s, ra_s, a[s] = krbgs(aw[0], kk)
                f[s] = sqab(xy, kk, pm_s, ra_s)[0]
            kk[i, j] = kr[i, j]
            bs[n - 1] = (f[-1] - f[1]) / (2 * dt)
            mt[n - 1, n - 1] = rx * (f[-1] + f[1] - 2 * f0) / (dt * dt)
            mt[n - 1, mn - 1] = mt[mn - 1, n - 1] = (a[-1] - a[1]) / (2 * dt)
        bs[mn - 1] = a0 - aw[1]
        mt[mn - 1, mn - 1] = 0.0
        try:
            dx = np.linalg.solve(mt, bs)
        except np.linalg.LinAlgError:
            dx = None
        if dx is None:  # matrice singulière
            break
        for n in range(1, mn):
            kr[1 + (n - 1) // 2, (n - 1) % 2] += dx[n - 1]
        pm, ra, a0 = krbgs(aw[0], kr)
        _, mr1, mr2 = sqab(xy, kr, pm, ra)
        if mr1 > meilleur[0]:  # l'écart maximal remonte : on garde le meilleur
            break
        meilleur = (mr1, mr2, aw[1] - a0)
        kp[:] = kr
    return meilleur


def approximer(profil, zz, qg, narcs):
    """Profil (liste de points du demi-pas) remplacé par des arcs de cercle tangents.

    `narcs` arcs sur la développante, plus l'arc de pied, l'arc de tête et, si `qg > 0`, le
    raccord de pied (comme KREBAX : NE(K)+3 arcs, +1 avec raccord). Les points de départ des
    arcs sont répartis régulièrement sur les points de la développante.
    """
    if not 1 <= narcs <= NARCS_MAX:
        raise ValueError(f"nombre d'arcs entre 1 et {NARCS_MAX} attendu")
    xy = np.asarray(profil, dtype=float)
    ni = len(xy)
    nk2 = ni // 10
    m = 1 if qg > 0 else 0  # raccord rectiligne de pied
    nk1 = nk2 - 2 * m
    kp = [1, nk1 + 1] + ([nk2 + 1] if m else [])  # indices (base 1) des points de transition
    kp += [nk2 + j * (ni - 2 * nk2) // narcs for j in range(1, narcs + 1)] + [ni]
    kr = xy[[k - 1 for k in kp]].copy()
    aw = (-math.pi / zz, 0.0)  # angles de la tangente au départ et à la fin
    ecart, indice, erreur = lagrmo(aw, xy, kr)
    pm, ra, _ = krbgs(aw[0], kr)
    return ProfilArcs(kr, pm, ra, ecart, indice, erreur)


def echantillonner(arcs, par_arc=40):
    """Points (n,2) répartis sur les arcs, du début à la fin du profil."""
    pts = []
    kr = arcs.points
    for i in range(len(arcs.rayons)):
        c, p, q = arcs.centres[i], kr[i], kr[i + 1]
        a0 = math.atan2(p[1] - c[1], p[0] - c[0])
        a1 = math.atan2(q[1] - c[1], q[0] - c[0])
        da = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi  # angle balayé, le plus court
        r = abs(arcs.rayons[i])
        t = np.linspace(0.0, 1.0, par_arc, endpoint=(i == len(arcs.rayons) - 1))
        pts.append(np.column_stack((c[0] + r * np.cos(a0 + t * da), c[1] + r * np.sin(a0 + t * da))))
    return np.vstack(pts)


def ecrire_eno(chemin, roues, arcs):
    """Écrit les profils en arcs de cercle au format .Eno de ENGCRE (Engcre.F90).

    Par roue : nom, en-tête, `module ZZ nombre de points`, puis les points de transition
    (x, y, 0 aux extrémités, 1 entre les deux) ; ensuite les données de Hertz des deux roues
    et la ligne des paramètres de simulation. Colonnes fixes, fins de ligne CRLF.
    """
    lignes = []
    for role, roue, a in zip(("m", "e"), roues, arcs):
        n = a.n_points
        lignes += [f"Developpante-{role}", "   Module |      ZZ |   Nbr.point|",
                   f"{roue.mo:10.6f}{roue.zz:10d}{n:10d}"]
        for i, (x, y) in enumerate(a.points):
            lignes.append(f"{x:12.7f}{y:12.7f}{0 if i in (0, n - 1) else 1:6d}")
    lignes.append("   E-Mod.  Coeff.Poiss. Epaiss.")
    for roue in roues:
        if None in (roue.e_mod, roue.poisson, roue.epaisseur):
            raise ValueError("données de Hertz manquantes (E-module, Poisson, épaisseur)")
        lignes.append(f"{roue.e_mod:10.1f}{roue.poisson:10.3f}{roue.epaisseur:10.3f}")
    lignes += ["    -dA   |   +dA   |    GW   |    ZV   | ", "     -30.0      30.0     120          1"]
    with open(chemin, "w", encoding="ascii", newline="") as f:
        f.write("\r\n".join(lignes) + "\r\n")
