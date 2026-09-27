"""Optimisation des profils en arcs de cercle : noyau rapide et fonction objectif.

Repris de Krbgs.F90 (KRBGS, KORZF, SWTE), Kontac.F90, Enganac.F90 et Resulto.F90 de ENGOPT
(K. Trümpy), écrits pour numba : une évaluation complète (3 entraxes, 101 positions) prend
environ 3 ms par cœur et les lots de profils se calculent en parallèle. Mêmes conventions que
contact.py (repère du dessin : point primitif en (0,0), menante de centre (0,-r1), menée de
centre (0,r2), flanc menant = profil en miroir) ; les résultats sont ceux de
`contact.Engrenement(..., mobiles=...)` à 1e-3 % près (voir `python optim.py`).

Variables d'optimisation (`Probleme`) : courbure et longueur de chaque arc, pas les
coordonnées des points de transition. Le profil est une chaîne d'arcs tangents : déplacer un
seul point de 1 µm fait tourner la tangente suivante de plusieurs degrés, et l'on ne peut donc
guère s'écarter de la solution de départ (1 à 3 µm), alors qu'avec la courbure et la longueur
(intégrées depuis le premier point) on reste dans le domaine admissible jusqu'à 10 à 30 µm.
Le dernier arc est refermé comme KORZF : son centre est sur l'axe de la dent, ce qui fixe la
pente finale à 0 et le point de tête (le rayon de tête n'est donc pas une variable).

Dents dissymétriques : chaque roue a un demi-profil menant et un demi-profil arrière (le même si
la dent est symétrique). Le premier sert au contact sous charge (couple, Hertz, rendement), le
second au jeu (contact des flancs arrière) ; les deux servent au contrôle d'interférence, qui
cherche les croisements d'arcs entre les contours complets des dents des deux roues (dents
voisines de la ligne des centres, 6 positions par période, à chaque entraxe) : un contact
tangent n'est pas une interférence, un recouvrement de plus de `EPS_INTERF` en est une.
"""

import math
from dataclasses import dataclass

import numpy as np
from numba import njit, prange

import arcs

TWO_PI = 2.0 * math.pi
IK = 101  # positions de la menante sur une période (IK d'Enganac pour NJ = 0)
N_JEU = 12  # positions où le jeu est évalué pendant l'optimisation
PENALITE = 1e3  # valeur de l'objectif d'un profil invalide (1D33 dans Resulto)

# colonnes de la sortie de `evaluer` (une ligne par entraxe)
VARC, REND, HERTZ, SYML, REND_MIN, JEU_MAX, JEU_MIN, COUPLE_MOY, INTERF, R_TETE, ECART_TETE = range(11)
NCOL = 11
EPS_INTERF = 1e-3  # mm : recouvrement de deux arcs au-dessous duquel on parle de contact, pas d'interférence
MESSAGES = {0: "ok", 1: "entraxe trop petit", 2: "entraxe trop grand", 3: "pas de période d'engrènement",
            4: "arc-boutement",
            5: "profil invalide (rayons de courbure) : essayez un autre nombre d'arcs dans l'onglet Création"}


# ---------------------------------------------------------------------------------------
# noyau numba
# ---------------------------------------------------------------------------------------

@njit(cache=True)
def _balayage(a):
    return (a + math.pi) % TWO_PI - math.pi


@njit(cache=True)
def krbgs(a0, xy, pm, ra):
    """KRBGS : centres `pm` et rayons signés `ra` des n-1 arcs tangents passant par les points
    `xy`, pente de départ a0. Retourne la pente de départ du dernier arc (A1)."""
    n = xy.shape[0]
    co, si = math.cos(a0), math.sin(a0)
    a1 = 0.0
    for i in range(n - 1):
        dx, dy = xy[i + 1, 0] - xy[i, 0], xy[i + 1, 1] - xy[i, 1]
        den = dx * si - dy * co
        if den == 0.0:
            den = 1e-300
        rr = 0.5 * (dx * dx + dy * dy) / den
        pm[i, 0], pm[i, 1] = xy[i, 0] + rr * si, xy[i, 1] - rr * co
        vx, vy = xy[i + 1, 0] - pm[i, 0], xy[i + 1, 1] - pm[i, 1]
        ra[i] = abs(rr) if vx * dy - vy * dx >= 0.0 else -abs(rr)
        co, si = vy / abs(rr), -vx / abs(rr)
        if i == n - 3:
            a1 = si / co
        if i > 0:  # arc de plus de 180 degrés : le signe du rayon suit celui des centres voisins
            dot = ((pm[i - 1, 0] - xy[i, 0]) * (pm[i, 0] - xy[i, 0])
                   + (pm[i - 1, 1] - xy[i, 1]) * (pm[i, 1] - xy[i, 1]))
            if (dot >= 0.0) != (ra[i] * ra[i - 1] >= 0.0):
                ra[i] = -ra[i]
    return a1


@njit(cache=True)
def korzf(a0, xy, pm, ra):
    """KORZF : comme `krbgs`, mais le dernier arc est refermé (centre sur l'axe de la dent,
    pente finale nulle) ; le dernier point (la tête) est recalculé dans `xy`."""
    n = xy.shape[0]
    a1 = krbgs(a0, xy, pm, ra)
    if a1 == 0.0:
        a1 = 1e-300
    dx = xy[n - 2, 0]
    dy = dx / a1
    pm[n - 2, 1] = xy[n - 2, 1] + dy
    pm[n - 2, 0] = 0.0
    rr = math.sqrt(dx * dx + dy * dy)
    ra[n - 2] = rr if ra[n - 2] >= 0.0 else -rr
    xy[n - 1, 1] = pm[n - 2, 1] + rr
    xy[n - 1, 0] = 0.0


@njit(cache=True)
def _swte(xy, ra, rk):
    """SWTE : 1 si le profil est valide (un seul changement de signe des rayons, rayon de tête
    au moins rk fois celui du pied), sinon une autre valeur."""
    n = xy.shape[0]
    s = 0
    for i in range(n - 2):
        if ra[i] * ra[i + 1] < 0.0:
            s += 1
    if xy[n - 2, 0] < rk * xy[0, 0]:
        s = 2
    return s


@njit(cache=True)
def _flanc(pts, pm, ra, s):
    """Flanc d'un demi-profil : centres, rayons, convexité, angle de départ et balayage de chaque
    arc, dans le repère de la roue. s = -1 : image miroir (flanc menant, x < 0) ; s = +1 : tel quel."""
    n = pts.shape[0] - 1
    c, r, a0, sw = np.empty((n, 2)), np.empty(n), np.empty(n), np.empty(n)
    cv = np.empty(n, dtype=np.bool_)
    for i in range(n):
        c[i, 0], c[i, 1] = s * pm[i, 0], pm[i, 1]
        r[i] = abs(ra[i])
        cv[i] = ra[i] > 0.0
        p0 = math.atan2(pts[i, 1] - c[i, 1], s * pts[i, 0] - c[i, 0])
        p1 = math.atan2(pts[i + 1, 1] - c[i, 1], s * pts[i + 1, 0] - c[i, 0])
        a0[i] = p0
        sw[i] = _balayage(p1 - p0)
    return c, r, cv, a0, sw


@njit(cache=True)
def _dans_arc(theta, a0, sw):
    t = _balayage(theta - a0)
    if sw < 0.0:
        t = -t
    elif sw == 0.0:
        t = 0.0
    return -1e-9 <= t <= abs(sw) + 1e-9


@njit(cache=True)
def _contact(phi, f1, f2, o1, o2, out):
    """KONTAC : plus petite rotation de la menée (autour de o2) qui rend tangents un arc du
    flanc menant tourné de phi et un arc du flanc mené (inf si aucune). Remplit
    out = [x, y, tangente x, tangente y, rayon 1, rayon 2] (rayons signés, convexe > 0)."""
    c1, r1, v1, a01, sw1 = f1
    c2, r2, v2, a02s, sw2 = f2
    cph, sph = math.cos(phi), math.sin(phi)
    best = math.inf
    for m in range(c1.shape[0]):
        cx = o1[0] + cph * c1[m, 0] - sph * c1[m, 1]
        cy = o1[1] + sph * c1[m, 0] + cph * c1[m, 1]
        a0m = a01[m] + phi
        wx, wy = cx - o2[0], cy - o2[1]
        w = math.hypot(wx, wy)
        aw = math.atan2(wy, wx)
        for n in range(c2.shape[0]):
            if v1[m] and v2[n]:
                dt = r1[m] + r2[n]
            elif v1[m]:
                dt = r2[n] - r1[m]
            elif v2[n]:
                dt = r1[m] - r2[n]
            else:
                continue  # deux arcs concaves ne se touchent pas
            if dt <= 1e-12:
                continue
            rho = math.hypot(c2[n, 0], c2[n, 1])
            if rho * w < 1e-12:
                continue
            h = (rho * rho + w * w - dt * dt) / (2.0 * rho * w)
            if abs(h) > 1.0 + 1e-12:
                continue
            ac = math.acos(min(1.0, max(-1.0, h)))
            a02 = math.atan2(c2[n, 1], c2[n, 0])
            for sg in (1.0, -1.0):
                d = (aw - a02 + sg * ac) % TWO_PI
                if d >= best:
                    continue
                c2x, c2y = o2[0] + rho * math.cos(a02 + d), o2[1] + rho * math.sin(a02 + d)
                ux, uy = (c2x - cx) / dt, (c2y - cy) / dt
                s = -1.0 if (v1[m] and not v2[n]) else 1.0
                px, py = cx + s * r1[m] * ux, cy + s * r1[m] * uy
                th1 = math.atan2(py - cy, px - cx)
                th2 = math.atan2(py - c2y, px - c2x) - d
                if _dans_arc(th1, a0m, sw1[m]) and _dans_arc(th2, a02s[n], sw2[n]):
                    best = d
                    out[0], out[1], out[2], out[3] = px, py, -uy, ux
                    out[4] = r1[m] if v1[m] else -r1[m]
                    out[5] = r2[n] if v2[n] else -r2[n]
    return best


@njit(cache=True)
def _efforts(o1, o2, px, py, gx, gy, mu, res):
    """res = [rapport w2/w1, rapport de couple, rendement, force normale] pour un couple moteur
    de 1 (voir contact.efforts) ; nan en cas d'arc-boutement."""
    r1x, r1y, r2x, r2y = px - o1[0], py - o1[1], px - o2[0], py - o2[1]
    nx, ny = -gy, gx
    p1x, p1y, p2x, p2y = -r1y, r1x, -r2y, r2x
    rho = (nx * p1x + ny * p1y) / (nx * p2x + ny * p2y)
    gl = (rho * p2x - p1x) * gx + (rho * p2y - p1y) * gy
    if r1x * ny - r1y * nx < 0.0:
        nx, ny = -nx, -ny
    sg = 1.0 if gl > 0.0 else (-1.0 if gl < 0.0 else 0.0)
    dx, dy = nx - mu * sg * gx, ny - mu * sg * gy
    m1 = r1x * dy - r1y * dx
    if m1 <= 0.0:
        res[0], res[1], res[2], res[3] = rho, np.nan, np.nan, np.nan
        return
    cou = (r2x * dy - r2y * dx) / m1
    res[0], res[1], res[2], res[3] = rho, cou, rho * cou, 1.0 / m1


@njit(cache=True)
def _c(phi, k, a1, a2, f1, f2, o1, o2, out):
    """Position de la menée imposée par la paire k (dent k de la menante)."""
    return _contact(phi + k * a1, f1, f2, o1, o2, out) + k * a2


@njit(cache=True)
def _ecart(phi, a1, a2, f1, f2, o1, o2, out):
    c0 = _c(phi, 0, a1, a2, f1, f2, o1, o2, out)
    if math.isinf(c0):
        return math.inf
    c1 = _c(phi, 1, a1, a2, f1, f2, o1, o2, out)
    if math.isinf(c1):
        return -math.inf
    return c0 - c1


@njit(cache=True)
def _debut(a1, a2, f1, f2, o1, o2, out):
    """Angle de la menante au début de menée, la paire 0 rejoignant la paire 1 (FW d'Enganac) ;
    nan s'il n'y a pas de période."""
    tol = 1e-9
    pas = a1 / 8.0
    if _ecart(-2.0 * a1, a1, a2, f1, f2, o1, o2, out) <= tol:
        return math.nan
    lo, hi = -2.0 * a1, math.nan
    for k in range(1, 17):
        h = (k - 16) * pas
        if _ecart(h, a1, a2, f1, f2, o1, o2, out) <= tol:
            hi = h
            break
        lo = h
    if math.isnan(hi):
        return math.nan
    while hi - lo > 1e-7:
        m = 0.5 * (lo + hi)
        if _ecart(m, a1, a2, f1, f2, o1, o2, out) > tol:
            lo = m
        else:
            hi = m
    if math.isinf(_c(hi, 0, a1, a2, f1, f2, o1, o2, out)) or math.isinf(_c(hi, 1, a1, a2, f1, f2, o1, o2, out)):
        return math.nan
    return hi


@njit(cache=True)
def _arcs_monde(fa, fb, o, axe0, pas, k0, kn, cx, cy, ra, a0, sw):
    """Arcs de `kn` dents (indices k0-kn..k0+kn, dent k d'axe axe0 + k pas) autour de o : flanc
    menant `fa` (miroir) et flanc arrière `fb` (tel quel). Retourne le nombre d'arcs."""
    idx = 0
    for k in range(k0 - kn, k0 + kn + 1):
        th = axe0 + k * pas
        c, s = math.cos(th), math.sin(th)
        for cote in range(2):
            f = fa if cote == 0 else fb
            for m in range(f[0].shape[0]):
                lx, ly = f[0][m, 0], f[0][m, 1]
                cx[idx], cy[idx] = o[0] + c * lx - s * ly, o[1] + s * lx + c * ly
                ra[idx], a0[idx], sw[idx] = f[1][m], f[3][m] + th, f[4][m]
                idx += 1
    return idx


@njit(cache=True)
def _interference(ff1, fd1, ff2, fd2, o1, o2, phi, dpos, a1, a2, kn1, kn2):
    """Recouvrement maximal (mm) des contours complets des dents des deux roues, la menante
    tournée de phi et la menée de dpos : deux arcs de contours différents qui se coupent
    (les cercles s'interpénètrent de plus de EPS_INTERF et le point d'intersection est sur les
    deux arcs). Un simple contact (tangence) n'est pas une interférence."""
    n1, n2 = ff1[0].shape[0], ff2[0].shape[0]
    m1, m2 = (2 * kn1 + 1) * 2 * n1, (2 * kn2 + 1) * 2 * n2
    cx1, cy1, r1, a01, s1 = np.empty(m1), np.empty(m1), np.empty(m1), np.empty(m1), np.empty(m1)
    cx2, cy2, r2, a02, s2 = np.empty(m2), np.empty(m2), np.empty(m2), np.empty(m2), np.empty(m2)
    k0 = int(round(-phi / a1))
    j0 = int(round((math.pi - dpos) / a2))
    _arcs_monde(ff1, fd1, o1, phi, a1, k0, kn1, cx1, cy1, r1, a01, s1)
    _arcs_monde(ff2, fd2, o2, dpos, a2, j0, kn2, cx2, cy2, r2, a02, s2)
    pire = 0.0
    for i in range(m1):
        for j in range(m2):
            dx, dy = cx2[j] - cx1[i], cy2[j] - cy1[i]
            d = math.hypot(dx, dy)
            rs, rd = r1[i] + r2[j], abs(r1[i] - r2[j])
            if d >= rs - EPS_INTERF or d <= rd + EPS_INTERF:
                continue
            a = (r1[i] * r1[i] - r2[j] * r2[j] + d * d) / (2.0 * d)
            h2 = r1[i] * r1[i] - a * a
            if h2 <= 0.0:
                continue
            h = math.sqrt(h2)
            ux, uy = dx / d, dy / d
            mx, my = cx1[i] + a * ux, cy1[i] + a * uy
            for sg in (1.0, -1.0):
                px, py = mx - sg * h * uy, my + sg * h * ux
                if (_dans_arc(math.atan2(py - cy1[i], px - cx1[i]), a01[i], s1[i])
                        and _dans_arc(math.atan2(py - cy2[j], px - cx2[j]), a02[j], s2[j])):
                    pire = max(pire, min(rs - d, d - rd))
                    break
    return pire


@njit(cache=True)
def _analyser(f1, f2, fb1, fb2, ff1, fd1, ff2, fd2, o1, o2, a1, a2, kn1, kn2, mu, ph, ik, nj, res):
    """ENGANAC : une période d'engrènement à un entraxe. `f` : flancs menants (miroir), `fb` :
    flancs arrière (miroir, pour le jeu), `ff`/`fd` : flancs menant (miroir) et arrière (tel
    quel) des contours complets, pour les interférences. Remplit `res` (colonnes VARC…).
    Retourne 0 si tout va bien, 3 sans période, 4 en cas d'arc-boutement."""
    out, ef = np.empty(6), np.empty(4)
    d0 = _debut(a1, a2, f1, f2, o1, o2, out)
    if math.isnan(d0):
        return 3
    em = ph[0] * ph[3] / (ph[0] + ph[3])
    pa = 0.5 * (ph[1] + ph[4])
    dm = min(ph[2], ph[5])
    ca, ci, cm, nm, rend, ha = -1e33, 1e33, 0.0, 1e33, 0.0, -1e33
    x_debut, x_fin, interf = 0.0, 0.0, 0.0
    pas_i = max(1, (ik - 1) // 5)
    sortie = False
    for i in range(ik):
        phi = d0 + i * a1 / (ik - 1)
        if i == 0:  # comme ENGANAC : 1e-6 d'écart pour départager les deux paires aux extrémités
            phi += 1e-6
        if i == ik - 1:
            phi -= 1e-6
        c0 = _contact(phi, f1, f2, o1, o2, out)
        if math.isinf(c0):
            return 3
        q0 = out.copy()
        q1, c1 = q0, math.inf
        if not sortie:
            d1 = _contact(phi + a1, f1, f2, o1, o2, out)
            q1, c1 = out.copy(), d1 + a2
            sortie = (not math.isfinite(c1)) or c1 > c0 + 1e-3
        pair1 = (not sortie) and c1 < c0
        q = q1 if pair1 else q0
        _efforts(o1, o2, q[0], q[1], q[2], q[3], mu, ef)
        if not math.isfinite(ef[1]):
            return 4
        ce = abs(ef[1])
        ca, ci, cm = max(ca, ce), min(ci, ce), cm + ce
        rend += ef[2]
        nm = min(nm, ef[2])
        s = q[4] + q[5]
        if abs(s) < 1e-12:
            return 4
        ha = max(ha, math.sqrt(ef[3] * em / (math.pi * dm * (1.0 - pa * pa) * abs(q[4] * q[5] / s))))
        if i == 0:
            x_debut = q[0]
        if i == ik - 1:
            x_fin = q[0]
        if i % pas_i == 0 or i == ik - 1:
            interf = max(interf, _interference(ff1, fd1, ff2, fd2, o1, o2, phi, c1 if pair1 else c0,
                                               a1, a2, kn1, kn2))
    res[VARC] = (ca - ci) / (0.5 * (ca + ci))
    res[REND] = rend / ik
    res[HERTZ] = ha
    res[REND_MIN] = nm
    res[COUPLE_MOY] = cm / ik
    res[INTERF] = interf
    qq = -x_debut / x_fin if x_fin != 0.0 else 0.0  # symétrie de l'engrènement (SymL)
    if qq > 0.0:
        res[SYML] = min(qq, 1.0 / qq)
    else:
        res[SYML] = max(qq, 1.0 / qq) if qq != 0.0 else -1e33
    res[JEU_MAX] = res[JEU_MIN] = math.nan
    if nj > 0:  # jeu de la menée : theta_max(phi) + theta_max_arrière(-phi), replié modulo un pas
        jmax, jmin = -1e33, 1e33
        for j in range(nj):
            phi = d0 + j * a1 / (nj - 1)
            tm = 1e33
            for k in (-1, 0, 1):
                tm = min(tm, _c(phi, k, a1, a2, f1, f2, o1, o2, out))
            tb = 1e33
            for k in (-1, 0, 1):
                tb = min(tb, _c(-phi, k, a1, a2, fb1, fb2, o1, o2, out))
            jj = math.degrees((tm + tb + a2 / 2.0) % a2 - a2 / 2.0)
            jmax, jmin = max(jmax, jj), min(jmin, jj)
        res[JEU_MAX], res[JEU_MIN] = jmax, jmin
    return 0


@njit(cache=True)
def evaluer(p1, p2, p1d, p2d, z1, z2, mo1, mo2, sds, mu, ph, ik, njs, rk, sortie):
    """Un couple de dents (demi-profils menants p1, p2 et demi-profils arrière p1d, p2d ; points
    de transition, refermés par KORZF sur une copie) à chacun des entraxes `sds` (mm) ; `njs` :
    positions de jeu par entraxe (0 = pas de jeu). `sortie` (len(sds), NCOL). Code de retour : 0
    ok, 1 entraxe trop petit, 2 trop grand, 3 pas de période, 4 arc-boutement, 5 profil invalide
    (SWTE avec le rayon de tête minimal `rk`)."""
    xy1, xy2, xd1, xd2 = p1.copy(), p2.copy(), p1d.copy(), p2d.copy()
    n1, n2 = xy1.shape[0], xy2.shape[0]
    pm1, ra1, pm2, ra2 = np.zeros((n1 - 1, 2)), np.zeros(n1 - 1), np.zeros((n2 - 1, 2)), np.zeros(n2 - 1)
    pd1, rd1, pd2, rd2 = np.zeros((n1 - 1, 2)), np.zeros(n1 - 1), np.zeros((n2 - 1, 2)), np.zeros(n2 - 1)
    korzf(-math.pi / z1, xy1, pm1, ra1)
    korzf(-math.pi / z2, xy2, pm2, ra2)
    korzf(-math.pi / z1, xd1, pd1, rd1)
    korzf(-math.pi / z2, xd2, pd2, rd2)
    if (_swte(xy1, ra1, rk) != 1 or _swte(xy2, ra2, rk) != 1
            or _swte(xd1, rd1, rk) != 1 or _swte(xd2, rd2, rk) != 1):
        return 5
    f1, f2 = _flanc(xy1, pm1, ra1, -1.0), _flanc(xy2, pm2, ra2, -1.0)
    fb1, fb2 = _flanc(xd1, pd1, rd1, -1.0), _flanc(xd2, pd2, rd2, -1.0)
    fd1, fd2 = _flanc(xd1, pd1, rd1, 1.0), _flanc(xd2, pd2, rd2, 1.0)
    a1, a2 = TWO_PI / z1, TWO_PI / z2
    pied1, pied2 = math.hypot(xy1[0, 0], xy1[0, 1]), math.hypot(xy2[0, 0], xy2[0, 1])
    t1, t2 = xy1[n1 - 1, 1], xy2[n2 - 1, 1]
    r_tete = min(abs(ra1[n1 - 2]), abs(ra2[n2 - 2]), abs(rd1[n1 - 2]), abs(rd2[n2 - 2]))
    ecart_tete = max(abs(xy1[n1 - 1, 1] - xd1[n1 - 1, 1]), abs(xy2[n2 - 1, 1] - xd2[n2 - 1, 1]))
    r1n, r2n = mo1 * z1 / 2.0, mo2 * z2 / 2.0
    for k in range(sds.shape[0]):
        sd = sds[k]
        if pied1 + t2 > sd or pied2 + t1 > sd:
            return 1
        h1, h2 = (sd * sd + t1 * t1 - t2 * t2) / (2.0 * sd * t1), (sd * sd + t2 * t2 - t1 * t1) / (2.0 * sd * t2)
        if abs(h1) > 1.0:
            return 2
        # nombre de dents de chaque roue de part et d'autre de la ligne des centres à contrôler
        kn1 = int(math.acos(min(1.0, abs(h1))) / a1) + 2
        kn2 = int(math.acos(min(1.0, abs(h2))) / a2) + 2
        cy = sd / (r1n + r2n)  # cercles primitifs de fonctionnement, tangents en (0,0)
        o1, o2 = np.array([0.0, -r1n * cy]), np.array([0.0, r2n * cy])
        fi = _analyser(f1, f2, fb1, fb2, f1, fd1, f2, fd2, o1, o2, a1, a2, kn1, kn2, mu, ph, ik, njs[k],
                       sortie[k])
        if fi != 0:
            return fi
        sortie[k, R_TETE], sortie[k, ECART_TETE] = r_tete, ecart_tete
    return 0


@njit(cache=True)
def _construire(p0, th0, k, lg, n, out):
    """Chaîne d'arcs de courbures `k` et de longueurs `lg` depuis le point p0, tangente th0 ;
    les n-2 premiers arcs seulement (le dernier est refermé par KORZF)."""
    out[0, 0], out[0, 1] = p0[0], p0[1]
    th = th0
    for j in range(k.shape[0]):
        dth = k[j] * lg[j]
        m = 1.0 if abs(dth) < 1e-12 else math.sin(0.5 * dth) / (0.5 * dth)
        out[j + 1, 0] = out[j, 0] + lg[j] * m * math.cos(th + 0.5 * dth)
        out[j + 1, 1] = out[j, 1] + lg[j] * m * math.sin(th + 0.5 * dth)
        th += dth
    out[n - 1, 0], out[n - 1, 1] = 0.0, 0.0  # recalculé par KORZF


@njit(cache=True, parallel=True)
def evaluer_lot(k1, l1, k2, l2, k1d, l1d, k2d, l2d, p0, th0, z1, z2, mo1, mo2, sds, mu, ph, ik, njs, rk,
                sortie, codes):
    """Lot de N couples de dents donnés par les courbures et longueurs d'arcs des demi-profils
    menants (k1, l1, k2, l2) et arrière (k1d, …), tableaux (N, n-2) ; `sortie`
    (N, len(sds), NCOL), `codes` (N,)."""
    n1, n2 = k1.shape[1] + 2, k2.shape[1] + 2
    for i in prange(k1.shape[0]):
        a, b, c, d = np.empty((n1, 2)), np.empty((n2, 2)), np.empty((n1, 2)), np.empty((n2, 2))
        _construire(p0[0], th0[0], k1[i], l1[i], n1, a)
        _construire(p0[1], th0[1], k2[i], l2[i], n2, b)
        _construire(p0[0], th0[0], k1d[i], l1d[i], n1, c)
        _construire(p0[1], th0[1], k2d[i], l2d[i], n2, d)
        codes[i] = evaluer(a, b, c, d, z1, z2, mo1, mo2, sds, mu, ph, ik, njs, rk, sortie[i])


# ---------------------------------------------------------------------------------------
# problème d'optimisation
# ---------------------------------------------------------------------------------------

CRITERES = ("couple", "jeu", "rendement", "hertz")  # clés de l'objectif principal / secondaire


@dataclass
class Objectif:
    """Fonction objectif et contraintes.

    J = f1/f1₀ + poids·f2/f2₀ où f1, f2 sont les critères principal et secondaire (`couple` :
    variation du couple, `jeu` : jeu maximal de la menée, `rendement` : perte de rendement,
    `hertz` : pression de Hertz maximale ; valeur quadratique moyenne sur les entraxes -dA, 0,
    +dA, comme Resulto) et f₀ leur valeur pour le profil de départ. J vaut donc 1 au départ.
    Les contraintes ajoutent des pénalités à J : jeu de la menée au moins `jeu_min` (degrés) à
    l'entraxe le plus petit (sinon blocage) et au plus `jeu_max` à l'entraxe le plus grand (None :
    pas de limite), rayon de l'arc de tête au moins `r_tete_min` (mm, durée de vie des
    fraises), pas d'interférence entre les dents, tête arrière égale à la tête avant.
    """

    principal: str = "couple"
    secondaire: str = None
    poids: float = 0.5
    jeu_min: float = 0.1
    jeu_max: float = None
    r_tete_min: float = 0.0
    p_jeu: float = 10.0  # par degré de dépassement
    p_tete: float = 10.0  # par fraction (1 = 100 %) de rayon en dessous du minimum
    p_interf: float = 1000.0  # par mm de recouvrement (1 µm pèse autant que J au départ)
    p_ecart: float = 1000.0  # par mm de différence de hauteur de tête entre les deux flancs
    rk: float = 0.0  # rayon de tête minimal / rayon de pied (SWTE de Resulto)


class Probleme:
    """Deux roues à optimiser : décodage des variables, évaluation par lots, objectif.

    `roues` : les deux `evolvente.Roue` (module, nombre de dents), `mobiles` : les deux
    `arcs.Mobile` de départ (demi-profil menant, éventuellement demi-profil arrière, et données de
    Hertz). `libres[k]` : pour la roue k (0 = menante), les indices des arcs variables parmi les
    n-2 premiers (tous par défaut) ; le dernier arc est toujours refermé. Un vecteur de variables x
    contient, roue par roue, la variation de courbure puis de longueur de chaque arc libre (unités
    `ECHELLE` fois la courbure médiane, fois la longueur de l'arc), pour le demi-profil menant puis, si
    `dos_libre`, pour le demi-profil arrière (variation par rapport au profil menant symétrique,
    ou au profil arrière de départ s'il est déjà dissymétrique).
    """

    ECHELLE = 0.03

    def __init__(self, roues, mobiles, mu=0.2, da_neg=-30.0, da_pos=30.0, libres=None, objectif=None,
                 dos_libre=False, avant_libre=True, calibrer=True):
        """`calibrer` : calcule tout de suite la référence `f0` sur le profil de départ (x=0),
        ValueError s'il est invalide (SWTE, entraxe, période) ; `calibrer=False` diffère ce calcul
        (`self.f0 = None`) à un appel explicite de `_criteres_ref`, pour un profil de départ dont
        on ne sait pas encore s'il est utilisable (voir `recherche.rechercher`, qui en cherche un
        valide à proximité plutôt que d'abandonner)."""
        self.roues, self.mobiles0, self.mu = roues, mobiles, mu
        self.da_um = (da_neg, 0.0, da_pos)
        self.objectif = objectif or Objectif()
        if any(m.e_mod is None or m.poisson is None or m.epaisseur is None for m in mobiles):
            raise ValueError("données de Hertz manquantes")
        self.ph = np.array([mobiles[0].e_mod, mobiles[0].poisson, mobiles[0].epaisseur,
                            mobiles[1].e_mod, mobiles[1].poisson, mobiles[1].epaisseur], dtype=float)
        self.z = tuple(r.zz for r in roues)
        self.mo = tuple(r.mo for r in roues)
        self.sd0 = 0.5 * (self.z[0] * self.mo[0] + self.z[1] * self.mo[1])
        self.sds = self.sd0 + np.array(self.da_um) * 1e-3
        self.ref = [self._fermer(m.points, z) for m, z in zip(mobiles, self.z)]
        self.p0 = np.array([r[0] for r in self.ref])
        self.n = [len(r) for r in self.ref]
        self.th0, self.k0, self.l0, self.k_ech = [], [], [], []
        for r, z in zip(self.ref, self.z):
            th0, k, lg = self._vers_kl(r, z)
            self.th0.append(th0)
            self.k0.append(k)
            self.l0.append(lg)
            self.k_ech.append(self.ECHELLE * float(np.median(np.abs(k))))
        self.th0 = np.array(self.th0)
        # profil arrière de départ : celui de chaque mobile s'il est dissymétrique, sinon le profil
        # menant (les variables arrière sont alors des écarts par rapport au profil menant courant)
        self.dos_relatif = all(m.points_dos is None for m in mobiles)
        self.k0d, self.l0d = [], []
        for m, z, th0, k, lg in zip(mobiles, self.z, self.th0, self.k0, self.l0):
            if m.points_dos is None:
                self.k0d.append(k)
                self.l0d.append(lg)
            else:
                self.k0d.append(self._vers_kl(self._fermer(m.points_dos, z), z, th0)[1])
                self.l0d.append(self._vers_kl(self._fermer(m.points_dos, z), z, th0)[2])
        n_arcs = [n - 2 for n in self.n]
        if libres is None:
            libres = [range(n) for n in n_arcs]
        self.libres = [np.array(sorted(set(li)), dtype=int) for li in libres]
        for li, n in zip(self.libres, n_arcs):
            if len(li) and (li[0] < 0 or li[-1] >= n):
                raise ValueError(f"arc libre hors de 0..{n - 1}")
        self.dos_libre, self.avant_libre = dos_libre, avant_libre
        self.ik, self.n_jeu = IK, N_JEU
        self.f0 = self._criteres_ref() if calibrer else None

    @property
    def n_var(self):
        return 2 * sum(len(li) for li in self.libres) * (int(self.avant_libre) + int(self.dos_libre))

    # --- profils -----------------------------------------------------------------------

    @staticmethod
    def _fermer(points, zz):
        """Points de transition refermés par KORZF (dernier arc centré sur l'axe de la dent)."""
        xy = np.array(points, dtype=float)
        korzf(-math.pi / zz, xy, np.zeros((len(xy) - 1, 2)), np.zeros(len(xy) - 1))
        return xy

    @staticmethod
    def _vers_kl(xy, zz, th0=None):
        """(tangente de départ dans le sens de parcours, courbures, longueurs) des n-2 premiers arcs."""
        if th0 is None:
            a0 = -math.pi / zz  # pente de KRBGS ; le sens de parcours est celui qui va vers la corde
            b = math.atan2(xy[1, 1] - xy[0, 1], xy[1, 0] - xy[0, 0])
            th0 = a0 if math.cos(b - a0) > 0 else a0 + math.pi
        th = th0
        k, lg = [], []
        for j in range(len(xy) - 2):
            dx, dy = xy[j + 1] - xy[j]
            corde = math.hypot(dx, dy)
            dth = (2.0 * (math.atan2(dy, dx) - th) + math.pi) % TWO_PI - math.pi
            long = corde if abs(dth) < 1e-12 else corde * (dth / 2.0) / math.sin(dth / 2.0)
            k.append(dth / long)
            lg.append(long)
            th += dth
        return th0, np.array(k), np.array(lg)

    def decoder(self, x):
        """Variables x -> courbures et longueurs (k1, l1, k2, l2, k1d, l1d, k2d, l2d), tableaux
        (N, n-2), pour N vecteurs de variables (x de forme (n_var,) ou (N, n_var))."""
        x = np.atleast_2d(np.asarray(x, dtype=float))
        avant, arriere, i = [], [], 0
        for w in range(2):
            k = np.tile(self.k0[w], (len(x), 1))
            lg = np.tile(self.l0[w], (len(x), 1))
            m = len(self.libres[w])
            if self.avant_libre:
                k[:, self.libres[w]] += x[:, i:i + m] * self.k_ech[w]
                lg[:, self.libres[w]] *= 1.0 + self.ECHELLE * x[:, i + m:i + 2 * m]
                i += 2 * m
            avant.append((k, lg))
        for w in range(2):
            if self.dos_relatif:
                k, lg = avant[w][0].copy(), avant[w][1].copy()
            else:
                k, lg = np.tile(self.k0d[w], (len(x), 1)), np.tile(self.l0d[w], (len(x), 1))
            if self.dos_libre:
                m = len(self.libres[w])
                k[:, self.libres[w]] += x[:, i:i + m] * self.k_ech[w]
                lg[:, self.libres[w]] *= 1.0 + self.ECHELLE * x[:, i + m:i + 2 * m]
                i += 2 * m
            arriere.append((k, lg))
        return [avant[0][0], avant[0][1], avant[1][0], avant[1][1],
                arriere[0][0], arriere[0][1], arriere[1][0], arriere[1][1]]

    def points(self, x):
        """Points de transition fermés (menant 1, menant 2, arrière 1, arrière 2) pour un vecteur."""
        d = self.decoder(x)
        res = []
        for (k, lg), w in zip(((d[0][0], d[1][0]), (d[2][0], d[3][0]), (d[4][0], d[5][0]), (d[6][0], d[7][0])),
                              (0, 1, 0, 1)):
            xy = np.empty((self.n[w], 2))
            _construire(self.p0[w], self.th0[w], k, lg, self.n[w], xy)
            res.append(self._fermer(xy, self.z[w]))
        return res

    def ecart_um(self, xa, xb=None):
        """Écart maximal (µm) entre les demi-profils des deux roues (menants, ou arrière si seules
        les variables arrière sont libres) pour deux vecteurs de variables (xb : profil de départ
        si None) : plus grande distance d'un point d'un profil à l'autre, échantillonné en arcs."""
        pa, pb_ = self.points(xa), self.points(np.zeros(self.n_var) if xb is None else xb)
        pire, decal = 0.0, 0 if self.avant_libre else 2  # variables arrière seules : on compare le dos
        for w in range(2):
            a = self._echantillon(pa[w + decal], w)
            b = self._echantillon(pb_[w + decal], w)
            d = np.hypot(a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1])
            pire = max(pire, d.min(axis=1).max(), d.min(axis=0).max())
        return 1e3 * pire

    def _echantillon(self, pts, w):
        pm, ra = np.zeros((len(pts) - 1, 2)), np.zeros(len(pts) - 1)
        krbgs(-math.pi / self.z[w], np.array(pts, dtype=float), pm, ra)
        return arcs.echantillonner(pts, pm, ra, 160)

    def echelle_um(self, n=24, graine=0):
        """Écart de profil (µm) produit en moyenne par des variables d'écart-type 1."""
        rng = np.random.default_rng(graine)
        val = []
        for _ in range(n):
            try:
                val.append(self.ecart_um(rng.normal(size=self.n_var)))
            except (ValueError, FloatingPointError):
                continue
        return float(np.median(val)) if val else 1.0

    def mobiles(self, x):
        """Les deux `arcs.Mobile` (points refermés, demi-profil arrière si différent) pour x."""
        p1, p2, d1, d2 = self.points(x)
        res = []
        for m, p, d in zip(self.mobiles0, (p1, p2), (d1, d2)):
            dos = None if np.abs(p - d).max() < 1e-9 else d
            res.append(arcs.Mobile(m.nom, m.mo, m.zz, p, list(m.drapeaux), m.e_mod, m.poisson, m.epaisseur,
                                   dos))
        return res

    # --- évaluation --------------------------------------------------------------------

    def evaluer(self, x, avec_jeu=True):
        """Évalue N vecteurs de variables : (S, codes), S (N, 3, NCOL) : colonnes VARC… par entraxe
        (-dA, 0, +dA). Le jeu n'est évalué qu'aux entraxes extrêmes (et pas si `avec_jeu` est faux)."""
        d = self.decoder(x)
        njs = np.array([self.n_jeu if avec_jeu else 0, 0, self.n_jeu if avec_jeu else 0], dtype=np.int64)
        s = np.zeros((len(d[0]), 3, NCOL))
        codes = np.zeros(len(d[0]), dtype=np.int64)
        evaluer_lot(*d, self.p0, self.th0, self.z[0], self.z[1], self.mo[0], self.mo[1],
                    self.sds, self.mu, self.ph, self.ik, njs, self.objectif.rk, s, codes)
        return s, codes

    @staticmethod
    def critere(s, cle):
        """Critère `cle` de chaque profil (valeur quadratique moyenne sur les 3 entraxes)."""
        if cle == "couple":
            v = s[:, :, VARC]
        elif cle == "jeu":  # jeu maximal aux entraxes extrêmes (il n'est pas évalué à l'entraxe nominal)
            return np.sqrt(np.mean(np.nan_to_num(s[:, ::2, JEU_MAX], nan=0.0) ** 2, axis=1))
        elif cle == "rendement":
            v = 1.0 - s[:, :, REND]
        elif cle == "hertz":
            v = s[:, :, HERTZ]
        else:
            raise ValueError(f"critère inconnu : {cle}")
        return np.sqrt(np.mean(np.nan_to_num(v, nan=0.0) ** 2, axis=1))

    def _criteres_ref(self, x=None):
        """Valeur des critères (référence f0) au vecteur `x` (le profil de départ si None) ;
        ValueError si ce profil est invalide."""
        s, codes = self.evaluer(np.zeros(self.n_var) if x is None else x)
        if codes[0]:
            raise ValueError(f"profil invalide : {MESSAGES[int(codes[0])]}")
        return {c: max(float(self.critere(s, c)[0]), 1e-12) for c in CRITERES}

    def violations(self, s):
        """Dépassements des contraintes de chaque profil : dictionnaire de tableaux (N,) — jeu
        minimal (°), jeu maximal (°), rayon de tête (fraction), interférence (mm), écart de tête (mm)."""
        ob = self.objectif
        n = len(s)
        v = {"jeu_min": np.zeros(n), "jeu_max": np.zeros(n), "r_tete": np.zeros(n),
             "interference": np.maximum(0.0, s[:, :, INTERF].max(axis=1)),
             "ecart_tete": np.maximum(0.0, s[:, 0, ECART_TETE] - 1e-3)}  # 1 µm de différence de tête toléré
        if ob.jeu_min is not None:
            v["jeu_min"] = np.maximum(0.0, ob.jeu_min - s[:, 0, JEU_MIN])
        if ob.jeu_max is not None:
            v["jeu_max"] = np.maximum(0.0, s[:, 2, JEU_MAX] - ob.jeu_max)
        if ob.r_tete_min > 0:
            v["r_tete"] = np.maximum(0.0, (ob.r_tete_min - s[:, 0, R_TETE]) / ob.r_tete_min)
        return v

    def valeur(self, s, codes):
        """Fonction objectif J de chaque profil, pénalités comprises (`PENALITE` si invalide)."""
        ob = self.objectif
        j = self.critere(s, ob.principal) / self.f0[ob.principal]
        if ob.secondaire:
            j = j + ob.poids * self.critere(s, ob.secondaire) / self.f0[ob.secondaire]
        v = self.violations(s)
        j = (j + ob.p_jeu * (np.nan_to_num(v["jeu_min"]) + np.nan_to_num(v["jeu_max"])) + ob.p_tete * v["r_tete"]
             + ob.p_interf * v["interference"] + ob.p_ecart * v["ecart_tete"])
        return np.where(codes == 0, j, PENALITE)

    def __call__(self, x):
        """Objectif de N vecteurs de variables (tableau (N,))."""
        s, codes = self.evaluer(x)
        return self.valeur(s, codes)

    def rapport(self, x=None, n_jeu=60):
        """Résultats détaillés (un dictionnaire par entraxe) d'un vecteur de variables, ou du
        profil de départ si x est None ; le jeu est évalué en `n_jeu` positions."""
        x = np.zeros(self.n_var) if x is None else x
        d = self.decoder(x)
        s = np.zeros((1, 3, NCOL))
        code = np.zeros(1, dtype=np.int64)
        evaluer_lot(*d, self.p0, self.th0, self.z[0], self.z[1], self.mo[0], self.mo[1],
                    self.sds, self.mu, self.ph, self.ik, np.full(3, n_jeu, dtype=np.int64),
                    self.objectif.rk, s, code)
        if code[0]:
            raise ValueError(MESSAGES[int(code[0])])
        return [{"da_um": da, "couple_var": 100 * r[VARC], "rendement": 100 * r[REND],
                 "hertz": r[HERTZ], "jeu_min": r[JEU_MIN], "jeu_max": r[JEU_MAX], "symetrie": r[SYML],
                 "interference_um": 1e3 * r[INTERF], "r_tete": r[R_TETE], "ecart_tete_um": 1e3 * r[ECART_TETE]}
                for da, r in zip(self.da_um, s[0])]


def profil_valide(mobile, rk=0.0):
    """True si le demi-profil menant de `mobile` passe le contrôle SWTE (une seule inversion du
    signe des rayons de courbure, tête pas trop pointue par rapport au pied) : condition
    nécessaire pour servir de profil de départ à l'optimisation (voir `Probleme._criteres_ref`)."""
    xy = Probleme._fermer(mobile.points, mobile.zz)
    n = len(xy)
    pm, ra = np.zeros((n - 1, 2)), np.zeros(n - 1)
    korzf(-math.pi / mobile.zz, xy, pm, ra)
    return _swte(xy, ra, rk) == 1


def main():
    """Contrôle : compare le noyau à `contact.Engrenement` sur test.Cev et mesure sa vitesse."""
    import sys
    import time

    import contact
    import evolvente

    fichier = sys.argv[1] if len(sys.argv) > 1 else "test.Cev"
    narcs = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    projet = evolvente.lire_projet(fichier)
    profils = [evolvente.evolve(r) for r in projet.roues]
    hertz = ((100000.0, 0.3, 0.3), (200000.0, 0.3, 0.28))
    mobiles = []
    for i, (r, p) in enumerate(zip(projet.roues, profils)):
        m = arcs.vers_mobile(r, arcs.approximer(p, r.zz, r.qg, narcs), "m" if i == 0 else "e")
        m.e_mod, m.poisson, m.epaisseur = hertz[i]
        mobiles.append(m)
    pb = Probleme(projet.roues, mobiles, objectif=Objectif(jeu_max=None))
    print(f"{fichier} : {pb.z[0]}/{pb.z[1]} dents, {narcs} arcs, {pb.n_var} variables")
    ecart = 0.0
    for r in pb.rapport(None, n_jeu=101):
        sim = contact.Engrenement(projet.roues, profils, r["da_um"] * 1e-3, mobiles=pb.mobiles(np.zeros(pb.n_var))
                                  ).simuler(pb.mu, IK, n_jeu=101)
        st = sim.statistiques()
        ref = (st["couple"][3], sim.hertz_max()[0], sim.jeu_min, sim.jeu_max)
        neuf = (r["couple_var"], r["hertz"], r["jeu_min"], r["jeu_max"])
        ecart = max(ecart, abs(ref[0] - neuf[0]) / max(ref[0], 1e-9), abs(ref[1] - neuf[1]) / ref[1])
        print(f"  dA={r['da_um']:+5.0f} µm  contact.py : varC {ref[0]:.3f} %  Hertz {ref[1]:.1f}  jeu {ref[2]:.3f}..{ref[3]:.3f}"
              f"   optim.py : varC {neuf[0]:.3f} %  Hertz {neuf[1]:.1f}  jeu {neuf[2]:.3f}..{neuf[3]:.3f}"
              f"  interférence {r['interference_um']:.1f} µm")
    print(f"  écart relatif max (varC, Hertz) : {ecart:.1e}")
    x = np.random.default_rng(0).normal(0, 2.0, (2000, pb.n_var))
    pb(x[:16])  # compilation
    t = time.perf_counter()
    f = pb(x)
    dt = time.perf_counter() - t
    print(f"  2000 profils voisins : {np.mean(f < PENALITE):.0%} valides, {dt:.2f} s ({2000 / dt:.0f} évaluations/s)")


if __name__ == "__main__":
    main()
