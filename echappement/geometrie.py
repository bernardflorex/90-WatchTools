"""Géométrie 2D des pièces de l'échappement et recherche des contacts.

Chaque pièce est un corps qui tourne autour d'un axe fixe. Son contour est
une suite de segments (les arcs du DXF sont échantillonnés finement), donnés
dans la POSE DESSINÉE ; la pièce tournée d'un angle a s'obtient par rotation
autour de son axe. Toutes les longueurs sont en mm, les angles en radians,
comptés dans le sens trigonométrique de la vue standard (balancier en haut,
roue en bas, palette d'entrée à gauche).

Le contact entre deux pièces est la distance minimale entre leurs contours :
nulle quand elles se touchent ou s'interpénètrent. C'est une approche
purement géométrique, sans raideur de contact ; elle sert au relevé du
tracé et à la cinématique quasi statique.
"""

from dataclasses import dataclass, field

import numpy as np

#: flèche maximale tolérée en échantillonnant un arc (mm) : 0,05 µm
FLECHE_MAX = 5e-5
#: pas angulaire maximal en échantillonnant un arc (rad). Le contour en
#: facettes ne sert plus qu'au dessin et au test de liberté : les contacts
#: sont cherchés sur les vrais arcs (premier_contact_exact)
PAS_ARC_MAX = np.radians(3.0)
#: deux contours plus proches que cette distance sont en contact (mm) : 1 nm
TOL_CONTACT = 1e-6


def rotation(points, centre, angle):
    """Tourne des points (..., 2) d'un angle autour de `centre`."""
    c, s = np.cos(angle), np.sin(angle)
    p = np.asarray(points, dtype=float) - centre
    return np.stack([c * p[..., 0] - s * p[..., 1],
                     s * p[..., 0] + c * p[..., 1]], axis=-1) + centre


def echantillonner_arc(centre, rayon, a0, a1):
    """Points d'un arc parcouru de a0 à a1 (rad, sens trigonométrique)."""
    while a1 <= a0:
        a1 += 2 * np.pi
    pas = 2.0 * np.sqrt(2.0 * FLECHE_MAX / rayon) if rayon > 2 * FLECHE_MAX else a1 - a0
    # pas angulaire plafonné : voir PAS_ARC_MAX
    pas = min(pas, PAS_ARC_MAX)
    n = max(2, int(np.ceil((a1 - a0) / pas)) + 1)
    t = np.linspace(a0, a1, n)
    return np.asarray(centre) + rayon * np.c_[np.cos(t), np.sin(t)]


@dataclass(eq=False)
class Element:
    """Élément de contour tel que dessiné : segment ou arc, avec une étiquette
    (« face de repos d'entrée »...). Les points sont échantillonnés."""
    points: np.ndarray            # (n, 2), dans la pose dessinée
    etiquette: str = ""
    centre: np.ndarray = None     # arcs seulement
    rayon: float = 0.0

    @property
    def est_arc(self):
        return self.centre is not None

    @property
    def longueur(self):
        return float(np.sum(np.hypot(*np.diff(self.points, axis=0).T)))


@dataclass(eq=False)
class Corps:
    """Pièce tournant autour d'un axe fixe (ou pièce fixe : butée)."""
    nom: str
    axe: np.ndarray
    elements: list = field(default_factory=list)
    fixe: bool = False

    def __post_init__(self):
        self.axe = np.asarray(self.axe, dtype=float)
        self._reconstruire()

    def _reconstruire(self):
        segs, ids = [], []
        for k, e in enumerate(self.elements):
            p = e.points
            segs.append(np.stack([p[:-1], p[1:]], axis=1))
            ids.append(np.full(len(p) - 1, k))
        self.segments = np.concatenate(segs) if segs else np.zeros((0, 2, 2))
        self.indices = np.concatenate(ids) if ids else np.zeros(0, int)

    def segments_a(self, angle):
        """Segments (n, 2, 2) de la pièce tournée de `angle`."""
        if angle == 0.0 or self.fixe:
            return self.segments
        return rotation(self.segments, self.axe, angle)

    def vers_local(self, point, angle):
        """Point du plan ramené dans la pose dessinée de la pièce."""
        return rotation(point, self.axe, -angle)

    def element(self, etiquette):
        return next(e for e in self.elements if e.etiquette == etiquette)


def _cr(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def _point_segment(p, a, b):
    """Distance et point le plus proche de p sur [a, b], vectorisé."""
    ab = b - a
    l2 = np.maximum(np.sum(ab * ab, axis=-1), 1e-300)
    t = np.clip(np.sum((p - a) * ab, axis=-1) / l2, 0.0, 1.0)
    q = a + t[..., None] * ab
    return np.linalg.norm(p - q, axis=-1), q


def distance(S1, S2, marge=0.05):
    """Distance minimale entre deux ensembles de segments.

    Retourne (d, i, j, P1, P2) : la distance, les indices des segments les plus
    proches et les points correspondants. d vaut 0 si deux segments se coupent.
    Seuls les couples de segments dont les boîtes, élargies de `marge`, se
    recouvrent sont examinés ; au-delà, retourne d = inf.
    """
    lo1, hi1 = S1.min(axis=1), S1.max(axis=1)
    lo2, hi2 = S2.min(axis=1), S2.max(axis=1)
    # présélection grossière sur la boîte de chaque ensemble
    g1 = np.all((hi1 >= lo2.min(0) - marge) & (lo1 <= hi2.max(0) + marge), axis=1)
    g2 = np.all((hi2 >= lo1.min(0) - marge) & (lo2 <= hi1.max(0) + marge), axis=1)
    i1, i2 = np.nonzero(g1)[0], np.nonzero(g2)[0]
    if len(i1) == 0 or len(i2) == 0:
        return np.inf, -1, -1, None, None
    A, B = S1[i1], S2[i2]
    ok = np.all((A.max(1)[:, None] >= B.min(1)[None] - marge)
                & (A.min(1)[:, None] <= B.max(1)[None] + marge), axis=2)
    ia, ib = np.nonzero(ok)
    if len(ia) == 0:
        return np.inf, -1, -1, None, None
    a0, a1, b0, b1 = A[ia, 0], A[ia, 1], B[ib, 0], B[ib, 1]
    r, s = a1 - a0, b1 - b0
    den = _cr(r, s)
    qp = b0 - a0
    with np.errstate(divide="ignore", invalid="ignore"):
        t, u = _cr(qp, s) / den, _cr(qp, r) / den
    coupe = (den != 0) & (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1)
    if coupe.any():
        k = int(np.nonzero(coupe)[0][0])
        P = a0[k] + t[k] * r[k]
        return 0.0, int(i1[ia[k]]), int(i2[ib[k]]), P, P
    cand = []
    d, q = _point_segment(a0, b0, b1); cand.append((d, a0, q))
    d, q = _point_segment(a1, b0, b1); cand.append((d, a1, q))
    d, q = _point_segment(b0, a0, a1); cand.append((d, q, b0))
    d, q = _point_segment(b1, a0, a1); cand.append((d, q, b1))
    D = np.stack([c[0] for c in cand])            # (4, n)
    m, k = np.unravel_index(np.argmin(D), D.shape)
    return float(D[m, k]), int(i1[ia[k]]), int(i2[ib[k]]), cand[m][1][k], cand[m][2][k]


def en_contact(c1, a1, c2, a2):
    return distance(c1.segments_a(a1), c2.segments_a(a2))[0] < TOL_CONTACT


def tourner_jusquau_contact(mobile, angle0, sens, obstacles, pas=np.radians(0.02),
                            course=np.radians(30.0), tol=1e-10):
    """Tourne `mobile` à partir de `angle0` dans le sens `sens` (+1/-1) jusqu'au
    premier contact avec l'un des `obstacles` [(corps, angle)].

    Si la position de départ est déjà en contact, recule d'abord jusqu'à la
    dégager : le résultat est alors la position de contact la plus avancée
    dans le sens demandé (cas du recul de la roue au dégagement).
    Retourne l'angle de contact, ou None si rien n'est touché sur `course`.
    """
    def touche(a):
        S = mobile.segments_a(a)
        return any(distance(S, c.segments_a(ac))[0] < TOL_CONTACT for c, ac in obstacles)

    a = angle0
    if touche(a):
        libre = None
        for k in range(1, int(course / pas) + 1):
            if not touche(angle0 - sens * k * pas):
                libre = angle0 - sens * k * pas
                break
        if libre is None:
            return None
        bas, haut = libre, libre + sens * pas
    else:
        haut = None
        for k in range(1, int(course / pas) + 1):
            if touche(angle0 + sens * k * pas):
                haut = angle0 + sens * k * pas
                break
        if haut is None:
            return None
        bas = haut - sens * pas
    while abs(haut - bas) > tol:
        m = 0.5 * (bas + haut)
        if touche(m):
            haut = m
        else:
            bas = m
    return bas


def point_de_contact(c1, a1, c2, a2):
    """Points les plus proches entre deux corps posés, et l'élément de chacun
    qui porte le contact : (P, element1, element2)."""
    d, i, j, P1, P2 = distance(c1.segments_a(a1), c2.segments_a(a2))
    if i < 0:
        return None, None, None
    return 0.5 * (P1 + P2), c1.elements[c1.indices[i]], c2.elements[c2.indices[j]]


# ---------------------------------------------------------------------------
# Contact analytique : rotation exacte jusqu'au premier contact
# ---------------------------------------------------------------------------

def _cercle_segments(P, O, A, B):
    """Intersections des cercles de centre O passant par les points P (k, 2)
    avec les segments [A, B] (k, 2), paire par paire : angles polaires (k, 2)
    des points d'intersection (nan s'il n'y en a pas)."""
    rho2 = np.sum((P - O) ** 2, axis=1)
    d = B - A
    w = A - O
    a = np.sum(d * d, axis=1)
    b = 2 * np.sum(d * w, axis=1)
    c = np.sum(w * w, axis=1) - rho2
    disc = b * b - 4 * a * c
    ok = disc >= 0
    r = np.sqrt(np.where(ok, disc, 0.0))
    out = np.full(c.shape + (2,), np.nan)
    for k, sg in enumerate((-1.0, 1.0)):
        u = (-b + sg * r) / (2 * a)
        dans = ok & (u >= 0.0) & (u <= 1.0)
        X = w[:, 0] + u * d[:, 0]
        Y = w[:, 1] + u * d[:, 1]
        out[:, k] = np.where(dans, np.arctan2(Y, X), np.nan)
    return out


def _paires(V, S, O):
    """Paires (sommet i, segment j) telles que le cercle de centre O passant
    par le sommet coupe le segment : le rayon du sommet doit tomber entre les
    distances min et max du segment à O. Les sommets triés par rayon, chaque
    segment ne retient que sa tranche ; avec des arcs finement échantillonnés,
    c'est une infime partie du produit cartésien."""
    rho = np.hypot(*(V - O).T)
    rmin, rmax = _rayons_segments(S, O)
    ordre = np.argsort(rho)
    rs = rho[ordre]
    lo = np.searchsorted(rs, rmin - 1e-12, "left")
    hi = np.searchsorted(rs, rmax + 1e-12, "right")
    nb = np.maximum(hi - lo, 0)
    total = int(nb.sum())
    j = np.repeat(np.arange(len(S)), nb)
    i = ordre[np.arange(total) - np.repeat(np.cumsum(nb) - nb, nb) + np.repeat(lo, nb)]
    return i, j


def _premier_sommet(V, S, O, sens):
    """Plus petite rotation t >= 0 des sommets V autour de O (sens `sens`)
    qui en amène un sur l'un des segments S. Rend (t, i, j), t = inf si
    aucun ne l'atteint."""
    i, j = _paires(V, S, O)
    if not len(i):
        return np.inf, -1, -1
    th = _cercle_segments(V[i], O, S[j, 0], S[j, 1])
    thP = np.arctan2(V[i, 1] - O[1], V[i, 0] - O[0])[:, None]
    t = np.mod(sens * (th - thP), 2 * np.pi)
    if not np.isfinite(t).any():
        return np.inf, -1, -1
    k = np.nanargmin(t)
    return t.flat[k], i[k // 2], j[k // 2]


def _rayons_segments(S, O):
    """Distances min et max au point O de chaque segment (m, 2, 2)."""
    A, B = S[:, 0] - O, S[:, 1] - O
    d = B - A
    l2 = np.maximum(np.sum(d * d, axis=1), 1e-300)
    u = np.clip(-np.sum(A * d, axis=1) / l2, 0.0, 1.0)
    rmin = np.hypot(*(A + u[:, None] * d).T)
    rmax = np.maximum(np.hypot(*A.T), np.hypot(*B.T))
    return rmin, rmax


def _filtrer(V, S, O, sens, course):
    """Ne garde que les sommets V (tournant autour de O dans le sens `sens`,
    sur au plus `course`) qui peuvent atteindre l'un des segments S, et ces
    segments-là. Filtre sur le rayon, puis sur l'angle si la course est
    inférieure à un demi-tour."""
    iv, js = np.arange(len(V)), np.arange(len(S))
    if not len(V) or not len(S):
        return V[:0], S[:0], iv[:0], js[:0]
    rho = np.hypot(*(V - O).T)
    rmin, rmax = _rayons_segments(S, O)
    ks = (rmin <= rho.max()) & (rmax >= rho.min())
    S, rmin, rmax, js = S[ks], rmin[ks], rmax[ks], js[ks]
    if not len(S):
        return V[:0], S, iv[:0], js
    kv = (rho >= rmin.min()) & (rho <= rmax.max())
    V, iv = V[kv], iv[kv]
    if course < np.pi and len(V):
        # secteur angulaire couvert par les segments, vu de O
        ang = np.arctan2(S[..., 1] - O[1], S[..., 0] - O[0]).ravel()
        ref = ang[0]
        rel = np.mod(ang - ref + np.pi, 2 * np.pi) - np.pi
        lo, hi = rel.min(), rel.max()
        if hi - lo < np.pi:
            th = np.mod(np.arctan2(V[:, 1] - O[1], V[:, 0] - O[0]) - ref + np.pi, 2 * np.pi) - np.pi
            # le sommet balaie [th, th + sens.course] : doit croiser [lo, hi]
            a0 = np.where(sens > 0, th, th - course)
            a1 = np.where(sens > 0, th + course, th)
            marge = 1e-3
            garde = np.zeros(len(V), bool)
            for k in (-1, 0, 1):
                garde |= (a1 + 2 * np.pi * k >= lo - marge) & (a0 + 2 * np.pi * k <= hi + marge)
            V, iv = V[garde], iv[garde]
        # et réciproquement : les segments que les sommets peuvent atteindre
        if len(V):
            th = np.mod(np.arctan2(V[:, 1] - O[1], V[:, 0] - O[0]) - ref + np.pi, 2 * np.pi) - np.pi
            v_lo, v_hi = th.min(), th.max()
            if v_hi - v_lo < np.pi:
                if sens > 0:
                    v_hi = v_hi + course
                else:
                    v_lo = v_lo - course
                s_ang = np.mod(np.arctan2(S[..., 1] - O[1], S[..., 0] - O[0]) - ref + np.pi,
                               2 * np.pi) - np.pi
                s_lo, s_hi = s_ang.min(axis=1), s_ang.max(axis=1)
                large = s_hi - s_lo > np.pi          # segment à cheval sur ±pi : gardé
                k = large | ((s_hi >= v_lo - 1e-3) & (s_lo <= v_hi + 1e-3))
                S, js = S[k], js[k]
    return V, S, iv, js


def sommets(S):
    """Sommets d'un ensemble de segments : le k-ième est l'extrémité de
    départ (k < m) ou d'arrivée (k >= m) du segment k mod m."""
    return np.concatenate([S[:, 0], S[:, 1]])


def premier_contact(S_mob, O, sens, S_fix, course=2 * np.pi):
    """Plus petite rotation t >= 0 (rad) d'un corps mobile autour de O, dans
    le sens `sens`, qui amène son contour au contact d'un contour fixe.

    S_* : segments (m, 2, 2), dans la pose courante. Deux familles de
    contacts, toutes deux exactes pour des contours polygonaux : un sommet du
    mobile qui vient sur un segment fixe, ou un sommet fixe sur lequel vient
    un segment du mobile. Le corps doit partir d'une position libre.

    Retourne (t, P, n, j_fixe, j_mobile) : la rotation, le point de contact
    (pose tournée), la normale unitaire dirigée du fixe vers le mobile et les
    indices des segments (fixe et mobile) qui portent le contact ; t = inf
    si rien n'est touché sur `course`.
    """
    m_mob, m_fix = len(S_mob), len(S_fix)
    V_mob, V_fix = sommets(S_mob), sommets(S_fix)
    meilleur = (np.inf, None, None, -1, -1)
    Vm, Sf, ivm, jsf = _filtrer(V_mob, S_fix, O, sens, course)
    Vf, Sm, ivf, jsm = _filtrer(V_fix, S_mob, O, -sens, course)
    if len(Vm) and len(Sf):
        tk, i, j = _premier_sommet(Vm, Sf, O, sens)
        if np.isfinite(tk):
            P = rotation(Vm[i], O, sens * tk)
            sg = Sf[j, 1] - Sf[j, 0]
            n = np.array([-sg[1], sg[0]]) / np.hypot(*sg)
            if (rotation(P, O, -sens * 1e-6) - P) @ n < 0:
                n = -n
            meilleur = (tk, P, n, int(jsf[j]), int(ivm[i] % m_mob))
    if len(Vf) and len(Sm):
        tk, i, j = _premier_sommet(Vf, Sm, O, -sens)
        if tk < meilleur[0]:
            P = Vf[i].copy()
            a0 = rotation(Sm[j, 0], O, sens * tk)
            sg = rotation(Sm[j, 1], O, sens * tk) - a0
            n = np.array([-sg[1], sg[0]]) / np.hypot(*sg)
            if (rotation(P, O, -sens * 1e-6) - P) @ n < 0:
                n = -n
            meilleur = (tk, P, n, int(ivf[i] % m_fix), int(jsm[j]))
    if meilleur[0] > course:
        return np.inf, None, None, -1, -1
    return meilleur


# ---------------------------------------------------------------------------
# Contact exact sur les vrais arcs
# ---------------------------------------------------------------------------
#
# Un contour échantillonné en facettes donne des positions de contact exactes
# à la flèche près (quelques picomètres), mais pas des rapports de transmission
# exacts : à chaque sommet franchi, la normale tourne d'un pas d'échantillonnage
# et le rapport saute. En dynamique, ces sauts deviennent des pointes de couple
# (le « bruit » de la courbe de couple) et obligent l'intégrateur à des pas de
# l'ordre de la microseconde. On cherche donc le premier contact sur la
# géométrie du DXF elle-même : segments, arcs et leurs extrémités.
#
# Le premier contact d'un corps qui tourne est une configuration où les deux
# contours se touchent en un point. Il en existe de deux familles : un SOMMET
# (extrémité d'élément) sur un segment ou un arc de l'autre corps, et une
# TANGENCE (arc sur segment, arc sur arc). Toutes se résolvent en fermé : le
# sommet, ou le centre de l'arc mobile, décrit un cercle autour de l'axe, qu'on
# coupe avec une droite (décalée du rayon pour une tangence) ou un cercle. Le
# point de contact doit appartenir aux deux éléments. Le plus petit angle sur
# l'ensemble des candidats est le premier contact ; une solution « fausse »
# (tangence d'un arc concave, par exemple) n'est atteignable qu'après une
# interpénétration, donc jamais la plus petite : il n'y a pas besoin de connaître
# la convexité des arcs.

DEUX_PI = 2 * np.pi


@dataclass(eq=False)
class Primitives:
    """Contour exact d'un ensemble d'éléments, dans une pose donnée.

    Segments [A, B] (js : élément porteur), arcs de centre C, rayon R, partant
    de l'angle a0 avec un balayage signé bal (ja), sommets V (jv)."""
    A: np.ndarray
    B: np.ndarray
    js: np.ndarray
    C: np.ndarray
    R: np.ndarray
    a0: np.ndarray
    bal: np.ndarray
    ja: np.ndarray
    V: np.ndarray
    jv: np.ndarray

    def tournees(self, O, angle):
        if angle == 0.0:
            return self
        return Primitives(rotation(self.A, O, angle), rotation(self.B, O, angle), self.js,
                          rotation(self.C, O, angle), self.R, self.a0 + angle, self.bal,
                          self.ja, rotation(self.V, O, angle), self.jv)


def primitives(elements):
    """Primitives exactes d'une liste d'éléments (indices dans cette liste)."""
    A, B, js, C, R, a0, bal, ja, V, jv = [], [], [], [], [], [], [], [], [], []
    for k, e in enumerate(elements):
        p = e.points
        V += [p[0], p[-1]]
        jv += [k, k]
        if e.est_arc:
            c = np.asarray(e.centre, float)
            ang = np.unwrap(np.arctan2(p[:, 1] - c[1], p[:, 0] - c[0]))
            b = ang[-1] - ang[0]
            if np.hypot(*(p[-1] - p[0])) < 1e-9 * max(e.rayon, 1e-9) and len(p) > 2:
                b = DEUX_PI * np.sign(b or 1.0)          # cercle complet
            C.append(c)
            R.append(e.rayon)
            a0.append(ang[0])
            bal.append(b)
            ja.append(k)
        else:
            for q0, q1 in zip(p[:-1], p[1:]):
                A.append(q0)
                B.append(q1)
                js.append(k)
    arr = lambda x, d=2: np.array(x, float).reshape(-1, d) if d else np.array(x, float)
    return Primitives(arr(A), arr(B), np.array(js, int), arr(C), np.array(R, float),
                      np.array(a0, float), np.array(bal, float), np.array(ja, int),
                      arr(V), np.array(jv, int))


def _ang(v):
    return np.arctan2(v[..., 1], v[..., 0])


def _sur_arc(theta, a0, bal, tol=1e-9):
    """L'angle theta tombe-t-il dans l'arc (a0, balayage signé bal) ?"""
    m = np.mod(theta - a0 - 0.5 * bal + np.pi, DEUX_PI) - np.pi
    return np.abs(m) <= 0.5 * np.abs(bal) + tol


def _cercle_droite(O, rho, A, u):
    """Intersections du cercle (O, rho) avec les droites A + s.u (u unitaire).
    rho (n, 1), A et u (1, m, 2) : rend X (2, n, m, 2) et le masque (2, n, m)."""
    w = A - O
    b = np.sum(w * u, axis=-1)
    c = np.sum(w * w, axis=-1) - rho ** 2
    disc = b * b - c
    ok = disc >= 0
    r = np.sqrt(np.where(ok, disc, 0.0))
    X = np.stack([A + (-b - r)[..., None] * u, A + (-b + r)[..., None] * u])
    return X, np.stack([ok, ok])


def _cercle_cercle(O, rho, C, d):
    """Intersections du cercle (O, rho) avec les cercles (C, d).
    rho (n, 1), C (1, m, 2), d (1, m) : rend X (2, n, m, 2) et le masque."""
    w = C - O
    L = np.hypot(w[..., 0], w[..., 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        x = (rho ** 2 - d ** 2 + L ** 2) / (2 * L)
        h2 = rho ** 2 - x ** 2
        ok = (L > 1e-15) & (h2 >= 0) & (d > 1e-12)
        h = np.sqrt(np.where(ok, h2, 0.0))
        e = w / np.where(L > 0, L, 1.0)[..., None]
    f = np.stack([-e[..., 1], e[..., 0]], axis=-1)
    base = O + x[..., None] * e
    X = np.stack([base - h[..., None] * f, base + h[..., None] * f])
    return np.nan_to_num(X), np.stack([ok, ok])


def premier_contact_exact(mob, O, sens, fix, course=DEUX_PI):
    """Comme `premier_contact`, sur les primitives exactes (`Primitives`) :
    plus petite rotation t >= 0 du corps mobile autour de O, dans le sens
    `sens`, qui amène son contour au contact du contour fixe.

    Retourne (t, P, n, élément fixe, élément mobile) ; t = inf si rien n'est
    touché sur `course`. n est dirigée du fixe vers le mobile."""
    T, P, N, JF, JM = [], [], [], [], []

    def garder(theta, X, Nrm, ok, jf, jm):
        # theta : rotation appliquée au mobile ; t : comptée dans le sens demandé
        t = np.mod(sens * theta, DEUX_PI)
        ok = ok & (t <= course)
        if ok.any():
            T.append(t[ok])
            P.append(X[ok])
            N.append(Nrm[ok])
            JF.append(np.broadcast_to(jf, ok.shape)[ok])
            JM.append(np.broadcast_to(jm, ok.shape)[ok])

    def seg_unit(A, B):
        d = B - A
        L = np.hypot(d[:, 0], d[:, 1])
        return d / L[:, None], L

    # (a) sommet mobile sur segment fixe / (c) sur arc fixe
    if len(mob.V):
        rho = np.hypot(*(mob.V - O).T)[:, None]
        thV = _ang(mob.V - O)[:, None]
        if len(fix.A):
            u, L = seg_unit(fix.A, fix.B)
            X, ok = _cercle_droite(O, rho, fix.A[None], u[None])
            s = np.sum((X - fix.A[None]) * u[None], axis=-1)
            ok &= (s >= -1e-12) & (s <= L[None] + 1e-12)
            nrm = np.broadcast_to(np.stack([-u[:, 1], u[:, 0]], -1)[None, None], X.shape)
            garder(_ang(X - O) - thV, X, nrm, ok, fix.js[None, None], mob.jv[None, :, None])
        if len(fix.C):
            X, ok = _cercle_cercle(O, rho, fix.C[None], fix.R[None])
            v = X - fix.C[None]
            ok &= _sur_arc(_ang(v), fix.a0[None], fix.bal[None])
            garder(_ang(X - O) - thV, X, v / fix.R[None, None, :, None], ok,
                   fix.ja[None, None], mob.jv[None, :, None])
    # (b) sommet fixe sur segment mobile / (d) sur arc mobile : dans le repère
    # du mobile, le sommet fixe tourne de -theta
    if len(fix.V):
        rho = np.hypot(*(fix.V - O).T)[:, None]
        thV = _ang(fix.V - O)[:, None]
        Vf = np.broadcast_to(fix.V[None, :, None], (2, len(fix.V), 1, 2))
        if len(mob.A):
            u, L = seg_unit(mob.A, mob.B)
            X, ok = _cercle_droite(O, rho, mob.A[None], u[None])
            s = np.sum((X - mob.A[None]) * u[None], axis=-1)
            ok &= (s >= -1e-12) & (s <= L[None] + 1e-12)
            theta = thV - _ang(X - O)
            n0 = np.stack([-u[:, 1], u[:, 0]], -1)[None, None]
            c, sn = np.cos(theta), np.sin(theta)
            nrm = np.stack([c * n0[..., 0] - sn * n0[..., 1], sn * n0[..., 0] + c * n0[..., 1]], -1)
            garder(theta, np.broadcast_to(Vf, X.shape), nrm, ok,
                   fix.jv[None, :, None], mob.js[None, None])
        if len(mob.C):
            X, ok = _cercle_cercle(O, rho, mob.C[None], mob.R[None])
            v = X - mob.C[None]
            ok &= _sur_arc(_ang(v), mob.a0[None], mob.bal[None])
            theta = thV - _ang(X - O)
            c, sn = np.cos(theta), np.sin(theta)
            v = v / mob.R[None, None, :, None]
            nrm = np.stack([c * v[..., 0] - sn * v[..., 1], sn * v[..., 0] + c * v[..., 1]], -1)
            garder(theta, np.broadcast_to(Vf, X.shape), nrm, ok,
                   fix.jv[None, :, None], mob.ja[None, None])
    # (e) arc mobile tangent à un segment fixe : son centre, qui tourne autour
    # de O, doit se trouver à la distance R de la droite du segment
    if len(mob.C) and len(fix.A):
        rho = np.hypot(*(mob.C - O).T)[:, None]
        thC = _ang(mob.C - O)[:, None]
        u, L = seg_unit(fix.A, fix.B)
        nu = np.stack([-u[:, 1], u[:, 0]], -1)
        for eps in (-1.0, 1.0):
            off = eps * mob.R[:, None, None] * nu[None]                     # (n, m, 2)
            X, ok = _cercle_droite(O, rho, fix.A[None] + off, u[None])
            Pc = X - off[None]
            s = np.sum((Pc - fix.A[None]) * u[None], axis=-1)
            ok &= (s >= -1e-12) & (s <= L[None] + 1e-12)
            theta = _ang(X - O) - thC
            ok &= _sur_arc(_ang(Pc - X) - theta, mob.a0[None, :, None], mob.bal[None, :, None])
            garder(theta, Pc, np.broadcast_to(nu[None, None], X.shape), ok,
                   fix.js[None, None], mob.ja[None, :, None])
    # (f) arc fixe tangent à un segment mobile (repère du mobile)
    if len(fix.C) and len(mob.A):
        rho = np.hypot(*(fix.C - O).T)[:, None]
        thC = _ang(fix.C - O)[:, None]
        u, L = seg_unit(mob.A, mob.B)
        nu = np.stack([-u[:, 1], u[:, 0]], -1)
        for eps in (-1.0, 1.0):
            off = eps * fix.R[:, None, None] * nu[None]
            X, ok = _cercle_droite(O, rho, mob.A[None] + off, u[None])
            Pm = X - off[None]
            s = np.sum((Pm - mob.A[None]) * u[None], axis=-1)
            ok &= (s >= -1e-12) & (s <= L[None] + 1e-12)
            theta = thC - _ang(X - O)
            Pw = _tourner(Pm, O, theta)
            ok &= _sur_arc(_ang(Pw - fix.C[None, :, None]), fix.a0[None, :, None],
                           fix.bal[None, :, None])
            garder(theta, Pw, _tourner(np.broadcast_to(nu[None, None], X.shape), 0.0, theta), ok,
                   fix.ja[None, :, None], mob.js[None, None])
    # (g) arc mobile tangent à un arc fixe : distance des centres R1 + R2
    # (extérieurs) ou |R1 - R2| (l'un dans l'autre)
    if len(mob.C) and len(fix.C):
        rho = np.hypot(*(mob.C - O).T)[:, None]
        thC = _ang(mob.C - O)[:, None]
        Rm, Rf = mob.R[:, None], fix.R[None, :]
        for d, signe in ((Rm + Rf, 1.0), (np.abs(Rm - Rf), np.where(Rf >= Rm, 1.0, -1.0))):
            X, ok = _cercle_cercle(O, rho, fix.C[None], d)
            u = (X - fix.C[None]) / np.where(d > 0, d, 1.0)[None, ..., None]
            Pc = fix.C[None] + (signe * Rf)[None, ..., None] * u
            theta = _ang(X - O) - thC
            ok &= _sur_arc(_ang(Pc - fix.C[None]), fix.a0[None, None], fix.bal[None, None])
            ok &= _sur_arc(_ang(Pc - X) - theta, mob.a0[None, :, None], mob.bal[None, :, None])
            garder(theta, Pc, u, ok, fix.ja[None, None], mob.ja[None, :, None])
    if not T:
        return np.inf, None, None, -1, -1
    T = np.concatenate(T)
    k = int(np.argmin(T))
    Pk = np.concatenate(P)[k]
    nk = np.concatenate(N)[k]
    nk = nk / np.hypot(*nk)
    # normale du fixe vers le mobile : le mobile qui recule s'en écarte
    if (rotation(Pk, O, -sens * 1e-6) - Pk) @ nk < 0:
        nk = -nk
    return float(T[k]), Pk, nk, int(np.concatenate(JF)[k]), int(np.concatenate(JM)[k])


def _tourner(p, O, theta):
    """Rotation de points p (..., 2) autour de O par des angles theta (...)."""
    c, s = np.cos(theta), np.sin(theta)
    q = p - O
    return np.stack([c * q[..., 0] - s * q[..., 1], s * q[..., 0] + c * q[..., 1]], -1) + O
