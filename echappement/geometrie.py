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

#: flèche maximale tolérée en échantillonnant un arc (mm) : 0,01 µm
FLECHE_MAX = 1e-5
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
