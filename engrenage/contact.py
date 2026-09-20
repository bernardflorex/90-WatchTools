"""Contact entre les flancs de deux roues, passage d'une paire de dents et efforts.

Repris de Kontap.F90 / Enganap.F90 / Draweng.F90 (K. Trümpy) : les flancs sont des
lignes brisées (affinées par une spline naturelle) ; on cherche de combien il faut
tourner la roue menée, depuis une position éloignée, pour qu'elle touche le flanc
menant. Repère du dessin (voir evolvente.geometrie_engrenage) : point primitif en
(0,0), roue 1 (menante) de centre (0,-rt1), roue 2 (menée) de centre (0,rt2),
angles positifs dans le sens trigonométrique.
"""

import math
from dataclasses import dataclass

import numpy as np

NA = 100  # points par flanc pour la recherche grossière du contact
NF = 3000  # points par flanc pour l'affinage local (Enganap.F90 : NA=5000 partout)
FEN = 40  # demi-largeur (en points fins) de la zone affinée autour du contact grossier
TOL = 3e-5  # tolérance (rad) sur l'égalité des positions de deux paires de dents
SEUIL_PAIRE = 5e-4  # écart (rad) à partir duquel la paire 1 remplace la paire 0 pour les efforts
N_POS = 150  # positions de la menante sur une période (IK dans Enganap.F90)


class ErreurContact(ValueError):
    pass


def spline_naturelle(pts, na=NA):
    """`na` points équidistants (abscisse curviligne) sur la spline cubique naturelle
    passant par `pts`, comme SPLIPH avec AW(1)=1D33."""
    p = np.asarray(pts, dtype=float)
    p = p[np.concatenate([[True], np.hypot(*np.diff(p, axis=0).T) > 1e-12])]
    t = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    h = np.diff(t)
    n = len(t)
    b = np.zeros((n, 2))
    for i in range(1, n - 1):
        b[i] = 6 * ((p[i + 1] - p[i]) / h[i] - (p[i] - p[i - 1]) / h[i - 1])
    # système tridiagonal h[i-1]*m[i-1] + 2(h[i-1]+h[i])*m[i] + h[i]*m[i+1] = b[i], m[0]=m[-1]=0
    diag = np.zeros(n)
    m = np.zeros((n, 2))
    for i in range(1, n - 1):
        diag[i] = 2 * (h[i - 1] + h[i]) - (h[i - 1] ** 2 / diag[i - 1] if i > 1 else 0.0)
        b[i] -= (h[i - 1] * b[i - 1] / diag[i - 1]) if i > 1 else 0.0
    for i in range(n - 2, 0, -1):  # dérivées secondes
        m[i] = (b[i] - (h[i] * m[i + 1] if i < n - 2 else 0.0)) / diag[i]
    ts = np.linspace(0.0, t[-1], na)
    k = np.clip(np.searchsorted(t, ts, side="right") - 1, 0, n - 2)
    hk, u, v = h[k][:, None], (ts - t[k])[:, None], (t[k + 1] - ts)[:, None]
    return (m[k] * v**3 + m[k + 1] * u**3) / (6 * hk) + (p[k] / hk - m[k] * hk / 6) * v + (
        p[k + 1] / hk - m[k + 1] * hk / 6
    ) * u


def rotation(pts, angle, centre=(0.0, 0.0)):
    c, s = math.cos(angle), math.sin(angle)
    p = np.asarray(pts, dtype=float)
    return np.column_stack((p[:, 0] * c - p[:, 1] * s + centre[0], p[:, 0] * s + p[:, 1] * c + centre[1]))


def _angle(p, o):
    return np.arctan2(p[..., 1] - o[1], p[..., 0] - o[0])


def _tangentes(p):
    """Tangentes unitaires aux sommets de la ligne p (différences centrales)."""
    t = np.gradient(p, axis=0)
    return t / np.hypot(*t.T)[:, None]


def _tangente_en(tang, seg, t):
    v = (1 - t) * tang[seg] + t * tang[seg + 1]
    return v / np.hypot(*v)


def _cercle_segments(o, rayons, a, d):
    """Intersections des cercles (centre o, `rayons`) avec les segments a + t*d, 0<=t<=1.

    Produit, pour la racine t- puis t+, (indice du cercle, indice du segment, t) des seules
    intersections qui existent.
    """
    ao = a - o
    aa = np.einsum("ij,ij->i", d, d)
    bb = np.einsum("ij,ij->i", ao, d)
    cc = np.einsum("ij,ij->i", ao, ao)[None, :] - (rayons**2)[:, None]
    disc = bb[None, :] ** 2 - aa[None, :] * cc
    racine = np.sqrt(np.where(disc >= 0, disc, 0.0))
    for signe in (-1.0, 1.0):
        t = (-bb[None, :] + signe * racine) / aa[None, :]
        ok = (disc >= 0) & (t >= -1e-12) & (t <= 1 + 1e-12)
        ir, js = np.nonzero(ok)
        yield ir, js, np.clip(t[ir, js], 0.0, 1.0)


def premier_contact(c1, c2, o2, sens=1):
    """Plus petite rotation d in [0, 2π[ de `c2` autour de `o2` (sens = +1: trigonométrique)
    qui amène `c2` au contact du flanc fixe `c1`, comme KONTAP.

    Le contact est soit un sommet de c2 sur un segment de c1, soit un sommet de c1 sur un
    segment de c2 (rotation inverse). Retourne (d, point de contact, tangente commune au
    contact) ou (None, None, None). La tangente est celle du flanc sur lequel tombe le
    contact, pas celle du sommet (qui peut être un coin).
    """
    c1, c2 = np.asarray(c1), np.asarray(c2)
    o2 = np.asarray(o2, dtype=float)
    deux_pi = 2 * math.pi
    meilleur, point, tangente = math.inf, None, None

    # sommets de c2 sur les segments de c1
    v = c2
    ang_v = _angle(v, o2)
    a, dd = c1[:-1], np.diff(c1, axis=0)
    for ir, js, t in _cercle_segments(o2, np.hypot(*(v - o2).T), a, dd):
        if not len(t):
            continue
        q = a[js] + t[:, None] * dd[js]
        d = (sens * (_angle(q, o2) - ang_v[ir])) % deux_pi
        k = int(np.argmin(d))
        if d[k] < meilleur:
            meilleur, point = d[k], q[k]
            tangente = _tangente_en(_tangentes(c1), js[k], t[k])

    # sommets de c1 sur les segments de c2 (c2 tourne de d <=> le sommet tourne de -d)
    u = c1
    ang_u = _angle(u, o2)
    a, dd = c2[:-1], np.diff(c2, axis=0)
    for ir, js, t in _cercle_segments(o2, np.hypot(*(u - o2).T), a, dd):
        if not len(t):
            continue
        q = a[js] + t[:, None] * dd[js]
        d = (sens * (ang_u[ir] - _angle(q, o2))) % deux_pi
        k = int(np.argmin(d))
        if d[k] < meilleur:
            meilleur, point = d[k], u[ir[k]]
            g = _tangente_en(_tangentes(c2), js[k], t[k])  # flanc mené avant rotation
            c, s_ = math.cos(sens * d[k]), math.sin(sens * d[k])
            tangente = np.array([g[0] * c - g[1] * s_, g[0] * s_ + g[1] * c])

    if point is None:
        return None, None, None
    return meilleur, point, tangente


def _croix(a, b):
    return a[0] * b[1] - a[1] * b[0]


def efforts(o1, o2, p, g, mu):
    """Rapport de transmission, rapport de couple, rendement et force normale pour un
    contact en `p` de tangente `g` (Draweng.F90 : dWg/dWt, Mg/Mt, produit des deux).

    Rapport de transmission w2/w1 (négatif : les roues tournent en sens opposés).
    Rapport de couple M2/M1 : couple exercé sur la menée / couple moteur, avec frottement
    de Coulomb `mu` opposé à la vitesse de glissement. Rendement = (w2/w1)*(M2/M1).
    Force normale pour un couple moteur de 1 (N par N.mm, donc N pour 1 mN.m).
    """
    r1, r2 = p - o1, p - o2
    n = np.array([-g[1], g[0]])
    perp1, perp2 = np.array([-r1[1], r1[0]]), np.array([-r2[1], r2[0]])
    rho = n @ perp1 / (n @ perp2)  # les vitesses normales de contact sont égales
    glissement = (rho * perp2 - perp1) @ g  # vitesse de la menée par rapport à la menante
    if _croix(r1, n) < 0:  # normale dirigée de la menante vers la menée : couple moteur > 0
        n = -n
    d = n - mu * np.sign(glissement) * g  # direction de la force exercée sur la menée
    m1 = _croix(r1, d)
    if m1 <= 0:  # arc-boutement
        return rho, math.nan, math.nan, math.nan
    couple = _croix(r2, d) / m1
    return rho, couple, rho * couple, 1.0 / m1


@dataclass
class Position:
    """Une position de la menante dans la période."""

    alpha1: float  # angle de la menante (rad) ; 0 = dent centrée sur la ligne des centres
    alpha2: float  # angle de l'axe de la dent 0 de la menée (rad) ; pi - pi/Z2 = creux centré
    contact: tuple  # point d'engrènement
    contact2: tuple = None  # second contact, pendant le recouvrement de deux paires de dents


@dataclass
class Simulation:
    """Une période d'engrènement : de l'entrée d'une paire de dents à celle de la suivante."""

    mu: float
    alpha1: np.ndarray  # (n,) angle de la menante
    alpha2: np.ndarray  # (n,) angle de la menée
    contact: np.ndarray  # (n,2) point d'engrènement
    contact2: np.ndarray  # (n,2) second contact (nan hors recouvrement)
    paire: np.ndarray  # (n,) paire de dents utilisée (0 ou 1)
    rapport: np.ndarray  # (n,) w2/w1
    couple: np.ndarray  # (n,) M2/M1
    rendement: np.ndarray  # (n,)
    normale: np.ndarray  # (n,) force normale pour un couple moteur de 1 (N par mN.m)
    ideal_rapport: float  # -Z1/Z2
    ideal_couple: float  # -Z2/Z1
    jeu_min: float = math.nan  # jeu de la menée sur la période (degrés)
    jeu_max: float = math.nan

    def __len__(self):
        return len(self.alpha1)

    def pourcent(self):
        """Les trois paramètres en % de leur valeur idéale (100 % = engrènement parfait)."""
        return {
            "rapport": 100.0 * self.rapport / self.ideal_rapport,
            "couple": 100.0 * self.couple / self.ideal_couple,
            "rendement": 100.0 * self.rendement,
        }

    def statistiques(self):
        """Par paramètre : (min, moyenne, max, variation) en %, variation = 2(max-min)/(max+min)."""
        stats = {}
        for nom, v in self.pourcent().items():
            v = v[np.isfinite(v)]
            mn, mx = float(v.min()), float(v.max())
            stats[nom] = (mn, float(v.mean()), mx, 200.0 * (mx - mn) / (mx + mn))
        return stats

    def position(self, i):
        c2 = None if np.isnan(self.contact2[i, 0]) else tuple(self.contact2[i])
        return Position(self.alpha1[i], self.alpha2[i], tuple(self.contact[i]), c2)


class Engrenement:
    def __init__(self, roues, profils, da=0.0, na=NA):
        """`da` : variation de l'entraxe par rapport à rt1 + rt2 (mm). Les cercles primitifs de
        fonctionnement (rayons proportionnels aux nombres de dents) sont tangents en (0,0)."""
        self.roues = roues
        z1, z2 = roues[0].zz, roues[1].zz
        r1, r2 = roues[0].mo * z1 / 2.0, roues[1].mo * z2 / 2.0
        entraxe = r1 + r2 + da
        self._tester_entraxe(profils, entraxe)
        cy = entraxe / (r1 + r2)  # rapport des cercles primitifs de fonctionnement aux nominaux
        r1, r2 = r1 * cy, r2 * cy
        self.a1, self.a2 = 2 * math.pi / z1, 2 * math.pi / z2
        self.o1, self.o2 = np.array([0.0, -r1]), np.array([0.0, r2])
        # flancs menants (côté x<0), en résolution grossière puis fine
        gauche = [[(-x, y) for x, y in prof] for prof in profils]
        self.g1, self.g2 = (spline_naturelle(g, na) for g in gauche)
        self.f1, self.f2 = (spline_naturelle(g, NF) for g in gauche)
        self.ref2 = rotation(self.g2, 0.0, self.o2)
        self.ref2f = rotation(self.f2, 0.0, self.o2)
        self._cache = {}

    @staticmethod
    def _tester_entraxe(profils, aa):
        """Entraxe trop petit (pied contre tête) ou trop grand (plus de contact), comme les
        tests FI=1 et FI=2 de Enganap.F90 : rayon de pied = 1er point, rayon de tête = dernier."""
        pied = [math.hypot(*p[0]) for p in profils]
        tete = [p[-1][1] for p in profils]  # la tête est sur l'axe de la dent (x = 0)
        if pied[0] + tete[1] > aa or pied[1] + tete[0] > aa:
            raise ErreurContact("entraxe trop petit (la tête d'une roue touche le pied de l'autre)")
        if abs((aa * aa + tete[0] ** 2 - tete[1] ** 2) / (2 * aa * tete[0])) > 1:
            raise ErreurContact("entraxe trop grand (les cercles de tête ne se coupent plus)")

    def psi(self, phi):
        """(rotation de la menée, point, tangente) pour la dent de la menante d'axe phi."""
        if phi not in self._cache:
            self._cache[phi] = self._contact(phi)
        return self._cache[phi]

    def _contact(self, phi):
        # recherche grossière, puis même recherche sur les seuls morceaux de flanc voisins
        # du contact, très finement échantillonnés : position et tangente du contact
        # bien plus précises qu'avec des cordes de 1e-3 (erreur en racine de la flèche)
        d0, q0, g0 = premier_contact(rotation(self.g1, phi, self.o1), self.ref2, self.o2, +1)
        if d0 is None:
            return math.inf, None, None
        f1 = rotation(self.f1, phi, self.o1)
        i1 = int(np.argmin(np.hypot(*(f1 - q0).T)))
        i2 = int(np.argmin(np.hypot(*(rotation(self.f2, d0, self.o2) - q0).T)))
        s1, s2 = slice(max(0, i1 - FEN), i1 + FEN + 1), slice(max(0, i2 - FEN), i2 + FEN + 1)
        d, q, g = premier_contact(f1[s1], self.ref2f[s2], self.o2, +1)
        if d is None or abs(d - d0) > 1e-3:  # contact réel hors de la zone affinée
            return d0, q0, g0
        return d, q, g

    def theta_max(self, phi):
        """Position la plus haute de la menée quand la menante est en phi (contact flanc menant)."""
        return min(self._c(phi, k) for k in (-1, 0, 1))

    def jeu(self, debut, n=10):
        """Jeu de la menée (degrés) sur une période : (min, max), comme AW(1) et AW(2).

        Par symétrie par rapport à la ligne des centres, la position la plus basse de la menée
        avec la menante en phi est l'image de la plus haute avec la menante en -phi ; l'écart
        entre les deux, replié modulo un pas de la menée, est le jeu.
        """
        jeux = []
        for j in range(n):
            phi = debut + j * self.a1 / (n - 1)
            somme = self.theta_max(phi) + self.theta_max(-phi)
            jeux.append(math.degrees((somme + self.a2 / 2) % self.a2 - self.a2 / 2))
        return min(jeux), max(jeux)

    def _c(self, phi, k):
        """Position de la menée imposée par la paire k (dent k de la menante)."""
        return self.psi(phi + k * self.a1)[0] + k * self.a2

    def debut(self):
        """Angle de la menante au début de menée, comme la recherche de FW dans Enganap.F90.

        Pour la paire k (dent k de la menante, dent -k de la menée) le flanc menant touche
        le flanc mené quand l'axe de la dent menée est à psi(alpha1 + k*A1) + k*A2, avec
        psi(phi) = premier contact en tournant la menée depuis une position éloignée. La
        menée se place sur la plus petite de ces valeurs. La paire 0 entre en contact
        quand sa valeur rejoint celle de la paire 1 (encore en prise).
        """

        def ecart(phi):  # > 0 : la paire 0 n'est pas encore en contact
            c0, c1 = self._c(phi, 0), self._c(phi, 1)
            if math.isinf(c0):
                return math.inf
            return -math.inf if math.isinf(c1) else c0 - c1

        # balayage croissant depuis -2 pas jusqu'à la première position où la paire 0 touche
        pas, n = self.a1 / 8, 8
        grille = [(k - 2 * n) * pas for k in range(2 * n + 1)]
        if ecart(grille[0]) <= TOL:
            raise ErreurContact("début de menée introuvable (paire 0 déjà en contact)")
        for lo, hi in zip(grille, grille[1:]):
            if ecart(hi) <= TOL:
                break
        else:
            raise ErreurContact("début de menée introuvable (pas de contact continu ?)")
        while hi - lo > 1e-7:
            m = (lo + hi) / 2
            lo, hi = (m, hi) if ecart(m) > TOL else (lo, m)

        # les deux paires doivent vraiment se toucher (sinon rapport de conduite < 1)
        theta2 = min(self._c(hi, 0), self._c(hi, 1))
        limite = 5e-3 * self.roues[0].mo
        for k in (0, 1):
            q = self.psi(hi + k * self.a1)[1]
            if _distance(q, rotation(self.g2, theta2 - k * self.a2, self.o2)) > limite:
                raise ErreurContact("deux paires de dents ne sont jamais en prise simultanément")
        return hi

    def simuler(self, mu=0.2, n_pos=N_POS):
        """`n_pos` positions de la menante, du début de menée à l'entrée de la paire suivante.

        Dans le recouvrement de deux paires on garde la paire 0 (ligne d'engrènement
        continue), sauf si la paire 1 impose une position de la menée nettement plus basse
        (profils qui ne se conjuguent pas) ; la menée reste sur la plus basse des deux.
        """
        debut = self.debut()
        n = n_pos
        res = {k: np.full(n, np.nan) for k in ("a1", "a2", "rap", "cou", "ren", "nor", "paire")}
        contact, contact2 = np.full((n, 2), np.nan), np.full((n, 2), np.nan)
        sortie = False  # la paire 1 a quitté le contact, elle n'y reviendra pas
        for i in range(n):
            phi = debut + i * self.a1 / (n - 1)
            c0, q0, g0 = self.psi(phi)
            c1 = math.inf
            if not sortie:
                d1, q1, g1 = self.psi(phi + self.a1)
                c1 = d1 + self.a2
                sortie = not math.isfinite(c1) or c1 > c0 + 1e-3
            paire = 1 if (not sortie and c1 < c0 - SEUIL_PAIRE) else 0
            q, g = (q1, g1) if paire else (q0, g0)
            contact[i] = q
            if not sortie and abs(c1 - c0) <= 2e-4:  # les deux paires se touchent
                contact2[i] = q0 if paire else q1
            res["a1"][i], res["a2"][i], res["paire"][i] = phi, min(c0, c1), paire
            res["rap"][i], res["cou"][i], res["ren"][i], res["nor"][i] = efforts(self.o1, self.o2, q, g, mu)
        z1, z2 = self.roues[0].zz, self.roues[1].zz
        jeu_min, jeu_max = self.jeu(debut)
        return Simulation(mu, res["a1"], res["a2"], contact, contact2, res["paire"],
                          res["rap"], res["cou"], res["ren"], res["nor"], -z1 / z2, -z2 / z1,
                          jeu_min, jeu_max)


def _distance(p, ligne):
    """Distance du point p à la ligne brisée `ligne`."""
    a, d = ligne[:-1], np.diff(ligne, axis=0)
    t = np.clip(((p - a) * d).sum(1) / np.maximum((d * d).sum(1), 1e-30), 0.0, 1.0)
    return float(np.min(np.hypot(*(a + t[:, None] * d - p).T)))


def simuler(roues, profils, mu=0.2, n_pos=N_POS, da=0.0):
    return Engrenement(roues, profils, da).simuler(mu, n_pos)
