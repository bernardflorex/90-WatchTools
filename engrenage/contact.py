"""Contact entre les flancs de deux roues, passage d'une paire de dents et efforts.

Repris de Kontap.F90 / Enganap.F90 / Draweng.F90 (K. Trümpy) : les flancs sont des
lignes brisées (affinées par une spline naturelle) ; on cherche de combien il faut
tourner la roue menée, depuis une position éloignée, pour qu'elle touche le flanc
menant. Repère du dessin (voir evolvente.geometrie_engrenage) : point primitif en
(0,0), roue 1 (menante) de centre (0,-rt1), roue 2 (menée) de centre (0,rt2),
angles positifs dans le sens trigonométrique.
"""

import math
from dataclasses import dataclass, replace

import numpy as np

import arcs as arcs_mod

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


DEMI_COURBURE = 25  # demi-fenêtre (en points fins) du cercle ajusté pour la courbure locale


def rayon_courbure(ligne, i, demi=DEMI_COURBURE):
    """Rayon de courbure signé de la ligne brisée `ligne` en son point `i` : positif si le flanc
    est convexe, négatif s'il est concave -- la convention de `flanc_arcs`, pour que la pression
    de Hertz se calcule de la même façon avec des points et avec des arcs.

    Cercle des moindres carrés (Kåsa) sur une fenêtre de +-`demi` points. Sur un flanc
    échantillonné à partir d'une chaîne d'arcs, il redonne le rayon exact de l'arc à 5 chiffres
    (vérifié sur les 12 arcs de `validation.json`, hors arcs trop courts pour la fenêtre).
    """
    p = np.asarray(ligne, dtype=float)[max(0, i - demi):i + demi + 1]
    if len(p) < 5:
        return math.nan
    x, y = p[:, 0] - p[:, 0].mean(), p[:, 1] - p[:, 1].mean()
    try:
        d, e, f = np.linalg.lstsq(np.column_stack((x, y, np.ones_like(x))), x * x + y * y,
                                  rcond=None)[0]
    except np.linalg.LinAlgError:
        return math.nan
    r = math.sqrt(max(f + (d * d + e * e) / 4.0, 0.0))
    if r < 1e-12:
        return math.nan
    j = len(p) // 2
    tx, ty = p[-1] - p[0]  # sens de parcours
    vx, vy = d / 2.0 - x[j], e / 2.0 - y[j]  # du point vers le centre du cercle
    return -math.copysign(r, tx * vy - ty * vx)


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


def _balayage(a):
    """Angle ramené dans ]-pi, pi]."""
    return (a + math.pi) % (2 * math.pi) - math.pi


def flanc_arcs(mobile, dos=False, miroir=-1.0):
    """Flanc menant d'un mobile défini par des arcs de cercle : dictionnaire de centres,
    rayons, convexité, angle de départ et balayage de chaque arc, dans le repère de la roue.

    Le flanc menant est l'image miroir (x -> -x) du profil. Un arc est convexe (matière du côté
    du centre) quand son rayon signé est positif -- au signe près de `arcs.convexite`, qui rétablit
    la convexité physique pour une couronne intérieure (matière à l'extérieur, flanc actif
    concave). `dos` : le demi-profil arrière (dent dissymétrique), traité de la même façon, pour
    le jeu.
    """
    pts = np.asarray(mobile.profil_dos() if dos else mobile.points, dtype=float)
    pm, ra = mobile.arcs_dos() if dos else mobile.arcs()
    ra = arcs_mod.convexite(mobile.zz, ra)
    # un miroir est une réflexion : il ne change pas la convexité, seulement les coordonnées
    m = np.array([miroir, 1.0])
    pts, pm = pts * m, pm * m
    a0 = np.arctan2(pts[:-1, 1] - pm[:, 1], pts[:-1, 0] - pm[:, 0])
    a1 = np.arctan2(pts[1:, 1] - pm[:, 1], pts[1:, 0] - pm[:, 0])
    return {"centres": pm, "rayons": np.abs(ra), "convexe": ra > 0, "a0": a0,
            "balayage": _balayage(a1 - a0)}


def _flanc_place(fl, angle, centre):
    """Le flanc tourné de `angle` autour de `centre` (repère du dessin)."""
    return {**fl, "centres": rotation(fl["centres"], angle, centre), "a0": fl["a0"] + angle}


def _dans_arc(theta, a0, balayage, tol=1e-9):
    t = _balayage(theta - a0) * np.sign(balayage)
    return (t >= -tol) & (t <= np.abs(balayage) + tol)


def premier_contact_arcs(c1, c2, o2, sens=1.0):
    """Comme `premier_contact`, mais entre deux flancs faits d'arcs de cercle (KONTA/KONTW du
    Fortran) : plus petite rotation d in [0, 2π[ de `c2` autour de `o2` (sens trigonométrique)
    qui met un arc de `c2` en tangence avec un arc du flanc fixe `c1`.

    `sens` : sens de la rotation cherchée, comme dans `premier_contact` (+1 trigonométrique,
    -1 horaire pour un engrènement intérieur). `d` est toujours rendu positif, compté dans ce
    sens-là ; la rotation réellement appliquée est `sens * d`.

    Deux arcs sont tangents quand la distance de leurs centres vaut r1 + r2 (deux arcs
    convexes) ou r_concave - r_convexe (un convexe dans un concave) ; deux arcs concaves ne se
    touchent pas -- une denture intérieure, dont le flanc actif est concave, passe donc par le
    cas « convexe dans concave », déjà prévu. Le centre de l'arc de `c2` décrit un cercle autour de `o2` : la distance
    voulue donne d par la loi des cosinus. Le point de tangence, sur la ligne des centres,
    doit se trouver sur les deux arcs. Retourne (d, point, tangente, (rayon 1, rayon 2)) ou
    (None, None, None, None) ; les rayons de courbure au contact sont positifs pour un arc
    convexe, négatifs pour un arc concave.
    """
    o2 = np.asarray(o2, dtype=float)
    c1c, r1, v1 = c1["centres"], c1["rayons"][:, None], c1["convexe"][:, None]
    c2c, r2, v2 = c2["centres"], c2["rayons"][None, :], c2["convexe"][None, :]
    dt = np.where(v1 & v2, r1 + r2, np.where(v1 & ~v2, r2 - r1, np.where(~v1 & v2, r1 - r2, -1.0)))
    w_vec = c1c - o2
    w, aw = np.hypot(*w_vec.T)[:, None], np.arctan2(w_vec[:, 1], w_vec[:, 0])[:, None]
    rho_vec = c2c - o2
    rho, a02 = np.hypot(*rho_vec.T)[None, :], np.arctan2(rho_vec[:, 1], rho_vec[:, 0])[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        h = (rho**2 + w**2 - dt**2) / (2 * rho * w)
    possible = (dt > 1e-12) & (np.abs(h) <= 1.0 + 1e-12)
    ac = np.arccos(np.clip(np.where(possible, h, 0.0), -1.0, 1.0))
    meilleur, resultat = math.inf, None
    for signe in (1.0, -1.0):
        d = (sens * (aw - a02 + signe * ac)) % (2 * math.pi)
        theta = sens * d  # rotation réellement appliquée à c2
        c2r = o2 + rho[..., None] * np.stack((np.cos(a02 + theta), np.sin(a02 + theta)), axis=-1)
        u = (c2r - c1c[:, None, :]) / np.where(possible, dt, 1.0)[..., None]
        s = np.where(v1 & ~v2, -1.0, 1.0)  # convexe dans concave : tangence du côté opposé
        p = c1c[:, None, :] + (s * r1)[..., None] * u
        th1 = np.arctan2(p[..., 1] - c1c[:, None, 1], p[..., 0] - c1c[:, None, 0])
        th2 = np.arctan2(p[..., 1] - c2r[..., 1], p[..., 0] - c2r[..., 0]) - theta  # avant rotation
        ok = (possible & _dans_arc(th1, c1["a0"][:, None], c1["balayage"][:, None])
              & _dans_arc(th2, c2["a0"][None, :], c2["balayage"][None, :]))
        if ok.any():
            dv = np.where(ok, d, np.inf)
            i = np.unravel_index(np.argmin(dv), dv.shape)
            if dv[i] < meilleur:
                meilleur = dv[i]
                rk1 = c1["rayons"][i[0]] * (1.0 if c1["convexe"][i[0]] else -1.0)
                rk2 = c2["rayons"][i[1]] * (1.0 if c2["convexe"][i[1]] else -1.0)
                resultat = (p[i], np.array([-u[i][1], u[i][0]]), (rk1, rk2))
    if resultat is None:
        return None, None, None, None
    return meilleur, resultat[0], resultat[1], resultat[2]


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
    hertz: np.ndarray = None  # (n,) pression de Hertz (N/mm²) pour un couple moteur de 1 mN.m
    hertz_local: tuple = None  # ((x, y) sur la menante, (x, y) sur la menée) de la pression maximale
    hertz_indice: int = -1  # position où elle est atteinte

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

    def blocage(self, seuil=0.0):
        """Mesures du sens inverse : `(rendement_max, fraction_bloquee, blocage_initial)`.

        Une position est **bloquée** quand le rendement y est `<= seuil` (0 par défaut) ou n'est pas
        fini : le couple de sortie est alors de signe contraire au couple d'entrée, ou le contact
        arc-boute (`efforts` rend nan) -- dans les deux cas la paire ne transmet pas.

        - `rendement_max` : le meilleur rendement de toute la période. C'est lui qu'il faut
          regarder pour « rendement faible **partout** » : une seule position à bon rendement
          suffit à rendre l'engrenage réversible.
        - `fraction_bloquee` : part des positions bloquées.
        - `blocage_initial` : part des positions bloquées **d'affilée depuis le début** de la
          période. Un blocage groupé en début de cycle vaut mieux qu'un blocage éparpillé : il
          empêche le mouvement de s'amorcer, alors que des positions bloquées isolées sont
          franchies sur l'élan.
        """
        r = np.asarray(self.rendement, dtype=float)
        n = len(r)
        bloque = ~np.isfinite(r) | (r <= seuil)
        fini = r[np.isfinite(r)]
        debut = int(np.argmin(bloque)) if bloque.any() else 0
        if bloque.all():
            debut = n
        elif not bloque[0]:
            debut = 0
        return (float(fini.max()) if fini.size else -math.inf,
                float(bloque.sum()) / n, float(debut) / n)

    def hertz_max(self):
        """(pression maximale en N/mm², indice de position), ou (nan, -1) sans données de Hertz."""
        if self.hertz is None or not np.isfinite(self.hertz).any():
            return math.nan, -1
        i = int(np.nanargmax(self.hertz))
        return float(self.hertz[i]), i

    def criteres(self):
        """Les quatre critères d'optimisation : variation de couple (%, à minimiser), jeu maximal
        (degrés, à minimiser), rendement moyen (%, à maximiser), pression de Hertz maximale
        (N/mm² pour 1 mN.m, à minimiser)."""
        st = self.statistiques()
        return {"couple": st["couple"][3], "jeu": self.jeu_max, "rendement": st["rendement"][1],
                "hertz": self.hertz_max()[0]}

    def position(self, i):
        c2 = None if np.isnan(self.contact2[i, 0]) else tuple(self.contact2[i])
        return Position(self.alpha1[i], self.alpha2[i], tuple(self.contact[i]), c2)


class Engrenement:
    def __init__(self, roues, profils, da=0.0, na=NA, mobiles=None):
        """`da` : variation de l'entraxe par rapport à rt1 + rt2 (mm). Les cercles primitifs de
        fonctionnement (rayons proportionnels aux nombres de dents) sont tangents en (0,0).

        `mobiles` (arcs.Mobile, menant et mené) : les profils sont des arcs de cercle et le
        contact est calculé analytiquement entre arcs (`premier_contact_arcs`), sans aucun
        point de profil ; `profils` sert alors seulement au contrôle d'entraxe. Sinon les
        profils sont des lignes brisées (points de evolvente.evolve), affinées par spline.
        """
        self.roues = roues
        z1, z2 = roues[0].zz, roues[1].zz
        r1, r2 = roues[0].mo * z1 / 2.0, roues[1].mo * z2 / 2.0
        # Entraxe signé : négatif si l'une des roues est une couronne intérieure (les deux
        # centres sont alors du même côté du point primitif). `da` s'ajoute *algébriquement*,
        # comme `evolvente.entraxe_signe` : voir là-bas pourquoi ce n'est pas la variation de
        # l'entraxe géométrique pour une couronne.
        entraxe = (r1 + r2) + da
        self._tester_entraxe(profils, entraxe, z1 * z2 < 0)
        cy = entraxe / (r1 + r2)  # rapport des cercles primitifs de fonctionnement aux nominaux
        r1, r2 = r1 * cy, r2 * cy
        # Pas angulaires *toujours positifs* : la menante avance dans le sens trigonométrique
        # quel que soit le signe de son nombre de dents, et les positions de la menée sont comptées
        # dans le sens de recherche (`sens`), où son avance par dent de la menante vaut 2pi/|z2|.
        self.a1, self.a2 = 2 * math.pi / abs(z1), 2 * math.pi / abs(z2)
        # `psi` cherche la rotation de la menée qui la ramène contre la menante en la faisant
        # tourner *à l'envers de son mouvement* : sens trigonométrique si elle tourne dans
        # l'autre sens que la menante (engrènement extérieur), horaire sinon (intérieur, les
        # deux roues tournent du même côté). Toutes les positions de la menée sont donc
        # comptées dans ce sens-là (`a2u`), et converties en angle réel par `self.sens`.
        self.sens = -1.0 if z1 * z2 < 0 else 1.0
        self.a2u = self.a2
        self.o1, self.o2 = np.array([0.0, -r1]), np.array([0.0, r2])
        self.interne = z1 * z2 < 0
        # Côté du demi-profil qui porte le flanc actif, dans le repère propre de chaque roue.
        # La menante tourne dans le sens trigonométrique, donc au point primitif elle pousse vers
        # -x si son centre est en dessous (zz > 0) et vers +x s'il est au-dessus (couronne) ; un
        # flanc pousse à l'opposé de sa matière, d'où `miroir1 = -signe(z1)`. La menée présente le
        # flanc opposé, et son repère est retourné de 180° quand son centre est de l'autre côté du
        # point primitif (z2 > 0), d'où `miroir2 = -signe(z1*z2)`.
        self.miroir1 = -1.0
        self.miroir2 = -math.copysign(1.0, z1 * z2)
        # les dents doivent pointer vers le point primitif : le repère propre de la menante est
        # retourné d'un demi-tour quand son centre est de l'autre côté (couronne menante). Celui
        # de la menée n'a pas besoin de l'être, la recherche trouve librement sa rotation.
        self.delta1 = math.pi if z1 < 0 else 0.0
        self.arcs = mobiles is not None
        self._cache = {}
        self._courbures = {}  # rayons de courbure signés au contact, par angle de menante (arcs)
        # données de Hertz : celles des mobiles en mode arcs, celles des roues en mode points
        # (les .Cev en ont un bloc, les .json une clé "hertz")
        self.hertz = None
        porteurs = mobiles if mobiles is not None else roues
        if all(getattr(p, "e_mod", None) and p.poisson and p.epaisseur for p in porteurs):
            self.hertz = [(p.e_mod, p.poisson, p.epaisseur) for p in porteurs]
        if self.arcs:
            self.tol, self.seuil_paire, self.seuil_contact2 = 1e-9, 0.0, 1e-3
            self.fl1 = flanc_arcs(mobiles[0], miroir=self.miroir1)
            fl2 = flanc_arcs(mobiles[1], miroir=self.miroir2)
            self.ref2 = _flanc_place(fl2, 0.0, self.o2)
            # flancs arrière (jeu) : les mêmes si les dents sont symétriques
            self.fb1 = flanc_arcs(mobiles[0], dos=True, miroir=self.miroir1)
            fb2 = flanc_arcs(mobiles[1], dos=True, miroir=self.miroir2)
            self.ref2b = _flanc_place(fb2, 0.0, self.o2)
            self._cache_dos = {}
            # flanc mené échantillonné, seulement pour vérifier les contacts trouvés
            self.g2 = mobiles[1].echantillonner() * np.array([self.miroir2, 1.0])
            self._reference()
            return
        self.tol, self.seuil_paire, self.seuil_contact2 = TOL, SEUIL_PAIRE, 2e-4
        # flancs menants : celui de la menante est toujours le demi-profil mirroré (x<0) ; celui
        # de la menée aussi tant que l'engrènement est extérieur (elle est alors tournée de ~180°,
        # ce qui l'amène face à la menante). Dans un engrènement intérieur la couronne n'est pas
        # retournée : c'est son demi-profil *non mirroré* qui présente sa matière du bon côté.
        gauche = [[(s * x, y) for x, y in prof] for s, prof in zip((self.miroir1, self.miroir2), profils)]
        self.g1, self.g2 = (spline_naturelle(g, na) for g in gauche)
        self.f1, self.f2 = (spline_naturelle(g, NF) for g in gauche)
        self.ref2 = rotation(self.g2, 0.0, self.o2)
        self.ref2f = rotation(self.f2, 0.0, self.o2)
        self._reference()

    @staticmethod
    def _tester_entraxe(profils, aa, interne=False):
        """Entraxe trop petit (pied contre tête) ou trop grand (plus de contact), comme les
        tests FI=1 et FI=2 de Enganap.F90 : rayon de pied = 1er point, rayon de tête = dernier.

        Pour une denture intérieure le pied est au *grand* rayon et la tête au petit, et l'entraxe
        signé est négatif : le pignon doit tenir dans le cercle de pied de la couronne, la tête de
        la couronne doit dégager le pied du pignon, et les deux cercles de tête doivent encore se
        couper pour qu'il y ait contact.
        """
        pied = [math.hypot(*p[0]) for p in profils]
        tete = [p[-1][1] for p in profils]  # la tête est sur l'axe de la dent (x = 0)
        a = abs(aa)
        if interne:
            i = 0 if pied[0] > tete[0] else 1  # la couronne : son pied est plus grand que sa tête
            j = 1 - i
            if a + tete[j] > pied[i]:
                raise ErreurContact("entraxe trop grand (la tête du pignon sort du pied de la couronne)")
            if tete[i] - a < pied[j]:
                raise ErreurContact("entraxe trop grand (la tête de la couronne touche le pied du pignon)")
            if a <= tete[i] - tete[j]:
                raise ErreurContact("entraxe trop petit (les cercles de tête ne se coupent plus)")
            return
        if pied[0] + tete[1] > a or pied[1] + tete[0] > a:
            raise ErreurContact("entraxe trop petit (la tête d'une roue touche le pied de l'autre)")
        if abs((a * a + tete[0] ** 2 - tete[1] ** 2) / (2 * a * tete[0])) > 1:
            raise ErreurContact("entraxe trop grand (les cercles de tête ne se coupent plus)")

    def _reference(self):
        """Position de la menée au voisinage de laquelle on déroule toutes les autres.

        `premier_contact` rend une rotation dans [0, 2pi[ : si la position de repos de la menée est
        voisine de 0, elle bascule entre ~0 et ~2pi d'une position à l'autre et `debut` ne voit plus
        qu'un saut de 2pi au lieu du croisement des deux paires. C'est exactement le cas d'un
        engrènement intérieur, où la recherche se fait dans l'autre sens et où le contact tombe donc
        juste de l'autre côté de zéro. On déroule donc tout autour d'une valeur de référence prise
        une fois pour toutes ; l'écart physique entre paires voisines (un pas de la menée) reste
        très petit devant pi, donc l'ordre entre paires est conservé."""
        self.u0 = 0.0
        d = self._contact(0.0)[0]
        self.u0 = 0.0 if not math.isfinite(d) else d

    def _derouler(self, d):
        """`d` ramené dans l'intervalle de largeur 2pi centré sur la référence."""
        if not math.isfinite(d):
            return d
        return self.u0 + (d - self.u0 + math.pi) % (2 * math.pi) - math.pi

    def psi(self, phi, dos=False):
        """(rotation de la menée, point, tangente) pour la dent de la menante d'axe phi ; `dos` :
        avec les flancs arrière des dents (seulement en mode arcs)."""
        cache = self._cache_dos if (dos and self.arcs) else self._cache
        if phi not in cache:
            cache[phi] = self._contact(phi, dos and self.arcs)
        return cache[phi]

    def _contact(self, phi, dos=False):
        if self.arcs:
            fl1 = _flanc_place(self.fb1 if dos else self.fl1, phi + self.delta1, self.o1)
            d, q, g, courbures = premier_contact_arcs(
                fl1, self.ref2b if dos else self.ref2, self.o2, self.sens)
            if not dos:
                self._courbures[phi] = courbures
            return self._derouler(math.inf if d is None else d), q, g
        # recherche grossière, puis même recherche sur les seuls morceaux de flanc voisins
        # du contact, très finement échantillonnés : position et tangente du contact
        # bien plus précises qu'avec des cordes de 1e-3 (erreur en racine de la flèche)
        d0, q0, g0 = premier_contact(rotation(self.g1, phi + self.delta1, self.o1), self.ref2, self.o2, self.sens)
        if d0 is None:
            self._courbures[phi] = (math.nan, math.nan)
            return math.inf, None, None
        d0 = self._derouler(d0)
        f1 = rotation(self.f1, phi + self.delta1, self.o1)
        i1 = int(np.argmin(np.hypot(*(f1 - q0).T)))
        i2 = int(np.argmin(np.hypot(*(rotation(self.f2, self.sens * d0, self.o2) - q0).T)))
        s1, s2 = slice(max(0, i1 - FEN), i1 + FEN + 1), slice(max(0, i2 - FEN), i2 + FEN + 1)
        d, q, g = premier_contact(f1[s1], self.ref2f[s2], self.o2, self.sens)
        # `d` est brut, `d0` est déjà déroulé : comparer sans dérouler rejetterait toujours
        # l'affinage (2pi d'écart) et ne laisserait que la recherche grossière
        d = None if d is None else self._derouler(d)
        if d is None or abs(d - d0) > 1e-3:  # contact réel hors de la zone affinée
            d, q, g = d0, q0, g0
        # courbures au contact, dans le repère propre de chaque roue (la courbure ne dépend pas
        # de la rotation) : l'équivalent, sur des points, des rayons des arcs en contact
        j1 = int(np.argmin(np.hypot(*(f1 - q).T)))
        j2 = int(np.argmin(np.hypot(*(rotation(self.f2, self.sens * d, self.o2) - q).T)))
        self._courbures[phi] = (rayon_courbure(self.f1, j1), rayon_courbure(self.f2, j2))
        return self._derouler(d), q, g

    def theta_max(self, phi, dos=False):
        """Position la plus haute de la menée quand la menante est en phi (contact flanc menant ;
        `dos` : contact des flancs arrière, qui limite la course de la menée dans l'autre sens)."""
        return min(self._c(phi, k, dos) for k in (-1, 0, 1))

    def jeu(self, debut, n=10):
        """Jeu de la menée (degrés) sur une période : (min, max), comme AW(1) et AW(2), évalué
        en `n` positions de la menante régulièrement réparties.

        Par symétrie par rapport à la ligne des centres, la position la plus basse de la menée
        avec la menante en phi est l'image de la plus haute avec la menante en -phi ; l'écart
        entre les deux, replié modulo un pas de la menée, est le jeu.
        """
        jeux = []
        for j in range(n):
            phi = debut + j * self.a1 / (n - 1)
            somme = self.theta_max(phi) + self.theta_max(-phi, dos=True)
            jeux.append(math.degrees((somme + self.a2u / 2) % self.a2u - self.a2u / 2))
        return min(jeux), max(jeux)

    def _c(self, phi, k, dos=False):
        """Position de la menée imposée par la paire k (dent k de la menante)."""
        return self.psi(phi + k * self.a1, dos)[0] + k * self.a2u

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
        if ecart(grille[0]) <= self.tol:
            raise ErreurContact("début de menée introuvable (paire 0 déjà en contact)")
        for lo, hi in zip(grille, grille[1:]):
            if ecart(hi) <= self.tol:
                break
        else:
            raise ErreurContact("début de menée introuvable (pas de contact continu ?)")
        while hi - lo > 1e-7:
            m = (lo + hi) / 2
            lo, hi = (m, hi) if ecart(m) > self.tol else (lo, m)

        # les deux paires doivent vraiment se toucher (sinon rapport de conduite < 1)
        theta2 = min(self._c(hi, 0), self._c(hi, 1))
        limite = 5e-3 * self.roues[0].mo
        for k in (0, 1):
            q = self.psi(hi + k * self.a1)[1]
            place = self.sens * (theta2 - k * self.a2u)
            if q is None or _distance(q, rotation(self.g2, place, self.o2)) > limite:
                raise ErreurContact("deux paires de dents ne sont jamais en prise simultanément")
        return hi

    def simuler(self, mu=0.2, n_pos=N_POS, n_jeu=None):
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
        hertz = np.full(n, np.nan) if self.hertz else None
        for i in range(n):
            phi = debut + i * self.a1 / (n - 1)
            c0, q0, g0 = self.psi(phi)
            c1 = math.inf
            if not sortie:
                d1, q1, g1 = self.psi(phi + self.a1)
                c1 = d1 + self.a2u
                sortie = not math.isfinite(c1) or c1 > c0 + 1e-3
            paire = 1 if (not sortie and c1 < c0 - self.seuil_paire) else 0
            q, g = (q1, g1) if paire else (q0, g0)
            contact[i] = q
            if not sortie and abs(c1 - c0) <= self.seuil_contact2:  # les deux paires se touchent
                contact2[i] = q0 if paire else q1
            res["a1"][i], res["a2"][i], res["paire"][i] = phi, min(c0, c1), paire
            res["rap"][i], res["cou"][i], res["ren"][i], res["nor"][i] = efforts(self.o1, self.o2, q, g, mu)
            if hertz is not None:
                hertz[i] = self._hertz(res["nor"][i], self._courbures[phi + paire * self.a1])
        z1, z2 = self.roues[0].zz, self.roues[1].zz
        # jeu : toutes les positions avec des arcs (moteur rapide), 10 avec des points
        jeu_min, jeu_max = self.jeu(debut, n_jeu or (n if self.arcs else 10))
        # `res["a2"]` est compté dans le sens de recherche ; la simulation expose l'angle réel
        sim = Simulation(mu, res["a1"], self.sens * res["a2"], contact, contact2, res["paire"],
                         res["rap"], res["cou"], res["ren"], res["nor"], -z1 / z2, -z2 / z1,
                         jeu_min, jeu_max, hertz)
        _, i = sim.hertz_max()
        if i >= 0:  # où se trouve la pression maximale sur chacun des deux profils
            k = int(res["paire"][i])
            local1 = rotation([contact[i] - self.o1], -(res["a1"][i] + k * self.a1))[0]
            local2 = rotation([contact[i] - self.o2], -self.sens * (res["a2"][i] - k * self.a2u))[0]
            sim.hertz_local, sim.hertz_indice = (tuple(local1), tuple(local2)), i
        return sim

    def _hertz(self, force, rayons):
        """Pression de Hertz maximale (N/mm²) d'un contact linéaire, pour la force normale `force`
        (N) : sqrt(F * E* / (pi * H * (1 - nu^2) * R)), avec E* = E1*E2/(E1+E2), nu moyen, H la
        plus petite épaisseur, R = |R1*R2/(R1+R2)| (rayons signés : concave négatif), comme
        ENGANAC (KW(17))."""
        (e1, n1, h1), (e2, n2, h2) = self.hertz
        somme = rayons[0] + rayons[1]
        if not math.isfinite(force) or abs(somme) < 1e-12:
            return math.nan
        r = abs(rayons[0] * rayons[1] / somme)
        e = e1 * e2 / (e1 + e2)
        nu = (n1 + n2) / 2.0
        return math.sqrt(force * e / (math.pi * min(h1, h2) * (1.0 - nu * nu) * r))


def _distance(p, ligne):
    """Distance du point p à la ligne brisée `ligne`."""
    a, d = ligne[:-1], np.diff(ligne, axis=0)
    t = np.clip(((p - a) * d).sum(1) / np.maximum((d * d).sum(1), 1e-30), 0.0, 1.0)
    return float(np.min(np.hypot(*(a + t[:, None] * d - p).T)))


def inverser(roues, profils, mobiles, n=400):
    """Paire « sens inverse » : quand la menée devient menante, l'effort change de signe et le
    contact bascule sur les **flancs arrière**. Rend `(roues, profils, mobiles)` à passer tels
    quels à `Engrenement`.

    Deux opérations, et rien d'autre : on échange les deux roues (la menée devient la roue 1) et on
    met le demi-profil arrière de chacune en position de flanc actif. Le résultat est l'image
    miroir de la configuration réelle -- une réflexion ne change ni les vitesses, ni le rendement,
    ni la pression de Hertz, ni le jeu -- donc tout le moteur direct s'applique sans retouche.

    Vérifié sur une paire symétrique : à mu = 0 les deux sens donnent exactement 100 % de
    rendement, le rapport de transmission s'inverse exactement (-7.415725 contre -0.134767) et à
    mu = 0.2 le sens inverse donne 95.37 % contre 95.80 % à l'aller -- l'écart attendu, approche et
    retraite échangeant leurs rôles.

    Le jeu est le **même jeu physique** des deux côtés, mais `Engrenement.jeu` le rend en degrés de
    la roue menée, et la menée n'est pas la même : les deux nombres sont dans le rapport des
    nombres de dents. C'est aussi pourquoi épaissir le flanc arrière pour le raidir mange le jeu,
    et qu'il faut garder la contrainte de jeu minimal pendant l'optimisation du sens inverse.
    """
    mi = []
    for m in mobiles:
        # les deux flancs **échangent leurs rôles**, l'arrière devient l'actif et l'avant le dos.
        # Collapser le dos (points_dos=None) serait plus court mais ferait calculer le jeu entre un
        # flanc et lui-même : `Engrenement.jeu` a besoin des deux flancs de la même dent.
        mi.append(replace(m, points=np.asarray(m.profil_dos(), dtype=float),
                          points_dos=np.asarray(m.points, dtype=float)))
    pr = [np.asarray(m.echantillonner(n), dtype=float) for m in mi]
    return [roues[1], roues[0]], [pr[1], pr[0]], [mi[1], mi[0]]


def simuler(roues, profils, mu=0.2, n_pos=N_POS, da=0.0):
    return Engrenement(roues, profils, da).simuler(mu, n_pos)
