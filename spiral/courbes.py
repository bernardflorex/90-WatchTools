# courbes.py
"""
Courbe de raccordement générique (Grossmann intérieure, côté collet, OU
Phillips extérieure, côté piton), représentée par une courbe de Bézier sur le
rayon r(u), u in [0,1], à nombre de points de contrôle configurable.

Rappel de théorie (Phillips, 1860 -- voir recherche du 2026-09-22) : la
condition de concentricité ne prescrit PAS de forme précise pour la courbe
terminale, seulement une condition sur la position du centre de gravité de
cette courbe (voir `condition_phillips` ci-dessous : d = R^2/L). Les formes
classiques (arc unique, double arc, double arc + segment droit) sont des
simplifications héritées du pliage MANUEL d'un ressort métallique. En
silicium usiné par DRIE cette contrainte de fabrication disparaît : on peut
laisser l'optimiseur explorer des formes à plusieurs points de contrôle
(inflexions possibles), sous les contraintes physiques réelles (jeu
inter-spires, cf. collision.py) et l'objectif réel (erreur de marche, cf.
isochronisme.py). La condition de Phillips reste calculée à titre de repère
de plausibilité, pas comme loi de construction.

Chaque zone a toujours :
  - une extrémité ANCRÉE, raccordée en C1 (position ET pente imposées) au
    pas constant du corps ;
  - une extrémité LIBRE, à rayon imposé (contrainte dure : r_min au collet,
    r_max au piton) mais pente et points intermédiaires libres -- c'est là
    que la forme se dessine.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import comb


def _bezier(u: np.ndarray, points) -> np.ndarray:
    n = len(points) - 1
    r = np.zeros_like(u, dtype=float)
    for i, p_i in enumerate(points):
        r = r + comb(n, i) * (1 - u) ** (n - i) * u ** i * p_i
    return r


@dataclass
class ZoneRaccordement:
    """Une zone de raccordement (Grossmann ou Phillips), avec ses points de
    contrôle de Bézier déjà résolus -- utile pour diagnostics (monotonie) et
    pour `condition_phillips`."""
    points_controle: np.ndarray
    ancre_a_u1: bool

    @property
    def n_points(self) -> int:
        return len(self.points_controle)

    def r(self, u: np.ndarray) -> np.ndarray:
        return _bezier(u, self.points_controle)

    def monotone(self, tol: float = 1e-9) -> bool:
        """True si les points de contrôle sont monotones non-décroissants
        (condition SUFFISANTE, via les propriétés des polynômes de
        Bernstein, pour que r(u) soit elle-même monotone -- un spiral qui
        ne l'est pas replie physiquement le ruban sur lui-même)."""
        return bool(np.all(np.diff(self.points_controle) >= -tol))


def construire_zone(n_points: int, r_ancre: float, pente_ancre_du: float,
                     r_libre: float, offsets_libres: np.ndarray,
                     ancre_a_u1: bool) -> ZoneRaccordement:
    """Construit les points de contrôle de la zone.

    n_points : nombre TOTAL de points de contrôle (>= 4). Degré du Bézier
        = n_points - 1. n_points - 3 d'entre eux sont libres (variables
        d'optimisation, voir `offsets_libres`) ; les 3 autres sont fixés par
        les contraintes dures (position à chaque bout) et par la continuité
        C1 (pente imposée à l'extrémité ancrée).
    r_ancre, pente_ancre_du : rayon et pente dr/du (PAS dr/dtheta -- la
        conversion par la largeur angulaire de la zone se fait chez
        l'appelant, voir parametrage.py) imposés à l'extrémité ancrée.
    r_libre : rayon imposé à l'extrémité libre (r_min ou r_max).
    offsets_libres : écarts (mm) des points de contrôle intermédiaires par
        rapport à une référence rectiligne entre les deux points fixés
        adjacents -- ce sont les VRAIES variables libres de l'optimisation
        pour cette zone (0 partout = raccord parfaitement rectiligne, la
        forme la plus neutre possible). len(offsets_libres) = n_points - 3.
    ancre_a_u1 : True si l'extrémité ancrée est en u=1 (courbe Grossmann,
        côté collet, l'extrémité libre r_min étant en u=0) ; False si
        l'extrémité ancrée est en u=0 (courbe Phillips, côté piton,
        l'extrémité libre r_max étant en u=1).
    """
    n_libres = n_points - 3
    if n_libres != len(offsets_libres):
        raise ValueError(
            f"n_points={n_points} implique {n_libres} points libres, "
            f"reçu {len(offsets_libres)} offsets."
        )
    degre = n_points - 1

    if ancre_a_u1:
        p0 = r_libre
        p_nm1 = r_ancre - pente_ancre_du / degre
        t = np.linspace(0.0, 1.0, n_libres + 2)[1:-1]
        ref = p0 + t * (p_nm1 - p0)
        points = [p0, *(ref + np.asarray(offsets_libres)), p_nm1, r_ancre]
    else:
        p1 = r_ancre + pente_ancre_du / degre
        t = np.linspace(0.0, 1.0, n_libres + 2)[1:-1]
        ref = p1 + t * (r_libre - p1)
        points = [r_ancre, p1, *(ref + np.asarray(offsets_libres)), r_libre]

    return ZoneRaccordement(points_controle=np.array(points, dtype=float), ancre_a_u1=ancre_a_u1)


def condition_phillips(x_zone: np.ndarray, y_zone: np.ndarray, ds_zone: np.ndarray,
                        x_centre_spiral: float, y_centre_spiral: float,
                        x_decrochement: float, y_decrochement: float) -> float:
    """Repère de plausibilité (PAS une contrainte imposée) : la condition
    classique de Phillips (1860) veut que le CG de la courbe terminale soit
    à une distance d = R^2/L du centre du spiral, sur le rayon perpendiculaire
    à la droite (centre spiral -> point de décrochement), R étant le rayon au
    point de décrochement et L la longueur de la courbe. Retourne le rapport
    d_mesure / d_attendue (proche de 1.0 pour une courbe classiquement
    "réglante" -- un écart n'est pas une erreur ici, puisque le vrai critère
    de validation du projet est l'orbite du CG calculée par l'elastica
    complet sur amplitude/position réelles, cf. isochronisme.py)."""
    L = float(np.sum(ds_zone))
    if L < 1e-9:
        return float("nan")
    R = np.hypot(x_decrochement - x_centre_spiral, y_decrochement - y_centre_spiral)
    x_cg = np.sum(x_zone * ds_zone) / L
    y_cg = np.sum(y_zone * ds_zone) / L
    d_mesure = np.hypot(x_cg - x_centre_spiral, y_cg - y_centre_spiral)
    d_attendue = (R ** 2) / L
    if d_attendue < 1e-9:
        return float("nan")
    return float(d_mesure / d_attendue)


def construire_zone_s(n_points: int, r_debut: float, pente_debut_du: float,
                       r_fin: float, pente_fin_du: float,
                       offsets_libres: np.ndarray) -> ZoneRaccordement:
    """Courbe de raccordement dont les DEUX extrémités sont ancrées en position
    ET en pente -- par opposition à `construire_zone`, où l'extrémité libre n'a
    que son rayon imposé.

    C'est ce qu'exige la transition vers un arc concentrique de raquette : il
    faut à la fois rejoindre la spirale d'Archimède avec sa pente (continuité
    C1 côté corps) et arriver sur l'arc avec une pente NULLE (l'arc est à rayon
    constant), sans quoi le spiral présenterait un coude à l'entrée de la zone
    balayée par les goupilles.

    Quatre conditions déterminent exactement une cubique : avec n_points=4 la
    courbe est entièrement fixée (interpolation d'Hermite cubique), et il n'y a
    aucune variable d'optimisation. Au-delà, les n_points-4 points de contrôle
    intermédiaires redeviennent libres, sans toucher aux tangences.

    La forme en S (inflexion) n'est pas imposée ici : elle apparaît d'elle-même
    dès que la pente moyenne exigée -- (r_fin - r_debut) / étendue angulaire --
    dépasse la pente d'entrée, ce qui est le cas courant (mesuré sur un calibre
    de référence : 117 à 700 µm/rad exigés contre 25 µm/rad de pas d'Archimède,
    selon l'angle de transition retenu).
    """
    n_libres = n_points - 4
    if n_points < 4:
        raise ValueError(f"n_points doit valoir au moins 4 (reçu {n_points})")
    if n_libres != len(offsets_libres):
        raise ValueError(
            f"n_points={n_points} implique {n_libres} points libres, "
            f"reçu {len(offsets_libres)} offsets."
        )
    degre = n_points - 1

    p0 = r_debut
    p1 = r_debut + pente_debut_du / degre
    pn = r_fin
    pn_1 = r_fin - pente_fin_du / degre

    if n_libres > 0:
        t = np.linspace(0.0, 1.0, n_libres + 2)[1:-1]
        ref = p1 + t * (pn_1 - p1)
        interieurs = list(ref + np.asarray(offsets_libres))
    else:
        interieurs = []

    points = [p0, p1, *interieurs, pn_1, pn]
    return ZoneRaccordement(points_controle=np.array(points, dtype=float), ancre_a_u1=False)
