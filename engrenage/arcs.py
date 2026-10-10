"""Approximation d'un profil par des arcs de cercle tangents, et export .Eno.

Repris de Krebax.F90 et Lagrmo.F90 (K. Trümpy) : le demi-pas du profil (points de
evolvente.evolve, `ZP` du Fortran) est remplacé par une suite d'arcs de cercle qui se
raccordent tangentiellement. Les arcs sont définis par leurs points de transition ; les
rayons et les centres s'en déduisent (`krbgs`) à partir de l'angle de la tangente au
départ (-pi/Z, dans le repère du profil). Les points de transition sont optimisés
(multiplicateur de Lagrange) pour minimiser la somme des carrés des distances du profil aux
arcs, avec une tangente horizontale (angle 0) à la fin, sur l'axe de la dent.
"""

import json
import math
from dataclasses import dataclass, replace
from pathlib import Path

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
    ecart_flanc: float = None  # la même, sur le flanc actif seul (None : pas de profil de référence)

    @property
    def n_points(self):
        return len(self.points)


def _ecart_angle(a, b):
    """a - b ramené dans (-pi, pi] : la cible peut être pi (couronne), où `krbgs` renvoie
    indifféremment +pi ou -pi, et deux évaluations voisines peuvent tomber de part et d'autre."""
    return (np.asarray(a) - b + math.pi) % (2 * math.pi) - math.pi


def angle_depart(zz):
    """Angle de la tangente au départ du demi-profil, dans la convention de `krbgs` (qui la prend
    dans le sens opposé au parcours) : `-pi/|Z|`, aussi bien pour une denture extérieure que pour
    une couronne intérieure -- la tangente de parcours vaut `pi - pi/|Z|` dans les deux cas."""
    return -math.pi / abs(zz)


def angle_fin(zz):
    """Angle visé par la tangente à l'arrivée, sur l'axe de la dent : 0 pour une denture
    extérieure, pi pour une couronne (le profil y arrive en descendant vers le centre)."""
    return 0.0 if zz > 0 else math.pi


SAGITTA_DROITE = 5e-5  # mm : en dessous, un arc est une droite et son signe ne veut rien dire


def inversions(kr, ra, sagitta=SAGITTA_DROITE):
    """Nombre de changements de signe des rayons, en ignorant les arcs rectilignes en pratique.

    Un raccord de pied droit (`qg > 0`) sort de `krbgs` avec un rayon énorme dont le signe ne
    tient qu'au bruit : le compter fait apparaître deux inversions fantômes et fait rejeter des
    profils parfaits (mesuré sur une couronne Z=-150 : le meilleur ajustement, à 0.0009 µm,
    était écarté au profit d'un autre à 1.6 µm). Un arc dont la flèche est sous `sagitta` est
    donc sauté.
    """
    ra = np.asarray(ra, dtype=float)
    kr = np.asarray(kr, dtype=float)
    demi = np.hypot(*(kr[1:] - kr[:-1]).T) / 2.0
    r = np.abs(ra)
    fleche = r - np.sqrt(np.maximum(r * r - demi * demi, 0.0))
    s = np.sign(ra[fleche >= sagitta])
    return int(np.count_nonzero(np.diff(s))) if len(s) > 1 else 0


def convexite(zz, ra):
    """Rayons signés ramenés à la convexité *physique* (positif = matière du côté du centre de
    l'arc). `krbgs` les donne relativement au sens de parcours ; pour une couronne, dont la
    matière est à l'extérieur, c'est l'opposé. Multiplier par un signe ne change pas le nombre
    d'inversions, donc le contrôle SWTE (une seule inversion) reste valable tel quel -- et un
    flanc actif de couronne, concave, y est légitime."""
    return ra if zz > 0 else -ra


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


def _krbgs_lot(a1, kr):
    """`krbgs` pour un lot de chaînes : `kr` (lot, points, 2) -> centres (lot, arcs, 2), rayons
    signés (lot, arcs), angle final (lot,). Même récurrence, vectorisée sur le lot : LAGRMO a
    besoin d'une chaîne par coordonnée perturbée et les évaluer ensemble est ce qui rend
    l'ajustement rapide."""
    co = np.full(len(kr), math.cos(a1))
    si = np.full(len(kr), math.sin(a1))
    n = kr.shape[1] - 1
    pm, ra = np.empty((len(kr), n, 2)), np.empty((len(kr), n))
    for i in range(n):
        dx, dy = kr[:, i + 1, 0] - kr[:, i, 0], kr[:, i + 1, 1] - kr[:, i, 1]
        rr = 0.5 * (dx * dx + dy * dy) / (dx * si - dy * co)
        pm[:, i, 0], pm[:, i, 1] = kr[:, i, 0] + rr * si, kr[:, i, 1] - rr * co
        vx, vy = kr[:, i + 1, 0] - pm[:, i, 0], kr[:, i + 1, 1] - pm[:, i, 1]
        ra[:, i] = np.copysign(rr, vx * dy - vy * dx)
        co, si = vy / np.abs(rr), -vx / np.abs(rr)
    return pm, ra, np.arctan2(-si, co)


def _dehors_ecart(xy, kr, pm, ra):
    """(point hors du secteur de l'arc, distance du point à l'arc) pour un lot de chaînes :
    tableaux (lot, arcs, points)."""
    a = xy[None, None, :, :] - pm[:, :, None, :]  # (lot, arcs, points, 2)
    c1 = kr[:, :-1, None, :] - pm[:, :, None, :]
    c2 = kr[:, 1:, None, :] - pm[:, :, None, :]

    def signe(c):
        """signe de l'angle (rayon de l'extrémité -> rayon du point) : celui du produit
        vectoriel, sauf sur l'alignement exact où l'angle vaut 0 ou pi (atan2 est inutile ici,
        et c'est le point chaud de tout l'ajustement)"""
        cr = a[..., 0] * c[..., 1] - a[..., 1] * c[..., 0]
        do = a[..., 0] * c[..., 0] + a[..., 1] * c[..., 1]
        return np.where(cr != 0.0, np.sign(cr), np.where(do < 0.0, 1.0, 0.0))

    dehors = signe(c1) * signe(c2) > 0
    return dehors, np.abs(np.abs(ra)[:, :, None] - np.hypot(a[..., 0], a[..., 1]))


def _affectation(dehors):
    """Arc affecté à chaque point, pour un lot : on passe à l'arc suivant dès qu'un point sort
    de l'arc courant et on reste sur le dernier arc ensuite (SQAB du Fortran). Le balayage se
    fait arc par arc, pas point par point."""
    nl, n_arcs, n = dehors.shape
    rang = np.arange(n)
    affect = np.zeros((nl, n), dtype=np.intp)
    depart = np.zeros(nl, dtype=np.intp)
    courant = np.zeros(nl, dtype=np.intp)
    for j in range(n_arcs):
        actif = courant == j
        if not actif.any():
            break
        # le point qui fait passer à l'arc suivant appartient à cet arc sans être retesté
        # (SQAB teste puis affecte) ; seul l'arc 0 est testé dès son premier point
        seuil = depart if j == 0 else depart + 1
        d = dehors[:, j, :] & (rang[None, :] >= seuil[:, None]) & actif[:, None]
        trouve = d.any(1)
        fin = np.where(trouve, d.argmax(1), n)
        if j == n_arcs - 1:  # plus d'arc suivant : tous les points restants y restent
            fin = np.full(nl, n)
            trouve = np.zeros(nl, dtype=bool)
        affect[actif[:, None] & (rang[None, :] >= depart[:, None])
               & (rang[None, :] < fin[:, None])] = j
        depart = np.where(actif, fin, depart)
        courant = np.where(actif & trouve, j + 1, courant)
    return affect


def _distances(xy, kr, pm, ra):
    """(écart de chaque point à chaque arc, arc affecté à chaque point) pour une seule chaîne."""
    dehors, ecart = _dehors_ecart(xy, kr[None], pm[None], ra[None])
    return ecart[0], _affectation(dehors)[0]


def sqab(xy, kr, pm, ra):
    """Somme des carrés des distances des points `xy` à leur arc, distance maximale, indice
    (SQAB du Fortran)."""
    ecart, affect = _distances(xy, kr, pm, ra)
    d = ecart[affect, np.arange(len(xy))]
    imax = int(np.argmax(d))
    return float(d @ d), float(d[imax]), imax


def _sqab_lot(xy, kr, pm, ra):
    """Somme des carrés des distances, pour un lot de chaînes : (lot,)."""
    dehors, ecart = _dehors_ecart(xy, kr, pm, ra)
    affect = _affectation(dehors)
    d = np.take_along_axis(ecart, affect[:, None, :], axis=1)[:, 0, :]
    return (d * d).sum(1)


GAIN_MINI = 0.002  # gain relatif d'écart maximal sous lequel une itération de LAGRMO ne sert plus


N_ITER_LAGRMO = 100  # garde-fou : le Fortran n'itère pas non plus indéfiniment


def lagrmo(aw, xy, kp):
    """Optimise les points de transition `kp` (les deux extrémités restent fixes) : LAGRMO.

    Le Fortran enchaîne les itérations jusqu'à ce que l'écart maximal remonte ; on s'arrête aussi
    dès que le gain relatif tombe sous `GAIN_MINI`, car la fin de la descente gagne des fractions
    de nanomètre en consommant l'essentiel du temps de calcul.
    """
    nk = len(kp)
    dt = 1e-4  # variation relative des coordonnées de transition
    rx = nk - 2  # relaxation
    mn = 2 * nk - 3
    kr = kp.copy()
    dt *= math.hypot(kr[-1, 0] - kr[0, 0], kr[-1, 1] - kr[0, 1])
    meilleur = (math.inf, 0, 0.0)  # (ecart max, indice, erreur d'angle)
    # inconnues : x, y des points de transition intérieurs ; la dernière ligne/colonne porte le
    # multiplicateur de Lagrange (tangente finale imposée)
    lig = np.arange(mn - 1)
    pt, coord = 1 + lig // 2, lig % 2
    for _ in range(N_ITER_LAGRMO):
        pm, ra, a0 = krbgs(aw[0], kr)
        f0 = sqab(xy, kr, pm, ra)[0]
        # une chaîne perturbée par (inconnue, sens) : évaluées toutes ensemble
        lot = np.repeat(kr[None], 2 * (mn - 1), axis=0)
        lot[2 * lig, pt, coord] -= dt
        lot[2 * lig + 1, pt, coord] += dt
        pm_l, ra_l, a_l = _krbgs_lot(aw[0], lot)
        f_l = _sqab_lot(xy, lot, pm_l, ra_l)
        fm, fp, am, ap = f_l[0::2], f_l[1::2], a_l[0::2], a_l[1::2]
        mt, bs = np.zeros((mn, mn)), np.zeros(mn)
        bs[:mn - 1] = (fm - fp) / (2 * dt)
        mt[lig, lig] = rx * (fm + fp - 2 * f0) / (dt * dt)
        mt[lig, mn - 1] = mt[mn - 1, lig] = _ecart_angle(am, ap) / (2 * dt)
        bs[mn - 1] = _ecart_angle(a0, aw[1])
        try:
            dx = np.linalg.solve(mt, bs)
        except np.linalg.LinAlgError:
            dx = None
        if dx is None:  # matrice singulière
            break
        kr[pt, coord] += dx[:mn - 1]  # (point, coordonnée) distincts : pas de cumul à gérer
        pm, ra, a0 = krbgs(aw[0], kr)
        _, mr1, mr2 = sqab(xy, kr, pm, ra)
        if mr1 > meilleur[0]:  # l'écart maximal remonte : on garde le meilleur
            break
        gain = meilleur[0] - mr1
        meilleur = (mr1, mr2, -_ecart_angle(a0, aw[1]))
        kp[:] = kr
        if gain < GAIN_MINI * mr1:  # la descente ne rapporte plus rien
            break
    return meilleur


def _mesures(xy):
    """(abscisse curviligne, rotation cumulée de la tangente) aux points du profil : deux façons
    de « parcourir » le profil indépendamment de la densité de son échantillonnage."""
    d = np.diff(xy, axis=0)
    s = np.concatenate(([0.0], np.cumsum(np.hypot(*d.T))))
    a = np.arctan2(d[:, 1], d[:, 0])
    da = np.abs((np.diff(a) + math.pi) % (2 * math.pi) - math.pi)
    return s, np.concatenate(([0.0, 0.0], np.cumsum(da)))


def _bornes(mesure, n_arcs):
    """Indices intérieurs découpant le profil en `n_arcs` morceaux de `mesure` égale."""
    cible = np.linspace(0.0, mesure[-1], n_arcs + 1)[1:-1]
    return np.clip(np.searchsorted(mesure, cible), 1, len(mesure) - 2)


def rayon_flanc(points, zz):
    """Rayon où finit le congé de pied d'une chaîne d'arcs : fin du dernier arc concave (rayon
    signé négatif). Au-dessous, le profil n'est pas une surface de travail -- rien n'y touche, et
    une fraise de taillage ne l'y reproduit pas (elle y creuse sa propre trochoïde). 0 s'il n'y a
    pas d'arc concave."""
    pts = np.asarray(points, dtype=float)
    ra = convexite(zz, krbgs(angle_depart(zz), pts)[1])
    negatifs = np.flatnonzero(ra < 0.0)
    return 0.0 if len(negatifs) == 0 else float(np.hypot(*pts[negatifs[-1] + 1]))


def _ajuster_points(xy, aw, kr, utile=None):
    """Ajuste (LAGRMO) la chaîne d'arcs passant par les points de transition `kr` (les deux
    extrémités restent fixes) et retourne (clé de comparaison, ProfilArcs, écart max par arc).
    `aw` : (tangente au départ, tangente visée à l'arrivée)."""
    a1 = aw[0]
    ecart, indice, erreur = lagrmo(aw, xy, kr)
    pm, ra, _ = krbgs(a1, kr)
    ec, affect = _distances(xy, kr, pm, ra)
    d = ec[affect, np.arange(len(xy))]
    par_arc = np.zeros(len(pm))
    np.maximum.at(par_arc, affect, d)
    flanc = float(d[utile].max()) if utile is not None else None
    return (max(0, inversions(kr, ra) - 1), ecart), \
        ProfilArcs(kr, pm, ra, ecart, indice, erreur, flanc), par_arc


def _ajuster(xy, aw, idx, utile=None):
    """Ajuste la chaîne d'arcs découpée aux indices intérieurs `idx` du profil.

    La clé classe d'abord par nombre d'inversions de signe des rayons au-delà d'une seule : un
    profil qui en a plusieurs est refusé plus loin par le contrôle SWTE (`optim.profil_valide`)
    et ne sert donc à rien, même avec un écart plus faible. `utile` (masque des points du flanc
    actif) ne sert qu'à renseigner `ecart_flanc` : classer là-dessus a été essayé et dégrade le
    résultat quand le nombre d'arcs demandé n'est pas celui de la référence.
    """
    return _ajuster_points(xy, aw, xy[np.concatenate(([0], idx, [len(xy) - 1]))].copy(), utile)


def _redistribuer(par_arc, s, idx, n_max):
    """Nouveau découpage égalisant l'écart maximal par arc : équirépartition en `écart^(1/3)`,
    l'exposant habituel pour une approximation d'ordre 3 (une borne se déplace vers l'arc qui
    porte le plus d'erreur). None si le découpage n'est pas exploitable."""
    poids = np.maximum(par_arc, 1e-12) ** (1.0 / 3.0)
    cum = np.concatenate(([0.0], np.cumsum(poids)))
    bornes = np.concatenate(([0.0], s[idx], [s[-1]]))
    cible = np.interp(np.linspace(0.0, cum[-1], len(idx) + 2)[1:-1], cum, bornes)
    neuf = np.clip(np.searchsorted(s, cible), 1, n_max)
    if len(np.unique(neuf)) != len(idx) or np.array_equal(neuf, idx):
        return None
    return neuf


N_REDISTRIB = 6  # redécoupages essayés par découpage de départ


def approximer(profil, zz, qg, narcs, reference=None):
    """Profil (liste de points du demi-pas) remplacé par des arcs de cercle tangents.

    `narcs` arcs sur la développante, plus l'arc de pied, l'arc de tête et, si `qg > 0`, le
    raccord de pied (comme KREBAX : NE(K)+3 arcs, +1 avec raccord).

    `reference` : points de transition du profil dont ce profil-ci est issu, quand on le sait --
    typiquement la roue d'origine d'une roue taillée par `fraise.tailler`. S'ils décrivent le
    même nombre d'arcs, ils fournissent un découpage de départ de plus, et c'est le meilleur qui
    soit puisque c'est la structure de la forme cherchée : sur `validation.json`, l'erreur dans
    la zone de contact de la roue de 12 dents tombe de 1.43 à 0.16 µm et la variation de couple
    après taillage de 3.2 % à 1.2 % (celle du profil d'origine est 1.35 %). Ils servent aussi à
    renseigner `ecart_flanc`, l'écart au-dessus du congé de pied (`rayon_flanc`) : le congé,
    creusé en trochoïde par la fraise, s'écarte légitimement de plusieurs microns et ne touche
    rien, donc l'écart maximal seul donne une fausse idée de la qualité de l'ajustement.

    Le découpage initial de KREBAX repère les arcs par des *indices* de points (les NI//10
    premiers pour le pied, les NI//10 derniers pour la tête), ce qui n'est juste que pour
    l'échantillonnage de `evolvente.evolve`. Sur un profil échantillonné autrement -- une roue
    taillée par `fraise.tailler`, un profil en arcs rééchantillonné -- ces indices tombent
    n'importe où, `lagrmo` ne déplace les points que par petits pas et l'ajustement reste bloqué
    à plusieurs microns avec des rayons qui changent de signe à chaque arc (profil refusé par
    SWTE). On part donc aussi de deux découpages indépendants de l'échantillonnage (abscisse
    curviligne et rotation de la tangente égales) et, depuis chacun, on enchaîne des
    redécoupages qui égalisent l'erreur par arc, en gardant le meilleur ajustement obtenu.
    """
    if not 1 <= narcs <= NARCS_MAX:
        raise ValueError(f"nombre d'arcs entre 1 et {NARCS_MAX} attendu")
    xy = np.asarray(profil, dtype=float)
    ni = len(xy)
    m = 1 if qg > 0 else 0  # raccord rectiligne de pied
    n_arcs = narcs + 2 + m
    aw = (angle_depart(zz), angle_fin(zz))  # tangentes au départ et sur l'axe de la dent
    zz = abs(zz)  # à partir d'ici seul le pas compte, le sens est passé par `aw`
    s, theta = _mesures(xy)
    nk2 = ni // 10
    nk1 = nk2 - 2 * m
    kp = [1, nk1 + 1] + ([nk2 + 1] if m else [])  # indices (base 1) des points de transition
    kp += [nk2 + j * (ni - 2 * nk2) // narcs for j in range(1, narcs + 1)] + [ni]
    departs = [np.array([k - 1 for k in kp[1:-1]], dtype=np.intp)]  # KREBAX
    departs += [_bornes(mesure, n_arcs) for mesure in (s, theta) if mesure[-1] > 0.0]
    utile = None
    if reference is not None:
        ref = np.asarray(reference, dtype=float)
        r = np.hypot(*xy.T) >= rayon_flanc(ref, zz)
        if r.sum() >= 10:  # sinon le congé occupe tout : écart du flanc non renseigné
            utile = r
        if len(ref) == n_arcs + 1:  # même nombre d'arcs : on repart de ses points de transition
            # les points eux-mêmes, pas le point échantillonné le plus proche : un arc court
            # (une tête, un congé) tombe sinon entre deux indices voisins, l'arc n'est plus
            # contraint que par un point et LAGRMO lui donne n'importe quelle courbure
            cle, profil_arcs, _ = _ajuster_points(
                xy, aw, np.vstack((xy[0], ref[1:-1], xy[-1])).copy(), utile)
            # on garde cet ajustement tel quel s'il est exploitable, sans le mettre en
            # concurrence : aucun écart scalaire ne permet de choisir à sa place (essayé, sur
            # l'écart maximal comme sur celui du flanc actif ; les deux retiennent
            # régulièrement une chaîne voisine tout aussi « précise » mais qui double ou triple
            # la variation de couple)
            if cle[0] == 0:
                return profil_arcs
    departs = [d for d in departs if len(d) == n_arcs - 1 and len(np.unique(d)) == n_arcs - 1]
    if not departs:
        raise ValueError("profil trop court pour ce nombre d'arcs")
    meilleur = None
    vus = set()
    for idx in departs:
        for _ in range(N_REDISTRIB):
            if tuple(idx) in vus:
                break
            vus.add(tuple(idx))
            # on poursuit la série même quand un redécoupage est moins bon : la suite du chemin
            # retombe souvent bien plus bas (s'arrêter à la première remontée coûte, mesuré sur
            # validation.json, 5 µm au lieu de 1 µm)
            cle, profil_arcs, par_arc = _ajuster(xy, aw, idx, utile)
            if meilleur is None or cle < meilleur[0]:
                meilleur = (cle, profil_arcs)
            suivant = _redistribuer(par_arc, s, idx, ni - 2)
            if suivant is None:
                break
            idx = suivant
    return meilleur[1]


def normales(points, centres):
    """Normales unitaires (du centre de l'arc vers le point) aux points de transition ; à un
    raccord tangentiel, les deux arcs voisins ont la même normale."""
    pts = np.asarray(points, dtype=float)
    idx = [min(i, len(pts) - 2) for i in range(len(pts))]
    v = pts - np.asarray(centres)[idx]
    return v / np.hypot(*v.T)[:, None]


def echantillonner(points, centres, rayons, n_total=400):
    """Points (n,2) répartis sur les arcs (à peu près proportionnellement à leur longueur),
    du début à la fin du profil."""
    n_arcs = len(rayons)
    arcs_ = []
    for i in range(n_arcs):
        c, p, q = centres[i], points[i], points[i + 1]
        a0 = math.atan2(p[1] - c[1], p[0] - c[0])
        a1 = math.atan2(q[1] - c[1], q[0] - c[0])
        da = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi  # angle balayé, le plus court
        arcs_.append((c, a0, da, abs(rayons[i])))
    longueurs = [r * abs(da) for _, _, da, r in arcs_]
    total = sum(longueurs) or 1.0
    pts = []
    for i, (c, a0, da, r) in enumerate(arcs_):
        n = max(8, round(n_total * longueurs[i] / total))
        t = np.linspace(0.0, 1.0, n, endpoint=(i == n_arcs - 1))
        pts.append(np.column_stack((c[0] + r * np.cos(a0 + t * da), c[1] + r * np.sin(a0 + t * da))))
    return np.vstack(pts)


def _echantillons(self, n_total=400):
    return echantillonner(self.points, self.centres, self.rayons, n_total)


ProfilArcs.echantillonner = _echantillons
ProfilArcs.normales = lambda self: normales(self.points, self.centres)


def _milieu_arc(points, centres, rayons, i):
    """Point (x, y) au milieu de l'arc `i` (même paramétrage que `echantillonner`) : krbgs
    reconstruit ensuite ce même arc en deux moitiés sur le même cercle, sans toucher aux autres
    (vérifié à 1e-12 mm), donc l'insérer ne change pas la forme du profil."""
    c, r = centres[i], abs(rayons[i])
    p, q = points[i], points[i + 1]
    a0 = math.atan2(p[1] - c[1], p[0] - c[0])
    a1 = math.atan2(q[1] - c[1], q[0] - c[0])
    da = (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi
    return c[0] + r * math.cos(a0 + da / 2.0), c[1] + r * math.sin(a0 + da / 2.0)


def ajouter_point(mobile, i):
    """Nouveau `Mobile` avec un point ajouté au milieu de l'arc `i` (0-based parmi les arcs
    actuels du demi-profil menant) : la forme du profil ne change pas (un arc de plus, sur le
    même cercle de part et d'autre du nouveau point), c'est donc toujours possible sans risque
    ; utile pour redonner une variable à l'optimisation sans repartir de zéro. Si la dent est
    dissymétrique, le même point (au même endroit relatif) est ajouté au demi-profil arrière,
    pour que les deux gardent le même nombre d'arcs (nécessaire à `optim.Probleme`)."""
    pts = np.asarray(mobile.points, dtype=float)
    if not 0 <= i < len(pts) - 1:
        raise ValueError(f"arc {i + 1} inexistant (1 à {len(pts) - 1})")
    pm, ra = mobile.arcs()
    p = _milieu_arc(pts, pm, ra, i)
    npts = np.insert(pts, i + 1, p, axis=0)
    dr = list(mobile.drapeaux[:i + 1]) + [1] + list(mobile.drapeaux[i + 1:])
    dos = None
    if mobile.points_dos is not None:
        pd = np.asarray(mobile.points_dos, dtype=float)
        pmd, rad = mobile.arcs_dos()
        dos = np.insert(pd, i + 1, _milieu_arc(pd, pmd, rad, i), axis=0)
    return replace(mobile, points=npts, drapeaux=dr, points_dos=dos)


def retirer_point(mobile, i):
    """Nouveau `Mobile` sans le point de transition qui commence l'arc `i` (fusionne cet arc
    avec le précédent ; pour l'arc 0, qui n'en a pas, on retire son point de fin à la place,
    fusionné avec le suivant). Contrairement à `ajouter_point`, la forme du profil change après
    ce point (l'arc fusionné ne passe généralement plus par lui) : voir `ecart_um` pour la
    mesurer avant de confirmer. Le même point est retiré du demi-profil arrière s'il existe."""
    n = len(mobile.points)
    idx = i if i >= 1 else 1
    if not 1 <= idx <= n - 2:
        raise ValueError("aucun point à retirer pour cet arc (il ne reste que les deux extrémités)")
    pts = np.delete(np.asarray(mobile.points, dtype=float), idx, axis=0)
    dr = [f for j, f in enumerate(mobile.drapeaux) if j != idx]
    dos = None
    if mobile.points_dos is not None:
        dos = np.delete(np.asarray(mobile.points_dos, dtype=float), idx, axis=0)
    return replace(mobile, points=pts, drapeaux=dr, points_dos=dos)


def ecart_um(mobile_a, mobile_b, n=400):
    """Écart maximal (µm) entre les demi-profils menants échantillonnés de deux mobiles de la
    même roue : plus grande distance d'un point de l'un au profil de l'autre, dans les deux sens."""
    a, b = mobile_a.echantillonner(n), mobile_b.echantillonner(n)
    d = np.hypot(a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1])
    return 1e3 * max(d.min(axis=1).max(), d.min(axis=0).max())


@dataclass
class Mobile:
    """Un mobile dont le profil est donné par ses points de transition d'arcs de cercle.

    Le profil est celui d'un demi-pas (repère de evolvente.evolve : départ sur l'axe du creux,
    fin sur l'axe de la dent). Les arcs se déduisent des points (`krbgs`, tangente de départ
    -pi/Z), comme dans les fichiers .Eno et .Eng.
    """

    nom: str
    mo: float
    zz: int
    points: np.ndarray  # (n,2)
    drapeaux: list  # 0 aux extrémités, 1 entre les deux (Eno)
    e_mod: float = None
    poisson: float = None
    epaisseur: float = None
    points_dos: np.ndarray = None  # demi-profil arrière si la dent est dissymétrique (sinon None : symétrique)

    @property
    def symetrique(self):
        return self.points_dos is None

    def arcs(self):
        """(centres, rayons signés) des arcs tangents passant par les points."""
        pm, ra, _ = krbgs(angle_depart(self.zz), np.asarray(self.points, dtype=float))
        return pm, ra

    def echantillonner(self, n_total=400):
        pm, ra = self.arcs()
        return echantillonner(self.points, pm, ra, n_total)

    def normales(self):
        return normales(self.points, self.arcs()[0])

    def profil_dos(self):
        """Points du demi-profil arrière (celui du profil menant si la dent est symétrique)."""
        return self.points if self.points_dos is None else self.points_dos

    def arcs_dos(self):
        pm, ra, _ = krbgs(angle_depart(self.zz), np.asarray(self.profil_dos(), dtype=float))
        return pm, ra

    def echantillonner_dos(self, n_total=400):
        pm, ra = self.arcs_dos()
        return echantillonner(self.profil_dos(), pm, ra, n_total)

    def normales_dos(self):
        return normales(self.profil_dos(), self.arcs_dos()[0])


def vers_mobile(roue, profil_arcs, role):
    """Mobile (`role` : "m" menante, "e" menée) d'une roue de evolvente et de ses arcs."""
    n = profil_arcs.n_points
    return Mobile(f"Developpante-{role}", roue.mo, roue.zz, profil_arcs.points.copy(),
                  [0 if i in (0, n - 1) else 1 for i in range(n)],
                  roue.e_mod, roue.poisson, roue.epaisseur)


# paramètres utilisés par la suite du programme, avec les valeurs de l'exemple profile.json
EXTRAS_DEFAUT = {"dA_neg": -30.0, "dA_pos": 30.0, "h_contact": 120.0,
                 "simulation": {"opti_1": "couple", "opti_2": "PHZ", "ponderation": 99999, "centrage": 0}}


CLE_DOS = "profile_arriere"  # demi-profil arrière d'une dent dissymétrique (sym: false), même forme que "profile"


def _flag(v):
    return 1 if v else 0


def _mobile_depuis_dict(d):
    rows = d["profile"]
    pts = np.array([[float(r[0]), float(r[1])] for r in rows])
    flags = [_flag(r[2]) if len(r) > 2 else 0 for r in rows]
    dos = None
    if not d.get("sym", True) and d.get(CLE_DOS):  # dent dissymétrique : deuxième demi-profil
        dos = np.array([[float(r[0]), float(r[1])] for r in d[CLE_DOS]])
    return Mobile(str(d.get("nom", "")), float(d["m"]), int(d["z"]), pts, flags,
                  d.get("E") or None, d.get("poisson") or None, d.get("epaisseur") or None, dos)


def lire_arcs(chemin):
    """Lit un fichier d'arcs de cercle : (mobile menant, mobile mené), extras).

    Formats : .json au schéma `profile.json` (liste `Engrenage` / `Mobiles` / `Simulation`),
    .Eng / .json de l'ancien schéma (`{"m": ..., "s": ..., "sim": ...}`) et .Eno (texte du
    Fortran). `extras` : plage de variation d'entraxe (dA_neg, dA_pos en µm), h_contact et
    paramètres de simulation, à conserver pour la suite du programme.
    """
    chemin = Path(chemin)
    extras = {**EXTRAS_DEFAUT, "simulation": dict(EXTRAS_DEFAUT["simulation"])}
    if chemin.suffix.lower() == ".eno":
        return _lire_eno(chemin, extras)
    doc = json.loads(chemin.read_text(encoding="utf-8"))
    if isinstance(doc, list):
        items = {it.get("type"): it.get("details", {}) for it in doc}
        vec = items["Mobiles"]["mobile_vec"]
        maitre = next(m for m in vec if m.get("is_master"))
        esclave = next(m for m in vec if not m.get("is_master"))
        for cle in ("dA_neg", "dA_pos", "h_contact"):
            if cle in items.get("Engrenage", {}):
                extras[cle] = float(items["Engrenage"][cle])
        extras["simulation"].update(items.get("Simulation", {}))
    elif isinstance(doc, dict) and "m" in doc and "s" in doc:
        maitre, esclave = doc["m"], doc["s"]
        sim = doc.get("sim", {})
        for cle in ("dA_neg", "dA_pos"):
            if cle in sim:
                extras[cle] = float(sim[cle])
        if "ponderation" in sim:
            extras["h_contact"] = float(sim["ponderation"])
        for cle in ("opti_1", "opti_2", "centrage"):
            if cle in sim:
                extras["simulation"][cle] = sim[cle]
    else:
        raise ValueError("fichier d'arcs de cercle non reconnu")
    return (_mobile_depuis_dict(maitre), _mobile_depuis_dict(esclave)), extras


def _lire_eno(chemin, extras):
    lignes = chemin.read_text(encoding="latin-1").splitlines()
    mobiles, k = [], 0
    for role in ("m", "e"):
        nom = lignes[k].strip()
        mo, zz, n = lignes[k + 2].split()[:3]
        rows = [l.split() for l in lignes[k + 3:k + 3 + int(n)]]
        pts = np.array([[float(r[0]), float(r[1])] for r in rows])
        mobiles.append(Mobile(nom, float(mo), int(zz), pts, [_flag(int(r[2])) for r in rows]))
        k += 3 + int(n)
    hertz = [[float(x) for x in l.split()[:3]] for l in lignes[k + 1:k + 3]]
    for m, h in zip(mobiles, hertz):
        m.e_mod, m.poisson, m.epaisseur = h
    try:
        dneg, dpos, gw = (float(x) for x in lignes[k + 4].split()[:3])
        extras.update(dA_neg=dneg, dA_pos=dpos, h_contact=gw)
    except (IndexError, ValueError):
        pass
    return tuple(mobiles), extras


def ecrire_json(chemin, mobiles, extras=None):
    """Écrit les mobiles en arcs de cercle au schéma de `profile.json` (liste `Engrenage`,
    `Mobiles`, `Simulation`), avec en plus les centres des arcs (`arc_centre`) comme dans les
    fichiers .Eng."""
    extras = extras or EXTRAS_DEFAUT
    vec = []
    for k, m in enumerate(mobiles):
        pm, ra = m.arcs()
        d = {"is_master": k == 0, "nom": m.nom, "z": int(m.zz), "m": float(m.mo), "sym": m.symetrique,
             "profile": [[float(x), float(y), int(f)] for (x, y), f in zip(m.points, m.drapeaux)],
             "arc_rayon": [float(r) for r in ra], "arc_centre": [[float(x), float(y)] for x, y in pm]}
        if not m.symetrique:  # deuxième demi-profil, avec ses propres rayons et centres
            pmd, rad = m.arcs_dos()
            d[CLE_DOS] = [[float(x), float(y), int(f)] for (x, y), f in zip(m.points_dos, m.drapeaux)]
            d["arc_rayon_arriere"] = [float(r) for r in rad]
            d["arc_centre_arriere"] = [[float(x), float(y)] for x, y in pmd]
        for cle, v in (("epaisseur", m.epaisseur), ("E", m.e_mod), ("poisson", m.poisson)):
            if v is not None:
                d[cle] = v
        d["ty"] = "gear"
        vec.append(d)
    doc = [
        {"id": 1, "type": "Engrenage", "details": {
            "dA_pos": extras["dA_pos"], "dA_neg": extras["dA_neg"], "h_contact": extras["h_contact"]}},
        {"id": 2, "type": "Mobiles", "details": {"mobile_vec": vec}},
        {"id": 3, "type": "Simulation", "details": dict(extras["simulation"])},
    ]
    Path(chemin).write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def ecrire_eno(chemin, mobiles, extras=None):
    """Écrit les mobiles en arcs de cercle au format .Eno de ENGCRE (Engcre.F90).

    Par mobile : nom, en-tête, `module ZZ nombre de points`, puis les points de transition
    (x, y, drapeau) ; ensuite les données de Hertz des deux mobiles et la ligne des paramètres
    de simulation. Colonnes fixes, fins de ligne CRLF.
    """
    extras = extras or EXTRAS_DEFAUT
    if not all(m.symetrique for m in mobiles):
        raise ValueError("dent dissymétrique : le format .Eno ne contient qu'un demi-profil, exportez en .json")
    lignes = []
    for m in mobiles:
        lignes += [m.nom, "   Module |      ZZ |   Nbr.point|",
                   f"{m.mo:10.6f}{m.zz:10d}{len(m.points):10d}"]
        lignes += [f"{x:12.7f}{y:12.7f}{f:6d}" for (x, y), f in zip(m.points, m.drapeaux)]
    lignes.append("   E-Mod.  Coeff.Poiss. Epaiss.")
    for m in mobiles:
        if None in (m.e_mod, m.poisson, m.epaisseur):
            raise ValueError("données de Hertz manquantes (E-module, Poisson, épaisseur)")
        lignes.append(f"{m.e_mod:10.1f}{m.poisson:10.3f}{m.epaisseur:10.3f}")
    lignes += ["    -dA   |   +dA   |    GW   |    ZV   | ",
               f"{extras['dA_neg']:10.1f}{extras['dA_pos']:10.1f}{extras['h_contact']:8.0f}{1:11d}"]
    with open(chemin, "w", encoding="ascii", errors="replace", newline="") as f:
        f.write("\r\n".join(lignes) + "\r\n")
