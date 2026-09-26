# physique.py
"""
Solveur elastica unique du programme (méthode de tir : on impose l'angle du
balancier au collet, on cherche par optimisation locale (Fx, Fy, M0) au piton
qui amène la fibre neutre intégrée depuis le piton jusqu'à la bonne position
et au bon angle au collet). Remplace les deux implémentations quasi
identiques qui existaient dans spiran_v1.py et spiropt_v1/physics.py.

Utilisé de façon identique que la GeometrieSpiral provienne d'un DXF importé
ou d'une génération paramétrique : la physique ne connaît que
(x, y, epaisseur, hauteur), jamais l'origine de la géométrie.
"""
from dataclasses import replace

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import minimize as scipy_minimize, root

from .constantes import E_SILICIUM, RHO_SILICIUM
from .modele import GeometrieSpiral
from .collision import jeu_min_spires

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _integrer_tir_numba(Fx, Fy, M0, phi0, Xm0, Ym0, ds, raideur_flexion,
                         x_piton, y_piton, phi_piton, x_cible, y_cible, phi_cible):
    """Intègre la lame du piton vers le collet pour un torseur d'effort
    (Fx, Fy, M0) imposé au piton.

    Le bras de levier du moment est pris sur la forme DÉFORMÉE, pas sur la
    forme au repos (correction du 2026-09-24) : l'équilibre statique s'écrit
    M(s) = M0 + Fx.(Y_def(s) - y_piton) - Fy.(X_def(s) - x_piton), avec les
    coordonnées du point tel qu'il est RÉELLEMENT placé sous charge. Utiliser
    les coordonnées au repos revient à supposer les déplacements petits devant
    les bras de levier -- hypothèse fausse ici, où le collet tourne jusqu'à
    340° : mesuré sur un spiral de référence, cela surestimait l'énergie
    élastique de 7e-5 à 40° mais de 2.2e-2 à 340°, et déplaçait le centre de
    gravité de 183 µm. Une erreur qui CROÎT avec l'amplitude est exactement ce
    qui se lit ensuite comme un défaut d'isochronisme : l'essentiel des
    ~650 s/j d'anisochronisme que donnait le modèle précédent venait de là.

    Les deux intégrations (angle et position) sont menées dans la MÊME boucle,
    ce qui rend le bras de levier disponible au fur et à mesure -- pas de
    point fixe à itérer, donc pas de surcoût : un prédicteur explicite au pas
    précédent suffit à estimer la position du pas courant, le correcteur
    recalcule le moment là où le point se trouve réellement. Vérifié contre
    une résolution par point fixe menée jusqu'à convergence (10 à 11
    itérations) : écart relatif résiduel sur l'énergie de l'ordre de 1e-9.
    """
    n = len(phi0)
    phi_def = np.zeros(n)
    X_def = np.zeros(n)
    Y_def = np.zeros(n)
    M_def = np.zeros(n)
    phi_def[-1] = phi_piton
    X_def[-1] = x_piton
    Y_def[-1] = y_piton
    M_def[-1] = M0

    for i in range(n - 2, -1, -1):
        dphi0 = phi0[i + 1] - phi0[i]
        h = ds[i] * 1e-3

        # Prédicteur : moment évalué là où se trouve le pas déjà intégré.
        M_pred = M0 + Fx * (Y_def[i + 1] - y_piton) - Fy * (X_def[i + 1] - x_piton)
        phi_pred = phi_def[i + 1] - dphi0 - (M_pred / raideur_flexion[i]) * h
        X_pred = X_def[i + 1] - ds[i] * np.cos(phi_pred)
        Y_pred = Y_def[i + 1] - ds[i] * np.sin(phi_pred)

        # Correcteur : moment recalculé à la position estimée du pas courant.
        M = M0 + Fx * (Y_pred - y_piton) - Fy * (X_pred - x_piton)
        M_def[i] = M
        phi_def[i] = phi_def[i + 1] - dphi0 - (M / raideur_flexion[i]) * h
        X_def[i] = X_def[i + 1] - ds[i] * np.cos(phi_def[i])
        Y_def[i] = Y_def[i + 1] - ds[i] * np.sin(phi_def[i])

    cout = (X_def[0] - x_cible) ** 2 + (Y_def[0] - y_cible) ** 2 + ((phi_def[0] - phi_cible) * 0.1) ** 2
    return cout, X_def, Y_def, phi_def, M_def


def resoudre_elastica(geom: GeometrieSpiral, angle_deg: float, r_min: float,
                       return_diagnostics: bool = False):
    """Résout l'elastica non-linéaire pour une rotation `angle_deg` du
    balancier (le collet, index 0, tourne de cet angle ; le piton, index -1,
    reste fixe). r_min sert à replacer le point cible du collet sur le cercle
    de rayon r_min (convention identique à l'origine : le collet est solidaire
    de l'axe du balancier).

    Retourne (Xe, Ye, Xi, Yi, x_cg, y_cg) : bords extérieur/intérieur déformés
    et centre de gravité (mm), plus un dict de diagnostics si demandé.

    Sur un spiral à RAQUETTE, seule la partie active (collet -> point de
    comptage) est résolue : au-delà, la lame est tenue entre les goupilles et
    le piton, ses deux extrémités sont donc fixes et elle ne se déforme pas.
    Cette portion morte est simplement recollée telle quelle aux bords rendus,
    pour que le contrôle de jeu inter-spires et l'affichage voient bien la
    pièce entière -- elle occupe de la place même si elle ne travaille pas.
    """
    geom_complet = geom
    geom = geom.partie_active()
    Xm0, Ym0 = geom.x, geom.y
    # Épaisseur FINIE (oxyde compris) pour tout ce qui est encombrement : c'est
    # la pièce oxydée qui occupe l'espace entre les spires.
    epaisseurs = geom.epaisseur_finie
    alpha = np.radians(angle_deg)
    dX0, dY0 = np.gradient(Xm0), np.gradient(Ym0)
    ds = np.sqrt(dX0 ** 2 + dY0 ** 2)
    masse_elementaire = geom.masse_elementaire(RHO_SILICIUM)

    phi0 = np.unwrap(np.arctan2(dY0, dX0))
    # Rigidité du COMPOSITE silicium/oxyde ; identique à E_silicium·I quand il
    # n'y a pas d'oxyde.
    raideur_flexion = geom.rigidite_flexion()

    theta_init = np.arctan2(Ym0[0], Xm0[0])
    x_cible = r_min * np.cos(theta_init + alpha)
    y_cible = r_min * np.sin(theta_init + alpha)
    phi_cible = phi0[0] + alpha
    x_piton, y_piton, phi_piton = Xm0[-1], Ym0[-1], phi0[-1]

    def integrer_tir(parametres):
        Fx, Fy, M0 = parametres[0] * 1e-3, parametres[1] * 1e-3, parametres[2] * 1e-6
        return _integrer_tir_numba(Fx, Fy, M0, phi0, Xm0, Ym0, ds, raideur_flexion,
                                    x_piton, y_piton, phi_piton, x_cible, y_cible, phi_cible)

    def residu_tir(parametres):
        """Les 3 conditions au collet : position (x, y) et angle de la
        tangente. Autant d'équations que d'inconnues (Fx, Fy, M0) -- le
        problème de tir est exactement déterminé, sa solution exacte annule
        les trois."""
        _, X, Y, phi, _ = integrer_tir(parametres)
        return np.array([X[0] - x_cible, Y[0] - y_cible, phi[0] - phi_cible])

    souplesse = np.sum((ds * 1e-3) / raideur_flexion)
    p0 = [0.0, 0.0, (alpha / souplesse) * 1e6 if abs(alpha) > 1e-5 else 0.0]

    # RÉSOLUTION du système (Powell hybride), et non minimisation d'un coût
    # agrégé (diagnostic 2026-09-24). La minimisation BFGS du coût
    # (dx^2 + dy^2 + (0.1.dphi)^2) utilisée auparavant s'arrêtait sur des
    # résidus de ~0.1 µm en position et jusqu'à 0.24° en ANGLE : la condition
    # angulaire, sous-pondérée d'un facteur 10 dans le coût, n'était en
    # pratique presque pas imposée. Sans conséquence visible sur l'orbite du
    # centre de gravité (effet du premier ordre), mais rédhibitoire pour
    # l'énergie élastique U(alpha), qui s'en trouvait faussée jusqu'à 2.4e-3
    # en relatif -- soit ~10x l'écart à l'isochronisme qu'on cherche
    # justement à mesurer dessus (8 s/j ~ 2e-4 en relatif sur U). La
    # résolution directe ramène les résidus à ~1e-11 µm et ~1e-13 degré, à
    # coût de calcul équivalent (mesuré : 3 à 7 ms par angle dans les deux
    # cas), ce qui rend l'isochronisme propre du spiral enfin exploitable.
    res_tir = root(residu_tir, p0, method="hybr")
    cout_final, X_f, Y_f, phi_f, M_f = integrer_tir(res_tir.x)

    # Garde-fou anti-divergence (diagnostic 2026-09-25, étude de convergence
    # de maillage) : le tir peut, pour un angle/maillage particulier, partir
    # sur une branche aberrante (résidu de plusieurs dizaines de µm, contre
    # une fraction de µm partout ailleurs) -- observé concrètement à un angle
    # donné pour un seul n de maillage sur ~35 testés, mais suffisant pour
    # corrompre la spline de l'orbite du CG et fausser silencieusement la
    # dispersion de marche calculée à partir d'elle (269 s/j au lieu de
    # ~13 s/j sur le cas ayant révélé le problème). Deux stratégies de repli
    # sont essayées et la meilleure est retenue : le même solveur reparti d'un
    # autre point (M0=0 au lieu de l'estimation linéaire), puis, en dernier
    # recours, la minimisation du coût agrégé -- moins précise, mais qui ne
    # peut pas diverger, là où une résolution de système le peut. Ne coûte
    # rien dans le cas normal (déclenché seulement si le premier essai a
    # clairement divergé).
    erreur_um = float(np.hypot(X_f[0] - x_cible, Y_f[0] - y_cible) * 1000.0)
    if erreur_um > 1.0:
        for essai in (root(residu_tir, [0.0, 0.0, 0.0], method="hybr"),
                      scipy_minimize(lambda ps: integrer_tir(ps)[0], p0, method="BFGS",
                                     options={"gtol": 1e-12, "maxiter": 300})):
            cout_alt, X_alt, Y_alt, phi_alt, M_alt = integrer_tir(essai.x)
            if cout_alt < cout_final:
                res_tir, cout_final = essai, cout_alt
                X_f, Y_f, phi_f, M_f = X_alt, Y_alt, phi_alt, M_alt

    # Énergie élastique de flexion emmagasinée, U(alpha) = somme M^2/(2.EI).ds,
    # avec le champ de moment issu de l'intégration elle-même (donc sur bras de
    # levier déformé, cf. _integrer_tir_numba). C'est la grandeur qui porte
    # TOUTE l'information d'isochronisme propre du spiral : la marche en
    # fonction de l'amplitude, dans TOUTES les positions (y compris
    # horizontales, où le terme de gravité est nul par construction), se déduit
    # de la seule courbe U(alpha) -- cf. isochronisme.
    # marche_isochronisme_intrinseque. M est en N.m par construction de
    # l'équation d'équilibre du tir (dphi = -(M/EI).ds avec ds en m), quelles
    # que soient les mises à l'échelle 1e-3/1e-6 appliquées aux paramètres
    # de tir Fx/Fy/M0 (qui ne sont que du préconditionnement numérique).
    energie_J = float(np.sum((M_f ** 2) / (2.0 * raideur_flexion) * (ds * 1e-3)))

    dX_d, dY_d = np.gradient(X_f), np.gradient(Y_f)
    norme_d = np.sqrt(dX_d ** 2 + dY_d ** 2)
    nX, nY = -dY_d / norme_d, dX_d / norme_d

    x_cg = np.sum(X_f * masse_elementaire) / np.sum(masse_elementaire)
    y_cg = np.sum(Y_f * masse_elementaire) / np.sum(masse_elementaire)

    Xe = X_f + (epaisseurs / 2.0) * nX
    Ye = Y_f + (epaisseurs / 2.0) * nY
    Xi = X_f - (epaisseurs / 2.0) * nX
    Yi = Y_f - (epaisseurs / 2.0) * nY

    if geom_complet.a_partie_morte:
        ext_morte, int_morte = geom_complet.contour_flancs()
        n_a = geom_complet.n_actif
        Xe = np.concatenate([Xe, ext_morte[n_a:, 0]])
        Ye = np.concatenate([Ye, ext_morte[n_a:, 1]])
        Xi = np.concatenate([Xi, int_morte[n_a:, 0]])
        Yi = np.concatenate([Yi, int_morte[n_a:, 1]])

    if not return_diagnostics:
        return Xe, Ye, Xi, Yi, x_cg, y_cg

    erreur_position_um = float(np.hypot(X_f[0] - x_cible, Y_f[0] - y_cible) * 1000.0)
    erreur_angle_deg = float(np.degrees(abs(((phi_f[0] - phi_cible) + np.pi) % (2 * np.pi) - np.pi)))
    diagnostics = {
        "angle_deg": float(angle_deg), "success": bool(res_tir.success),
        "message": str(res_tir.message),
        "n_iter": int(getattr(res_tir, "nit", getattr(res_tir, "nfev", 0))),
        "cout_final": float(cout_final), "erreur_position_um": erreur_position_um,
        "erreur_angle_deg": erreur_angle_deg, "energie_J": energie_J,
        # Résultante que le spiral applique au piton, en newtons. Le collet
        # subit l'opposée : c'est elle qui pousse le pivot du balancier contre
        # sa pierre, et sa DIRECTION par rapport à la ligne des centres de
        # l'échappement est ce que vise la règle classique du nombre de spires
        # entier. Les paramètres de tir sont en newtons par construction (le
        # moment M = F.bras vaut des N.m avec le bras en mm et le facteur
        # 1e-3 appliqué au paramètre).
        "force_piton_N": (float(res_tir.x[0]), float(res_tir.x[1])),
    }
    return Xe, Ye, Xi, Yi, x_cg, y_cg, diagnostics


def balayage_angulaire(geom: GeometrieSpiral, r_min: float, amplitude_max_deg: float,
                        pas_angulaire_deg: float = 10.0):
    """Balaie l'angle du balancier sur [-amplitude_max_deg, +amplitude_max_deg]
    et résout l'elastica à CHAQUE angle -- une seule fois, réutilisé à la fois
    pour l'orbite du CG (isochronisme.py) et pour le jeu inter-spires
    (optimisation.py), qui avant cette fusion étaient calculés par deux
    balayages séparés : l'orbite sur toute la plage (pas_angulaire_deg) mais
    le jeu à seulement 3 angles (0, +max, -max, cf. optimisation._jeu_a_angle,
    supprimée). Un spiral peut atteindre son jeu minimal à un angle
    INTERMÉDIAIRE (l'elastica n'est pas monotone en angle près du minimum de
    jeu) -- ne vérifier que les 3 extrêmes pouvait donc laisser passer une
    collision réelle non détectée par l'optimiseur (diagnostic 2026-09-23).

    Retourne (angles_deg, x_cg_mm, y_cg_mm, jeu_mm, energie_J) -- jeu_mm est
    le jeu inter-spires (mm) à CHAQUE angle balayé, permettant de retrouver
    ensuite l'angle de jeu minimal exact via angles_deg[np.argmin(jeu_mm)] ;
    energie_J est l'énergie élastique emmagasinée à chaque angle, dont se
    déduit l'isochronisme propre du spiral (marche en fonction de l'amplitude,
    valable dans TOUTES les positions y compris horizontales) -- même
    mutualisation que pour le jeu : ce sont trois exploitations différentes du
    MÊME balayage elastica, qui est la seule partie coûteuse.

    Nombre de points ajusté pour être impair (theta=0 exactement au centre,
    requis par les métriques de symétrie ET par le calcul d'isochronisme, qui
    a besoin du minimum d'énergie U(0)=0 comme référence exacte)."""
    n_side = max(1, round(amplitude_max_deg / pas_angulaire_deg))
    angles_deg = np.linspace(-amplitude_max_deg, amplitude_max_deg, 2 * n_side + 1)
    x_cg = np.zeros_like(angles_deg)
    y_cg = np.zeros_like(angles_deg)
    jeu_mm = np.zeros_like(angles_deg)
    energie_J = np.zeros_like(angles_deg)
    for i, ang in enumerate(angles_deg):
        Xe, Ye, Xi, Yi, x_c, y_c, diag = resoudre_elastica(geom, ang, r_min, return_diagnostics=True)
        x_cg[i] = x_c
        y_cg[i] = y_c
        # Épaisseur FINIE : c'est la lame oxydée qui doit passer entre les
        # spires, pas la cote de gravure. L'oxyde élargit la lame de ~4 µm,
        # soit un tiers du jeu minimal exigé -- l'oublier rendait le contrôle
        # de collision optimiste d'autant.
        jeu_mm[i] = jeu_min_spires((Xe + Xi) / 2.0, (Ye + Yi) / 2.0, geom.epaisseur_finie)
        energie_J[i] = diag["energie_J"]
    return angles_deg, x_cg, y_cg, jeu_mm, energie_J


def orbite_cg(geom: GeometrieSpiral, r_min: float, amplitude_max_deg: float,
              pas_angulaire_deg: float = 10.0):
    """Compatibilité : identique à balayage_angulaire sans le jeu ni l'énergie,
    pour les appelants qui n'ont besoin que de l'orbite du CG."""
    angles_deg, x_cg, y_cg, _, _ = balayage_angulaire(geom, r_min, amplitude_max_deg, pas_angulaire_deg)
    return angles_deg, x_cg, y_cg


def souplesse_integree(s_m: np.ndarray, e_m: np.ndarray) -> float:
    """Souplesse intégrée S = integrale( ds / e³ ), en m⁻², par la règle du
    trapèze sur l'abscisse curviligne.

    Quadrature du trapèze et non somme de np.gradient(s) : la seconde suppose
    implicitement un pas régulier et devient fausse là où la densité du
    maillage change brusquement. C'est exactement le cas d'un spiral à
    raquette, dont l'arc concentrique est discrétisé bien plus finement que le
    corps (mesuré : 5 µm par élément sur l'arc contre 64 µm dans le corps), et
    dont le point de comptage tombe justement dans cette zone. L'erreur de
    quadrature y atteignait 12% sur la plage de réglage calculée
    (diagnostic 2026-09-25). Le trapèze est exact sur un pas quelconque et se
    tronque proprement en tout point, ce qui est indispensable pour comparer
    deux longueurs actives différentes."""
    return float(trapezoid(1.0 / np.asarray(e_m) ** 3, x=np.asarray(s_m)))


def frequence_propre(geom: GeometrieSpiral, inertie_balancier: float) -> float:
    """Fréquence propre (Hz) de l'oscillateur balancier-spiral pour CETTE
    géométrie, quelle que soit son origine (DXF importé, génération
    paramétrique, solution optimisée) -- même formule (poutre linéaire,
    couple = E.h / (12.souplesse_int)) que la mise à l'échelle automatique de
    parametrage.construire_geometrie, mais réutilisable après coup sur une
    géométrie déjà figée, pour afficher où l'oscillateur atterrit réellement
    par rapport à la fréquence cible visée.

    Sur un spiral à raquette, seule la partie active compte : la portion
    au-delà du point de comptage n'apporte aucune raideur. La rigidité est
    celle du composite silicium/oxyde quand la lame est oxydée."""
    couple = 1.0 / geom.souplesse()
    return (1.0 / (2.0 * np.pi)) * np.sqrt(couple / inertie_balancier)


def marche_raquette(geom: GeometrieSpiral, inertie_balancier: float, n_pas_arc: int):
    """Plage de marche (s/j) réellement couverte par la raquette, mesurée en
    déplaçant le point de comptage d'un bout à l'autre de l'arc concentrique.

    Retourne (marche_rapide, marche_lente, ecart_total) : la marche obtenue
    quand les goupilles sont poussées au plus près du corps (spiral actif le
    plus court, montre la plus rapide), quand elles sont au plus près du piton
    (le plus long, la plus lente), et l'étendue totale entre les deux. Les deux
    valeurs sont rapportées au réglage NOMINAL (goupilles au milieu de l'arc).

    C'est une VÉRIFICATION, pas une construction : l'angle de l'arc a été
    dimensionné analytiquement par parametrage.angle_arc_raquette, et cette
    fonction contrôle ce que la géométrie obtenue donne réellement. Même
    logique que la mise à l'échelle de l'épaisseur, calculée analytiquement
    puis vérifiée sur la géométrie finale."""
    if not geom.a_partie_morte or n_pas_arc < 3:
        return 0.0, 0.0, 0.0
    # Le parcours des goupilles est centré sur le point de comptage, et il ne
    # va PAS jusqu'au bout de la géométrie : l'arc se prolonge au-delà jusqu'au
    # piton (portion que la raquette n'atteint jamais). Les bornes se déduisent
    # donc de n_actif, pas de la fin du tableau.
    n_arc = n_pas_arc | 1
    i_debut = max(geom.n_actif - 1 - n_arc // 2, 4)
    i_fin = min(i_debut + n_arc - 1, geom.n - 1)
    f_nominal = frequence_propre(geom, inertie_balancier)
    f_rapide = frequence_propre(replace(geom, n_actif=i_debut + 1), inertie_balancier)
    f_lente = frequence_propre(replace(geom, n_actif=i_fin + 1), inertie_balancier)
    marche_rapide = 86400.0 * (f_rapide / f_nominal - 1.0)
    marche_lente = 86400.0 * (f_lente / f_nominal - 1.0)
    return float(marche_rapide), float(marche_lente), float(marche_rapide - marche_lente)


def axe_force_pivot(geom: GeometrieSpiral, r_min: float, amplitude_max_deg: float):
    """Résultante que le spiral applique au pivot du balancier : module maximal
    (µN) et orientation de son AXE (deg, modulo 180°).

    Cette force n'intervient pas dans la marche calculée ici -- elle agit par
    le frottement du pivot dans sa pierre, qui n'est pas modélisé. Elle est
    exposée parce qu'elle est la justification physique de la règle horlogère
    du nombre entier de spires : mesuré sur un calibre de référence, son axe
    suit la droite collet -> point de comptage à 1 à 4° près, et son module
    atteint treize fois le poids du spiral. Placer cette droite sur la ligne
    des centres de l'échappement fixe donc l'orientation de cette force par
    rapport à la réaction de l'ancre.

    L'axe est donné modulo 180° parce que la force s'inverse à chaque
    alternance : ce qui compte est sa direction, pas son sens.
    """
    angles = np.array([-amplitude_max_deg, 0.0, amplitude_max_deg])
    forces = []
    for a in angles:
        *_, diag = resoudre_elastica(geom, float(a), r_min, return_diagnostics=True)
        forces.append(diag["force_piton_N"])
    forces = np.array(forces)
    modules = np.hypot(forces[:, 0], forces[:, 1])
    i = int(np.argmax(modules))
    axe = float(np.degrees(np.arctan2(forces[i, 1], forces[i, 0])) % 180.0)
    return float(modules[i] * 1e6), axe


def angle_point_comptage(geom: GeometrieSpiral) -> float:
    """Position angulaire du point de comptage (deg, modulo 180°), pour la
    comparer à l'axe de la force : les deux doivent coïncider."""
    i = geom.n_actif - 1
    return float(np.degrees(np.arctan2(geom.y[i], geom.x[i])) % 180.0)
