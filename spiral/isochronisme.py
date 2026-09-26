# isochronisme.py
"""
Marche diurne d'une géométrie de spiral. TROIS contributions, de natures
différentes, calculées séparément puis sommées (analyser_marches) :

  1. GRAVITÉ -- déplacement du centre de gravité du spiral sous son propre
     poids, via la formule d'Airy (perturbation au premier ordre d'un
     oscillateur harmonique par un couple perturbateur périodique).
     Dépend de la POSITION : maximal dans les quatre positions verticales,
     rigoureusement NUL à plat (la gravité y est perpendiculaire au plan du
     spiral). C'est ce terme, et lui seul, qui crée la dispersion de marche
     entre positions -- l'objectif historique de ce projet.

  2. ISOCHRONISME PROPRE -- non-linéarité du couple de rappel du spiral en
     grande déformation, qui fait varier sa période avec l'amplitude.
     Indépendant de la position : agit à l'identique dans les six positions.
     Calculé par intégration de la période exacte de l'oscillateur
     conservatif à partir de la courbe d'énergie élastique U(alpha), PAS par
     perturbation -- l'effet cherché est justement l'anharmonicité.

  3. ÉCHAPPEMENT -- retard aux faibles amplitudes (cf. echappement.py).
     Également indépendant de la position.

La séparation n'est pas cosmétique : comme (1) est nul à plat et que (3) ne
l'est pas, l'erreur d'échappement ne peut être compensée QUE par (2). Les
deux objectifs d'optimisation qui en découlent sont donc indépendants, et
non deux lectures de la même grandeur.

x_cg(phi)/y_cg(phi) et U(alpha) viennent du VRAI solveur elastica
(physique.balayage_angulaire), pas d'une approximation à courbure uniforme
(invalide dès que l'épaisseur varie le long de la lame, cf. note historique
du projet : écart mesuré ~50x sur la dispersion de l'orbite CG entre les deux
approches).
"""
import numpy as np
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from .constantes import G_TERRE, RHO_SILICIUM
from .echappement import ReglageEchappement, marche_echappement
from .modele import GeometrieSpiral
from .physique import balayage_angulaire
from .raquette import analyser_jeu_goupilles

POSITIONS_VERTICALES_DEFAUT = {
    "3H": np.radians(180.0),
    "6H": np.radians(90.0),
    "9H": np.radians(0.0),
    "12H": np.radians(270.0),
}

# Les deux positions horizontales (cadran haut / cadran bas) sont, pour ce
# modèle, strictement identiques : la gravité y est perpendiculaire au plan du
# spiral, son couple perturbateur autour de l'axe du balancier est nul, donc
# le déplacement du centre de gravité n'y produit AUCUNE erreur de marche.
# Elles sont nommées séparément quand même, parce que ce qui s'y passe (erreur
# d'échappement + isochronisme propre du spiral) est précisément ce qu'on
# cherche à régler, et qu'un chronocomparateur les mesure séparément.
POSITIONS_HORIZONTALES = ("CH", "CB")


def marche_positions_verticales(geom: GeometrieSpiral, r_min: float, freq_cible: float,
                                 inertie_balancier: float, amplitudes_deg: np.ndarray,
                                 positions: dict = None, pas_balayage_deg: float = 10.0,
                                 return_diagnostics: bool = False):
    """Calcule la marche diurne (s/j) aux amplitudes et positions verticales
    données. Retourne dict{nom_position: array marches (s/j), même longueur
    que amplitudes_deg}, et si demandé un dict de diagnostics (angles_deg,
    x_cg_mm, y_cg_mm, jeu_mm du balayage déjà calculé) -- pour réutiliser ce
    balayage ailleurs (décentrage moyen et jeu inter-spires réel côté
    optimisation.py) sans ré-exécuter l'elastica complet une deuxième fois."""
    if positions is None:
        positions = POSITIONS_VERTICALES_DEFAUT

    amplitude_max = float(np.max(amplitudes_deg))
    angles_deg, x_cg_mm, y_cg_mm, jeu_mm, energie_J = balayage_angulaire(
        geom, r_min, amplitude_max, pas_balayage_deg)
    x_cg_m, y_cg_m = x_cg_mm * 1e-3, y_cg_mm * 1e-3
    masse_spiral = geom.masse_totale(RHO_SILICIUM)

    phi_rad_ech = np.radians(angles_deg)
    cs_x, cs_y = CubicSpline(phi_rad_ech, x_cg_m), CubicSpline(phi_rad_ech, y_cg_m)
    dcs_x, dcs_y = cs_x.derivative(), cs_y.derivative()

    omega = 2.0 * np.pi * freq_cible
    psi_sample = np.linspace(0, 2 * np.pi, 2000)
    amplitudes_rad = np.radians(amplitudes_deg)

    marches = {}
    for nom_pos, angle_g in positions.items():
        gx, gy = G_TERRE * np.cos(angle_g), G_TERRE * np.sin(angle_g)
        valeurs = []
        for Theta_0 in amplitudes_rad:
            phi_sample = Theta_0 * np.sin(psi_sample)
            Mg = masse_spiral * (gx * dcs_x(phi_sample) + gy * dcs_y(phi_sample))
            int_val = trapezoid(Mg * np.sin(psi_sample), psi_sample)
            delta_f_f = -1.0 / (2.0 * np.pi * (omega ** 2) * inertie_balancier * Theta_0) * int_val
            valeurs.append(-86400.0 * delta_f_f)
        marches[nom_pos] = np.array(valeurs)

    if not return_diagnostics:
        return marches

    diagnostics = {
        "angles_deg": angles_deg, "x_cg_mm": x_cg_mm, "y_cg_mm": y_cg_mm, "jeu_mm": jeu_mm,
        "energie_J": energie_J,
        "angle_jeu_min_deg": float(angles_deg[np.argmin(jeu_mm)]), "jeu_min_mm": float(np.min(jeu_mm)),
    }
    return marches, diagnostics


def _quart_periode(E: float, a_tournant: float, spline_U, derivee_U,
                    inertie_balancier: float, u: np.ndarray) -> float:
    """Quart de période : temps mis par le balancier pour aller de l'équilibre
    (theta=0) au point de rebroussement a_tournant, sur UN côté, pour le
    niveau d'énergie E.

        dt = sqrt(J/2) . d(theta) / sqrt(E - U(theta))

    L'intégrande diverge en 1/sqrt au point de rebroussement (divergence
    intégrable). Le changement de variable theta = a.sin(u), u dans [0, pi/2],
    la régularise complètement : le cos(u) du jacobien s'annule exactement à
    la même vitesse que le sqrt(E - U) au dénominateur. La limite au point de
    rebroussement vaut sqrt(2.a / U'(a)) et est imposée analytiquement (0/0 en
    flottant sinon). Pour un ressort LINÉAIRE l'intégrande devient exactement
    constant, égal à sqrt(2/k) : le quart de période vaut alors
    (pi/2).sqrt(J/k), indépendant de l'amplitude -- c'est le cas isochrone
    exact, retrouvé sans aucune erreur de quadrature, ce qui est précisément
    la propriété recherchée ici (on ne mesure que l'ÉCART à l'isochronisme,
    de l'ordre de 1e-4 en relatif : une quadrature qui ne serait pas exacte
    sur le cas linéaire noierait le signal utile dans son propre biais)."""
    theta = a_tournant * np.sin(u)
    ecart = E - spline_U(theta)
    limite = np.sqrt(2.0 * a_tournant / derivee_U(a_tournant))
    with np.errstate(divide="ignore", invalid="ignore"):
        integrande = a_tournant * np.cos(u) / np.sqrt(ecart)
    # Point de rebroussement + éventuels points où l'erreur d'interpolation
    # rend `ecart` négatif juste avant lui : on y impose la limite analytique.
    integrande[~np.isfinite(integrande)] = limite
    integrande[ecart <= 0.0] = limite
    integrande[-1] = limite
    return float(np.sqrt(inertie_balancier / 2.0) * simpson(integrande, x=u))


def marche_isochronisme_intrinseque(angles_deg: np.ndarray, energie_J: np.ndarray,
                                     inertie_balancier: float, amplitudes_deg: np.ndarray,
                                     amplitude_reference_deg: float,
                                     n_quadrature: int = 401):
    """Marche diurne (s/j) due au seul ISOCHRONISME PROPRE du spiral, c'est-à-
    dire à la non-linéarité de son couple de rappel en grande déformation.

    Contrairement au terme de gravité (marche_positions_verticales), ce terme
    NE DÉPEND PAS de l'orientation de la montre : il agit identiquement dans
    les quatre positions verticales ET les deux positions horizontales. C'est
    donc le SEUL levier dont dispose la conception du spiral pour compenser
    l'erreur d'échappement (cf. echappement.py), qui a la même propriété.

    Méthode -- énergie, pas perturbation. On ne peut pas utiliser ici la
    formule d'Airy du terme de gravité : celle-ci est une perturbation au
    premier ordre d'un oscillateur supposé harmonique, alors que l'effet
    cherché EST l'anharmonicité elle-même. On intègre donc la période exacte
    de l'oscillateur conservatif à partir de la courbe d'énergie élastique
    U(alpha) mesurée sur le balayage elastica :

        T(Theta_0) = 2 . [ quart de période côté + , + côté - ]

    Le spiral n'est PAS symétrique en armage/désarmage : U(+alpha) != U(-alpha)
    (une spire qui se serre ne se raidit pas comme une spire qui s'ouvre). Les
    deux points de rebroussement d'une oscillation ne sont donc pas symétriques
    -- ils sont à ÉNERGIE égale, pas à angle égal. Les deux côtés sont traités
    séparément à partir du même niveau d'énergie E, et l'amplitude rapportée
    est la demi-course (a+ + a-)/2, qui est ce que mesure un chronocomparateur.
    Symétriser U aurait supprimé cette asymétrie, qui contribue pourtant au
    défaut d'isochronisme.

    `amplitude_reference_deg` fixe le zéro de la courbe (la marche moyenne de
    la montre est de toute façon rattrapée au réglage : seule la VARIATION
    avec l'amplitude a un sens physique ici).

    Retourne (marches_s_par_j, diagnostics). Si la courbe d'énergie n'est pas
    exploitable (non strictement croissante de part et d'autre de l'équilibre
    -- géométrie dégénérée, ou balayage trop grossier), retourne des NaN
    plutôt qu'un chiffre faux, et le signale dans les diagnostics.
    """
    angles_deg = np.asarray(angles_deg, dtype=float)
    U = np.asarray(energie_J, dtype=float)
    theta = np.radians(angles_deg)
    amplitudes_deg = np.atleast_1d(np.asarray(amplitudes_deg, dtype=float))
    nan = np.full(amplitudes_deg.shape, np.nan)

    i0 = int(np.argmin(np.abs(theta)))
    U = U - U[i0]

    # Deux branches monotones de part et d'autre de l'équilibre, ramenées
    # toutes deux à une abscisse positive a = |theta|.
    a_p, U_p = theta[i0:], U[i0:]
    a_n, U_n = -theta[i0::-1], U[i0::-1]
    if len(a_p) < 4 or len(a_n) < 4:
        return nan, {"valide": False, "raison": "balayage trop grossier"}
    if np.any(np.diff(U_p) <= 0.0) or np.any(np.diff(U_n) <= 0.0):
        return nan, {"valide": False, "raison": "énergie non monotone"}

    spl_p, spl_n = CubicSpline(a_p, U_p), CubicSpline(a_n, U_n)
    der_p, der_n = spl_p.derivative(), spl_n.derivative()

    # Les deux points de rebroussement d'une oscillation sont à ÉNERGIE égale.
    # On les apparie en RÉSOLVANT les équations, pas en interpolant une table :
    # une table interpolée introduisait un plancher de bruit de ~0.2 à 0.4 s/j
    # (mesuré sur un ressort linéaire, dont la marche exacte vaut 0 par
    # construction) -- soit 3 à 5% d'une cible de compensation de 8 s/j.
    # CubicSpline.solve donne les racines exactes du polynôme par morceaux,
    # sans grille : le même test retombe alors à 4e-9 s/j.
    def rebroussement_negatif(E: float) -> float:
        racines = spl_n.solve(E, extrapolate=False)
        racines = racines[racines > 0.0]
        return float(racines[0]) if len(racines) else float("nan")

    def amplitude_mesuree(ap: float) -> float:
        """Demi-course (a+ + a-)/2 pour un rebroussement positif ap donné --
        c'est cette demi-course que mesure un chronocomparateur."""
        an = rebroussement_negatif(float(spl_p(ap)))
        return 0.5 * (ap + an)

    # Plage utilisable : bornée par la branche la moins étendue en énergie. La
    # marge relative évite que le point de rebroussement conjugué tombe
    # EXACTEMENT sur le dernier noeud, où `solve(extrapolate=False)` peut le
    # rejeter pour un ulp -- ce qui faisait échouer tout le calcul au seul
    # voisinage de l'amplitude maximale.
    E_max_commun = min(float(spl_p(a_p[-1])), float(spl_n(a_n[-1]))) * (1.0 - 1e-9)
    racines_max = spl_p.solve(E_max_commun, extrapolate=False)
    racines_max = racines_max[racines_max > 0.0]
    if len(racines_max) == 0:
        return nan, {"valide": False, "raison": "points de rebroussement introuvables"}
    ap_max = float(racines_max[-1])
    amplitude_max_atteignable = float(np.degrees(amplitude_mesuree(ap_max)))
    if not np.isfinite(amplitude_max_atteignable):
        return nan, {"valide": False, "raison": "points de rebroussement introuvables"}

    # Borne basse du bracketing. PAS une valeur infinitésimale : au voisinage
    # immédiat de l'équilibre, U est du second ordre alors que le résidu du
    # solveur de tir sur U(0) ne l'est pas (mesuré : U(0) vaut ~0.1% de
    # U(20°) avant recalage), et la spline peut y passer légèrement sous zéro
    # -- le point de rebroussement conjugué devient alors introuvable et tout
    # le calcul retournait NaN. On part du premier noeud du balayage, très en
    # dessous des amplitudes d'intérêt (200° et plus), en remontant tant que
    # l'appariement des deux points de rebroussement n'est pas défini.
    ap_min = a_p[1]
    for noeud in a_p[1:]:
        if np.isfinite(amplitude_mesuree(float(noeud))):
            ap_min = float(noeud)
            break
    else:
        return nan, {"valide": False, "raison": "points de rebroussement introuvables"}

    A_min, A_max = amplitude_mesuree(ap_min), amplitude_mesuree(ap_max)

    def rebroussement_positif(amplitude_deg: float) -> float:
        """L'amplitude mesurée est monotone croissante en ap : brentq est donc
        garanti de converger DÈS LORS que la cible est dans l'intervalle. Les
        deux bornes sont traitées à part, sans quoi une cible égale à la borne
        au bit près (cas courant : l'amplitude d'analyse est justement écrêtée
        sur l'amplitude maximale atteignable) ne présente plus de changement
        de signe et brentq échoue."""
        cible = np.radians(amplitude_deg)
        if cible <= A_min:
            return ap_min
        if cible >= A_max:
            return ap_max
        try:
            return float(brentq(lambda ap: amplitude_mesuree(ap) - cible,
                                ap_min, ap_max, xtol=1e-14, rtol=8.9e-16))
        except (ValueError, IndexError):
            return float("nan")

    u = np.linspace(0.0, np.pi / 2.0, n_quadrature | 1)  # impair, requis par Simpson

    def periode(amplitude_deg: float) -> float:
        ap = rebroussement_positif(amplitude_deg)
        if not np.isfinite(ap):
            return float("nan")
        E = float(spl_p(ap))
        an = rebroussement_negatif(E)
        return 2.0 * (_quart_periode(E, ap, spl_p, der_p, inertie_balancier, u)
                      + _quart_periode(E, an, spl_n, der_n, inertie_balancier, u))

    amp_ref = min(float(amplitude_reference_deg), amplitude_max_atteignable)
    amplitudes_evaluees = np.minimum(amplitudes_deg, amplitude_max_atteignable)

    T_ref = periode(amp_ref)
    if not np.isfinite(T_ref) or T_ref <= 0.0:
        return nan, {"valide": False, "raison": "période de référence non calculable"}

    T = np.array([periode(a) for a in amplitudes_evaluees])
    marches = 86400.0 * (T_ref - T) / T_ref

    # Fréquence des petites oscillations : lue sur la COURBURE de U à
    # l'équilibre (k = U''(0)), pas sur une période intégrée à amplitude
    # quasi nulle, où la quadrature perd toute précision relative.
    k0 = 0.5 * (float(spl_p.derivative(2)(0.0)) + float(spl_n.derivative(2)(0.0)))
    diagnostics = {
        "valide": bool(np.all(np.isfinite(marches))),
        "frequence_petites_oscillations_hz": (
            float(np.sqrt(k0 / inertie_balancier) / (2.0 * np.pi)) if k0 > 0 else float("nan")),
        "amplitude_reference_deg": amp_ref,
        "amplitude_max_atteignable_deg": amplitude_max_atteignable,
        "periode_reference_s": T_ref,
    }
    return marches, diagnostics


def composer_marches(marches_gravite: dict, marche_iso: np.ndarray,
                      marche_ech: np.ndarray, marche_jeu: np.ndarray = None) -> dict:
    """Assemble la marche TOTALE réellement observée au chronocomparateur :

        verticales  : gravité(position, amplitude) + isochronisme + échappement
        horizontales: 0                            + isochronisme + échappement

    C'est cette somme, et elle seule, qui correspond à ce qu'on mesure. Les
    trois termes sont conservés séparément par ailleurs (affichage, réglage) :
    seul le premier dépend de la position, seuls les deux derniers existent à
    plat -- et c'est exactement pour ça que l'échappement ne peut être
    compensé que par l'isochronisme propre."""
    commun = marche_iso + marche_ech
    if marche_jeu is not None:
        commun = commun + marche_jeu
    totales = {nom: valeurs + commun for nom, valeurs in marches_gravite.items()}
    for nom in POSITIONS_HORIZONTALES:
        totales[nom] = commun
    return totales


def analyser_marches(geom: GeometrieSpiral, r_min: float, freq_cible: float,
                      inertie_balancier: float, amplitudes_deg: np.ndarray,
                      reglage_echappement: ReglageEchappement = None,
                      positions: dict = None, pas_balayage_deg: float = 10.0,
                      jeu_goupilles_um: float = 0.0) -> dict:
    """Analyse complète d'une géométrie en UN SEUL balayage elastica : les
    trois contributions à la marche diurne, séparées, plus leur somme par
    position et les grandeurs dérivées dont l'optimiseur a besoin.

    C'est le point d'entrée unique à utiliser partout (interface, optimiseur,
    étude de convergence) : le balayage elastica est la seule partie coûteuse,
    et il porte à lui seul l'orbite du centre de gravité, le jeu inter-spires
    ET la courbe d'énergie élastique.

    Clés du dictionnaire retourné :
        amplitudes_deg     grille d'amplitudes de toutes les courbes ci-dessous
        gravite            {position verticale: marche (s/j)} -- terme de CG
        isochronisme       marche propre du spiral (s/j), toutes positions
        echappement        marche due à l'échappement (s/j), toutes positions
        jeu_goupilles      marche due au jeu de raquette (s/j), toutes positions
        jeu_amplitude_basse  ce même retard à l'amplitude basse : c'est l'ÉCART
                           TOTAL entre un montage à piton fixe et un montage à
                           raquette, donc le coût du compromis entre les deux
        totales            {position: somme}, verticales ET horizontales
        dispersion_gravite dispersion entre positions verticales (s/j)
        ecart_compensation écart à la cible de compensation (s/j, signé)
        iso_amplitude_basse marche propre à l'amplitude basse (s/j)
        diag_orbite        diagnostics du balayage (jeu, orbite CG, énergie)
        diag_isochronisme  diagnostics du calcul de période non linéaire
    """
    reglage = reglage_echappement if reglage_echappement is not None else ReglageEchappement()
    amplitudes_deg = np.asarray(amplitudes_deg, dtype=float)

    marches_grav, diag_orbite = marche_positions_verticales(
        geom, r_min, freq_cible, inertie_balancier, amplitudes_deg,
        positions=positions, pas_balayage_deg=pas_balayage_deg, return_diagnostics=True,
    )

    # L'amplitude basse est évaluée EXPLICITEMENT, pas interpolée depuis la
    # grille d'affichage : rien ne garantit qu'elle y tombe (l'utilisateur
    # choisit librement [amp_min, amp_max]), et une extrapolation silencieuse
    # par bord de np.interp donnerait une compensation fausse sans le dire.
    amp_iso = np.append(amplitudes_deg, float(reglage.amplitude_basse_deg))
    iso_tout, diag_iso = marche_isochronisme_intrinseque(
        diag_orbite["angles_deg"], diag_orbite["energie_J"], inertie_balancier,
        amp_iso, reglage.amplitude_reference_deg,
    )
    marche_iso, iso_basse = iso_tout[:-1], float(iso_tout[-1])

    # Jeu de goupilles : quatrième contribution, elle aussi indépendante de la
    # position. Nulle si jeu_goupilles_um vaut 0, cas qui correspond EXACTEMENT
    # au piton fixe au point de comptage (cf. raquette.py).
    jeu_tout, diag_jeu = analyser_jeu_goupilles(
        geom, r_min, inertie_balancier, jeu_goupilles_um, amp_iso,
        reglage.amplitude_reference_deg)
    marche_jeu, jeu_basse = jeu_tout[:-1], float(jeu_tout[-1])

    # Le balayage elastica ne va que jusqu'à max(amplitudes_deg) : au-delà,
    # l'isochronisme serait écrêté silencieusement sur cette borne, et la
    # compensation calculée à une amplitude qui n'est pas celle demandée.
    amplitude_balayee = float(np.max(amplitudes_deg))
    diag_iso = dict(diag_iso)
    diag_iso["amplitude_basse_hors_plage"] = bool(
        reglage.amplitude_basse_deg > amplitude_balayee + 1e-9
        or reglage.amplitude_reference_deg > amplitude_balayee + 1e-9)
    if diag_iso["amplitude_basse_hors_plage"]:
        diag_iso["valide"] = False
        diag_iso["raison"] = (
            f"amplitude basse ({reglage.amplitude_basse_deg:.0f}°) ou de référence "
            f"({reglage.amplitude_reference_deg:.0f}°) au-delà de la plage analysée "
            f"({amplitude_balayee:.0f}°)")
        iso_basse = float("nan")
        marche_iso = np.full_like(marche_iso, np.nan)
    marche_ech = marche_echappement(amplitudes_deg, reglage)

    # Part du retard de jeu que le spiral doit compenser, selon le montage visé
    # (cf. ReglageEchappement.variante_visee) : 0 pour un piton fixe, 1 pour un
    # montage à raquette, 1/2 pour répartir l'erreur symétriquement entre les
    # deux. Les trois coïncident exactement quand le jeu est nul.
    fraction = reglage.fraction_jeu()
    ecart_compensation = (iso_basse - float(reglage.retard_s_par_j)
                          + fraction * jeu_basse)

    return {
        "amplitudes_deg": amplitudes_deg,
        "gravite": marches_grav,
        "isochronisme": marche_iso,
        "echappement": marche_ech,
        "jeu_goupilles": marche_jeu,
        "totales": composer_marches(marches_grav, marche_iso, marche_ech, marche_jeu),
        "dispersion_gravite": critere_dispersion_marche(marches_grav),
        "iso_amplitude_basse": iso_basse,
        "jeu_amplitude_basse": jeu_basse,
        "ecart_compensation": ecart_compensation,
        "fraction_jeu_visee": fraction,
        "reglage_echappement": reglage,
        "diag_orbite": diag_orbite,
        "diag_isochronisme": diag_iso,
        "diag_jeu_goupilles": diag_jeu,
    }


def critere_dispersion_marche(marches: dict) -> float:
    """Critère scalaire pour l'optimisation : écart max entre positions
    verticales, à chaque amplitude, maximisé sur les amplitudes (pire cas).
    C'est CE critère (et non un proxy géométrique) qui doit être minimisé
    par l'optimiseur pour répondre directement à l'objectif du cahier des
    charges (erreur de marche en positions verticales en fonction de
    l'amplitude)."""
    valeurs = np.array(list(marches.values()))  # (n_positions, n_amplitudes)
    ecart_par_amplitude = valeurs.max(axis=0) - valeurs.min(axis=0)
    return float(np.max(ecart_par_amplitude))
