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
import numpy as np
from scipy.optimize import minimize as scipy_minimize

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
    n = len(phi0)
    phi_def = np.zeros(n)
    phi_def[-1] = phi_piton
    for i in range(n - 2, -1, -1):
        M = M0 + Fx * (Ym0[i] - y_piton) - Fy * (Xm0[i] - x_piton)
        phi_def[i] = phi_def[i + 1] - (phi0[i + 1] - phi0[i]) - (M / raideur_flexion[i]) * (ds[i] * 1e-3)

    X_def = np.zeros(n)
    Y_def = np.zeros(n)
    X_def[-1] = x_piton
    Y_def[-1] = y_piton
    for i in range(n - 2, -1, -1):
        X_def[i] = X_def[i + 1] - ds[i] * np.cos(phi_def[i])
        Y_def[i] = Y_def[i + 1] - ds[i] * np.sin(phi_def[i])

    cout = (X_def[0] - x_cible) ** 2 + (Y_def[0] - y_cible) ** 2 + ((phi_def[0] - phi_cible) * 0.1) ** 2
    return cout, X_def, Y_def, phi_def


def resoudre_elastica(geom: GeometrieSpiral, angle_deg: float, r_min: float,
                       return_diagnostics: bool = False):
    """Résout l'elastica non-linéaire pour une rotation `angle_deg` du
    balancier (le collet, index 0, tourne de cet angle ; le piton, index -1,
    reste fixe). r_min sert à replacer le point cible du collet sur le cercle
    de rayon r_min (convention identique à l'origine : le collet est solidaire
    de l'axe du balancier).

    Retourne (Xe, Ye, Xi, Yi, x_cg, y_cg) : bords extérieur/intérieur déformés
    et centre de gravité (mm), plus un dict de diagnostics si demandé.
    """
    Xm0, Ym0, epaisseurs = geom.x, geom.y, geom.epaisseur
    alpha = np.radians(angle_deg)
    dX0, dY0 = np.gradient(Xm0), np.gradient(Ym0)
    ds = np.sqrt(dX0 ** 2 + dY0 ** 2)
    masse_elementaire = RHO_SILICIUM * (geom.hauteur * 1e-3) * (epaisseurs * 1e-3) * (ds * 1e-3)

    phi0 = np.unwrap(np.arctan2(dY0, dX0))
    I_local = ((geom.hauteur * 1e-3) * (epaisseurs * 1e-3) ** 3) / 12.0
    raideur_flexion = E_SILICIUM * I_local

    theta_init = np.arctan2(Ym0[0], Xm0[0])
    x_cible = r_min * np.cos(theta_init + alpha)
    y_cible = r_min * np.sin(theta_init + alpha)
    phi_cible = phi0[0] + alpha
    x_piton, y_piton, phi_piton = Xm0[-1], Ym0[-1], phi0[-1]

    def integrer_tir(parametres):
        Fx, Fy, M0 = parametres[0] * 1e-3, parametres[1] * 1e-3, parametres[2] * 1e-6
        return _integrer_tir_numba(Fx, Fy, M0, phi0, Xm0, Ym0, ds, raideur_flexion,
                                    x_piton, y_piton, phi_piton, x_cible, y_cible, phi_cible)

    souplesse = np.sum((ds * 1e-3) / raideur_flexion)
    p0 = [0.0, 0.0, (alpha / souplesse) * 1e6 if abs(alpha) > 1e-5 else 0.0]

    res_tir = scipy_minimize(lambda ps: integrer_tir(ps)[0], p0, method="BFGS",
                              options={"gtol": 1e-6, "maxiter": 60})
    cout_final, X_f, Y_f, phi_f = integrer_tir(res_tir.x)

    # Garde-fou anti-divergence (diagnostic 2026-09-25, étude de convergence
    # de maillage) : le tir BFGS peut, pour un angle/maillage particulier,
    # converger vers un minimum local aberrant (résidu de plusieurs dizaines
    # de µm, contre <0.3 µm partout ailleurs) -- observé concrètement à un
    # angle donné pour un seul n de maillage sur ~35 testés, mais suffisant
    # pour corrompre la spline de l'orbite du CG et fausser silencieusement
    # la dispersion de marche calculée à partir d'elle (269 s/j au lieu de
    # ~13 s/j sur le cas ayant révélé le problème). Un second essai, parti
    # d'un point de départ différent (M0=0 au lieu de l'estimation linéaire),
    # suffit à s'en sortir quand ça se produit ; ne coûte rien dans le cas
    # normal (déclenché seulement si le premier essai a clairement divergé).
    erreur_um = float(np.hypot(X_f[0] - x_cible, Y_f[0] - y_cible) * 1000.0)
    if erreur_um > 1.0:
        res_tir_alt = scipy_minimize(lambda ps: integrer_tir(ps)[0], [0.0, 0.0, 0.0],
                                      method="BFGS", options={"gtol": 1e-6, "maxiter": 150})
        cout_alt, X_alt, Y_alt, phi_alt = integrer_tir(res_tir_alt.x)
        if cout_alt < cout_final:
            res_tir, cout_final, X_f, Y_f, phi_f = res_tir_alt, cout_alt, X_alt, Y_alt, phi_alt

    dX_d, dY_d = np.gradient(X_f), np.gradient(Y_f)
    norme_d = np.sqrt(dX_d ** 2 + dY_d ** 2)
    nX, nY = -dY_d / norme_d, dX_d / norme_d

    x_cg = np.sum(X_f * masse_elementaire) / np.sum(masse_elementaire)
    y_cg = np.sum(Y_f * masse_elementaire) / np.sum(masse_elementaire)

    Xe = X_f + (epaisseurs / 2.0) * nX
    Ye = Y_f + (epaisseurs / 2.0) * nY
    Xi = X_f - (epaisseurs / 2.0) * nX
    Yi = Y_f - (epaisseurs / 2.0) * nY

    if not return_diagnostics:
        return Xe, Ye, Xi, Yi, x_cg, y_cg

    erreur_position_um = float(np.hypot(X_f[0] - x_cible, Y_f[0] - y_cible) * 1000.0)
    erreur_angle_deg = float(np.degrees(abs(((phi_f[0] - phi_cible) + np.pi) % (2 * np.pi) - np.pi)))
    diagnostics = {
        "angle_deg": float(angle_deg), "success": bool(res_tir.success),
        "message": str(res_tir.message), "n_iter": int(res_tir.nit),
        "cout_final": float(cout_final), "erreur_position_um": erreur_position_um,
        "erreur_angle_deg": erreur_angle_deg,
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

    Retourne (angles_deg, x_cg_mm, y_cg_mm, jeu_mm) -- jeu_mm est le jeu
    inter-spires (mm) à CHAQUE angle balayé, permettant de retrouver ensuite
    l'angle de jeu minimal exact via angles_deg[np.argmin(jeu_mm)].

    Nombre de points ajusté pour être impair (theta=0 exactement au centre,
    requis par les métriques de symétrie)."""
    n_side = max(1, round(amplitude_max_deg / pas_angulaire_deg))
    angles_deg = np.linspace(-amplitude_max_deg, amplitude_max_deg, 2 * n_side + 1)
    x_cg = np.zeros_like(angles_deg)
    y_cg = np.zeros_like(angles_deg)
    jeu_mm = np.zeros_like(angles_deg)
    for i, ang in enumerate(angles_deg):
        Xe, Ye, Xi, Yi, x_c, y_c = resoudre_elastica(geom, ang, r_min)
        x_cg[i] = x_c
        y_cg[i] = y_c
        jeu_mm[i] = jeu_min_spires((Xe + Xi) / 2.0, (Ye + Yi) / 2.0, geom.epaisseur)
    return angles_deg, x_cg, y_cg, jeu_mm


def orbite_cg(geom: GeometrieSpiral, r_min: float, amplitude_max_deg: float,
              pas_angulaire_deg: float = 10.0):
    """Compatibilité : identique à balayage_angulaire sans le jeu, pour les
    appelants qui n'ont besoin que de l'orbite du CG."""
    angles_deg, x_cg, y_cg, _ = balayage_angulaire(geom, r_min, amplitude_max_deg, pas_angulaire_deg)
    return angles_deg, x_cg, y_cg


def frequence_propre(geom: GeometrieSpiral, inertie_balancier: float) -> float:
    """Fréquence propre (Hz) de l'oscillateur balancier-spiral pour CETTE
    géométrie, quelle que soit son origine (DXF importé, génération
    paramétrique, solution optimisée) -- même formule (poutre linéaire,
    couple = E.h / (12.souplesse_int)) que la mise à l'échelle automatique de
    parametrage.construire_geometrie, mais réutilisable après coup sur une
    géométrie déjà figée, pour afficher où l'oscillateur atterrit réellement
    par rapport à la fréquence cible visée."""
    ds_m = np.gradient(geom.s) * 1e-3
    e_m = np.clip(geom.epaisseur * 1e-3, 1e-9, None)
    souplesse_int = np.sum(ds_m / e_m ** 3)
    couple = (E_SILICIUM * (geom.hauteur * 1e-3)) / (12.0 * souplesse_int)
    return (1.0 / (2.0 * np.pi)) * np.sqrt(couple / inertie_balancier)
