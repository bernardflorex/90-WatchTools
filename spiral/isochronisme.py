# isochronisme.py
"""
Marche diurne induite par le déplacement du centre de gravité (CG) du spiral
sous l'effet de la gravité, dans les positions verticales de la montre,
via la formule d'Airy (perturbation au premier ordre d'un oscillateur
harmonique par un couple perturbateur périodique).

x_cg(phi)/y_cg(phi) viennent du VRAI solveur elastica (physique.balayage_angulaire),
pas d'une approximation à courbure uniforme (invalide dès que l'épaisseur
varie le long de la lame, cf. note historique du projet : écart mesuré ~50x
sur la dispersion de l'orbite CG entre les deux approches).
"""
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline

from .constantes import G_TERRE, RHO_SILICIUM
from .modele import GeometrieSpiral
from .physique import balayage_angulaire

POSITIONS_VERTICALES_DEFAUT = {
    "3H": np.radians(180.0),
    "6H": np.radians(90.0),
    "9H": np.radians(0.0),
    "12H": np.radians(270.0),
}


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
    angles_deg, x_cg_mm, y_cg_mm, jeu_mm = balayage_angulaire(geom, r_min, amplitude_max, pas_balayage_deg)
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
        "angle_jeu_min_deg": float(angles_deg[np.argmin(jeu_mm)]), "jeu_min_mm": float(np.min(jeu_mm)),
    }
    return marches, diagnostics


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
