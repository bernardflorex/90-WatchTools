# collision.py
"""
Jeu inter-spires minimal, généralisé au cas où le rayon n'est plus une
fonction affine de l'angle (pas constant) : on ne suppose plus que le "tour
suivant" est à un pas fixe, on va le chercher par interpolation sur la
géométrie réellement déformée.

Hypothèse conservée (doit être respectée par les bornes de la paramétrisation,
voir parametrage.py) : r(theta) reste une fonction univoque de theta, i.e. le
spiral ne "replie" jamais sur lui-même localement — seul le PAS peut varier,
pas le sens de parcours. C'est une hypothèse physique de toute façon
nécessaire pour qu'un spiral plat soit valide.
"""
import numpy as np


def jeu_min_spires(X_def: np.ndarray, Y_def: np.ndarray, epaisseurs: np.ndarray) -> float:
    """Jeu minimal (mm) entre spires adjacentes, sur une géométrie déformée
    (fibre neutre). Retourne +inf si le spiral fait moins d'un tour complet
    (pas de spire adjacente à comparer)."""
    R_m = np.hypot(X_def, Y_def)
    T_m = np.unwrap(np.arctan2(Y_def, X_def))

    idx_sort = np.argsort(T_m)
    T_s, R_s, e_s = T_m[idx_sort], R_m[idx_sort], epaisseurs[idx_sort]

    min_space = float("inf")
    for i in range(len(T_s)):
        t_cible = T_s[i] + 2.0 * np.pi
        if t_cible <= T_s[-1]:
            r_suivant = np.interp(t_cible, T_s, R_s)
            e_suivant = np.interp(t_cible, T_s, e_s)
            clearance = r_suivant - R_s[i] - (e_s[i] + e_suivant) / 2.0
            if clearance < min_space:
                min_space = clearance
    return min_space
