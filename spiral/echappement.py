# echappement.py
"""
Erreur de marche introduite par l'ÉCHAPPEMENT (et non par le spiral).

Pourquoi c'est un module séparé, et pourquoi ça change la structure de
l'optimisation :

L'échappement à ancre suisse perturbe l'oscillateur par un couple parasite
quasi constant en module (frottement de repos/tirage, impulsion, dégagement).
Un couple parasite CONSTANT appliqué à un oscillateur d'amplitude Theta_0
produit, au premier ordre, une erreur de marche proportionnelle à 1/Theta_0 :
l'effet est d'autant plus grand que l'amplitude est faible. C'est le
« retard aux faibles amplitudes » bien connu, et c'est cette dépendance que
l'on modélise ici.

Deux propriétés déterminantes pour la conception du spiral :

  1. L'erreur d'échappement dépend de l'AMPLITUDE, pas de l'ORIENTATION de la
     montre : elle est présente à l'identique dans les positions horizontales
     (cadran haut / cadran bas) comme dans les positions verticales.
  2. Le terme de gravité du spiral (déplacement du centre de gravité,
     isochronisme.marche_positions_verticales) est au contraire NUL en
     position horizontale -- la gravité y est perpendiculaire au plan du
     spiral, son couple perturbateur autour de l'axe du balancier est nul.

Conséquence directe : le terme de gravité ne peut PAS compenser l'échappement.
La compensation ne peut venir que de l'ISOCHRONISME PROPRE du spiral (la
non-linéarité de son couple de rappel en grande déformation, qui fait varier
la période avec l'amplitude indépendamment de toute orientation) -- c'est
isochronisme.marche_isochronisme_intrinseque, calculé à partir de l'énergie
élastique U(alpha) du balayage. On demande donc au spiral de gagner
`retard_s_par_j` à l'amplitude basse pour annuler la perte de l'échappement :
c'est la seule composante du modèle qui agisse dans les quatre positions
verticales ET les deux horizontales à la fois.

Ce module ne modélise PAS un échappement particulier à partir de sa géométrie
(angles de levée, de repos, de tirage, couple moteur) : il prend la valeur du
retard comme une DONNÉE mesurée ou spécifiée par l'utilisateur, et n'en
modélise que la loi de variation avec l'amplitude.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class ReglageEchappement:
    """Cible de compensation de l'erreur d'échappement.

    retard_s_par_j : retard (en valeur absolue, s/j) induit par l'échappement
        à `amplitude_basse_deg` PAR RAPPORT à `amplitude_reference_deg`.
        C'est donc l'avance que le spiral doit produire de lui-même à
        l'amplitude basse pour que le total soit plat. 8 s/j à 200° est
        l'ordre de grandeur usuel d'un échappement à ancre suisse réglé.
    amplitude_basse_deg : amplitude basse de référence (fin de réserve de
        marche), là où l'erreur d'échappement est la plus pénalisante.
    amplitude_reference_deg : amplitude « pleine », à laquelle la marche est
        conventionnellement mise à zéro (le réglage moyen de la montre absorbe
        un décalage constant, seule la VARIATION avec l'amplitude compte).
    tolerance_s_par_j : largeur de la bande acceptée autour de la cible pour
        la contrainte d'optimisation. Une cible exacte (tolérance nulle) est
        inatteignable et rendrait toute la population infaisable ; une bande
        laisse NSGA-II arbitrer entre compensation et dispersion de gravité.
        Le défaut de 2 s/j est calé sur l'INCERTITUDE DE MAILLAGE mesurée sur
        cette grandeur (étude de convergence du 2026-09-24 : l'isochronisme
        propre à 200° passe de 20.8 à 18.4 s/j entre 1670 et 13360 éléments,
        soit 1 à 2.4 s/j d'écart aux maillages de travail courants).
        Descendre sous cette valeur sans raffiner le maillage revient à
        contraindre du bruit de discrétisation -- le bouton « Vérifier à
        haute résolution » le signale explicitement quand c'est le cas.
    actif : si False, aucune compensation n'est demandée à l'optimiseur (la
        marche d'échappement reste calculée et affichée, à titre informatif).
    variante_visee : sur quel MONTAGE la compensation est réglée, quand le
        spiral est à raquette. Un même spiral ne peut pas être parfait dans les
        deux, puisque le jeu de goupilles n'existe que dans l'un (cf.
        raquette.py) et que l'écart entre les deux est une pente, pas un
        décalage constant -- la raquette elle-même ne peut donc pas le
        rattraper.
          "milieu"  : réglé à mi-chemin, chaque montage est faux de la moitié ;
          "raquette": réglé pour le montage à raquette, jeu compris ;
          "piton"   : réglé pour le piton fixe au point de comptage, comme si
                      le jeu n'existait pas.
        Sans effet si le spiral n'a pas de goupilles ou si le jeu est nul :
        les trois réglages coïncident alors exactement.
    """
    retard_s_par_j: float = 8.0
    amplitude_basse_deg: float = 200.0
    amplitude_reference_deg: float = 300.0
    tolerance_s_par_j: float = 2.0
    actif: bool = True
    variante_visee: str = "milieu"   # "milieu" | "raquette" | "piton"

    def fraction_jeu(self) -> float:
        """Part du retard de jeu que le spiral doit compenser de lui-même."""
        return {"piton": 0.0, "milieu": 0.5, "raquette": 1.0}.get(self.variante_visee, 0.5)


def marche_echappement(amplitudes_deg, reglage: ReglageEchappement) -> np.ndarray:
    """Marche (s/j) due au seul échappement, en fonction de l'amplitude,
    rapportée à `amplitude_reference_deg` (donc nulle à cette amplitude, et
    égale à -retard_s_par_j à `amplitude_basse_deg`).

    Loi en 1/Theta_0 : erreur de marche d'un oscillateur harmonique soumis à
    un couple parasite constant, au premier ordre. C'est le seul degré de
    liberté du modèle ; deux valeurs (le retard et l'amplitude à laquelle il
    est spécifié) suffisent à le caler entièrement, sans avoir à décrire la
    géométrie de l'échappement.

    Identique dans TOUTES les positions (horizontales comprises) : aucune
    dépendance à l'orientation, contrairement au terme de gravité.
    """
    amplitudes_deg = np.asarray(amplitudes_deg, dtype=float)
    theta_b = float(reglage.amplitude_basse_deg)
    theta_ref = float(reglage.amplitude_reference_deg)
    if abs(1.0 / theta_b - 1.0 / theta_ref) < 1e-12:
        return np.zeros_like(amplitudes_deg)
    K = -float(reglage.retard_s_par_j) / (1.0 / theta_b - 1.0 / theta_ref)
    return K / amplitudes_deg - K / theta_ref
