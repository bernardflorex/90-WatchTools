# raquette.py
"""
Erreur de marche due au JEU ENTRE LES GOUPILLES de la raquette.

Mécanisme. La lame passe entre deux goupilles avec un jeu libre. Tant qu'elle
est plaquée contre l'une d'elles, le spiral est effectivement encastré AUX
GOUPILLES : la longueur active va du collet aux goupilles, et la raideur est
`k_goupille`. Dès que la lame décolle et traverse le jeu, plus rien ne la tient
là : le spiral se comporte alors comme s'il était encastré AU PITON, donc plus
long et plus souple, `k_libre < k_goupille`. La montre retarde pendant ces
traversées.

Ce qui commande l'effet, ce sont les deux grandeurs que l'expérience
horlogère désigne depuis toujours :

  - le JEU : plus les goupilles sont écartées, plus la lame met de temps à
    traverser, donc plus la phase souple dure longtemps ;
  - la LONGUEUR APRÈS LES GOUPILLES : elle fixe à la fois l'ampleur du saut de
    raideur (plus il reste de lame morte, plus la perte de raideur est grande)
    ET la vitesse de traversée, puisque c'est la flexion de ce tronçon libre
    qui déplace latéralement la lame au droit des goupilles. Mesuré sur un
    calibre de référence : l'excursion latérale varie comme le CARRÉ de cette
    longueur (0.9 µm pour 124 µm de lame morte, 7.6 µm pour 371 µm, 30 µm pour
    742 µm, sur ±340° de balancier).

Ces deux effets vont dans le même sens : jeu ou longueur morte plus grands
=> retard plus marqué, et d'autant plus marqué que l'amplitude est faible.

Méthode. La raideur est constante par morceaux en fonction de l'angle du
balancier : `k_libre` dans la bande morte |phi| < phi_c, `k_goupille` au-delà.
Ce n'est pas une non-linéarité lisse, donc ni la formule d'Airy (perturbation
d'un oscillateur harmonique) ni l'intégration exacte de la période (qui exige
une énergie régulière et monotone) ne conviennent. L'outil adapté est la
FONCTION DESCRIPTIVE de la zone morte (équilibrage harmonique), qui donne la
raideur équivalente vue par une oscillation sinusoïdale d'amplitude Theta_0.

Cas limite important : jeu = 0 donne phi_c = 0, donc un facteur de zone morte
égal à 1, donc k_eq = k_goupille EXACTEMENT, à toute amplitude, donc une marche
identiquement nulle. C'est rigoureusement la configuration « piton fixe au
point de comptage ». Le spiral libre et le spiral à raquette ne sont donc pas
deux modèles à réconcilier : ce sont les deux bouts d'un seul paramètre
continu, et le passage de l'un à l'autre se fait sans discontinuité.
"""
from dataclasses import replace
import numpy as np

from .modele import GeometrieSpiral
from .physique import frequence_propre, resoudre_elastica


def facteur_zone_morte(epsilon) -> np.ndarray:
    """Fonction descriptive d'une zone morte de demi-largeur relative epsilon
    (= phi_c / Theta_0) : fraction de la raideur additionnelle effectivement
    ressentie par une oscillation sinusoïdale d'amplitude Theta_0.

        N(0) = 1  -> bande morte négligeable devant l'amplitude : tout se passe
                     comme si la lame était toujours plaquée (pas de jeu).
        N(1) = 0  -> l'amplitude ne dépasse pas la bande morte : la lame ne
                     touche jamais, la raquette n'a aucune autorité.
    """
    e = np.clip(np.asarray(epsilon, dtype=float), 0.0, 1.0)
    return 1.0 - (2.0 / np.pi) * (np.arcsin(e) + e * np.sqrt(np.maximum(1.0 - e * e, 0.0)))


def excursion_laterale_um_par_rad(geom: GeometrieSpiral, r_min: float,
                                   amplitude_max_deg: float = 340.0) -> float:
    """Vitesse de déplacement latéral de la lame AU DROIT DES GOUPILLES, en µm
    par radian de rotation du balancier, dans l'état libre (encastré au piton
    seul) -- c'est cet état qui décide du moment où le contact se produit.

    Le point de mesure est FIXE DANS L'ESPACE (les goupilles ne bougent pas
    pendant la marche, c'est la lame qui coulisse entre elles) : on cherche
    donc le point de la lame déformée le plus proche du centre des goupilles,
    et on projette l'écart sur la normale radiale. Suivre un point MATÉRIEL
    donnerait une autre grandeur, sans rapport avec le contact.

    Retourne 0.0 si la géométrie n'a pas de portion morte (spiral libre).
    """
    if not geom.a_partie_morte:
        return 0.0
    i = geom.n_actif - 1
    Px, Py = float(geom.x[i]), float(geom.y[i])
    norme = float(np.hypot(Px, Py))
    if norme < 1e-12:
        return 0.0
    ux, uy = Px / norme, Py / norme

    # Recherche limitée à la dernière spire : une spire intérieure pourrait
    # sinon être plus proche que la lame qui passe réellement aux goupilles.
    j0 = int(np.searchsorted(geom.theta, geom.theta[i] - 2.0 * np.pi))
    geom_libre = replace(geom, n_actif=geom.n)

    angles = np.array([-amplitude_max_deg, 0.0, amplitude_max_deg])
    ecarts = np.empty(len(angles))
    for n, a in enumerate(angles):
        Xe, Ye, Xi, Yi, _, _ = resoudre_elastica(geom_libre, float(a), r_min)
        xm, ym = (Xe[j0:] + Xi[j0:]) / 2.0, (Ye[j0:] + Yi[j0:]) / 2.0
        k = int(np.argmin((xm - Px) ** 2 + (ym - Py) ** 2))
        ecarts[n] = ((xm[k] - Px) * ux + (ym[k] - Py) * uy) * 1000.0
    pente = float(np.polyfit(np.radians(angles), ecarts, 1)[0])
    return abs(pente)


def analyser_jeu_goupilles(geom: GeometrieSpiral, r_min: float, inertie_balancier: float,
                            jeu_um: float, amplitudes_deg: np.ndarray,
                            amplitude_reference_deg: float):
    """Marche (s/j) due au seul jeu de goupilles, rapportée à l'amplitude de
    référence, plus un dict de diagnostics.

    Indépendante de la position de la montre, comme l'erreur d'échappement et
    l'isochronisme propre : elle s'ajoute donc à ce que la conception du spiral
    doit compenser, et non à la dispersion entre positions verticales.
    """
    amplitudes_deg = np.atleast_1d(np.asarray(amplitudes_deg, dtype=float))
    nul = np.zeros_like(amplitudes_deg)
    if not geom.a_partie_morte or jeu_um <= 0.0:
        return nul, {"actif": False, "raison": "piton fixe au point de comptage (jeu nul)"}

    f_goupille = frequence_propre(geom, inertie_balancier)
    f_libre = frequence_propre(replace(geom, n_actif=geom.n), inertie_balancier)
    excursion = excursion_laterale_um_par_rad(geom, r_min, float(np.max(amplitudes_deg)))

    if excursion <= 1e-12:
        # La lame ne bouge pas latéralement : elle ne touchera jamais une
        # goupille, quel que soit le jeu. La raquette est sans effet et le
        # spiral se comporte comme encastré au piton.
        phi_c = np.inf
    else:
        phi_c = (jeu_um / 2.0) / excursion  # rad de balancier

    rapport_k = (f_goupille / f_libre) ** 2

    def raideur_relative(theta_deg):
        theta = np.radians(np.asarray(theta_deg, dtype=float))
        with np.errstate(divide="ignore", invalid="ignore"):
            eps = np.where(theta > 0, phi_c / np.maximum(theta, 1e-12), 1.0)
        return 1.0 + (rapport_k - 1.0) * facteur_zone_morte(eps)

    k_ref = float(raideur_relative(amplitude_reference_deg))
    marches = 86400.0 * (np.sqrt(raideur_relative(amplitudes_deg) / k_ref) - 1.0)

    diagnostics = {
        "actif": True,
        "excursion_um_par_rad": float(excursion),
        "amplitude_contact_deg": float(np.degrees(phi_c)),
        "ecart_raideur_s_par_j": float(86400.0 * (f_goupille / f_libre - 1.0)),
        "autorite": float(facteur_zone_morte(phi_c / np.radians(float(np.max(amplitudes_deg))))),
    }
    return marches, diagnostics
