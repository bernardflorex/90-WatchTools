# parametrage.py
"""
Génération paramétrique d'un spiral, en remplacement des "secteurs longs" de
spiropt_v1 (5 noeuds corps + 3 noeuds Grossmann interpolés, un seul scalaire
de forme "cambrure" par extrémité) par :

  - une courbe de GROSSMANN (intérieure, côté collet) et une courbe de
    PHILLIPS (extérieure, côté piton) -- terminologie historique correcte
    (voir courbes.py) -- chacune paramétrée par :
      * `*_angle_deg`  : étendue ANGULAIRE de la zone de raccordement.
      * `*_n_pas`      : nombre d'ÉLÉMENTS DISCRETS réservés à cette zone
                         dans le maillage physique (densité locale, utile
                         car c'est là que la courbure est la plus forte).
      * `*_n_points`   : nombre de points de contrôle de Bézier de la
                         courbe elle-même (>= 4). Remplace le scalaire
                         unique "cambrure" par une forme réglable, avec
                         inflexions possibles si souhaité -- justifié par la
                         théorie (Phillips ne prescrit aucune forme précise,
                         seulement une condition sur son centre de gravité,
                         voir courbes.condition_phillips).
      * `*_offsets`    : les n_points-3 variables de forme libres
                         elles-mêmes (écarts en mm par rapport au raccord
                         rectiligne neutre), longueur = n_points - 3.
  - un profil d'épaisseur e(s) et un écart radial delta_r(s) libres, chacun
    porté par un petit nombre de points de contrôle (PCHIP, préserve la
    monotonie locale, pas de dépassement/oscillation contrairement à une
    spline cubique classique) échantillonnés sur l'abscisse curviligne
    normalisée du spiral déjà construit.

Principe de construction (généralisation du panneau 2 de spiropt_v1 aux deux
extrémités, cf. courbes.py pour le détail Bézier) : une ligne de référence à
pas constant `pas_theorique = (r_max - r_min) / (n_tours_fictifs * 2*pi)`
sert UNIQUEMENT à définir la pente imposée (continuité C1) à la jonction
entre le corps et chaque courbe. Le corps est par construction sur cette
ligne (donc sa jonction est toujours exactement à la bonne pente, sans
écrêtage). Côté piton, r_max est une contrainte dure indépendante de cette
ligne : c'est le désaccord entre "où la ligne de référence arriverait" et "où
il faut vraiment arriver" que la forme de la courbe de Phillips doit
absorber -- sa vraie raison d'être historique. Côté collet, r_min est
également une contrainte dure, mais elle est PAR CONSTRUCTION sur la ligne de
référence (qui part de r_min à theta=0) : la courbe de Grossmann n'a donc pas
de désaccord à absorber, ses points libres ne pilotent que la forme
(tangente au collet), pas une fermeture géométrique. Asymétrie voulue,
physiquement cohérente (le piton est une contrainte d'encombrement dure, le
collet un choix de forme) -- cf. spiral_lab/courbes.py pour les sources.
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.interpolate import PchipInterpolator

from .modele import GeometrieSpiral
from .courbes import construire_zone, condition_phillips


def epaisseur_theorique(r_min: float, r_max: float, n_tours_reels: float, hauteur_h: float,
                         freq_cible: float, inertie_balancier: float) -> float:
    """Estimation analytique grossière (poutre à épaisseur constante,
    longueur = rayon moyen x nombre de tours) de l'épaisseur nécessaire pour
    atteindre freq_cible -- sert uniquement à initialiser `epaisseur_noeuds`
    à une valeur raisonnable (au lieu d'une valeur arbitraire) avant sondage/
    optimisation. Reprend calculer_epaisseur_theorique de spiropt_v1/physics.py."""
    from .constantes import E_SILICIUM
    omega = 2.0 * np.pi * freq_cible
    k_cible = (omega ** 2) * inertie_balancier
    r_moyen = (r_min + r_max) / 2.0
    longueur_m = r_moyen * n_tours_reels * 2.0 * np.pi * 1e-3
    hauteur_m = hauteur_h * 1e-3
    e_m = ((12.0 * longueur_m * k_cible) / (E_SILICIUM * hauteur_m)) ** (1.0 / 3.0)
    return e_m * 1000.0


@dataclass
class ParametresSpiral:
    # --- Contexte fixe (encombrement, cible chronométrique) ---
    r_min: float
    r_max: float
    n_tours_reels: float
    n_tours_fictifs: float
    hauteur_h: float
    freq_cible: float
    inertie_balancier: float
    n_elements_corps: int = 1400

    # --- Courbe de Grossmann (intérieure, côté collet) ---
    grossmann_angle_deg: float = 90.0
    grossmann_n_pas: int = 120
    grossmann_n_points: int = 4
    grossmann_offsets: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # --- Courbe de Phillips (extérieure, côté piton) ---
    phillips_angle_deg: float = 120.0
    phillips_n_pas: int = 150
    phillips_n_points: int = 4
    phillips_offsets: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # --- Profils libres (valeurs aux points de contrôle, PCHIP sur s/L in [0,1]) ---
    epaisseur_noeuds: np.ndarray = field(default_factory=lambda: np.full(20, 0.030))
    delta_r_noeuds: np.ndarray = field(default_factory=lambda: np.zeros(12))

    @property
    def n_elements(self) -> int:
        return self.n_elements_corps + self.grossmann_n_pas + self.phillips_n_pas


def rayons_jonction(p: ParametresSpiral):
    """Rayons de jonction corps<->Grossmann et corps<->Phillips sur la ligne
    de référence à pas constant, et le pas lui-même -- calcul partagé entre
    construire_geometrie et optimisation.construire_bornes, pour que les
    bornes des offsets Grossmann/Phillips puissent être dimensionnées en
    proportion de la largeur radiale RÉELLE de chaque zone (qui dépend de
    l'angle et du pas donc du calibre) plutôt qu'en une valeur fixe en mm
    -- une largeur fixe s'est révélée beaucoup trop grande pour un calibre
    à pas serré (diagnostic 2026-09-22 : ~55% des tirages aléatoires
    devenaient non-monotones avec spread_grossmann_mm=0.030 sur une zone de
    ~39 µm de large seulement)."""
    theta_total = p.n_tours_reels * 2.0 * np.pi
    pas_theorique = (p.r_max - p.r_min) / (p.n_tours_fictifs * 2.0 * np.pi)
    theta_gross = np.radians(p.grossmann_angle_deg)
    theta_phil = np.radians(p.phillips_angle_deg)
    theta_corps = max(theta_total - theta_gross - theta_phil, 1e-6)
    r_body_start = p.r_min + pas_theorique * theta_gross
    r_body_end = r_body_start + pas_theorique * theta_corps
    return pas_theorique, r_body_start, r_body_end


def construire_geometrie(p: ParametresSpiral, return_diagnostics: bool = False):
    """Construit la géométrie complète : squelette collet -> Grossmann ->
    corps -> Phillips -> piton, puis applique épaisseur et écart radial
    libres. Retourne une GeometrieSpiral (et un dict de diagnostics si
    demandé -- validité géométrique, monotonie des deux courbes, repère de
    plausibilité Phillips -- à utiliser comme pénalité/contrainte côté
    optimiseur plutôt que de planter)."""
    if len(p.grossmann_offsets) != p.grossmann_n_points - 3:
        raise ValueError(f"grossmann_n_points={p.grossmann_n_points} implique "
                          f"{p.grossmann_n_points - 3} offsets, reçu {len(p.grossmann_offsets)}")
    if len(p.phillips_offsets) != p.phillips_n_points - 3:
        raise ValueError(f"phillips_n_points={p.phillips_n_points} implique "
                          f"{p.phillips_n_points - 3} offsets, reçu {len(p.phillips_offsets)}")

    theta_total = p.n_tours_reels * 2.0 * np.pi
    theta_gross = np.radians(p.grossmann_angle_deg)
    theta_phil = np.radians(p.phillips_angle_deg)
    theta_corps = theta_total - theta_gross - theta_phil
    geometrie_valide = theta_corps > 1e-6

    pas_theorique, r_body_start, r_body_end = rayons_jonction(p)
    theta_corps = max(theta_corps, 1e-6)

    # --- Zone de Grossmann : ancree en u=1 (jonction corps), libre en u=0 (collet, r_min) ---
    zone_g = construire_zone(
        p.grossmann_n_points, r_ancre=r_body_start, pente_ancre_du=pas_theorique * theta_gross,
        r_libre=p.r_min, offsets_libres=p.grossmann_offsets, ancre_a_u1=True,
    )
    u_g = np.linspace(0.0, 1.0, p.grossmann_n_pas, endpoint=False)
    theta_g_arr = u_g * theta_gross
    r_g = zone_g.r(u_g)

    # --- Zone corps : pas constant, exactement sur la ligne de reference ---
    theta_b_arr = theta_gross + np.linspace(0.0, theta_corps, p.n_elements_corps, endpoint=False)
    r_b = p.r_min + pas_theorique * theta_b_arr

    # --- Zone de Phillips : ancree en u=0 (jonction corps), libre en u=1 (piton, r_max) ---
    zone_ph = construire_zone(
        p.phillips_n_points, r_ancre=r_body_end, pente_ancre_du=pas_theorique * theta_phil,
        r_libre=p.r_max, offsets_libres=p.phillips_offsets, ancre_a_u1=False,
    )
    u_ph = np.linspace(0.0, 1.0, p.phillips_n_pas)
    theta_ph_arr = theta_gross + theta_corps + u_ph * theta_phil
    r_ph = zone_ph.r(u_ph)

    theta_final = np.concatenate([theta_g_arr, theta_b_arr, theta_ph_arr])
    r_ref = np.concatenate([r_g, r_b, r_ph])

    # --- Ecart radial libre, PCHIP sur l'indice normalise (proxy de s/L avant deformation) ---
    u_idx = np.linspace(0.0, 1.0, len(theta_final))
    if len(p.delta_r_noeuds) >= 2:
        u_noeuds_r = np.linspace(0.0, 1.0, len(p.delta_r_noeuds))
        delta_r = PchipInterpolator(u_noeuds_r, p.delta_r_noeuds)(u_idx)
    else:
        delta_r = np.zeros_like(theta_final)
    r_final = r_ref + delta_r

    x = r_final * np.cos(theta_final)
    y = r_final * np.sin(theta_final)

    # --- Epaisseur, PCHIP sur l'abscisse curviligne normalisee reelle ---
    ds = np.hypot(np.diff(x, prepend=x[0]), np.diff(y, prepend=y[0]))
    ds[0] = 0.0
    s_cum = np.cumsum(ds)
    s_norm = s_cum / s_cum[-1]
    u_noeuds_e = np.linspace(0.0, 1.0, len(p.epaisseur_noeuds))
    epaisseur = PchipInterpolator(u_noeuds_e, p.epaisseur_noeuds)(s_norm)

    # --- Mise a l'echelle analytique de l'epaisseur pour tomber pile sur freq_cible ---
    from .constantes import E_SILICIUM
    souplesse_int = np.sum((ds * 1e-3) / np.clip(epaisseur * 1e-3, 1e-6, None) ** 3)
    couple = (E_SILICIUM * (p.hauteur_h * 1e-3)) / (12.0 * souplesse_int)
    freq_actuelle = (1.0 / (2.0 * np.pi)) * np.sqrt(couple / p.inertie_balancier)
    beta = (p.freq_cible / freq_actuelle) ** (2.0 / 3.0)
    epaisseur = epaisseur * beta

    geom = GeometrieSpiral(x=x, y=y, epaisseur=epaisseur, hauteur=p.hauteur_h)

    if not return_diagnostics:
        return geom

    n_g = len(theta_g_arr)
    n_b = len(theta_b_arr)
    ratio_phillips = condition_phillips(
        x[n_g + n_b:], y[n_g + n_b:], ds[n_g + n_b:],
        x_centre_spiral=0.0, y_centre_spiral=0.0,
        x_decrochement=x[-1], y_decrochement=y[-1],
    )
    diagnostics = {
        "geometrie_valide": geometrie_valide,
        "theta_corps_deg": float(np.degrees(theta_corps)),
        "r_body_start_mm": float(r_body_start),
        "r_body_end_mm": float(r_body_end),
        "grossmann_monotone": zone_g.monotone(),
        "phillips_monotone": zone_ph.monotone(),
        "phillips_ratio_plausibilite": ratio_phillips,
        "freq_avant_scaling_hz": float(freq_actuelle),
        "beta_scaling": float(beta),
    }
    return geom, diagnostics
