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
from dataclasses import dataclass, field, replace
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .modele import GeometrieSpiral
from .courbes import construire_zone, construire_zone_s, condition_phillips


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
    # Orientation du spiral dans le mouvement (voir SpecMouvement) : angle de
    # la ligne d'échappement par rapport à l'axe de tige, et sens de
    # développement. Ils ne changent RIEN à la physique interne du spiral --
    # seulement son orientation par rapport aux positions de contrôle, dont
    # dépend entièrement le critère DVH.
    angle_echappement_deg: float = -6.2
    sens_trigonometrique: bool = True
    collet_cote_ancre: bool = True
    # Wafer : plan cristallin et orientation du dessin sur ce plan. Seul le
    # (111) est isotrope dans le plan ; voir constantes.module_young_silicium.
    plan_cristallin: str = "111"
    orientation_cristal_deg: float = 0.0

    @property
    def angle_collet_deg(self) -> float:
        """Direction de l'attache au collet dans le repère du tracé, en degrés.

        Deux informations s'y combinent, et il faut les deux : la ligne des
        centres est une DROITE, son angle ne la définit qu'à 180° près.

          1. `angle_echappement_deg` : l'angle qu'elle fait avec l'axe de tige.
             La tige sort à 3H, soit à gauche (180°) en vue de fond avec 12H en
             haut, d'où le premier terme.
          2. `collet_cote_ancre` : lequel de ses deux bouts porte le collet.
             Vrai (défaut) = celui qui pointe vers l'ancre, donc à l'OPPOSÉ de
             la tige, ce qui ramène le collet à angle_echappement_deg lui-même
             (−6.2° pour la valeur par défaut). Faux = l'autre bout, à 173.8°.

        Correction du 2026-09-29 : seul le terme 1 était pris en compte, et le
        spiral partait de 173.8° au lieu de −6.2°, soit un demi-tour d'écart
        avec le mouvement réel. Une rotation ne change pas le MODULE de la
        sensibilité à la gravité, mais elle en inverse la phase de 180°, donc
        le signe du DVH -- l'erreur était invisible sur la dispersion.
        """
        angle = self.angle_ancre_deg
        if not self.collet_cote_ancre:
            angle += 180.0
        return (angle + 180.0) % 360.0 - 180.0

    @property
    def angle_ancre_deg(self) -> float:
        """Direction, depuis l'axe du balancier, du DEMI-AXE de la ligne des
        centres qui porte l'ancre et la roue d'échappement.

        Le tracé est en vue côté rouage, 6H en haut : le repère 3H, d'où sort
        la tige, est alors à DROITE, c'est-à-dire à 0°. L'angle de la ligne des
        centres étant compté depuis l'axe de tige et son demi-axe utile
        pointant du côté du 3H, la direction cherchée est `angle_echappement_deg`
        lui-même (−6.2° par défaut : juste sous l'horizontale, vers la droite).

        Cette direction ne dépend PAS de `collet_cote_ancre` : celui-ci ne dit
        que sur lequel des deux demi-axes le spiral vient s'attacher, pas où se
        trouve l'échappement.
        """
        return (self.angle_echappement_deg + 180.0) % 360.0 - 180.0

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

    # --- Terminaison extérieure à RAQUETTE (alternative à la courbe de Phillips) ---
    # Quand raquette_active vaut True, la courbe de Phillips est remplacée par
    # une courbe en S suivie d'un ARC CONCENTRIQUE au rayon r_max. L'arc est ce
    # qu'exige une raquette : ses goupilles y glissent en tournant autour de
    # l'axe du balancier sans pousser le spiral radialement, de sorte que seule
    # la LONGUEUR ACTIVE change. Sur une courbe de Phillips (ou une spirale
    # d'Archimède), déplacer la raquette déformerait le spiral en plus d'en
    # changer la longueur.
    raquette_active: bool = False
    raquette_plage_s_par_j: float = 60.0        # ± cette valeur autour du nominal
    raquette_transition_angle_deg: float = 120.0
    raquette_n_points: int = 4                  # points de contrôle de la courbe en S
    raquette_offsets: np.ndarray = field(default_factory=lambda: np.zeros(0))
    raquette_transition_n_pas: int = 150
    # Étendue angulaire de l'arc concentrique AVANT le point de comptage,
    # c'est-à-dire entre la fin du raccord en S et les goupilles. Par défaut
    # égale à l'angle du piton au-delà, ce qui place le point de comptage au
    # MILIEU de la partie circulaire terminale. Cette portion fait partie du
    # spiral ACTIF : l'allonger remplace de la spirale d'Archimède par de l'arc
    # concentrique dans la partie qui travaille.
    raquette_arc_avant_deg: float = 30.0
    raquette_arc_avant_n_pas: int = 100
    # Étendue angulaire de l'arc AU-DELÀ du parcours des goupilles, jusqu'au
    # piton. Indépendante de la plage de réglage (diagnostic 2026-09-25) : le
    # piton n'est PAS à l'extrémité du parcours des goupilles. C'est cette
    # longueur morte qui permet à la lame de fléchir latéralement au droit des
    # goupilles, donc de venir en contact avec l'une d'elles -- sans quoi la
    # raquette n'aurait aucune autorité. Mesuré : l'excursion latérale varie
    # comme le CARRÉ de cette longueur (124 -> 371 µm de longueur morte fait
    # passer l'excursion de 0.9 à 7.6 µm sur ±340°).
    raquette_piton_angle_deg: float = 30.0
    raquette_piton_n_pas: int = 100
    # Jeu LIBRE total entre les goupilles, épaisseur de lame déjà déduite : la
    # lame peut donc s'écarter de ±jeu/2 avant de toucher. C'est ce que contrôle
    # un horloger à la loupe.
    raquette_jeu_um: float = 20.0

    # --- Compensation thermique par oxydation ---
    # L'épaisseur d'oxyde est RECALCULÉE pour chaque candidat : elle dépend du
    # profil d'épaisseur, que l'optimiseur fait justement varier (la part de
    # rigidité portée par l'oxyde passe de 48% sur une section de 20 µm à 22%
    # sur une section de 50 µm). `oxyde_impose_um` > 0 fige au contraire
    # l'épaisseur, pour étudier ce que donne une valeur que le procédé sait
    # tenir.
    oxydation_active: bool = False
    oxyde_impose_um: float = 0.0

    # --- Profils libres (valeurs aux points de contrôle, PCHIP sur s/L in [0,1]) ---
    epaisseur_noeuds: np.ndarray = field(default_factory=lambda: np.full(20, 0.030))
    delta_r_noeuds: np.ndarray = field(default_factory=lambda: np.zeros(12))

    @property
    def n_elements(self) -> int:
        n = self.n_elements_corps + self.grossmann_n_pas
        if self.raquette_active:
            return (n + self.raquette_transition_n_pas + self.n_pas_arc_avant
                    + self.n_pas_piton)
        return n + self.phillips_n_pas

    @property
    def n_pas_arc_avant(self) -> int:
        """Éléments de l'arc concentrique avant le point de comptage. Au moins
        un : le point de comptage est ce dernier élément."""
        if not self.raquette_active or self.raquette_arc_avant_deg <= 1e-9:
            return 1
        return max(1, self.raquette_arc_avant_n_pas)

    @property
    def n_pas_piton(self) -> int:
        """Éléments du prolongement d'arc jusqu'au piton. Nul si l'angle l'est :
        sinon on créerait des éléments de longueur nulle, qui font diverger la
        souplesse intégrée (division par ds nul) sans rien représenter."""
        if not self.raquette_active or self.raquette_piton_angle_deg <= 1e-9:
            return 0
        return self.raquette_piton_n_pas

    @property
    def angle_exterieur_deg(self) -> float:
        """Étendue angulaire de la zone terminale extérieure, hors arc de
        raquette (dont l'angle est DÉDUIT de la plage de réglage demandée,
        pas choisi -- voir angle_arc_raquette)."""
        return (self.raquette_transition_angle_deg if self.raquette_active
                else self.phillips_angle_deg)


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
    # n_tours_reels compte du collet au POINT DE COMPTAGE (piton sur un spiral
    # libre, goupilles sur un spiral à raquette), et non jusqu'au bout de la
    # matière. Seule la PREMIÈRE MOITIÉ du parcours des goupilles entre donc
    # dans le budget angulaire ; la seconde moitié et l'arc qui mène au piton
    # viennent en plus. C'est la convention horlogère : un nombre entier de
    # tours place le point de comptage sur le même rayon que l'attache à la
    # virole, ce qui aligne sur ce rayon l'axe de la force que le spiral
    # exerce sur le pivot du balancier (vérifié : l'axe suit la position
    # angulaire du point de comptage à 1-4° près).
    theta_ext = np.radians(p.angle_exterieur_deg)
    if p.raquette_active:
        theta_ext = theta_ext + np.radians(p.raquette_arc_avant_deg)
    theta_corps = max(theta_total - theta_gross - theta_ext, 1e-6)
    r_body_start = p.r_min + pas_theorique * theta_gross
    r_body_end = r_body_start + pas_theorique * theta_corps
    return pas_theorique, r_body_start, r_body_end


def angle_arc_raquette(geom: GeometrieSpiral, p: ParametresSpiral) -> float:
    """Étendue angulaire (rad) de l'arc concentrique donnant la plage de
    réglage demandée, déduite de la géométrie fournie.

    Déplacer la raquette le long de l'arc raccourcit ou rallonge la partie
    ACTIVE du spiral. Pour une poutre en flexion, la raideur varie comme
    l'inverse de la souplesse intégrée S = somme(ds / e³), donc :

        marche (s/j) = 86400 · ½ · ΔS / S       avec  ΔS = L_arc / e_piton³

    La plage totale vaut deux fois raquette_plage_s_par_j (réglage de part et
    d'autre du nominal, point de comptage au MILIEU de l'arc), d'où :

        L_arc = 4 · plage · S · e_piton³ / 86400        θ_arc = L_arc / r_max

    La sensibilité varie comme 1/e³ à l'extrémité : un piton aminci rend la
    raquette beaucoup plus sensible et raccourcit d'autant l'arc nécessaire.
    C'est pourquoi l'angle est déduit de la géométrie réelle et non d'une
    estimation à épaisseur uniforme."""
    # Souplesse de la PARTIE ACTIVE : c'est elle qui fixe la fréquence, donc
    # elle qui détermine de combien un déplacement des goupilles la change.
    # Utiliser la souplesse totale surestimait la plage de la fraction que
    # représente la longueur morte -- imperceptible tant que le piton était
    # collé aux goupilles, mais 1.4% d'erreur dès qu'on l'en éloigne de 30°
    # (diagnostic 2026-09-25). Exprimée en rigidité E·I et non en e³, pour
    # rester valable sur une lame oxydée, dont la rigidité n'est plus
    # proportionnelle au cube de l'épaisseur.
    actif = geom.partie_active()
    souplesse = actif.souplesse()
    rigidite_piton = float(actif.rigidite_flexion()[-1])
    L_arc_m = 4.0 * float(p.raquette_plage_s_par_j) * souplesse * rigidite_piton / 86400.0
    return (L_arc_m * 1e3) / p.r_max


def construire_geometrie(p: ParametresSpiral, return_diagnostics: bool = False,
                          reglage_thermique=None):
    """Construit la géométrie complète et applique épaisseur et écart radial
    libres. Retourne une GeometrieSpiral (et un dict de diagnostics si
    demandé -- validité géométrique, monotonie des courbes, repère de
    plausibilité Phillips -- à utiliser comme pénalité/contrainte côté
    optimiseur plutôt que de planter).

    Deux terminaisons extérieures possibles :
      - spiral LIBRE  : collet -> Grossmann -> corps -> Phillips -> piton ;
      - spiral à RAQUETTE : collet -> Grossmann -> corps -> courbe en S ->
        arc concentrique -> piton, le point de comptage (goupilles) étant au
        milieu de l'arc.

    En mode raquette la construction se fait en DEUX passes : l'angle de l'arc
    dépend de la souplesse intégrée de la lame, qui dépend elle-même de la
    géométrie. Une première passe sans arc donne cette souplesse, la seconde
    construit la géométrie définitive. L'arc ne représentant que quelques
    dixièmes de pour-cent de l'étendue angulaire, une seule itération suffit
    largement ; le coût est négligeable devant le balayage elastica qui suit."""
    resultat = _construire(p, return_diagnostics, reglage_thermique)
    if not p.raquette_active:
        return resultat

    # Le parcours des goupilles se DÉDUIT de la géométrie (il dépend de la
    # souplesse et de l'épaisseur au point de comptage) mais ne la modifie
    # plus : l'arc concentrique a une étendue fixée par ses propres
    # paramètres. Une seule passe suffit donc désormais, là où l'ancienne
    # découpe -- où le parcours occupait l'arc -- imposait un point fixe.
    geom = resultat[0] if return_diagnostics else resultat
    demi_parcours = 0.5 * angle_arc_raquette(geom, p) * p.r_max
    geom = replace(geom, demi_parcours_raquette=demi_parcours)
    if not return_diagnostics:
        return geom

    diagnostics = resultat[1]
    diagnostics["raquette_arc_deg"] = float(np.degrees(2.0 * demi_parcours / p.r_max))
    diagnostics["raquette_arc_mm"] = float(2.0 * demi_parcours)
    # Les goupilles doivent rester sur l'arc concentrique de part et d'autre du
    # point de comptage : au-delà, elles quitteraient la portion circulaire et
    # déformeraient le spiral radialement au lieu de n'en changer que la
    # longueur active.
    marge = min(p.raquette_arc_avant_deg, p.raquette_piton_angle_deg)
    diagnostics["raquette_parcours_tient_dans_arc"] = bool(
        np.degrees(demi_parcours / p.r_max) <= marge + 1e-9)
    diagnostics["raquette_demi_parcours_deg"] = float(np.degrees(demi_parcours / p.r_max))
    return geom, diagnostics


def _construire(p: ParametresSpiral, return_diagnostics: bool, reglage_thermique=None):
    if len(p.grossmann_offsets) != p.grossmann_n_points - 3:
        raise ValueError(f"grossmann_n_points={p.grossmann_n_points} implique "
                          f"{p.grossmann_n_points - 3} offsets, reçu {len(p.grossmann_offsets)}")
    if p.raquette_active:
        if len(p.raquette_offsets) != p.raquette_n_points - 4:
            raise ValueError(f"raquette_n_points={p.raquette_n_points} implique "
                              f"{p.raquette_n_points - 4} offsets, reçu {len(p.raquette_offsets)}")
    elif len(p.phillips_offsets) != p.phillips_n_points - 3:
        raise ValueError(f"phillips_n_points={p.phillips_n_points} implique "
                          f"{p.phillips_n_points - 3} offsets, reçu {len(p.phillips_offsets)}")

    theta_total = p.n_tours_reels * 2.0 * np.pi
    theta_gross = np.radians(p.grossmann_angle_deg)
    # n_tours_reels compte du collet au POINT DE COMPTAGE (piton sur un spiral
    # libre, goupilles sur un spiral à raquette), et non jusqu'au bout de la
    # matière. Seule la PREMIÈRE MOITIÉ du parcours des goupilles entre donc
    # dans le budget angulaire ; la seconde moitié et l'arc qui mène au piton
    # viennent en plus. C'est la convention horlogère : un nombre entier de
    # tours place le point de comptage sur le même rayon que l'attache à la
    # virole, ce qui aligne sur ce rayon l'axe de la force que le spiral
    # exerce sur le pivot du balancier (vérifié : l'axe suit la position
    # angulaire du point de comptage à 1-4° près).
    theta_ext = np.radians(p.angle_exterieur_deg)
    if p.raquette_active:
        theta_ext = theta_ext + np.radians(p.raquette_arc_avant_deg)
    theta_corps = theta_total - theta_gross - theta_ext
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

    # --- Zone terminale exterieure ---
    theta_jonction = theta_gross + theta_corps
    theta_trans = np.radians(p.angle_exterieur_deg)
    if p.raquette_active:
        # Courbe en S : ancree en pente AUX DEUX BOUTS (pas d'Archimede cote
        # corps, pente NULLE cote arc), suivie de l'arc concentrique a r_max.
        zone_ext = construire_zone_s(
            p.raquette_n_points, r_debut=r_body_end,
            pente_debut_du=pas_theorique * theta_trans,
            r_fin=p.r_max, pente_fin_du=0.0, offsets_libres=p.raquette_offsets,
        )
        u_ext = np.linspace(0.0, 1.0, p.raquette_transition_n_pas, endpoint=False)
        theta_ext_arr = theta_jonction + u_ext * theta_trans
        r_ext = zone_ext.r(u_ext)

        # Arc concentrique, en DEUX portions d'étendue fixée : avant le point
        # de comptage (partie active) et au-delà jusqu'au piton (partie morte).
        # Le parcours des goupilles n'entre plus dans cette découpe -- il est
        # déduit de la plage de réglage et se contente de devoir tenir DANS
        # l'arc, ce qui supprime le couplage circulaire entre la plage demandée
        # et le budget angulaire du spiral.
        theta_debut_arc = theta_jonction + theta_trans
        theta_avant = np.radians(p.raquette_arc_avant_deg)
        n_avant = p.n_pas_arc_avant
        theta_arc_arr = theta_debut_arc + np.linspace(
            theta_avant / n_avant, theta_avant, n_avant)
        r_arc = np.full(n_avant, p.r_max)

        n_pit = p.n_pas_piton
        if n_pit > 0:
            theta_piton = np.radians(p.raquette_piton_angle_deg)
            theta_pit_arr = theta_debut_arc + theta_avant + np.linspace(
                theta_piton / n_pit, theta_piton, n_pit)
            theta_arc_arr = np.concatenate([theta_arc_arr, theta_pit_arr])
            r_arc = np.concatenate([r_arc, np.full(n_pit, p.r_max)])

        theta_final = np.concatenate([theta_g_arr, theta_b_arr, theta_ext_arr, theta_arc_arr])
        r_ref = np.concatenate([r_g, r_b, r_ext, r_arc])

        # Point de comptage = dernier élément de l'arc AVANT, donc à
        # raquette_arc_avant_deg de la fin du raccord en S.
        i_arc_debut = len(theta_g_arr) + len(theta_b_arr) + len(theta_ext_arr)
        i_comptage = i_arc_debut + n_avant - 1
        n_actif = i_comptage + 1
    else:
        zone_ext = construire_zone(
            p.phillips_n_points, r_ancre=r_body_end, pente_ancre_du=pas_theorique * theta_trans,
            r_libre=p.r_max, offsets_libres=p.phillips_offsets, ancre_a_u1=False,
        )
        u_ext = np.linspace(0.0, 1.0, p.phillips_n_pas)
        theta_ext_arr = theta_jonction + u_ext * theta_trans
        r_ext = zone_ext.r(u_ext)

        theta_final = np.concatenate([theta_g_arr, theta_b_arr, theta_ext_arr])
        r_ref = np.concatenate([r_g, r_b, r_ext])
        i_arc_debut = None
        i_comptage = len(theta_final) - 1
        n_actif = len(theta_final)
    zone_ph = zone_ext

    # --- Ecart radial libre, PCHIP sur l'indice normalise (proxy de s/L avant deformation) ---
    u_idx = np.linspace(0.0, 1.0, len(theta_final))
    if len(p.delta_r_noeuds) >= 2:
        u_noeuds_r = np.linspace(0.0, 1.0, len(p.delta_r_noeuds))
        delta_r = PchipInterpolator(u_noeuds_r, p.delta_r_noeuds)(u_idx)
    else:
        delta_r = np.zeros_like(theta_final)
    r_final = r_ref + delta_r

    # Orientation dans le mouvement : le collet est posé sur la ligne
    # d'échappement, et le spiral se développe dans le sens choisi. Le nombre
    # de tours étant entier, le point de comptage retombe sur cette même ligne.
    sens = 1.0 if p.sens_trigonometrique else -1.0
    theta_final = np.radians(p.angle_collet_deg) + sens * theta_final

    x = r_final * np.cos(theta_final)
    y = r_final * np.sin(theta_final)

    # --- Epaisseur, PCHIP sur l'abscisse curviligne normalisee reelle ---
    ds = np.hypot(np.diff(x, prepend=x[0]), np.diff(y, prepend=y[0]))
    ds[0] = 0.0
    s_cum = np.cumsum(ds)
    s_norm = s_cum / s_cum[-1]
    u_noeuds_e = np.linspace(0.0, 1.0, len(p.epaisseur_noeuds))
    epaisseur = PchipInterpolator(u_noeuds_e, p.epaisseur_noeuds)(s_norm)

    # --- Mise a l'echelle de l'epaisseur pour tomber pile sur freq_cible ---
    # Sur la partie ACTIVE seulement : au-dela du point de comptage, la lame est
    # tenue entre les goupilles et le piton, elle ne participe pas au rappel.
    from .physique import frequence_propre

    def _geom(profil, t_mm):
        return GeometrieSpiral(x=x, y=y, epaisseur=profil, hauteur=p.hauteur_h,
                                n_actif=n_actif, epaisseur_oxyde=t_mm,
                                plan_cristallin=p.plan_cristallin,
                                orientation_cristal_deg=p.orientation_cristal_deg)

    # La fréquence de départ est MESURÉE sur la géométrie, pas recalculée par
    # une formule à module constant : sur un wafer anisotrope, E varie le long
    # de la lame et aucune expression fermée ne donne la souplesse intégrée.
    # La loi d'échelle beta^(3/2), elle, reste exacte — E ne dépend pas de
    # l'épaisseur.
    freq_actuelle = frequence_propre(_geom(epaisseur, 0.0), p.inertie_balancier)
    beta = (p.freq_cible / freq_actuelle) ** (2.0 / 3.0)
    t_oxyde_mm = 0.0

    if p.oxydation_active:
        # Épaisseur d'oxyde et mise à l'échelle sont COUPLÉES : l'oxyde change
        # la rigidité (donc la fréquence, donc beta), et beta change
        # l'épaisseur de lame (donc l'oxyde requis, qui dépend du rapport
        # oxyde/silicium). Point fixe sur les deux -- il converge en deux ou
        # trois tours, l'oxyde ne pesant que quelques pour-cent de la rigidité.
        # La formule analytique beta = (f_cible/f)^(2/3) ne vaut plus ici : la
        # rigidité d'un composite n'est pas proportionnelle au cube de
        # l'épaisseur. On résout donc numériquement.
        from .thermique import ReglageThermique, epaisseur_oxyde_optimale
        reglage_th = reglage_thermique or ReglageThermique()
        for _ in range(4):
            # ORDRE IMPORTANT : l'épaisseur d'oxyde d'abord, la mise à
            # l'échelle ensuite. beta est ainsi toujours résolu avec l'oxyde
            # définitif, donc la fréquence tombe exactement sur la cible --
            # l'ordre inverse laissait un résidu (et jusqu'à 1% d'écart quand
            # l'épaisseur d'oxyde est imposée).
            t_precedent = t_oxyde_mm
            if p.oxyde_impose_um > 0.0:
                t_oxyde_mm = p.oxyde_impose_um * 1e-3
            else:
                t_um = epaisseur_oxyde_optimale(_geom(epaisseur * beta, t_oxyde_mm),
                                                 p.inertie_balancier, reglage_th)
                t_oxyde_mm = t_um * 1e-3 if np.isfinite(t_um) else 0.0

            def ecart_frequence(facteur):
                return (frequence_propre(_geom(epaisseur * facteur, t_oxyde_mm),
                                         p.inertie_balancier) - p.freq_cible)
            try:
                beta = float(brentq(ecart_frequence, 0.2 * beta, 5.0 * beta,
                                    xtol=1e-12, rtol=8.9e-16))
            except ValueError:
                break
            # Convergé dès que l'épaisseur d'oxyde ne bouge plus d'un nanomètre.
            if p.oxyde_impose_um > 0.0 or abs(t_oxyde_mm - t_precedent) < 1e-6:
                break

    epaisseur = epaisseur * beta
    geom = _geom(epaisseur, t_oxyde_mm)

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
        "raquette_active": bool(p.raquette_active),
        "raquette_arc_avant_deg": float(p.raquette_arc_avant_deg),
        "raquette_piton_angle_deg": float(p.raquette_piton_angle_deg),
        "i_comptage": int(i_comptage),
        "n_actif": int(n_actif),
        "longueur_morte_mm": float(s_cum[-1] - s_cum[n_actif - 1]),
        "oxydation_active": bool(p.oxydation_active),
        "epaisseur_oxyde_um": float(t_oxyde_mm * 1000.0),
        "croissance_lame_um": float(2.0 * 0.56 * t_oxyde_mm * 1000.0),
    }
    return geom, diagnostics
