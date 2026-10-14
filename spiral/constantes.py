# constantes.py
"""Constantes matériau. Le silicium est pris en orientation <111>, qui est
transversalement isotrope dans le plan : le module d'Young y est le même
quelle que soit la direction de flexion dans le plan du spiral, ce qui
autorise un scalaire au lieu d'un tenseur."""

# E_SILICIUM est DÉFINI PLUS BAS, à partir des souplesses du cristal : c'est le
# module dans le plan d'un wafer (111), seul plan isotrope dans le plan. Il
# valait 169e9 en dur ; la valeur exacte est 168,92 GPa, et les deux écritures
# ont fini par se contredire de 0,05 % (détecté par le contrôle de fréquence
# absolue, 2026-10-08). Une seule source de vérité : les souplesses.
RHO_SILICIUM = 2330.0   # Masse volumique du silicium (kg/m^3)
G_TERRE = 9.81          # m/s^2

# --- Oxyde thermique de silicium (compensation thermique) ---
# Le silicium a un coefficient thermoélastique fortement NÉGATIF : il
# s'assouplit en chauffant, et un spiral nu retarde d'environ 3.8 s/j/°C.
# L'oxyde thermique a le coefficient inverse, largement POSITIF. Une couche de
# quelques microns, placée là où elle travaille le plus (la périphérie de la
# section, la plus éloignée de la fibre neutre), suffit donc à annuler l'effet.
E_OXYDE = 70e9          # Module de Young de l'oxyde thermique (Pa)
RHO_OXYDE = 2200.0      # Masse volumique de l'oxyde (kg/m^3)

# Coefficients thermoélastiques (1/E)·dE/dT, en 1/K. Valeurs de littérature :
# à confirmer par mesure sur le procédé réel, ce sont elles qui fixent
# l'épaisseur d'oxyde à viser.
TCE_SILICIUM = -60e-6
TCE_OXYDE = +215e-6

# Termes du SECOND ordre (1/E)·d²E/dT², en 1/K². Nuls par défaut faute de
# données fiables -- et c'est une limite assumée du modèle : avec des modules
# strictement linéaires en température, le résidu après compensation ressort
# quasi nul (mesuré 0.008 s/j crête-à-crête sur 8-38 °C), ce qui est
# certainement optimiste. C'est la courbure réelle de E(T) qui limite les
# chronomètres. À renseigner dès que mesurée.
TCE2_SILICIUM = 0.0
TCE2_OXYDE = 0.0

# Dilatation thermique (1/K)
CTE_SILICIUM = 2.6e-6
CTE_OXYDE = 0.5e-6

# Géométrie de l'oxydation thermique : croître 1 µm d'oxyde consomme 0.44 µm
# de silicium et dépasse de 0.56 µm vers l'extérieur (rapport volumique 2.27).
CONSOMMATION_OXYDATION = 0.44
CROISSANCE_OXYDATION = 0.56


# --- Anisotropie élastique du silicium -------------------------------------
# E_SILICIUM ci-dessus (169 GPa) n'est PAS une propriété du silicium : c'est sa
# valeur dans le plan d'un wafer (111), le seul plan où le module soit isotrope
# dans le plan. Sur un wafer (001), le module varie de 130 à 169 GPa selon la
# direction, avec une période de 90° — 30 % d'écart, soit 13 % sur la raideur
# intégrée d'un spiral qui balaie toutes les directions. L'ignorer est une
# erreur de modèle, pas une approximation (remonté par l'utilisateur,
# 2026-10-08).
#
# Souplesses du silicium monocristallin, en Pa^-1. Ce sont elles les vraies
# constantes matériau ; tous les modules s'en déduisent.
S11_SILICIUM = 7.68e-12
S12_SILICIUM = -2.14e-12
S44_SILICIUM = 12.6e-12

# Base orthonormée DANS LE PLAN de chaque wafer, exprimée dans le repère
# cristallin. Le premier vecteur est l'axe de référence : c'est par rapport à
# lui qu'est compté l'angle d'orientation du dessin sur le wafer.
_BASES_PLAN = {
    #        axe de référence           axe à 90° de celui-ci
    "001": ((1.0, 0.0, 0.0),            (0.0, 1.0, 0.0)),            # réf. <100>
    "111": ((1.0, -1.0, 0.0),           (1.0, 1.0, -2.0)),           # réf. <1-10>
    "110": ((0.0, 0.0, 1.0),            (1.0, -1.0, 0.0)),           # réf. <001>
}
PLANS_CRISTALLINS = tuple(_BASES_PLAN)


def module_young_silicium(plan: str, angle_rad):
    """Module de Young du silicium (Pa) dans le PLAN du wafer `plan`, pour une
    direction faisant `angle_rad` avec l'axe de référence de ce plan.

    Loi standard pour un cristal cubique, en fonction des cosinus directeurs
    (l1, l2, l3) de la direction de traction :

        1/E = s11 - 2·(s11 - s12 - s44/2)·(l1²l2² + l2²l3² + l3²l1²)

    C'est bien le module en TRACTION-COMPRESSION uniaxiale le long de la lame
    qui intervient : en flexion plane, les fibres s'allongent et se raccourcissent
    dans la direction de la poutre, pas transversalement.

    Valeurs de contrôle : E<100> = 130,2 GPa, E<110> = 168,9 GPa,
    E<111> = 187,5 GPa ; dans le plan (111), E = 168,9 GPa quelle que soit la
    direction.
    """
    try:
        u, v = _BASES_PLAN[plan]
    except KeyError:
        raise ValueError(
            f"Plan cristallin inconnu : {plan!r}. Attendu l'un de "
            f"{', '.join(PLANS_CRISTALLINS)}.") from None
    import numpy as _np
    u = _np.asarray(u, dtype=float)
    v = _np.asarray(v, dtype=float)
    u /= _np.linalg.norm(u)
    v /= _np.linalg.norm(v)
    a = _np.asarray(angle_rad, dtype=float)
    # l a pour forme (..., 3) : une direction unitaire par angle.
    l = _np.cos(a)[..., None] * u + _np.sin(a)[..., None] * v
    l2 = l ** 2
    gamma = l2[..., 0] * l2[..., 1] + l2[..., 1] * l2[..., 2] + l2[..., 2] * l2[..., 0]
    s_aniso = S11_SILICIUM - S12_SILICIUM - S44_SILICIUM / 2.0
    return 1.0 / (S11_SILICIUM - 2.0 * s_aniso * gamma)


# Module dans le plan d'un wafer (111) : isotrope, c'est la valeur de référence
# du projet et le défaut partout où un module scalaire suffit. Déduite des
# souplesses ci-dessus, donc toujours cohérente avec module_young_silicium.
E_SILICIUM = float(module_young_silicium("111", 0.0))
