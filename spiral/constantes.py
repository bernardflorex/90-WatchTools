# constantes.py
"""Constantes matériau. Le silicium est pris en orientation <111>, qui est
transversalement isotrope dans le plan : le module d'Young y est le même
quelle que soit la direction de flexion dans le plan du spiral, ce qui
autorise un scalaire au lieu d'un tenseur."""

E_SILICIUM = 169e9      # Module de Young du silicium <111> (Pa)
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
