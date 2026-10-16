# Échappement à ancre suisse

Simulation 2D d'un échappement à ancre suisse (roue de 20 dents, pointage
20.3) : rendement de la roue au balancier, en statique et en dynamique, pour
la palette d'entrée et la palette de sortie ; retard (marche) induit par
l'échappement ; optimisation NSGA-II du tracé. Règles communes à tous les
outils (environnement, sessions parallèles, git, interface) : `CLAUDE.md` de
la racine.

## Lancement

    python main.py [tracé.dxf]               depuis ce dossier (ou « Run » / « Échappement » dans VSCode)
    python -m echappement.main [tracé.dxf]   depuis la racine
    python -m echappement.dynamique          une oscillation en console (terminal VSCode de ce dossier)

Tracé par défaut : `simulation.dxf` (dessin simplifié : roue, ancre et
fourchette, deux butées, ellipse du balancier).

## Modules

- `geometrie.py` : `Element` (segment ou arc, étiqueté), `Corps` (pièce qui
  tourne autour d'un axe), distance entre contours en facettes, et surtout
  `premier_contact_exact` : plus petite rotation d'un corps qui l'amène au
  contact d'un autre, sur la géométrie EXACTE (segments, vrais arcs et leurs
  extrémités ; sommet sur segment ou arc, tangence arc-segment et arc-arc,
  toutes en fermé). Les facettes (`PAS_ARC_MAX` = 3°) ne servent plus qu'au
  dessin et au test de liberté.
- `dxf_import.py` : lecture, regroupement par connexité, étiquetage
  automatique (face de repos, plan d'impulsion, arêtes, dos de palette,
  flancs, cornes, baguette, ellipse, butées), orientation en vue standard.
- `tables.py` : tables de contact en fonction de l'angle d'ancre (pas 0,02°,
  grille débordant de 0,5° au-delà des butées) : position de la roue au
  contact de chaque palette, bornes du balancier imposées par la fourchette.
  Toute la mécanique lit ces tables (environ 8 s par tracé).
- `cinematique.py` : relevé (`Releve`) : fin du dégagement, impulsion, chute,
  chemin perdu, repos, recul, levée, angles par rapport à la ligne des centres.
- `statique.py` : efforts quasi statiques avec frottement de Coulomb,
  rendement, retard par la formule d'Airy.
- `parametres.py` : paramètres physiques (`Parametres`).
- `dynamique.py` : dynamique à événements, 3 degrés de liberté (roue, ancre,
  balancier), contacts unilatéraux exprimés en jeux angulaires lus dans les
  tables (splines cubiques d'interpolation), multiplicateurs de Lagrange +
  Baumgarte, frottement régularisé, chocs de Newton avec contacts persistants,
  règle de Zénon, bilan d'énergie fermé (contrôle), marche par la phase de
  l'oscillateur. Exploitation : `couple_fonction`, `decollement`,
  `fenetres_action`, `amplitude_equilibre`. Environ 1,2 s par oscillation.
- `trace_parametrique.py` : 13 variables qui modifient le tracé de référence
  (entraxes, ellipse, fourchette, course, et par palette position,
  inclinaison et longueur du plan d'impulsion, tirage).
- `optimisation.py` : NSGA-II (pymoo), évaluation statique (~4 s) ou
  dynamique (~5 s, ~10 s si l'amplitude d'équilibre est un objectif).
- `app.py` : fenêtre ; à gauche Tracé / Calculs / Optimisation / Résultats ;
  à droite Échappement (animation 100 images, statique ou dynamique), Couple
  au balancier, Cinématique, Énergie, Retard, Pareto.

## Décisions de l'utilisateur

- Vue standard : balancier en haut, roue en bas, palette d'entrée à gauche ;
  la roue tourne dans le sens horaire.
- Chemin perdu : course (µm) de la dent sur la face de repos entre le premier
  contact de chute et le repos, ancre sur la butée.
- Données : couple sur la roue 1,13 µN·m ; µ = 0,15 ; 3,5 Hz ; amplitude
  visée 310° ; facteur de qualité 220 (horizontal, 310° à 250°, frottement de
  l'air compris : pas de modèle d'air séparé). Inerties actuelles pour des
  pièces acier/rubis, pas encore silicium : roue 0,021 + pignon 0,0002 g·mm²
  (on ne les change pas pour « arranger » un résultat), rouage ramené
  0,0007 g·mm², ancre avec dard 0,00814 g·mm², balancier 8 mg·cm² (comme le
  programme spiral).
- Restitution 0 par défaut : aucun rebond n'est visible à la caméra rapide.
- Toutes les pièces sont rigides vis-à-vis des couples en jeu (pas de modèle
  élastique du rouage).
- Départ de la simulation : roue et ancre au repos sur la palette d'entrée ;
  l'ellipse arrive sur la fourchette avec la vitesse que lui donne son
  amplitude, et c'est elle qui lance le dégagement.
- Couple au balancier : positif quand l'ancre est motrice, négatif quand le
  balancier fournit le dégagement, nul hors contact.
- Le cycle (animation, tableaux, colonnes E/S) commence par la palette
  d'entrée.

## Résultats de référence (`simulation.dxf`, paramètres par défaut)

Statique : rendement 47,0 %, retard −6,65 s/j à 310°. Dynamique : rendement
37,1 % (l'utilisateur mesure environ 36 % sur l'échappement réel), marche
−8,88 s/j, amplitude 310° → 309,17° en une oscillation ; amplitude
d'équilibre 284,7° (rendement 38,6 %, marche −10,9 s/j). Décollement de la
dent après le dégagement : 6,2 µm à l'entrée (reprise à 13 % du plan
d'impulsion), 1,75 µm à la sortie. Levées 48,67° (entrée) et 49,38° (sortie) :
l'écart vient du dessin (la palette de sortie lâche sa dent 0,17° d'ancre plus
tard). Toute modification du moteur de calcul doit retrouver ces chiffres.

## Pièges connus

- `Releve.alternances` liste la palette de SORTIE en premier. Toujours
  indexer par `alt.palette` ou trier (`key=lambda a: a.palette != "entrée"`) :
  l'ordre brut a déjà inversé l'animation statique et les colonnes E/S des
  résultats d'optimisation.
- Facettes : un contour échantillonné donne des positions exactes mais des
  rapports de transmission en escalier ; en dynamique cela faisait des
  pointes de couple de 750 nN·m et des pas d'intégration de 1 µs (100 s par
  oscillation). D'où `premier_contact_exact`. Ne pas revenir aux facettes
  pour les contacts.
- Lisser les tables (spline de lissage, validation croisée) a donné des
  dérivées de mauvais signe : on garde des splines d'interpolation.
- Fin d'oscillation : l'événement de retour n'est armé que trois quarts de
  période après le départ ou le retour précédent, sinon le segment qui repart
  sur le point de départ le compte deux fois.
- Tolérances d'intégration 1e-8 / 1e-11 : 1e-10 ne change rien à 4 décimales.
- Il n'y a pas de suite de tests : vérifier par les chiffres ci-dessus et par
  l'interface sans écran (recette dans le `CLAUDE.md` de la racine).

## Reste à faire (quand l'utilisateur le demandera)

Inerties silicium ; calibrage du tracé sur le rendement mesuré ; couplage de
la courbe de marche de l'échappement avec le module spiral.
