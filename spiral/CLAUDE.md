# Spiral

Laboratoire de spiral : géométrie (import DXF ou génération paramétrique,
courbes terminales), physique (élastica), isochronisme et marche diurne,
compensation thermique, optimisation NSGA-II. Règles communes à tous les
outils (environnement, sessions parallèles, git, interface) : `CLAUDE.md` de
la racine.

## Lancement et tests

    python main.py                          depuis ce dossier (ou « Run » / « Spiral » dans VSCode)
    python -m spiral.main                   depuis la racine
    python tests/test_spiral_lab.py         ~200 contrôles, quelques minutes ; marche de n'importe où
                                            (le script met la racine dans sys.path)

## Modules (d'après leurs en-têtes, à lire avant de modifier)

`modele.py` (représentation canonique unique d'un spiral), `parametrage.py`
(génération paramétrique), `courbes.py` (courbes de raccordement),
`physique.py` (solveur élastica unique, méthode de tir), `isochronisme.py`
(marche diurne, trois contributions), `echappement.py` (erreur de marche due à
l'échappement), `raquette.py` (jeu entre goupilles), `thermique.py`
(compensation par couche d'oxyde), `collision.py` (jeu inter-spires),
`constantes.py` (matériaux ; silicium <111>), `optimisation.py` (NSGA-II
pymoo), `presets.py` + `presets/` (réglages et bornes en JSON),
`sauvegarde.py` / `rapport.py` (fichiers produits, dossier de projet dans
`Results/`), `dxf_io.py`, `ancien_format.py`, `app.py` (fenêtre),
`threads.py`, `rendu.py`, `convergence.py`.

## Points récents

- `ancien_format.py` relit les « expériences » JSON des anciens programmes
  (v4 à v7, par exemple `Results/spiralv4_1.json`). Les champs absents de
  l'ancien fichier sont marqués « à valider » dans l'interface
  (`app._marquer_a_valider`) jusqu'à ce que l'utilisateur les confirme.
- `physique.frequence_propre` corrige deux biais de la formule de flexion
  pure, mesurés le 14.10.2026 sur les spirals de référence : la raideur
  inclut la réaction du pivot (`raideur_lineaire`, +80 à +125 s/j) et
  l'inertie du spiral s'ajoute à celle du balancier
  (`inertie_equivalente_spiral`, méthode de Rayleigh, −40 à −80 s/j).
- La fréquence propre ne concorde toujours pas avec le programme d'un
  collègue de l'utilisateur ; on attend ses constantes du silicium et de
  l'oxyde de silicium. Son traitement en non-linéarité géométrique
  (Green-Lagrange) a été jugé inutile ici : l'élastica est déjà exacte en
  grandes rotations, et les déformations du spiral restent petites.
- Le module `echappement/` (outil séparé) calcule la marche induite par un
  échappement réel ; son couplage avec ce module est prévu, pas fait.
