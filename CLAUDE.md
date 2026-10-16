# WatchTools

Famille d'outils de calcul horloger : `came/` (coeur de chronographe et marteau),
`engrenage/` (profils de denture, engrènement, optimisation, fraise), `spiral/`
(spiral : courbes terminales, isochronisme, optimisation), `echappement/`
(échappement à ancre suisse : tracé, rendement statique et dynamique, retard,
optimisation). Chaque outil est un paquet Python et un programme à part ; ils
partagent `commun/` (thème sombre, widgets de panneau, navigation souris des
vues, icônes, lanceur). Le `main.py` de la racine est le lanceur : une icône
par outil, chaque outil démarré dans son propre processus
(`python -m <outil>.main`). Ajouter un outil : une ligne dans `OUTILS` de
`main.py` et une fonction de dessin dans `commun/icones.py`.

Chaque dossier d'outil a son propre `CLAUDE.md` ; celui-ci ne contient que ce
qui est commun à tous. L'utilisateur est francophone : interface, commentaires
et documentation en français (le `CLAUDE.md` d'engrenage est resté en anglais).

## Environnement

Interpréteur : `/home/florian/DEV/virtual_env/.watchtools/bin/python`
(PySide6, numpy, scipy, numba, pymoo, joblib, matplotlib, python-docx, ezdxf,
cma). Le `python3` du système n'a pas ces dépendances.

## Lancement

Depuis la racine :

    python main.py                  le lanceur
    python main.py <outil>          un outil sans passer par le lanceur
    python -m <outil>.main          un outil seul
    python -m <outil>.<module>      un module en ligne de commande

Depuis le dossier d'un outil : `python main.py` (chaque `main.py` met
lui-même la racine dans `sys.path` quand il est lancé comme script, et
importe ses modules en absolu). Les autres modules d'un outil utilisent des
imports relatifs et se lancent en `python -m <outil>.<module>` avec la racine
dans le chemin : c'est le cas dans le terminal VSCode d'un dossier d'outil.

Chaque dossier d'outil a un `.vscode/` pour être ouvert seul : `settings.json`
(interpréteur `.watchtools`, `python.analysis.extraPaths: [".."]` pour que
Pylance trouve `commun`, `PYTHONPATH=..` dans le terminal intégré) et
`launch.json` (l'outil, et « module courant » = `python -m <outil>.<fichier
ouvert>`, plus quelques entrées propres à l'outil). Les chemins par défaut
des outils sont tous ancrés sur `__file__`, jamais sur le répertoire courant.

## Sessions en parallèle

L'utilisateur ouvre un VSCode, et donc une session Claude, par dossier
d'outil, pour faire avancer plusieurs outils en même temps. Il n'y a qu'un
dépôt git et un seul arbre de travail pour tous :

- chaque session voit les modifications non commitées des autres ; ne jamais
  « nettoyer », reformater ou annuler ce qui est hors de son dossier ;
- commiter seulement son dossier, quand l'utilisateur le demande :
  `git add <outil>/` puis `git commit`, ou `git commit -- <outil>/` ; jamais
  `git add -A`, `git commit -a`, `git stash` ni changement de branche ;
- `commun/` et le `main.py` de la racine servent à tous les outils : n'y
  toucher que si c'est nécessaire, le dire à l'utilisateur, et vérifier que
  les autres outils démarrent encore ;
- deux optimisations lancées en même temps se partagent le processeur : chaque
  outil a un nombre de processus réglable.

## Interface commune

Panneau de réglages à gauche (onglets), vues matplotlib à droite. Une figure
remplit sa vue : `Figure(layout="constrained")` et, pour un dessin à l'échelle,
`set_aspect("equal", adjustable="datalim")` (étend les données au lieu de
rétrécir le cadre). Navigation souris de type CAO sur chaque canevas :
`commun.navigation.activer(canvas)` (molette = zoom sous le curseur, glisser =
déplacer, double-clic = vue d'ensemble). Couleurs par `commun.theme`
(`theme.encre()`, `theme.gris(n)`, `ACC`…`ACC4`), jamais `"k"` ni un gris fixe.
L'utilisateur veut des infobulles sur tous les contrôles :
`commun.widgets.aide(widget, texte)` (un texte long va en infobulle, pas dans
le panneau).

## Tester l'interface sans écran

`QT_QPA_PLATFORM=offscreen`, la racine dans `PYTHONPATH`. Neutraliser
`QMessageBox.critical/information/warning/question` et `QFileDialog` (une boîte
non neutralisée bloque sans rien dire). Les calculs tournent dans des
`QThread` : attendre réellement (`app.processEvents()` + `time.sleep` dans une
boucle) avant de lire le résultat, puis `fenetre.grab().save("x.png")` et
regarder l'image. Écrire ces scripts hors du dépôt.

## Git

Dépôt `https://github.com/bernardflorex/90-WatchTools.git`. Les historiques
des dépôts d'origine (015_eng_cre, 040-Spiral branche histoire-ancienne,
050-came) ont été raccordés par une fusion (code seulement) ; le commit de
copie 3734500 leur est greffé localement par `git replace` : sans les refs
`refs/replace/*`, `log` et `blame` des fichiers repris s'arrêtent à cette
copie. Le dossier `home/` de la racine est un environnement créé par erreur
(chemin relatif), non suivi par git (son contenu tombe sous le `.gitignore`) ;
ne pas le supprimer sans l'accord de l'utilisateur.
