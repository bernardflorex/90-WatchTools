# Dimensionnement d'un coeur de chronographe — contexte de travail

Document de reprise. Il est ecrit pour qu'une session reparte d'ici et non de
zero : ce que le programme fait, comment le lancer et le verifier, les
conventions du modele, et surtout **les erreurs deja commises** avec leur
cause, pour ne pas les refaire.

Les echanges se font en francais. Tout chiffre cite ici vient d'une mesure
reelle, pas d'une estimation.

---

## 1. Ce que fait le programme

Il dimensionne un **coeur de chronographe** remis a zero par un marteau, d'un
bout a l'autre, dans un seul panneau. Deux phases :

1. **Le coeur**, pan de marteau suppose infini. Donnees : frottement,
   distance du plat a l'axe (nominal 0,6 mm), encombrement maxi, couple de
   rappel voulu, pivot et pente du pan. Sorties : la forme du coeur **et le
   ressort** (prearmage, raideur), calcules pour atteindre le couple voulu.
2. **Le pan reel** : longueur, decalage le long de lui-meme, rayons de coin.
   Le couple du marteau vient de la phase 1 ; on voit ce que le **vrai**
   marteau delivre.

La came simple et le galet ont ete **abandonnes** sur decision de
l'utilisateur. Ne pas les reintroduire.

Projet frere : `/home/florian/DEV/Python/040-Spiral` (Spiral Lab). L'objectif
annonce est une **famille d'outils d'optimisation horlogere en open source**,
avec une IHM commune. Un programme d'engrenages est prevu.

---

## 2. Environnement et lancement

Depuis la fusion dans WatchTools, ce dossier est le paquet `came` (l'ancien
sous-dossier `came/came/` a ete remonte d'un niveau) et l'environnement est
**`/home/florian/DEV/virtual_env/.watchtools`**. Les lanceurs `came.sh` /
`came.bat` n'ont pas ete repris.

Depuis la racine de WatchTools :

    python main.py                   l'accueil : une icone par outil, chacun dans son processus
    python main.py came              cet outil, sans passer par l'accueil
    python -m came.main              cet outil seul
    python -m came.main --chrono     calcul en console
    python came/verifier.py          versions + calcul d'essai, ne lance rien

Depuis CE dossier (l'utilisateur ouvre un VSCode et une session Claude par
outil pour travailler sur plusieurs en parallele ; regles communes dans le
`CLAUDE.md` de la racine) :

    python main.py [--chrono]        main.py met lui-meme la racine dans sys.path
    python verifier.py
    python -m came.<module>          dans le terminal VSCode (PYTHONPATH=.. via .vscode/settings.json)

`.vscode/launch.json` : « Came », « module courant », « verifier
l'installation », « synthese en console ».

Le theme (couleurs, feuille de style) et les widgets de panneau
(`OngletsCompacts`, `Etiquette`, `soft`, `aide`) sont dans `commun/`, partages
avec les autres outils.

Dependances : **PySide6** (obligatoire), numpy + numba (acceleration), ezdxf
(export DXF). Qt est **PySide6** et non PyQt6 — choix de licence, pour
diffuser en LGPL.

### Tester l'interface sans ecran

Recette qui marche, a reprendre telle quelle :

```python
import os; os.environ['QT_QPA_PLATFORM'] = 'offscreen'
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox
b = []
for n in ('critical', 'information', 'warning', 'question'):
    setattr(QMessageBox, n, staticmethod(lambda *a, **k: b.append(a[1])))
QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: ('fichier.dxf', ''))
QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: ('sortie.dxf', ''))
```

**Neutraliser les quatre boites est indispensable** : une seule non stubbee
bloque le test indefiniment, sans message.

**Attendre reellement.** La synthese tourne dans un `QThread` et prend 7 a
13 s. Une boucle de `processEvents()` sans `sleep` ne suffit pas :

```python
def att(s):
    t = time.time()
    while time.time() - t < s:
        app.processEvents(); time.sleep(0.03)
```

Ne jamais conclure « le coeur ne se construit pas » sans avoir attendu.

---

## 3. Architecture

| fichier | role |
|---|---|
| `main.py` | CLI : `--chrono` (= `--cli`), `--nsga`, `--marteau`, `--face`, `--figer`, `--procs`, `--cellules`, `--population`, `--generations`, `--scan` |
| `core.py` | `Params`, `FlatFace`/`ArcFace`, `contour_plan()`, `synthesize()` (came simple, conservee mais hors flux), `pose()` |
| `chrono.py` | `ChronoParams`, `synthesize_chrono()` : deux flancs, plat, raccords, pointe. `support()` donne (bet, Hs, ss) |
| `profil.py` | `Profil` (courbure a bande limitee **+ plat explicite**), `evaluer()`, `_ouverture_nez()`, `dimensionner_ressort()`, `longueur_plat()`, `balayage_marteau()`, `ProfilImporte` |
| `noyau.py` | noyaux numba + repli Python pur, valides l'un contre l'autre. `DISPONIBLE` |
| `nsga2.py` | NSGA-II ecrit a la main, generique : `evaluer_x(x) -> (objectifs, violation)`. `MARTEAU_LIBRE`, `MARTEAU_UTILE` |
| `arcs.py` | approximation par arcs tangents, `mesurer_coeur()` |
| `dxfio.py` | lecture/ecriture DXF (ezdxf), `aligner_contour()` |
| `projet.py` | etudes `.json` |
| `ui.py` | la fenetre, ~2700 lignes |
| `verifier.py` | controle d'installation |

---

## 4. Conventions du modele — a ne pas melanger

**Fonction d'appui.** Le coeur est decrit par `H(beta)`, distance de l'axe a
la tangente de normale sortante `beta`. Alors `s = dH/dbeta` est l'abscisse
du point de contact, `rho = H + H''` le rayon de courbure, et le point de
contact vaut `X = H.u(beta) + s.t(beta)`.

**Couple avec frottement**, identique pour un contact sur le plan ou sur un
arrondi — seules les deux grandeurs geometriques changent :

    q = s - t(gamma).P        e = H - u(gamma).P
    T = Cm (s + sigma.mu.H) / (q + sigma.mu.e)
    N = |Cm / (q + sigma.mu.e)|

`sigma = +/-1` est le sens de glissement : les deux sont essayes, on retient
celui qui rend le mouvement possible ; si aucun ne le rend possible la
position est **arc-boutee** et compte dans la plage morte.

**Pente du pan.** `angle_plan` est mesure **de la direction O->P a la normale
du plan, au repos**. Elle est donc independante de l'armage. D'ou
`gamma0 = atan2(P) + radians(angle_plan)`.

**Le plat du zero est un Dirac de courbure explicite**, pas une serie. Une
serie a bande limitee n'en rend que 0,31 mm sur 1,46, et 0,67 encore a 96
harmoniques. La fonction d'appui qui le porte est

    H_plat(beta) = L.sin(x/2),   x = (beta - beta0) mod 2.pi

`s` y saute de `L` — c'est exactement un segment droit. **Le demi-argument
est essentiel** : avec un argument entier, `cos(0) = cos(2.pi)` et rien ne
saute. La courbure lissee ajoutee, `(3L/4).sin(x/2)`, est positive ; et
`H_plat` etant periodique, sa courbure a un rang 1 nul, donc le contour se
referme sans liberer le rang 1 de la serie.

**`beta0` doit valoir `gamma0`.** Changer la pente ou le pivot deplace la
normale de repos : le plat du coeur **doit suivre**, sinon le coeur ne part
plus plat contre plat. Vrai pour l'edition manuelle, pour les solutions du
front et pour le balayage.

**Non-penetration du nez du marteau.** Le nez est le segment des deux centres
de coin dilate du rayon `rho`. Le coeur, vu du contact, est son enveloppe
convexe, donc une intersection de demi-plans. Un segment en est hors **si et
seulement si** un appui separe ses **deux** extremites :

    min(u_j.A, u_j.B) >= H_j + rho

L'ouverture cherchee est le minimum sur j. Resolu en forme close.

**L'encombrement est une limite, pas un but.** `r_max` est mesure sur tout le
profil. Le rapport annonce celui qui est atteint et la monotonie du rayon.

---

## 5. Erreurs deja commises — ne pas les refaire

Chacune a coute une ou plusieurs iterations. Symptome, puis cause.

**« Le marteau quitte le coeur. »** Quatre causes distinctes, cumulees :
1. la branche « coin » acceptait un `phi` de signe oppose a l'ouverture, donc
   le marteau passait de l'autre cote du coeur ;
2. le contact arrondi/facette etait structurellement manque — sur une facette
   tous les points partagent **une** normale, donc la portion droite de la
   courbe offset se reduit a un echantillon et la racine passe a cote ;
3. le bras de levier etait pris sur un sommet quelconque de la facette, alors
   qu'un arrondi touche au **pied de la perpendiculaire** ;
4. prendre, par caracteristique, la premiere ouverture qui la degage puis le
   **maximum** des trois est faux : l'ensemble admissible d'un coin est un
   **arc**, pas un intervalle croissant. D'ou le critere unique du §4.

**Le centre de l'arrondi etait pose sur le plan** au lieu d'etre rentre de
`rc` dans la matiere : l'arc depassait devant le plan et le coin touchait a
chaque position.

**La plage morte se comptait en echantillons.** `max(runs) * 360/n` donnait
3 deg a n=120 et 1 deg a n=360 — un critere d'optimisation qui variait avec
la finesse du calcul. Elle est maintenant mesuree par interpolation de la
marge `|s| - mu.H`.

**Le plancher `2.atan(mu)` n'est pas incontournable.** Il ne vaut que pour une
**pointe vive en appui plan**. Sur l'arrondi du marteau, le bras de levier ne
s'annule plus : 14,1 deg avec un pan infini, **1,05 deg** avec un pan de
1,4 mm. C'est le levier principal, et l'utilisateur l'avait vu avant moi.

**`reparer` brouillait le vecteur.** `profil.reparer` attend le vecteur
**complet**, cx et cy compris. Lui passer `x[:-nq]` les retirait deux fois :
les harmoniques se decalaient d'un rang et la graine ressortait irrecevable.

**Un depassement de 0,9 um disqualifiait la graine.** En domination
contrainte, toute solution recevable battait alors la meilleure.
`TOL_ENCOMBREMENT = 0,005 mm`.

**Les processus d'evaluation reconstruisent le probleme.** Il faut leur
transmettre `plat_mini`, `gabarit` et `marteau`, sinon ils evaluent un autre
coeur et toute la generation sort irrecevable (« 0 faisables »).

**`pivot_t` s'accumulait.** C'est un **deplacement** : il doit etre compte
depuis le pivot du **lancement**, pas depuis le pivot courant. Sinon il se
reporte a chaque clic et le pivot tourne autour du coeur sans revenir.

**L'etude enregistrait le mauvais coeur** : elle prenait le profil de l'onglet
actif. Le sujet est le contour importe ou optimise, quel que soit l'affichage.

**Qt : deux pieges de mise en page**, tous deux responsables de centaines de
pixels de vide :
1. un `QLabel` replie en politique horizontale `Ignored` annonce une hauteur
   **mini** calculee a largeur minimale — 280 px pour 110 occupes. D'ou
   `Etiquette(QLabel)` qui redefinit `minimumSizeHint`;
2. `heightForWidth` d'un `QTabWidget` porte sur **toutes** ses pages.

**PySide6 n'appelle pas une methode virtuelle posee sur l'instance.**
`lbl.minimumSizeHint = lambda: ...` marchait en PyQt6, pas ici : il faut une
sous-classe.

**Un DXF sans table de calques** est mal relu par les logiciels de CAO : il en
sort des points et des droites, pas un coeur. D'ou ezdxf.

**Le zoom prenait le centre de la fenetre pour origine.** `_m` a pour point
fixe `center`, projete au milieu de la fenetre decale de `offset` : ne
multiplier que `scale` zoomait donc sur un point qui, apres un deplacement de
la vue, n'est pas forcement sur la partie affichee. Il faut garder le point du
dessin **sous le curseur** :

    offset' = a - (a - offset) . k        a = souris - milieu de la fenetre

exact et symetrique en x et y malgre l'axe des ordonnees inverse. Verifie :
l'ecart du point sous le curseur reste a 2e-15 mm, y compris apres
deplacement de la vue. Le zoom est borne a [0,5 ; 20000] px/mm et la molette
horizontale est ignoree.

**Piege Batch** dans `came.bat` : une variable posee et lue dans le meme bloc
entre parentheses est developpee **avant** d'etre ecrite. D'ou les
sous-routines appelees par `call`.

---

## 6. A verifier apres toute modification

Le script de regression vit dans le scratchpad de session ; le recreer au
besoin. Ce qu'il doit couvrir :

- **invariants physiques** sur tous les chemins (profil importe, reconstruit,
  synthese) : contact effectif a chaque position, marteau rigide, un seul cote
  du coeur, effort normal fini, au plus deux changements de signe du couple,
  continuite de `phi` ;
- **accord numba / Python pur** : il doit rester a ~2e-15 sur phi, N et T ;
- **reference par force brute** : fermer le marteau par dichotomie jusqu'a la
  premiere intersection de son contour reel avec le polygone du coeur, et
  comparer l'ouverture. L'ecart doit rester au millieme de degre. **Tester
  l'intersection des segments, pas seulement les sommets** — un plan peut
  traverser le coeur sans qu'aucune extremite soit dedans ;
- `./came.sh --verifier`, `main.py --chrono`, `main.py --nsga`, et l'interface
  complete : construction, simulation, arcs, import/export DXF, aller-retour
  JSON, balayage, NSGA, selection de solutions.

Valeurs de reference actuelles (ChronoParams par defaut, pan 1,4 mm) :

    plat porte 1,4541 mm (impose 1,4613)   plat mesure 1,5589
    ressort C0 = 6,2993 N.mm   k = 8,5963 N.mm/rad   residu 0,7 %
    couple moyen 0,28999 pour 0,3000 vise
    plage morte 1,05 deg

---

## 7. Ce qui reste ouvert

**La synthese suppose le pan infini.** Elle cherche la tangence sans borne
d'abscisse. La phase 2 montre apres coup ce que le pan reel delivre, elle ne
l'anticipe pas. Un coeur dimensionne **pour** un pan court demanderait de
reprendre l'integration inverse avec la borne et les arrondis.

**NSGA-II plafonne la ou un balayage direct fait mieux.** A 32 individus sur
30 generations il atteint 23,9 % d'ondulation, quand la grille
`balayage_marteau()` en trouve 12,5 en 2,6 s. C'est normal : il repartit son
effort sur quatre objectifs et sur tout l'espace des formes. Ne pas vendre
NSGA-II comme la bonne reponse a un probleme a deux variables.

**Le levier du pivot est sa position angulaire, pas sa distance.** Mesure sur
+/- 4 mm : radialement l'ondulation varie de 1,6 point, perpendiculairement a
O->P de 68,8 points. Facteur 43. D'ou la variable `pivot_t`.

---

## 8. Le projet Spiral — ce que son solveur est vraiment

Lu dans `Spiral_Lab_Theorie_et_Structure.docx`. **A savoir avant d'envisager
de le reutiliser.**

Ce n'est **pas** un solveur d'elements finis 2-D. C'est une **elastica
non lineaire 1-D** : poutre d'Euler-Bernoulli inextensible et non cisaillable
en grandes rotations, resolue par **methode de tir** — trois inconnues
`(Fx, Fy, M0)` au piton, trois conditions au collet, systeme exactement
determine, resolu par `scipy.optimize.root(method='hybr')` (Powell hybride).
L'integration est un **predicteur-correcteur** compile Numba, en une passe du
piton vers le collet, qui resout le couplage implicite entre le moment et la
position deformee.

Deux enseignements de methode, transposables :

- resoudre le **systeme** plutot que minimiser un cout agrege : le passage de
  BFGS a Powell a fait tomber le residu angulaire de 0,24 deg a 1e-9 deg, et
  l'erreur sur l'energie elastique de 2,4e-3 — dix fois l'effet
  d'isochronisme qu'on mesurait dessus. Quatre fois plus rapide, en plus ;
- la quadrature de la periode est regularisee par `theta = a.sin u`, choisie
  pour etre **exacte sur le cas lineaire** : sans cela son propre biais
  noierait un signal de 1e-4.

**Conclusion pour la simulation 2-D demandee** : le moteur du spiral ne se
branche pas. Physique differente (poutre 1-D contre milieu continu 2-D),
dimensionnalite differente, et il ne traite aucun contact. Ce qui se reprend
de `spiral_lab`, en revanche : l'architecture (fenetre PySide6 unique,
panneaux + onglets, calculs en fil de fond), l'import DXF qui reconstruit une
fibre neutre et une epaisseur locale depuis deux flancs, la configuration JSON
a compatibilite ascendante, le rapport Word, et l'etude de convergence
automatisee.

Pour la chaine « DXF -> maillage -> conditions aux limites -> contraintes »,
mon avis : **Gmsh** pour le maillage (pip, API Python, maillage de surface
contraint), un solveur elastique lineaire a ecrire, et le **contact traite en
dernier** comme un sujet a part entiere — c'est un probleme non lineaire et
non regulier, et c'est precisement la ou j'ai fait le plus d'erreurs ici,
avec un seul point de contact et une solution analytique.

---

## 9. Facon de travailler attendue

- **Mesurer, pas affirmer.** Chaque correction est verifiee par un chiffre, et
  le chiffre est donne. Plusieurs de mes erreurs venaient d'un script de test
  faux, pas du programme : verifier le test avant d'accuser le code.
- **Dire ce qui echoue.** Une tentative annulee se raconte, avec sa cause. Le
  plat explicite a demande deux essais faux avant le bon.
- **L'utilisateur connait son metier.** Il a eu raison sur la plage morte, sur
  la lecture du DXF, sur l'incoherence des parametres et sur le pivot — a
  chaque fois avant moi. Ses remarques sont des mesures, pas des impressions.
- **Ne pas elargir le perimetre.** La came simple et le galet sont abandonnes.
- Commentaires et libelles **en francais sans accents** dans le code, pour
  rester lisible partout.
