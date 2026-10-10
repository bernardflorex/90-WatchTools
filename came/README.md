# Dimensionnement d'un coeur de chronographe

L'application traite le **coeur de chronographe**, d'un bout a l'autre, dans
un seul panneau. La came simple est remise a plus tard, et le galet est
abandonne.

## Le pivot derivait a chaque solution selectionnee

Signale : le pivot change de cadran et ne revient plus en place quand on
remonte dans les solutions, le depart n'est plus plat sur plat, et on ne sait
plus si la simulation est correcte. Deux defauts, tous deux reproduits.

**`pivot_t` s'accumulait.** C'est un deplacement, et il etait compte depuis le
pivot **courant** au lieu du pivot avec lequel l'optimisation a ete lancee.
Chaque clic le reportait donc sur le precedent :

| clic | pivot | \|OP\| | ecart plat / plan |
|---|---|---|---|
| depart | (16.00, 5.00) | 16.76 | 0.000 deg |
| ligne 0 | (15.60, 6.28) | 16.82 | 4.574 deg |
| ligne 1 | (15.10, 7.53) | 16.87 | 9.168 deg |
| ligne 0 | (14.50, 8.73) | 16.92 | 13.713 deg |
| ligne 1 | (13.80, 9.89) | 16.98 | 18.279 deg |
| ligne 0 | (13.02, 10.98) | 17.03 | 22.795 deg |

La distance a l'axe restant constante, le pivot **tournait** autour du coeur -
exactement « change de cadran ». La selection est maintenant idempotente :
ligne 0 redonne toujours (15.60, 6.28), ligne 1 toujours (15.60, 6.29).

**Le plat du coeur ne suivait pas.** Le profil gardait la normale de repos du
lancement, alors que la solution deplace la pente ou le pivot : d'ou l'ecart
qui grandit dans la colonne de droite. Le gabarit est maintenant reconstruit
avec la nouvelle normale, et l'ecart revient a 0.000 deg.

### « Je ne sais pas si la simulation est correcte »

Deux reponses. D'abord le controle : pour chaque solution, les quatre
objectifs affiches ont ete compares a ceux que l'optimiseur avait evalues, sur
trois lignes du front - ils coincident **chiffre pour chiffre** (48.333 /
3.677 / 62.396 / 0.713, etc.).

Ensuite, c'est desormais visible dans le panneau. Les *Reperes* portent une
ligne **depart de la simulation** : en vert « plat sur plan (0.004 deg) », en
rouge « decale de X deg » au-dela de 0.2 deg. Elle se met a jour apres une
construction, apres chaque solution selectionnee et apres un balayage.

## Le pivot : une seule direction porte

Hypothese posee : deplacer le pivot en **distance** n'apporte pas grand
chose, le faire varier en x et y dans un carre en apporterait plus. Verifiee,
et avec une marge nette. Coeur fixe, pan de 1.4 mm, ressort redimensionne a
chaque position pour que le couple vise soit toujours atteint :

| direction, sur +/- 4 mm | ondulation | etendue |
|---|---|---|
| radiale (change \|OP\|) | 41.3 a 42.9 % | **1.6 point** |
| perpendiculaire a O->P | 13.5 a 82.3 % | **68.8 points** |

Un facteur **43** entre les deux. La distance a l'axe n'apporte quasiment
rien ; c'est la position **angulaire** du pivot qui compte.

### La grille carree le montre directement

Balayage de (x, y) par pas de 2 mm, pente relative tenue a 60 deg,
ondulation en % :

```
          x= 12.0  x= 14.0  x= 16.0  x= 18.0  x= 20.0
  y=  1.0     12.9     12.5     12.5     12.7     13.3
  y=  3.0     34.5     26.0     24.9     22.9     18.9
  y=  5.0     54.3     47.5     42.0     36.6     33.9
  y=  7.0     80.3     69.3     60.2     52.0     47.8
  y=  9.0    104.6     89.6     77.2     68.1     60.5
```

Les lignes sont plates, les colonnes non : +/- 4 mm en x ne change rien
(12.9 contre 13.3 a y = 1), -4 mm en y fait passer l'ondulation de 42.0 a
12.5 % et le residu du ressort de 11.0 a 2.7 %. Ici x est a peu pres la
direction radiale et y la perpendiculaire : le tableau dit la meme chose que
la mesure ci-dessus.

### Pivot et pente sont couples

La pente seule, a pivot fixe, descend l'ondulation de 42.0 a 28.5 % (optimum
50 deg). Le couple (position du pivot, pente) descend a **12.5 %**. La carte
montre une vallee diagonale : l'optimum de pente se deplace avec la position
du pivot - 45 deg a y = 5, 60 deg a y = 1, 65 deg a y = 0. Les deux variables
ne sont donc pas interchangeables.

### Ce qui en a ete tire

`pivot_t` remplace `pivot_x` et `pivot_y` dans les variables conseillees :
c'est le deplacement du pivot **perpendiculairement a O->P**, une seule
variable au lieu de deux dont une inutile. X et Y restent disponibles pour qui
veut le carre complet.

Les directions utiles sont desormais **semees** dans la population de depart,
comme l'etait la courbure de face : sans cela, dix generations plafonnaient a
35 % d'ondulation. Avec le semis, 23 %.

### Et un constat qui vaut d'etre dit

Meme a 32 individus sur 30 generations, NSGA-II plafonne a **23.9 %**
d'ondulation la ou un balayage direct de ces deux nombres en trouve **12.5**.
C'est normal : il repartit son effort sur quatre objectifs et sur tout
l'espace des formes de coeur, alors qu'on cherche ici deux valeurs.

D'ou le bouton **Balayer pivot x pente** : `profil.balayage_marteau()` parcourt
la grille a coeur fixe, en 2.6 s pour 61 points, applique le meilleur point et
l'explique dans le rapport. Le critere se choisit - ondulation, effort, plage
morte, tour sous la consigne - et chacun donne un point different :

| critere | pivot retenu | pente | ondulation | effort | plage morte |
|---|---|---|---|---|---|
| depart | (16.00, 5.00) | 60.0 | 38.2 % | 0.647 N | 1.05 deg |
| ondulation | (17.50, 0.20) | 60.0 | **15.6 %** | 0.608 N | 1.01 deg |
| effort | (17.50, 0.20) | 53.3 | 27.2 % | **0.581 N** | 1.06 deg |

### Un defaut trouve en chemin

Le profil porte son plat a la normale de repos `gamma0`. Changer la pente ou
le pivot la deplace, et le profil ne suivait pas : le balayage annoncait 15.4 %
d'ondulation, la simulation en affichait 38.8 apres application. Le meme
defaut touchait l'edition **manuelle** de la pente ou du pivot sur un coeur
synthetise. Le profil est maintenant reconstruit avec la nouvelle normale, et
les deux chiffres concordent a 0.2 point.

## Les leviers de l'optimisation

Question posee : la plage morte et l'effort de contact sont des objectifs,
mais sur quoi l'optimisation a-t-elle prise ? Elle ne faisait varier que la
courbure du coeur, harmonique par harmonique - or le profil du coeur n'a
presque **aucune** prise sur ces deux grandeurs. Les leviers sont la longueur
du pan, ses arrondis, sa pente et la position du pivot. Les laisser fixes
revenait a demander d'ameliorer ce qui ne pouvait pas bouger.

Six grandeurs du marteau peuvent desormais entrer dans les variables, case par
case : longueur du pan, decalage, rayons de coin, pente du pan, pivot X et
pivot Y. Les bornes sont une fenetre autour de la valeur de depart - pan de
0.15 a 1.4 fois l'encombrement, pente a +/- 28 deg, pivot a +/- 25 % de sa
distance a l'axe - pour rester dans un dimensionnement voisin plutot que
d'explorer des geometries ou la synthese n'a plus de sens.

La pente et le pivot deplacent la normale de repos : le **plat du coeur suit**,
sinon le coeur ne partirait plus plat contre plat. Et la solution retenue
reporte ses valeurs de marteau dans les champs du panneau, sinon la simulation
affichee ne serait pas celle qui a ete evaluee.

### Ce que cela apporte, mesure

Meme graine, meme budget (28 individus, 26 generations), gain sur chaque
objectif par rapport a la graine :

| variables | sous consigne | plage morte | ondulation | effort |
|---|---|---|---|---|
| coeur seul | -27.5 % | -0.50 deg | -3.92 % | -0.19 N |
| coeur + pan | -23.3 % | -0.17 deg | **-6.18 %** | -0.12 N |
| coeur + pan + pivot | -24.2 % | -0.43 deg | -3.58 % | -0.11 N |

La lecture honnete : chaque case ajoute une variable, donc agrandit l'espace
de recherche, et a budget egal la convergence est plus lente. L'**ondulation**
y gagne tout de suite - c'est le seul objectif ou le marteau apporte plus que
le coeur. Pour les autres, il faut plus de generations avant que le gain
compense le cout.

Autrement dit : les leviers existent maintenant, mais ils ne sont pas
gratuits. Ouvrir les six cases d'emblee n'est pas la bonne facon de s'en
servir.

## Le paragraphe Reperes annoncait un plancher, pas une mesure

Sa ligne « plage morte au sommet » donnait `2.atan(mu)`, le plancher d'une
pointe vive en appui plan. On lisait donc 13.7 deg quand la simulation en
donnait 1.05, et la valeur ne bougeait pas en raccourcissant le pan : ce
n'etait pas la meme grandeur. Les deux sont desormais distinguees, et la
mesure suit la simulation - 14.14 deg avec un pan infini, 1.05 a 1.4 mm de
pan, 1.50 a 2.5 mm.

## Le coeur construit se simule tout de suite

Le profil discretise etait le produit d'un **import** ou d'une
**optimisation** : il fallait donc lancer NSGA-II pour connaitre le couple
reellement delivre, la plage morte, les contacts sur coin ou l'approximation
par arcs d'un coeur que l'on venait de construire. Il est maintenant deduit de
la synthese, avec son plat exact, des que le coeur est calcule.

C'est la **simulation** qui est affichee, et non la synthese. Celle-ci donne
un couple constant par construction - une droite a la consigne - et une plage
morte de pan infini : les courbes montraient donc 0.300 N.mm tout du long et
13 deg de plage morte, et le curseur placait le plat contre plat **au milieu**
de la course, puisque sur cet onglet l'abscisse allait de `-course_b` a
`+course_a`.

Le trace, les courbes et le curseur passent donc sur le profil discretise et
son pan reel. `theta = 0` y est la position de repos, plat contre plat, au
**debut** de la course, et le couple trace est celui que le marteau delivre
vraiment : de -0.315 a +0.379 N.mm sur l'exemple ci-dessous, au lieu d'une
droite.

Le rapport porte un bloc *Simulation du coeur construit*. Exemple, pan de
1.4 mm avec des coins de 0.1 mm :

    plat du zero (a 2 um pres)       1.5589 mm
    couple moyen sur le coeur       0.30361 N.mm (vise 0.3000)
    ondulation                        38.19 %
    effort de contact maxi           0.6465 N
    plage morte                        1.05 deg
    contacts sur un coin du pan         220 / 360
    residu de l'ajustement affine      11.1 %

La plage morte affichee tombe ainsi a **1.05 deg** au lieu de 14.05, et les
trois champs du pan la font bouger : 1.07 deg a 2.2 mm de pan, 1.36 a 0.2 mm
de rayon de coin, 1.38 avec 0.3 mm de decalage.

L'approximation par arcs, l'export DXF du profil et le dimensionnement du
ressort suivent de la meme facon, sans optimisation. Verifie : 20 arcs et un
DXF sans erreur d'audit sur un coeur fraichement construit.

Le rapport garde les blocs de construction - courses, encombrement, couples
requis, plat impose par la geometrie - au-dessus de la simulation : on voit
donc d'un coup ce que la synthese a decide et ce que la piece delivre.

`is_optim` ne veut plus dire « un profil existe » mais « un profil est
**affiche** » : la vue du chronographe reste a l'ecran apres une
construction, parce qu'elle montre le plat, la corde et le creux de pointe,
alors qu'un import ou une solution du front basculent sur la vue du contour.

## Le galet a disparu du panneau

Ses deux grandeurs - rayon d'axe et frottement de palier - restaient
affichees alors que le galet est abandonne et que plus rien ne les lisait.
Elles sortent des champs. `profil.mu_galet` demeure, documente, si le sujet
revient un jour.

## Le plat du zero est maintenant exact

C'etait le verrou : une serie a bande limitee ne peut pas porter un plat - il
demande un Dirac de courbure - et la graine n'en rendait que 0.31 mm sur 1.46.
Le coeur optimise perdait donc sa position nulle.

### La reparametrisation

La courbure n'est plus une seule serie mais **un Dirac plus une serie** :

    rho(beta) = L.delta(beta - beta0) + rho_serie(beta)

La fonction d'appui qui porte le Dirac, periodique et continue, est

    H_plat(beta) = L.sin(x/2),   x = (beta - beta0) mod 2.pi

`s = dH/dbeta` y vaut `(L/2).cos(x/2)` et saute donc de `-L/2` a `+L/2`, soit
`L` : c'est exactement un segment droit de longueur L. Un argument entier ne
sauterait pas, `cos(0) = cos(2.pi)` - c'est l'erreur que j'avais faite au
premier essai. La courbure lissee qu'elle ajoute, `(3L/4).sin(x/2)`, est
positive, donc convexe ; et comme `H_plat` est periodique, sa courbure a un
rang 1 nul : le contour se referme, sans avoir a liberer le rang 1 de la
serie.

### Trois consequences dans le code

**La contrainte d'usinabilite porte sur la courbure totale.** La serie doit
rester au-dessus de `mini - (3L/4).sin(x/2)`, un plancher qui depend de beta.
`reparer` et `noyau.reparer_coeffs` recoivent donc un plancher **par
echantillon** ; un plancher scalaire gonflait la serie la ou le plat apportait
deja de la courbure, et le profil sortait a 0.24..1.41 mm au lieu de
0.30..3.00.

**Le Dirac doit etre retire des donnees, pas de leur interpolation.** Dans la
fonction d'appui analytique, le plat est un saut de `s` repere a
`beta = gamma0` et de hauteur exactement la longueur du plat. La
discontinuite est repartie sur un intervalle d'echantillonnage : une
soustraction analytique laissait une pointe de `L/dbeta`, plus de deux cents,
dans une cellule de courbure. On soustrait donc le plat des echantillons
eux-memes.

**Le centre se cale sur la cible d'origine.** Le caler contre la cible privee
de son plat laissait un rang 1 de 0.6 mm dans l'ecart, c'est-a-dire un centre
de rotation deplace d'autant.

### Et deux dans l'optimisation

Le plat est un **gabarit**, pas une variable : il est repris de la synthese et
transmis a `probleme`, de sorte que toutes les solutions gardent le meme plat,
exact. Il faut aussi le transmettre aux **processus d'evaluation** : sans cela
les fils reconstruisaient le probleme sans le plat, evaluaient donc un autre
coeur, et toute la generation sortait irrecevable (« 0 faisables »).

### Ce que cela donne

| | avant | apres |
|---|---|---|
| plat rendu (impose 1.46 mm) | 0.31 mm | **1.56 mm** |
| ecart de la fonction d'appui | - | 0.029 mm |
| residu de l'ajustement du ressort | 103 % d'ondulation | **0.7 %** |
| couple moyen pour 0.300 vise | - | 0.289 N.mm |

Et la graine, qui est desormais un vrai coeur a plat, part de bien plus haut :
54.2 % du tour sous la consigne au lieu de 62.8, 66.6 % d'ondulation au lieu
de 103.3. Le front la contient et l'ameliore sur les quatre objectifs, en
gardant le plat : la solution affichee en mesure 1.5611 mm.

## Qt : PySide6, et non PyQt6

L'interface est passee a **PySide6** (6.11.2). Les deux sont le meme Qt ; la
difference est la licence : PyQt6 est GPL ou commercial, PySide6 est LGPL -
donc diffusable en logiciel libre sans contrainte sur ce qui l'utilise.

Le portage a ete court : huit `pyqtSignal` devenus `Signal`, trois lignes
d'import. Les enumerations qualifiees (`Qt.PenStyle.DashLine`,
`QSizePolicy.Policy.Ignored`) marchent telles quelles.

**Un seul idiome n'a pas survecu**, et c'etait le correctif de mise en page le
plus important : poser `minimumSizeHint` sur l'**instance** d'un label.

    lbl.minimumSizeHint = lambda: QSize(0, 0)      # PyQt6 seulement

PySide6 n'appelle pas une methode virtuelle attachee ainsi : il faut une
sous-classe. D'ou `Etiquette(QLabel)`, qui redefinit `minimumSizeHint`
proprement. Sans elle, les centaines de pixels orphelins du panneau
revenaient. Verification apres portage : les ecarts entre blocs valent tous
6 px, l'espacement normal.

Tout le reste a ete repasse sous PySide6 : construction du coeur, simulation,
arcs, import et export DXF, aller-retour JSON, balayage, NSGA-II et selection
des solutions.

## Lancer le programme

Hors VSCode, sans rien connaitre de Python :

| systeme | fichier | double-clic |
|---|---|---|
| Linux, Kubuntu | `lancer.sh` | oui, apres `chmod +x lancer.sh` |
| Windows | `lancer.bat` | oui |

Au premier demarrage le lanceur cree l'environnement et installe les
dependances de `requirements.txt`. Ensuite il le retrouve seul.

    ./lancer.sh                 lance l'interface
    ./lancer.sh --verifier      controle l'installation, ne lance rien
    ./lancer.sh --raccourci     pose une entree dans le menu KDE / GNOME
    ./lancer.sh --chrono        calcul en console

L'environnement est choisi dans cet ordre : la variable `CAME_VENV`, puis
`.venv` ou `.came` a cote du programme, puis `~/DEV/virtual_env/.came` et
`~/.venvs/came`. Un environnement **complet** passe avant un environnement
incomplet - sans cette regle, un `.venv` de projet auquel il manque une
dependance etait choisi avant celui de l'utilisateur, qui les a toutes. Un
environnement incomplet est complete au lieu d'etre rejete.

`--verifier` dit ce qui repond et ce qui manque, puis fait un calcul d'essai :

    environnement : /home/florian/DEV/virtual_env/.came
    python        : 3.14.4
    PyQt6         : 6.11.0
    numpy         : 2.5.3
    numba         : 0.68.0
    ezdxf         : 1.4.4
    noyau compile : oui
    essai de calcul : ok, couple 0.30000 N.mm, plat 1.4613 mm

Le chemin de premiere installation a ete essaye a blanc, dans un dossier
vierge et sans environnement joignable : creation, installation, lancement.
Le `.bat` n'a pas pu etre execute ici, faute de Windows ; sa structure a ete
controlee ligne a ligne. Il evite notamment le piege classique de Batch - une
variable posee et lue dans le meme bloc entre parentheses est developpee
avant d'etre ecrite - en passant par des sous-routines appelees par `call`.

## Les deux phases

**1. Le coeur, pan de marteau suppose infini.** On donne le coefficient de
frottement, la distance du plat a l'axe (nominal 0.6 mm), l'encombrement maxi,
le couple de rappel voulu, le pivot du marteau et la pente de son pan. La
synthese inverse en deduit la forme du coeur et le **ressort** : prearmage et
raideur sont des sorties, calculees pour que le couple voulu soit atteint.
Elles s'affichent dans le panneau, elles ne se saisissent pas.

**2. Le pan reel.** Le plat du coeur etant obtenu, on pose le pan tel qu'il
sera : longueur, decalage le long de lui-meme, rayons de coin. Le couple du
marteau est celui de la phase 1 ; on voit alors le couple que le **vrai**
marteau delivre sur l'axe du coeur. Exemple mesure : couple voulu 0.300 N.mm,
couple moyen delivre par un pan de 1.4 mm avec des coins de 0.1 mm :
0.260 N.mm, plage morte 1.79 deg au lieu de 14.05.

Tout champ agit et le resultat suit. La phase 1 reconstruit le coeur - quelques
secondes, en tache de fond, avec une temporisation pour ne pas relancer le
calcul a chaque cran de molette - la phase 2 ne refait que la simulation.

## L'encombrement est une limite, pas un but

Le rayon maxi est mesure sur **tout** le profil, pas a la pointe, et le
rapport annonce celui qui est reellement atteint pour celui qui est admis. Il
dit aussi si le rayon est monotone croissant le long des flancs : un coeur non
monotone peut tenir dans moins de place que son rayon de pointe ne le
suggererait.

## Repartir d'un dessin

A l'import d'un DXF, les grandeurs de construction sont **mesurees sur le
contour** et portees dans les champs : distance du plat a l'axe, encombrement,
longueur du plat, rayon des raccords concaves, rayon et corde de la pointe, et
bien sur le pivot, la pente et le pan du marteau. Les mesures viennent de
l'approximation par arcs du dessin, elles sont donc approchees et le journal le
dit.

Deux grandeurs ne sont dans aucun contour : le **frottement**, qui est une
donnee physique, et le **couple de rappel voulu**, qui est une specification.
Leurs libelles sont mis en evidence dans le panneau, avec l'explication en
infobulle. Verification sur `Sketches - minute.dxf` : distance plat-axe
mesuree a 0.601 mm pour un nominal de 0.6.

## L'optimisation est passee en second plan

Elle reste disponible, avec son avertissement : tant que la representation a
bande limitee ne porte pas un plat exact, le front n'est pas exploitable pour
un coeur a plat. Voir plus bas pourquoi, et ce qu'il faudrait pour le corriger.

## Sauvegarde

Les etudes `.json` portent les parametres, le contour, le ressort calcule et
les reglages ; l'export DXF est ecrit par ezdxf, avec un seul marteau a la
position de repos par defaut.

---

# Dimensionnement de cames horlogeres

Deux familles de cames sont prevues ; la seconde reste a faire.

1. **Coeur simple** : remise en position d'un mobile par un marteau, a couple
   constant. Onglet *Coeur simple*.
2. **Coeur de chronographe** : coeur a deux flancs, remise a zero dans un sens
   ou dans l'autre, plat du zero et rayons de raccord. Onglet *Coeur de
   chronographe*.
3. *(a venir)* cames de commande imposant un mouvement defini a un palpeur.

## Lancement

```bash
python3 main.py              # interface graphique (PyQt6)
python3 main.py --cli        # synthese en console
python3 main.py --cli --angle-plan 65 --mu 0.15 --export coeur.dxf
```

## Installation

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Seule l'interface graphique a une dependance (PyQt6) : le calcul est en Python
pur. Le mode `--cli` fonctionne donc sans rien installer.

## Modele

Le coeur tourne autour de O = (0,0), le marteau autour de P. Le marteau porte un
plan (ou un arc) de contact et recoit un couple `Cm = C0 + k.|psi|` : `C0` est le
couple mini (coeur en place), `k` la **pente du marteau**, `psi` son ouverture —
ce qui remplace un ressort de rappel separe.

Le coeur est decrit par sa **fonction d'appui** `H(beta)` (distance de O a la
tangente de normale sortante). Le point de contact a pour abscisse tangentielle
`s = dH/dbeta`. Le bilan des moments au contact donne, avec frottement :

```
N  = Cm / (q + sigma.mu.e)                q = s - t.P   (bras de N au pivot)
Tc = Cm (s + sigma.mu.H) / (q + sigma.mu.e)   e = H - u.P   (bras du frottement)
```

`sigma = -signe(v_glissement)` distingue le frottement **rentrant** (retour) du
frottement **sortant** (armement) : les deux sont calcules et traces.

**Synthese inverse** : on impose `Tc = constante`, ce qui donne `s` explicitement,
puis la coherence `dH/dtheta = s.(dbeta/dtheta)` fournit une EDO sur `phi(theta)`
integree en RK4. Le profil se reconstruit par `X(beta) = H.u(beta) + s.t(beta)`.
La valeur du couple constant est determinee par recherche de racine pour que le
profil atteigne exactement l'encombrement maxi admissible.

Verifications automatiques : convexite (`H + H'' > 0`), monotonie du rayon,
respect de `r_min` / `r_max`, absence d'arc-boutement (changement de signe de
`q + sigma.mu.e`), course du marteau.

## Coeur de chronographe

Le **marteau est moteur** : son couple `Cm = C0 + k.|psi|` actionne le coeur,
qui restitue un couple constant sur l'axe du mobile. Depuis n'importe quelle
position le coeur part d'un cote ou de l'autre : il a deux flancs.

### Dissymetrie et courses

D'un flanc a l'autre le coeur tourne en sens inverse, donc le glissement au
contact s'inverse : frottement **rentrant** d'un cote, **sortant** de l'autre.
Les deux courses ne sont pas des entrees, elles sont **calculees** par deux
conditions :

1. la pointe est commune aux deux flancs et sur le diametre maxi ;
2. les deux sens de remise a zero donnent le meme couple sur le coeur.

La course jusqu'a la pointe ne vaut donc pas 180 deg et depend du coefficient
de frottement. Avec r_min = 0.3, r_max = 3, mu = 0.12 : **171.7 deg d'un cote,
150.1 deg de l'autre**, la pointe occupant les 38.2 deg restants.

### Niveau des couples

Le profil ne depend que des **rapports** T_coeur/C0 et k/C0 : le niveau absolu
se regle a posteriori. On saisit donc le couple vise sur le coeur et le
logiciel en deduit le couple marteau necessaire (0.3 N.mm sur le coeur demande
C0 = 5.04 N.mm avec la geometrie par defaut).

### Le plat du zero n'est pas reglable

Au zero le plan du marteau est tangent au coeur, et les deux flancs demarrent
tous deux **sur cette tangente**, aux abscisses -s_A(0) et +s_B(0). Le segment
qui les separe EST le plat du zero :

    L_plat = |s_A(0)| + s_B(0)       avec  s = dH/dbeta

Comme s vaut en moyenne (r_max - r_min)/course, l'ordre de grandeur est

    L_plat ~ 2 (r_max - r_min) / course[rad]

soit 1.78 mm pour r_min = 0.3 et r_max = 3 mm. Un plat de 0.6 mm exigerait
r_max ~ 1.1 mm. Le logiciel affiche la longueur imposee, la compare a celle
souhaitee et donne le r_max correspondant.

Raccourcir le plat a l'usinage ne raccourcit pas le plat **vu par le marteau** :
le plan enjambe le creux et porte sur les deux aretes. Le contact saute alors
d'une arete a l'autre, et c'est la que la resultante peut passer par l'axe.

### Non-blocage

Le bras de levier de la resultante par rapport a l'axe du coeur vaut

    d = |s + sigma.mu.H|

Il s'annule quand `s = -sigma.mu.H` : la resultante passe par l'axe, le couple
de rappel disparait, l'aiguille reste coincee a un petit angle du zero. Au bord
du plat `s = +/-L` et `H = r_min`, d'ou la condition

    **L > mu . r_min**

et le bras minimal au bord du plat vaut `L - mu.r_min`. Avec les cotes par
defaut : 0.80 mm pour un minimum de 0.036 mm, soit 22 fois la marge. Le bras
est trace sur 360 deg avec la limite `mu.r_min`.

### Raccords plat / flancs

La tangente est commune au plat et au flanc a leur jonction : le raccord est
**naturel**, aucun rayon n'est necessaire pour la continuite. Un rayon
d'usinage reste possible tant qu'il reste inferieur au rayon de courbure du
flanc a son depart (0.16 et 0.35 mm ici) ; au-dela il mordrait dans le profil,
ce que le logiciel signale.

### Pointe et plage morte

Au sommet, la normale de contact balaye la direction radiale : le bras de
levier `s` passe par zero, le couple s'annule, et le frottement elargit ce
point en une **plage morte** ou le coeur ne part d'aucun cote.

Cette plage est inevitable avec un marteau a face plane sur une came convexe.
Le blocage a lieu pour `|s| < mu.H`, et la vitesse a laquelle `s` traverse zero
vaut `ds/dbeta = rho_c - H`, bornee par `-H` puisque le rayon de courbure
`rho_c` est positif. La largeur minimale est donc

    2 . atan(mu)

atteinte pour `rho_c = 0`, c'est-a-dire une **pointe vive**. Mesures a
r_min 0.3 / r_max 3 / mu 0.12 (plancher 13.69 deg) :

| corde de pointe [mm] | 0.05 | 0.40 | 0.72 | 1.20 | 1.80 | 2.40 |
|---|---|---|---|---|---|---|
| plage morte [deg] | 14.10 | 20.56 | 27.14 | 35.60 | 47.12 | 62.76 |

Tout creux ou arrondi au sommet **elargit** la plage. La pente du plan du
marteau n'a pratiquement pas d'effet (35.6 a 36.0 deg entre 55 et 82 deg de
pente). Le seul levier est le coefficient de frottement. En pratique c'est
l'energie du marteau qui franchit cette plage.

Le parametre `corde_pointe` permet quand meme de tronquer le sommet par un
creux (pointe moins fragile) : le logiciel chiffre le cout en plage morte.


## Influence de la pente du plan du marteau

Etude parametrique a couple coeur impose (0.3 N.mm), r_min 0.3, r_max 3 :

| pente [deg] | C0 requis [N.mm] | N maxi [N] | course marteau [deg] | plat [mm] | courbure flanc [mm] | courses |
|---|---|---|---|---|---|---|
| 46 | 3.65 | 0.733 | 18.4 | 2.19 | 0.051 | 191/165 |
| 50 | 4.56 | 0.629 | 13.5 | 1.88 | 0.143 | 194/166 |
| 55 | 5.03 | 0.581 | 12.3 | 1.78 | 0.155 | 188/172 |
| 60 | 5.43 | 0.548 | 11.4 | 1.72 | 0.168 | 186/174 |
| 65 | 5.77 | 0.524 | 10.7 | 1.67 | 0.178 | 184/176 |
| 70 | 6.06 | 0.506 | 10.2 | 1.64 | 0.357 | 183/177 |
| 80 | 6.48 | 0.480 | 9.6 | 1.60 | 0.365 | 180/180 |

- la pente fixe le bras de levier au pivot, `q = |OP|.sin(pente)` (10.8 mm a
  40 deg, 16.7 mm a 85 deg) : le **couple marteau requis croit avec la
  pente**, c'est le prix a payer ;
- en echange l'effort de contact, la course du marteau et la longueur du plat
  diminuent, et les deux courses se rapprochent de 180 deg ;
- **en dessous de 48 deg le calcul change de branche** : le rayon de courbure
  du flanc a son depart tombe a 0.05 mm (sous le rayon de raccord demande,
  donc non usinable) et les courses sortent de 180 deg +/- 20 %. Borne
  inferieure dure ;
- au-dessus de 70 deg tout plafonne (moins de 5 % de gain jusqu'a 82 deg)
  alors que C0 continue de monter.

**Pente conseillee : 65 a 70 deg.** Prendre 50 a 55 deg si le couple
disponible au marteau est le facteur limitant. La plage morte au sommet, elle,
ne depend pas de la pente (13.6 a 14.3 deg sur toute la plage).

La pente de couple `k` a un effet secondaire : de k/C0 = 0 a 2, le C0 requis
baisse de 14 % mais le couple marteau maxi monte de 19 % et l'effort de 17 %.
La garder faible.

## Acceleration du calcul

Deux leviers, dans cet ordre.

### 1. Compilation du noyau (numba)

[came/noyau.py](came/noyau.py) compile les boucles numeriques quand numpy et
numba sont presents, et retombe sinon en Python pur. Les deux chemins sont
**valides l'un contre l'autre** : ecart nul sur la reparation du profil,
1e-9 sur la plage morte, le couple, l'effort et la courbure.

Le piege : compiler l'evaluation mecanique ne servait presque a rien, parce
qu'elle ne pesait que **1 %** du temps d'une optimisation. Le gros du cout
etait la **reparation du profil** (projections alternees pour rho >= r_outil
et fermeture du contour), appelee pour chaque individu a chaque generation et
restee en Python. Mesure par individu :

| | Python pur | numba | apres tables trigo |
|---|---|---|---|
| reparation du profil | 150.5 ms | 14.3 ms | **0.97 ms** |
| evaluation mecanique | 1.2 ms (deja) | 1.25 ms | 1.25 ms |
| **total** | **151.7 ms** | 15.5 ms | **2.22 ms** (68x) |

Les tables de cosinus et sinus precalculees une fois par reparation font le
dernier facteur 15 : sans elles les appels trigonometriques dominaient.

### 2. Repartition sur les coeurs (multiprocessing)

Le noyau compile est **mono-thread** : un seul coeur travaille, ce qui sur
une machine a 16 coeurs plafonne a **7 % d'occupation processeur**. Les
individus d'une generation etant independants, c'est la qu'il faut
paralleliser. `pool_lot()` evalue toute une generation sur un pool de
processus (stdlib, aucune dependance).

Optimisation de reference (48 individus, 20 generations, face libre) :

| | duree |
|---|---|
| Python pur, 1 coeur | ~150 s |
| numba, 1 coeur | 95 s puis 15 s apres correction de la reparation |
| numba + 15 processus | **2 s** |

Le nombre de processus se regle dans l'interface (par defaut : coeurs - 1).

### Environnement

Le Python du systeme n'a ni pip ni numpy : il fait tourner le programme en
Python pur, sans acceleration. L'environnement virtuel du projet est

    /home/florian/DEV/virtual_env/.came

Il contient PyQt6, numpy et numba. **Il faut les trois** : sans PyQt6 dans le
venv il faudrait choisir entre l'interface et l'acceleration.

    /home/florian/DEV/virtual_env/.came/bin/python main.py

Pour le recreer ailleurs :

    python3 -m venv .came
    .came/bin/pip install -r requirements.txt

Sans cela tout fonctionne a l'identique en Python pur, simplement plus
lentement ; l'interface indique en vert ou en orange si l'acceleration est
active.

## Sens du calcul : le ressort est une sortie

C'est le sens naturel du probleme. On connait :

- l'**encombrement** du coeur (r_min, r_max) ;
- le **pivot du marteau** ;
- la **pente du plan** du marteau ;
- le **couple de rappel voulu** sur le coeur ;
- le frottement.

On en deduit la **geometrie du coeur**, puis le **ressort du marteau** :
prearmage et raideur.

### Comment le ressort se calcule

Le rapport de transmission `T / Cm` ne depend que de la geometrie et du
frottement, pas du niveau de couple. On le mesure donc une fois avec un
couple marteau unitaire, puis on en deduit le couple marteau **exige** en
chaque position :

    Cm_exige(theta) = couple_vise / rapport(theta)

Un ressort reel est affine en l'ouverture du marteau. On ajuste donc
`Cm = C0 + k.psi` au sens des moindres carres sur la partie utile du
parcours : **C0 est le prearmage, k la raideur**. Le resultat est ensuite mis
a l'echelle pour que le couple **mini** atteigne la consigne (l'ajustement
passe au milieu du nuage, il laisserait le minimum en dessous), et le residu
de l'ajustement est l'ondulation qu'un ressort lineaire ne peut pas eviter.

Toutes les solutions delivrant alors le meme couple, on les compare sur ce
qu'il en coute : **ondulation, plage morte, effort de contact, couple marteau
maxi** (qui dimensionne le ressort).

### Exemple complet

Entrees : coeur 0.30 a 3.00 mm, pivot (16, 5), pente 65 deg, couple vise
0.300 N.mm, mu 0.12, galet au bout du marteau.

    Ressort de marteau (sortie)
      prearmage C0                 4.856 N.mm
      raideur k                   17.304 N.mm/rad   (0.302 N.mm/deg)
      couple a pleine ouverture    7.859 N.mm
      course du marteau            9.95 deg
      residu de l'ajustement       11.1 %

    Couple obtenu
      couple mini garanti          0.300 N.mm sur 77 % du tour
      ondulation                   42 %
      effort de contact maxi       0.530 N
      plage morte                  1.0 deg

    Profil
      galet R = 6.26 mm, mu au contact 0.0037
      rayon de courbure 0.089 a 11.36 mm (raccord usinable 0.08 impose)
      rayon du profil 0.407 a 2.998 mm

### Ou lire le ressort dans l'interface

## Le plat du zero disparaissait du profil optimise

Signale : l'optimisation fait n'importe quoi, le plat du coeur n'existe plus,
donc pas de position nulle.

### Le modele ne peut pas porter un plat exact

Mesure faite : la **graine** elle-meme n'a qu'un plat de 0.31 mm la ou la
synthese en impose 1.72, et cela ne vient pas de l'optimisation. Un plat
demande un Dirac de courbure, et la convergence d'une serie a bande limitee y
est trop lente pour qu'un nombre d'harmoniques y suffise :

| harmoniques | plat rendu (tol 2 um) |
|---|---|
| 8 | 0.27 mm |
| 12 | 0.31 mm |
| 24 | 0.42 mm |
| 40 | 0.52 mm |
| 96 | 0.67 mm |

J'ai essaye de le porter explicitement : `H = H_plat + H_bande` avec
`H_plat = L.sin(x/2)`, dont la derivee saute bien de L - le plat mesure est
alors exact. Mais la decomposition **n'est pas valide** : elle exigerait que
le coeur contienne la forme du plat comme sommande de Minkowski, ce que rien
ne garantit. Le profil sortait a 0.24..1.41 mm au lieu de 0.30..3.00, et sa
courbure totale devenait negative. L'essai est annule, et la raison est
consignee dans `Profil`.

Le porter correctement demande de **reparametrer l'espace d'optimisation** :
`rho = Dirac + partie lissee`, avec le rang 1 de la partie lissee libre pour
compenser celui du Dirac - alors que c'est aujourd'hui l'absence de rang 1 qui
assure la fermeture du contour. C'est un chantier a part.

### Ce qui est fait

`profil.longueur_plat()` mesure le plat sur la geometrie : l'etendue, le long
de la tangente au repos, des points qui s'en ecartent de moins de 2 um. Le
rapport l'affiche, avec un avertissement quand la serie ne rend pas ce que la
synthese visait. Et `nsga2.probleme` recoit un `plat_mini`, fixe par
l'appelant a 90 % du plat de la graine : l'optimisation **ne peut plus le
degrader**. Verifie sur douze perturbations : les solutions recevables gardent
0.296 a 0.306 mm, contre 0.309 pour la graine.

## Enregistrer une solution de l'optimisation levait une TypeError

`projet.sauver` lisait les points dans `profil.pts`. Un profil importe en a,
un profil **optimise** est une suite de coefficients et n'en a pas : `pts`
valait alors l'objet `Profil` lui-meme, d'ou
`TypeError: 'Profil' object is not iterable`. `_points()` echantillonne
desormais le contour quand l'objet n'a pas de points, et accepte aussi une
simple liste.

## Le marteau se placait sur la pointe, pas sur le plat

Sur l'onglet du chronographe, le curseur de position part a zero - et zero y
est le **debut de la course**, soit `-course_b`, c'est-a-dire la pointe. Le
marteau apparaissait donc toujours sur la pointe a l'ouverture de l'onglet. Le
curseur est maintenant amene a la position de repos, `x = 0`, ou le plat du
coeur est contre le plan du marteau. Verifie : curseur a 48 %, x = -1.2 deg,
rayon de contact 0.998 mm pour un coeur de 0.3 a 3.0 - donc sur le plat.

## Le couple sur un arrondi de coin

Question posee : qu'est-ce qui est calcule au niveau du couple quand le coeur
travaille sur un des rayons des coins du plan.

Rien de particulier : **la meme formule que sur le plan**. Ce qui change, ce
sont les deux grandeurs geometriques qui y entrent.

La normale de contact est celle du **coeur** au point touche, `beta_j`, et non
celle du plan. Le point de contact est le pied de la perpendiculaire issue du
centre de l'arrondi : si `A` est ce centre et `rho` son rayon,

    H = u_j.A - rho          s = t_j.A

ce qui est exact, le point de contact et le centre partageant la normale. `s`
est le bras de levier autour de l'axe du coeur - c'est lui qui ne s'annule
plus a la pointe, et c'est tout le benefice de l'arrondi.

Ensuite, a l'identique du contact sur le plan :

    q = s - t(gamma).P       e = H - u(gamma).P
    T = Cm (s + sigma.mu.H) / (q + sigma.mu.e)
    N = |Cm / (q + sigma.mu.e)|

avec `sigma = +/-1` le sens du glissement. Les deux sens sont essayes et l'on
retient celui qui rend le mouvement possible ; si aucun ne le rend possible,
la position est **arc-boutee** et compte dans la plage morte. C'est le meme
`_couple()` pour les deux types de contact.

## A quoi servent les trois onglets

La question a ete posee, et l'interface ne la repondait pas. Chaque onglet
porte maintenant une ligne qui dit son role, et une infobulle qui la detaille.

| onglet | ce qu'il fait | plan du marteau |
|---|---|---|
| Coeur simple | **construit** un coeur a un flanc a partir du couple voulu | suppose assez long |
| Coeur de chronographe | **construit** le coeur a deux flancs : plat du zero, raccords, pointe | suppose assez long |
| Optimisation du profil | **simule et optimise** un coeur qui existe deja | **modelise tel qu'il est** |

Les deux premiers font une synthese inverse : on impose le couple et on en
deduit la forme, en integrant le long du contact. Ils cherchent la tangence
**sans borne d'abscisse**, donc ils supposent le plan du marteau assez long.
Le troisieme ne construit rien : il prend le coeur qui existe - importe,
relu d'une etude, ou issu d'une des deux syntheses - et il est le **seul** a
modeliser le plan reel, longueur finie, arrondis de bout et decalage. C'est
donc la que les chiffres decrivent la piece.

### Les parametres ne reagissaient pas tous de la meme facon

Sur l'onglet d'optimisation, changer l'angle du plan ou le pivot relancait le
calcul - c'est le recalage plat sur plan - mais changer la **longueur du
plat**, le rayon des coins, le decalage, le frottement ou le couple vise ne
mettait rien a jour. L'interface repondait a certains champs et pas a
d'autres, sans que rien ne le laisse voir.

Ces cinq-la ne changent pas la forme du coeur, seulement sa simulation :
`_resimuler_profil()` la rejoue desormais sans toucher au profil. Verifie sur
le coeur de minuterie, chaque champ deplace bien le resultat :

| champ modifie | plage morte | effort maxi |
|---|---|---|
| (depart) | 1.247 deg | 1.012 N |
| longueur du plan 2.0 mm | 1.101 deg | 0.956 N |
| rayon des coins 0.25 mm | 1.919 deg | 0.940 N |
| decalage 0.4 mm | 2.150 deg | 0.983 N |
| frottement 0.20 | 3.469 deg | 1.576 N |

### « Le coeur ne se construit plus » apres un import

C'est vrai, et c'est geometrique. Un pivot et une pente lus dans un dessin ne
viennent pas d'une synthese : rien ne garantit qu'un coeur existe pour eux.
Avec un plan plat, le marteau ne peut ecarter son plan de l'axe que de
`|OP|(1 - cos a) + r_min` ; si cela ne couvre pas l'encombrement demande,
aucun profil n'existe, quelle que soit sa forme. Sur `test.dxf` : 3.055 mm
atteignables pour 3.500 demandes.

L'onglet du coeur simple le disait deja. Celui du chronographe repondait
« pas de solution : verifier pivot, pente et rayons », ce qui n'aide pas :
`_pourquoi_pas()` donne maintenant le meme diagnostic chiffre, avec les trois
leviers - augmenter la pente, eloigner le pivot, ou reduire l'encombrement.
L'onglet d'optimisation, lui, continue de simuler sans probleme : il n'a pas
besoin qu'un coeur soit **synthetisable**, seulement qu'il existe.

## L'optimisation rendait des solutions pires que sa graine

Signale : les solutions sont nettement moins bonnes que le coeur construit par
l'onglet du chronographe, quelle que soit la combinaison d'objectifs. C'etait
exact, et pour deux raisons qui se cumulaient.

### La reparation brouillait le vecteur

`profil.reparer` attend le vecteur **complet**, `cx` et `cy` compris : il les
lit en queue et en deduit le nombre d'harmoniques. `nsga2.reparer_x` lui
passait `x[:-nq]`, donc deja prive de ces deux composantes, et il en retirait
deux autres : le nombre d'harmoniques tombait de 11 a 10 et les series `an` et
`bn` se trouvaient **decalees d'un rang**.

Effet mesure sur la graine analytique : `viol` passait de 0 a 0.032 et les
quatre objectifs de (62.8 ; 14.76 ; 103.3 ; 0.674) a
(64.4 ; 15.64 ; 107.2 ; 0.689) - exactement la solution que l'optimisation
finissait par rendre. La reparation est maintenant neutre sur un profil deja
admissible, verifie au bit pres.

### Un depassement de 0.9 um disqualifiait la graine

La contrainte d'encombrement etait stricte. La graine, reconstruite a bande
limitee, sort a `r_max = 3.0009` mm pour 3.0000 demandes : 0.9 um, residu de
reconstruction et non defaut de conception. En domination contrainte, **toute**
solution recevable battait alors la meilleure, aussi mauvaise soit-elle.
`TOL_ENCOMBREMENT = 0.005` mm donne la tolerance qui manquait.

### Resultat

La graine est a nouveau recevable et figure dans le front, qui ne peut donc
plus etre globalement plus mauvais qu'elle. Front obtenu en 8 generations,
plan suppose infini :

| | graine | front |
|---|---|---|
| tour sous la consigne | 62.8 % | 57.2 a 64.4 % |
| plage morte | 14.76 deg | 14.31 a 15.52 deg |
| ondulation | 103.3 % | 102.6 a 146.3 % |
| effort de contact | 0.674 N | 0.672 a 0.796 N |

### Ce que l'optimisation cherche, et ou est le levier

Le rapport l'affiche maintenant en clair, avec le gain objectif par objectif
par rapport a la graine. Le point a retenir : **le profil du coeur n'a presque
pas de prise sur la plage morte**. Avec un plan de marteau suppose infini elle
reste au plancher d'arc-boutement, et les quatre objectifs se disputent des
dixiemes. En renseignant une longueur de plan de 1.4 mm, le meme optimiseur
descend a **1.52 deg** de plage morte et 66 % d'ondulation, contre 14.76 et
103 %. Le levier est le plan du marteau, pas le coeur - l'avertissement est
desormais ecrit dans le rapport quand la longueur de plan est laissee a zero.

## Le profil approche par des arcs de cercle tangents

Un dessin d'atelier ne se cote pas sur une polyligne de 720 points mais sur
une suite d'arcs, chacun par son centre et son rayon. `came/arcs.py` les
calcule, avec l'**ecart maxi** comme donnee : le nombre d'arcs en decoule.
C'est le bon sens de la demande - un nombre d'arcs impose ne dirait rien de la
precision obtenue.

### La tangence est automatique

Le calcul se fait dans la fonction d'appui. Un cercle de centre C et de rayon
R y vaut `H(beta) = C.u(beta) + R`. Deux arcs voisins qui se raccordent a la
normale `beta_k` partagent le point de contact si

    C_k = C_{k-1} + (R_{k-1} - R_k) u(beta_k)

et ils y partagent alors la **meme normale**, donc la meme tangente : le
raccordement est G1 par construction, sans biarc ni point de jonction libre.
Toute la chaine se deduit de `C_0` et des rayons, donc `H_arcs(beta)` est
**lineaire** en ces inconnues : l'ajustement est un simple moindres carres.

L'ecart se lit directement : `H(beta) - H_arcs(beta)` est la distance entre
les deux tangentes de meme normale, soit l'ecart mesure le long de la normale.
C'est lui que l'on borne.

La fermeture du contour donne deux equations lineaires sur les rayons,
`somme_k (R_{k-1} - R_k) u(beta_k) = 0`, posees en contrainte de l'ajustement
et non corrigees apres coup.

### Resultats sur le coeur de minuterie

| ecart vise | arcs | ecart obtenu | residu de fermeture |
|---|---|---|---|
| 0.05 mm | 8 | 0.043 mm | 5e-09 mm |
| 0.02 mm | 14 | 0.018 mm | 1e-09 mm |
| 0.01 mm | 18 | 0.0086 mm | 7e-10 mm |
| 0.005 mm | 23 | 0.0047 mm | 4e-10 mm |
| 0.002 mm | 56 | 0.0020 mm | 1e-11 mm |

Discontinuite entre arcs voisins : 7e-10 mm. Quatre des dix-huit arcs ont un
rayon **negatif** : ce sont les raccords concaves entre les flancs et le plat
du zero. `Arc.pour_dxf()` les rend avec un rayon positif et les angles
decales de 180 deg, ce qu'attend une entite ARC.

### Dans l'interface

Un groupe *Profil en arcs de cercle* porte l'ecart maxi et une case
*montrer le profil en arcs*. La chaine est superposee au contour en pointille
bleu, avec un point a chaque raccord, pour voir ou l'approximation derive. Le
rapport donne le nombre d'arcs, l'ecart obtenu, le residu de fermeture et la
plage des rayons. Les arcs sont ecrits dans le DXF sur le calque `COEUR_ARCS`,
en vraies entites ARC.

## L'export ecrivait douze marteaux

Un seul marteau est desormais ecrit, **a la position de repos**, plat sur
plat. Une case *exporter les 12 positions du marteau* ajoute les onze autres
pour qui veut verifier la cinematique relative.

## L'etude enregistree ne contenait pas toujours le bon coeur

Signale : une etude sauvee apres un changement d'angle du plan montrait une
plage morte de 12 deg et perdait la position zero, mais la meme etude relue
simulait correctement.

`_on_sauver` prenait le profil de l'onglet **actif**. Apres un import, revenir
sur l'onglet d'une synthese y affiche le coeur **synthetise** - dont la plage
morte a la pointe est le plancher d'arc-boutement, soit les 12 deg observes -
et c'est lui qui etait enregistre, a la place des 720 points importes. L'etude
relue repassait alors par le chemin du profil importe, d'ou une simulation
correcte mais sur une autre geometrie.

Le sujet de l'etude est maintenant le contour importe ou optimise, quel que
soit l'onglet affiche, et le commentaire du fichier dit lequel des deux a ete
ecrit. Verifie : 720 points enregistres depuis chacun des trois onglets.

## L'export DXF est ecrit par ezdxf

Signale : le fichier produit donnait, dans un logiciel de CAO, des points, de
petits cercles et des droites, mais rien qui ressemble a un coeur ou a un
marteau.

L'ecriture etait faite a la main, avec une section HEADER **vide** et
**aucune table de calques**. Beaucoup de logiciels ne relisent alors que
partiellement les entites. `dxfio.Dessin` accumule maintenant les entites et
les ecrit par **ezdxf** en R2010 complet - en-tete, tables, poignees - avec
les calques declares et colores, `$INSUNITS` en millimetres, et des contours
qui sont de vraies polylignes fermees. Les trois ecrivains (coeur simple,
chronographe, profil optimise ou importe) passent par la.

Verification : `ezdxf.audit()` ne releve **aucune erreur** sur les fichiers
produits, le coeur sort en une polyligne fermee de 720 sommets, et les
positions de cinematique en polylignes de 22 points sur leur propre calque.

Le repli en texte brut reste en place si ezdxf manque - il ecrit desormais une
table de calques - mais `ezdxf>=1.1` est ajoute a `requirements.txt` : c'est
la voie normale.

## Le recalage plat sur plan ne suivait pas les changements d'angle

Le recalage a l'import vaut pour l'angle de plan et le pivot de ce moment-la.
Or `gamma0 = angle(O->P) + angle_plan` : changer l'angle du marteau ou son
pivot deplace la normale du plan, et le coeur n'etait plus a sa position de
repos - le marteau recommencait sa course deja ouvert.

`dxfio.aligner_contour()` est maintenant une fonction publique, idempotente,
que l'interface rappelle a chaque changement de `angle_plan`, `pivot_x` ou
`pivot_y`. Mesure sur `test.dxf`, en balayant l'angle du plan :

| angle du plan | phi(0) | ecart plat / plan |
|---|---|---|
| 43.6 deg | -0.008 deg | 0.03 deg |
| 50 deg | -0.017 deg | -0.06 deg |
| 60 deg | -0.014 deg | 0.19 deg |
| 30 deg | -0.030 deg | -0.06 deg |
| -40 deg | +0.012 deg | -0.06 deg |

Le residu est la resolution de la grille de la fonction d'appui (0.25 deg).

Le **decalage du plan** ne change pas `gamma0` : il fait glisser le plan le
long de lui-meme, le plat du coeur reste parallele, et `phi(0)` ne bouge pas -
ce que la mesure confirme. Il ne deplace le contact au repos hors du plat que
si le plan est trop court pour le couvrir, ce qui est un fait de geometrie et
non un defaut de recalage.

## La selection du contact etait fausse : le marteau penetrait le coeur

Signale sur `test.dxf` apres le recalage de la position de repos : perte de
contact vers 239 deg. Mesure contre une reference par force brute, le solveur
**sous-ouvrait de 26 deg** a partir de theta = 209 : ouverture exacte 31.06
deg, ouverture retenue 4.63. Le marteau penetrait le coeur de 0.34 mm, et phi
sautait de 26 deg d'une position a la suivante - impossible pour un solide
rigide.

### L'erreur etait conceptuelle

Je calculais, pour chacune des trois caracteristiques du nez - le plan et ses
deux coins - la **premiere** ouverture qui la degage, puis je prenais la plus
grande des trois. Or l'ensemble des ouvertures admissibles d'un coin est un
**arc**, pas un intervalle croissant : a `psi = 0` le coin cote -t etait
degage, a `psi = 4.6` il etait rentre dans le coeur, et la regle du maximum
retenait 4.6. Un maximum de premieres-solutions n'est pas une solution
commune.

### Le bon critere est unique et exact

Le nez du marteau est le segment des deux centres de coin, dilate du rayon
`rayon_coin`. Le coeur, vu du contact, est son enveloppe convexe, donc une
intersection de demi-plans `u_j.x <= H_j`. Un segment est hors d'une telle
intersection **si et seulement si** un appui le separe en entier :

    min(u_j.A, u_j.B) >= H_j + rho

L'ensemble admissible est exactement l'union sur j de ces conditions, et
l'ouverture cherchee en est le minimum - sans hypothese de monotonie, et sans
avoir a combiner trois candidats. `profil._ouverture_nez()` et
`noyau.scan_plan_fini` l'appliquent, en forme close appui par appui (chaque
appui donne l'intersection de deux arcs en psi).

Accord avec la force brute sur la zone fautive : **0.002 deg** au lieu de
26 deg. Plus aucun saut de phi superieur a 0.29 deg, les deux moteurs
d'accord a 3.6e-15.

### Plan ou coin : par la projection, pas par la tangence

Premiere version du correctif : le contact etait toujours classe « coin » et
le bras de levier pris au centre du coin, meme quand c'est le plan qui porte -
jusqu'a 2 mm d'erreur sur le bras, et sept changements de signe du couple sur
un tour. Le classement se fait maintenant sur la **projection** du point de
contact : si le point d'appui du coeur se projette a l'interieur du segment,
c'est le plan qui porte et le bras de levier est celui de ce point ; sinon
c'est un coin, et c'est le pied de la perpendiculaire issue de son centre. Sur
le coeur de minuterie : 181 contacts sur coin et 179 sur le plan, deux
changements de signe du couple.

## Le coeur importe ne partait pas de sa position de repos

Signale sur `test.dxf` : le couple transmis tombait a zero sur une quinzaine
de degres **en plein profil actif**, loin de la pointe. C'etait bien la pointe,
decalee.

Le modele suppose que `theta = 0` **est** la position stable : le marteau y est
a `phi = 0` et touche le plat du zero. L'orientation du coeur dans un DXF est
en revanche quelconque. Sur `test.dxf`, la normale du plat etait a 129.50 deg
et celle du plan du marteau a 146.63 deg : **17.13 deg d'ecart**. Le marteau
demarrait donc deja ouvert (`phi(0) = -2.09 deg`, contact a H = 1.004 mm au
lieu de 0.701), et tout le parcours se trouvait decale d'autant - la plage
morte de la pointe apparaissait en plein profil.

`dxfio._aligner_sur_le_plat()` tourne desormais le contour pour amener la
normale de son plat sur celle du plan du marteau. La normale du plat est la
direction ou la fonction d'appui est minimale, c'est-a-dire exactement la
reference que le moteur prend au repos. Apres recalage : `phi(0) = -0.005 deg`,
et la plage de couple nul revient a theta = 179-192 deg, au rayon 3.500 mm,
soit r_max : la pointe.

### Le plan lu etait le bras du marteau

`test.dxf` passe par la lecture de secours, qui prenait la plus longue droite
du calque. Sur ce fichier c'est le bras pivot -> contact : **6.44 mm** au lieu
de 4.00 mm, une pente de 77.18 deg au lieu de 43.59, et aucun rayon de coin.

`_plan_depuis_polyligne()` lit maintenant le contour **materiel** du marteau
quand il est dessine comme une polyligne - ce que fait l'export : son plus
long segment est la partie droite, et les trois points de chaque bout donnent
le rayon de l'arrondi par cercle circonscrit. Avec la lecture correcte, la
plage de couple nul **disparait** : plage morte 0.12 deg, aucune position a
couple nul.

Les quatre fichiers d'essai s'importent desormais de facon coherente, tous a
`phi(0)` nul a 0.005 deg pres.

## Fichier d'etude .json

Reprendre une etude demandait de re-saisir une vingtaine de champs puis de
re-importer un DXF. `came/projet.py` ecrit et relit tout d'un coup :
parametres, contour du coeur, ressort calcule et reglages de l'optimisation.
Deux boutons dans le bloc *Fichiers* du panneau.

Le format est plat et lisible - une ligne par point, pour qu'un contour de
720 points tienne en 764 lignes au lieu de 2880 - avec les angles en degres et
les longueurs en millimetres comme dans l'interface. Un champ inconnu de la
version qui relit est signale au journal plutot qu'ignore en silence.
Aller-retour verifie exact sur les 24 parametres et les 720 points.

## Le plan reel du marteau n'entrait pas dans le moteur de contact

Les deux syntheses analytiques cherchent la tangence **sans borne
d'abscisse** : elles supposent le plan infini et n'ont jamais lu
`longueur_plan`, `rayon_coin` ni `decalage_plan`. Des que l'utilisateur les
renseignait, la cinematique affichee ne decrivait plus la piece et le marteau
quittait le coeur :

| | ecart marteau / coeur | positions decollees |
|---|---|---|
| Coeur simple, plan 4.0 mm | 0.367 mm | 22 / 61 |
| Coeur simple, plan 1.4 mm | **0.973 mm** | **50 / 61** |
| Chronographe, plan 1.4 mm | 0.650 mm | 20 / 35 |

`ui.simuler_plan_fini()` rejoue donc le contour obtenu dans le moteur de
contact a plan fini - celui qui gere la tangence bornee et les deux arrondis
de bout, verifie contre une reference par force brute. C'est cette simulation
qui sert au trace, au point de contact, aux courbes, au rapport et a l'export
DXF des deux onglets analytiques. L'ecart tombe a **0.006 mm**.

### Correction : ma regle sur la longueur du plan etait fausse

J'avais annonce `longueur_plan_maxi()` comme une regle de construction : la
plus grande longueur de plan telle que l'arrondi porte a la pointe, cherchee
par dichotomie. Deux defauts.

**La fenetre de pointe n'est pas monotone en longueur de plan** - 80 % a
3.0 mm, 40 % a 1.6, 45 % a 1.4 sur le coeur synthetise. Une dichotomie sur un
predicat non monotone n'a pas de sens.

**Et elle ignorait le secteur sans contact.** Elle conseillait 1.054 mm sur le
coeur synthetise, ou le marteau lache le coeur sur **62 deg**, alors que
3.0 mm donne la meme plage morte avec 2 deg sans contact. Le critere etait un
proxy de la plage morte, qu'on sait maintenant mesurer directement.

`balayage_plan()` balaye donc les longueurs et rend le tableau complet
(fenetre de pointe, plage morte, secteur sans contact, effort), et
`plan_recommande()` prend la plus petite plage morte parmi les longueurs dont
le secteur sans contact reste sous 5 deg :

| coeur | conseille | plage morte | sans contact | N max |
|---|---|---|---|---|
| chronographe synthetise | 2.655 mm | 1.37 deg | 2 deg | 0.337 N |
| importe (minute) | 1.564 mm | 1.00 deg | 0 deg | 0.462 N |

La piece reelle de l'utilisateur fait 1.119 mm, pour 1.14 deg de plage morte
et 0 deg sans contact : elle est dans la plage utile, et sur ce coeur le
compromis n'existe meme pas - le secteur sans contact reste nul de 0.7 a
3.0 mm.

Quand **aucune** longueur ne tient la contrainte, il n'y a pas de bon choix a
annoncer : sur le coeur simple par defaut la marche fait 4 mm de profondeur et
le minimum atteignable est 24 deg sans contact. Le rapport rend alors le front
non domine sur (plage morte, secteur sans contact), marque d'une fleche, et
dit que l'arbitrage revient au concepteur - avec l'avertissement que la plage
morte perd son sens quand une grande part du tour est sans contact.

La contrainte de fenetre de pointe a aussi ete **retiree de NSGA-II** : la
plage morte y est deja le second objectif, le proxy etait non monotone, et le
profil optimise est de toute facon convexe (rho >= rayon d'outil), donc sans
marche a enjamber.

### L'onglet chronographe ignorait les changements de plan

Toute modification de parametre y mettait `cres = None` en attendant un
nouveau clic sur *Calculer le coeur*, mais la **vue** gardait sa simulation
precedente : changer la longueur du plan ne changeait rien a l'affichage. Or
`longueur_plan`, `rayon_coin` et `decalage_plan` ne decrivent que le marteau
et ne changent pas la geometrie du coeur : la synthese reste valable, seule la
re-simulation est a refaire. C'est ce que fait desormais `PLAN_MARTEAU`.
Ecart mesure apres correction : 0.00033 mm, et le secteur sans contact suit la
longueur du plan (1 deg a 3.0 mm, 29 deg a 1.4 mm, 78 deg a 0.8 mm) - un plan
trop court finit par lacher le coeur, ce que le rapport annonce.

### Le reste est une concavite, pas un bug

Un plan, et meme un coin arrondi, ne peut toucher que l'**enveloppe convexe**
du coeur. La marche qui ramene la pointe au plat du zero est une concavite :
sur le coeur simple par defaut elle fait 3.86 mm de creux, et sur ce secteur
le marteau franchit la corde de l'enveloppe sans rien toucher. Mesure : le
secteur va de theta = 45 a 84 deg, ou le rayon de contact tombe de 2.71 a
1.10 mm, soit exactement la marche.

Ce n'est pas un defaut de calcul, c'est une chute libre - le modele quasi
statique n'y a pas de sens. `profil.secteur_marche()` le detecte (point de
contact hors du dessin), le marteau est trace **en orange** sur ce secteur
avec la mention correspondante, et le rapport en donne la largeur en degres :
43 deg pour le coeur simple a plan de 1.4 mm, 30 deg pour le chronographe.

### Le trace suit la fonction d'appui, pas les points

`ProfilImporte.contour()` renvoyait le contour reconstruit par la fonction
d'appui. Il coupe les angles entre deux echantillons de normale, jusqu'a
14 um au sommet d'un coeur synthetise. Renvoyer les points eux-memes donnait
0.3 mm d'ecart apparent, parce que le calcul travaille sur l'enveloppe
convexe : le trace reste donc la reconstruction, qui est ce que le contact
voit, et `contour_reconstruit()` reste disponible.

## L'arrondi du marteau est le levier sur la plage morte

Constat de l'utilisateur sur sa piece reelle (`Sketches - minute.dxf`) : le
point mort au passage de la pointe n'y vaut que **1,2 deg**, parce que le
contact s'y fait sur le **rayon du bord du marteau** et non sur l'arete du
plan. C'est exact, et cela invalide le caractere incontournable du plancher
`2.atan(mu)` que j'avais annonce : ce plancher ne vaut que pour une pointe vive
en appui **plan**. Sur un arrondi, le bras de levier ne s'annule plus.

Mesure sur le coeur synthetise, en raccourcissant le seul plan du marteau :

| longueur du plan | hors arrondi dans +/-10 deg | plage morte |
|---|---|---|
| 8.0 mm (ce que la synthese exige) | 100 % | 15.85 deg |
| 4.0 mm | 65 % | 11.96 deg |
| 3.0 mm | 43 % | 3.07 deg |
| 2.0 mm | 38 % | 2.73 deg |
| 1.2 mm | **0 %** | **2.73 deg** |

Influence du rayon de coin, a 1.2 mm de plan : 2.32 deg a rc = 0.02 mm,
2.58 a 0.10, 3.51 a 0.30. Un petit rayon reduit la plage morte mais augmente
la pression de contact : c'est le compromis a regler.

### Regle de construction

`profil.longueur_plan_maxi(p, pr)` renvoie par dichotomie la plus grande
longueur de plan telle que l'arrondi porte sur toute la fenetre de +/- 10 deg
autour de la pointe. Verification croisee :

* coeur synthetise : 1.553 mm, pour une plage morte de 2.48 deg ;
* coeur importe : 1.942 mm -- et la piece reelle de l'utilisateur fait
  1.119 mm, donc largement du bon cote.

`profil.fenetre_pointe(e, demi)` donne le pourcentage de cette fenetre ou
l'arrondi ne porte **pas** : zero est le but. Le rapport l'affiche avec la
longueur maxi admise, et `nsga2.probleme` en fait une contrainte (penalite
proportionnelle a ce pourcentage). Avec un plan de 1.4 mm, la graine
analytique passe de 14.84 a **2.11 deg** de plage morte.

### Nouveau parametre : le decalage du plan

Le milieu du plan etait implicitement au pied de la normale issue de l'axe du
coeur, au repos. `decalage_plan` le deplace le long du plan, depuis ce meme
pied, quand le systeme est dans sa position stable (plat du coeur contre plat
du marteau) : c'est lui qui decide **lequel des deux arrondis** verra passer la
pointe. Il agit sur le trace, sur la fenetre de tangence du plan et sur la
position des deux centres de coin, dans les deux moteurs (accord 1.8e-15).

## Panneau de gauche : deux vides de 237 px

Le bloc d'onglets descendait de plus de 200 px au passage sur le
chronographe, et le resume du bas s'en detachait d'autant. Deux causes
superposees.

1. **Le resume replie annoncait un minimum de 280 px** pour 110 px reellement
   occupes. Avec une politique horizontale `Ignored`, Qt calcule la hauteur
   mini d'un `QLabel` replie a largeur minimale. Ce minimum superieur a la
   place prise faisait distribuer a Qt des centaines de pixels orphelins.
   `soft()` annule donc ce minimum.
2. **`heightForWidth` d'un QTabWidget porte sur toutes ses pages.** La zone
   defilante dimensionne le panneau par `heightForWidth` des qu'un texte y est
   replie : elle recevait la hauteur de l'onglet d'optimisation (985 px) sur
   l'onglet du chronographe (511 px), soit 474 px de trop, repartis en deux
   vides de 237 px. `OngletsCompacts` redefinit `hasHeightForWidth` et
   `heightForWidth` sur la **page affichee**.

Tous les ecarts entre blocs sont maintenant de 6 px, l'espacement normal, et
le bloc d'onglets ne bouge plus d'un onglet a l'autre.

## Infobulles sur les champs

Les 24 champs de parametres et les 5 reglages de l'optimisation portent une
infobulle, sur la case **et** sur son libelle : ce que la valeur represente,
son unite et sa plage admise. Les textes sont dans `ui.AIDES` et
`ui.AIDES_OPTIM`.

## Revue du pipeline de contact

Symptome signale : le marteau perd souvent le contact avec le coeur. Mesure
faite sur le coeur importe de `Sketches - minute.dxf`, distance reelle entre
l'enveloppe convexe du coeur et la matiere du marteau, 120 positions :

| etat | positions sans contact | jeu maxi | ecart au calcul exact |
|---|---|---|---|
| avant | 74 / 120 | 0.063 mm | 0.50 deg |
| apres | 0 / 120 | 0.0007 mm | 0.0042 deg |

La reference est une force brute independante : on ferme le marteau par
dichotomie jusqu'a la premiere intersection entre son contour reel et le
polygone du coeur. Quatre defauts distincts.

### Le marteau passait de l'autre cote du coeur

La branche « coin » acceptait `abs(ph) <= psi_max`, donc un `phi` de signe
oppose a l'ouverture. Entre 171 et 198 deg, `phi` valait -5.9 deg contre +6.0
juste apres : le marteau traversait le coeur et y revenait d'un bond. Le
contact y demandait une normale a 125 deg du plan, impossible sur un quart de
cercle. On impose maintenant `0 <= eps*phi <= psi_max`.

### Le contact coin / facette plane etait structurellement manque

La condition etait cherchee par racine sur la courbe offset du coeur,
`Y(beta) = (H+rho).u + s.t`. Or sur une **facette droite** tous les points
partagent une seule normale beta : la portion droite de l'offset se reduit a
un echantillon isole, et la racine passait a cote. Le marteau penetrait alors
le coeur jusqu'a 0.3 deg.

La condition est desormais exacte et en forme close. Le nez du marteau est le
voisinage de rayon `rho` du segment de contact : le coin ne penetre pas si son
centre A reste hors du coeur dilate de rho, c'est-a-dire s'il existe un appui
j tel que `u_j.A >= H_j + rho`. Comme A tourne autour du pivot,

    u_j.A(psi) = u_j.P + L.cos(beta_j - alpha_c - eps.psi)

chaque appui donne un intervalle d'angles admissibles, et l'ouverture requise
est le plus petit psi appartenant a l'un d'eux. Plus de recherche de racine,
plus de pas de balayage, plus de racine manquee. Facettes et sommets sont
traites de la meme facon.

### Le bras de levier etait pris sur un sommet quelconque

Pour un contact sur coin, `(H, s)` etait lu dans la table du profil a la
normale retenue. Sur une facette, `s` est la coordonnee tangentielle d'un
sommet arbitraire de cette facette : le rapport de transmission sautait d'une
position a l'autre et `N_max` passait de 7.96 a 42.27 N selon le pas de
calcul. Un arrondi touche une facette au **pied de la perpendiculaire** issue
de son centre, d'ou

    H = u_j.A - rho        s = t_j.A

Le rapport de transmission est redevenu lisse et monotone par flanc, et
`N_max` ne depend plus du pas (0.448 N a n = 90, 120, 180 et 360).

### La plage morte se mesurait en nombre d'echantillons

`max(runs) * 360/n` : un echantillon bloque isole pesait 3 deg a n=120 et
1 deg a n=360, ce qui faisait varier un critere d'optimisation avec la finesse
du calcul (0 / 2.9 / 0 / 1.5 deg selon n). On interpole maintenant les
passages par zero de la marge d'arc-boutement `|s| - mu.H`, evaluee avec le H
et le s **exacts** du profil sur une grille fine en theta, la normale visitee
etant obtenue en deroulant `beta(theta)`. Resultat stable et convergent :
1.34 / 1.14 / 1.12 / 1.12 deg pour n = 90 / 180 / 360 / 720.

A noter : une pointe vive en appui plan a pour plancher `2.atan(mu)`, mais un
**coin arrondi passe sous cette valeur** - c'est le mecanisme par lequel le
rayon de coin reduit la plage morte.

### Les deux onglets analytiques ignoraient la longueur du plan

`synthesize` et `synthesize_chrono` n'utilisent ni `longueur_plan` ni
`rayon_coin` : ils cherchent la tangence sans borne d'abscisse, donc supposent
le plan assez long. Le trace, lui, utilisait `0.6*r_max + 0.5`, sans rapport
avec le plan reel : sur l'onglet *Coeur simple* le contact sortait de 0.45 mm
du plan dessine sur 28 positions sur 361, et le marteau quittait visiblement
le coeur.

Les deux syntheses calculent maintenant `demi_requis`, la demi-longueur de
plan que le profil exige reellement. Le trace prend le plan reel quand il est
donne, sinon une longueur qui couvre la tangence. Et le rapport annonce la
longueur exigee, avec un avertissement explicite quand le plan reel est plus
court : le couple calcule ne vaut alors plus, et seul l'onglet d'optimisation
modelise le plan fini.

### `--nsga` lancait l'interface graphique

L'option n'est implementee que dans `cli_chrono`. Sans `--chrono` elle tombait
sur `run()` : la fenetre s'ouvrait, aucun calcul ne demarrait, et le processus
restait a 0 % de CPU. `--nsga` implique desormais `--chrono`.

### L'export DXF de l'onglet d'optimisation ecrivait le mauvais coeur

`_on_export` ne connaissait que les onglets *Coeur simple* et *Chronographe* :
sur l'onglet d'optimisation il tombait sur `write_dxf(self.res, ...)` et
ecrivait le coeur de la **synthese analytique**, pas le profil affiche. Le
fichier sortait avec les calques `COEUR`/`BRAS`/`PIVOT_MARTEAU` de la synthese
au lieu des 12 positions de cinematique. Le bouton contextuel aiguille
maintenant vers `_on_export_profil`.

Controle sur le fichier produit : 13 contours de marteau de 22 points chacun
(plan plus les deux arrondis de 0.1 mm), aucun point de marteau a l'interieur
du coeur, distance mini marteau / coeur 5e-6 mm.

### Invariants verifies

Un controle systematique passe maintenant sur tous les chemins de calcul
(profil importe, reconstruction a bande limitee, synthese de chronographe,
plan infini, synthese analytique), moteur numba et Python pur :

* contact effectif a chaque position, sans penetration ;
* marteau rigide (rayon pivot -> face constant a 3.6e-15 mm) ;
* un seul cote du coeur sur tout le tour ;
* effort normal fini et positif ;
* au plus deux changements de signe du couple sur un tour ;
* accord numba / Python pur a 1.8e-15 sur phi, N et T.

### Le centre des arrondis de bout etait pose sur le plan

`scan_plan_fini` et `_evaluer_plan_fini` placaient le centre du quart de
cercle de coin **sur** le plan de contact. L'arc avancait donc de `rc` devant
le plan, vers le coeur : le coin touchait a chaque position et le plan ne
pouvait jamais porter. Le centre est maintenant rentre de `rc` dans la matiere
(+u0), l'arrondi est tangent au plan et ne le depasse pas. Les deux moteurs
(numba et Python pur) donnent des resultats identiques au bit pres.

`core.contour_plan()` fournit le contour **materiel** du plan, arrondis
compris ; il sert au trace (`FlatFace.segment`, `profil.pose`) et a l'export
DXF, qui reduisait auparavant la face a une corde. Convention : `longueur_plan`
est la longueur de la **partie droite**, chaque bout ajoutant `rayon_coin`.

### Le bloc *Fichiers DXF* est en tete du panneau

Il y etait en bas, sous les onglets. Or le `QTabWidget` des modes annoncait la
hauteur de sa page la plus chargee (l'optimisation, 985 px) quel que soit
l'onglet affiche : le panneau atteignait 2240 px, et les boutons *Importer* /
*Exporter* tombaient a y = 1610 px, hors de la fenetre. Deux corrections :

* `OngletsCompacts` redefinit `sizeHint` et `minimumSizeHint` pour ne retenir
  que la page affichee ;
* le bloc *Fichiers DXF* est insere en position 0, donc toujours visible.

Les champs *Couple marteau C0* et *Pente k* du panneau du haut changent de
role selon l'onglet :

| onglet | role de C0 et k |
|---|---|
| Coeur simple | **donnee**, champs actifs |
| Coeur de chronographe | **sortie**, champs grises, prefixe `=` |
| Optimisation du profil | **sortie**, champs grises, prefixe `=` |

En mode sortie ce ne sont que des afficheurs : la valeur n'est pas ecrite
dans les parametres. Sans quoi la synthese analytique repartirait du ressort
calcule au tour suivant (son profil depend de C0) et le resultat deriverait a
chaque relance. En revenant sur l'onglet *Coeur simple*, les valeurs saisies
par l'utilisateur sont restituees.

Le ressort complet se lit dans le rapport, section *Ressort de marteau
(sortie)*.

## Comment le coeur est calcule

Deux methodes coexistent, et elles ne repondent pas a la meme question.

### 1. Synthese analytique inverse (core.py, chrono.py)

Le profil est une **sortie**, pas une entree. On impose `couple = constante`,
ce qui rend `s` explicite, puis la coherence `dH/dtheta = s.(dbeta/dtheta)`
donne une equation differentielle sur l'angle du marteau que l'on integre en
RK4. Le profil se reconstruit par `X(beta) = H.u(beta) + s.t(beta)`.

Il ne reste que 4 ou 5 parametres (pente du plan, pivot, pente de couple,
niveau de couple). C'est exact, rapide, et le couple constant est garanti par
construction. Mais on ne peut **pas arbitrer** entre des criteres
contradictoires : on ne peut que rejeter les solutions qui ne conviennent pas.

### 2. Profil discretise + NSGA-II (profil.py, nsga2.py)

Le profil devient la variable. La representation n'est ni le rayon polaire ni
une suite de points, mais le **rayon de courbure en fonction de l'angle de
normale** :

    rho(beta) = H(beta) + H''(beta)

parce que les deux proprietes dures du probleme y deviennent triviales :

- **convexite** : `rho >= 0`, une simple borne sur chaque variable. C'est
  indispensable ici, le marteau a une face plane donc il ne touche que
  l'enveloppe convexe ;
- **fermeture** du contour : les harmoniques de rang 1 de rho sont nulles. Il
  suffit de ne pas les mettre dans le jeu de variables (une projection les
  retire apres chaque croisement ou mutation).

`rho` est une **serie de Fourier a bande limitee**, et c'est exactement la
courbure de la courbe reconstruite : la contrainte `rho >= r_outil` porte
donc sur la geometrie que l'on usinera.

Une premiere version stockait `rho` par cellules constantes par morceaux. Les
valeurs de cellules etaient bien positives, mais la reconstruction a bande
limitee sonnait autour du Dirac que represente un plat du profil : la
courbure reelle descendait a **-1.2 mm**, donc les profils n'etaient pas
convexes. Un plat exact n'existe pas dans un espace a bande limitee, et c'est
de toute facon ce qu'interdit la contrainte d'usinabilite. La graine
analytique, qui a un plat au zero, en est donc l'approximation la plus proche
qui soit usinable.

La reparation d'un jeu de coefficients se fait par **projections alternees**
entre l'ensemble convexe `{rho >= r_outil}` et le sous-espace a bande limitee
sans harmonique de rang 1, avec un relevement final du terme constant qui
garantit la borne.

Les deux harmoniques de rang 1 de H ne changent pas la forme : elles
**deplacent le centre de rotation** dans le profil. Ce sont donc exactement
les deux degres de liberte qui comptent pour une came, et ils sont separes de
la forme. Jeu de variables : `J` courbures + `cx, cy`.

L'analyse est **directe** : pour chaque position du coeur on cherche la
tangence du plan du marteau, on en deduit le bras de levier, le couple dans
les deux sens de glissement et le blocage eventuel. Environ 15 ms par profil.

La population de depart est semee avec la solution analytique, discretisee :
l'optimiseur part d'un coeur valide et l'ameliore, au lieu de chercher la
convexite et la fermeture au hasard.

### Le galet : le seul vrai levier sur la plage morte

La plage morte au sommet vaut `2.atan(mu)` **ou mu est le frottement au
contact**, et la mesure confirme qu'elle ne depend presque que de lui :

- le rapport `dtheta/dbeta` mesure **0.95** : la fenetre de frottement en
  normale se transporte telle quelle en angle de coeur ;
- la distance au pivot n'y change rien (15 a 16 deg de |OP| = 8 a 40 mm) ;
- la forme du coeur non plus : 2400 formes libres explorees par NSGA-II ne
  descendent pas sous ce plancher ;
- une face de marteau convexe la divise au mieux par deux (voir ci-dessous).

Un **galet monte fou** au bout du marteau remplace le glissement par du
roulement. Il ne reste au contact que le frottement de l'axe du galet, ramene
par le rapport des rayons :

    mu_contact = mu_palier . r_axe / R_galet  +  f_roulement / R_galet

Mesures (axe de galet 0.15 mm, frottement de palier 0.12, 1440 pas) :

| bout du marteau | mu au contact | plage morte | couple coeur | N maxi |
|---|---|---|---|---|
| plat glissant | 0.120 | **15.5 deg** | 0.189 | 0.357 |
| galet R = 1.2 | 0.0192 | 1.00 deg | 0.210 | 0.356 |
| galet R = 1.6 | 0.0144 | 0.75 deg | 0.214 | 0.332 |
| galet R = 2.0 | 0.0115 | 0.75 deg | 0.217 | 0.316 |
| galet R = 3.0 | 0.0077 | 0.75 deg | 0.222 | 0.294 |

**La plage morte est divisee par 15 a 30**, et au passage le couple sur le
coeur augmente de 15 % pendant que l'effort de contact baisse de 18 %. Tous
les controles de coherence passent et la course du marteau reste a 9.5 deg.

Elle reste proportionnelle au frottement du palier du galet :

| mu palier | 0.05 | 0.10 | 0.15 | 0.25 |
|---|---|---|---|---|
| plage morte [deg] | 0.25 | 0.75 | 0.75 | 1.25 |

Case **Galet monte fou** dans l'onglet d'optimisation (elle libere
automatiquement la courbure de la face, puisque le rayon du galet devient une
variable et que `mu_contact` en depend). En ligne de commande : `--galet`.

### Le plat du marteau comme variable

La borne `2.atan(mu)` demontree plus haut porte sur la **normale de contact**.
Avec une face plane, la normale est rigidement liee a l'angle du marteau :
`dtheta/dbeta` vaut environ 1, donc la fenetre de frottement se transporte
telle quelle en angle de coeur. **Une face courbe decouple les deux.**

L'analyse se fait sur la **developpee** du coeur : le centre du cercle de la
face est en `Y(beta) = (H + R).u + s.t` dans le repere du coeur, et comme il
est a distance fixe du pivot, chaque position donne l'equation scalaire
`|Y(beta) - P_theta| = L`. Parmi ses racines, le ressort arrete le marteau au
**premier contact**, donc a la plus petite ouverture sans interpenetration.

Mesures a **couple marteau fige** (C0 = 4 N.mm, k = 2 N.mm/rad), profil issu
de la synthese analytique, r_min 0.3, r_max 3, mu 0.12 :

| R face [mm] | plage morte [deg] | couple coeur [N.mm] | N maxi [N] | changements de cote |
|---|---|---|---|---|
| +0.8 (convexe) | 5.0 | 0.158 | 0.352 | **53 — invalide** |
| +1.0 | 4.0 | 0.165 | 0.334 | **16 — invalide** |
| **+1.3** | **6.0** | 0.167 | 0.314 | 0 |
| +2.0 | 7.0 | 0.170 | 0.288 | 0 |
| +3.0 | 9.0 | 0.174 | 0.287 | 0 |
| +6.0 | 11.0 | 0.179 | 0.299 | 0 |
| plane | 16.0 | 0.189 | 0.357 | 0 |
| -10 (concave) | 21.0 | 0.211 | 0.589 | 0 |
| -30 | 17.0 | 0.196 | 0.397 | 0 |

**Une face convexe de rayon de l'ordre de r_max divise la plage morte par
deux** (16 deg a plat, 6 a 7 deg a R = 1.3 a 2 mm), en reduisant au passage
l'effort de contact de 20 %. Le prix est environ 10 % de couple coeur.

En dessous de R = 1.3 mm le marteau **change de cote du coeur** en cours de
course : la plage morte affichee descend a 4 deg mais la solution n'est pas
mecaniquement realisable, et le logiciel la rejette.

Une face **concave** ne sert a rien : elle degrade la plage morte et double
l'effort.

### Controles de coherence mecanique

Une premiere version de ce calcul donnait des faces concaves a plage morte
nulle. C'etaient des artefacts, elimines par trois controles, verifies a
chaque evaluation :

- **le ressort du marteau ne peut pas s'inverser.** Le couple applique garde
  un sens fixe sur tout le tour ; une version intermediaire le liait au signe
  de l'angle du marteau, ce qui faisait changer le couple de sens en pleine
  descente de flanc ;
- **le marteau ne change pas de cote du coeur.** L'angle du marteau doit
  garder un signe (seuil relatif a la course, parce qu'il passe par zero au
  repos). Selectionner la racine de plus grande ouverture faisait sauter le
  contact de l'autre cote du coeur : le marteau restait ouvert de 26 deg
  meme au zero ;
- **le couple ne change de sens que deux fois par tour**, au zero (equilibre
  stable) et au sommet (instable). Un changement de sens en pleine descente
  de flanc voudrait dire que le coeur repart en arriere, impossible s'il doit
  ramener son plat sur celui du marteau.

S'y ajoutent le controle d'interpenetration de la face et la continuite de la
position du marteau d'un pas a l'autre.

### Objectifs et contraintes

### Deux etapes : ressort calcule, puis ressort fige

1. **Ressort calcule.** Le coeur initial (synthese analytique) donne son
   ressort : prearmage et raideur, par ajustement affine du couple marteau
   exige. C'est le dimensionnement.
2. **Ressort fige.** On garde ce C0 et ce k, et on fait varier le coeur et
   le marteau pour ameliorer la version initiale. Toutes les solutions sont
   alors comparees **a armement identique**, seul moyen de savoir laquelle
   est reellement meilleure.

Case *Figer le ressort sur la version initiale* dans l'onglet
d'optimisation, ou `--figer` en ligne de commande. Resultat type
(C0 = 3.669 N.mm, k = 32.01 N.mm/rad figes, galet) :

| | part du tour sous la consigne | plage morte | effort |
|---|---|---|---|
| version initiale | 49.2 % | 16.00 deg | 0.758 N |
| **apres optimisation** | **26.4 %** | **1.00 deg** | 0.566 N |

### Les criteres de couple portent sur la partie utile

Le couple s'annule **forcement** aux deux equilibres : le zero (stable) et le
sommet (instable). Les inclure dans les criteres les rend trompeurs :

- une **moyenne** correcte peut cacher un couple qui s'effondre sur une bonne
  partie du parcours, et l'optimiseur exploite exactement cette faille : il
  rend le coeur large et presque circulaire, ce qui baisse l'effort de
  contact sans rien couter sur la moyenne ;
- l'**ondulation** rapportee a la moyenne vaut mecaniquement 150 a 200 %
  quel que soit le profil, puisque le minimum est nul par construction.

On ecarte donc une fenetre de 20 deg autour de chaque changement de signe, et
sur le reste (environ 75 % du tour) on mesure le **couple mini garanti** :
c'est lui qui doit valoir le couple vise. L'ondulation calculee sur cette
meme partie tombe a 65 a 80 %, et devient un nombre qui veut dire quelque
chose.

Le front est classe par manque de couple d'abord : une solution sans plage
morte mais sans couple n'est pas utilisable, et elle ne doit pas se presenter
en premier.

### Le couple vise doit etre atteignable

Le rapport de transmission coeur/marteau vaut ici environ **0.04**. Avec
C0 = 4 N.mm le coeur ne peut delivrer que 0.15 a 0.18 N.mm, quel que soit le
profil : tout le front affiche alors 48 % de manque. Ce n'est pas un defaut
du profil mais du **ressort du marteau**.

| C0 [N.mm] | couple mini garanti | manque |
|---|---|---|
| 4 | 0.182 | 39 % |
| 6 | 0.273 | 9 % |
| **7** | **0.309** | **0 %** |
| 8 | 0.364 | 0 % |

Le rapport signale le facteur a appliquer a C0 quand le couple est
insuffisant. Avec C0 = 7 N.mm, k = 3.5 et un galet de 4.3 mm, la solution
retenue donne **0.309 N.mm garantis sur 77 % du tour, 0 deg de plage morte**,
65 % d'ondulation, le couple ne tombant a zero qu'aux deux equilibres.

### Le couple marteau est une donnee

Le **couple du marteau est une donnee** : C0 et sa pente k sont identiques
d'une simulation a l'autre, sinon les solutions ne sont pas comparables. Le
couple obtenu sur le coeur est donc une **sortie**, et c'est son ecart a la
valeur visee qui est un objectif. (Une version precedente remettait les
couples a l'echelle solution par solution, ce qui rendait la comparaison
vide : le couple coeur valait toujours la cible par construction.)

**Objectifs** (4) : **part du tour sous la consigne**, plage morte,
ondulation du couple, effort de contact.

Le premier a remplace un critere fonde sur le couple moyen puis sur un
minimum fenetre. Les deux etaient piegeux :

- le **couple moyen** masque un couple qui s'effondre : l'optimiseur rendait
  le coeur large et presque circulaire, ce qui baisse l'effort de contact
  sans rien couter sur la moyenne ;
- le minimum calcule **hors des equilibres** ecartait du critere *tout*
  echantillon a couple nul. Plus la plage morte etait longue, plus elle etait
  cachee : l'optimiseur etait litteralement recompense de l'allonger. C'est
  le defaut qui faisait deriver les solutions vers de longues plages a couple
  nul.

La part du tour sous la consigne n'a aucune fenetre arbitraire : une longue
plage a couple faible la fait monter, point.

**Contraintes** : encombrement, rayon mini, exactement deux changements de
signe du couple sur un tour (le zero stable et le sommet instable), et
**rayons de raccord usinables**.

Cette derniere est gratuite dans cet espace de conception : `rho_j` *est* le
rayon de courbure local, et un angle vif c'est `rho = 0`. Imposer un rayon
d'outil revient donc a relever la **borne inferieure des variables** de 0 a
`r_outil`. Le profil etant convexe il n'a pas d'angle rentrant, la seule
geometrie non usinable est l'angle vif, et elle est interdite par
construction. La projection de fermeture respecte cette borne.

Cout de la contrainte, mesure (graine analytique, r_outil impose) :

| r_outil [mm] | rho mini obtenu | plage morte [deg] | C0 [N.mm] |
|---|---|---|---|
| 0 (libre) | 0.005 | 13.50 | 6.098 |
| 0.05 | 0.050 | 13.50 | 6.115 |
| 0.08 | 0.080 | 13.50 | 6.127 |
| 0.15 | 0.150 | 13.50 | 6.154 |
| 0.30 | 0.300 | 16.50 | 6.159 |

Jusqu'a 0.15 mm elle ne coute presque rien (+0.9 % de couple marteau). A
0.30 mm la pointe doit etre trop arrondie et la plage morte augmente.
Le rayon est pris dans le parametre `Rayon de raccord du plat`
(`--rayon-raccord`).

Pilotage depuis l'interface, onglet **Optimisation du profil** :

- reglages : nombre de cellules de courbure, population, generations, pas
  d'analyse ;
- les quatre objectifs sont cochables (au moins deux) ; les contraintes
  restent toujours actives ;
- une case **Optimiser aussi le plat du marteau** ajoute la courbure de la
  face aux variables (evaluation environ 3 fois plus lente) ; la vue trace
  alors l'arc de la face au lieu du segment ;
- le calcul tourne dans un fil, l'avancement s'affiche dans la barre d'etat ;
- le **front de Pareto** s'affiche en liste et en nuage de points (onglet
  *Pareto* des graphiques). Cliquer une ligne ou un point affiche le profil
  correspondant dans la vue cinematique, avec ses courbes de couple, d'effort
  et de courbure, et son rapport ;
- la solution choisie s'exporte en DXF.

En ligne de commande :

    python3 main.py --chrono --angle-plan 65 --nsga

Resultat (48 cellules, population 60, 35 generations, 31 s, r_outil 0.08 mm) :

| plage morte [deg] | ondulation [%] | C0 [N.mm] | N [N] |
|---|---|---|---|
| 13.50 | 113.6 | 6.030 | 0.543 |
| 15.00 | 105.1 | 6.049 | 0.544 |
| 21.00 | 223.2 | 5.396 | 0.482 |
| 24.00 | 220.1 | 5.271 | 0.460 |
| 30.00 | 211.7 | 5.108 | 0.456 |

On peut donc baisser le couple marteau de 15 % et l'effort de contact de 16 %
en acceptant une plage morte deux fois plus large. Et **aucune des formes
libres explorees ne descend sous 13.5 deg de plage morte**, ce qui confirme
numeriquement la borne `2.atan(mu) = 13.7 deg` etablie analytiquement, cette
fois sans hypothese de forme.

Verification d'une solution du front : rho mini 0.0800 mm (exactement la
borne), rayon du profil 0.3017 a 2.9974 mm (dans les limites 0.30 / 3.00),
2 changements de signe, defaut de fermeture 3e-13, convexite 48/48 cellules.

## Parametres

| Parametre | Role |
|---|---|
| Pivot marteau X / Y | centre de rotation du marteau |
| Rayon mini du coeur | distance mini au centre de rotation (position de repos) |
| Encombrement maxi | rayon maxi admissible du coeur |
| Couple marteau C0 | couple mini, coeur en place |
| Pente marteau k | accroissement du couple par radian d'ouverture |
| Frottement dynamique | coefficient mu |
| Pente du plan | angle de la droite O-P a la **normale** de contact, au repos |
| Longueur du plan | longueur du plan de contact ; 0 = deduite de l'encombrement |
| Rayon du plan | 0 = plan ; > 0 = face convexe |
| Rotation synthetisee | course du coeur simple (350 deg pour un coeur a un lobe) |
| Nombre de lobes | repetition du profil (coeur simple) |

Parametres propres au coeur de chronographe :

| Parametre | Role |
|---|---|
| Couple vise sur le coeur | fixe le niveau des couples ; le C0 marteau en decoule |
| Corde de la pointe | troncature du sommet ; 0.05 = pointe vive (optimal) |
| Rayon des pointes | arrondi d'usinage des deux sommets |
| Rayon du creux de pointe | arc concave entre les deux sommets (trace) |
| Plat souhaite | valeur de controle ; la longueur reelle est imposee |
| Marge anti-blocage | facteur de securite souhaite sur L / (mu.r_min) |
| Diedre des deux pans | angle entre les deux pans du marteau |

Les courses des deux flancs ne figurent pas dans ces listes : elles sont
calculees.

## Optimisation (coeur simple)

Nelder-Mead multi-departs sur les variables cochees (pente du plan, pente du
marteau, rayon du plan, position du pivot), avec bornes reglables. Critere
pondere : augmentation du couple marteau, couple utile du coeur, ecart entre
frottement rentrant et sortant, effort de contact, course du marteau ;
penalite sur les profils non valides et sur l'ecart a l'encombrement cible.

## Comment se mesure la pente du marteau

L'angle saisi est mesure **depuis la direction O->P** (du centre du coeur
vers le pivot du marteau) **jusqu'a la normale** du plan de contact, et
**au repos**, coeur au zero.

Trois choix, chacun pour une raison :

- **rapporte a OP** et non a l'horizontale ou a la verticale : ainsi il ne
  depend pas de l'endroit ou l'on place le pivot autour du coeur. C'est lui
  qui fixe le bras de levier au pivot, `q = |OP|.sin(pente)`, donc le rapport
  de transmission ;
- **sur la normale** et non sur le plan : toute la mecanique du contact
  s'ecrit sur la normale. Le plan lui-meme est a 90 deg de plus ;
- **au repos** : en cours de course le marteau tourne de `phi`, la normale
  devient `gamma0 + phi` et l'orientation absolue du plan change d'autant.
  C'est exactement ce qu'on observe quand le coeur s'arme.

Exemple avec le pivot en (16, 5) et une pente de 65 deg :

| | deg / horizontale |
|---|---|
| direction O->P | 17.35 |
| normale au repos | 82.35 |
| **plan au repos** | **172.35** |
| plan a pleine ouverture | 161.94 |
| rotation du plan sur un tour | 10.42 |

Le rapport donne ces quatre valeurs : ce sont elles qu'on reporte en CAO, la
pente saisie n'etant qu'un parametre interne commode.

## Le plan du marteau a une longueur finie

Un plan infini est une commodite de calcul. Si le point de tangence sort du
plan reel, ce n'est plus lui qui touche le coeur mais son **coin**,
eventuellement arrondi au rayon `rayon_coin`. Le marteau se ferme jusqu'au
premier obstacle : parmi les trois contacts possibles - le plan et chacun de
ses deux coins - c'est celui qui exige la **plus grande ouverture** qui
s'impose.

Le contact sur un coin est un contact ponctuel sur le flanc du coeur : la
normale est celle du coeur au point touche, pas celle du plan. Un coin
arrondi se traite comme un galet de petit rayon.

Effet, sur le coeur analytique de reference :

| longueur du plan | contacts sur coin | plage morte | N maxi |
|---|---|---|---|
| infini | 0 / 360 | 16.0 deg | 0.36 N |
| 6 mm | 30 / 360 | 16.0 deg | 7.49 N |
| 4 mm | 87 / 360 | 16.0 deg | 6.34 N |
| 3 mm | 172 / 360 | 7.0 deg | 2.53 N |
| 2 mm | 220 / 360 | 0.0 deg | 8.62 N |

Et sur un coeur de minuterie importe, dont le plan mesure 1.119 mm :

| | contacts sur coin | plage morte | N maxi | C0 requis |
|---|---|---|---|---|
| plan suppose infini | 0 / 360 | 14.0 deg | 1.07 N | 7.50 N.mm |
| **plan reel 1.119 mm** | **253 / 360** | 2.0 deg | 6.76 N | 9.21 N.mm |
| plan reel, coins r = 0.05 | 285 / 360 | 1.0 deg | 5.39 N | 9.42 N.mm |

**70 % du tour se fait sur un coin** : l'hypothese du plan infini etait
fausse pour cette geometrie, et elle sous-estimait l'effort de contact d'un
facteur 6. Le rapport indique le nombre de positions concernees.

Le **rayon des coins est essentiel**. Un coin vif concentre la charge ; un
rayon de 0.1 mm la repartit. Sur le coeur de minuterie importe :

| | contacts sur coin | plage morte | N maxi | C0 requis |
|---|---|---|---|---|
| plan suppose infini | 0 / 360 | 14.0 deg | 1.07 N | 7.50 N.mm |
| plan reel, coin vif | 253 / 360 | 2.0 deg | **6.76 N** | 9.21 N.mm |
| **plan reel, coin r = 0.1** | **304 / 360** | **0.0 deg** | **1.05 N** | 9.41 N.mm |

Le coin vif multiplie l'effort de contact par 6.4 ; le rayon de 0.1 mm le
ramene au niveau du plan infini tout en supprimant la plage morte.

Le champ *Longueur du plan* vaut 0 par defaut, ce qui conserve l'hypothese
du plan toujours assez long ; *Rayon des coins* vaut 0.1 mm. L'import DXF
renseigne la longueur a partir du dessin.

Ce calcul est **compile** (voir plus haut) : 2.5 ms par evaluation au lieu de
285 ms en Python pur, resultats identiques. Une optimisation complete avec
plan fini prend 3 s au lieu de 1 s. Le profil importe en beneficie aussi :
seul son echantillonnage reste en Python, le contact passe par le noyau.

## Cinematique du marteau

Le marteau est un **solide** : sa face est a distance fixe du pivot et ne
peut que tourner autour. Le tracé la construisait a partir du pied de la
perpendiculaire issue de O, qui **glisse** le long de la face quand le
marteau tourne : les extremites dessinees n'etaient donc pas des points
materiels, et le marteau semblait se translater. La face est maintenant
definie une fois au repos puis tournee rigidement autour du pivot ; la
distance de ses points au pivot est constante a 7e-15 mm pres sur un tour.

Le coeur, lui, tourne autour de O : le spectre de ses rayons est rigoureusement
invariant sur un tour. Un coeur etant excentre, sa masse balaie beaucoup en
tournant, mais il n'y a pas de translation.

## DXF

- **Import** : le dessin est lu tel qu'un CAO l'exporte. Le coeur et le
  marteau y sont des **suites d'arcs et de segments**, pas des polylignes, et
  les cercles pleins reperent les axes. Le lecteur rechaine donc les contours
  par proximite de leurs extremites, puis :

  - le **coeur** est le contour qui entoure son axe de rotation en restant a
    distance de lui de tous les cotes. Un cercle plein marque souvent cet
    axe ; **a defaut c'est l'origine** qui en tient lieu, suivant la
    convention du dessin ;
  - le **marteau** est l'autre contour. Son pivot est le cercle plein
    eloigne, ou, s'il n'y en a pas, **le point du marteau le plus eloigne du
    coeur**, la ou aboutissent ses bras ;
  - le **plan de contact** est relu sur les **entites** du marteau et non
    sur les points discretises, sinon un arc y produit de minuscules
    segments. On ecarte les droites qui aboutissent au pivot - ce sont les
    bras - et on retient la plus proche du coeur. Sa **longueur est reprise
    du dessin** ;
  - les **arcs qui encadrent ce plan dans la chaine** donnent le **rayon des
    coins**, qui est determinant pour l'effort de contact ;
  - tout est ramene dans le repere du coeur, axe a l'origine. L'axe n'a donc
    pas besoin d'etre en (0, 0) dans le dessin.

  **Le contour importe est simule tel quel** : sa fonction d'appui est
  calculee directement sur les points, sans passer par la reconstruction a
  bande limitee. On obtient donc le comportement de la geometrie dessinee
  (couple, plage morte, efforts) et le ressort qu'elle reclame. La
  reconstruction ne sert plus que de **graine a l'optimisation**, qui a
  besoin d'un vecteur de coefficients ; le rapport donne l'ecart entre les
  deux.

  Les rayons mini et maxi du coeur, le pivot, la pente et la longueur du
  plan sont repris du dessin et remplacent les valeurs courantes.

  Le contour lu est aussi **affiche en rose par-dessus** le profil simule.
  Le marteau ayant une face plane, il ne voit que l'enveloppe convexe du
  dessin : un coeur tronque ou non convexe s'en ecarte, et c'est visible
  directement dans la vue.

  Deux dessins du meme mecanisme, l'un avec des cercles d'axe et l'autre
  sans mais avec les rayons de coin, donnent maintenant **le meme
  resultat** : rayons 0.556 a 2.318 mm, pivot a 14.72 mm, plan de 1.119 mm
  a -46.59 deg, coins de 0.100 mm. Avec son
  ressort ajuste (C0 = 7.50 N.mm, k = 0) il donne 0.305 N.mm de couple
  moyen, 14.0 deg de plage morte et 1.07 N d'effort, et 49 % du tour reste
  sous la consigne de 0.300 N.mm.

  A defaut d'arcs, l'ancienne lecture sert de secours : deux cercles
  concentriques pour les rayons, un cercle eloigne pour le pivot, une droite
  pour le plan.

- **Export** : l'extension `.dxf` est ajoutee si elle manque dans le nom
  saisi. Le coeur sort dans **son** repere (position d'usinage) et le marteau
  y est ramene : sans cela le plan de contact traversait le coeur de
  1.7 mm sur le plan, puisque le coeur etait dessine a theta = 0 et le
  marteau a la position du curseur. La tangence est maintenant exacte
  (penetration 0.000000 mm a toute position).

  Calques : `COEUR`, `CONSTRUCTION`, `PIVOT_MARTEAU`, `PLAN_CONTACT`, `BRAS`,
  et `MARTEAU_CINEMATIQUE` qui trace douze positions du marteau sur un tour.
  Dans le repere du coeur le pivot du marteau decrit un cercle autour de O :
  c'est bien le mouvement **relatif** du marteau vu du coeur. Pour le coeur
  de chronographe s'y ajoutent `PLAT_ZERO`, `CORDE_POINTE` et
  `CREUX_POINTE`.

## Organisation

```
main.py            point d'entree (GUI / CLI)
came/core.py       modele, frottement, synthese d'un flanc a couple constant
came/chrono.py     coeur de chronographe : deux flancs, plat, raccords, pointe
came/profil.py     profil discretise (courbure par cellules), analyse directe
came/nsga2.py      NSGA-II et le probleme multi-objectif du coeur
came/optimize.py   critere pondere et Nelder-Mead (coeur simple)
came/dxfio.py      lecture / ecriture DXF
came/ui.py         interface PyQt6 (vue cinematique, courbes, rapport)
```
