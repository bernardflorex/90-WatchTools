# outils/contenu_theorie.py
"""
Contenu du document de théorie, sous forme de données. Voir generer_theorie.py
pour le rendu Word et pour la raison d'être de cette séparation.

Chaque bloc est un tuple dont le premier élément est le type :
    ("h1", titre) ("h2", titre) ("p", texte[, amorce_en_gras])
    ("puce", texte) ("eq", [lignes]) ("tab", entetes, lignes[, poids])
    ("saut",)

RÈGLE DE MISE À JOUR : toute valeur numérique citée ici provient d'une
exécution réelle du programme, jamais d'une estimation. Quand le modèle change,
c'est cette règle qui coûte cher à tenir -- et c'est elle qui fait l'intérêt du
document pour celui qui doit recouper les résultats.
"""

_PARTIE_1 = [

    # ================================================================ 1
    ("h1", "1. Objet du document"),
    ("p", "Ce document décrit le modèle physique et les méthodes numériques du "
          "programme Spiral Lab, qui optimise la géométrie d'un spiral plat en silicium "
          "usiné par DRIE. Il est écrit pour permettre la reconstruction indépendante du "
          "modèle dans un autre outil de calcul, et la comparaison poste par poste des "
          "résultats."),
    ("p", "Il est organisé en quatre parties : la théorie physique (§2 et §3), les "
          "outils mathématiques employés (§4), les choix de modélisation et leur "
          "justification chiffrée (§5 et §6), et la structure du programme fichier par "
          "fichier (§7). Le §8 propose une liste ordonnée de points de comparaison pour "
          "la validation croisée, du plus élémentaire au plus intégré."),
    ("p", "Toutes les valeurs numériques citées comme « mesurées » proviennent "
          "d'exécutions réelles du programme, et non d'estimations. Les conditions de "
          "mesure sont précisées à chaque fois."),
    ("p", "Le critère final visé est le DVH du contrôle COSC — l'écart de marche entre "
          "la position 6H et la position horizontale. Il dépend de l'ORIENTATION du "
          "spiral dans le mouvement, et pas seulement de sa forme : le §3.3 fixe donc le "
          "repère de façon explicite. Quiconque veut recouper un chiffre de DVH doit "
          "d'abord s'assurer qu'il travaille dans le même repère, faute de quoi il "
          "obtiendra la bonne valeur au signe près.",
          "Avertissement d'orientation. "),

    # ================================================================ 2
    ("h1", "2. Le problème physique"),
    ("h2", "2.1 L'oscillateur"),
    ("p", "Un balancier d'inertie J est rappelé par un spiral plat. L'ensemble forme un "
          "oscillateur mécanique dont la fréquence propre vise ici 3,5 Hz. Le spiral est "
          "une lame de silicium de hauteur hors-plan h constante (100 µm, imposée par la "
          "gravure DRIE), d'épaisseur dans le plan e(s) librement variable le long de "
          "l'abscisse curviligne s, encastrée au piton (extrémité extérieure, fixe par "
          "rapport au bâti) et au collet (extrémité intérieure, solidaire de l'axe du "
          "balancier)."),
    ("p", "La montre mécanique est un oscillateur de précision : une erreur de marche de "
          "1 s/j correspond à un écart relatif de période de 1,16·10⁻⁵. Toute "
          "simplification du modèle doit être jugée à cette échelle, ce qui est le fil "
          "conducteur des choix décrits au §5."),

    ("h2", "2.2 Les sources d'erreur de marche"),
    ("p", "La marche diurne observée au chronocomparateur est la somme de contributions "
          "de natures différentes. Leur séparation n'est pas un choix de présentation : "
          "elle commande toute l'architecture du programme, parce qu'elles ne se "
          "corrigent pas avec les mêmes leviers. Les quatre premières dépendent de "
          "l'amplitude ; la cinquième, de la température."),
    ("tab",
     ["Contribution", "Origine physique", "Dépend de la position ?", "Dépend de l'amplitude ?"],
     [["Gravité",
       "Déplacement du centre de gravité du spiral sous son propre poids pendant l'oscillation",
       "OUI — maximale en position verticale, rigoureusement NULLE à plat", "Oui"],
      ["Isochronisme propre",
       "Non-linéarité du couple de rappel du spiral en grande déformation",
       "NON — identique dans les six positions", "Oui"],
      ["Échappement",
       "Couple parasite de l'échappement (repos, tirage, impulsion, dégagement)",
       "NON — identique dans les six positions",
       "Oui, d'autant plus fort que l'amplitude est faible"],
      ["Jeu de raquette",
       "La lame quitte l'appui d'une goupille et traverse le jeu : le spiral actif "
       "s'allonge alors jusqu'au piton",
       "NON — identique dans les six positions",
       "Oui, d'autant plus fort que l'amplitude est faible"],
      ["Température", "Le module d'Young du silicium chute de 60 ppm par kelvin", "NON",
       "Non — dépend de la température, pas de l'amplitude"]],
     [16, 34, 26, 24]),
    ("p", "En position horizontale (cadran haut ou cadran bas), la gravité est "
          "perpendiculaire au plan du spiral : son couple perturbateur autour de l'axe du "
          "balancier est nul. Le terme de gravité y disparaît exactement, alors que les "
          "deux autres subsistent.",
          "Conséquence décisive. "),
    ("p", "Il en découle que l'erreur d'échappement ne peut pas être compensée par la "
          "géométrie qui règle le déplacement du centre de gravité : celle-ci n'a aucun "
          "effet à plat, là où l'échappement agit pourtant pleinement. La compensation ne "
          "peut venir que de l'isochronisme propre du spiral, seul terme à agir dans les "
          "six positions. Le jeu de raquette a exactement la même propriété et s'ajoute "
          "donc à ce qu'il faut compenser. Les deux objectifs d'optimisation — minimiser "
          "la dispersion entre positions verticales, et compenser échappement plus jeu — "
          "portent sur des grandeurs physiquement indépendantes."),
    ("p", "Cette indépendance a été vérifiée numériquement : faire varier le retard "
          "d'échappement spécifié de 0 à 30 s/j laisse la dispersion entre positions "
          "verticales inchangée à moins de 10⁻⁹ s/j près."),
    ("p", "La même propriété donne au DVH sa signification exacte. Le DVH étant la "
          "différence entre une position verticale et une position à plat, tous les "
          "termes indépendants de la position s'y annulent deux à deux : il ne reste que "
          "le terme de gravité en position 6H. Vérifié à 5,6·10⁻¹⁶ s/j près. Optimiser le "
          "DVH, c'est donc optimiser l'orbite du centre de gravité, et rien d'autre — "
          "l'échappement, l'isochronisme propre et le jeu de goupilles n'y entrent pas.",
          "Le DVH en découle directement. "),

    ("h2", "2.3 Ce que le modèle ne traite pas"),
    ("p", "Sont explicitement hors du périmètre actuel, par décision de projet :"),
    ("puce", "l'erreur d'échappement n'est pas calculée depuis la géométrie de "
             "l'échappement (angles de levée, de repos, de tirage, couple moteur) : elle "
             "est prise comme donnée mesurée ou spécifiée, dont seule la loi de variation "
             "avec l'amplitude est modélisée ;"),
    ("puce", "le FROTTEMENT DU PIVOT du balancier dans sa pierre. C'est la seule "
             "omission qui ait une conséquence de conception : le spiral presse le pivot "
             "avec une force atteignant treize fois son propre poids, et c'est elle qui "
             "justifie la règle du nombre entier de spires (voir §3.8). Le programme "
             "calcule et affiche cette force, mais n'en déduit aucun frottement ;"),
    ("puce", "le scalloping résiduel de la gravure DRIE et les post-traitements ;"),
    ("puce", "la dynamique du choc d'impulsion et les effets de l'air ;"),
    ("puce", "l'anisotropie résiduelle du silicium : l'orientation cristalline ⟨111⟩ "
             "retenue est transversalement isotrope dans le plan, ce qui justifie un "
             "module d'Young scalaire constant E = 169 GPa."),

    # ================================================================ 3
    ("saut",),
    ("h1", "3. Théorie"),
    ("h2", "3.1 Élastica non linéaire du spiral"),
    ("p", "Le spiral est traité comme une poutre plane inextensible et non cisaillable "
          "(Euler-Bernoulli) en grandes rotations. La lame étant très élancée (épaisseur "
          "~35 µm pour une longueur développée de plusieurs dizaines de millimètres), "
          "l'allongement axial et le cisaillement transverse sont négligeables devant la "
          "flexion."),
    ("p", "Pour une rotation α du balancier, le spiral est chargé au piton par un "
          "torseur plan réduit à deux composantes de force et un moment, (Fx, Fy, M₀). "
          "L'équilibre statique d'un tronçon s'écrit, en tout point d'abscisse "
          "curviligne s :"),
    ("eq", ["M(s) = M₀ + Fx · [ Y(s) − y_piton ] − Fy · [ X(s) − x_piton ]",
            "X(s), Y(s) sont les coordonnées du point DANS LA CONFIGURATION DÉFORMÉE"]),
    ("p", "La loi de comportement relie la variation de courbure au moment fléchissant :"),
    ("eq", ["dφ/ds  =  dφ₀/ds  −  M(s) / [ E · I(s) ]        avec   I(s) = h · e(s)³ / 12"]),
    ("p", "et la cinématique du repère de Frenet ferme le système :"),
    ("eq", ["dX/ds = cos φ(s)          dY/ds = sin φ(s)",
            "φ₀(s) est l'orientation de la tangente dans la configuration au repos"]),
    ("p", "φ₀(s) est l'orientation de la tangente dans la configuration au repos : le "
          "spiral n'est pas une poutre droite, et c'est bien l'écart de courbure par "
          "rapport à sa forme libre qui produit le moment."),
    ("p", "Les conditions aux limites sont :"),
    ("puce", "au piton (s = L, extrémité extérieure) : position et orientation figées, "
             "X = x_piton, Y = y_piton, φ = φ₀(L) ;"),
    ("puce", "au collet (s = 0, extrémité intérieure) : le point est solidaire de l'axe "
             "du balancier, donc astreint au cercle de rayon r_min et tourné de α — "
             "X = r_min·cos(θ₀+α), Y = r_min·sin(θ₀+α), φ = φ₀(0) + α."),
    ("p", "Trois conditions au collet, trois inconnues (Fx, Fy, M₀) : le problème aux "
          "limites est exactement déterminé. C'est ce qui permet de le résoudre comme un "
          "système d'équations plutôt que comme une minimisation (§4.1)."),

    ("h2", "3.2 Erreur de marche due à la gravité — formule d'Airy"),
    ("p", "Le spiral se déforme pendant l'oscillation ; son centre de gravité décrit donc "
          "une orbite x_cg(φ), y_cg(φ) dans le plan, fonction de l'angle instantané φ du "
          "balancier. En position verticale, le poids du spiral exerce alors sur l'axe un "
          "couple perturbateur qui dépend de φ, donc un couple qui varie au cours de "
          "l'oscillation."),
    ("p", "Ce couple étant très faible devant le couple de rappel, son effet sur la "
          "période se traite par perturbation au premier ordre d'un oscillateur "
          "harmonique (méthode d'Airy). Pour une amplitude Θ₀ et une orientation de la "
          "gravité (g_x, g_y) :"),
    ("eq", ["M_g(φ) = m_spiral · [ g_x · dx_cg/dφ  +  g_y · dy_cg/dφ ]",
            "δf/f  =  − 1 / (2π · ω² · J · Θ₀)  ·  I",
            "avec  I = INTEGRALE[ψ = 0 → 2π]  M_g( Θ₀ · sin ψ ) · sin ψ · dψ",
            "ψ est la phase de l'oscillation ; marche (s/j) = − 86400 · δf/f"]),
    ("p", "Quatre positions verticales sont évaluées (3H, 6H, 9H, 12H), correspondant à "
          "quatre orientations du vecteur gravité dans le plan du spiral. Le premier "
          "critère minimisé par l'optimiseur est l'écart maximal entre positions, pris "
          "sur l'amplitude la plus défavorable :"),
    ("eq", ["dispersion = max sur Θ₀ [ max sur positions (marche) − min sur positions (marche) ]"]),
    ("p", "L'orbite du centre de gravité provient du solveur élastica complet (§3.1), et "
          "non d'une approximation à courbure uniforme. Ce point est important : "
          "l'approximation à courbure uniforme devient fausse dès que l'épaisseur varie "
          "le long de la lame — ce qui est précisément le degré de liberté exploité ici — "
          "et l'écart mesuré sur la dispersion de l'orbite atteignait un facteur ~50 "
          "entre les deux approches."),

    ("h2", "3.3 Repère du mouvement, positions COSC et critère DVH"),
    ("p", "Cette section fixe le repère. Le DVH change de SIGNE si l'on se trompe de "
          "demi-tour, sans qu'aucune autre grandeur du calcul ne bouge : c'est une erreur "
          "muette, et la seule façon de s'en prémunir est d'expliciter chaque convention."),

    ("p", "La marche due à la gravité en position verticale n'est pas une collection de "
          "quatre nombres indépendants : c'est la PROJECTION D'UN VECTEUR UNIQUE sur la "
          "direction de la pesanteur. La formule d'Airy étant linéaire en g, la marche "
          "obéit exactement à",
     "Le vecteur de sensibilité. "),
    ("eq", ["marche(θ) = A · cos( θ − φ )",
            "θ = direction de g dans le plan du spiral",
            "A = module de sensibilité (s/j)   φ = phase (deg)"]),
    ("p", "Vérifié sur 72 directions de g réparties sur le tour complet : résidu de "
          "7·10⁻¹⁶ s/j. Deux positions orthogonales suffisent donc à déterminer "
          "entièrement le comportement en gravité, et le couple (A, φ) dit d'un coup "
          "d'œil ce que l'orientation du spiral peut et ne peut pas faire :"),
    ("puce", "le MODULE A ne dépend que de la forme du spiral. Il est rigoureusement "
             "invariant par rotation — mesuré de 0,760261463 à 0,760261484 s/j sur quatre "
             "orientations, soit 3·10⁻⁸ en relatif. Tourner le spiral dans le mouvement "
             "ne peut donc jamais réduire sa sensibilité à la gravité ;"),
    ("puce", "la PHASE φ, elle, suit l'orientation degré pour degré. Elle décide de la "
             "répartition de l'effet entre les quatre positions, donc de la valeur du "
             "DVH, mais uniquement en le reportant d'une position sur une autre."),
    ("p", "Réduire A est le seul gain inconditionnel. Régler φ est un arbitrage : on "
          "choisit quelles positions payent."),

    ("p", "Le tracé est fait en vue CÔTÉ ROUAGE — le spiral est sous le coq, c'est de ce "
          "côté qu'on le voit — et orienté le repère 6H EN HAUT, c'est-à-dire dans "
          "l'attitude même de la position COSC 6H. La tige de remontoir, qui sort à 3H, "
          "est alors À DROITE du dessin.",
     "Repère du tracé. "),
    ("tab",
     ["Repère horaire", "Direction dans le tracé", "Position dans le tracé"],
     [["6H", "90°", "en haut"],
      ["3H", "0°", "à droite — c'est de ce côté que sort la tige"],
      ["12H", "270°", "en bas"],
      ["9H", "180°", "à gauche"]],
     [24, 26, 50]),
    ("p", "La chiralité est celle d'une vue de dos : 12-3-6-9 y tourne dans le sens "
          "ANTIHORAIRE, comme tout cadran regardé par derrière. Se tromper de côté "
          "ÉCHANGE 3H et 9H — la dispersion est inchangée, puisque c'est une permutation, "
          "mais les valeurs sont attribuées aux mauvaises positions dès qu'on les compare "
          "à un relevé de chronocomparateur."),

    ("p", "Une position COSC porte le nom du repère horaire qui se trouve EN HAUT. En "
          "position « 6H », le 6 est donc en haut et la pesanteur pointe vers le 12. La "
          "direction de g est par conséquent TOUJOURS celle du repère opposé :",
     "Nommage des positions. "),
    ("eq", ["direction de g  =  direction du repère horaire  +  180°"]),
    ("p", "Confondre les deux décale tout l'étiquetage d'un demi-tour. Dans le programme, "
          "les directions de g ne sont pas écrites à la main : elles sont DÉDUITES des "
          "positions des repères horaires par la relation ci-dessus. Tenir deux tables à "
          "jour séparément est exactement ce qui avait laissé passer un demi-tour lors de "
          "la mise au point."),

    ("p", "Le spiral n'est pas posé arbitrairement dans le mouvement. Son "
          "attache au collet se fait sur la LIGNE DES CENTRES — la droite qui joint l'axe "
          "du balancier à celui de l'ancre, dite aussi ligne droite de l'échappement. "
          "Trois données la décrivent entièrement :",
     "Orientation du spiral. "),
    ("tab",
     ["Donnée", "Rôle", "Valeur de référence"],
     [["angle_echappement_deg",
       "Angle de la ligne des centres avec l'axe de la tige de remontoir, compté dans le "
       "sens trigonométrique tel qu'on le voit côté rouage. La tige sortant à 3H, donc à "
       "0° dans le tracé, cet angle EST la direction cherchée.", "−6,2°"],
      ["collet_cote_ancre",
       "Lequel des DEUX BOUTS de la ligne des centres porte le collet. Une droite n'est "
       "définie qu'à 180° près par son angle : sans cette donnée, le spiral peut être "
       "placé en demi-tour. Vrai = le bout qui pointe vers l'ancre et la roue "
       "d'échappement.", "vrai"],
      ["sens_trigonometrique",
       "Sens de développement du spiral depuis le collet, vu côté rouage. Le sens inverse "
       "donne la pièce MIROIR, qui n'a pas la même phase φ.", "vrai"]],
     [24, 56, 20]),
    ("p", "Sur ces valeurs, le collet part à −6,2° et le point de comptage retombe à "
          "173,8° lorsque le nombre de tours est entier (§3.8). La direction de l'ancre "
          "ne dépend PAS de collet_cote_ancre : celui-ci ne dit que sur lequel des deux "
          "demi-axes le spiral vient s'attacher, pas où se trouve l'échappement."),

    ("p", "Le DVH est l'écart de marche entre la position COSC 6H et la position "
          "horizontale. Tous les termes indépendants de la position s'annulant dans cette "
          "différence (§2.2), il se réduit exactement au terme de gravité en 6H, "
          "c'est-à-dire à la projection du vecteur de sensibilité sur la direction de g "
          "en position 6H :",
     "Le critère DVH. "),
    ("eq", ["DVH  =  marche_6H  −  marche_horizontale  =  A · cos( 90° − φ )"]),
    ("p", "Il en découle une INÉGALITÉ DE STRUCTURE, utile comme contrôle : les positions "
          "6H et 12H valent ±A·sin φ et figurent toutes deux dans l'éventail des quatre "
          "verticales, donc"),
    ("eq", ["dispersion  ≥  2 · | DVH |"]),
    ("p", "Sur la configuration de référence, l'égalité est atteinte (1,9587 s/j de "
          "dispersion pour 0,9794 s/j de DVH au pire cas) : la phase tombe presque "
          "exactement sur l'axe 6H/12H, donc tout l'effet de gravité atterrit sur la "
          "position que le COSC note. C'est le cas le plus défavorable possible pour ce "
          "critère, et c'est ce que l'optimisation cherche à défaire."),
    ("p", "La dispersion ne voit que le module A ; le DVH voit AUSSI la phase. Ils ne "
          "sont donc pas redondants malgré leur parenté : à dispersion égale, minimiser "
          "aussi le DVH fait tourner φ vers 0° ou 180°, ce qui reporte l'effet résiduel "
          "de gravité sur les positions 3H et 9H plutôt que sur la position notée. C'est "
          "exactement le degré de liberté que la dispersion seule ignorait, et la raison "
          "pour laquelle les deux figurent dans le même objectif (§4.9)."),

    ("h2", "3.4 Isochronisme propre — période exacte de l'oscillateur conservatif"),
    ("p", "La formule d'Airy ne convient pas pour cette contribution : c'est une "
          "perturbation d'un oscillateur supposé harmonique, alors que l'effet recherché "
          "EST l'anharmonicité elle-même. On calcule donc la période exacte."),
    ("p", "L'énergie élastique emmagasinée pour une rotation α du balancier se déduit du "
          "champ de moment de la solution élastica :"),
    ("eq", ["U(α)  =  INTEGRALE[s = 0 → L]  M(s)² / [ 2 · E · I(s) ] · ds"]),
    ("p", "Le couple de rappel vaut C(α) = dU/dα ; si U était exactement quadratique, "
          "l'oscillateur serait harmonique et parfaitement isochrone. L'écart de U à la "
          "parabole est donc exactement le défaut d'isochronisme, et c'est cette courbe "
          "U(α) qui porte toute l'information."),
    ("p", "La conservation de l'énergie mécanique de l'oscillateur libre s'écrit :"),
    ("eq", ["½ · J · (dθ/dt)²  +  U(θ)  =  E",
            "dt  =  √(J/2) · dθ / √( E − U(θ) )"]),
    ("p", "Le spiral n'est pas symétrique en armage et en désarmage : une spire qui se "
          "resserre ne se raidit pas comme une spire qui s'ouvre, donc U(+θ) ≠ U(−θ). Les "
          "deux points de rebroussement d'une oscillation ne sont par conséquent pas "
          "symétriques en angle — ils sont à ÉNERGIE égale. On note a₊ et a₋ les deux "
          "amplitudes de rebroussement, solutions de U(a₊) = U(−a₋) = E, et la période "
          "complète vaut :",
     "Asymétrie. "),
    ("eq", ["T(E)  =  2 · [ T₊ + T₋ ]",
            "T±  =  √(J/2) · INTEGRALE[θ = 0 → a±]  dθ / √( E − U(±θ) )"]),
    ("p", "L'amplitude rapportée est la demi-course Θ₀ = (a₊ + a₋)/2, qui est la grandeur "
          "que mesure un chronocomparateur. Symétriser U aurait supprimé cette asymétrie, "
          "qui contribue pourtant réellement au défaut d'isochronisme (contribution au "
          "second ordre, vérifiée contre la solution analytique de Duffing, §4.3)."),
    ("p", "La marche correspondante est rapportée à une amplitude de référence — seule la "
          "variation avec l'amplitude a un sens physique, un décalage constant étant de "
          "toute façon rattrapé au réglage de la raquette :"),
    ("eq", ["marche_iso(Θ₀) = 86400 · [ T(Θ_ref) − T(Θ₀) ] / T(Θ_ref)",
            "positive = la montre avance à cette amplitude par rapport à la référence"]),

    ("h2", "3.5 Erreur d'échappement"),
    ("p", "L'échappement à ancre suisse perturbe l'oscillateur par un couple parasite de "
          "module quasi constant (frottement de repos et de tirage, impulsion, "
          "dégagement). Un couple parasite constant appliqué à un oscillateur d'amplitude "
          "Θ₀ produit, au premier ordre, une erreur de marche proportionnelle à 1/Θ₀ : "
          "l'effet est d'autant plus grand que l'amplitude est faible. C'est le retard aux "
          "faibles amplitudes."),
    ("p", "Le modèle ne décrit pas la géométrie de l'échappement. Il prend le retard "
          "comme une donnée — sa valeur à une amplitude basse de référence — et n'en "
          "modélise que la loi de variation :"),
    ("eq", ["E(Θ) = K/Θ − K/Θ_ref        avec  K = − retard / ( 1/Θ_basse − 1/Θ_ref )"]),
    ("p", "Deux nombres (le retard et l'amplitude à laquelle il est spécifié) calent "
          "entièrement le modèle. Réglage par défaut : 8 s/j de retard à 200°, référence "
          "à 300°."),
    ("p", "La cible d'optimisation en découle directement : pour que le total soit plat, "
          "le spiral doit produire de lui-même, par son isochronisme propre, une AVANCE "
          "de +8 s/j à 200°. L'écart à cette cible est la grandeur contrainte :"),
    ("eq", ["écart_compensation = marche_iso(Θ_basse) − retard_échappement"]),
]

_PARTIE_2 = [

    ("h2", "3.6 Courbes terminales — Grossmann (intérieure) et Phillips (extérieure)"),
    ("p", "La terminologie retenue suit l'usage historique : la courbe de Grossmann est "
          "la courbe terminale INTÉRIEURE (côté collet), la courbe de Phillips la courbe "
          "terminale EXTÉRIEURE (côté piton)."),
    ("p", "La théorie de Phillips (1860) établit une condition de concentricité du "
          "développement du spiral. Point essentiel pour ce projet : elle ne prescrit "
          "AUCUNE forme précise pour la courbe terminale, seulement une condition sur la "
          "position du centre de gravité de cette courbe — celui-ci doit se trouver à une "
          "distance d du centre du spiral, sur le rayon perpendiculaire à la droite "
          "joignant le centre au point de décrochement :"),
    ("eq", ["d = R² / L",
            "R = rayon au point de décrochement, L = longueur développée de la courbe"]),
    ("p", "Les formes classiques (arc unique, double arc, double arc plus segment droit) "
          "sont des simplifications héritées du pliage MANUEL d'un ressort métallique. En "
          "silicium usiné par DRIE, cette contrainte de fabrication disparaît "
          "entièrement."),
    ("p", "Le programme exploite cette liberté : chaque courbe terminale est une "
          "courbe de Bézier sur le rayon r(u), à nombre de points de contrôle "
          "configurable, dont les points intermédiaires sont des variables "
          "d'optimisation. Des inflexions sont donc possibles. La condition de Phillips "
          "est calculée et affichée comme repère de plausibilité, mais n'est PAS imposée "
          "comme contrainte : le vrai critère de validation est l'orbite du centre de "
          "gravité calculée par l'élastica complet sur les amplitudes et positions "
          "réelles.",
     "Choix. "),
    ("p", "Chaque zone terminale possède une extrémité ANCRÉE, raccordée en continuité C1 "
          "(position et pente imposées) au pas constant du corps du spiral, et une "
          "extrémité LIBRE, à rayon imposé (r_min au collet, r_max au piton) mais de "
          "pente libre."),
    ("p", "L'asymétrie entre les deux courbes est voulue et physiquement cohérente. Côté "
          "piton, r_max est une contrainte d'encombrement dure, indépendante de la ligne "
          "de référence à pas constant : c'est le désaccord entre « où la ligne de "
          "référence arriverait » et « où il faut vraiment arriver » que la courbe de "
          "Phillips doit absorber — sa raison d'être historique. Côté collet, r_min est "
          "par construction sur la ligne de référence : la courbe de Grossmann n'a aucun "
          "désaccord géométrique à rattraper, ses points libres ne pilotent que la forme."),

    ("h2", "3.7 Terminaison à raquette — arc concentrique et courbe en S"),
    ("p", "Un spiral destiné à recevoir une raquette ne peut pas se terminer par une "
          "courbe de Phillips. Les goupilles de la raquette pivotent autour de l'axe du "
          "balancier et glissent le long du spiral : pour que leur déplacement ne change "
          "QUE la longueur active, la portion qu'elles balaient doit être un ARC CENTRÉ "
          "SUR CET AXE. Sur toute autre forme, déplacer la raquette pousserait aussi le "
          "spiral radialement."),
    ("p", "La terminaison devient donc : corps d'Archimède, puis courbe en S, puis arc "
          "concentrique au rayon r_max, puis piton. La courbe en S doit être tangente au "
          "corps d'un côté (continuité C1 avec le pas d'Archimède) et à l'arc de l'autre "
          "(pente NULLE, l'arc étant à rayon constant). Quatre conditions déterminent "
          "exactement une cubique ; au-delà de quatre points de contrôle, les points "
          "supplémentaires redeviennent des variables d'optimisation sans toucher aux "
          "tangences."),
    ("p", "L'arc concentrique est décrit par DEUX angles indépendants, comptés de part et "
          "d'autre du point de comptage : l'arc AVANT (entre la fin de la courbe en S et "
          "le point de comptage) et l'arc APRÈS (entre le point de comptage et le piton). "
          "Les prendre égaux place le point de comptage au milieu de la portion "
          "circulaire terminale, ce qui laisse à la raquette la même course dans les deux "
          "sens ; les dissocier permet de loger le piton plus loin sans déplacer le point "
          "de comptage.",
     "Découpage de l'arc. "),
    ("p", "La forme en S n'est pas imposée, elle apparaît d'elle-même : la pente "
          "moyenne exigée pour gagner le rayon manquant dépasse toujours largement le pas "
          "d'Archimède (mesuré : 117 à 700 µm/rad selon l'angle de transition, contre "
          "25 µm/rad de pas), donc la courbe doit se redresser puis s'aplatir.",
     "Pourquoi un S. "),
    ("p", "La longueur de l'arc n'est pas un réglage : elle se DÉDUIT de la plage de "
          "réglage exigée. La raideur variant comme l'inverse de la souplesse intégrée :"),
    ("eq", ["plage (s/j) = 86400 · ½ · ΔS / S        avec  S = INTEGRALE( ds / E·I )",
            "L_arc = 4 · plage · S · (E·I)_piton / 86400        θ_arc = L_arc / r_max"]),
    ("p", "La sensibilité varie comme le cube de l'épaisseur au piton : un piton aminci "
          "rend la raquette beaucoup plus sensible et raccourcit d'autant l'arc "
          "nécessaire. La plage obtenue est ensuite VÉRIFIÉE sur la géométrie construite "
          "(mesuré : ±30 et ±60 s/j demandés, obtenus à mieux que 0,1 %), et non supposée "
          "acquise du fait du dimensionnement analytique."),

    ("h2", "3.8 Point de comptage, nombre entier de spires et force sur le pivot"),
    ("p", "Sur un spiral à raquette, le spiral ACTIF s'arrête aux goupilles, pas au "
          "piton : la portion au-delà est tenue à ses deux extrémités, elle ne travaille "
          "pas et ne se déforme pas (vérifié : déplacement nul à la précision machine "
          "entre 0 et 300° de balancier). Elle occupe pourtant de la place et doit être "
          "usinée, donc elle compte pour le jeu entre spires et pour l'export."),
    ("p", "Le nombre de tours réels et fictifs compte du collet au POINT DE COMPTAGE, "
          "jamais jusqu'au bout de la matière. Un nombre entier place ainsi le point de "
          "comptage sur la ligne des centres, comme l'attache à la virole, conformément à "
          "la règle horlogère classique. Ce n'est toutefois pas une contrainte du "
          "programme : le nombre de tours peut être fractionnaire, et le point de comptage "
          "tombe alors où il tombe."),
    ("p", "Trois hypothèses ont été testées sur ce que cette règle apporte ; deux "
          "sont fausses.",
     "Ce que la règle apporte réellement. "),
    ("p", "Le décentrage statique varie de 52,7 à 58,3 µm quand on parcourt un "
          "tour entier de partie fractionnaire, sans minimum aux tours entiers. Le "
          "décentrage d'une spirale d'Archimède vaut pas/π indépendamment du nombre de "
          "tours (calcul : 50,0 µm, mesure : 56,4 µm) : chaque spire complète contribue "
          "dans le même sens.",
     "Ce n'est PAS une question de poise. "),
    ("p", "À géométrie identique, tourner le spiral "
          "dans le boîtier fait varier la dispersion de 1,774 à 2,135 s/j avec une "
          "période de exactement 90° — celle du jeu de positions 3H/6H/9H/12H. C'est un "
          "effet d'échantillonnage, pas une propriété du spiral.",
     "Ce n'est PAS l'orientation vis-à-vis des positions. "),
    ("p", "Le spiral applique au pivot du balancier une "
          "résultante atteignant 99,4 µN, soit treize fois son propre poids, dont l'axe "
          "bascule de 180° à chaque alternance. Mesuré : cet axe suit la droite collet → "
          "point de comptage à quelques degrés près. Un nombre entier de spires aligne "
          "donc cette force sur la ligne des centres, ce qui fixe son orientation par "
          "rapport à la réaction de l'ancre.",
     "C'est l'orientation de la FORCE SUR LE PIVOT. "),
    ("p", "La règle reste justifiée, non par ce que ce "
          "programme optimise, mais par le frottement du pivot, qui n'est pas modélisé. "
          "Elle ne coûte rien — c'est un choix d'entier au moment de définir le calibre — "
          "et protège d'un effet que le calcul ne voit pas. Le programme affiche le module "
          "et l'axe de cette force pour que l'alignement soit vérifiable d'un coup d'œil. "
          "L'axe est rendu MODULO 180°, puisque c'est une droite et non un vecteur : le "
          "comparer à un angle rendu modulo 360° conclurait à un demi-tour inexistant.",
     "Conséquence pour la conception. "),

    ("h2", "3.9 Jeu entre les goupilles de la raquette"),
    ("p", "Tant que la lame est plaquée contre une goupille, le spiral est encastré AUX "
          "GOUPILLES : longueur active jusqu'à elles, raideur k_goupille. Dès qu'elle "
          "décolle et traverse le jeu, plus rien ne la tient là : le spiral se comporte "
          "comme encastré AU PITON, donc plus long et plus souple, k_libre. La montre "
          "retarde pendant ces traversées, et d'autant plus que l'amplitude est faible — "
          "les traversées occupent alors une part plus grande de l'oscillation."),
    ("p", "La raideur est donc constante par morceaux en fonction de l'angle du "
          "balancier : k_libre dans la bande morte |φ| < φ_c, k_goupille au-delà. La "
          "demi-largeur de bande vaut :"),
    ("eq", ["φ_c = (jeu / 2) / c",
            "c = excursion latérale de la lame aux goupilles, en µm/rad"]),
    ("p", "c se calcule par le solveur elastica dans l'état libre, en mesurant l'écart de "
          "la lame déformée au point FIXE des goupilles — un point matériel donnerait une "
          "autre grandeur, sans rapport avec le contact, puisque la lame coulisse entre "
          "les goupilles. Le jeu est évalué sur l'épaisseur FINIE de la lame, oxyde "
          "compris, et non sur la cote de gravure."),
    ("p", "Deux dépendances, toutes deux vérifiées numériquement, et "
          "toutes deux conformes à l'expérience horlogère :",
     "Ce qui commande l'effet. "),
    ("tab",
     ["Longueur après les goupilles", "Excursion latérale sur ±340°",
      "Contact pour 20 µm de jeu", "Saut de marche si le maintien est perdu"],
     [["124 µm", "0,9 µm", "jamais atteint", "50 s/j"],
      ["371 µm", "7,6 µm", "jamais atteint", "150 s/j"],
      ["742 µm", "30 µm", "227° — dans la plage", "301 s/j"],
      ["1856 µm", "175 µm", "39°", "758 s/j"]],
     [26, 26, 24, 24]),
    ("p", "L'excursion varie comme le CARRÉ de la longueur après les goupilles. Trop "
          "courte, la lame flotte en permanence et la raquette n'a aucune autorité ; trop "
          "longue, le retard devient important. Le pire cas est intermédiaire, là où la "
          "bande morte est comparable à l'amplitude de travail."),

    ("h2", "3.10 Compensation thermique par oxyde de silicium"),
    ("p", "Le silicium s'assouplit fortement en chauffant : son coefficient "
          "thermoélastique vaut environ −60 ppm/K. Un spiral nu retarde donc de "
          "3,85 s/j/°C, soit 115 s/j entre 8 et 38 °C — six fois le critère COSC. L'oxyde "
          "thermique a le comportement inverse (+215 ppm/K)."),
    ("p", "En flexion, la contribution d'une couche varie comme le cube "
          "de sa distance à la fibre neutre. L'oxyde étant à la périphérie de la section, "
          "il pèse très au-dessus de sa part volumique : 3,7 µm d'oxyde sur une lame de "
          "31 µm portent déjà un tiers de la rigidité. La section oxydée étant le "
          "rectangle extérieur moins le cœur de silicium :",
     "Pourquoi si peu suffit. "),
    ("eq", ["E·I = E_oxyde · ( I_ext − I_coeur )  +  E_silicium · I_coeur"]),
    ("p", "L'oxydation thermique consomme 0,44 µm de silicium par micron d'oxyde crû et "
          "déborde de 0,56 µm vers l'extérieur. La cote de GRAVURE et la cote FINIE "
          "diffèrent donc de 2 × 0,56 × t, soit typiquement 4 µm — un tiers du jeu minimal "
          "exigé entre spires. Le programme exporte les deux fichiers DXF séparément."),
    ("p", "L'épaisseur d'oxyde est RECALCULÉE pour chaque "
          "candidat, et c'est nécessaire : la part de rigidité portée par l'oxyde passe "
          "de 48 % sur une section de 20 µm à 22 % sur une section de 50 µm, or "
          "l'optimiseur fait varier l'épaisseur de ±35 %. Épaisseur d'oxyde et mise à "
          "l'échelle en fréquence sont de plus couplées (l'oxyde change la raideur, "
          "l'échelle change le rapport oxyde/silicium) : les deux sont résolues par point "
          "fixe.",
     "Couplage avec le profil d'épaisseur. "),
    ("p", "Modeler le profil d'épaisseur n'apporte AUCUN levier "
          "supplémentaire sur la courbure thermique. Le coefficient global reste "
          "⟨β⟩·α_oxyde + (1−⟨β⟩)·α_silicium, une famille à un seul paramètre — exactement "
          "ce que règle déjà l'épaisseur d'oxyde. On ne peut pas annuler à la fois la "
          "pente et la courbure par la géométrie seule.",
     "Un résultat négatif utile. "),
]

_PARTIE_3 = [

    # ================================================================ 4
    ("saut",),
    ("h1", "4. Outils mathématiques"),
    ("h2", "4.1 Résolution du problème de tir"),
    ("p", "Le problème aux limites du §3.1 est résolu par méthode de tir : on se donne "
          "(Fx, Fy, M₀) au piton, on intègre la poutre jusqu'au collet, et on ajuste ces "
          "trois paramètres pour que les trois conditions au collet soient satisfaites."),
    ("p", "Le système étant exactement déterminé (3 équations, 3 inconnues), il est résolu "
          "comme un système d'équations non linéaires par la méthode hybride de Powell "
          "(scipy.optimize.root, méthode « hybr »), et non par minimisation d'un coût "
          "agrégé."),
    ("p", "Ce point n'est pas cosmétique. La version antérieure "
          "minimisait par BFGS le coût dx² + dy² + (0,1·dφ)², dans lequel la condition "
          "angulaire était sous-pondérée d'un facteur 10 : la minimisation s'arrêtait sur "
          "des résidus atteignant 0,24° d'angle au collet. L'effet sur l'orbite du centre "
          "de gravité était invisible (grandeur du premier ordre), mais l'énergie "
          "élastique s'en trouvait faussée jusqu'à 2,4·10⁻³ en relatif — soit environ dix "
          "fois l'effet d'isochronisme que l'on cherche justement à mesurer dessus "
          "(8 s/j correspond à ≈ 2·10⁻⁴ en relatif sur U).",
     "Justification mesurée. "),
    ("tab",
     ["Méthode", "Résidu de position", "Résidu d'angle", "Erreur relative sur U", "Coût"],
     [["Minimisation BFGS du coût agrégé", "~0,1 µm", "jusqu'à 0,24°", "jusqu'à 2,4·10⁻³",
       "2,1 s / 35 angles"],
      ["Résolution du système (Powell hybride)", "~10⁻⁶ µm", "~10⁻⁹ °", "référence",
       "0,51 s / 35 angles"]],
     [30, 17, 15, 20, 18]),
    ("p", "Le résidu de la résolution directe -- 10⁻⁶ µm, soit un picomètre sur une "
          "pièce de 2,7 mm -- est celui de l'arithmétique flottante et non de la "
          "convergence : les coordonnées passent par les sinus et cosinus de "
          "l'orientation dans le mouvement (§3.3), qui coûtent quelques unités dans le "
          "dernier chiffre."),
    ("p", "La résolution directe est donc à la fois beaucoup plus précise et quatre fois "
          "plus rapide. Un garde-fou anti-divergence subsiste : si le résidu de position "
          "dépasse 1 µm, deux stratégies de repli sont essayées (même solveur reparti d'un "
          "autre point, puis minimisation du coût agrégé, qui ne peut pas diverger) et la "
          "meilleure est retenue."),

    ("h2", "4.2 Intégration de la poutre — schéma prédicteur-correcteur"),
    ("p", "L'équilibre du §3.1 est implicite : le moment en un point dépend de la position "
          "DÉFORMÉE de ce point, qui dépend elle-même du moment. L'intégration résout ce "
          "couplage en une seule passe, du piton vers le collet, en menant conjointement "
          "l'intégration de l'angle et celle de la position :"),
    ("eq", ["pour i de n-2 à 0 :",
            "    M_pred = M₀ + Fx·(Y[i+1] − y_p) − Fy·(X[i+1] − x_p)   # prédicteur",
            "    φ_pred = φ[i+1] − dφ₀ − (M_pred / EI[i])·ds",
            "    X_pred = X[i+1] − ds·cos(φ_pred) ;  Y_pred = Y[i+1] − ds·sin(φ_pred)",
            "    M      = M₀ + Fx·(Y_pred − y_p) − Fy·(X_pred − x_p)   # correcteur",
            "    φ[i]   = φ[i+1] − dφ₀ − (M / EI[i])·ds",
            "    X[i]   = X[i+1] − ds·cos(φ[i]) ;  Y[i] = Y[i+1] − ds·sin(φ[i])"]),
    ("p", "Le prédicteur explicite estime la position du pas courant, le correcteur y "
          "recalcule le moment. Le résultat a été vérifié contre une résolution par point "
          "fixe menée jusqu'à convergence (10 à 11 itérations de la boucle complète) : "
          "l'écart relatif résiduel sur l'énergie est de l'ordre de 10⁻¹², pour un coût de "
          "calcul inchangé."),
    ("p", "La boucle est compilée par Numba (njit) ; elle est appelée quelques dizaines de "
          "fois par résolution de tir, elle-même appelée à chaque angle de balayage, pour "
          "chaque candidat de chaque génération de l'optimiseur. C'est le seul point "
          "réellement coûteux du programme."),

    ("h2", "4.3 Quadrature de l'intégrale de période et régularisation"),
    ("p", "L'intégrande de la période diverge en 1/√ au point de rebroussement (divergence "
          "intégrable). Le changement de variable θ = a·sin u, avec u parcourant [0, π/2], "
          "la régularise complètement : le cos u du jacobien s'annule exactement à la même "
          "vitesse que le √(E − U) du dénominateur."),
    ("eq", ["T±  =  √(J/2) · INTEGRALE[u = 0 → π/2]  a·cos u / √( E − U(a·sin u) ) · du",
            "limite quand u → π/2   =   √( 2a / U'(a) )"]),
    ("p", "Pour un ressort LINÉAIRE, l'intégrande devient exactement "
          "constant, égal à √(2/k) : la quadrature est alors exacte, indépendamment du "
          "nombre de points. C'est une propriété essentielle ici, puisqu'on ne mesure "
          "qu'un ÉCART à l'isochronisme de l'ordre de 10⁻⁴ en relatif — une quadrature qui "
          "ne serait pas exacte sur le cas linéaire noierait le signal utile dans son "
          "propre biais.",
     "Propriété recherchée. "),
    ("p", "L'intégration numérique utilise la règle de Simpson sur 401 points. Les points "
          "de rebroussement conjugués sont appariés par RÉSOLUTION (racines exactes du "
          "polynôme par morceaux de la spline, via CubicSpline.solve, et recherche de zéro "
          "par la méthode de Brent pour l'amplitude cible), et non par interpolation d'une "
          "table : une table interpolée introduisait un plancher de bruit de 0,2 à "
          "0,4 s/j, soit 3 à 5 % d'une cible de 8 s/j."),
    ("p", "Le calcul a été validé sur trois cas de référence analytiques :",
     "Validation. "),
    ("tab",
     ["Cas de référence", "Valeur attendue", "Valeur calculée", "Accord"],
     [["Ressort linéaire (isochrone exact par construction)", "0 s/j", "3,6·10⁻⁹ s/j",
       "précision machine"],
      ["Ressort durcissant U = ½kθ² + ¼bθ⁴, b/k = 10⁻⁴ (Duffing, 1er ordre : "
       "δω/ω = 3bA²/8k)", "−49,348 s/j", "−49,283 s/j", "0,13 %"],
      ["Asymétrie U = ½kθ² + (a₃/3)θ³, a₃/k = 10⁻² (Duffing quadratique, 2e ordre)",
       "54,83 s/j", "55,00 s/j", "3 chiffres significatifs"],
      ["Sensibilité au pas du balayage angulaire (30° → 5°)", "—",
       "0,005 s/j d'écart", "converge"]],
     [42, 17, 17, 24]),

    ("h2", "4.4 Représentation des géométries"),
    ("p", "Courbes de Bézier sur le rayon r(u), en base de "
          "Bernstein :",
     "Courbes terminales. "),
    ("eq", ["r(u)  =  SOMME[i = 0 → n]  C(n,i) · (1−u)^(n−i) · u^i · P_i",
            "u dans [0,1],   n = n_points − 1"]),
    ("p", "Sur n_points points de contrôle, trois sont fixés par les contraintes dures "
          "(rayon imposé à chaque extrémité) et par la continuité C1 à l'extrémité "
          "ancrée ; les n_points − 3 autres sont libres et constituent les variables de "
          "forme. Ils sont paramétrés en ÉCART par rapport à un raccord rectiligne, de "
          "sorte que la valeur nulle corresponde à la forme la plus neutre possible. La "
          "courbe en S de la terminaison à raquette obéit aux mêmes règles, avec une "
          "condition de tangence de plus (pente nulle à l'arrivée sur l'arc), donc un "
          "point de contrôle fixé de plus."),
    ("p", "La monotonie de r(u) est vérifiée via une propriété des polynômes de "
          "Bernstein : si les points de contrôle sont monotones non décroissants, la "
          "courbe l'est aussi (condition suffisante). Un spiral non monotone replierait "
          "physiquement le ruban sur lui-même ; les candidats concernés sont rejetés "
          "d'emblée par l'optimiseur."),
    ("p", "Interpolation PCHIP "
          "(Piecewise Cubic Hermite Interpolating Polynomial) sur un petit nombre de "
          "nœuds de contrôle. Le choix de PCHIP plutôt qu'une spline cubique classique est "
          "délibéré : PCHIP préserve la monotonie locale et ne produit ni dépassement ni "
          "oscillation entre les nœuds, ce qui garantit qu'une épaisseur reste positive et "
          "physiquement plausible quelles que soient les valeurs tirées par l'optimiseur.",
     "Profils d'épaisseur e(s) et d'écart radial δr(s). "),
    ("p", "Splines cubiques (CubicSpline), pour deux "
          "usages distincts : l'orbite du centre de gravité, dont la formule d'Airy exige "
          "la DÉRIVÉE par rapport à l'angle, et la courbe d'énergie U(θ), dont le calcul "
          "de période exige les racines exactes. Le choix de la spline cubique est ici "
          "justifié par le fait qu'elle reproduit exactement les polynômes de degré ≤ 3 : "
          "l'erreur d'interpolation sur la partie quasi quadratique de U — la partie "
          "dominante — est donc nulle, et ne subsiste que sur l'écart, lui-même petit.",
     "Interpolation des grandeurs balayées. "),

    ("p", "Le maillage de calcul et la "
          "description du contour dans un fichier de fabrication répondent à deux "
          "exigences sans rapport. Le premier est un compromis entre précision et temps, "
          "et il densifie VOLONTAIREMENT les zones terminales, là où la courbure varie "
          "vite ; le second doit seulement décrire la pièce assez finement pour la "
          "machine, uniformément. Le programme dissocie donc les deux : l'export DXF "
          "accepte un nombre de points par flanc indépendant du maillage.",
     "Ré-échantillonnage pour l'export. "),
    ("p", "Le ré-échantillonnage interpole à ABSCISSE CURVILIGNE uniforme, et deux "
          "précautions le rendent sûr. La partie active et la portion morte sont traitées "
          "séparément, pour que le point de comptage — frontière physique, là où les "
          "goupilles saisissent la lame — reste exactement un point du maillage au lieu "
          "de tomber entre deux. Et l'épaisseur est interpolée par PCHIP et non par "
          "spline cubique : une spline peut dépasser entre deux nœuds et produire "
          "localement une cote hors de la plage d'origine, ce qui n'a aucun sens sur une "
          "cote de fabrication."),
    ("p", "La longueur développée diminue légèrement quand on enlève des points — une "
          "polyligne à cordes coupe les virages. Mesuré : 0,09 % d'écart à 600 points, "
          "0,006 % à 5 000, sur un maillage de calcul de 1 670. C'est le seul effet du "
          "ré-échantillonnage sur la géométrie, et il décroît comme attendu."),

    ("h2", "4.5 Mise à l'échelle analytique de l'épaisseur"),
    ("p", "La fréquence propre cible est atteinte exactement, sans itération, par une mise "
          "à l'échelle homothétique du profil d'épaisseur. La raideur d'une poutre de "
          "section rectangulaire vaut C = E·h / (12 · Σ ds/e³) ; multiplier tout le profil "
          "par β multiplie donc la raideur par β³ et la fréquence par β^(3/2). D'où :"),
    ("eq", ["C  =  E · h / ( 12 · SOMME[ ds / e³ ] )",
            "β  =  ( f_cible / f_obtenue )^(2/3)          puis          e(s) ← β · e(s)"]),
    ("p", "Cette correction est appliquée systématiquement à la construction de toute "
          "géométrie paramétrique. Elle libère l'optimiseur de la contrainte de "
          "fréquence : il n'explore que des FORMES, la fréquence étant garantie par "
          "construction. Mesuré : 3,5000000 Hz dans toutes les configurations testées, "
          "oxyde imposé et terminaison à raquette compris."),

    ("h2", "4.6 Jeu inter-spires"),
    ("p", "Le jeu entre spires adjacentes ne peut pas être calculé en supposant un pas "
          "constant, puisque le pas est précisément une variable d'optimisation. Il est "
          "obtenu en recherchant, pour chaque point de la fibre neutre déformée, le point "
          "situé un tour plus loin en angle polaire déroulé, par interpolation sur la "
          "géométrie réellement déformée :"),
    ("eq", ["jeu(i)  =  r(θ_i + 2π) − r(θ_i)  −  [ e(θ_i) + e(θ_i + 2π) ] / 2"]),
    ("p", "Le jeu retenu est le minimum sur TOUS les angles du balayage, et non seulement "
          "aux amplitudes extrêmes : l'élastica n'est pas monotone en angle au voisinage "
          "du minimum de jeu, et une vérification aux seuls extrêmes laissait passer des "
          "collisions réelles. L'épaisseur employée est la cote FINIE, oxyde compris."),

    ("h2", "4.7 Fonction descriptive d'une zone morte (jeu de raquette)"),
    ("p", "Le jeu de goupilles introduit une non-linéarité PAR MORCEAUX : la raideur "
          "change brutalement quand la lame quitte ou retrouve l'appui. Ni la formule "
          "d'Airy (perturbation d'un oscillateur supposé harmonique) ni l'intégration "
          "exacte de la période (qui exige une énergie régulière et monotone) ne "
          "s'appliquent. L'outil adapté est l'équilibrage harmonique : on cherche la "
          "raideur qu'une oscillation sinusoïdale d'amplitude Θ₀ ressent en moyenne."),
    ("eq", ["k_eq(Θ₀) = k_libre + ( k_goupille − k_libre ) · N(ε)     ε = min( φ_c/Θ₀ , 1 )",
            "N(ε) = 1 − (2/π) · [ arcsin(ε) + ε·√(1−ε²) ]",
            "fonction descriptive classique d'une zone morte"]),
    ("p", "Deux cas limites donnent tout le sens de la formule :"),
    ("puce", "ε = 0 (jeu nul) : N = 1, donc k_eq = k_goupille EXACTEMENT, à toute "
             "amplitude, donc une marche identiquement nulle. C'est rigoureusement la "
             "configuration à piton fixe. Le spiral libre et le spiral à raquette ne sont "
             "donc pas deux modèles à réconcilier : ce sont les deux bouts d'un seul "
             "paramètre continu."),
    ("puce", "ε = 1 (amplitude égale à la bande morte) : N = 0, k_eq = k_libre. La lame ne "
             "touche jamais, la raquette n'a aucune autorité."),
    ("p", "Conséquence de conception : un même spiral ne peut pas être parfait dans les "
          "deux montages, et l'écart entre eux est une PENTE en fonction de l'amplitude, "
          "pas un décalage constant — la raquette, qui ne corrige qu'un décalage, ne peut "
          "donc pas le rattraper. Le programme laisse choisir le montage visé (piton fixe, "
          "raquette, ou mi-chemin) et affiche l'erreur résiduelle dans chacun."),

    ("h2", "4.8 Section composite silicium/oxyde"),
    ("p", "La rigidité de la lame oxydée se calcule comme rectangle extérieur moins cœur "
          "de silicium, ce qui est exact pour une couche d'épaisseur uniforme sur une "
          "section rectangulaire. Elle se réduit EXACTEMENT à E_silicium·I quand "
          "l'épaisseur d'oxyde est nulle (vérifié à la précision machine) : la lame nue "
          "n'est pas un cas particulier traité à part."),
    ("p", "Toute la physique utilise ensuite la rigidité du composite et l'encombrement de "
          "la lame FINIE : le solveur elastica, la fréquence propre, la masse (l'oxyde est "
          "un peu moins dense que le silicium qu'il remplace) et surtout le contrôle de "
          "jeu entre spires. L'épaisseur d'oxyde qui annule la variation de marche est "
          "cherchée par la méthode de Brent, bornée par l'épaisseur qui consommerait "
          "entièrement la section la plus mince."),
    ("p", "Ne pas coupler le calcul et ajouter l'oxyde APRÈS "
          "l'optimisation introduirait +1449 s/j d'erreur de fréquence (mesuré) : la lame "
          "gravée doit être amincie pour compenser le raidissement apporté par l'oxyde.",
     "Ordre de grandeur de l'enjeu. "),

    ("h2", "4.9 Optimisation multi-objectifs — NSGA-II"),
    ("p", "Le vecteur de conception est de longueur fixe, dans l'ordre suivant :"),
    ("eq", ["[ epaisseur_noeuds (Me) , delta_r_noeuds (Mr) ,",
            "  grossmann_angle_deg , grossmann_offsets (Ng) ,",
            "  angle_exterieur_deg , offsets_exterieurs (Np) ]"]),
    ("p", "La zone terminale extérieure est décrite par les mêmes cases du vecteur, qu'il "
          "s'agisse d'une courbe de Phillips (spiral libre) ou de la courbe en S d'une "
          "terminaison à raquette : le mode est figé pour la durée d'une optimisation, le "
          "vecteur garde donc bien une longueur constante, comme l'exige pymoo."),
    ("p", "Me, Mr, Ng, Np sont des choix de DISCRÉTISATION DE LA CONCEPTION (nombre de "
          "degrés de liberté de forme), figés pour la durée d'une optimisation. La densité "
          "de maillage physique et le degré des Bézier n'en font pas partie, pour la même "
          "raison."),
    ("p", "L'algorithme est NSGA-II (Non-dominated Sorting Genetic Algorithm II, "
          "implémentation pymoo), choisi pour un espace de recherche de grande dimension, "
          "non convexe, sans gradient disponible et parsemé de candidats invalides."),
    ("tab",
     ["Terme", "Expression", "Nature"],
     [["F₀", "dispersion de marche entre positions verticales  +  |DVH|  (s/j)",
       "objectif, à minimiser"],
      ["F₁", "− min( jeu inter-spires , plafond )", "objectif, à minimiser"],
      ["G₀", "jeu_min_requis − jeu mesuré ≤ 0", "contrainte"],
      ["G₁", "|écart de compensation| − tolérance ≤ 0", "contrainte"],
      ["G₂", "|oxyde requis − centre de la fenêtre| − demi-largeur ≤ 0",
       "contrainte, facultative"]],
     [10, 60, 30]),
    ("p", "Les deux se lisent en s/j et vont dans le même "
          "sens : aucun arbitrage à exposer, donc rien à gagner à diluer le front de "
          "Pareto d'une dimension de plus. Ils ne sont pourtant pas redondants (§3.3) — la "
          "dispersion ne voit que le module A du vecteur de sensibilité, le DVH voit aussi "
          "sa phase. À dispersion égale, la somme préfère donc les géométries dont la "
          "phase reporte l'effet résiduel sur 3H et 9H plutôt que sur la position 6H, "
          "celle que le COSC note. Le DVH est pris au PIRE CAS sur la plage d'amplitude, "
          "par cohérence avec la dispersion. Le poids reste exposé dans le code, à zéro "
          "pour retrouver l'objectif antérieur.",
     "Pourquoi le DVH s'AJOUTE à la dispersion au lieu de former un 3e objectif. "),
    ("p", "Ce n'est pas une "
          "grandeur à minimiser mais une cible à atteindre : au-delà de la tolérance, un "
          "spiral qui surcompense n'est pas « meilleur », il est faux dans l'autre sens. "
          "En faire un objectif aurait porté le front de Pareto à trois dimensions, le "
          "diluant fortement à population constante, pour un arbitrage qui n'a pas lieu "
          "d'être. NSGA-II minimisant d'abord la violation de contrainte, la population "
          "converge vers la cible avant d'arbitrer dispersion contre jeu — le comportement "
          "recherché. Si la cible se révèle inatteignable, l'algorithme dégrade proprement "
          "en s'en approchant au mieux au lieu de tout rejeter.",
     "Pourquoi la compensation d'échappement est une CONTRAINTE et non un objectif. "),
    ("p", "Le plafonnement du jeu dans F₁ évite une incitation à agrandir indéfiniment un "
          "jeu déjà suffisant au détriment de l'objectif chronométrique. Le décentrage "
          "moyen du centre de gravité est calculé mais délibérément non ajouté comme "
          "objectif : sur un spiral correctement réglé en Grossmann et Phillips, il est "
          "très fortement corrélé à la dispersion de marche — les deux découlent de la "
          "même orbite — et l'ajouter diluerait le front sans apporter d'arbitrage "
          "indépendant."),
    ("p", "La contrainte G₂ n'est active que si une fenêtre d'épaisseur d'oxyde est "
          "donnée : elle n'a de sens que lorsqu'une même oxydation doit servir toute une "
          "série. Sans fenêtre, l'épaisseur requise se lit simplement sur le front de "
          "Pareto, candidat par candidat."),
    ("p", "Un sondage aléatoire préalable classe des tirages "
          "uniformes dans les bornes et fournit le point de départ du warm start. La "
          "géométrie de départ elle-même y est incluse : les bornes étant centrées sur "
          "elle, un tirage uniforme s'en écarte forcément dans toutes les directions à la "
          "fois, et il est courant qu'aucun tirage ne fasse mieux que le point central. "
          "Mesuré sur un cas réel : la géométrie de départ offrait 21,7 µm de jeu, le "
          "meilleur de 60 tirages aléatoires −224 µm (spires croisées).",
     "Sondage et warm start. "),
]

_PARTIE_4 = [

    ("h2", "4.10 Reconstitution d'un contour importé"),
    ("p", "Un DXF ne livre pas un spiral : il livre un tas d'entités géométriques, dans "
          "un ordre quelconque, parfois mêlées à des accessoires de dessin. Trois étapes "
          "l'en extraient, et chacune a son piège, tous les trois rencontrés sur des "
          "fichiers réels."),
    ("p", "Les extrémités sont indexées dans un ARBRE K-D construit une "
          "fois, et chaque rattachement est une interrogation de cet arbre. Une recherche "
          "par balayage du pool restant coûte un temps CUBIQUE en nombre d'entités : "
          "invisible sur un dessin fait de quelques polylignes, rédhibitoire sur un export "
          "CAO en segments élémentaires. Mesuré sur un fichier de 6 087 LINE : plusieurs "
          "minutes contre 0,9 s.",
     "Chaînage des entités. "),
    ("p", "Une tolérance fixe ne convient pas. Trop serrée, "
          "elle laisse le contour en morceaux ; trop lâche, elle recolle au spiral des "
          "entités qui n'en font pas partie — une ligne de construction partant du centre "
          "ramenait ainsi le rayon de collet reconstruit de 0,468 à 0,249 mm, sans aucun "
          "signalement. Le programme essaie donc du plus serré au plus lâche et retient la "
          "première tolérance qui produise deux flancs de longueurs comparables, critère "
          "physique : les deux flancs d'un ruban ne diffèrent que de l'épaisseur intégrée "
          "sur la rotation totale, soit quelques pour cent. Si aucune ne convient, l'import "
          "échoue en détaillant ce qu'il a trouvé, plutôt que de livrer une géométrie fausse.",
     "Tolérance de recollement. "),
    ("p", "Le contour parcourt un flanc, "
          "contourne le bout du ruban, et revient par l'autre. Les deux plis se cherchent "
          "aux extrema de l'ANGLE POLAIRE DÉROULÉ, et non du rayon : θ croît strictement le "
          "long d'un flanc et décroît le long de l'autre, son maximum est donc au bout par "
          "construction. Le rayon, lui, est maximal sur le flanc EXTÉRIEUR, plus loin du "
          "centre : en contournant la pointe il continue de croître. Sur un spiral dont la "
          "courbe terminale rentre légèrement, le maximum de rayon tombait 102° APRÈS la "
          "pointe, la découpe était décalée d'autant, et l'intersection angulaire des deux "
          "flancs amputait le spiral de 102° au piton — silencieusement.",
     "Séparation des deux flancs. "),
    ("p", "L'appariement interpole r(θ) linéairement entre "
          "points consécutifs, or r(θ) n'est pas linéaire le long d'une corde : il passe "
          "par un minimum au pied de la perpendiculaire issue du centre. Là où le dessin "
          "est grossier — une portion rectiligne de courbe terminale dessinée en deux "
          "cordes de 948 et 593 µm quand le reste est échantillonné tous les 41 µm — "
          "l'épaisseur, différence de deux rayons interpolés, devient la différence de "
          "deux erreurs : elle tombait à 10 µm au lieu de 36. Les flancs sont donc "
          "densifiés LE LONG DE LEURS PROPRES SEGMENTS avant appariement. Des points posés "
          "sur une droite ne changent rien à la géométrie, et la flèche d'une corde de "
          "50 µm à r = 2,4 mm vaut 0,13 µm.",
     "Densification préalable. "),
    ("p", "Effet cumulé de ces trois corrections sur la boucle fermée export → import, "
          "mesurée sur une géométrie dont on connaît la vérité : l'erreur de rayon au "
          "collet passe de 11,80 à 0,04 µm. Les 12 µm étaient documentés depuis des "
          "semaines comme une « limite connue du pipeline » ; ils n'étaient que la "
          "tolérance de recollement."),

    # ================================================================ 5
    ("saut",),
    ("h1", "5. Choix de modélisation et leur justification"),
    ("p", "Le tableau ci-dessous récapitule les décisions structurantes. Chaque ligne "
          "indique l'alternative écartée et la mesure qui a motivé le choix."),
    ("tab",
     ["Choix retenu", "Alternative écartée", "Justification mesurée"],
     [["Bras de levier du moment pris sur la configuration DÉFORMÉE",
       "Bras de levier sur la configuration au repos (petits déplacements)",
       "L'erreur croît avec l'amplitude : 7·10⁻⁵ sur l'énergie à 40°, mais 2,2·10⁻² à "
       "340°, avec un déplacement du centre de gravité atteignant 183 µm. Une erreur "
       "croissante avec l'amplitude se lit exactement comme un défaut d'isochronisme : "
       "elle donnait −646 s/j à 200°, contre +19,3 s/j avec l'équilibre correct."],
      ["Résolution du système de tir (Powell hybride)",
       "Minimisation BFGS d'un coût agrégé",
       "Résidu angulaire ramené de 0,24° à 10⁻⁹ °, erreur sur l'énergie de 2,4·10⁻³ à la "
       "précision machine, et calcul 4× plus rapide. Le résidu résiduel (10⁻⁶ µm, soit un "
       "picomètre) est celui de la virgule flottante, pas de la convergence."],
      ["Période exacte par intégration de l'énergie",
       "Perturbation au premier ordre (Airy) appliquée à l'isochronisme",
       "Incohérent par construction : Airy suppose l'oscillateur harmonique, or l'effet "
       "cherché est précisément son anharmonicité."],
      ["Points de rebroussement à énergie égale, deux branches traitées séparément",
       "Symétrisation de U(θ)",
       "La symétrisation supprime la contribution de l'asymétrie armage/désarmage, "
       "vérifiée non nulle et conforme au second ordre de Duffing."],
      ["Appariement des rebroussements par résolution exacte",
       "Interpolation d'une table énergie → amplitude",
       "Plancher de bruit ramené de 0,2–0,4 s/j à 3,6·10⁻⁹ s/j sur le cas linéaire de "
       "référence."],
      ["Tracé côté rouage, repère 6H EN HAUT et tige à DROITE",
       "Tracé 12H en haut, tige à gauche",
       "C'est l'attitude dans laquelle la pièce est examinée, et celle de la position "
       "COSC 6H : la vue par défaut est donc directement celle du critère noté. La "
       "chiralité est la même dans les deux cas — l'erreur porte sur un demi-tour, "
       "invisible sur la dispersion, et n'inverse que le SIGNE du DVH."],
      ["Directions de g DÉDUITES des positions des repères horaires (g = repère + 180°)",
       "Deux tables tenues à jour séparément",
       "Le nom d'une position COSC désigne le repère qui est EN HAUT, donc la direction "
       "opposée à g. Tenir les deux tables à la main est précisément ce qui avait laissé "
       "passer un demi-tour d'étiquetage."],
      ["Bout de la ligne des centres portant le collet donné explicitement",
       "Collet déduit du seul angle de la ligne des centres",
       "Une droite n'est définie qu'à 180° près par son angle. Sans cette donnée, le "
       "spiral peut être posé en demi-tour : module de sensibilité et dispersion "
       "inchangés, DVH de signe opposé."],
      ["|DVH| ajouté à la dispersion dans le MÊME objectif",
       "Troisième objectif du front de Pareto, ou DVH ignoré",
       "Les deux vont dans le même sens, donc aucun arbitrage à exposer. Mais la "
       "dispersion ne voit que le module du vecteur de sensibilité, le DVH voit aussi sa "
       "phase : la somme fait tourner la phase pour reporter l'effet sur 3H/9H plutôt que "
       "sur la position notée."],
      ["Ligne des centres tracée en DEMI-DROITE, du centre vers l'ancre",
       "Diamètre entier ; repère du point de comptage affiché",
       "De l'autre côté du balancier il n'y a rien : un diamètre suggère une symétrie "
       "inexistante. Le point de comptage n'est plus tracé, le nombre de tours n'étant "
       "pas forcément entier — sa position angulaire n'a alors rien d'imposé."],
      ["Compensation d'échappement traitée en contrainte",
       "Troisième objectif du front de Pareto",
       "C'est une cible, pas une grandeur à minimiser ; un front 3D se diluerait "
       "fortement à population constante. Mesuré : la population converge sur la cible en "
       "15 à 20 générations."],
      ["Courbes terminales en Bézier à points de contrôle libres",
       "Formes classiques imposées (arc, double arc)",
       "La théorie de Phillips ne prescrit aucune forme, seulement une condition sur le "
       "centre de gravité. Les formes classiques sont des contraintes de pliage manuel, "
       "sans objet en silicium DRIE."],
      ["Jeu évalué sur tout le balayage angulaire",
       "Jeu évalué aux seules amplitudes extrêmes (0, ±max)",
       "L'élastica n'est pas monotone en angle près du minimum de jeu : des collisions "
       "réelles n'étaient pas détectées."],
      ["Grille d'amplitude interne fixe à 5° pour l'objectif",
       "Grille d'affichage choisie par l'utilisateur",
       "Une grille grossière peut manquer le pire cas entre deux points — écart constaté "
       "entre la dispersion annoncée par l'optimiseur et celle relue sur le graphe."],
      ["Profils portés par PCHIP", "Spline cubique classique",
       "PCHIP ne produit ni dépassement ni oscillation entre nœuds : une épaisseur reste "
       "positive et plausible pour tout tirage de l'optimiseur."],
      ["Mise à l'échelle analytique de l'épaisseur (β)",
       "Fréquence traitée comme contrainte d'optimisation",
       "La fréquence cible est atteinte par construction (3,5000000 Hz mesuré), ce qui "
       "libère entièrement l'optimiseur pour la recherche de forme."],
      ["Module d'Young scalaire E = 169 GPa", "Tenseur d'élasticité anisotrope complet",
       "L'orientation cristalline ⟨111⟩ est transversalement isotrope dans le plan ; le "
       "module y est constant quelle que soit la direction dans le plan."],
      ["Arc concentrique + courbe en S pour la terminaison à raquette",
       "Courbe de Phillips sous les goupilles",
       "Sur toute forme non concentrique, déplacer la raquette pousse le spiral "
       "radialement au lieu de ne changer que sa longueur active."],
      ["Arc décrit par deux angles indépendants de part et d'autre du point de comptage",
       "Un seul angle, point de comptage à une extrémité de l'arc",
       "Les prendre égaux centre le point de comptage sur l'arc, donc égalise la course "
       "de la raquette dans les deux sens ; les dissocier permet d'éloigner le piton sans "
       "déplacer le point de comptage."],
      ["Piton placé au-delà du parcours des goupilles, par un paramètre indépendant",
       "Piton à l'extrémité du parcours des goupilles",
       "Avec 124 µm de lame après les goupilles, l'excursion latérale n'est que de "
       "0,9 µm : la lame ne touche jamais les goupilles et la raquette n'a aucune "
       "autorité. L'excursion varie comme le carré de cette longueur."],
      ["Condition limite au point de comptage, portion morte exclue du calcul",
       "Encastrement au piton",
       "La portion au-delà des goupilles ne travaille pas (déplacement nul vérifié à la "
       "précision machine). L'inclure fausserait la fréquence de ~0,2 %, soit 150 s/j. "
       "L'exclure du centre de gravité est EXACT et non approché : la formule d'Airy "
       "n'utilise que la dérivée du moment statique, et celui de la portion morte est "
       "constant."],
      ["Jeu de goupilles par fonction descriptive de zone morte",
       "Airy, ou intégration exacte de la période",
       "La non-linéarité est par morceaux : Airy suppose un oscillateur harmonique, "
       "l'intégration de période exige une énergie régulière. L'équilibrage harmonique "
       "est l'outil adapté, et il redonne exactement le piton fixe quand le jeu est nul."],
      ["Épaisseur d'oxyde recalculée pour chaque candidat, couplée à la mise à l'échelle",
       "Oxyde ajouté après l'optimisation",
       "Mesuré : +1449 s/j d'erreur de fréquence si l'oxyde est ajouté après coup, et "
       "4 µm de jeu entre spires perdus sans que le calcul le voie."],
      ["Plis du ruban cherchés aux extrema de l'ANGLE POLAIRE",
       "Extrema de rayon ; ou pic de courbure",
       "θ croît strictement le long d'un flanc : son maximum est au bout du ruban par "
       "construction. Le rayon est maximal sur le flanc extérieur, donc APRÈS la pointe — "
       "mesuré 102° après, d'où autant de spiral amputé. Le pic de courbure, lui, se "
       "laissait tromper par les accessoires du dessin."],
      ["Tolérance de recollement adaptative, validée par l'égalité des longueurs de flancs",
       "Tolérance fixe de 0,1 mm",
       "0,1 mm vaut trois fois l'épaisseur de la lame : assez pour souder au spiral une "
       "ligne de construction et ramener le rayon de collet de 0,468 à 0,249 mm sans rien "
       "signaler. L'aller-retour export/import passe de 11,80 à 0,04 µm d'erreur."],
      ["Flancs densifiés le long de leurs segments avant appariement angulaire",
       "Appariement sur le dessin tel quel",
       "r(θ) n'est pas linéaire le long d'une corde. Sur une portion rectiligne dessinée "
       "grossièrement, l'épaisseur — différence de deux rayons interpolés — tombait à "
       "10 µm pour 36 réels. Les points ajoutés sont sur la droite : la géométrie est "
       "inchangée."],
      ["Chaînage des entités par arbre k-d", "Balayage du pool restant",
       "Coût cubique en nombre d'entités : plusieurs minutes contre 0,9 s sur un export "
       "CAO de 6 087 segments élémentaires, l'application paraissant figée."],
      ["Densité de points des DXF indépendante du maillage de calcul",
       "DXF hérité du découpage de calcul",
       "Les deux répondent à des exigences sans rapport : le maillage de calcul est un "
       "compromis précision/temps qui densifie les zones terminales, un fichier de "
       "gravure doit décrire le contour uniformément. Ré-échantillonner ne change pas la "
       "courbe (0,006 % de longueur à 5 000 points), seulement sa description."],
      ["Deux fichiers DXF : cote de gravure et pièce finie", "Un seul fichier",
       "Les deux diffèrent de 2 × 0,56 × épaisseur d'oxyde, soit ~4 µm, un tiers du jeu "
       "minimal exigé. Les confondre donnerait une lame trop épaisse d'autant."],
      ["Nombre de tours compté jusqu'au point de comptage",
       "Compté jusqu'au bout de la matière",
       "Convention horlogère : un nombre entier place le point de comptage sur la ligne "
       "des centres, ce qui y aligne l'axe de la force exercée sur le pivot. La règle "
       "reste facultative, le nombre de tours pouvant être fractionnaire."]],
     [26, 24, 50]),

    # ================================================================ 6
    ("saut",),
    ("h1", "6. Domaine de validité et incertitudes numériques"),
    ("p", "Cette section est la plus importante pour une validation croisée : elle indique "
          "à quel niveau les chiffres produits sont significatifs, et à partir de quand ils "
          "ne le sont plus."),

    ("h2", "6.1 Convergence en maillage"),
    ("p", "Le maillage physique est le nombre d'éléments discrets le long de la lame. "
          "L'étude de convergence a été conduite sur un spiral d'Archimède muni de ses "
          "deux courbes terminales, avec un pas de balayage angulaire de 20°."),
    ("tab",
     ["Nombre d'éléments", "Fréquence (Hz)", "Jeu min. (µm)", "Dispersion (s/j)",
      "DVH (s/j)", "Isochronisme à 200° (s/j)"],
     [["1 670", "3,5000000", "21,72", "1,929", "−0,560", "+20,755"],
      ["3 340", "3,5000000", "22,58", "1,853", "−0,457", "+19,317"],
      ["6 680", "3,5000000", "22,99", "1,815", "−0,406", "+18,682"],
      ["13 360", "3,5000000", "23,19", "1,796", "−0,380", "+18,384"]],
     [16, 17, 14, 16, 13, 24]),
    ("p", "Les écarts successifs sur l'isochronisme (1,44 puis 0,64 puis 0,30 s/j) se "
          "réduisent de moitié à chaque doublement : la convergence est du premier ordre, "
          "et l'extrapolation situe la limite vers +18,1 s/j. Au seuil de 1 % de la valeur "
          "la plus fine, seule la fréquence est acquise dès 1 670 éléments ; le jeu "
          "minimal l'est à 6 680, et la dispersion, le DVH et l'isochronisme propre ne le "
          "sont pas encore à 13 360."),
    ("p", "La fréquence propre ne converge pas : elle est EXACTE à tout "
          "maillage (3,5000000 Hz), parce que la mise à l'échelle analytique de "
          "l'épaisseur (§4.5) la cale par construction sur la géométrie discrétisée telle "
          "qu'elle est. Elle ne renseigne donc en rien sur la finesse du maillage, et ne "
          "doit pas être lue comme un indicateur de convergence. La dispersion de gravité "
          "demande davantage. L'isochronisme propre est le plus exigeant des trois : il se "
          "lit sur la COURBURE de la courbe d'énergie, donc sur une différence de "
          "différences. Il faut compter environ 3 000 éléments pour qu'il soit "
          "exploitable, contre environ 1 500 pour la dispersion.",
     "Hiérarchie de sensibilité. "),
    ("p", "Le DVH subit la MÊME erreur absolue que la "
          "dispersion — les deux se lisent sur la même orbite de centre de gravité — mais "
          "il est au moins deux fois plus petit qu'elle par construction (§3.3). La même "
          "erreur y devient donc bien plus lourde en relatif. Mesuré sur la configuration "
          "de référence, entre 1 670 et 6 680 éléments :",
     "Le DVH est le plus exposé de tous. "),
    ("tab",
     ["Grandeur", "1 670 éléments", "6 680 éléments", "Écart absolu", "Écart relatif"],
     [["Dispersion", "1,929 s/j", "1,815 s/j", "0,114 s/j", "6 %"],
      ["DVH (COSC)", "−0,560 s/j", "−0,406 s/j", "0,154 s/j", "38 %"],
      ["Isochronisme propre à 200°", "+20,755 s/j", "+18,682 s/j", "2,073 s/j", "11 %"]],
     [30, 18, 18, 17, 17]),
    ("p", "Vérifier la dispersion ne suffit donc PAS à valider le DVH : un maillage qui "
          "résout honnêtement la première peut laisser le second largement dans le bruit. "
          "À 1 670 éléments, la valeur du DVH bouge de plus d'un tiers — son signe tient, "
          "mais sa valeur n'est pas exploitable."),
    ("p", "Ce constat serait "
          "trompeur s'il s'arrêtait là. Le DVH n'est pas un nombre mais une COURBE en "
          "fonction de l'amplitude, et les deux ne convergent pas au même rythme. Mesuré "
          "sur la configuration de référence, entre 3 270 et 7 998 éléments :",
     "Mais la FORME converge bien plus vite que le NIVEAU. "),
    ("tab",
     ["Grandeur lue sur la courbe de DVH", "Valeur", "Incertitude de maillage"],
     [["Valeur à 200°", "−0,581 s/j", "0,088 s/j"],
      ["Valeur à 340°", "+0,170 s/j", "0,096 s/j"],
      ["Minimum de la courbe", "−0,929 s/j, vers 240°", "position stable"],
      ["Changement de signe", "vers 330°", "position stable"],
      ["Variation sur toute la plage", "1,099 s/j", "0,128 s/j au pire point"]],
     [42, 30, 28]),
    ("p", "L'incertitude de maillage (0,13 s/j) ne vaut que le huitième de la variation du "
          "DVH sur la plage (1,10 s/j) : les deux courbes se déplacent presque en bloc. "
          "L'amplitude du minimum, celle du changement de signe et la pente générale sont "
          "donc lisibles au maillage de travail, alors que la valeur à une amplitude "
          "donnée ne l'est pas."),
    ("p", "Conséquence pratique, et c'est la raison d'être de la vue dédiée dans le "
          "programme : un DVH doit se juger sur sa COURBE, pas sur un point. Le diagnostic "
          "compare l'incertitude de maillage à la VARIATION du DVH sur la plage et non à "
          "sa valeur, et avertit quand l'évolution elle-même cesse d'être lisible. Une "
          "validation croisée a tout intérêt à comparer d'abord la forme (position du "
          "minimum, amplitude du changement de signe), qui discrimine bien avant que les "
          "niveaux absolus ne soient convergés."),
    ("p", "Le programme rend cette incertitude explicite de deux façons. Un bouton "
          "recalcule fréquence, jeu, dispersion, DVH et isochronisme propre à maillage "
          "quadruplé, et affiche pour chacun l'écart absolu ET relatif ; trois "
          "avertissements se déclenchent automatiquement — si l'écart de maillage est du "
          "même ordre que la dispersion, que le DVH, ou s'il dépasse la tolérance de "
          "compensation demandée. Un second bouton retrace la marche en position 6H et à "
          "plat sur toute la plage d'amplitude, aux deux maillages superposés, avec leur "
          "différence et la bande d'incertitude : c'est là que se juge la forme."),

    ("h2", "6.2 Convergence en pas de balayage angulaire"),
    ("p", "Le pas de balayage est l'échantillonnage en angle du balancier, chaque angle "
          "nécessitant une résolution élastica complète. Après la correction du solveur de "
          "tir, la sensibilité est négligeable : l'isochronisme propre à 200° varie de "
          "0,005 s/j entre un pas de 30° et un pas de 5°. Avant cette correction, la "
          "dispersion était de ±10 s/j sur la même plage."),

    ("h2", "6.3 Biais de sélection de l'algorithme génétique"),
    ("p", "Un point à connaître pour interpréter les résultats. L'optimiseur évalue des "
          "milliers de candidats et retient le MINIMUM d'un critère bruité : il "
          "sélectionne donc préférentiellement les candidats dont l'erreur numérique est "
          "la plus favorable, et cet avantage disparaît dès qu'on recalcule dans d'autres "
          "conditions. La valeur annoncée par l'optimiseur est systématiquement optimiste."),
    ("p", "Règle pratique : ne comparer deux solutions que si leur écart dépasse "
          "l'incertitude de maillage affichée par la vérification haute résolution. En "
          "deçà, on départage du bruit numérique. Cette règle est d'autant plus stricte "
          "pour le DVH, dont l'incertitude relative est la plus forte (§6.1)."),

    ("h2", "6.4 Incertitude sur les coefficients thermiques"),
    ("p", "Après compensation, le modèle donne un résidu de 0,008 s/j crête-à-crête sur "
          "8–38 °C, c'est-à-dire zéro. IL NE FAUT PAS LE CROIRE : les modules d'Young sont "
          "pris linéaires en température, faute de données sur leur courbure, et c'est "
          "justement cette courbure qui limite les chronomètres réels. Les deux termes du "
          "second ordre sont exposés dans l'interface, à zéro, et signalés comme non "
          "renseignés."),
    ("p", "La vraie limite est ailleurs : dans l'incertitude sur les coefficients "
          "eux-mêmes, qui sont des valeurs de littérature. Chaque coefficient est perturbé "
          "seul, l'épaisseur d'oxyde étant recalculée (procédé recalé) ou non :"),
    ("tab",
     ["Coefficient", "Incertitude", "Épaisseur d'oxyde requise",
      "Pente résiduelle si non recalée"],
     [["Coefficient thermoélastique de l'oxyde", "±20 ppm/K", "3,38 – 3,95 µm",
       "±0,285 s/j/°C"],
      ["Module d'Young de l'oxyde", "±5 GPa", "3,48 – 3,83 µm", "±0,181 s/j/°C"],
      ["Coefficient thermoélastique du silicium", "±5 ppm/K", "3,50 – 3,79 µm",
       "±0,145 s/j/°C"],
      ["Dilatation du balancier", "±1 ppm/K", "3,56 – 3,73 µm", "±0,086 s/j/°C"],
      ["Dilatation du silicium", "±0,1 ppm/K", "3,64 – 3,65 µm", "±0,002 s/j/°C"]],
     [34, 16, 24, 26]),
    ("p", "Cumulées en quadrature : 0,377 s/j/°C si l'épaisseur calculée est appliquée "
          "sans être recalée. C'est sous le critère COSC mais très au-dessus d'un objectif "
          "serré. Deux conclusions : l'incertitude sur OÙ VISER (±0,29 µm) vaut six fois "
          "le contrôle de procédé nécessaire pour atteindre 0,05 s/j/°C (±0,05 µm), donc "
          "l'épaisseur doit être calée expérimentalement ; et le coefficient de l'oxyde "
          "domine tout le reste, c'est lui qu'il faut mesurer en priorité."),

    ("h2", "6.5 Hypothèses de modélisation restant à vérifier"),
    ("puce", "Poutre d'Euler-Bernoulli : allongement axial et cisaillement transverse "
             "négligés."),
    ("puce", "Comportement plan : aucun déversement hors-plan ni torsion. Pour une lame de "
             "100 µm de hauteur et 35 µm d'épaisseur, le rapport d'aspect rend l'hypothèse "
             "raisonnable, mais un calcul tridimensionnel la confirmerait utilement."),
    ("puce", "Le collet est modélisé comme un encastrement ponctuel astreint au cercle de "
             "rayon r_min ; la virole réelle a une raideur finie et une longueur "
             "d'encastrement non nulle."),
    ("puce", "Le piton est un encastrement parfait, sans jeu ni raideur de fixation."),
    ("puce", "r(θ) reste une fonction univoque de θ : le spiral ne se replie jamais "
             "localement sur lui-même. Cette hypothèse est vérifiée explicitement à chaque "
             "construction de géométrie."),
    ("puce", "Aucun effet dynamique : le calcul est quasi statique, l'inertie propre du "
             "spiral est négligée devant celle du balancier."),
    ("puce", "Le frottement du pivot du balancier dans sa pierre n'est pas modélisé, alors "
             "que le spiral l'y presse avec jusqu'à treize fois son propre poids (§3.8). "
             "C'est la seule omission qui commande une décision de conception."),
    ("puce", "Les coefficients thermoélastiques sont des valeurs de littérature, et leurs "
             "termes du second ordre sont pris nuls (§6.4)."),
    ("puce", "L'orientation du spiral dans le mouvement (§3.3) est une DONNÉE d'entrée, "
             "non une mesure : l'angle de la ligne des centres, le bout qui porte le "
             "collet et le sens de développement doivent être relevés sur le calibre "
             "réel. Une erreur de demi-tour y inverse le signe du DVH sans rien changer "
             "d'autre."),
]

_PARTIE_5 = [

    # ================================================================ 7
    ("saut",),
    ("h1", "7. Structure du programme"),
    ("p", "Le programme est un package Python unique, spiral_lab, sans dépendance entre "
          "l'analyse et l'optimisation : une géométrie importée d'un DXF et une géométrie "
          "produite par l'optimiseur empruntent exactement le même chemin de calcul. "
          "L'interface est une fenêtre Qt (PySide6) unique, les calculs longs s'exécutant "
          "en arrière-plan."),
    ("tab",
     ["Fichier", "Lignes", "Rôle"],
     [["main.py", "29",
       "Point d'entrée. Lance l'application Qt. Gère l'import absolu pour que le fichier "
       "soit exécutable directement depuis un IDE."],
      ["app.py", "2609",
       "Fenêtre unique de l'application : cinq panneaux de réglage (origine avec repère "
       "du mouvement, forme avec terminaison à raquette et compensation thermique, "
       "sondage, optimisation avec cible d'échappement, résultats) et huit onglets "
       "graphiques (géométrie, marche diurne, front de Pareto, épaisseur, jeu minimal, "
       "DVH, front de Pareto, épaisseur, jeu minimal, orbite du centre de gravité, "
       "thermique, animation). Contient aussi les deux vérifications à haute "
       "résolution — métriques à maillage quadruplé, et marche 6H / à plat sur toute la "
       "plage d'amplitude — les garde-fous de cohérence et l'enregistrement du projet."],
      ["modele.py", "251",
       "Représentation canonique unique d'un spiral (GeometrieSpiral) : fibre neutre "
       "x, y, épaisseur de gravure, hauteur, épaisseur d'oxyde, et nombre d'éléments "
       "élastiquement actifs. En dérive l'épaisseur finie, la rigidité du composite, la "
       "souplesse intégrée et les deux contours (gravure et pièce finie). Tout le reste "
       "du programme ne manipule que cet objet, quelle que soit l'origine de la "
       "géométrie. Porte aussi le ré-échantillonnage (§4.4), qui redécrit la même "
       "courbe avec un autre nombre de points sans déplacer le point de comptage. "
       "Convention : tout en millimètres, conversion SI locale au point d'usage."],
      ["constantes.py", "42",
       "Constantes matériau : silicium ⟨111⟩ (169 GPa, 2330 kg/m³) et oxyde thermique "
       "(70 GPa, 2200 kg/m³), coefficients thermoélastiques et de dilatation, rapports de "
       "consommation et de croissance de l'oxydation, pesanteur."],
      ["dxf_io.py", "491",
       "Import et export DXF. À l'import, rechaîne les entités du dessin (arbre k-d), "
       "isole le contour du spiral, le scinde en deux flancs aux extrema d'angle "
       "polaire, puis reconstruit la fibre neutre et l'épaisseur locale par mise en "
       "correspondance des deux flancs sur cet angle (§4.10). À l'export, régénère les "
       "deux flancs par décalage de ±e/2 le "
       "long de la normale locale, à une densité de points indépendante du maillage de "
       "calcul, en DEUX fichiers quand la lame est oxydée : cote de gravure pour la "
       "DRIE, pièce finie pour la vérification."],
      ["courbes.py", "182",
       "Courbes terminales de Grossmann et de Phillips : construction des points de "
       "contrôle de Bézier sous contrainte de continuité C1, test de monotonie par la "
       "propriété des polynômes de Bernstein, et calcul de la condition de Phillips "
       "d = R²/L comme repère de plausibilité."],
      ["parametrage.py", "566",
       "Génération paramétrique d'une géométrie complète : assemblage de la courbe de "
       "Grossmann, du corps à pas constant et de la terminaison extérieure (Phillips ou "
       "courbe en S + arc de raquette), application des profils d'épaisseur et d'écart "
       "radial (PCHIP), puis mise à l'échelle analytique de l'épaisseur. Porte aussi "
       "l'ORIENTATION dans le mouvement : angle de la ligne des centres, bout qui porte "
       "le collet, sens de développement (§3.3). Définit le vecteur de conception "
       "(ParametresSpiral)."],
      ["physique.py", "411",
       "Cœur du calcul. Solveur élastica non linéaire par méthode de tir (§3.1, §4.1, "
       "§4.2), incluant l'intégration prédicteur-correcteur compilée Numba et le calcul "
       "de l'énergie élastique. Fournit le balayage angulaire, qui produit d'un seul "
       "passage l'orbite du centre de gravité, le jeu inter-spires et la courbe "
       "d'énergie ; plus la plage de réglage de la raquette et l'axe de la force exercée "
       "sur le pivot (§3.8)."],
      ["collision.py", "36",
       "Jeu minimal entre spires adjacentes sur une géométrie déformée, sans hypothèse de "
       "pas constant (§4.6)."],
      ["isochronisme.py", "555",
       "Les quatre contributions à la marche diurne et leur composition : gravité par la "
       "formule d'Airy, isochronisme propre par intégration de la période exacte, "
       "échappement, jeu de raquette ; puis somme par position (quatre verticales, deux "
       "horizontales). Porte aussi le repère du mouvement — positions des repères "
       "horaires, directions de g qui s'en déduisent, DVH et vecteur de sensibilité "
       "(§3.3). Point d'entrée unique d'analyse d'une géométrie."],
      ["echappement.py", "117",
       "Modèle de l'erreur d'échappement : loi en 1/Θ calée sur un retard spécifié à une "
       "amplitude basse, et réglage de la cible de compensation, y compris le choix du "
       "montage visé (§3.5)."],
      ["raquette.py", "149",
       "Erreur de marche due au jeu entre les goupilles : fonction descriptive de la zone "
       "morte, et mesure de l'excursion latérale de la lame au droit des goupilles dans "
       "l'état libre (§3.9, §4.7)."],
      ["thermique.py", "268",
       "Compensation thermique : rigidité du composite silicium/oxyde en fonction de la "
       "température, recherche de l'épaisseur d'oxyde annulant la variation de marche, "
       "marche résiduelle, et étude de sensibilité aux coefficients matériau "
       "(§3.10, §6.4)."],
      ["optimisation.py", "450",
       "Branchement NSGA-II : construction des bornes, évaluation d'un candidat "
       "(parallélisée par joblib), définition du problème pymoo à deux objectifs et deux "
       "à trois contraintes, sondage aléatoire préalable et construction de la population "
       "initiale du warm start (§4.9)."],
      ["presets.py", "204",
       "Dataclasses de configuration, séparées par nature : spécification du mouvement "
       "(y compris son repère), bornes de recherche, hyperparamètres de l'algorithme "
       "génétique. Sauvegarde et rechargement JSON de la configuration complète, avec "
       "compatibilité ascendante."],
      ["rapport.py", "747",
       "Rapport de simulation Word et enregistrement du projet complet : configuration, "
       "résultats numériques, quatre états caractéristiques du spiral dans le repère du "
       "mouvement, quatre vues d'attitude COSC, orbite du centre de gravité, marche, "
       "épaisseur, température. Écrit dans un dossier unique le rapport, l'animation GIF, "
       "les DXF et la configuration JSON."],
      ["rendu.py", "202",
       "Tracés communs à l'interface et au rapport : repère du mouvement (repères "
       "horaires, axe de tige, demi-droite de la ligne des centres, collet), vecteur de "
       "gravité, calcul des images du déploiement et encodage GIF. Séparé de physique.py, "
       "qui doit rester utilisable sans dépendance graphique."],
      ["convergence.py", "105",
       "Étude automatisée de convergence du maillage sur les métriques de marche, DVH "
       "compris, avec rapport de convergence à seuil paramétrable (§6.1)."],
      ["threads.py", "106",
       "QThread pour le sondage, l'optimisation et le calcul d'animation, afin de ne pas "
       "geler l'interface. Transmet la progression génération par génération."],
      ["tests/test_spiral_lab.py", "781",
       "Suite de vérification consolidée, exécutable seule : références analytiques "
       "(ressort linéaire, Duffing), composition des marches, repère du mouvement et DVH, "
       "solveur elastica, terminaison à raquette, thermique, DXF, configuration, "
       "optimisation, import et export DXF, ré-échantillonnage, rapport. Environ 175 "
       "contrôles."]],
     [20, 9, 71]),
    ("p", "Total : environ 7 500 lignes réparties sur 19 modules, plus une suite de "
          "tests de 780 lignes. Un dossier outils/ porte le générateur du présent "
          "document, dont le contenu est versionné comme des données."),

    ("h2", "7.1 Chaîne de calcul"),
    ("p", "Le flux de données, de la définition d'une forme jusqu'à un chiffre de marche :"),
    ("eq", ["ParametresSpiral   (vecteur de conception + repère du mouvement)",
            "      ↓  parametrage.construire_geometrie",
            "         Bézier + PCHIP, orientation sur la ligne des centres,",
            "         puis DEUX points fixes imbriqués :",
            "           arc de raquette <-> souplesse    et    oxyde <-> échelle β",
            "GeometrieSpiral  (x, y, e_gravure, h, oxyde, n_actif)  — OU dxf_io.charger_dxf",
            "      ↓  physique.balayage_angulaire     — une résolution élastica par angle,",
            "                                           sur la partie ACTIVE, E·I composite",
            "orbite du CG  +  jeu(θ)  +  énergie U(θ)  — un seul passage, trois exploitations",
            "      ↓  isochronisme.analyser_marches",
            "gravité (Airy) + isochronisme propre (période exacte)",
            "               + échappement (1/Θ) + jeu de goupilles (zone morte)",
            "               → marche par position, DVH, vecteur de sensibilité (A, φ)",
            "      ↓  optimisation.evaluer_candidat",
            "dispersion , jeu minimal , décentrage , écart de compensation ,",
            "             épaisseur d'oxyde , DVH (pire cas)",
            "      ↓  optimisation.ProblemeSpiral  →  NSGA-II",
            "front de Pareto"]),
    ("p", "Le balayage angulaire est le seul point coûteux ; il est mutualisé entre les "
          "exploitations qui en dépendent, ce qui évite plusieurs parcours redondants du "
          "même calcul."),

    # ================================================================ 8
    ("saut",),
    ("h1", "8. Points de comparaison recommandés pour la validation croisée"),
    ("p", "Les comparaisons sont classées du plus élémentaire au plus intégré. Il est "
          "fortement conseillé de les mener dans cet ordre : un écart apparu à une étape "
          "rend inexploitables toutes les suivantes."),
    ("tab",
     ["Étape", "Grandeur à comparer", "Valeur de référence Spiral Lab", "Remarque"],
     [["1", "Longueur développée, masse totale, moment quadratique local", "Géométrie seule",
       "Vérifie d'abord que les deux outils décrivent bien la même pièce."],
      ["2", "Fréquence propre en petites oscillations",
       "3,5000000 Hz à tout maillage",
       "Exacte par construction (§4.5) : elle vérifie que les deux outils décrivent la "
       "même raideur, mais ne dit RIEN sur la convergence du maillage."],
      ["3", "Forme déformée à un angle donné (ex. +200°, +340°)",
       "Résidu de fermeture 3,4·10⁻⁷ µm à 280° (3 340 éléments)",
       "Comparer la position du collet et l'allure générale, pas seulement le déplacement "
       "maximal."],
      ["4", "Énergie élastique U(α) sur toute la plage d'angle",
       "6,83·10⁻⁶ J à 340° (3 340 éléments)",
       "C'est LA grandeur critique. Un écart de 10⁻³ en relatif suffit à fausser "
       "complètement l'isochronisme."],
      ["5", "Position du centre de gravité en fonction de l'angle",
       "Orbite complète sur ±340°",
       "Sensible au bras de levier déformé : écart de 183 µm à 340° entre les deux "
       "formulations."],
      ["6", "Couple de rappel C(α) = dU/dα", "Déduit de U(α)",
       "Un outil qui donne directement le couple permet de court-circuiter l'étape 4."],
      ["7", "Isochronisme propre : marche en fonction de l'amplitude, à plat",
       "+18,4 s/j à 200° (13 360 éléments), référence 300°",
       "Se compare sans aucun modèle de gravité ni d'échappement. Le test le plus "
       "discriminant."],
      ["8", "Vecteur de sensibilité à la gravité : module A et phase φ",
       "A = 0,648 s/j, φ = 324,1° à 300° d'amplitude (13 360 éléments ; 0,754 s/j et "
       "312,1° à 1 670 éléments — grandeur lente à converger)",
       "Se compare SANS convention de position : A et φ décrivent le comportement pour "
       "toute direction de g. Vérifier d'abord que marche(θ) = A·cos(θ−φ) est bien "
       "reproduit, ce qui valide la linéarité du modèle de gravité."],
      ["9", "Dispersion de marche entre positions verticales",
       "1,796 s/j (13 360 éléments)",
       "Dépend en plus de la formule d'Airy et du maillage. Insensible au repère : c'est "
       "une permutation des quatre positions."],
      ["10", "DVH : marche en position COSC 6H moins marche à plat",
       "−0,380 s/j à 300° (13 360 éléments ; −0,560 à 1 670)",
       "À ne comparer qu'APRÈS avoir aligné le repère du §3.3 : le signe s'inverse d'un "
       "demi-tour. Vérifier au passage l'inégalité dispersion ≥ 2·|DVH|, qui ne dépend "
       "d'aucune convention."],
      ["11", "Forme de la courbe de DVH en fonction de l'amplitude",
       "Minimum −0,929 s/j vers 240°, changement de signe vers 330°, variation 1,10 s/j "
       "sur 200–340°",
       "Converge BEAUCOUP plus vite que le niveau (incertitude de maillage 0,13 s/j pour "
       "1,10 s/j de variation) : à comparer avant l'étape 10, elle discrimine alors que "
       "les valeurs absolues ne sont pas encore acquises."],
      ["12", "Jeu inter-spires minimal et angle auquel il survient", "Selon géométrie",
       "Vérifier que le minimum est cherché sur tout le balayage, et sur la lame FINIE "
       "(oxyde compris) et non sur la cote de gravure."],
      ["13", "Marche en fonction de la température, lame nue",
       "−3,848 s/j/°C, 115,4 s/j sur 8–38 °C",
       "Ne dépend que du coefficient thermoélastique et des dilatations : se compare sans "
       "aucun modèle d'oxyde."],
      ["14", "Épaisseur d'oxyde annulant la variation de marche",
       "3,683 µm (lame +4,13 µm, silicium consommé 3,21 µm)",
       "Dépend en plus du modèle de section composite. Comparer aussi la part de rigidité "
       "portée par l'oxyde : 33 % sur une lame de 31 µm."],
      ["15", "Résultante exercée sur le pivot du balancier",
       "99,4 µN à ±340°, axe à quelques degrés de la ligne des centres",
       "Sous-produit du solveur de tir. L'axe est une DROITE, à comparer modulo 180°. "
       "Utile pour alimenter un modèle de frottement de pivot que ce programme n'a pas."]],
     [8, 30, 26, 36]),
    ("p", "Les étapes 7 et 8 sont les plus utiles. La 7 ne fait "
          "intervenir ni la gravité, ni l'échappement, ni la formule d'Airy : seulement la "
          "mécanique de la poutre et l'intégration de la période — un désaccord y isole "
          "immédiatement un écart sur le modèle d'élastica. La 8 teste tout le modèle de "
          "gravité SANS engager la moindre convention de position ni de repère, ce qui "
          "évite de confondre une erreur de physique avec une erreur d'orientation. Ce "
          "n'est qu'une fois ces deux-là d'accord que les étapes 9 à 11 ont un sens — et "
          "parmi elles, la forme de la courbe de DVH (étape 11) se compare utilement bien "
          "avant son niveau (étape 10).",
     "Recommandation. "),

    ("h2", "8.1 Configuration de référence"),
    ("p", "Les valeurs ci-dessus correspondent à la configuration par défaut du "
          "programme :"),
    ("tab",
     ["Paramètre", "Valeur"],
     [["Rayon au collet r_min", "0,47 mm"],
      ["Rayon au piton r_max", "2,67 mm"],
      ["Nombre de tours réels", "12"],
      ["Nombre de tours fictifs (ligne de référence)", "14"],
      ["Hauteur hors-plan h", "0,100 mm (constante)"],
      ["Fréquence cible", "3,5 Hz"],
      ["Inertie du balancier J", "8·10⁻¹⁰ kg·m² (affichée 8 mg·cm² dans l'interface)"],
      ["Jeu inter-spires minimal requis", "0,012 mm"],
      ["Amplitude maximale", "340°"],
      ["Côté de vue", "rouage (le spiral est sous le coq)"],
      ["Orientation du tracé", "repère 6H en haut, tige de remontoir à droite"],
      ["Angle de la ligne des centres avec l'axe de tige", "−6,2°"],
      ["Bout de la ligne des centres portant le collet", "celui qui pointe vers l'ancre"],
      ["Sens de développement depuis le collet", "trigonométrique (vu côté rouage)"],
      ["Module d'Young du silicium ⟨111⟩", "169 GPa"],
      ["Masse volumique du silicium", "2330 kg/m³"],
      ["Retard d'échappement spécifié", "8 s/j à 200°, référence 300°"],
      ["Oxyde thermique", "70 GPa, 2200 kg/m³, +215 ppm/K, 0,5 ppm/K"],
      ["Coefficient thermoélastique du silicium", "−60 ppm/K"],
      ["Oxydation : consommation / croissance", "0,44 / 0,56 µm par µm d'oxyde"],
      ["Dilatation du balancier", "16,7 ppm/K (CuBe2)"],
      ["Plage de température visée", "8 à 38 °C, référence 23 °C"]],
     [52, 48]),
    ("p", "La configuration complète d'un spiral — spécification du mouvement et de son "
          "repère, vecteur de conception, bornes de recherche et hyperparamètres — se "
          "sauvegarde dans un fichier JSON unique, ce qui permet de transmettre exactement "
          "la géométrie évaluée ET son orientation. La géométrie elle-même s'exporte en "
          "DXF (deux flancs en LWPOLYLINE) pour reprise directe dans un mailleur."),

    ("h2", "8.2 Vérifications intégrées au programme"),
    ("p", "Certaines de ces comparaisons se font sans outil externe, le programme se "
          "contrôlant lui-même. Elles sont regroupées dans une suite exécutable "
          "(tests/test_spiral_lab.py, environ 175 contrôles) :"),
    ("puce", "le calcul de période non linéaire retombe à 3,6·10⁻⁹ s/j sur un ressort "
             "strictement linéaire, dont l'isochronisme est exact par construction ;"),
    ("puce", "il s'accorde à 0,13 % avec la formule de Duffing au premier ordre sur un "
             "ressort durcissant, et à trois chiffres avec le terme du second ordre sur un "
             "ressort asymétrique ;"),
    ("puce", "l'intégration en une passe du solveur elastica s'accorde à 10⁻¹² avec une "
             "résolution par point fixe menée jusqu'à convergence ;"),
    ("puce", "la section composite redonne E_silicium·I exactement quand l'épaisseur "
             "d'oxyde est nulle ;"),
    ("puce", "la marche de jeu de goupilles est identiquement nulle quand le jeu est nul, "
             "c'est-à-dire que le montage à piton fixe est retrouvé exactement ;"),
    ("puce", "la fréquence propre tombe sur 3,5000000 Hz dans toutes les configurations, "
             "oxyde imposé et raquette compris ;"),
    ("puce", "le réimport d'un DXF exporté redonne la même épaisseur à 0,01 µm près, pour "
             "chacune des deux cotes ;"),
    ("puce", "sur le repère du mouvement : chaque direction de g est diamétralement "
             "opposée à son repère horaire, le module de sensibilité est invariant par "
             "rotation à 10⁻⁶ près alors que le DVH, lui, varie, un demi-tour du spiral "
             "inverse exactement le signe du DVH sans toucher au module, et l'inégalité "
             "dispersion ≥ 2·|DVH| est vérifiée sur les trois variantes de géométrie ;"),
    ("puce", "un DXF fait de milliers de segments élémentaires, portant en plus un cercle "
             "de virole et une ligne de construction vers le centre, s'importe en moins "
             "d'une seconde sans que ces accessoires soient pris pour le spiral ; un ruban "
             "dont la terminaison rentre sur 120° est importé sans troncature angulaire, "
             "et une portion rectiligne dessinée en cordes grossières n'y creuse plus "
             "l'épaisseur ;"),
    ("puce", "le ré-échantillonnage conserve la longueur développée à 0,006 % près à "
             "5 000 points, laisse le point de comptage exactement sur le maillage, et "
             "ne fait jamais sortir l'épaisseur de sa plage d'origine."),

    ("h2", "8.3 Résultat d'optimisation obtenu"),
    ("p", "À titre indicatif, le résultat d'une optimisation complète (60 individus, 150 "
          "générations, tolérance de compensation 3 s/j, départ au spiral neutre) : 32 "
          "solutions non dominées, dont les meilleures présentent une dispersion de "
          "0,37 s/j, un jeu minimal de 55 à 64 µm et un isochronisme propre de +5,4 s/j à "
          "200°. Ce chiffre est antérieur à l'ajout du DVH dans l'objectif : il donne un "
          "ordre de grandeur de ce que la recherche de forme sait gagner, non la valeur "
          "que produirait le critère actuel."),
    ("p", "Ces chiffres appellent deux réserves explicites. D'abord, les solutions se "
          "tiennent toutes au bord de la bande de tolérance (−2,45 à −2,80 s/j pour une "
          "bande de ±3,00) : la contrainte de compensation est ACTIVE, l'optimiseur "
          "préfère sous-compenser parce que cela lui achète de la dispersion et du jeu. Il "
          "existe donc un arbitrage réel entre compenser l'échappement et minimiser la "
          "dispersion de gravité. Ensuite, une dispersion de 0,37 s/j doit être confrontée "
          "à l'incertitude de maillage avant d'être retenue (§6.1 et §6.3) — et le DVH qui "
          "l'accompagne davantage encore, son incertitude relative étant environ cinq fois "
          "plus forte."),
]

CONTENU = _PARTIE_1 + _PARTIE_2 + _PARTIE_3 + _PARTIE_4 + _PARTIE_5
