# app.py
"""
Fenêtre unique du laboratoire de spiral -- remplace spiran_v1 (analyse DXF)
ET spiropt_v1 (optimisation), qui étaient deux programmes séparés. Un seul
flux, toujours visible (pas d'onglets à deviner pour les réglages, cf. note
de conception de spiropt_v1/gui.py -- principe conservé) :

    1. Origine       : charger un DXF existant OU initialiser un spiral
                       d'Archimède depuis une spec de mouvement.
    2. Paramètres    : structure des courbes de Grossmann/Phillips (angle,
                       n_pas, n_points) et des profils libres (épaisseur,
                       écart radial) -- génère un aperçu.
    3. Sondage       : recherche d'un point de départ faisable.
    4. Optimisation  : NSGA-II, objectif = vraie dispersion de marche, sous
                       contrainte de jeu inter-spires ET de compensation de
                       l'erreur d'échappement.
    5. Résultats     : marche diurne par position (verticales ET horizontales,
                       avec les trois contributions séparées), export DXF,
                       validation par réimport (boucle fermée).

Les tracés (géométrie / marche diurne / front de Pareto) sont dans des
onglets À DROITE -- ce sont des VUES, pas des réglages, donc les onglets n'y
posent pas le problème que la note de spiropt_v1/gui.py visait à éviter.
"""
import json
import os
from dataclasses import asdict, replace

import numpy as np
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QDoubleSpinBox, QSpinBox, QScrollArea, QSplitter,
    QTabWidget, QTableWidget, QTableWidgetItem, QFileDialog, QMessageBox,
    QAbstractItemView, QSlider, QApplication, QHeaderView, QCheckBox, QComboBox,
    QInputDialog, QLineEdit,
)
from PySide6.QtCore import Qt, QTimer

from commun import navigation, theme
from .parametrage import ParametresSpiral, construire_geometrie, epaisseur_theorique
from .presets import (
    SpecMouvement, BornesOptimisation, HyperparametresGA,
    sauvegarder_configuration, charger_configuration,
)
from .physique import (
    resoudre_elastica, balayage_angulaire, frequence_propre, marche_raquette,
    axe_force_pivot, angle_point_comptage,
)
from .constantes import RHO_SILICIUM
from .isochronisme import (analyser_marches, POSITIONS_HORIZONTALES,
                            positions_verticales)
from .echappement import ReglageEchappement
from .thermique import ReglageThermique, analyser_thermique, etude_sensibilite
from . import ancien_format, dxf_io, rapport, sauvegarde
from .optimisation import vecteur_depuis_params, params_depuis_vecteur, PAS_DISPERSION_DEG_DEFAUT
from .threads import SondageThread, OptimisationThread, AnimationThread
from .rendu import (sauvegarder_gif_deploiement, indices_aller_retour,
                     tracer_repere_mouvement)

_LARGEUR_MAX_CHAMP = 110

# Racine des projets enregistrés, dans le dossier du paquet : un emplacement
# fixe et prévisible, plutôt que le dossier courant, qui dépend de l'endroit
# d'où l'application a été lancée.
_RACINE_RESULTS = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "Results")

# L'inertie du balancier vaut ~8e-10 kg.m^2 : illisible dans un champ de saisie,
# et facile à saisir de travers d'un facteur dix. Elle est donc AFFICHÉE en
# mg.cm^2, unité d'usage en horlogerie, où le calibre de référence tombe
# exactement sur 8 — un nombre qu'on relit sans le compter.
#
# Tout le reste du programme — spec, physique, rapport, JSON — reste en kg.m^2,
# l'unité SI : la conversion n'existe qu'ICI, aux deux fonctions ci-dessous, et
# nulle part ailleurs. C'est la seule façon de ne pas semer des facteurs 1e10
# dans le code. L'aller-retour est exact à la résolution du champ (3 décimales
# en mg.cm^2 = 1e-13 kg.m^2, soit celle de l'ancien champ à 13 décimales).
_J_KGM2_VERS_MGCM2 = 1.0e10    # 1 kg.m^2 = 1e6 mg x (100 cm)^2 = 1e10 mg.cm^2


def _inertie_affichee(j_kg_m2: float) -> float:
    """kg.m^2 (interne) -> mg.cm^2 (affiché)."""
    return j_kg_m2 * _J_KGM2_VERS_MGCM2


def _inertie_interne(j_mg_cm2: float) -> float:
    """mg.cm^2 (affiché) -> kg.m^2 (interne). Seule porte d'entrée."""
    return j_mg_cm2 / _J_KGM2_VERS_MGCM2


def _spinbox(minimum, maximum, valeur, pas=1.0, decimales=3):
    sb = QDoubleSpinBox()
    sb.setRange(minimum, maximum)
    sb.setDecimals(decimales)
    sb.setSingleStep(pas)
    sb.setValue(valeur)
    sb.setMaximumWidth(_LARGEUR_MAX_CHAMP)
    return sb


def _intbox(minimum, maximum, valeur):
    sb = QSpinBox()
    sb.setRange(minimum, maximum)
    sb.setValue(valeur)
    sb.setMaximumWidth(_LARGEUR_MAX_CHAMP)
    return sb


def _texte_ecourte(etiquette: QLabel, texte: str, maxi: int = 40):
    """Affiche un texte long en le RACCOURCISSANT par le milieu, avec le texte
    entier en info-bulle.

    Le retour à la ligne ne coupe que sur les espaces : un nom de fichier, qui
    n'en contient pas, reste un seul « mot » et impose sa largeur au panneau
    entier -- un chemin de 70 caractères le faisait passer de 374 à 418 px et
    faisait apparaître une barre de défilement horizontale. Raccourcir par le
    MILIEU garde le début et l'extension, les deux parties qui identifient le
    fichier."""
    mots = texte.split()
    if all(len(mot) <= maxi for mot in mots):
        etiquette.setToolTip("")
        etiquette.setText(texte)
        return
    # Seuls les MOTS trop longs sont raccourcis, pas la phrase entière : le
    # retour à la ligne se charge du reste, et tronquer tout le texte rendait
    # illisible un message dont seul le chemin posait problème.
    garde = (maxi - 3) // 2
    etiquette.setToolTip(texte)
    etiquette.setText(" ".join(
        mot if len(mot) <= maxi else f"{mot[:garde]}...{mot[-garde:]}" for mot in mots))


def _relatif(ecart: float, reference: float) -> str:
    """« Δ 0.15 s/j, 37 % » : l'écart de maillage, absolu puis rapporté à la
    valeur convergée. Sans le pourcentage, deux nombres voisins à l'oeil
    peuvent cacher une dérive qui dépasse la grandeur mesurée."""
    if not np.isfinite(reference) or abs(reference) < 1e-12:
        return f"Δ {abs(ecart):.2f}"
    return f"Δ {abs(ecart):.2f}, {100.0 * abs(ecart) / abs(reference):.0f} %"


def _relatif(ecart: float, reference: float) -> str:
    """« Δ 0.15 s/j, 37 % » : l'écart de maillage, absolu puis rapporté à la
    valeur convergée. Sans le pourcentage, deux nombres voisins à l'oeil
    peuvent cacher une dérive qui dépasse la grandeur mesurée."""
    if not np.isfinite(reference) or abs(reference) < 1e-12:
        return f"Δ {abs(ecart):.2f}"
    return f"Δ {abs(ecart):.2f}, {100.0 * abs(ecart) / abs(reference):.0f} %"


def _assurer_extension(chemin: str, extension: str) -> str:
    """Complète l'extension d'un chemin choisi via QFileDialog si elle est
    absente. Certaines boîtes de dialogue (selon l'OS/le gestionnaire de
    fichiers) ne complètent pas automatiquement l'extension du filtre choisi
    -- constaté concrètement sur l'export GIF (PIL levait une ValueError
    "unknown file extension"), le même risque existe pour tout fichier
    sauvé (DXF, JSON) : mieux vaut le corriger partout, une fois, en
    fonction du contexte (extension attendue) plutôt qu'au cas par cas."""
    if not chemin.lower().endswith(extension.lower()):
        return chemin + extension
    return chemin


class LaboSpiralApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Laboratoire Spiral -- Grossmann/Phillips, discrétisé, NSGA-II")
        self.resize(1650, 950)

        self.spec = SpecMouvement()
        self.bornes = BornesOptimisation()
        self.hyperparams = HyperparametresGA()
        self.base: ParametresSpiral | None = None
        self.geom = None
        # Ancienne expérience (format v4-v7) reprise : son fichier et la
        # géométrie qui en est issue, pour la reposer si l'orientation change.
        self._ancienne = None
        self._a_valider = {}        # widget -> libellé, champs encore à valider
        self.sondage_resultats = None
        self.solutions_pareto = []  # liste de dicts {x, dispersion, jeu}
        self._thread = None

        splitter = QSplitter(Qt.Horizontal)
        self.scroll_gauche = None
        splitter.addWidget(self._construire_panneau_gauche())
        splitter.addWidget(self._construire_panneau_droit())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([420, 1230])  # largeur initiale, librement ajustable à la souris
        splitter.setChildrenCollapsible(False)
        self.setCentralWidget(splitter)

        self._rafraichir_base_depuis_spec()

    # ------------------------------------------------------------------
    # Construction de l'interface
    # ------------------------------------------------------------------
    def _construire_panneau_gauche(self) -> QWidget:
        # Les réglages sont répartis en TROIS ONGLETS, un par étape du travail.
        # Empilés dans un seul panneau déroulant, ils faisaient plus de deux
        # mille pixels de haut : on passait son temps à faire défiler, et les
        # commandes d'analyse étaient toujours hors de vue de la table des
        # solutions.
        #
        # Les outils d'analyse sont construits DEUX FOIS, dans l'onglet
        # « Calculs » et dans l'onglet « Optimisation » : analyser une solution
        # du front de Pareto ne doit pas obliger à changer d'onglet, puis à
        # revenir pour en sélectionner une autre. Les boutons dupliqués
        # appellent les mêmes gestionnaires et écrivent dans les mêmes
        # résultats (voir _enregistrer_commande).
        self._commandes = {}
        self._affichages = {}

        self.onglets_reglages = QTabWidget()
        self.onglets_reglages.setDocumentMode(True)

        forme = QWidget()
        v_forme = QVBoxLayout(forme)
        v_forme.setContentsMargins(2, 4, 2, 2)
        v_forme.addWidget(self._panneau_origine())
        v_forme.addWidget(self._panneau_parametres())
        v_forme.addWidget(self._panneau_sauvegarde_forme())
        v_forme.addStretch(1)
        self.onglets_reglages.addTab(forme, "1. Forme")

        calculs = QWidget()
        v_calc = QVBoxLayout(calculs)
        v_calc.setContentsMargins(2, 4, 2, 2)
        v_calc.addWidget(self._panneau_echappement())
        v_calc.addWidget(self._boite_analyse("calc"))
        v_calc.addWidget(self._panneau_sauvegarde_calculs())
        v_calc.addStretch(1)
        self.onglets_reglages.addTab(calculs, "2. Calculs")

        optim = QWidget()
        v_opt = QVBoxLayout(optim)
        v_opt.setContentsMargins(2, 4, 2, 2)
        v_opt.addWidget(self._panneau_sondage())
        v_opt.addWidget(self._panneau_optimisation())
        v_opt.addWidget(self._boite_analyse("optim"))
        v_opt.addWidget(self._panneau_sauvegarde_optimisation())
        v_opt.addStretch(1)
        self.onglets_reglages.addTab(optim, "3. Optimisation")

        contenu = QWidget()
        layout = QVBoxLayout(contenu)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.addWidget(self.onglets_reglages)
        layout.addStretch(0)

        # Rendre le panneau COMPRESSIBLE. Sans cela son contenu imposait une
        # largeur minimale de 496 px alors que la zone visible en fait 404 :
        # tout débordait, tableaux compris (mesuré 2026-09-28 -- les tableaux
        # eux-mêmes n'étaient pas en cause, leur minimum est de 70 px).
        # Deux réglages, appliqués ici en une fois plutôt que formulaire par
        # formulaire, pour qu'ils valent aussi pour les panneaux ajoutés plus
        # tard :
        #   - WrapLongRows : quand la ligne ne tient pas, l'étiquette passe
        #     AU-DESSUS du champ. Le minimum devient le plus grand des deux au
        #     lieu de leur somme.
        #   - les listes déroulantes se dimensionnent sur une longueur minimale
        #     et non sur leur item le plus long (218 px pour « Piton fixe au
        #     point de comptage ») ; le texte complet reste lisible déroulé.
        for formulaire in contenu.findChildren(QFormLayout):
            formulaire.setRowWrapPolicy(QFormLayout.WrapLongRows)
        # Retour à la ligne sur TOUTES les étiquettes du panneau. Sans cela, une
        # étiquette qui reçoit un texte long en cours de route impose sa largeur
        # au panneau entier, et le splitter l'élargit pour la satisfaire : le
        # commentaire du sondage exigeait à lui seul 908 px, faisant passer le
        # minimum du contenu de 374 à 950 px (mesuré 2026-09-28). Appliqué en
        # une fois plutôt qu'étiquette par étiquette, pour valoir aussi pour
        # celles ajoutées plus tard.
        for etiquette in contenu.findChildren(QLabel):
            etiquette.setWordWrap(True)
        for liste in contenu.findChildren(QComboBox):
            liste.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
            liste.setMinimumContentsLength(10)

        # Pas de largeur MAXIMALE ici : elle empêchait d'élargir le panneau en
        # tirant la poignée du splitter, donc de lire un texte qui dépasse
        # (tableaux, résultats numériques). Seule une largeur minimale est
        # imposée, la largeur initiale étant fixée par le splitter lui-même
        # (voir __init__ : setSizes).
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(contenu)
        # Largeur minimale DÉDUITE du contenu plutôt que fixée en dur : en
        # dessous, la poignée du splitter pourrait amener le panneau dans un
        # état où une barre de défilement HORIZONTALE apparaît, c'est-à-dire
        # exactement le débordement qu'on cherche à éviter. Recalculée à la
        # construction, elle suit les panneaux ajoutés plus tard.
        marge = scroll.verticalScrollBar().sizeHint().width() + 8
        # minimumSizeHint et NON sizeHint : un QLabel à retour à la ligne
        # annonce comme taille SOUHAITÉE la largeur de son texte sur une seule
        # ligne (Qt ne sait pas d'avance à quelle largeur il sera replié), ce
        # qui fait bondir le sizeHint du panneau de 374 à 616 px dès qu'un
        # résultat s'affiche. Le minimumSizeHint, lui, tient compte du repli :
        # c'est la seule largeur en dessous de laquelle une barre de défilement
        # horizontale apparaîtrait vraiment. Vérifié à 388/434/504/684 px :
        # aucune barre horizontale, le texte se replie (2026-09-29).
        scroll.setMinimumWidth(contenu.minimumSizeHint().width() + marge)
        self.scroll_gauche = scroll
        return scroll

    def _onglet_avec_barre(self, canvas: FigureCanvas) -> QWidget:
        """Enveloppe un canvas avec la barre d'outils matplotlib standard
        (zoom rectangle, pan, retour vue précédente, sauvegarde image...)."""
        conteneur = QWidget()
        v = QVBoxLayout(conteneur)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(NavigationToolbar(canvas, conteneur), 0)
        navigation.activer(canvas)  # molette, glisser, double-clic : comme la vue de la came
        v.addWidget(canvas, 1)
        return conteneur

    def _construire_panneau_droit(self) -> QWidget:
        self.tabs = QTabWidget()

        # Géométrie : vue d'ensemble + 2 zooms dédiés (Grossmann/collet et
        # Phillips/piton), car ces deux courbes ne couvrent que quelques
        # dizaines à quelques centaines de microns sur un spiral de plusieurs
        # mm de rayon -- invisibles sur la vue d'ensemble seule.
        self.fig_geom = Figure(figsize=(9, 6))
        gs = self.fig_geom.add_gridspec(2, 2, width_ratios=[2, 1])
        self.ax_geom = self.fig_geom.add_subplot(gs[:, 0])
        self.ax_zoom_grossmann = self.fig_geom.add_subplot(gs[0, 1])
        self.ax_zoom_phillips = self.fig_geom.add_subplot(gs[1, 1])
        self.canvas_geom = FigureCanvas(self.fig_geom)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_geom), "Géométrie")

        self.fig_marche = Figure(figsize=(6, 5))
        self.ax_marche = self.fig_marche.add_subplot(111)
        self.canvas_marche = FigureCanvas(self.fig_marche)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_marche), "Marche diurne")

        # Onglet DVH : la position 6H, la position à plat et leur différence, en
        # fonction de l'amplitude. Le DVH étant la grandeur la plus lente à
        # converger (même erreur absolue que la dispersion, pour une valeur au
        # moins deux fois plus petite), il lui faut sa propre vue, où le
        # maillage de travail et un maillage raffiné se superposent : c'est
        # l'écart entre les deux courbes qui dit si la forme du DVH est résolue,
        # pas la valeur à une seule amplitude.
        self.fig_dvh = Figure(figsize=(6, 6))
        gs_dvh = self.fig_dvh.add_gridspec(2, 1, height_ratios=[1.3, 1])
        self.ax_dvh_positions = self.fig_dvh.add_subplot(gs_dvh[0])
        self.ax_dvh_ecart = self.fig_dvh.add_subplot(gs_dvh[1],
                                                      sharex=self.ax_dvh_positions)
        self.canvas_dvh = FigureCanvas(self.fig_dvh)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_dvh), "DVH")

        self.fig_pareto = Figure(figsize=(6, 5))
        self.ax_pareto = self.fig_pareto.add_subplot(111)
        self.canvas_pareto = FigureCanvas(self.fig_pareto)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_pareto), "Front de Pareto")

        self.fig_epaisseur = Figure(figsize=(6, 4))
        self.ax_epaisseur = self.fig_epaisseur.add_subplot(111)
        self.canvas_epaisseur = FigureCanvas(self.fig_epaisseur)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_epaisseur), "Épaisseur")

        self.fig_jeu_min = Figure(figsize=(6, 6))
        self.ax_jeu_min = self.fig_jeu_min.add_subplot(111)
        self.canvas_jeu_min = FigureCanvas(self.fig_jeu_min)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_jeu_min), "Jeu minimal")

        # Orbite du centre de gravité : la grandeur que l'optimisation minimise
        # réellement. Deux vues — la trajectoire dans le plan, dont c'est
        # l'ÉTENDUE qui crée l'erreur de marche en position verticale, et les
        # deux composantes en fonction de l'angle du balancier, où se lisent
        # les accidents locaux que la trajectoire seule masque.
        # Disposition en DEUX LIGNES et non deux colonnes : l'orbite est très
        # allongée (typiquement 270 µm en x pour 60 µm en y) et doit être
        # tracée à aspect égal pour ne pas mentir sur sa forme. Côte à côte,
        # elle se réduisait à un trait écrasé.
        self.fig_cg = Figure(figsize=(8, 7))
        gs_cg = self.fig_cg.add_gridspec(2, 2, height_ratios=[1, 1.5],
                                          width_ratios=[40, 1])
        self.ax_cg_xy = self.fig_cg.add_subplot(gs_cg[0, 0])
        # Axes DÉDIÉ à la barre de couleur, créé une fois : la recréer à chaque
        # tracé l'empilait et rétrécissait le graphe un peu plus à chaque
        # calcul, et la détacher après un clear() de son parent lève une erreur.
        self.ax_cg_barre = self.fig_cg.add_subplot(gs_cg[0, 1])
        self.ax_cg_angle = self.fig_cg.add_subplot(gs_cg[1, :])
        self.canvas_cg = FigureCanvas(self.fig_cg)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_cg), "Orbite du CG")

        self.fig_thermique = Figure(figsize=(6, 5))
        self.ax_thermique = self.fig_thermique.add_subplot(111)
        self.canvas_thermique = FigureCanvas(self.fig_thermique)
        self.tabs.addTab(self._onglet_avec_barre(self.canvas_thermique), "Thermique")

        self.tabs.addTab(self._onglet_animation(), "Animation")

        return self.tabs

    def _onglet_animation(self) -> QWidget:
        """Lecteur interactif du déploiement du spiral -- calcule une fois
        (bouton "Calculer l'animation" au panneau 5), puis se pilote comme un
        lecteur vidéo : curseur de position, lecture/pause en va-et-vient,
        vitesse réglable, pas-à-pas manuel avant/arrière. L'export GIF
        (panneau 5) réutilise les mêmes frames sans recalcul."""
        self._anim_angles = None
        self._anim_frames = None
        self._anim_sequence = []
        self._anim_pos = 0  # position dans _anim_sequence (pas dans frames directement)
        self._anim_timer = QTimer(self)
        self._anim_timer.timeout.connect(self._anim_avancer_auto)

        conteneur = QWidget()
        v = QVBoxLayout(conteneur)
        v.setContentsMargins(0, 0, 0, 0)

        self.fig_anim = Figure(figsize=(6, 6))
        self.ax_anim = self.fig_anim.add_subplot(111)
        self.canvas_anim = FigureCanvas(self.fig_anim)
        barre_anim = NavigationToolbar(self.canvas_anim, conteneur)
        navigation.activer(self.canvas_anim)
        v.addWidget(barre_anim, 0)
        v.addWidget(self.canvas_anim, 1)

        ligne_slider = QHBoxLayout()
        self.slider_anim = QSlider(Qt.Horizontal)
        self.slider_anim.setEnabled(False)
        self.slider_anim.valueChanged.connect(self._anim_afficher_position)
        ligne_slider.addWidget(self.slider_anim, 1)
        self.label_anim_angle = QLabel("angle = --")
        self.label_anim_angle.setMinimumWidth(110)
        ligne_slider.addWidget(self.label_anim_angle)
        v.addLayout(ligne_slider, 0)

        ligne_boutons = QHBoxLayout()
        self.btn_anim_prev = QPushButton("◀◀ Précédent")
        self.btn_anim_prev.setToolTip("Recule d'une frame (pas-à-pas manuel).")
        self.btn_anim_prev.setEnabled(False)
        self.btn_anim_prev.clicked.connect(lambda: self._anim_pas_manuel(-1))
        ligne_boutons.addWidget(self.btn_anim_prev)

        self.btn_anim_play = QPushButton("▶ Lecture")
        self.btn_anim_play.setToolTip("Lecture automatique en va-et-vient (comme l'oscillation réelle du balancier).")
        self.btn_anim_play.setEnabled(False)
        self.btn_anim_play.clicked.connect(self._anim_basculer_lecture)
        ligne_boutons.addWidget(self.btn_anim_play)

        self.btn_anim_next = QPushButton("Suivant ▶▶")
        self.btn_anim_next.setToolTip("Avance d'une frame (pas-à-pas manuel).")
        self.btn_anim_next.setEnabled(False)
        self.btn_anim_next.clicked.connect(lambda: self._anim_pas_manuel(+1))
        ligne_boutons.addWidget(self.btn_anim_next)

        ligne_boutons.addSpacing(20)
        ligne_boutons.addWidget(QLabel("Vitesse (im/s) :"))
        self.sb_anim_fps = _intbox(1, 30, 12)
        self.sb_anim_fps.setToolTip("Images par seconde en lecture automatique.")
        self.sb_anim_fps.valueChanged.connect(self._anim_mettre_a_jour_intervalle)
        ligne_boutons.addWidget(self.sb_anim_fps)
        v.addLayout(ligne_boutons, 0)

        return conteneur

    # --- 1. Origine -----------------------------------------------------
    def _panneau_origine(self) -> QGroupBox:
        box = QGroupBox("1. Origine du spiral")
        layout = QVBoxLayout(box)

        ligne_dxf = QHBoxLayout()
        btn_dxf = QPushButton("Charger un DXF...")
        btn_dxf.setToolTip("Importe un spiral déjà dessiné (DXF machine) : extraction de la fibre neutre et de l'épaisseur locale par ray-casting.")
        btn_dxf.clicked.connect(self.on_charger_dxf)
        ligne_dxf.addWidget(btn_dxf)
        self.label_dxf = QLabel("(aucun)")
        ligne_dxf.addWidget(self.label_dxf, 1)
        layout.addLayout(ligne_dxf)

        # Charger une configuration est une ORIGINE au même titre qu'un DXF :
        # sa place est en tête de l'onglet, à côté du chargement de dessin, et
        # non en bas parmi les sauvegardes.
        btn_cfg = QPushButton("Charger une configuration (.json)...")
        btn_cfg.setToolTip(
            "Recharge une configuration complète — spec du mouvement, repère, wafer et "
            "vecteur de conception — et régénère la géométrie correspondante.")
        btn_cfg.clicked.connect(self.on_charger_configuration)
        layout.addWidget(btn_cfg)

        form = QFormLayout()
        self.sb_r_min = _spinbox(0.05, 5.0, self.spec.r_min, 0.01, 4)
        self.sb_r_min.setToolTip("Rayon de la virole/collet (mm) -- point d'attache central, solidaire de l'axe du balancier.")
        self.sb_r_max = _spinbox(0.1, 10.0, self.spec.r_max, 0.01, 4)
        self.sb_r_max.setToolTip("Rayon maximal (mm) au piton -- limite d'encombrement du spiral (ne doit pas toucher le pont/coq).")
        self.sb_n_tours = _spinbox(1.0, 30.0, self.spec.n_tours_reels, 0.5, 2)
        self.sb_n_tours.setToolTip("Nombre de tours RÉELS du spiral (virole -> piton, incluant les courbes Grossmann et Phillips).")
        self.sb_n_tours_fictifs = _spinbox(1.0, 40.0, self.spec.n_tours_fictifs, 0.5, 2)
        self.sb_n_tours_fictifs.setToolTip(
            "Nombre de tours \"fictifs\" servant à définir le pas théorique de référence "
            "= (r_max - r_min) / (n_tours_fictifs x 360°). Sert UNIQUEMENT à fixer la pente "
            "imposée à la jonction avec les courbes Grossmann/Phillips (continuité C1) -- "
            "ce n'est pas un nombre de tours physiquement présent sur le spiral."
        )
        self.sb_hauteur = _spinbox(0.01, 1.0, self.spec.hauteur_h, 0.01, 3)
        self.sb_hauteur.setToolTip("Hauteur hors-plan de la lame (mm), constante sur tout le spiral (usinage DRIE).")
        self.sb_freq = _spinbox(1.0, 10.0, self.spec.freq_cible, 0.1, 3)
        self.sb_freq.setToolTip("Fréquence chronométrique cible (Hz). L'épaisseur est automatiquement mise à l'échelle pour l'atteindre.")
        self.sb_inertie = _spinbox(
            _inertie_affichee(1e-11), _inertie_affichee(1e-7),
            _inertie_affichee(self.spec.inertie_balancier), 0.1, 3)
        self.sb_inertie.setToolTip(
            "Moment d'inertie du balancier, en mg·cm² — utilisé avec freq_cible pour "
            "calculer la raideur, donc l'épaisseur requise.\n\n"
            "Unité d'affichage seulement : le calcul, la sauvegarde JSON et le rapport "
            "restent en kg·m² (SI). 8 mg·cm² = 8·10⁻¹⁰ kg·m².\n\n"
            "La résolution du champ (10⁻³ mg·cm²) est celle de l'ancien champ en kg·m² : "
            "rien n'est perdu à la conversion."
        )
        self.sb_jeu_min = _spinbox(0.001, 0.1, self.spec.jeu_min_spires, 0.001, 4)
        self.sb_jeu_min.setToolTip("Jeu minimal requis entre spires adjacentes (mm), sous déformation à amplitude max. Contrainte dure de l'Optimisation (pas du Sondage).")
        self.sb_amplitude_max = _spinbox(100.0, 360.0, self.spec.amplitude_max, 5.0, 1)
        self.sb_amplitude_max.setToolTip("Amplitude d'oscillation maximale du balancier (deg) -- pire cas utilisé pour le jeu inter-spires, le GIF et la vue \"jeu minimal\".")
        self.ed_projet = QLineEdit("K1")
        self.ed_projet.setToolTip(
            "Nom du calibre ou de l'étude. Il préfixe TOUS les fichiers enregistrés, "
            "et nomme leur dossier dans Results/ : « K2_20261010-1432 ».\n\n"
            "Le nom d'abord, l'horodatage ensuite : les études d'un même calibre se "
            "regroupent, et se trient chronologiquement à l'intérieur. L'ordre inverse "
            "mêlerait tous les calibres."
        )
        form.addRow("nom du projet", self.ed_projet)
        form.addRow("r_min (mm)", self.sb_r_min)
        form.addRow("r_max (mm)", self.sb_r_max)
        form.addRow("n_tours_reels", self.sb_n_tours)
        form.addRow("n_tours_fictifs", self.sb_n_tours_fictifs)
        form.addRow("hauteur_h (mm)", self.sb_hauteur)
        form.addRow("freq_cible (Hz)", self.sb_freq)
        form.addRow("inertie_balancier (mg.cm²)", self.sb_inertie)

        self.cb_plan = QComboBox()
        for _code, _texte in (("111", "(111) — isotrope dans le plan, 169 GPa"),
                              ("001", "(001) — anisotrope, 130 à 169 GPa"),
                              ("110", "(110) — anisotrope, 130 à 187 GPa")):
            self.cb_plan.addItem(_texte, _code)
        self.cb_plan.setToolTip(
            "Plan cristallin du wafer dans lequel le spiral est gravé.\n\n"
            "Le silicium n'a un module de Young ISOTROPE DANS LE PLAN que sur un wafer "
            "(111), où il vaut 169 GPa quelle que soit la direction. Sur un (001) il "
            "varie de 130 à 169 GPa avec une période de 90°, et sur un (110) de 130 à "
            "187 GPa.\n\n"
            "La lame d'un spiral balaie forcément toutes les directions du plan : sur un "
            "wafer anisotrope, son module change donc tout au long du ruban. En moyenne "
            "harmonique, un (001) donne 147 GPa, soit 13 % de raideur en moins qu'un "
            "(111) — ce n'est pas une approximation qu'on peut négliger."
        )
        self.cb_plan.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.cb_plan.setMinimumContentsLength(10)
        self.cb_plan.currentIndexChanged.connect(self._maj_orientation_cristal)
        form.addRow("plan cristallin du wafer", self.cb_plan)

        self.sb_orientation_cristal = _spinbox(-180.0, 180.0, 0.0, 5.0, 1)
        self.sb_orientation_cristal.setToolTip(
            "Angle, dans le repère du dessin, de l'axe cristallin de référence du plan "
            "— ⟨100⟩ pour un (001), ⟨1-10⟩ pour un (110).\n\n"
            "Il décide de l'orientation des directions molles par rapport au spiral. "
            "Sans effet sur un wafer (111), isotrope dans le plan."
        )
        form.addRow("orientation du cristal (deg)", self.sb_orientation_cristal)
        self._maj_orientation_cristal()   # grise le champ si le plan est isotrope
        form.addRow("jeu_min_spires (mm)", self.sb_jeu_min)
        form.addRow("amplitude_max (deg)", self.sb_amplitude_max)
        layout.addLayout(form)

        self.chk_cote_rouage = QCheckBox("Vue côté rouage du mouvement")
        self.chk_cote_rouage.setChecked(True)
        self.chk_cote_rouage.setToolTip(
            "Côté depuis lequel le spiral est vu, et depuis lequel les quatre positions "
            "verticales sont nommées.\n\n"
            "Le spiral est sous le coq, côté rouage : c'est de ce côté qu'on le voit et "
            "qu'il est dessiné. Les heures 12-3-6-9 y tournent dans le sens ANTIHORAIRE, "
            "comme tout cadran vu de dos. Décoché, on se place côté cadran et le sens "
            "s'inverse.\n\n"
            "Le tracé est orienté 6H EN HAUT -- l'attitude de la position COSC 6H, celle "
            "qui sert au DVH. La tige, qui sort à 3H, est donc à DROITE.\n\n"
            "Se tromper de côté ÉCHANGE 3H et 9H. La dispersion est inchangée (c'est une "
            "permutation) et l'optimisation aussi, mais les valeurs seraient attribuées "
            "aux mauvaises positions dès qu'on les compare à un relevé de "
            "chronocomparateur."
        )
        self.chk_cote_rouage.stateChanged.connect(self._tracer_geometrie)
        layout.addWidget(self.chk_cote_rouage)

        form_repere = QFormLayout()
        form_repere.setRowWrapPolicy(QFormLayout.WrapLongRows)
        self.sb_angle_echappement = _spinbox(-180.0, 180.0,
                                              getattr(self.spec, "angle_echappement_deg", -6.2),
                                              0.1, 2)
        self.sb_angle_echappement.setToolTip(
            "Angle que fait la LIGNE DES CENTRES — celle qui joint l'axe du balancier à "
            "celui de l'ancre, dite aussi ligne droite de l'échappement — avec l'axe de "
            "la tige de remontoir, comptée dans le sens trigonométrique TEL QU'ON LE "
            "VOIT côté rouage.\n\n"
            "La tige sort à 3H, donc à droite du tracé (6H en haut) : à −6.2°, la ligne "
            "des centres pointe juste sous l'horizontale, vers la droite. C'est sur elle "
            "que le spiral se fixe au collet ; "
            "le nombre de tours étant entier, c'est aussi sur elle que retombe le point "
            "de comptage (piton ou goupilles), à l'opposé.\n\n"
            "Cet angle ORIENTE tout le spiral dans le mouvement. Il ne change pas le "
            "MODULE de la sensibilité à la gravité (une rotation ne modifie pas la forme), "
            "mais il change sa PHASE — donc le DVH, qui ne regarde que la position 6H."
        )
        self.sb_angle_echappement.valueChanged.connect(self._tracer_geometrie)
        form_repere.addRow("angle échappement (deg)", self.sb_angle_echappement)

        self.chk_sens_trigo = QCheckBox("Développement trigonométrique (vu côté rouage)")
        self.chk_sens_trigo.setChecked(getattr(self.spec, "sens_trigonometrique", True))
        self.chk_sens_trigo.setToolTip(
            "Sens dans lequel le spiral se développe depuis le collet, vu côté rouage. "
            "Coché : sens antihoraire. Décoché : sens horaire — le spiral est le MIROIR "
            "du précédent.\n\n"
            "Ce n'est pas une convention d'affichage : le miroir change la phase du "
            "vecteur de sensibilité à la gravité, donc le DVH."
        )
        self.chk_sens_trigo.stateChanged.connect(self._tracer_geometrie)

        self.chk_collet_cote_ancre = QCheckBox("Collet du côté de l'ancre")
        self.chk_collet_cote_ancre.setChecked(
            getattr(self.spec, "collet_cote_ancre", True))
        self.chk_collet_cote_ancre.setToolTip(
            "LEQUEL DES DEUX BOUTS de la ligne des centres porte l'attache au collet.\n\n"
            "Une droite n'est définie qu'à 180° près par son angle : l'angle ci-dessus "
            "ne suffit donc pas à orienter le spiral dans le mouvement.\n\n"
            "Coché (défaut) : le bout qui pointe vers l'ancre, à l'opposé de la tige — "
            "le collet part alors à l'angle saisi ci-dessus (−6.2° par défaut), et le "
            "point de comptage retombe à l'opposé (173.8°).\n"
            "Décoché : l'inverse, tout le spiral fait un demi-tour.\n\n"
            "Se tromper de bout ne change NI la dispersion NI le module de sensibilité à "
            "la gravité, mais inverse le SIGNE du DVH : l'erreur est muette sur les "
            "graphes et ne se voit qu'à la comparaison avec un relevé réel."
        )
        self.chk_collet_cote_ancre.stateChanged.connect(self._tracer_geometrie)
        for signal in (self.sb_angle_echappement.valueChanged, self.chk_sens_trigo.stateChanged,
                       self.chk_collet_cote_ancre.stateChanged, self.cb_plan.currentIndexChanged,
                       self.sb_orientation_cristal.valueChanged):
            signal.connect(self._reposer_ancienne_experience)

        layout.addLayout(form_repere)
        layout.addWidget(self.chk_sens_trigo)
        layout.addWidget(self.chk_collet_cote_ancre)

        ligne_btn = QHBoxLayout()
        btn_init = QPushButton("Initialiser Archimède")
        btn_init.setToolTip("Construit un spiral neuf à partir de la spec ci-dessus + des paramètres de forme du panneau 2.")
        btn_init.clicked.connect(self.on_init_archimede)
        ligne_btn.addWidget(btn_init)
        layout.addLayout(ligne_btn)

        return box

    # --- 2. Paramètres ----------------------------------------------------
    def _panneau_parametres(self) -> QGroupBox:
        box = QGroupBox("2. Paramètres de forme")
        box.setToolTip(
            "Forme des deux terminaisons — courbe de Grossmann côté collet, courbe de "
            "Phillips ou terminaison à raquette côté piton — et finesse des profils "
            "libres d'épaisseur et d'écart radial. Ces valeurs démarrent neutres : ce "
            "sont le Sondage puis l'Optimisation qui les font évoluer."
        )
        layout = QVBoxLayout(box)
        form = QFormLayout()

        self.sb_gross_angle = _spinbox(5.0, 300.0, 90.0, 5.0, 1)
        self.sb_gross_angle.setToolTip(
            "Étendue angulaire de la courbe de Grossmann (deg), depuis le collet. "
            "IMPORTANT : avec des offsets à 0 (valeur de départ), cette courbe reste "
            "très proche d'une simple continuation du pas théorique -- changer l'angle "
            "seul ne produit donc qu'un effet discret, visible surtout dans l'onglet "
            "Géométrie > zoom \"collet\". Le vrai levier de forme, ce sont les offsets, "
            "réglés par le Sondage/l'Optimisation, pas directement ici."
        )
        self.sb_gross_n_pas = _intbox(10, 2000, 120)
        self.sb_gross_n_pas.setToolTip("Nombre d'éléments discrets du maillage physique réservés à cette zone (densité locale, pas une variable de forme).")
        self.sb_gross_n_points = _intbox(4, 20, 4)
        self.sb_gross_n_points.setToolTip("Nombre de points de contrôle de la courbe de Bézier (>=4). Plus il y en a, plus la forme explorable par l'optimisation est riche (inflexions possibles).")
        self.sb_phil_angle = _spinbox(5.0, 300.0, 120.0, 5.0, 1)
        self.sb_phil_angle.setToolTip(
            "Étendue angulaire de la courbe de Phillips (deg), avant le piton. "
            "Même remarque que Grossmann : avec des offsets à 0, l'effet visuel est "
            "discret (voir onglet Géométrie > zoom \"piton\")."
        )
        self.sb_phil_n_pas = _intbox(10, 2000, 150)
        self.sb_phil_n_pas.setToolTip("Nombre d'éléments discrets du maillage physique réservés à cette zone (densité locale, pas une variable de forme).")
        self.sb_phil_n_points = _intbox(4, 20, 4)
        self.sb_phil_n_points.setToolTip("Nombre de points de contrôle de la courbe de Bézier (>=4). Plus il y en a, plus la forme explorable par l'optimisation est riche (inflexions possibles).")
        self.sb_n_corps = _intbox(50, 10000, 1400)
        self.sb_n_corps.setToolTip(
            "Nombre d'éléments discrets du CORPS seul (spirale à pas constant, entre les "
            "deux courbes terminales).\n\n"
            "Ce n'est PAS le maillage total : les zones terminales ont leur propre densité "
            "(Grossmann n_pas, Phillips ou raquette n_pas). Le total est affiché sous le "
            "bouton « Générer géométrie » et c'est lui que rapporte la vérification à "
            "haute résolution."
        )
        self.sb_n_noeuds_ep = _intbox(2, 60, 20)
        self.sb_n_noeuds_ep.setToolTip("Nombre de points de contrôle du profil d'épaisseur libre le long du spiral (plus il y en a, plus le profil peut varier localement).")
        self.sb_n_points_dxf = _intbox(0, 50000, 0)
        self.sb_n_points_dxf.setToolTip(
            "Nombre de points par flanc dans les fichiers DXF exportés. 0 = reprendre le "
            "maillage de calcul tel quel.\n\n"
            "Les deux n'ont aucune raison d'être égaux : le maillage de calcul est un "
            "compromis entre précision et temps, et il densifie volontairement les zones "
            "terminales. Un fichier de gravure, lui, doit seulement décrire le contour "
            "assez finement pour la machine.\n\n"
            "Ré-échantillonner ne change pas la courbe, seulement sa description. Le point "
            "de comptage reste un point du maillage."
        )
        self.sb_n_noeuds_dr = _intbox(0, 40, 12)
        self.sb_n_noeuds_dr.setToolTip("Nombre de points de contrôle de l'écart radial libre par rapport à la spirale de référence (0 = rayon strictement imposé par le pas théorique).")

        form.addRow("Grossmann -- angle (deg)", self.sb_gross_angle)
        form.addRow("Grossmann -- n_pas (maillage)", self.sb_gross_n_pas)
        form.addRow("Grossmann -- n_points (Bézier)", self.sb_gross_n_points)
        form.addRow("Phillips -- angle (deg)", self.sb_phil_angle)
        form.addRow("Phillips -- n_pas (maillage)", self.sb_phil_n_pas)
        form.addRow("Phillips -- n_points (Bézier)", self.sb_phil_n_points)
        form.addRow("Corps -- n_elements", self.sb_n_corps)
        form.addRow("Profil épaisseur -- n_noeuds", self.sb_n_noeuds_ep)
        form.addRow("Profil écart radial -- n_noeuds", self.sb_n_noeuds_dr)
        form.addRow("Points par flanc à l'export DXF", self.sb_n_points_dxf)
        layout.addLayout(form)

        layout.addWidget(self._panneau_raquette())
        layout.addWidget(self._panneau_thermique())

        btn = QPushButton("Générer géométrie")
        btn.setToolTip(
            "Les valeurs de forme (offsets, épaisseurs) démarrent neutres (raccord "
            "rectiligne / épaisseur théorique uniforme) : c'est le Sondage puis "
            "l'Optimisation (panneaux 3-4 ci-dessous) qui les font évoluer. Changer "
            "seulement l'angle ici ne produit qu'un effet visuel discret -- voir les "
            "zooms dédiés dans l'onglet Géométrie."
        )
        btn.clicked.connect(self.on_generer_geometrie)
        layout.addWidget(btn)

        self.label_diag = QLabel("")
        self.label_diag.setStyleSheet("font-family: monospace;")
        self.label_diag.setWordWrap(True)
        layout.addWidget(self.label_diag)

        return box

    def _panneau_raquette(self) -> QGroupBox:
        """Terminaison extérieure à raquette, en remplacement de la courbe de
        Phillips. Placé dans le panneau de FORME (et non dans celui de
        l'optimisation, comme l'échappement) parce qu'il change réellement la
        géométrie construite, pas seulement le critère optimisé."""
        box = QGroupBox("Terminaison à raquette")
        box.setToolTip(
            "Pour poser une raquette, la portion du spiral balayée par les goupilles "
            "doit être un ARC CENTRÉ SUR L'AXE du balancier : les goupilles y glissent "
            "sans pousser le spiral radialement, donc seule la LONGUEUR ACTIVE change. "
            "Une courbe de Phillips (ou une spirale d'Archimède) à cet endroit serait "
            "déformée par le déplacement de la raquette.\n\n"
            "Activer cette option remplace la courbe de Phillips par : une courbe en S "
            "(tangente au corps d'un côté, à l'arc de l'autre) suivie de l'arc "
            "concentrique au rayon r_max. Le point de comptage (goupilles) se place au "
            "milieu de l'arc."
        )
        v = QVBoxLayout(box)

        self.chk_raquette = QCheckBox("Spiral à raquette (remplace la courbe de Phillips)")
        self.chk_raquette.setToolTip(
            "Décoché : spiral LIBRE, terminaison par courbe de Phillips (comportement "
            "d'origine). Coché : terminaison courbe en S + arc concentrique."
        )
        self.chk_raquette.stateChanged.connect(self._maj_etat_raquette)
        v.addWidget(self.chk_raquette)

        form = QFormLayout()
        self.sb_raq_plage = _spinbox(1.0, 1000.0, 60.0, 5.0, 1)
        self.sb_raq_plage.setToolTip(
            "Plage de réglage exigée de la raquette, DE PART ET D'AUTRE du nominal "
            "(saisir 60 donne ±60 s/j). C'est cette exigence qui DÉTERMINE la longueur "
            "de l'arc : le programme la calcule, vous ne la saisissez pas. La "
            "sensibilité varie comme 1/épaisseur³ au piton, donc un piton aminci "
            "raccourcit beaucoup l'arc nécessaire. La plage réellement obtenue est "
            "vérifiée et affichée dans les résultats numériques."
        )
        self.sb_raq_trans_angle = _spinbox(5.0, 300.0, 120.0, 5.0, 1)
        self.sb_raq_trans_angle.setToolTip(
            "Étendue angulaire (deg) laissée au spiral pour passer du corps à l'arc. "
            "Plus elle est courte, plus la courbe en S doit se redresser vite : sur un "
            "calibre de référence, la pente exigée passe de 117 µm/rad (transition sur "
            "180°) à 700 µm/rad (sur 30°), contre 25 µm/rad pour le pas d'Archimède. "
            "L'inflexion en S est donc automatique quel que soit l'angle ; c'est la "
            "douceur du raccord qui est en jeu."
        )
        self.sb_raq_n_points = _intbox(4, 20, 4)
        self.sb_raq_n_points.setToolTip(
            "Points de contrôle de la courbe en S (>=4). À 4 points la courbe est "
            "entièrement déterminée par ses deux tangences (Hermite cubique) : aucune "
            "variable d'optimisation. Au-delà, chaque point supplémentaire ajoute un "
            "degré de liberté de forme sans toucher aux tangences. Pour comparer à "
            "armes égales avec une courbe de Phillips à N points, prendre N+1."
        )
        self.sb_raq_trans_n_pas = _intbox(10, 2000, 150)
        self.sb_raq_trans_n_pas.setToolTip("Nombre d'éléments discrets du maillage réservés à la courbe en S.")
        self.sb_raq_arc_avant = _spinbox(0.0, 270.0, 30.0, 5.0, 1)
        self.sb_raq_arc_avant.setToolTip(
            "Étendue angulaire (deg) de l'arc concentrique AVANT le point de comptage, "
            "c'est-à-dire entre la fin du raccord en S et les goupilles.\n\n"
            "Prendre la même valeur que l'angle du piton au-delà place le point de "
            "comptage au MILIEU de la partie circulaire terminale.\n\n"
            "Cette portion fait partie du spiral ACTIF : l'allonger remplace de la "
            "spirale d'Archimède par de l'arc concentrique dans la partie qui travaille. "
            "Elle doit aussi être au moins aussi longue que la demi-course des goupilles, "
            "sans quoi celles-ci quitteraient l'arc."
        )
        self.sb_raq_n_pas_arc = _intbox(1, 2000, 100)
        self.sb_raq_n_pas_arc.setToolTip(
            "Nombre d'éléments discrets du maillage sur l'arc avant le point de comptage."
        )
        form.addRow("Plage de réglage ± (s/j)", self.sb_raq_plage)
        form.addRow("Transition -- angle (deg)", self.sb_raq_trans_angle)
        form.addRow("Transition -- n_points (Bézier)", self.sb_raq_n_points)
        form.addRow("Transition -- n_pas (maillage)", self.sb_raq_trans_n_pas)
        self.sb_raq_piton_angle = _spinbox(0.0, 300.0, 30.0, 5.0, 1)
        self.sb_raq_piton_angle.setToolTip(
            "Étendue angulaire de l'arc APRÈS le point de comptage, jusqu'au piton.\n\n"
            "Cette portion est tenue à ses deux bouts : elle ne travaille pas et ne se "
            "déforme pas. Elle occupe néanmoins de la place (jeu entre spires) et doit "
            "être usinée.\n\n"
            "Lui donner la même valeur qu'à l'arc avant place le point de comptage au "
            "MILIEU de la portion circulaire, ce qui laisse à la raquette la même course "
            "dans les deux sens. Une lame trop courte au-delà des goupilles prive la "
            "raquette d'autorité : l'excursion latérale varie comme le CARRÉ de cette "
            "longueur (0,9 µm pour 124 µm de lame, 175 µm pour 1856)."
        )
        self.sb_raq_piton_n_pas = _intbox(0, 2000, 100)
        self.sb_raq_piton_n_pas.setToolTip(
            "Nombre d'éléments discrets du maillage sur l'arc après le point de comptage."
        )
        self.sb_raq_jeu = _spinbox(0.0, 200.0, 20.0, 1.0, 1)
        self.sb_raq_jeu.setToolTip(
            "Jeu LIBRE entre les deux goupilles de la raquette, au-delà de l'épaisseur "
            "de la lame (µm).\n\n"
            "Tant que la lame est plaquée contre une goupille, le spiral est encastré "
            "là ; dès qu'elle décolle et traverse le jeu, il se comporte comme encastré "
            "au piton, donc plus long et plus souple. La montre retarde pendant ces "
            "traversées, d'autant plus que l'amplitude est faible.\n\n"
            "Zéro redonne EXACTEMENT le montage à piton fixe."
        )
        form.addRow("Arc avant comptage -- angle (deg)", self.sb_raq_arc_avant)
        form.addRow("Arc avant comptage -- n_pas", self.sb_raq_n_pas_arc)
        form.addRow("Arc comptage → piton -- angle (deg)", self.sb_raq_piton_angle)
        form.addRow("Arc comptage → piton -- n_pas", self.sb_raq_piton_n_pas)
        form.addRow("Jeu entre goupilles (µm)", self.sb_raq_jeu)
        v.addLayout(form)
        self._champs_raquette = [self.sb_raq_plage, self.sb_raq_trans_angle,
                                  self.sb_raq_n_points, self.sb_raq_trans_n_pas,
                                  self.sb_raq_n_pas_arc, self.sb_raq_arc_avant,
                                  self.sb_raq_piton_angle, self.sb_raq_piton_n_pas,
                                  self.sb_raq_jeu]
        self._champs_phillips = [self.sb_phil_angle, self.sb_phil_n_pas, self.sb_phil_n_points]
        self._maj_etat_raquette()
        return box

    def _maj_etat_raquette(self):
        """Les deux terminaisons s'excluent : on grise celle qui ne sert pas,
        pour qu'il soit visible d'un coup d'oeil laquelle pilote la géométrie."""
        actif = self.chk_raquette.isChecked()
        for w in getattr(self, "_champs_raquette", []):
            w.setEnabled(actif)
        for w in getattr(self, "_champs_phillips", []):
            w.setEnabled(not actif)

    def _panneau_thermique(self) -> QGroupBox:
        """Compensation thermique par couche d'oxyde de silicium."""
        box = QGroupBox("Compensation thermique")
        box.setToolTip(
            "Le silicium s'assouplit fortement en chauffant : un spiral nu retarde "
            "d'environ 3.8 s/j/°C, soit 115 s/j entre 8 et 38 °C -- six fois le "
            "critère COSC. L'oxyde thermique a le comportement inverse (+215 ppm/K "
            "contre -60) et, placé à la périphérie de la section où il travaille le "
            "plus en flexion, quelques microns suffisent à annuler l'effet.\n\n"
            "L'épaisseur d'oxyde est RECALCULÉE pour chaque candidat : elle dépend du "
            "profil d'épaisseur, que l'optimiseur fait varier. La lame gravée doit "
            "donc être plus mince que la lame finie -- l'export DXF produit les deux "
            "fichiers."
        )
        v = QVBoxLayout(box)

        self.chk_oxydation = QCheckBox("Lame oxydée (compensation thermique)")
        self.chk_oxydation.setToolTip(
            "Décoché : silicium nu, non compensé. L'analyse thermique reste "
            "disponible et montre ce que coûte l'absence de compensation."
        )
        self.chk_oxydation.stateChanged.connect(self._maj_etat_thermique)
        v.addWidget(self.chk_oxydation)

        form = QFormLayout()
        self.sb_ox_impose = _spinbox(0.0, 30.0, 0.0, 0.1, 3)
        self.sb_ox_impose.setToolTip(
            "Épaisseur d'oxyde imposée (µm). 0 = calculée automatiquement pour annuler "
            "la variation de marche sur la plage visée. Imposer une valeur sert à "
            "étudier ce que donne une épaisseur que le procédé sait effectivement "
            "tenir : la marche résiduelle est alors affichée."
        )
        self.sb_T_min = _spinbox(-40.0, 60.0, 8.0, 1.0, 1)
        self.sb_T_min.setToolTip("Borne basse de la plage de température visée (°C).")
        self.sb_T_max = _spinbox(-40.0, 120.0, 38.0, 1.0, 1)
        self.sb_T_max.setToolTip("Borne haute de la plage de température visée (°C).")
        self.sb_T_ref = _spinbox(-40.0, 60.0, 23.0, 1.0, 1)
        self.sb_T_ref.setToolTip(
            "Température de référence (°C) : celle à laquelle la marche est mise à zéro "
            "et à laquelle la fréquence est réglée. Seule la VARIATION avec la "
            "température a un sens ici."
        )
        self.sb_cte_bal = _spinbox(0.0, 40.0, 16.7, 0.1, 2)
        self.sb_cte_bal.setToolTip(
            "Dilatation du BALANCIER (ppm/K) -- 16.7 pour du CuBe2 (Glucydur). Elle "
            "compte autant que le spiral : son inertie croît comme le carré du rayon, "
            "ce qui pèse directement sur la fréquence. Une erreur de 1 ppm/K déplace "
            "l'épaisseur d'oxyde requise de 0.09 µm."
        )
        self.sb_tol_therm = _spinbox(0.001, 5.0, 0.1, 0.01, 3)
        self.sb_tol_therm.setToolTip(
            "Écart de marche thermique acceptable (s/j/°C). Repère : le critère COSC "
            "est 0.6 s/j/°C."
        )
        form.addRow("Épaisseur d'oxyde imposée (µm)", self.sb_ox_impose)
        form.addRow("Température min (°C)", self.sb_T_min)
        form.addRow("Température max (°C)", self.sb_T_max)
        form.addRow("Température de référence (°C)", self.sb_T_ref)
        form.addRow("Dilatation balancier (ppm/K)", self.sb_cte_bal)
        form.addRow("Tolérance (s/j/°C)", self.sb_tol_therm)

        self.sb_tce_si = _spinbox(-300.0, 0.0, -60.0, 1.0, 2)
        self.sb_tce_si.setToolTip("Coefficient thermoélastique du silicium (ppm/K), (1/E)·dE/dT. Valeur de littérature : -60.")
        self.sb_tce_ox = _spinbox(0.0, 500.0, 215.0, 5.0, 2)
        self.sb_tce_ox.setToolTip(
            "Coefficient thermoélastique de l'oxyde thermique (ppm/K). Valeur de "
            "littérature : +215. C'est LE coefficient le plus incertain et le plus "
            "influent : ±20 ppm/K déplacent l'épaisseur requise de 3.42 à 3.99 µm."
        )
        self.sb_E_ox = _spinbox(10.0, 200.0, 70.0, 1.0, 1)
        self.sb_E_ox.setToolTip("Module d'Young de l'oxyde thermique (GPa). Valeur de littérature : 70.")
        self.sb_tce2_si = _spinbox(-10.0, 10.0, 0.0, 0.01, 4)
        self.sb_tce2_si.setToolTip(
            "Terme du SECOND ordre du silicium (ppm/K²), (1/E)·d²E/dT². Nul par défaut "
            "faute de données. C'est une limite assumée : avec des modules strictement "
            "linéaires, le résidu après compensation ressort quasi nul, ce qui est "
            "certainement optimiste. C'est la courbure réelle de E(T) qui limite les "
            "chronomètres -- à renseigner dès que mesurée."
        )
        self.sb_tce2_ox = _spinbox(-10.0, 10.0, 0.0, 0.01, 4)
        self.sb_tce2_ox.setToolTip("Terme du second ordre de l'oxyde (ppm/K²). Nul par défaut, même remarque.")
        form.addRow("TCE silicium (ppm/K)", self.sb_tce_si)
        form.addRow("TCE oxyde (ppm/K)", self.sb_tce_ox)
        form.addRow("Module de l'oxyde (GPa)", self.sb_E_ox)
        form.addRow("TCE² silicium (ppm/K²)", self.sb_tce2_si)
        form.addRow("TCE² oxyde (ppm/K²)", self.sb_tce2_ox)
        v.addLayout(form)

        self.btn_thermique = QPushButton("Analyser en température")
        self.btn_thermique.setToolTip(
            "Calcule la marche en fonction de la température, l'épaisseur d'oxyde "
            "requise, et de combien le résultat bouge dans l'incertitude de chaque "
            "coefficient matériau. Voir l'onglet \"Thermique\"."
        )
        self.btn_thermique.clicked.connect(self.on_analyser_thermique)
        v.addWidget(self.btn_thermique)

        self._champs_thermiques = [self.sb_ox_impose]
        self._maj_etat_thermique()
        return box

    def _maj_etat_thermique(self):
        for w in getattr(self, "_champs_thermiques", []):
            w.setEnabled(self.chk_oxydation.isChecked())

    def _reglage_thermique(self) -> ReglageThermique:
        return ReglageThermique(
            temperature_min=self.sb_T_min.value(),
            temperature_max=self.sb_T_max.value(),
            temperature_reference=self.sb_T_ref.value(),
            cte_balancier=self.sb_cte_bal.value() * 1e-6,
            tolerance_s_par_j_par_degc=self.sb_tol_therm.value(),
            actif=self.chk_oxydation.isChecked(),
            tce_silicium=self.sb_tce_si.value() * 1e-6,
            tce_oxyde=self.sb_tce_ox.value() * 1e-6,
            e_oxyde=self.sb_E_ox.value() * 1e9,
            tce2_silicium=self.sb_tce2_si.value() * 1e-6,
            tce2_oxyde=self.sb_tce2_ox.value() * 1e-6,
        )

    # --- 3. Sondage ------------------------------------------------------
    def _panneau_sondage(self) -> QGroupBox:
        box = QGroupBox("3. Sondage")
        box.setToolTip(
            "Tirage aléatoire dans les bornes, classé du moins mauvais au pire. "
            "Le jeu minimal requis N'EST PAS appliqué ici comme un filtre -- c'est "
            "l'étape suivante (Optimisation) qui doit combler l'écart ; le sondage "
            "sert juste à démarrer d'un point correct plutôt que d'un point arbitraire."
        )
        layout = QVBoxLayout(box)
        form = QFormLayout()
        self.sb_sondage_n = _intbox(10, 5000, 200)
        self.sb_sondage_n.setToolTip("Nombre de tirages aléatoires à évaluer.")
        self.sb_sondage_seed = _intbox(0, 999999, 0)
        self.sb_sondage_seed.setToolTip("Graine aléatoire (reproductibilité) : la même graine redonne le même tirage.")
        self.sb_amp_min = _spinbox(50.0, 360.0, 200.0, 5.0, 0)
        self.sb_amp_min.setToolTip("Borne basse de la plage d'amplitude analysée pour la dispersion de marche.")
        self.sb_amp_max = _spinbox(50.0, 360.0, self.spec.amplitude_max, 5.0, 0)
        self.sb_amp_max.setToolTip(
            "Borne haute de la plage d'amplitude analysée. Utilisée à la fois pour "
            "l'affichage et pour l'objectif d'optimisation (grille interne fixe à 5°, "
            "voir note du panneau)."
        )
        self.sb_pas_balayage = _spinbox(5.0, 90.0, 20.0, 5.0, 0)
        self.sb_pas_balayage.setToolTip(
            "Pas angulaire (deg) du balayage elastica servant à construire l'orbite du "
            "centre de gravité (coûteux : un vrai calcul elastica par angle). Plus fin "
            "= plus précis mais plus lent."
        )
        form.addRow("Échantillons", self.sb_sondage_n)
        form.addRow("Graine (seed)", self.sb_sondage_seed)
        form.addRow("Amplitude analyse min (deg)", self.sb_amp_min)
        form.addRow("Amplitude analyse max (deg)", self.sb_amp_max)
        form.addRow("Pas balayage orbite (deg)", self.sb_pas_balayage)
        layout.addLayout(form)

        self.btn_sondage = QPushButton("Lancer le sondage")
        self.btn_sondage.setToolTip(
            "La dispersion de marche est toujours évaluée finement (pas de 5°) en "
            "interne, quelle que soit la résolution d'affichage -- évite qu'un pire cas "
            "passe entre deux points grossiers."
        )
        self.btn_sondage.clicked.connect(self.on_lancer_sondage)
        layout.addWidget(self.btn_sondage)

        self.label_sondage = QLabel("(pas encore lancé)")
        layout.addWidget(self.label_sondage)

        self.table_sondage = QTableWidget(0, 3)
        self.table_sondage.setHorizontalHeaderLabels(
            ["Dispersion (s/j)", "Jeu min (µm)", "Compens. (s/j)"])
        self.table_sondage.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_sondage.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_sondage.itemSelectionChanged.connect(self.on_selectionner_sondage)
        self.table_sondage.setMaximumHeight(160)
        layout.addWidget(self.table_sondage)

        return box

    # --- 4. Optimisation ---------------------------------------------------
    def _panneau_optimisation(self) -> QGroupBox:
        box = QGroupBox("4. Optimisation NSGA-II")
        layout = QVBoxLayout(box)
        form = QFormLayout()
        self.sb_pop = _intbox(4, 500, 60)
        self.sb_pop.setToolTip("Taille de la population NSGA-II. Plus grand = explore plus large par génération, mais plus lent.")
        self.sb_gen = _intbox(1, 500, 40)
        self.sb_gen.setToolTip("Nombre de générations NSGA-II. Le croisement génétique peut trouver des solutions faisables qu'un sondage aléatoire seul ne trouve pas -- privilégier plus de générations si le sondage ne donne rien de bon.")
        self.sb_ga_seed = _intbox(0, 999999, 42)
        self.sb_ga_seed.setToolTip("Graine aléatoire NSGA-II (reproductibilité).")
        form.addRow("Taille population", self.sb_pop)
        form.addRow("Nombre de générations", self.sb_gen)
        form.addRow("Graine (seed)", self.sb_ga_seed)
        layout.addLayout(form)


        self.btn_optimisation = QPushButton("Lancer l'optimisation")
        self.btn_optimisation.setToolTip(
            "Objectifs minimisés : dispersion de marche entre positions verticales, et "
            "jeu inter-spires (maximisé). Contraintes : jeu >= jeu_min_spires (panneau 1) "
            "et compensation de l'échappement dans sa tolérance (ci-dessus). Démarre "
            "depuis le meilleur point du Sondage si disponible, sinon depuis la "
            "géométrie active."
        )
        self.btn_optimisation.clicked.connect(self.on_lancer_optimisation)
        layout.addWidget(self.btn_optimisation)

        self.label_optim = QLabel("(pas encore lancée)")
        layout.addWidget(self.label_optim)

        # Une colonne de cases a cocher EN PREMIER : elle sert a choisir les
        # geometries a exporter, independamment de la selection courante. Sans
        # elle, exporter plusieurs solutions obligerait a les charger une par
        # une et a relancer un enregistrement a chaque fois.
        self.table_pareto = QTableWidget(0, 4)
        self.table_pareto.setToolTip(
            "Solutions non-dominées du front de Pareto. Cliquer une ligne la charge "
            "comme géométrie active ; cocher la case la retient pour l'export DXF.")
        self.table_pareto.setHorizontalHeaderLabels(
            ["Export", "Dispersion (s/j)", "Jeu min (µm)", "Compens. (s/j)"])
        self.table_pareto.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_pareto.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_pareto.itemSelectionChanged.connect(self.on_selectionner_solution)
        self.table_pareto.setMaximumHeight(160)
        layout.addWidget(self.table_pareto)

        return box

    def _panneau_echappement(self) -> QGroupBox:
        """Cible de compensation de l'erreur d'échappement.

        Ce réglage ne décrit PAS le spiral mais l'échappement auquel il sera
        associé : c'est une donnée d'entrée du réglage chronométrique, au même
        titre que l'inertie du balancier. Il est placé ici, dans le panneau
        d'optimisation, parce qu'il n'agit que là -- il ajoute une contrainte,
        il ne change pas la géométrie construite."""
        box = QGroupBox("Erreur d'échappement")
        box.setToolTip(
            "L'échappement fait RETARDER la montre aux faibles amplitudes, dans "
            "TOUTES les positions (horizontales comprises). Le déplacement du centre "
            "de gravité du spiral ne peut pas compenser ce retard : son effet est nul "
            "à plat. Seul l'ISOCHRONISME PROPRE du spiral (variation de sa période "
            "avec l'amplitude, indépendante de l'orientation) le peut. On demande donc "
            "ici au spiral d'AVANCER de la valeur ci-dessous à l'amplitude basse."
        )
        v = QVBoxLayout(box)
        form = QFormLayout()

        self.chk_echappement = QCheckBox("Compenser à l'optimisation")
        self.chk_echappement.setChecked(True)
        self.chk_echappement.setToolTip(
            "Décoché : la marche d'échappement reste calculée et affichée, mais "
            "l'optimiseur ne cherche plus à la compenser (une contrainte de moins)."
        )
        v.addWidget(self.chk_echappement)

        self.sb_ech_retard = _spinbox(0.0, 120.0, 8.0, 0.5, 2)
        self.sb_ech_retard.setToolTip(
            "Retard (s/j) induit par l'échappement à l'amplitude basse, par rapport à "
            "l'amplitude de référence. Valeur MESURÉE ou spécifiée pour l'échappement "
            "retenu, pas calculée ici. C'est exactement l'avance que le spiral doit "
            "produire de lui-même pour que le total soit plat."
        )
        self.sb_ech_amp_basse = _spinbox(50.0, 359.0, 200.0, 5.0, 0)
        self.sb_ech_amp_basse.setToolTip(
            "Amplitude basse de référence (deg), typiquement la fin de réserve de "
            "marche -- là où l'erreur d'échappement est la plus pénalisante. Doit "
            "rester dans la plage d'amplitude analysée (panneau 3)."
        )
        self.sb_ech_amp_ref = _spinbox(50.0, 360.0, 300.0, 5.0, 0)
        self.sb_ech_amp_ref.setToolTip(
            "Amplitude de référence (deg) à laquelle la marche est conventionnellement "
            "mise à zéro. Seule la VARIATION de marche avec l'amplitude a un sens ici : "
            "un décalage constant est rattrapé au réglage de la raquette."
        )
        self.sb_ech_tolerance = _spinbox(0.05, 20.0, ReglageEchappement().tolerance_s_par_j, 0.1, 2)
        self.sb_ech_tolerance.setToolTip(
            "Largeur de la bande acceptée autour de la cible (s/j). Une cible exacte "
            "serait inatteignable et rendrait toute la population infaisable ; la bande "
            "laisse l'optimiseur arbitrer entre compensation et dispersion de gravité. "
            "Ne pas descendre sous l'incertitude de maillage (bouton \"Vérifier à haute "
            "résolution\" au panneau 5)."
        )
        self.cb_variante = QComboBox()
        self.cb_variante.addItem("Milieu des deux montages", "milieu")
        self.cb_variante.addItem("Montage à raquette (jeu compris)", "raquette")
        self.cb_variante.addItem("Piton fixe au point de comptage", "piton")
        self.cb_variante.setToolTip(
            "Sur quel MONTAGE la compensation est réglée. Un même spiral ne peut pas "
            "être parfait dans les deux : le jeu de goupilles n'existe que dans l'un, "
            "et l'écart entre les deux est une PENTE en fonction de l'amplitude, pas "
            "un décalage constant -- la raquette elle-même ne peut donc pas le "
            "rattraper.\n\n"
            "Milieu : chaque montage est faux de la moitié de l'écart, symétriquement.\n"
            "Raquette : parfait à raquette, la version piton fixe avancera à basse "
            "amplitude.\n"
            "Piton fixe : l'inverse.\n\n"
            "Les trois réglages coïncident exactement si le jeu est nul ou si le "
            "spiral n'a pas de goupilles."
        )
        form.addRow("Montage visé", self.cb_variante)
        form.addRow("Retard échappement (s/j)", self.sb_ech_retard)
        form.addRow("à l'amplitude (deg)", self.sb_ech_amp_basse)
        form.addRow("Amplitude référence (deg)", self.sb_ech_amp_ref)
        form.addRow("Tolérance (s/j)", self.sb_ech_tolerance)
        v.addLayout(form)
        return box

    def _verifier_coherence_echappement(self) -> bool:
        """L'isochronisme propre se lit sur la courbe d'énergie du balayage,
        qui ne monte que jusqu'à l'amplitude d'analyse maximale : demander une
        compensation au-delà rendrait TOUS les candidats invalides, sans que
        l'optimiseur puisse dire pourquoi (il ne voit qu'une pénalité). On le
        signale ici, avant de lancer, plutôt que de laisser tourner un calcul
        long pour rien. Retourne False si l'utilisateur renonce."""
        reglage = self._reglage_echappement()
        amp_max = self.sb_amp_max.value()
        depassement = max(reglage.amplitude_basse_deg, reglage.amplitude_reference_deg)
        if not reglage.actif or depassement <= amp_max:
            return True
        reponse = QMessageBox.warning(
            self, "Compensation d'échappement",
            f"L'amplitude basse ({reglage.amplitude_basse_deg:.0f}°) ou de référence "
            f"({reglage.amplitude_reference_deg:.0f}°) dépasse l'amplitude d'analyse "
            f"maximale ({amp_max:.0f}°, panneau 3).\n\n"
            f"L'isochronisme propre ne peut pas y être évalué : tous les candidats "
            f"seront rejetés.\n\nAugmenter l'amplitude d'analyse maximale, ou "
            f"décocher la compensation.",
            QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Cancel,
        )
        return reponse == QMessageBox.Ok

    def _jeu_goupilles_um(self) -> float:
        """Jeu appliqué aux calculs. Nul si la géométrie active n'a pas de
        portion morte (spiral libre ou DXF importé) : il n'y aurait alors
        aucune goupille où l'appliquer."""
        if self.geom is None or not self.geom.a_partie_morte:
            return 0.0
        return self.sb_raq_jeu.value() if hasattr(self, "sb_raq_jeu") else 0.0

    def _reglage_echappement(self) -> ReglageEchappement:
        return ReglageEchappement(
            retard_s_par_j=self.sb_ech_retard.value(),
            amplitude_basse_deg=self.sb_ech_amp_basse.value(),
            amplitude_reference_deg=self.sb_ech_amp_ref.value(),
            tolerance_s_par_j=self.sb_ech_tolerance.value(),
            actif=self.chk_echappement.isChecked(),
            variante_visee=self.cb_variante.currentData(),
        )

    # --- 5. Résultats / export -------------------------------------------
    # ------------------------------------------------------------------
    # Outils d'analyse, construits une fois par onglet qui en a besoin
    # ------------------------------------------------------------------
    def _dire_ecourte(self, cle: str, texte: str, maxi: int = 34):
        """Comme _dire, mais raccourcit les mots trop longs (chemins) par le
        milieu et met le texte entier en info-bulle."""
        for e in self._affichages.get(cle, ()):
            _texte_ecourte(e, texte, maxi)

    def _infobulle(self, cle: str, texte: str):
        for e in self._affichages.get(cle, ()):
            e.setToolTip(texte)

    def _enregistrer_commande(self, cle: str, widget):
        """Mémorise un widget dupliqué entre onglets, pour que les mises à jour
        s'appliquent à toutes ses copies. Sans ce registre, un bouton désactivé
        pendant un calcul ne le serait que dans l'onglet visible, et le même
        calcul pourrait être relancé depuis l'autre."""
        self._commandes.setdefault(cle, []).append(widget)
        return widget

    def _bouton(self, cle: str, texte: str, gestionnaire, infobulle: str = ""):
        b = QPushButton(texte)
        if infobulle:
            b.setToolTip(infobulle)
        b.clicked.connect(gestionnaire)
        return self._enregistrer_commande(cle, b)

    def _affichage(self, cle: str, hauteur_min: int = 0, monospace: bool = False):
        e = QLabel("")
        e.setWordWrap(True)
        if hauteur_min:
            e.setMinimumHeight(hauteur_min)
        if monospace:
            e.setStyleSheet("font-family: monospace;")
        self._affichages.setdefault(cle, []).append(e)
        return e

    def _dire(self, cle: str, texte: str, infobulle: str = None):
        """Écrit dans toutes les copies d'un affichage."""
        for e in self._affichages.get(cle, ()):
            e.setText(texte)
            if infobulle is not None:
                e.setToolTip(infobulle)

    def _occuper(self, cle: str, actif: bool, texte_actif: str = None):
        """Même rôle que _definir_bouton_actif, mais sur toutes les copies d'une
        commande : un calcul lancé depuis l'onglet Calculs doit aussi griser le
        bouton jumeau de l'onglet Optimisation, sinon on peut le relancer en
        parallèle et écraser le thread en cours."""
        for b in self._commandes.get(cle, ()):
            self._definir_bouton_actif(b, actif, texte_actif)

    def _boite_analyse(self, suffixe: str) -> QGroupBox:
        """Commandes d'analyse de la GÉOMÉTRIE ACTIVE, quelle que soit son
        origine : DXF importé, forme générée, ou solution retenue sur le front
        de Pareto. Construite dans l'onglet « Calculs » ET dans l'onglet
        « Optimisation » : comparer deux solutions ne doit pas obliger à faire
        l'aller-retour entre les deux."""
        box = QGroupBox("Analyse de la géométrie active")
        box.setToolTip(
            "Ces commandes portent toujours sur la géométrie ACTIVE — la dernière "
            "générée, importée ou sélectionnée dans une table de solutions. Les tracés "
            "correspondants apparaissent dans les onglets de droite."
        )
        layout = QVBoxLayout(box)
        layout.addWidget(self._bouton(
            "marche", "Calculer la marche diurne et le DVH", self.on_calculer_marche,
            "Marche par position en fonction de l'amplitude, avec les contributions "
            "séparées : gravité, isochronisme propre, échappement, jeu de goupilles.\n\n"
            "Enchaîne ensuite le DVH à HAUTE RÉSOLUTION — marche en position COSC 6H et "
            "à plat sur toute la plage, au maillage de travail et à un maillage 4x plus "
            "fin. Le DVH subit la même erreur de discrétisation absolue que la "
            "dispersion mais vaut au moins deux fois moins qu'elle : lire l'un sans "
            "l'autre donnerait une confiance que le chiffre ne mérite pas.\n\n"
            "Calcul long : trois balayages elastica complets."))
        layout.addWidget(self._affichage("dvh"))
        layout.addWidget(self._bouton(
            "anim", "Calculer l'animation et le jeu minimal", self.on_calculer_animation,
            "Calcule 41 formes déformées entre -amplitude_max et +amplitude_max, puis "
            "cherche l'angle où le jeu inter-spires est le plus faible et affiche la "
            "forme déformée à cet angle.\n\n"
            "Les deux vont ensemble : le jeu minimal est un instantané du déploiement "
            "qu'on vient d'animer, et il se lit sur le même balayage."))
        layout.addWidget(self._affichage("gif"))

        box_num = QGroupBox("Résultats numériques")
        v_num = QVBoxLayout(box_num)
        resume = self._affichage("resultats", hauteur_min=110)
        resume.setText("(aucune géométrie)")
        v_num.addWidget(resume)
        v_num.addWidget(self._bouton(
            "haute_res", "Vérifier à haute résolution", self.on_verifier_haute_resolution,
            "Recalcule fréquence, jeu minimal, dispersion, DVH et isochronisme propre "
            "avec un maillage 4x plus fin, et affiche l'écart absolu ET relatif. Un "
            "maillage grossier surestime le risque : c'est le sens conservateur, mais "
            "à vérifier avant de retenir un candidat."))
        v_num.addWidget(self._affichage("haute_res", hauteur_min=90))
        v_num.addWidget(self._affichage("thermique", hauteur_min=60))
        layout.addWidget(box_num)
        return box

    # ------------------------------------------------------------------
    # Sauvegardes, une par onglet
    # ------------------------------------------------------------------
    def _panneau_sauvegarde_forme(self) -> QGroupBox:
        box = QGroupBox("Enregistrer la forme")
        box.setToolTip(
            "Écrit la CONFIGURATION seule : spec du mouvement, repère, wafer et vecteur "
            "de conception. De quoi reconstruire exactement cette forme, sans les "
            "résultats de calcul."
        )
        v = QVBoxLayout(box)
        v.addWidget(self._bouton("sauv_forme", "Enregistrer la configuration (.json)",
                                 self.on_enregistrer_forme))
        v.addWidget(self._affichage("sauv_forme"))
        return box

    def _panneau_sauvegarde_calculs(self) -> QGroupBox:
        box = QGroupBox("Enregistrer les calculs")
        box.setToolTip(
            "Écrit la configuration, les grandeurs calculées et la géométrie en DXF "
            "dans un dossier horodaté de Results/. Deux DXF si la lame est oxydée : "
            "cote de gravure et pièce finie, qui diffèrent d'environ 4 µm."
        )
        v = QVBoxLayout(box)
        v.addWidget(self._bouton(
            "sauv_calc", "Tout exporter (config, calculs, DXF, rapport, GIF)",
            self.on_enregistrer_calculs,
            "Écrit TOUT dans un dossier horodaté de Results/, au nom du projet :\n"
            "  • la configuration, de quoi reconstruire la forme à l'identique ;\n"
            "  • les grandeurs calculées, en JSON, pour comparer d'une étude à l'autre ;\n"
            "  • la géométrie en DXF — deux fichiers si la lame est oxydée, cote de "
            "gravure et pièce finie, qui diffèrent d'environ 4 µm ;\n"
            "  • le rapport de simulation Word et l'animation GIF.\n\n"
            "Un seul dossier par session, pour qu'elle se transmette ou s'archive d'un "
            "bloc plutôt qu'en fichiers dispersés."))
        v.addWidget(self._affichage("sauv_calc"))
        v.addWidget(self._affichage("projet"))
        return box

    def _panneau_sauvegarde_optimisation(self) -> QGroupBox:
        box = QGroupBox("Enregistrer l'optimisation")
        box.setToolTip(
            "Écrit le spiral de départ, les réglages NSGA-II, et les solutions COCHÉES "
            "dans la table ci-dessus — chacune en DXF (gravure et fini) et en JSON, "
            "plus un index.csv qui les récapitule avec leurs métriques."
        )
        v = QVBoxLayout(box)
        v.addWidget(self._bouton("sauv_optim", "Enregistrer NSGA + solutions cochées",
                                 self.on_enregistrer_optimisation))
        v.addWidget(self._affichage("sauv_optim"))
        return box

    # ------------------------------------------------------------------
    # Enregistrements : un dossier horodaté par session, dans Results/
    # ------------------------------------------------------------------
    def _reinitialiser_analyses(self):
        """Efface TOUT ce qui découle de la géométrie précédente : résultats
        numériques, tracés, animation, solutions d'optimisation, et le dossier
        de session en cours.

        Appelé quand on change d'origine — chargement d'une configuration ou
        d'un DXF. Laisser en place la marche, le front de Pareto ou l'animation
        de la pièce précédente est pire qu'un affichage vide : rien ne signale
        qu'ils ne décrivent plus ce qu'on a sous les yeux, et on les lit comme
        s'ils étaient à jour.
        """
        for cle in self._affichages:
            self._dire(cle, "")
        self._dire("resultats", "(aucune géométrie)")

        for nom in ("_dernier_jeu_min_um", "_dernier_jeu_min_angle", "_dernier_dispersion",
                    "_dernier_dvh", "_dernier_iso_basse", "_dernier_ecart_compensation",
                    "_dernier_diag_iso", "_dernier_jeu_basse", "_dernier_diag_jeu",
                    "_dernier_force_pivot"):
            setattr(self, nom, None)

        # Un nouveau point de départ, c'est une nouvelle étude : elle mérite
        # son propre dossier horodaté plutôt que de venir s'ajouter au
        # précédent.
        self._dossier_courant = None

        self.solutions_pareto = []
        self.sondage_resultats = None
        for table in (self.table_pareto, self.table_sondage):
            table.setRowCount(0)

        self._anim_timer.stop()
        self._anim_angles = None
        self._anim_frames = None
        self._anim_sequence = []
        self._anim_pos = 0
        self.slider_anim.setEnabled(False)
        for b in (self.btn_anim_play, self.btn_anim_prev, self.btn_anim_next):
            b.setEnabled(False)
        self.label_anim_angle.setText("angle = --")
        self.ax_anim.clear()
        self.canvas_anim.draw_idle()

        for axes, canvas in ((("ax_marche",), "canvas_marche"),
                              (("ax_dvh_positions", "ax_dvh_ecart"), "canvas_dvh"),
                              (("ax_pareto",), "canvas_pareto"),
                              (("ax_epaisseur",), "canvas_epaisseur"),
                              (("ax_jeu_min",), "canvas_jeu_min"),
                              (("ax_cg_xy", "ax_cg_angle", "ax_cg_barre"), "canvas_cg"),
                              (("ax_thermique",), "canvas_thermique")):
            for nom in axes:
                getattr(self, nom).clear()
            getattr(self, canvas).draw_idle()

    def _origine_geometrie(self) -> str:
        """D'où vient la géométrie active, pour l'en-tête du rapport. Un
        rapport sans cette ligne oblige à deviner si les chiffres décrivent un
        dessin importé ou une forme calculée."""
        texte = self.label_dxf.text() if hasattr(self, "label_dxf") else ""
        return texte.strip("() ") if texte else ""

    def _n_points_dxf(self):
        """Densité demandée pour les DXF, ou None pour reprendre le maillage de
        calcul. Zéro dans l'interface veut dire « ne rien changer » : plus
        parlant qu'une case à cocher séparée, et sans état incohérent possible
        (case décochée mais valeur saisie)."""
        n = self.sb_n_points_dxf.value() if hasattr(self, "sb_n_points_dxf") else 0
        return n if n >= 8 else None

    def _dossier_session(self, cle_affichage: str):
        """Crée (une fois) le dossier de la session courante. Toutes les
        sauvegardes d'une même session y atterrissent, pour qu'elles se
        transmettent ou s'archivent d'un bloc."""
        if getattr(self, "_dossier_courant", None) and os.path.isdir(self._dossier_courant):
            return self._dossier_courant
        try:
            self._dossier_courant = sauvegarde.dossier_session(
                _RACINE_RESULTS, self.ed_projet.text())
        except ValueError as e:
            QMessageBox.warning(self, "Enregistrement", str(e))
            return None
        return self._dossier_courant

    def _annoncer(self, cle: str, dossier: str, chemins):
        noms = [os.path.basename(c) for c in chemins]
        self._dire_ecourte(cle, f"{len(noms)} fichier(s) dans "
                                f"{os.path.basename(dossier)} : {', '.join(noms)}")
        self._infobulle(cle, dossier + "\n" + "\n".join(noms))

    def on_enregistrer_forme(self):
        """Onglet 1 : la configuration seule."""
        self._synchroniser_spec()
        if self.base is None:
            QMessageBox.information(self, "Enregistrement",
                                     "Générez d'abord une géométrie (panneau 2).")
            return
        dossier = self._dossier_session("sauv_forme")
        if dossier is None:
            return
        chemin = sauvegarde.enregistrer_configuration(
            dossier, self.spec, self.base, self.bornes, self.hyperparams,
            self._reglage_echappement())
        self._annoncer("sauv_forme", dossier, [chemin])

    def on_enregistrer_calculs(self):
        """Onglet 2 : configuration, grandeurs calculées et DXF."""
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Enregistrement",
                                     "Générez ou chargez d'abord une géométrie.")
            return
        dossier = self._dossier_session("sauv_calc")
        if dossier is None:
            return
        self._dire("sauv_calc", "Calcul en cours...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            r_min = float(self.geom.rayon()[0])
            resultat = analyser_marches(
                self.geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
                amplitudes_deg=self._amplitudes_analyse(),
                reglage_echappement=self._reglage_echappement(),
                positions=positions_verticales(self.chk_cote_rouage.isChecked()),
                pas_balayage_deg=self.sb_pas_balayage.value(),
                jeu_goupilles_um=self._jeu_goupilles_um())
            diag = resultat["diag_orbite"]
            # On n'enregistre PAS le résultat brut : il porte des tableaux de
            # travail (orbite complète, balayage) qui pèsent lourd et ne se
            # relisent pas. Seules les grandeurs qu'on compare d'une étude à
            # l'autre sont écrites.
            analyse = {
                "projet": self.ed_projet.text(),
                "frequence_hz": frequence_propre(self.geom, self.spec.inertie_balancier),
                "frequence_cible_hz": self.spec.freq_cible,
                "masse_ug": self.geom.masse_totale(RHO_SILICIUM) * 1e9,
                "longueur_mm": self.geom.longueur,
                "longueur_active_mm": self.geom.partie_active().longueur,
                "epaisseur_min_um": self.geom.epaisseur_finie.min() * 1000.0,
                "epaisseur_max_um": self.geom.epaisseur_finie.max() * 1000.0,
                "epaisseur_oxyde_um": self.geom.epaisseur_oxyde * 1000.0,
                "plan_cristallin": self.spec.plan_cristallin,
                "orientation_cristal_deg": self.spec.orientation_cristal_deg,
                "jeu_min_um": float(diag["jeu_min_mm"]) * 1000.0,
                "jeu_min_angle_deg": float(diag["angle_jeu_min_deg"]),
                "dispersion_gravite_s_j": resultat["dispersion_gravite"],
                "dvh_reference_s_j": resultat["dvh_reference"],
                "dvh_max_s_j": resultat["dvh_max"],
                "sensibilite_module_s_j": resultat["sensibilite_module_reference"],
                "sensibilite_phase_deg": resultat["sensibilite_phase_deg"],
                "iso_amplitude_basse_s_j": resultat["iso_amplitude_basse"],
                "ecart_compensation_s_j": resultat["ecart_compensation"],
                "amplitudes_deg": resultat["amplitudes_deg"],
                "marches_par_position_s_j": {k: v for k, v in resultat["totales"].items()},
            }
            chemins = [sauvegarde.enregistrer_analyse(dossier, analyse)]
        finally:
            QApplication.restoreOverrideCursor()

        # Le rapport écrit lui-même la configuration, les DXF et le GIF : on ne
        # les double pas ici. Il reste séparé parce qu'il est long, et que sa
        # progression doit rester visible.
        chemins += self._ecrire_rapport(dossier)
        self._annoncer("sauv_calc", dossier, [c for c in chemins if c])

    def _ecrire_rapport(self, dossier: str) -> list:
        """Rapport Word, GIF, DXF et configuration, dans le dossier de session.
        Retourne les chemins écrits, ou une liste vide en cas d'échec."""
        frames = None
        if self._anim_frames is not None and self._anim_angles is not None:
            frames = (self._anim_angles, self._anim_frames)

        def progression(texte):
            self._dire("projet", texte)
            QApplication.processEvents()

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            fichiers = rapport.enregistrer_projet(
                os.path.dirname(dossier), os.path.basename(dossier),
                self.geom, self.spec, self.base,
                bornes=self.bornes, hyperparametres=self.hyperparams,
                reglage_echappement=self._reglage_echappement(),
                reglage_thermique=self._reglage_thermique(),
                amplitudes_deg=self._amplitudes_analyse(),
                pas_balayage_deg=self.sb_pas_balayage.value(),
                jeu_goupilles_um=self._jeu_goupilles_um(),
                r_min=float(self.geom.rayon()[0]),
                frames_animation=frames, n_frames_gif=31,
                origine=self._origine_geometrie() if hasattr(self, "_origine_geometrie") else "",
                ecraser=True,          # le dossier de session nous appartient
                n_points_dxf=self._n_points_dxf(), progression=progression)
        except Exception as e:
            QMessageBox.critical(self, "Rapport", str(e))
            self._dire("projet", "(rapport : échec)")
            return []
        finally:
            QApplication.restoreOverrideCursor()
        self._dire("projet", "")
        chemins = []
        for v in fichiers.values():
            chemins += v if isinstance(v, (list, tuple)) else [v]
        return [c for c in chemins if isinstance(c, str)]

    def _solutions_cochees(self) -> list:
        coches = []
        for i in range(self.table_pareto.rowCount()):
            case = self.table_pareto.item(i, 0)
            if case is not None and case.checkState() == Qt.Checked:
                coches.append(i)
        return coches

    def on_enregistrer_optimisation(self):
        """Onglet 3 : réglages NSGA-II, spiral de départ, et les solutions
        cochées — chacune en DXF et en JSON, plus un index.csv."""
        self._synchroniser_spec()
        if not self.solutions_pareto:
            QMessageBox.information(self, "Enregistrement",
                                     "Lancez d'abord une optimisation.")
            return
        rangs = self._solutions_cochees()
        if not rangs:
            QMessageBox.information(
                self, "Enregistrement",
                "Aucune solution cochée. Cochez dans la première colonne de la table "
                "celles dont vous voulez les fichiers de gravure.")
            return
        dossier = self._dossier_session("sauv_optim")
        if dossier is None:
            return
        self._dire("sauv_optim", f"Export de {len(rangs)} solution(s)...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            chemins = [sauvegarde.enregistrer_configuration(
                dossier, self.spec, self.base, self.bornes, self.hyperparams,
                self._reglage_echappement())]
            chemins.append(sauvegarde.enregistrer_analyse(dossier, {
                "projet": self.ed_projet.text(),
                "hyperparametres_ga": asdict(self.hyperparams),
                "bornes": asdict(self.bornes),
                "echappement": asdict(self._reglage_echappement()),
                "n_solutions": len(self.solutions_pareto),
                "solutions_exportees": [r + 1 for r in rangs],
            }, nom="nsga"))

            solutions = []
            for r in rangs:
                sol = self.solutions_pareto[r]
                params = params_depuis_vecteur(sol["x"], self.base)
                geom = construire_geometrie(
                    params, reglage_thermique=self._reglage_thermique())
                solutions.append({
                    "rang": r + 1,
                    "geometrie": geom,
                    "dispersion": sol.get("dispersion"),
                    "dvh": sol.get("dvh"),
                    "jeu_um": (sol.get("jeu") or 0.0) * 1000.0,
                    "ecart_compensation": sol.get("ecart_compensation"),
                    "oxyde_um": geom.epaisseur_oxyde * 1000.0,
                    "frequence_hz": frequence_propre(geom, self.spec.inertie_balancier),
                })
            ecrits = sauvegarde.enregistrer_solutions(
                dossier, solutions, n_points_dxf=self._n_points_dxf())
            for liste in ecrits.values():
                chemins += liste
        finally:
            QApplication.restoreOverrideCursor()
        self._annoncer("sauv_optim", dossier, chemins)

    def _maj_orientation_cristal(self):
        """L'orientation n'a de sens que sur un plan anisotrope : sur un (111)
        le module est le même dans toutes les directions, le champ est donc
        grisé plutôt que laissé actif et sans effet."""
        isotrope = self.cb_plan.currentData() == "111"
        self.sb_orientation_cristal.setEnabled(not isotrope)

    def _lire_spec_depuis_formulaire(self) -> SpecMouvement:
        return SpecMouvement(
            r_min=self.sb_r_min.value(), r_max=self.sb_r_max.value(),
            n_tours_reels=self.sb_n_tours.value(), n_tours_fictifs=self.sb_n_tours_fictifs.value(),
            hauteur_h=self.sb_hauteur.value(), freq_cible=self.sb_freq.value(),
            inertie_balancier=_inertie_interne(self.sb_inertie.value()),
            jeu_min_spires=self.sb_jeu_min.value(),
            amplitude_max=self.sb_amplitude_max.value(),
            vue_cote_rouage=(self.chk_cote_rouage.isChecked()
                             if hasattr(self, "chk_cote_rouage") else True),
            angle_echappement_deg=(self.sb_angle_echappement.value()
                                   if hasattr(self, "sb_angle_echappement") else -6.2),
            sens_trigonometrique=(self.chk_sens_trigo.isChecked()
                                  if hasattr(self, "chk_sens_trigo") else True),
            collet_cote_ancre=(self.chk_collet_cote_ancre.isChecked()
                               if hasattr(self, "chk_collet_cote_ancre") else True),
            plan_cristallin=(self.cb_plan.currentData()
                             if hasattr(self, "cb_plan") else "111"),
            orientation_cristal_deg=(self.sb_orientation_cristal.value()
                                     if hasattr(self, "sb_orientation_cristal") else 0.0),
        )

    def _synchroniser_spec(self):
        """Recharge self.spec depuis les champs actuellement affichés, SANS
        reconstruire la géométrie -- à appeler en tête de toute action qui lit
        self.spec (marche, GIF, jeu minimal, validation, sondage,
        optimisation...), pour ne jamais utiliser une spec périmée. Avant
        cette correction, self.spec n'était resynchronisé que par certaines
        actions (Initialiser Archimède, Générer géométrie, Charger spec...) :
        partir d'un DXF importé sans jamais passer par ces actions-là pouvait
        laisser self.spec sur ses valeurs par défaut, indépendamment de ce
        qui était affiché à l'écran (diagnostic 2026-09-24)."""
        self.spec = self._lire_spec_depuis_formulaire()

    def _rafraichir_base_depuis_spec(self):
        self.spec = self._lire_spec_depuis_formulaire()
        e_theo = epaisseur_theorique(
            self.spec.r_min, self.spec.r_max, self.spec.n_tours_reels,
            self.spec.hauteur_h, self.spec.freq_cible, self.spec.inertie_balancier,
        )
        n_e = max(2, self.sb_n_noeuds_ep.value()) if hasattr(self, "sb_n_noeuds_ep") else 20
        n_dr = self.sb_n_noeuds_dr.value() if hasattr(self, "sb_n_noeuds_dr") else 12
        n_g_pts = self.sb_gross_n_points.value() if hasattr(self, "sb_gross_n_points") else 4
        n_p_pts = self.sb_phil_n_points.value() if hasattr(self, "sb_phil_n_points") else 4

        self.base = ParametresSpiral(
            r_min=self.spec.r_min, r_max=self.spec.r_max,
            n_tours_reels=self.spec.n_tours_reels, n_tours_fictifs=self.spec.n_tours_fictifs,
            hauteur_h=self.spec.hauteur_h, freq_cible=self.spec.freq_cible,
            inertie_balancier=self.spec.inertie_balancier,
            # Repère du mouvement : c'est lui qui ORIENTE le spiral, donc qui
            # fixe la phase du vecteur de sensibilité à la gravité et le DVH.
            angle_echappement_deg=self.spec.angle_echappement_deg,
            sens_trigonometrique=self.spec.sens_trigonometrique,
            collet_cote_ancre=self.spec.collet_cote_ancre,
            plan_cristallin=self.spec.plan_cristallin,
            orientation_cristal_deg=self.spec.orientation_cristal_deg,
            n_elements_corps=self.sb_n_corps.value() if hasattr(self, "sb_n_corps") else 1400,
            grossmann_angle_deg=self.sb_gross_angle.value() if hasattr(self, "sb_gross_angle") else 90.0,
            grossmann_n_pas=self.sb_gross_n_pas.value() if hasattr(self, "sb_gross_n_pas") else 120,
            grossmann_n_points=n_g_pts, grossmann_offsets=np.zeros(n_g_pts - 3),
            phillips_angle_deg=self.sb_phil_angle.value() if hasattr(self, "sb_phil_angle") else 120.0,
            phillips_n_pas=self.sb_phil_n_pas.value() if hasattr(self, "sb_phil_n_pas") else 150,
            phillips_n_points=n_p_pts, phillips_offsets=np.zeros(n_p_pts - 3),
            epaisseur_noeuds=np.full(n_e, e_theo),
            delta_r_noeuds=np.zeros(max(n_dr, 2)) if n_dr > 0 else np.zeros(2),
            **self._parametres_raquette(),
        )

    def _parametres_raquette(self) -> dict:
        """Champs de terminaison à raquette lus depuis l'interface. Les widgets
        peuvent ne pas exister encore au tout premier rafraîchissement (le
        panneau est construit après) -- on retombe alors sur les valeurs par
        défaut de la dataclass, c'est-à-dire spiral libre."""
        if not hasattr(self, "chk_raquette"):
            return {}
        n_pts = self.sb_raq_n_points.value()
        return {
            "raquette_active": self.chk_raquette.isChecked(),
            "raquette_plage_s_par_j": self.sb_raq_plage.value(),
            "raquette_transition_angle_deg": self.sb_raq_trans_angle.value(),
            "raquette_n_points": n_pts,
            "raquette_offsets": np.zeros(n_pts - 4),
            "raquette_transition_n_pas": self.sb_raq_trans_n_pas.value(),
            "raquette_arc_avant_deg": self.sb_raq_arc_avant.value(),
            "raquette_arc_avant_n_pas": self.sb_raq_n_pas_arc.value(),
            "raquette_piton_angle_deg": self.sb_raq_piton_angle.value(),
            "raquette_piton_n_pas": self.sb_raq_piton_n_pas.value(),
            "raquette_jeu_um": self.sb_raq_jeu.value(),
            "oxydation_active": self.chk_oxydation.isChecked() if hasattr(self, "chk_oxydation") else False,
            "oxyde_impose_um": self.sb_ox_impose.value() if hasattr(self, "sb_ox_impose") else 0.0,
        }

    def on_init_archimede(self):
        self._rafraichir_base_depuis_spec()
        try:
            geom, diag = construire_geometrie(self.base, return_diagnostics=True)
        except ValueError as e:
            QMessageBox.warning(self, "Erreur", str(e))
            return
        self.geom = geom
        self.label_dxf.setText("(spiral d'Archimède généré)")
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()

    def on_charger_dxf(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Charger un DXF", "", "DXF (*.dxf)")
        if not chemin:
            return
        try:
            # n_elements=3000 : l'étude de convergence du 2026-09-25 montre un
            # biais encore ~10% sur le jeu minimal à 1500 éléments (toujours
            # dans le sens conservateur -- surestime le risque, pas l'inverse),
            # ramené à ~5% à 3000 pour un coût de calcul encore négligeable
            # pour un import ponctuel (voir aussi "Vérifier à haute
            # résolution" au panneau 5 pour une geometrie generee).
            geom = dxf_io.charger_dxf(
                chemin, n_elements=3000, hauteur_h=self.sb_hauteur.value(),
                plan_cristallin=self.cb_plan.currentData(),
                orientation_cristal_deg=self.sb_orientation_cristal.value())
        except Exception as e:
            QMessageBox.critical(self, "Erreur d'import", str(e))
            return
        self._reinitialiser_analyses()
        self.geom = geom
        _texte_ecourte(self.label_dxf, chemin.split("/")[-1])
        r = geom.rayon()
        self.sb_r_min.setValue(float(r.min()))
        self.sb_r_max.setValue(float(r.max()))
        self._tracer_geometrie()
        QMessageBox.information(
            self, "DXF chargé",
            "Géométrie chargée. Pour l'analyser (marche diurne, jeu minimal, GIF, "
            "export), utilisez directement le panneau \"5. Résultats et export\" "
            "ci-dessous -- les panneaux 2 à 4 ne servent qu'à CONSTRUIRE ou OPTIMISER "
            "une géométrie, pas à en analyser une déjà importée.",
        )

    def on_generer_geometrie(self):
        self._rafraichir_base_depuis_spec()
        try:
            geom, diag = construire_geometrie(self.base, return_diagnostics=True)
        except ValueError as e:
            QMessageBox.warning(self, "Erreur", str(e))
            return
        self.geom = geom
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()

    def _afficher_diagnostics(self, diag: dict):
        raquette = diag.get("raquette_active", False)
        nom_ext = "courbe en S" if raquette else "Phillips"
        # Le budget est détaillé et non résumé : le champ « Corps » ne pilote que
        # le corps, et la différence avec le total a déjà prêté à confusion.
        detail = f"{self.base.n_elements_corps} corps + {self.base.grossmann_n_pas} Grossmann"
        if raquette:
            detail += (f" + {self.base.raquette_transition_n_pas} S"
                       f" + {self.base.n_pas_arc_avant} arc"
                       f" + {self.base.n_pas_piton} piton")
        else:
            detail += f" + {self.base.phillips_n_pas} Phillips"
        lignes = [
            f"maillage total        : {self.base.n_elements}  ({detail})",
            f"geometrie_valide      : {diag['geometrie_valide']}",
            f"grossmann_monotone    : {diag['grossmann_monotone']}",
            f"{nom_ext + ' monotone':22s}: {diag['phillips_monotone']}",
            f"phillips_plausibilite : {diag['phillips_ratio_plausibilite']:.2f}  (repère ~1.0, pas une contrainte)",
            f"freq avant scaling    : {diag['freq_avant_scaling_hz']:.3f} Hz  (beta={diag['beta_scaling']:.3f})",
        ]
        if raquette:
            lignes.append(
                f"arc de raquette       : parcours {diag['raquette_arc_deg']:.2f}° "
                f"({diag['raquette_arc_mm']*1000:.0f} µm) dans un arc de "
                f"{diag['raquette_arc_avant_deg']:.0f}°+{diag['raquette_piton_angle_deg']:.0f}°, "
                f"longueur morte {diag['longueur_morte_mm']*1000:.0f} µm"
            )
            if not diag.get("raquette_parcours_tient_dans_arc", True):
                lignes.append(
                    f"ATTENTION : la demi-course des goupilles ({diag['raquette_demi_parcours_deg']:.1f}°) "
                    f"dépasse l'arc concentrique disponible -- elles en sortiraient et "
                    f"déformeraient le spiral radialement. Allonger l'arc ou réduire la plage."
                )
        self.label_diag.setText("\n".join(lignes))

    def _tracer_geometrie(self):
        # Une nouvelle géométrie invalide le jeu minimal et la dispersion de
        # marche affichés (calculés pour l'ANCIENNE géométrie) -- on les
        # efface plutôt que de laisser afficher une valeur périmée.
        self._dernier_jeu_min_um = None
        self._dernier_dispersion = None
        self._dernier_dvh = None
        self._dernier_iso_basse = None
        self._dernier_ecart_compensation = None
        self._dernier_diag_iso = None
        self._dernier_jeu_basse = None
        self._dernier_diag_jeu = None
        self._dernier_force_pivot = None

        self.ax_geom.clear()
        self.ax_zoom_grossmann.clear()
        self.ax_zoom_phillips.clear()
        if self.geom is not None:
            ext, inte = self.geom.contour_flancs()
            self.ax_geom.plot(ext[:, 0], ext[:, 1], "-", color=theme.encre(), lw=0.8)
            self.ax_geom.plot(inte[:, 0], inte[:, 1], "-", color=theme.encre(), lw=0.8)
            self.ax_geom.plot(0, 0, "+", color=theme.encre(), ms=10)
            tracer_repere_mouvement(
                self.ax_geom, float(self.geom.rayon().max()),
                self.chk_cote_rouage.isChecked(),
                angle_collet_deg=self.base.angle_collet_deg if self.base else None,
                angle_ancre_deg=self.base.angle_ancre_deg if self.base else None,
                rayon_min=float(self.geom.rayon().min()))
            self.ax_geom.set_aspect("equal")
            self.ax_geom.set_xlabel("x (mm)")
            self.ax_geom.set_ylabel("y (mm)")

            # Zooms dédiés collet/Grossmann et piton/Phillips : ces courbes ne
            # couvrent que quelques dizaines à quelques centaines de microns,
            # invisibles sur la vue d'ensemble (cf. échange du 2026-09-24).
            # Utilise les n_pas exacts si la géométrie active vient bien de
            # `self.base` (générée paramétriquement) ; sinon (DXF importé),
            # se rabat sur une fraction de la longueur totale.
            if self.base is not None and self.geom.n == self.base.n_elements:
                n_g = self.base.grossmann_n_pas
                n_p = (self.base.raquette_transition_n_pas + self.base.n_pas_arc_avant
                       + self.base.n_pas_piton
                       if self.base.raquette_active else self.base.phillips_n_pas)
            else:
                n_g = n_p = max(5, self.geom.n // 20)
            marge_g = n_g // 2
            marge_p = n_p // 2

            i_fin_g = min(n_g + marge_g, self.geom.n)
            self.ax_zoom_grossmann.plot(ext[:i_fin_g, 0], ext[:i_fin_g, 1], ".-", color=theme.encre(), ms=2, lw=0.8)
            self.ax_zoom_grossmann.plot(inte[:i_fin_g, 0], inte[:i_fin_g, 1], ".-", color=theme.encre(), ms=2, lw=0.8)
            self.ax_zoom_grossmann.set_aspect("equal")
            self.ax_zoom_grossmann.set_title("Zoom collet (Grossmann)", fontsize=9)
            self.ax_zoom_grossmann.tick_params(labelsize=7)

            i_debut_p = max(self.geom.n - n_p - marge_p, 0)
            self.ax_zoom_phillips.plot(ext[i_debut_p:, 0], ext[i_debut_p:, 1], ".-", color=theme.encre(), ms=2, lw=0.8)
            self.ax_zoom_phillips.plot(inte[i_debut_p:, 0], inte[i_debut_p:, 1], ".-", color=theme.encre(), ms=2, lw=0.8)
            self.ax_zoom_phillips.set_aspect("equal")
            self.ax_zoom_phillips.tick_params(labelsize=7)

            if self.geom.a_partie_morte:
                # La portion morte (au-delà des goupilles) est tracée en gris :
                # elle est bien usinée et occupe de la place, mais ne travaille
                # pas. Sans ce repère, rien ne distingue à l'écran la longueur
                # active de la longueur totale.
                na = self.geom.n_actif
                for ax, lw in ((self.ax_geom, 1.6), (self.ax_zoom_phillips, 2.2)):
                    ax.plot(ext[na - 1:, 0], ext[na - 1:, 1], "-", color=theme.gris(0.65), lw=lw)
                    ax.plot(inte[na - 1:, 0], inte[na - 1:, 1], "-", color=theme.gris(0.65), lw=lw)
                    ax.plot(self.geom.x[na - 1], self.geom.y[na - 1], "o",
                            color="tab:red", ms=6, mfc="none", mew=1.5)
                self.ax_zoom_phillips.set_title(
                    "Zoom piton -- arc de raquette\n"
                    "rouge = point de comptage, gris = portion morte", fontsize=8)
            else:
                self.ax_zoom_phillips.set_title("Zoom piton (Phillips)", fontsize=9)

        self.fig_geom.tight_layout()
        self.canvas_geom.draw_idle()
        self._tracer_epaisseur()
        self._afficher_resultats_numeriques()

    def _afficher_resultats_numeriques(self):
        """Résumé numérique de la géométrie active -- en particulier la
        fréquence propre atteinte (absente jusqu'ici : rien n'indiquait où
        l'oscillateur atterrissait réellement par rapport à freq_cible après
        optimisation). Se met à jour automatiquement à chaque changement de
        géométrie ; le jeu minimal et la dispersion de marche, plus coûteux à
        calculer, restent affichés depuis leur dernier calcul explicite
        (boutons dédiés) plutôt que recalculés à chaque tracé."""
        self._synchroniser_spec()
        if self.geom is None:
            self._dire("resultats", "(aucune géométrie)")
            return

        freq = frequence_propre(self.geom, self.spec.inertie_balancier)
        ecart_pct = 100.0 * (freq - self.spec.freq_cible) / self.spec.freq_cible
        masse_ug = self.geom.masse_totale(RHO_SILICIUM) * 1e9

        lignes = [
            f"Fréquence propre : {freq:.4f} Hz",
            f"Cible : {self.spec.freq_cible:.4f} Hz (écart {ecart_pct:+.2f}%)",
            f"Masse : {masse_ug:.2f} µg",
        ]
        if self.geom.a_partie_morte:
            # Plage VÉRIFIÉE sur la géométrie obtenue, pas la plage demandée :
            # l'angle de l'arc est dimensionné analytiquement, il faut contrôler
            # ce qu'il donne réellement (même logique que la mise à l'échelle de
            # l'épaisseur, calculée puis vérifiée).
            rapide, lente, etendue = marche_raquette(self.geom, self.spec.inertie_balancier)
            demande = self.sb_raq_plage.value() if hasattr(self, "sb_raq_plage") else 0.0
            lignes.append(
                f"Raquette : {rapide:+.1f} / {lente:+.1f} s/j "
                f"(étendue {etendue:.1f} s/j pour ±{demande:.0f} demandés)")
            lignes.append(
                f"Longueur active : {self.geom.partie_active().longueur:.2f} mm "
                f"sur {self.geom.longueur:.2f} mm")
            if self._dernier_jeu_basse is not None:
                d = self._dernier_diag_jeu or {}
                contact = d.get("amplitude_contact_deg", float("inf"))
                lignes.append(
                    f"Jeu de goupilles : {self._dernier_jeu_basse:+.2f} s/j à "
                    f"{self._reglage_echappement().amplitude_basse_deg:.0f}° "
                    f"(contact dès {contact:.0f}° de balancier)")
                fr = self._reglage_echappement().fraction_jeu()
                reste_raquette = abs(self._dernier_jeu_basse) * (1.0 - fr)
                reste_piton = abs(self._dernier_jeu_basse) * fr
                lignes.append(
                    f"Montage visé : {self.cb_variante.currentText().lower()} -- "
                    f"erreur résiduelle à l'amplitude basse : "
                    f"{reste_raquette:.2f} s/j à raquette, {reste_piton:.2f} s/j en piton fixe")
                if not np.isfinite(contact) or contact > self.spec.amplitude_max:
                    lignes.append(
                        "ATTENTION : la lame n'atteint jamais les goupilles dans la plage "
                        "d'amplitude -- la raquette n'a aucune autorité. Allonger l'arc "
                        "au-delà des goupilles, ou resserrer le jeu.")
        if self._dernier_jeu_min_um is not None:
            lignes.append(f"Jeu minimal : {self._dernier_jeu_min_um:.1f} µm à {self._dernier_jeu_min_angle:+.0f}°")
        else:
            lignes.append("Jeu minimal : non calculé")
        if self._dernier_dispersion is not None:
            lignes.append(f"Dispersion de marche (gravité, verticales) : {self._dernier_dispersion:.2f} s/j")
        else:
            lignes.append("Dispersion de marche : non calculée")
        if self._dernier_dvh is not None:
            dvh, dvh_max, module, phase, amp_dvh = self._dernier_dvh
            lignes.append(
                f"DVH (COSC, 6H moins horizontale) : {dvh:+.2f} s/j à {amp_dvh:.0f}° "
                f"(pire cas sur la plage : {dvh_max:.2f} s/j)")
            # La marche en verticale suit exactement A.cos(theta - phi). Le
            # module ne dépend que de la FORME, la phase de l'ORIENTATION dans
            # le mouvement : afficher les deux dit d'un coup d'oeil ce qu'une
            # rotation du spiral peut encore gagner, et ce qu'elle ne peut pas.
            lignes.append(
                f"Sensibilité à la gravité : module {module:.2f} s/j, phase {phase:.1f}° "
                f"(marche = A·cos(θ − φ) ; une rotation change φ, pas A)")

        if self._dernier_force_pivot is not None:
            # Treize fois le poids du spiral : c'est cette force, et non la
            # poise, qui justifie la règle du nombre entier de spires. L'axe
            # est une DROITE, rendu modulo 180° comme le point de comptage.
            module, axe = self._dernier_force_pivot
            comptage = angle_point_comptage(self.geom)
            ecart = (axe - comptage + 90.0) % 180.0 - 90.0
            poids = self.geom.masse_totale(RHO_SILICIUM) * 9.81 * 1e6
            lignes.append(
                f"Force sur le pivot : {module:.1f} µN max ({module/poids:.1f}x le poids), "
                f"axe {axe:.1f}° (point de comptage {comptage:.1f}°, écart {ecart:+.1f}°)")

        reglage = self._reglage_echappement()
        if self._dernier_iso_basse is not None:
            lignes.append(
                f"Isochronisme propre à {reglage.amplitude_basse_deg:.0f}° : "
                f"{self._dernier_iso_basse:+.2f} s/j "
                f"(cible {reglage.retard_s_par_j:+.2f} pour annuler l'échappement)")
            etat = ("dans la tolérance"
                    if abs(self._dernier_ecart_compensation) <= reglage.tolerance_s_par_j
                    else "HORS tolérance")
            lignes.append(
                f"Écart de compensation : {self._dernier_ecart_compensation:+.2f} s/j "
                f"(± {reglage.tolerance_s_par_j:.2f} -- {etat})")
            if self._dernier_diag_iso is not None and not self._dernier_diag_iso.get("valide", False):
                lignes.append(f"ATTENTION isochronisme non calculable : "
                              f"{self._dernier_diag_iso.get('raison', 'inconnue')}")
        else:
            lignes.append("Compensation d'échappement : non calculée")
        self._dire("resultats", "\n".join(lignes))

    def _tracer_epaisseur(self):
        self.ax_epaisseur.clear()
        if self.geom is not None:
            self.ax_epaisseur.plot(self.geom.s, self.geom.epaisseur * 1000.0, "-")
            self.ax_epaisseur.set_xlabel("Longueur depuis la virole/collet (mm)")
            self.ax_epaisseur.set_ylabel("Épaisseur (µm)")
            self.ax_epaisseur.grid(True, linestyle=":", alpha=0.6)
        self.canvas_epaisseur.draw_idle()

    # ------------------------------------------------------------------
    # Logique -- Sondage
    # ------------------------------------------------------------------
    def _amplitudes_analyse(self) -> np.ndarray:
        """Grille d'amplitudes pour l'affichage ET pour l'objectif d'optimisation
        -- pas fixe (PAS_DISPERSION_DEG_DEFAUT, cf. optimisation.py), pour que
        le graphe affiché corresponde toujours exactement à ce que l'optimiseur
        évalue en interne (plus de "amp_pas" séparé, source de confusion)."""
        return np.arange(self.sb_amp_min.value(), self.sb_amp_max.value() + 1e-6, PAS_DISPERSION_DEG_DEFAUT)

    def _definir_bouton_actif(self, bouton: QPushButton, actif: bool, texte_actif: str = None):
        """Retour visuel qu'un traitement en arrière-plan est en cours --
        avant cette correction, un bouton restait identique visuellement
        pendant tout le calcul, seul le texte d'un label séparé changeait
        (facile à manquer). Désactiver le bouton empêche aussi un
        double-lancement qui écraserait self._thread pendant qu'il tourne
        encore (référence perdue sur le QThread en cours, comportement
        indéfini côté PySide6)."""
        if actif:
            bouton.setProperty("_texte_repos", bouton.text())
            if texte_actif:
                bouton.setText(texte_actif)
            bouton.setEnabled(False)
            bouton.setStyleSheet(theme.ACTIF)
        else:
            texte_repos = bouton.property("_texte_repos")
            if texte_repos:
                bouton.setText(texte_repos)
            bouton.setEnabled(True)
            bouton.setStyleSheet("")

    def on_lancer_sondage(self):
        self._synchroniser_spec()
        if not self._verifier_coherence_echappement():
            return
        self._definir_bouton_actif(self.btn_sondage, True, "Sondage en cours...")
        self.label_sondage.setText("Sondage en cours...")
        self._thread = SondageThread(
            self.base, self.bornes, self._amplitudes_analyse(),
            n_echantillons=self.sb_sondage_n.value(), seed=self.sb_sondage_seed.value(),
            pas_balayage_deg=self.sb_pas_balayage.value(),
            reglage_echappement=self._reglage_echappement(),
        )
        self._thread.signal_termine.connect(self.on_sondage_termine)
        self._thread.start()

    def on_sondage_termine(self, resultats: dict):
        self._definir_bouton_actif(self.btn_sondage, False)
        self.sondage_resultats = resultats
        n_ok = int(np.sum(resultats["jeu_mm"] >= self.spec.jeu_min_spires))
        self.label_sondage.setText(
            f"{len(resultats['jeu_mm'])} tirages, {n_ok} au jeu requis "
            f"({self.spec.jeu_min_spires*1000:.1f} µm). Les "
            f"{len(resultats['classement'])} meilleurs ci-dessous."
        )
        self.label_sondage.setToolTip(
            "Le jeu requis n'est PAS un filtre ici : le sondage classe les tirages du "
            "moins mauvais au pire pour fournir un point de départ à l'Optimisation, "
            "qui applique elle la vraie contrainte. Sélectionner une ligne en fait la "
            "géométrie active ; sans sélection, l'Optimisation part du meilleur."
        )
        self.table_sondage.setRowCount(len(resultats["classement"]))
        for i, sol in enumerate(resultats["classement"]):
            self.table_sondage.setItem(i, 0, QTableWidgetItem(f"{sol['dispersion']:.2f}"))
            self.table_sondage.setItem(i, 1, QTableWidgetItem(f"{sol['jeu_mm']*1000:.1f}"))
            ecart = sol.get("ecart_compensation")
            self.table_sondage.setItem(
                i, 2, QTableWidgetItem("--" if ecart is None else f"{ecart:+.2f}"))

    def on_selectionner_sondage(self):
        lignes = self.table_sondage.selectionModel().selectedRows()
        if not lignes or not self.sondage_resultats:
            return
        idx = lignes[0].row()
        x = self.sondage_resultats["classement"][idx]["x"]
        self.base = params_depuis_vecteur(x, self.base)
        geom, diag = construire_geometrie(self.base, return_diagnostics=True)
        self.geom = geom
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()

    # ------------------------------------------------------------------
    # Logique -- Optimisation
    # ------------------------------------------------------------------
    def on_lancer_optimisation(self):
        self._synchroniser_spec()
        if not self._verifier_coherence_echappement():
            return
        self._definir_bouton_actif(self.btn_optimisation, True, "Optimisation en cours...")
        self.hyperparams = HyperparametresGA(
            pop_size=self.sb_pop.value(), n_gen=self.sb_gen.value(), seed=self.sb_ga_seed.value(),
        )
        x0 = self.sondage_resultats["x_meilleur"] if self.sondage_resultats is not None else None
        self.label_optim.setText("Optimisation en cours... (génération 0)")
        self._thread = OptimisationThread(
            self.base, self.bornes, self.hyperparams, self.spec.jeu_min_spires,
            self._amplitudes_analyse(), x0=x0,
        )
        self._thread.signal_progres.connect(self.on_optimisation_progres)
        self._thread.signal_termine.connect(self.on_optimisation_terminee)
        self._thread.start()

    def on_optimisation_progres(self, gen: int, dispersion: float, jeu: float):
        self.label_optim.setText(
            f"Gén. {gen}/{self.hyperparams.n_gen} -- dispersion {dispersion:.2f} s/j, "
            f"jeu {jeu*1000:.1f} µm"
        )

    def on_optimisation_terminee(self, res):
        self._definir_bouton_actif(self.btn_optimisation, False)
        if res.X is None:
            self.label_optim.setText("Terminée -- AUCUNE solution faisable dans le front final.")
            self.solutions_pareto = []
        else:
            X = np.atleast_2d(res.X)
            F = np.atleast_2d(res.F)
            ordre = np.argsort(F[:, 0])
            # L'écart de compensation n'est pas dans F (c'est une contrainte,
            # pas un objectif) : on le relit dans G, si la contrainte est
            # active. G[:,1] = |écart| - tolérance, donc |écart| = G + tolérance.
            G = np.atleast_2d(res.G) if res.G is not None else None
            tol = self._reglage_echappement().tolerance_s_par_j
            def _ecart(i):
                if G is None or G.shape[1] < 2:
                    return None
                return float(G[i, 1]) + tol
            self.solutions_pareto = [
                {"x": X[i], "dispersion": float(F[i, 0]), "jeu": float(-F[i, 1]),
                 "ecart_compensation": _ecart(i)}
                for i in ordre
            ]
            self.label_optim.setText(f"Terminée -- {len(self.solutions_pareto)} solution(s).")

        self.table_pareto.setRowCount(len(self.solutions_pareto))
        for i, sol in enumerate(self.solutions_pareto):
            case = QTableWidgetItem()
            case.setFlags((case.flags() | Qt.ItemIsUserCheckable) & ~Qt.ItemIsEditable)
            case.setCheckState(Qt.Unchecked)
            case.setToolTip("Cocher pour exporter cette géométrie en DXF.")
            self.table_pareto.setItem(i, 0, case)
            self.table_pareto.setItem(i, 1, QTableWidgetItem(f"{sol['dispersion']:.2f}"))
            self.table_pareto.setItem(i, 2, QTableWidgetItem(f"{sol['jeu']*1000:.1f}"))
            ecart = sol.get("ecart_compensation")
            self.table_pareto.setItem(
                i, 3, QTableWidgetItem("--" if ecart is None else f"|{ecart:.2f}|"))
        self.table_pareto.resizeColumnToContents(0)
        self._tracer_pareto()

    def _tracer_pareto(self):
        self.ax_pareto.clear()
        if self.solutions_pareto:
            disp = [s["dispersion"] for s in self.solutions_pareto]
            jeu_um = [s["jeu"] * 1000 for s in self.solutions_pareto]
            self.ax_pareto.scatter(disp, jeu_um)
            self.ax_pareto.set_xlabel("Dispersion de marche (s/j)")
            self.ax_pareto.set_ylabel("Jeu min inter-spires (µm)")
            self.ax_pareto.axhline(self.spec.jeu_min_spires * 1000, color="tab:red", ls="--", lw=0.8, label="jeu min requis")
            self.ax_pareto.legend()
        self.canvas_pareto.draw_idle()

    def on_selectionner_solution(self):
        lignes = self.table_pareto.selectionModel().selectedRows()
        if not lignes or not self.solutions_pareto:
            return
        idx = lignes[0].row()
        sol = self.solutions_pareto[idx]
        self.base = params_depuis_vecteur(sol["x"], self.base)
        self.geom = construire_geometrie(self.base)
        self._tracer_geometrie()

    # ------------------------------------------------------------------
    # Logique -- Résultats / export
    # ------------------------------------------------------------------
    def on_calculer_marche(self):
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Marche diurne", "Générez ou chargez d'abord une géométrie.")
            return
        r_min = float(self.geom.rayon()[0])
        amplitudes = self._amplitudes_analyse()
        reglage = self._reglage_echappement()
        resultat = analyser_marches(
            self.geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
            amplitudes_deg=amplitudes, reglage_echappement=reglage,
            positions=positions_verticales(self.chk_cote_rouage.isChecked()),
            pas_balayage_deg=self.sb_pas_balayage.value(),
            jeu_goupilles_um=self._jeu_goupilles_um(),
        )

        # Deux vues superposées : la marche TOTALE par position (ce que
        # mesurerait un chronocomparateur, seule grandeur comparable à un
        # relevé réel) et, en pointillés, les deux contributions qui ne
        # dépendent PAS de la position -- l'isochronisme propre du spiral et
        # l'erreur d'échappement. Elles se lisent l'une contre l'autre : c'est
        # leur somme, et non la dispersion entre positions, qui dit si le
        # retard d'échappement est compensé.
        self.ax_marche.clear()
        for nom, valeurs in resultat["totales"].items():
            style = "s--" if nom in POSITIONS_HORIZONTALES else "o-"
            largeur = 2.0 if nom in POSITIONS_HORIZONTALES else 1.0
            self.ax_marche.plot(amplitudes, valeurs, style, lw=largeur, ms=4, label=nom)
        self.ax_marche.plot(amplitudes, resultat["isochronisme"], ":", color=theme.gris(0.35),
                            lw=1.5, label="isochronisme propre")
        self.ax_marche.plot(amplitudes, resultat["echappement"], ":", color="tab:red",
                            lw=1.5, label="échappement")
        if resultat["diag_jeu_goupilles"].get("actif"):
            self.ax_marche.plot(amplitudes, resultat["jeu_goupilles"], ":", color="tab:purple",
                                lw=1.5, label="jeu de goupilles")
        self.ax_marche.axvline(reglage.amplitude_basse_deg, color=theme.gris(0.7), ls="-.", lw=0.8)
        self.ax_marche.set_xlabel("Amplitude (deg)")
        self.ax_marche.set_ylabel("Marche diurne (s/j)")
        self.ax_marche.set_title(
            f"Compensation à {reglage.amplitude_basse_deg:.0f}° : "
            f"spiral {resultat['iso_amplitude_basse']:+.2f} s/j "
            f"contre {-reglage.retard_s_par_j:+.2f} s/j d'échappement "
            f"(écart {resultat['ecart_compensation']:+.2f} s/j)", fontsize=9)
        self.ax_marche.axhline(0, color=theme.encre(), lw=0.8)
        self.ax_marche.grid(True, linestyle=":", alpha=0.5)
        self.ax_marche.legend(fontsize=8, ncol=2)
        self.fig_marche.tight_layout()
        self.canvas_marche.draw_idle()
        self.tabs.setCurrentWidget(self.canvas_marche.parentWidget())

        self._dernier_dispersion = resultat["dispersion_gravite"]
        self._dernier_dvh = (resultat["dvh_reference"], resultat["dvh_max"],
                             resultat["sensibilite_module_reference"],
                             resultat["sensibilite_phase_deg"],
                             resultat["amplitude_dvh_deg"])
        self._dernier_iso_basse = resultat["iso_amplitude_basse"]
        self._dernier_ecart_compensation = resultat["ecart_compensation"]
        self._dernier_diag_iso = resultat["diag_isochronisme"]
        self._dernier_jeu_basse = resultat["jeu_amplitude_basse"]
        self._dernier_diag_jeu = resultat["diag_jeu_goupilles"]
        self._tracer_orbite_cg(resultat["diag_orbite"])
        self._dernier_force_pivot = axe_force_pivot(
            self.geom, r_min, self.spec.amplitude_max)
        self._afficher_resultats_numeriques()

        # Le DVH suit automatiquement : son incertitude de maillage est
        # beaucoup plus lourde en relatif que celle de la dispersion, et lire
        # l'une sans l'autre donne une confiance que le chiffre ne mérite pas.
        # Silencieux : sur une géométrie importée, le raffinement n'est pas
        # possible et cela ne doit pas interrompre le calcul de marche par une
        # boîte de dialogue.
        self.on_marche_haute_resolution(silencieux=True)
        self.tabs.setCurrentWidget(self.canvas_marche.parentWidget())

    def on_analyser_thermique(self):
        """Tenue en température de la géométrie active, plus l'étude de
        sensibilité aux coefficients matériau.

        La sensibilité n'est pas un supplément d'agrément : les coefficients
        sont des valeurs de littérature, pas des mesures, et sans elle le
        chiffre de résidu affiché inspirerait une confiance qu'il ne mérite
        pas."""
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Thermique", "Générez ou chargez d'abord une géométrie.")
            return
        reglage = self._reglage_thermique()
        self._definir_bouton_actif(self.btn_thermique, True, "Analyse en cours...")
        QApplication.processEvents()
        try:
            t_impose = self.geom.epaisseur_oxyde * 1000.0
            a_nu = analyser_thermique(self.geom, self.spec.inertie_balancier, reglage,
                                       t_oxyde_um=0.0)
            a_ox = analyser_thermique(self.geom, self.spec.inertie_balancier, reglage,
                                       t_oxyde_um=t_impose if t_impose > 0 else None)
            sensibilite = etude_sensibilite(self.geom, self.spec.inertie_balancier, reglage)
        finally:
            self._definir_bouton_actif(self.btn_thermique, False)

        self.ax_thermique.clear()
        self.ax_thermique.plot(a_nu["temperatures"], a_nu["marches"], "o-", color="tab:red",
                                label=f"silicium nu ({a_nu['pente_s_par_j_par_degc']:+.2f} s/j/°C)")
        if a_ox["valide"]:
            self.ax_thermique.plot(a_ox["temperatures"], a_ox["marches"], "s-", color="tab:blue",
                                    label=f"oxydé {a_ox['epaisseur_oxyde_um']:.2f} µm "
                                          f"({a_ox['pente_s_par_j_par_degc']:+.3f} s/j/°C)")
        # Repère COSC : la bande dans laquelle la marche doit rester pour tenir
        # 0.6 s/j/°C sur la plage visée.
        etendue = reglage.temperature_max - reglage.temperature_min
        self.ax_thermique.axhspan(-0.6 * etendue / 2, 0.6 * etendue / 2, color=theme.gris(0.9),
                                   zorder=0, label="bande COSC (0.6 s/j/°C)")
        self.ax_thermique.axhline(0, color=theme.encre(), lw=0.8)
        self.ax_thermique.set_xlabel("Température (°C)")
        self.ax_thermique.set_ylabel("Marche (s/j)")
        self.ax_thermique.grid(True, linestyle=":", alpha=0.5)
        self.ax_thermique.legend(fontsize=8)
        self.ax_thermique.set_title("Marche en fonction de la température", fontsize=9)
        self.fig_thermique.tight_layout()
        self.canvas_thermique.draw_idle()
        self.tabs.setCurrentWidget(self.canvas_thermique.parentWidget())

        lignes = [
            f"Silicium nu : {a_nu['pente_s_par_j_par_degc']:+.3f} s/j/°C, "
            f"{a_nu['amplitude_s_par_j']:.1f} s/j sur {reglage.temperature_min:.0f}-"
            f"{reglage.temperature_max:.0f} °C",
        ]
        if a_ox["valide"]:
            lignes += [
                f"Oxyde requis : {a_ox['epaisseur_oxyde_um']:.3f} µm "
                f"(lame +{a_ox['croissance_lame_um']:.2f} µm, silicium consommé "
                f"{a_ox['silicium_consomme_um']:.2f} µm)",
                f"Après compensation : {a_ox['pente_s_par_j_par_degc']:+.4f} s/j/°C, "
                f"{a_ox['amplitude_s_par_j']:.3f} s/j crête-à-crête",
                f"Part de rigidité portée par l'oxyde : "
                f"{a_ox['part_oxyde_min']*100:.1f} à {a_ox['part_oxyde_max']*100:.1f} % "
                f"le long de la lame",
            ]
            if abs(a_ox["pente_s_par_j_par_degc"]) > reglage.tolerance_s_par_j_par_degc:
                lignes.append(f"HORS tolérance ({reglage.tolerance_s_par_j_par_degc:.3f} s/j/°C)")
        else:
            lignes.append(f"Compensation impossible : {a_ox.get('raison', 'inconnue')}")

        lignes.append("")
        lignes.append("Sensibilité aux coefficients (épaisseur requise | pente si non recalée) :")
        for l in sensibilite:
            lignes.append(
                f"  {l['coefficient']:<14s} {l['t_moins_um']:.3f}..{l['t_plus_um']:.3f} µm | "
                f"{l['pente_non_recalee_moins']:+.3f}..{l['pente_non_recalee_plus']:+.3f} s/j/°C")
        quadrature = float(np.hypot.reduce(
            [max(abs(l["pente_non_recalee_moins"]), abs(l["pente_non_recalee_plus"]))
             for l in sensibilite]))
        lignes.append(
            f"  -> cumulées en quadrature : {quadrature:.3f} s/j/°C si l'épaisseur "
            f"d'oxyde n'est PAS recalée sur le procédé réel.")
        lignes.append(
            "  Les coefficients étant de littérature, l'épaisseur doit être calée "
            "par la mesure : le calcul donne la cible, pas la certitude.")
        self._dire("thermique", "\n".join(lignes))

    def _tracer_orbite_cg(self, diag_orbite: dict):
        """Orbite du centre de gravité : la grandeur que l'optimisation
        minimise réellement.

        Deux vues plutôt qu'une. La TRAJECTOIRE dans le plan dit l'essentiel —
        c'est son étendue, non sa position moyenne, qui crée l'erreur de marche
        en position verticale. Les deux COMPOSANTES en fonction de l'angle du
        balancier montrent les accidents locaux que la trajectoire seule
        masque : un repli, un point de rebroussement, une portion où le CG
        repart en arrière.
        """
        x = np.asarray(diag_orbite["x_cg_mm"], dtype=float) * 1000.0   # µm
        y = np.asarray(diag_orbite["y_cg_mm"], dtype=float) * 1000.0
        angles = np.asarray(diag_orbite["angles_deg"], dtype=float)

        self.ax_cg_xy.clear()
        self.ax_cg_barre.clear()
        points = self.ax_cg_xy.scatter(x, y, c=angles, cmap="coolwarm", s=8, zorder=3)
        self.ax_cg_xy.plot(x, y, "-", color=theme.gris(0.6), lw=0.6, zorder=2)
        self.fig_cg.colorbar(points, cax=self.ax_cg_barre, label="angle du balancier (°)")
        # Aspect égal mais limites ÉTIRABLES, et cadrage carré centré sur
        # l'orbite : une orbite quasi circulaire laissée à l'aspect strict se
        # réduisait à une bande de 274 px de large, le reste de la figure étant
        # perdu en marges.
        self.ax_cg_xy.set_aspect("equal", adjustable="datalim")
        cx, cy = 0.5 * (x.min() + x.max()), 0.5 * (y.min() + y.max())
        demi = max(x.max() - x.min(), y.max() - y.min(), 1e-6) * 0.6
        self.ax_cg_xy.set_xlim(cx - demi, cx + demi)
        self.ax_cg_xy.set_ylim(cy - demi, cy + demi)
        self.ax_cg_xy.set_xlabel("x du centre de gravité (µm)")
        self.ax_cg_xy.set_ylabel("y (µm)")
        self.ax_cg_xy.grid(True, linestyle=":", alpha=0.5)
        self.ax_cg_xy.set_title(
            f"Orbite du centre de gravité — étendue "
            f"{x.max()-x.min():.1f} x {y.max()-y.min():.1f} µm", fontsize=9)

        self.ax_cg_angle.clear()
        self.ax_cg_angle.plot(angles, x, "-", color="tab:blue", lw=1.0, label="x")
        self.ax_cg_angle.plot(angles, y, "-", color="tab:red", lw=1.0, label="y")
        self.ax_cg_angle.axvline(0, color=theme.encre(), lw=0.8)
        self.ax_cg_angle.set_xlabel("angle du balancier (°)")
        self.ax_cg_angle.set_ylabel("position du CG (µm)")
        self.ax_cg_angle.grid(True, linestyle=":", alpha=0.5)
        self.ax_cg_angle.legend(fontsize=8)
        self.fig_cg.tight_layout()
        self.canvas_cg.draw_idle()

    def _metriques_geometrie(self, geom, r_min: float = None):
        """(freq_hz, jeu_min_um, dispersion_s_j, iso_amplitude_basse_s_j,
        dvh_s_j) d'une géométrie, en UN seul balayage elastica (marche, jeu,
        isochronisme propre et DVH partagent le même).

        r_min est passé EXPLICITEMENT par l'appelant plutôt que déduit de la
        géométrie : deux géométries comparées entre elles (maillage de travail
        contre maillage fin) doivent partir de la même condition limite, sinon
        on compare deux calculs qui ne diffèrent pas que par le maillage.

        Le DVH est rendu à l'amplitude de référence, signé : c'est le chiffre
        que le COSC note, et il doit être vérifié convergé au même titre que la
        dispersion — il est plus petit qu'elle d'un facteur deux au moins, donc
        plus exposé au bruit de discrétisation."""
        if r_min is None:
            r_min = float(geom.rayon()[0])
        resultat = analyser_marches(
            geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
            amplitudes_deg=self._amplitudes_analyse(),
            reglage_echappement=self._reglage_echappement(),
            positions=positions_verticales(self.chk_cote_rouage.isChecked()),
            pas_balayage_deg=self.sb_pas_balayage.value(),
            jeu_goupilles_um=self._jeu_goupilles_um(),
        )
        return (frequence_propre(geom, self.spec.inertie_balancier),
                float(resultat["diag_orbite"]["jeu_mm"].min()) * 1000.0,
                resultat["dispersion_gravite"], resultat["iso_amplitude_basse"],
                resultat["dvh_reference"])

    def _geometrie_generee(self, titre: str) -> bool:
        """Le raffinement de maillage suppose une géométrie PARAMÉTRIQUE : un
        DXF importé n'a pas de vecteur de conception à reconstruire plus fin.
        Le test `geom.n != base.n_elements` détecte le cas où la géométrie
        affichée ne vient plus des paramètres courants."""
        if self.geom is None or self.base is None or self.geom.n != self.base.n_elements:
            QMessageBox.information(
                self, titre,
                "Disponible pour une géométrie générée (panneaux 2 à 4). Pour un DXF "
                "importé, l'import se fait déjà à 3000 éléments.",
            )
            return False
        return True

    def _base_raffinee(self, facteur_max: float = 4.0, plafond: int = 8000):
        """Mêmes paramètres de forme, maillage multiplié. Toutes les densités
        sont multipliées ENSEMBLE : ne raffiner qu'une zone changerait la
        répartition des éléments, donc la géométrie discrétisée, et l'écart
        mesuré ne serait plus imputable au seul maillage."""
        facteur = min(facteur_max, plafond / self.base.n_elements)
        return replace(
            self.base,
            n_elements_corps=int(self.base.n_elements_corps * facteur),
            grossmann_n_pas=int(self.base.grossmann_n_pas * facteur),
            phillips_n_pas=int(self.base.phillips_n_pas * facteur),
            raquette_transition_n_pas=int(self.base.raquette_transition_n_pas * facteur),
            raquette_arc_avant_n_pas=int(self.base.raquette_arc_avant_n_pas * facteur),
            raquette_piton_n_pas=int(self.base.raquette_piton_n_pas * facteur),
        )

    def _marches_6h_et_plat(self, geom, r_min, amplitudes):
        """(marche 6H, marche à plat, DVH) sur toute la plage d'amplitude.

        La position à plat n'a aucun terme de gravité : elle porte donc
        exactement les contributions indépendantes de la position
        (isochronisme propre, échappement, jeu de goupilles). Leur différence
        ne laisse subsister que la gravité en 6H — c'est le DVH, et c'est
        pourquoi il ne dépend que de l'orbite du centre de gravité."""
        resultat = analyser_marches(
            geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
            amplitudes_deg=amplitudes,
            reglage_echappement=self._reglage_echappement(),
            positions=positions_verticales(self.chk_cote_rouage.isChecked()),
            pas_balayage_deg=self.sb_pas_balayage.value(),
            jeu_goupilles_um=self._jeu_goupilles_um(),
        )
        totales = resultat["totales"]
        nom_plat = next(iter(POSITIONS_HORIZONTALES))
        return totales["6H"], totales[nom_plat], resultat["dvh"], resultat

    def on_marche_haute_resolution(self, silencieux: bool = False):
        """Marche en 6H et à plat sur toute la plage d'amplitude, au maillage
        de travail et à un maillage raffiné, plus leur différence (le DVH).

        Une seule amplitude ne suffit pas à juger le DVH : ce qui intéresse le
        réglage, c'est son ÉVOLUTION avec l'amplitude, et c'est justement elle
        que le bruit de discrétisation déforme en premier."""
        self._synchroniser_spec()
        if self.geom is None or self.base is None or self.geom.n != self.base.n_elements:
            # Géométrie importée : pas de vecteur de conception à reconstruire
            # plus fin. On le dit, sans interrompre un enchaînement.
            self._dire("dvh", "DVH à haute résolution : indisponible sur une géométrie "
                              "importée (pas de paramétrage à raffiner).")
            if not silencieux:
                self._geometrie_generee("Marche à haute résolution")
            return

        base_fine = self._base_raffinee()
        # Grille d'amplitude DENSIFIÉE à 10° : la grille d'affichage (souvent
        # 20°) suffit à situer un niveau, pas à lire une évolution.
        amplitudes = self._amplitudes_analyse()
        amplitudes = np.arange(float(amplitudes.min()), float(amplitudes.max()) + 1e-6, 10.0)
        r_min = float(self.geom.rayon()[0])

        self._dire("dvh", f"Calcul à {self.base.n_elements} puis "
                                f"{base_fine.n_elements} éléments...")
        self._occuper("dvh", True, "Calcul...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            six0, plat0, dvh0, res0 = self._marches_6h_et_plat(self.geom, r_min, amplitudes)
            geom_fine = construire_geometrie(base_fine,
                                             reglage_thermique=self._reglage_thermique())
            six1, plat1, dvh1, res1 = self._marches_6h_et_plat(geom_fine, r_min, amplitudes)
        except Exception as e:
            QMessageBox.critical(self, "Marche à haute résolution", str(e))
            return
        finally:
            QApplication.restoreOverrideCursor()
            self._occuper("dvh", False,
                                        "Marche 6H / à plat à haute résolution")

        n0, n1 = self.base.n_elements, base_fine.n_elements
        ax = self.ax_dvh_positions
        ax.clear()
        ax.plot(amplitudes, six0, "--", color="tab:red", lw=1.0, alpha=0.6,
                label=f"6H — {n0} éléments")
        ax.plot(amplitudes, six1, "-", color="tab:red", lw=1.6,
                label=f"6H — {n1} éléments")
        ax.plot(amplitudes, plat0, "--", color="tab:blue", lw=1.0, alpha=0.6,
                label=f"à plat — {n0} éléments")
        ax.plot(amplitudes, plat1, "-", color="tab:blue", lw=1.6,
                label=f"à plat — {n1} éléments")
        ax.set_ylabel("marche (s/j)")
        ax.set_title("Marche en position COSC 6H et à plat", fontsize=10)
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend(fontsize=8, ncol=2)

        ax2 = self.ax_dvh_ecart
        ax2.clear()
        ax2.plot(amplitudes, dvh0, "--", color=theme.gris(0.35), lw=1.0, alpha=0.7,
                 label=f"DVH — {n0} éléments")
        ax2.plot(amplitudes, dvh1, "-", color=theme.encre(), lw=1.8, label=f"DVH — {n1} éléments")
        ax2.fill_between(amplitudes, dvh0, dvh1, color="tab:orange", alpha=0.25,
                         label="incertitude de maillage")
        ax2.axhline(0, color=theme.encre(), lw=0.8)
        ax2.set_xlabel("amplitude (deg)")
        ax2.set_ylabel("DVH (s/j)")
        ax2.set_title("DVH = 6H − à plat", fontsize=10)
        ax2.grid(True, linestyle=":", alpha=0.5)
        ax2.legend(fontsize=8)
        self.fig_dvh.tight_layout()
        self.canvas_dvh.draw_idle()
        self.tabs.setCurrentWidget(self.canvas_dvh.parentWidget())

        # Le diagnostic compare l'incertitude de maillage à la VARIATION du DVH
        # sur la plage, et non à sa valeur : c'est l'évolution qu'on vient lire
        # sur ce graphe, et c'est elle qui peut être noyée alors même que le
        # niveau moyen paraît solide.
        etendue = float(np.max(dvh1) - np.min(dvh1))
        incertitude = float(np.max(np.abs(dvh1 - dvh0)))
        lignes = [
            f"DVH à {amplitudes[0]:.0f}° : {dvh0[0]:+.3f} → {dvh1[0]:+.3f} s/j ; "
            f"à {amplitudes[-1]:.0f}° : {dvh0[-1]:+.3f} → {dvh1[-1]:+.3f} s/j",
            f"Variation sur la plage : {etendue:.3f} s/j "
            f"(de {np.min(dvh1):+.3f} à {np.max(dvh1):+.3f})",
            f"Incertitude de maillage : {incertitude:.3f} s/j au pire point",
        ]
        if incertitude > etendue:
            lignes.append(
                "ATTENTION : l'incertitude de maillage DÉPASSE la variation du DVH sur "
                "toute la plage — la forme de la courbe n'est pas résolue, son évolution "
                "avec l'amplitude ne veut rien dire à ce maillage.")
        elif incertitude > 0.3 * etendue:
            lignes.append(
                "L'incertitude de maillage vaut plus du tiers de la variation : le sens "
                "de l'évolution est lisible, pas son ampleur.")
        else:
            lignes.append("La forme de la courbe est résolue à ce maillage.")
        self._dire("dvh", "\n".join(lignes))

    def on_verifier_haute_resolution(self):
        """Compare les métriques au maillage de travail et à un maillage 4x
        plus fin (plafonné à 8000 éléments), pour quantifier directement le
        biais de discrétisation sur LA géométrie active (cf. convergence.py)."""
        self._synchroniser_spec()
        if not self._geometrie_generee("Haute résolution"):
            return

        base_fine = self._base_raffinee()
        self._dire("haute_res", f"Calcul à {base_fine.n_elements} éléments...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            f0, j0, d0, i0, v0 = self._metriques_geometrie(self.geom, r_min=self.base.r_min)
            f1, j1, d1, i1, v1 = self._metriques_geometrie(
                construire_geometrie(base_fine, reglage_thermique=self._reglage_thermique()),
                r_min=self.base.r_min)
        finally:
            QApplication.restoreOverrideCursor()

        reglage = self._reglage_echappement()
        lignes = [
            f"Maillage de travail ({self.base.n_elements}) -> fin ({base_fine.n_elements}) :",
            f"Fréquence : {f0:.4f} -> {f1:.4f} Hz",
            f"Jeu minimal : {j0:.1f} -> {j1:.1f} µm",
            # L'écart est donné EN RELATIF : les trois grandeurs de marche
            # subissent à peu près la même erreur absolue, mais elles n'ont pas
            # le même ordre de grandeur. Le DVH, au moins deux fois plus petit
            # que la dispersion, est le plus exposé.
            f"Dispersion : {d0:.2f} -> {d1:.2f} s/j ({_relatif(d1 - d0, d1)})",
            f"DVH (COSC) : {v0:+.2f} -> {v1:+.2f} s/j ({_relatif(v1 - v0, v1)})",
            f"Isochronisme propre à {reglage.amplitude_basse_deg:.0f}° : "
            f"{i0:+.2f} -> {i1:+.2f} s/j ({_relatif(i1 - i0, i1)}, "
            f"cible {reglage.retard_s_par_j:+.2f})",
        ]
        if abs(d1 - d0) > 0.5 * abs(d1):
            lignes.append(
                f"ATTENTION : l'écart de maillage ({abs(d1-d0):.2f} s/j) est du même "
                f"ordre que la dispersion elle-même — ce chiffre n'est pas résolu.")
        if abs(v1 - v0) > 0.5 * abs(v1):
            lignes.append(
                f"ATTENTION : l'écart de maillage sur le DVH ({abs(v1-v0):.2f} s/j) est "
                f"du même ordre que le DVH lui-même — son signe comme sa valeur sont "
                f"douteux à ce maillage.")
        if abs(i1 - i0) > reglage.tolerance_s_par_j:
            lignes.append(
                f"ATTENTION : l'isochronisme propre bouge de {abs(i1-i0):.2f} s/j avec "
                f"le maillage, soit plus que la tolérance de compensation demandée "
                f"({reglage.tolerance_s_par_j:.2f} s/j) — l'optimiseur contraindrait "
                f"du bruit.")
        self._dire("haute_res", "\n".join(lignes))

    def on_afficher_jeu_minimal(self):
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Jeu minimal", "Générez ou chargez d'abord une géométrie.")
            return
        r_min = float(self.geom.rayon()[0])
        angles_deg, _, _, jeu_mm, _ = balayage_angulaire(
            self.geom, r_min, self.spec.amplitude_max, self.sb_pas_balayage.value()
        )
        idx = int(np.argmin(jeu_mm))
        angle_pire = float(angles_deg[idx])
        Xe, Ye, Xi, Yi, _, _ = resoudre_elastica(self.geom, angle_pire, r_min)

        self.ax_jeu_min.clear()
        ext0, inte0 = self.geom.contour_flancs()
        self.ax_jeu_min.plot(ext0[:, 0], ext0[:, 1], "-", color=theme.gris(0.83), lw=0.7, label="repos (0°)")
        self.ax_jeu_min.plot(inte0[:, 0], inte0[:, 1], "-", color=theme.gris(0.83), lw=0.7)
        self.ax_jeu_min.plot(Xe, Ye, "-", color="tab:red", lw=0.9, label=f"déformé à {angle_pire:+.0f}°")
        self.ax_jeu_min.plot(Xi, Yi, "-", color="tab:red", lw=0.9)
        self.ax_jeu_min.plot(0, 0, "+", color=theme.encre(), ms=10)
        self.ax_jeu_min.set_aspect("equal")
        self.ax_jeu_min.set_xlabel("x (mm)")
        self.ax_jeu_min.set_ylabel("y (mm)")
        self.ax_jeu_min.set_title(f"Jeu minimal = {jeu_mm[idx]*1000:.1f} µm à l'angle {angle_pire:+.0f}°")
        self.ax_jeu_min.legend()
        self.canvas_jeu_min.draw_idle()
        self.tabs.setCurrentWidget(self.canvas_jeu_min)

        self._dernier_jeu_min_um = float(jeu_mm[idx] * 1000.0)
        self._dernier_jeu_min_angle = angle_pire
        self._afficher_resultats_numeriques()

    def on_calculer_animation(self):
        self._lancer_calcul_animation(chemin_gif_a_exporter=None)

    def _lancer_calcul_animation(self, chemin_gif_a_exporter):
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Animation", "Générez ou chargez d'abord une géométrie.")
            return
        self._chemin_gif_a_exporter = chemin_gif_a_exporter
        self._occuper("anim", True, "Calcul en cours...")
        self._dire("gif", f"Calcul des frames en cours (amplitude ±{self.spec.amplitude_max:.0f}°)...")
        self._thread = AnimationThread(self.geom, float(self.geom.rayon()[0]), self.spec.amplitude_max, n_frames=41)
        self._thread.signal_termine.connect(self.on_animation_calculee)
        self._thread.start()

    def on_animation_calculee(self, angles_deg, frames):
        self._occuper("anim", False)
        self._anim_charger(angles_deg, frames)
        self._dire("gif", f"Animation calculée ({len(frames)} frames) -- voir l'onglet \"Animation\".")
        # Le jeu minimal se lit sur le même déploiement : l'enchaîner évite de
        # refaire un balayage pour une vue qui en découle.
        self.on_afficher_jeu_minimal()
        self.tabs.setCurrentWidget(self.canvas_anim.parentWidget())
        if getattr(self, "_chemin_gif_a_exporter", None):
            self._exporter_gif_depuis_frames(self._chemin_gif_a_exporter)
            self._chemin_gif_a_exporter = None

    def _exporter_gif_depuis_frames(self, chemin: str):
        self._dire("gif", "Encodage du GIF...")
        try:
            sauvegarder_gif_deploiement(self._anim_angles, self._anim_frames, chemin, aller_retour=True)
        except Exception as e:
            self._dire("gif", "(échec de l'enregistrement)")
            QMessageBox.critical(self, "Erreur d'enregistrement du GIF", str(e))
            return
        self._dire_ecourte("gif", f"Enregistré : {chemin}")

    # ------------------------------------------------------------------
    # Logique -- Onglet Animation (lecteur interactif)
    # ------------------------------------------------------------------
    def _anim_charger(self, angles_deg, frames):
        self._anim_angles = angles_deg
        self._anim_frames = frames
        self._anim_sequence = indices_aller_retour(len(frames))

        # Cadrage FIXE, calculé une fois sur l'enveloppe de TOUTES les frames
        # (le spiral le plus déployé, à l'amplitude max) -- sans ça, chaque
        # frame se recadre sur sa propre étendue (ax.clear() efface aussi les
        # limites), ce qui donne l'illusion trompeuse que le centre de
        # rotation se déplace pendant la lecture (diagnostic 2026-09-25).
        toutes_x = np.concatenate([np.concatenate([f[0], f[2]]) for f in frames])
        toutes_y = np.concatenate([np.concatenate([f[1], f[3]]) for f in frames])
        marge = 0.05 * max(toutes_x.max() - toutes_x.min(), toutes_y.max() - toutes_y.min(), 1e-6)
        self._anim_xlim = (toutes_x.min() - marge, toutes_x.max() + marge)
        self._anim_ylim = (toutes_y.min() - marge, toutes_y.max() + marge)

        self.slider_anim.setEnabled(True)
        self.slider_anim.setRange(0, len(self._anim_sequence) - 1)
        self.btn_anim_play.setEnabled(True)
        self.btn_anim_prev.setEnabled(True)
        self.btn_anim_next.setEnabled(True)
        self.slider_anim.setValue(0)
        self._anim_afficher_position(0)

    def _anim_afficher_position(self, pos: int):
        if not self._anim_sequence:
            return
        i = self._anim_sequence[pos]
        Xe, Ye, Xi, Yi = self._anim_frames[i]
        # le zoom de l'utilisateur tient pendant la lecture
        xlim, ylim = (self.canvas_anim.navigation.vue_utilisateur(self.ax_anim)
                      or (self._anim_xlim, self._anim_ylim))
        self.ax_anim.clear()
        ext0, inte0 = self.geom.contour_flancs() if self.geom is not None else (None, None)
        if ext0 is not None:
            self.ax_anim.plot(ext0[:, 0], ext0[:, 1], "-", color=theme.gris(0.83), lw=0.6)
            self.ax_anim.plot(inte0[:, 0], inte0[:, 1], "-", color=theme.gris(0.83), lw=0.6)
        self.ax_anim.plot(Xe, Ye, "-", color="tab:blue", lw=0.9)
        self.ax_anim.plot(Xi, Yi, "-", color="tab:blue", lw=0.9)
        self.ax_anim.plot(0, 0, "+", color=theme.encre(), ms=10)
        self.ax_anim.set_xlim(*xlim)
        self.ax_anim.set_ylim(*ylim)
        self.ax_anim.set_aspect("equal")
        self.ax_anim.set_xlabel("x (mm)")
        self.ax_anim.set_ylabel("y (mm)")
        self.canvas_anim.draw_idle()
        self.label_anim_angle.setText(f"angle = {self._anim_angles[i]:+.0f}°")

    def _anim_pas_manuel(self, delta: int):
        if self._anim_timer.isActive():
            self._anim_basculer_lecture()  # un pas manuel interrompt la lecture auto
        nouvelle_pos = int(np.clip(self.slider_anim.value() + delta, 0, self.slider_anim.maximum()))
        self.slider_anim.setValue(nouvelle_pos)

    def _anim_avancer_auto(self):
        pos = (self.slider_anim.value() + 1) % (self.slider_anim.maximum() + 1)
        self.slider_anim.setValue(pos)

    def _anim_basculer_lecture(self):
        if self._anim_timer.isActive():
            self._anim_timer.stop()
            self.btn_anim_play.setText("▶ Lecture")
        else:
            self._anim_mettre_a_jour_intervalle()
            self._anim_timer.start()
            self.btn_anim_play.setText("⏸ Pause")

    def _anim_mettre_a_jour_intervalle(self):
        self._anim_timer.setInterval(int(1000 / self.sb_anim_fps.value()))

    # --- Anciennes expériences (format v4-v7) ---------------------------
    def _champ_widget(self, cle):
        """Widget de l'interface qui porte un champ de la spec, s'il y en a un."""
        return {"angle_echappement_deg": self.sb_angle_echappement,
                "collet_cote_ancre": self.chk_collet_cote_ancre,
                "sens_trigonometrique": self.chk_sens_trigo,
                "vue_cote_rouage": self.chk_cote_rouage,
                "plan_cristallin": self.cb_plan,
                "orientation_cristal_deg": self.sb_orientation_cristal}.get(cle)

    def _marquer_a_valider(self, widget, libelle):
        """Encadre un champ tant que l'utilisateur ne l'a pas confirmé.

        Le confirmer, c'est y toucher : changer la valeur, valider la saisie
        (Entrée ou sortie du champ), cliquer la case, choisir dans la liste --
        même pour garder la valeur proposée."""
        self._a_valider[widget] = libelle
        widget.setStyleSheet(f"border: 2px solid {theme.ACC2.name()};")
        widget.setToolTip(f"À valider : {libelle} — absent de l'ancien fichier, "
                          "valeur proposée par défaut.")
        confirmer = lambda *_a, w=widget: self._confirmer_champ(w)
        if isinstance(widget, QCheckBox):
            widget.clicked.connect(confirmer)
        elif isinstance(widget, QComboBox):
            widget.activated.connect(confirmer)
        else:
            widget.valueChanged.connect(confirmer)
            widget.editingFinished.connect(confirmer)

    def _confirmer_champ(self, widget):
        if self._a_valider.pop(widget, None) is None:
            return
        widget.setStyleSheet("")
        widget.setToolTip("")
        if self._a_valider:
            self.statusBar().showMessage(
                "Encore à valider : " + ", ".join(self._a_valider.values()))
        else:
            self.statusBar().showMessage("Tous les champs repris de l'ancien fichier sont validés.")

    def _effacer_marques(self):
        for w in list(self._a_valider):
            w.setStyleSheet("")
            w.setToolTip("")
        self._a_valider.clear()

    def _charger_ancienne_experience(self, chemin):
        """Reprend une expérience des anciens programmes (v4-v7) : le spiral
        est reconstruit par l'ancien générateur et repris comme une géométrie
        importée ; les champs que l'ancien format ignorait sont mis en
        évidence jusqu'à ce que l'utilisateur les ait validés."""
        try:
            spec, geom, champs, diag = ancien_format.reprendre_experience(
                chemin, self._lire_spec_depuis_formulaire())
        except Exception as e:
            QMessageBox.critical(self, "Ancien fichier illisible", str(e))
            return
        self._effacer_marques()
        self._reinitialiser_analyses()
        self.spec = spec
        for sb, v in ((self.sb_r_min, spec.r_min), (self.sb_r_max, spec.r_max),
                      (self.sb_n_tours, spec.n_tours_reels),
                      (self.sb_n_tours_fictifs, spec.n_tours_fictifs),
                      (self.sb_hauteur, spec.hauteur_h), (self.sb_freq, spec.freq_cible),
                      (self.sb_inertie, _inertie_affichee(spec.inertie_balancier)),
                      (self.sb_jeu_min, spec.jeu_min_spires),
                      (self.sb_amplitude_max, spec.amplitude_max),
                      (self.sb_amp_max, spec.amplitude_max)):
            sb.setValue(v)
        self.geom = geom
        self._ancienne = (chemin, geom)
        _texte_ecourte(self.label_dxf, "ancien format : " + chemin.split("/")[-1])
        self._tracer_geometrie()
        a_valider = [c for c in champs if c.statut == ancien_format.A_VALIDER]
        for c in a_valider:
            w = self._champ_widget(c.cle)
            # un champ grisé (orientation sur un wafer isotrope) est sans effet :
            # il n'y a rien à y valider
            if w is not None and w.isEnabled():
                self._marquer_a_valider(w, c.libelle)
        f = frequence_propre(geom, spec.inertie_balancier)
        QMessageBox.information(self, "Ancien fichier repris",
                                self._texte_reprise(chemin, champs, diag, f, spec))
        self.statusBar().showMessage("Encore à valider : " + ", ".join(self._a_valider.values()))

    @staticmethod
    def _texte_reprise(chemin, champs, diag, freq, spec):
        def ligne(c):
            v = c.valeur
            if isinstance(v, float):
                v = f"{v:g}"
            elif isinstance(v, bool):
                v = "oui" if v else "non"
            elif isinstance(v, dict):
                v = "…"
            note = f" <i>— {c.note}</i>" if c.note else ""
            return f"<li><b>{c.libelle}</b> : {v}{note}</li>"
        orange = theme.ACC2.name()
        blocs = []
        for statut, titre in ((ancien_format.A_VALIDER,
                               f"<span style='color:{orange}'>À valider</span> "
                               "(absents de l'ancien format, encadrés dans le panneau)"),
                              (ancien_format.CONVERTI, "Convertis"),
                              (ancien_format.REPRIS, "Repris tels quels"),
                              (ancien_format.IGNORE, "Non repris")):
            lignes = [ligne(c) for c in champs if c.statut == statut]
            if lignes:
                blocs.append(f"<p><b>{titre}</b></p><ul>{''.join(lignes)}</ul>")
        controle = (f"<p>Contrôle : fréquence recalculée par le modèle actuel "
                    f"{freq:.4f} Hz pour une cible de {spec.freq_cible:g} Hz ; "
                    f"épaisseurs mises à l'échelle par l'ancien programme "
                    f"(facteur {diag['facteur_epaisseur']:.4f}).</p>")
        fin = ("<p>Changer l'angle d'échappement, le côté du collet, le sens ou le wafer "
               "repose aussitôt le spiral. Le spiral repris s'analyse au panneau Résultats ; "
               "il ne se ré-optimise pas avec les panneaux de forme.</p>")
        return (f"<p><b>{chemin.split('/')[-1]}</b> est un fichier des anciens programmes "
                f"(v4 à v7) : le spiral a été reconstruit exactement comme ils le dessinaient.</p>"
                + "".join(blocs) + controle + fin)

    def _reposer_ancienne_experience(self, *_a):
        """Repose le spiral repris quand son orientation ou le wafer change :
        l'ancien fichier ne les fixait pas, ils viennent de l'interface."""
        if self._ancienne is None or self.geom is not self._ancienne[1]:
            return
        try:
            _spec, geom, _c, _d = ancien_format.reprendre_experience(
                self._ancienne[0], self._lire_spec_depuis_formulaire())
        except Exception:
            return
        self._reinitialiser_analyses()
        self.geom = geom
        self._ancienne = (self._ancienne[0], geom)
        self._tracer_geometrie()

    def on_charger_configuration(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Charger une configuration complète", "", "JSON (*.json)")
        if not chemin:
            return
        try:
            with open(chemin, encoding="utf-8") as f:
                ancien = ancien_format.est_experience(json.load(f))
        except Exception:
            ancien = False
        if ancien:
            self._charger_ancienne_experience(chemin)
            return
        self._effacer_marques()
        try:
            spec, base, bornes, hyperparams, echappement = charger_configuration(chemin)
        except Exception as e:
            QMessageBox.critical(self, "Erreur de chargement", str(e))
            return
        self.spec, self.base, self.bornes, self.hyperparams = spec, base, bornes, hyperparams

        self.sb_r_min.setValue(spec.r_min); self.sb_r_max.setValue(spec.r_max)
        self.sb_n_tours.setValue(spec.n_tours_reels); self.sb_n_tours_fictifs.setValue(spec.n_tours_fictifs)
        self.sb_hauteur.setValue(spec.hauteur_h); self.sb_freq.setValue(spec.freq_cible)
        self.sb_inertie.setValue(_inertie_affichee(spec.inertie_balancier))
        self.sb_jeu_min.setValue(spec.jeu_min_spires)
        self.cb_plan.setCurrentIndex(max(0, self.cb_plan.findData(spec.plan_cristallin)))
        self.sb_orientation_cristal.setValue(spec.orientation_cristal_deg)
        self._maj_orientation_cristal()
        self.sb_amplitude_max.setValue(spec.amplitude_max); self.sb_amp_max.setValue(spec.amplitude_max)
        self.chk_cote_rouage.setChecked(spec.vue_cote_rouage)

        self.sb_gross_angle.setValue(base.grossmann_angle_deg); self.sb_gross_n_pas.setValue(base.grossmann_n_pas)
        self.sb_gross_n_points.setValue(base.grossmann_n_points)
        self.sb_phil_angle.setValue(base.phillips_angle_deg); self.sb_phil_n_pas.setValue(base.phillips_n_pas)
        self.sb_phil_n_points.setValue(base.phillips_n_points)
        self.sb_n_corps.setValue(base.n_elements_corps)
        self.sb_n_noeuds_ep.setValue(len(base.epaisseur_noeuds)); self.sb_n_noeuds_dr.setValue(len(base.delta_r_noeuds))

        self.sb_pop.setValue(hyperparams.pop_size); self.sb_gen.setValue(hyperparams.n_gen)
        self.sb_ga_seed.setValue(hyperparams.seed)

        self.sb_ech_retard.setValue(echappement.retard_s_par_j)
        self.sb_ech_amp_basse.setValue(echappement.amplitude_basse_deg)
        self.sb_ech_amp_ref.setValue(echappement.amplitude_reference_deg)
        self.sb_ech_tolerance.setValue(echappement.tolerance_s_par_j)
        self.chk_echappement.setChecked(echappement.actif)
        i = self.cb_variante.findData(echappement.variante_visee)
        if i >= 0:
            self.cb_variante.setCurrentIndex(i)

        self._reinitialiser_analyses()
        geom, diag = construire_geometrie(
            self.base, return_diagnostics=True,
            reglage_thermique=self._reglage_thermique())
        self.geom = geom
        self.label_dxf.setText("(configuration complète chargée)")
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()
