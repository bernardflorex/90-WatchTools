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
from dataclasses import replace
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
)
from PySide6.QtCore import Qt, QTimer

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
from .isochronisme import analyser_marches, POSITIONS_HORIZONTALES
from .echappement import ReglageEchappement
from .thermique import ReglageThermique, analyser_thermique, etude_sensibilite
from . import dxf_io
from .optimisation import vecteur_depuis_params, params_depuis_vecteur, PAS_DISPERSION_DEG_DEFAUT
from .threads import SondageThread, OptimisationThread, AnimationThread
from .rendu import sauvegarder_gif_deploiement, indices_aller_retour

_LARGEUR_MAX_CHAMP = 110


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
        self.sondage_resultats = None
        self.solutions_pareto = []  # liste de dicts {x, dispersion, jeu}
        self._thread = None

        splitter = QSplitter(Qt.Horizontal)
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
        contenu = QWidget()
        layout = QVBoxLayout(contenu)
        layout.addWidget(self._panneau_origine())
        layout.addWidget(self._panneau_parametres())
        layout.addWidget(self._panneau_sondage())
        layout.addWidget(self._panneau_optimisation())
        layout.addWidget(self._panneau_resultats())
        layout.addStretch(1)

        # Pas de largeur MAXIMALE ici : elle empêchait d'élargir le panneau en
        # tirant la poignée du splitter, donc de lire un texte qui dépasse
        # (tableaux, résultats numériques). Seule une largeur minimale est
        # imposée, la largeur initiale étant fixée par le splitter lui-même
        # (voir __init__ : setSizes).
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(contenu)
        scroll.setMinimumWidth(300)
        return scroll

    def _onglet_avec_barre(self, canvas: FigureCanvas) -> QWidget:
        """Enveloppe un canvas avec la barre d'outils matplotlib standard
        (zoom rectangle, pan, retour vue précédente, sauvegarde image...)."""
        conteneur = QWidget()
        v = QVBoxLayout(conteneur)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(NavigationToolbar(canvas, conteneur), 0)
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
        # réellement. Deux vues -- la trajectoire dans le plan (c'est son
        # ÉTENDUE qui crée l'erreur de marche en position verticale) et les
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
        # calcul, et la détacher après un clear() de son axes parent lève une
        # erreur dans matplotlib.
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
        self.sb_inertie = _spinbox(1e-11, 1e-7, self.spec.inertie_balancier, 1e-11, 13)
        self.sb_inertie.setToolTip("Moment d'inertie du balancier (kg.m^2) -- utilisé avec freq_cible pour calculer la raideur/épaisseur requise.")
        self.sb_jeu_min = _spinbox(0.001, 0.1, self.spec.jeu_min_spires, 0.001, 4)
        self.sb_jeu_min.setToolTip("Jeu minimal requis entre spires adjacentes (mm), sous déformation à amplitude max. Contrainte dure de l'Optimisation (pas du Sondage).")
        self.sb_amplitude_max = _spinbox(100.0, 360.0, self.spec.amplitude_max, 5.0, 1)
        self.sb_amplitude_max.setToolTip("Amplitude d'oscillation maximale du balancier (deg) -- pire cas utilisé pour le jeu inter-spires, le GIF et la vue \"jeu minimal\".")
        form.addRow("r_min (mm)", self.sb_r_min)
        form.addRow("r_max (mm)", self.sb_r_max)
        form.addRow("n_tours_reels", self.sb_n_tours)
        form.addRow("n_tours_fictifs", self.sb_n_tours_fictifs)
        form.addRow("hauteur_h (mm)", self.sb_hauteur)
        form.addRow("freq_cible (Hz)", self.sb_freq)
        form.addRow("inertie_balancier (kg.m^2)", self.sb_inertie)
        form.addRow("jeu_min_spires (mm)", self.sb_jeu_min)
        form.addRow("amplitude_max (deg)", self.sb_amplitude_max)
        layout.addLayout(form)

        ligne_btn = QHBoxLayout()
        btn_init = QPushButton("Initialiser Archimède + Grossmann/Phillips")
        btn_init.setToolTip("Construit un spiral neuf à partir de la spec ci-dessus + des paramètres de forme du panneau 2.")
        btn_init.clicked.connect(self.on_init_archimede)
        ligne_btn.addWidget(btn_init)
        btn_save = QPushButton("Sauver spec...")
        btn_save.setToolTip("Sauve UNIQUEMENT la spec du mouvement ci-dessus (pas les paramètres de forme -- voir \"Sauver configuration complète\" au panneau 5).")
        btn_save.clicked.connect(self.on_sauver_spec)
        ligne_btn.addWidget(btn_save)
        btn_load = QPushButton("Charger spec...")
        btn_load.setToolTip("Charge une spec de mouvement sauvée précédemment.")
        btn_load.clicked.connect(self.on_charger_spec)
        ligne_btn.addWidget(btn_load)
        layout.addLayout(ligne_btn)

        return box

    # --- 2. Paramètres ----------------------------------------------------
    def _panneau_parametres(self) -> QGroupBox:
        box = QGroupBox("2. Paramètres de forme (terminaisons intérieure et extérieure)")
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
        self.sb_n_corps.setToolTip("Nombre d'éléments discrets du corps actif (spirale à pas constant, entre les deux courbes terminales).")
        self.sb_n_noeuds_ep = _intbox(2, 60, 20)
        self.sb_n_noeuds_ep.setToolTip("Nombre de points de contrôle du profil d'épaisseur libre le long du spiral (plus il y en a, plus le profil peut varier localement).")
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
        box = QGroupBox("Terminaison extérieure à raquette")
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
        self.sb_raq_n_pas_arc = _intbox(5, 2000, 61)
        self.sb_raq_n_pas_arc.setToolTip(
            "Nombre d'éléments discrets sur l'arc. Ramené au nombre IMPAIR immédiatement "
            "supérieur s'il est pair, pour que le point de comptage tombe exactement au "
            "milieu et que la raquette dispose du même parcours des deux côtés. Il fixe "
            "la finesse de positionnement de la raquette, pas la précision du calcul."
        )
        self.sb_raq_piton_angle = _spinbox(0.0, 270.0, 30.0, 5.0, 1)
        self.sb_raq_piton_angle.setToolTip(
            "Étendue angulaire (deg) de l'arc AU-DELÀ du parcours des goupilles, "
            "jusqu'au piton. C'est le paramètre le plus important du montage à "
            "raquette : c'est la flexion de ce tronçon libre qui déplace la lame "
            "latéralement au droit des goupilles, donc qui lui permet de venir en "
            "contact avec l'une d'elles. Trop court, la lame flotte en permanence "
            "entre les goupilles et la raquette n'a AUCUNE autorité. Trop long, le "
            "retard dû au jeu devient important. L'excursion latérale varie comme le "
            "CARRÉ de cette longueur."
        )
        self.sb_raq_piton_n_pas = _intbox(0, 2000, 100)
        self.sb_raq_piton_n_pas.setToolTip("Nombre d'éléments discrets du maillage sur l'arc entre goupilles et piton.")
        self.sb_raq_jeu = _spinbox(0.0, 200.0, 20.0, 1.0, 1)
        self.sb_raq_jeu.setToolTip(
            "Jeu LIBRE total entre les goupilles, épaisseur de lame déjà déduite : la "
            "lame peut donc s'écarter de ±jeu/2 avant de toucher.\n\n"
            "0 = piton fixe au point de comptage : aucune erreur de jeu, aucune "
            "dépendance en amplitude. Ce n'est pas un cas particulier traité à part, "
            "c'est la limite continue du même modèle.\n\n"
            "Au-delà, la lame traverse le jeu à chaque alternance ; pendant la "
            "traversée le spiral n'est plus tenu aux goupilles mais au piton, donc "
            "plus long et plus souple : la montre retarde, d'autant plus que "
            "l'amplitude est faible."
        )
        form.addRow("Plage de réglage ± (s/j)", self.sb_raq_plage)
        form.addRow("Piton -- angle au-delà (deg)", self.sb_raq_piton_angle)
        form.addRow("Piton -- n_pas (maillage)", self.sb_raq_piton_n_pas)
        form.addRow("Jeu libre entre goupilles (µm)", self.sb_raq_jeu)
        form.addRow("Transition -- angle (deg)", self.sb_raq_trans_angle)
        form.addRow("Transition -- n_points (Bézier)", self.sb_raq_n_points)
        form.addRow("Transition -- n_pas (maillage)", self.sb_raq_trans_n_pas)
        form.addRow("Arc -- n_pas (maillage)", self.sb_raq_n_pas_arc)
        v.addLayout(form)
        self._champs_raquette = [self.sb_raq_plage, self.sb_raq_piton_angle,
                                  self.sb_raq_piton_n_pas, self.sb_raq_jeu,
                                  self.sb_raq_trans_angle, self.sb_raq_n_points,
                                  self.sb_raq_trans_n_pas, self.sb_raq_n_pas_arc]
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
        box = QGroupBox("Compensation thermique (oxyde de silicium)")
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

        self.btn_thermique = QPushButton("Analyser la tenue en température")
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
        box = QGroupBox("3. Sondage (point de départ le moins mauvais)")
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
        self.table_sondage.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
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

        layout.addWidget(self._panneau_echappement())

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

        self.table_pareto = QTableWidget(0, 3)
        self.table_pareto.setToolTip("Solutions non-dominées du front de Pareto. Cliquer une ligne charge cette solution comme géométrie active.")
        self.table_pareto.setHorizontalHeaderLabels(["Dispersion (s/j)", "Jeu min (µm)", "Compens. (s/j)"])
        self.table_pareto.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_pareto.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_pareto.itemSelectionChanged.connect(self.on_selectionner_solution)
        self.table_pareto.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
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
        box = QGroupBox("Compensation de l'erreur d'échappement")
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
    def _panneau_resultats(self) -> QGroupBox:
        box = QGroupBox("5. Résultats et export")
        box.setToolTip(
            "Toutes les actions ci-dessous travaillent sur la GÉOMÉTRIE ACTIVE, quelle "
            "que soit son origine : DXF chargé au panneau 1, géométrie générée au "
            "panneau 2, ou solution sélectionnée au panneau 3 ou 4. Après un simple "
            "chargement DXF, c'est ICI qu'il faut venir directement -- pas besoin des "
            "panneaux 1-4, qui ne servent qu'à CONSTRUIRE ou OPTIMISER une géométrie."
        )
        layout = QVBoxLayout(box)

        btn_marche = QPushButton("Calculer la marche diurne (géométrie active)")
        btn_marche.setToolTip(
            "Fonctionne quelle que soit l'origine de la géométrie active (DXF importé, "
            "génération paramétrique, ou solution issue du Sondage/de l'Optimisation) -- "
            "pas besoin de repasser par les panneaux 1-4."
        )
        btn_marche.clicked.connect(self.on_calculer_marche)
        layout.addWidget(btn_marche)

        btn_jeu_min = QPushButton("Afficher la géométrie au jeu inter-spires minimal")
        btn_jeu_min.setToolTip("Cherche l'angle (dans ±amplitude_max) où le jeu inter-spires est le plus faible, et affiche la forme déformée à cet angle.")
        btn_jeu_min.clicked.connect(self.on_afficher_jeu_minimal)
        layout.addWidget(btn_jeu_min)

        self.btn_calculer_anim = QPushButton("Calculer l'animation du déploiement (mini <-> maxi)")
        self.btn_calculer_anim.setToolTip(
            "Calcule 41 formes déformées entre -amplitude_max et +amplitude_max "
            "(quelques secondes à quelques dizaines de secondes selon la finesse du "
            "maillage). Alimente le lecteur interactif de l'onglet \"Animation\" "
            "(curseur, lecture en va-et-vient, vitesse réglable, pas-à-pas) -- "
            "l'export GIF ci-dessous réutilise ces frames sans recalcul."
        )
        self.btn_calculer_anim.clicked.connect(self.on_calculer_animation)
        layout.addWidget(self.btn_calculer_anim)

        self.btn_exporter_gif = QPushButton("Exporter l'animation en GIF...")
        self.btn_exporter_gif.setToolTip(
            "Calcule l'animation si ce n'est pas déjà fait, puis l'encode en GIF "
            "animé en va-et-vient (comme l'oscillation réelle du balancier, pas un "
            "aller simple qui boucle en sautant)."
        )
        self.btn_exporter_gif.clicked.connect(self.on_exporter_gif)
        layout.addWidget(self.btn_exporter_gif)
        self.label_gif = QLabel("")
        layout.addWidget(self.label_gif)

        ligne = QHBoxLayout()
        btn_export = QPushButton("Exporter en DXF...")
        btn_export.setToolTip("Exporte les deux flancs (offset ±épaisseur/2) de la géométrie active en LWPOLYLINE.")
        btn_export.clicked.connect(self.on_exporter_dxf)
        ligne.addWidget(btn_export)
        btn_valider = QPushButton("Valider l'export (réimport DXF)")
        btn_valider.setToolTip(
            "Réexporte la géométrie active en DXF, la RÉIMPORTE comme un fichier "
            "externe, et compare : détecte une perte de précision à l'export/import "
            "(arrondi, ré-échantillonnage) -- PAS une validation de la conception "
            "elle-même (déjà couverte par Sondage/Optimisation). Nécessite d'avoir "
            "cliqué \"Exporter en DXF...\" d'abord."
        )
        btn_valider.clicked.connect(self.on_valider_boucle_fermee)
        ligne.addWidget(btn_valider)
        layout.addLayout(ligne)

        self.label_validation = QLabel("")
        self.label_validation.setWordWrap(True)
        self.label_validation.setStyleSheet("font-family: monospace;")
        layout.addWidget(self.label_validation)

        box_resultats_num = QGroupBox("Résultats numériques (géométrie active)")
        v_num = QVBoxLayout(box_resultats_num)
        self.label_resultats_num = QLabel("(aucune géométrie)")
        # Pas de police monospace ici (plus large par caractère qu'une police
        # proportionnelle -> forçait le retour à la ligne même sur des lignes
        # courtes dans ce panneau étroit) + hauteur minimale garantie : le
        # texte se retrouvait tronqué (diagnostic 2026-09-25).
        self.label_resultats_num.setWordWrap(True)
        self.label_resultats_num.setMinimumHeight(110)
        v_num.addWidget(self.label_resultats_num)

        self.btn_verif_haute_res = QPushButton("Vérifier à haute résolution (maillage x4)")
        self.btn_verif_haute_res.setToolTip(
            "Recalcule fréquence, jeu minimal et dispersion de marche de la géométrie "
            "générée (panneau 2/3/4) avec un maillage 4x plus fin, et affiche l'écart "
            "avec le maillage de travail. D'après l'étude de convergence (2026-09-25), "
            "un maillage grossier surestime le risque (jeu plus pessimiste, dispersion "
            "plus élevée que la valeur convergée) : c'est le sens conservateur, mais "
            "à vérifier avant de retenir un candidat final. Pour un DXF importé, "
            "l'import se fait déjà à 3000 éléments."
        )
        self.btn_verif_haute_res.clicked.connect(self.on_verifier_haute_resolution)
        v_num.addWidget(self.btn_verif_haute_res)
        self.label_verif_haute_res = QLabel("")
        self.label_verif_haute_res.setWordWrap(True)
        self.label_verif_haute_res.setMinimumHeight(90)
        v_num.addWidget(self.label_verif_haute_res)
        layout.addWidget(box_resultats_num)

        ligne_config = QHBoxLayout()
        btn_save_config = QPushButton("Sauver configuration complète...")
        btn_save_config.setToolTip("Sauve en UN fichier JSON : spec du mouvement + tous les paramètres de forme (Grossmann/Phillips, profils) + bornes + hyperparamètres GA.")
        btn_save_config.clicked.connect(self.on_sauver_configuration)
        ligne_config.addWidget(btn_save_config)
        btn_load_config = QPushButton("Charger configuration complète...")
        btn_load_config.setToolTip("Recharge une configuration complète sauvée précédemment et régénère la géométrie correspondante.")
        btn_load_config.clicked.connect(self.on_charger_configuration)
        ligne_config.addWidget(btn_load_config)
        layout.addLayout(ligne_config)

        return box

    # ------------------------------------------------------------------
    # Logique -- Origine
    # ------------------------------------------------------------------
    def _lire_spec_depuis_formulaire(self) -> SpecMouvement:
        return SpecMouvement(
            r_min=self.sb_r_min.value(), r_max=self.sb_r_max.value(),
            n_tours_reels=self.sb_n_tours.value(), n_tours_fictifs=self.sb_n_tours_fictifs.value(),
            hauteur_h=self.sb_hauteur.value(), freq_cible=self.sb_freq.value(),
            inertie_balancier=self.sb_inertie.value(), jeu_min_spires=self.sb_jeu_min.value(),
            amplitude_max=self.sb_amplitude_max.value(),
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
            "raquette_n_pas_arc": self.sb_raq_n_pas_arc.value(),
            "raquette_piton_angle_deg": self.sb_raq_piton_angle.value(),
            "raquette_piton_n_pas": self.sb_raq_piton_n_pas.value(),
            "raquette_jeu_um": self.sb_raq_jeu.value(),
            "oxydation_active": self.chk_oxydation.isChecked() if hasattr(self, "chk_oxydation") else False,
            "oxyde_impose_um": self.sb_ox_impose.value() if hasattr(self, "sb_ox_impose") else 0.0,
        }

    def on_init_archimede(self):
        self._rafraichir_base_depuis_spec()
        try:
            geom, diag = construire_geometrie(self.base, return_diagnostics=True,
                                              reglage_thermique=self._reglage_thermique())
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
            geom = dxf_io.charger_dxf(chemin, n_elements=3000, hauteur_h=self.sb_hauteur.value())
        except Exception as e:
            QMessageBox.critical(self, "Erreur d'import", str(e))
            return
        self.geom = geom
        self.label_dxf.setText(chemin.split("/")[-1])
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

    def on_sauver_spec(self):
        chemin, _ = QFileDialog.getSaveFileName(self, "Sauver la spec", "", "JSON (*.json)")
        if chemin:
            self._lire_spec_depuis_formulaire().sauvegarder(_assurer_extension(chemin, ".json"))

    def on_charger_spec(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Charger une spec", "", "JSON (*.json)")
        if not chemin:
            return
        self.spec = SpecMouvement.charger(chemin)
        self.sb_r_min.setValue(self.spec.r_min)
        self.sb_r_max.setValue(self.spec.r_max)
        self.sb_n_tours.setValue(self.spec.n_tours_reels)
        self.sb_n_tours_fictifs.setValue(self.spec.n_tours_fictifs)
        self.sb_hauteur.setValue(self.spec.hauteur_h)
        self.sb_freq.setValue(self.spec.freq_cible)
        self.sb_inertie.setValue(self.spec.inertie_balancier)
        self.sb_jeu_min.setValue(self.spec.jeu_min_spires)
        self.sb_amplitude_max.setValue(self.spec.amplitude_max)
        # Synchronise la fenêtre d'analyse (sondage/optimisation) sur la spec
        # chargée -- sinon on peut se retrouver à optimiser pour une
        # amplitude différente de celle réellement visée par le mouvement.
        self.sb_amp_max.setValue(self.spec.amplitude_max)

    # ------------------------------------------------------------------
    # Logique -- Paramètres
    # ------------------------------------------------------------------
    def on_generer_geometrie(self):
        self._rafraichir_base_depuis_spec()
        try:
            geom, diag = construire_geometrie(self.base, return_diagnostics=True,
                                              reglage_thermique=self._reglage_thermique())
        except ValueError as e:
            QMessageBox.warning(self, "Erreur", str(e))
            return
        self.geom = geom
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()

    def _afficher_diagnostics(self, diag: dict):
        raquette = diag.get("raquette_active", False)
        nom_ext = "courbe en S" if raquette else "Phillips"
        lignes = [
            f"geometrie_valide      : {diag['geometrie_valide']}",
            f"grossmann_monotone    : {diag['grossmann_monotone']}",
            f"{nom_ext + ' monotone':22s}: {diag['phillips_monotone']}",
            f"phillips_plausibilite : {diag['phillips_ratio_plausibilite']:.2f}  (repère ~1.0, pas une contrainte)",
            f"freq avant scaling    : {diag['freq_avant_scaling_hz']:.3f} Hz  (beta={diag['beta_scaling']:.3f})",
        ]
        if diag.get("oxydation_active"):
            lignes.append(
                f"oxyde thermique       : {diag['epaisseur_oxyde_um']:.3f} µm  "
                f"(lame +{diag['croissance_lame_um']:.2f} µm à la gravure)")
        if raquette:
            lignes.append(
                f"arc de raquette       : {diag['raquette_arc_deg']:.2f}°  "
                f"({diag['raquette_arc_mm']*1000:.0f} µm), longueur morte "
                f"{diag['longueur_morte_mm']*1000:.0f} µm"
            )
        self.label_diag.setText("\n".join(lignes))

    def _tracer_geometrie(self):
        # Une nouvelle géométrie invalide le jeu minimal et la dispersion de
        # marche affichés (calculés pour l'ANCIENNE géométrie) -- on les
        # efface plutôt que de laisser afficher une valeur périmée.
        self._dernier_jeu_min_um = None
        self._dernier_dispersion = None
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
            self.ax_geom.plot(ext[:, 0], ext[:, 1], "k-", lw=0.8)
            self.ax_geom.plot(inte[:, 0], inte[:, 1], "k-", lw=0.8)
            self.ax_geom.plot(0, 0, "k+", ms=10)
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
                n_p = (self.base.raquette_transition_n_pas + (self.base.raquette_n_pas_arc | 1)
                       if self.base.raquette_active else self.base.phillips_n_pas)
            else:
                n_g = n_p = max(5, self.geom.n // 20)
            marge_g = n_g // 2
            marge_p = n_p // 2

            i_fin_g = min(n_g + marge_g, self.geom.n)
            self.ax_zoom_grossmann.plot(ext[:i_fin_g, 0], ext[:i_fin_g, 1], "k.-", ms=2, lw=0.8)
            self.ax_zoom_grossmann.plot(inte[:i_fin_g, 0], inte[:i_fin_g, 1], "k.-", ms=2, lw=0.8)
            self.ax_zoom_grossmann.set_aspect("equal")
            self.ax_zoom_grossmann.set_title("Zoom collet (Grossmann)", fontsize=9)
            self.ax_zoom_grossmann.tick_params(labelsize=7)

            i_debut_p = max(self.geom.n - n_p - marge_p, 0)
            self.ax_zoom_phillips.plot(ext[i_debut_p:, 0], ext[i_debut_p:, 1], "k.-", ms=2, lw=0.8)
            self.ax_zoom_phillips.plot(inte[i_debut_p:, 0], inte[i_debut_p:, 1], "k.-", ms=2, lw=0.8)
            self.ax_zoom_phillips.set_aspect("equal")
            self.ax_zoom_phillips.tick_params(labelsize=7)

            if self.geom.a_partie_morte:
                # La portion morte (au-delà des goupilles) est tracée en gris :
                # elle est bien usinée et occupe de la place, mais ne travaille
                # pas. Sans ce repère, rien ne distingue à l'écran la longueur
                # active de la longueur totale.
                na = self.geom.n_actif
                for ax, lw in ((self.ax_geom, 1.6), (self.ax_zoom_phillips, 2.2)):
                    ax.plot(ext[na - 1:, 0], ext[na - 1:, 1], "-", color="0.65", lw=lw)
                    ax.plot(inte[na - 1:, 0], inte[na - 1:, 1], "-", color="0.65", lw=lw)
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
            self.label_resultats_num.setText("(aucune géométrie)")
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
            n_arc = self.sb_raq_n_pas_arc.value() if hasattr(self, "sb_raq_n_pas_arc") else 61
            rapide, lente, etendue = marche_raquette(
                self.geom, self.spec.inertie_balancier, n_arc)
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

        if self._dernier_force_pivot is not None:
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
        self.label_resultats_num.setText("\n".join(lignes))

    def _tracer_orbite_cg(self, diag_orbite: dict):
        """Trajectoire du centre de gravité pendant le déploiement.

        C'est la grandeur que l'optimisation minimise : en position verticale,
        le poids du spiral exerce sur l'axe un couple proportionnel au
        DÉPLACEMENT du centre de gravité, et c'est sa variation au cours de
        l'oscillation -- pas sa position moyenne -- qui fait l'erreur de
        marche. Un spiral parfaitement réglé aurait une orbite réduite à un
        point."""
        angles = diag_orbite["angles_deg"]
        xc, yc = diag_orbite["x_cg_mm"] * 1000.0, diag_orbite["y_cg_mm"] * 1000.0

        self.ax_cg_xy.clear()
        self.ax_cg_angle.clear()
        self.ax_cg_xy.set_anchor("W")
        trace = self.ax_cg_xy.scatter(xc, yc, c=angles, cmap="coolwarm", s=14, zorder=3)
        self.ax_cg_xy.plot(xc, yc, "-", color="0.6", lw=0.8, zorder=2)
        self.ax_cg_xy.plot(xc.mean(), yc.mean(), "k+", ms=12, zorder=4)
        etendue = float(np.hypot(xc - xc.mean(), yc - yc.mean()).max())
        self.ax_cg_xy.set_aspect("equal")
        self.ax_cg_xy.set_xlabel("x du CG (µm)")
        self.ax_cg_xy.set_ylabel("y du CG (µm)")
        self.ax_cg_xy.grid(True, linestyle=":", alpha=0.5)
        self.ax_cg_xy.set_title(
            f"Orbite du CG -- rayon max {etendue:.1f} µm\n"
            f"décentrage moyen {np.hypot(xc.mean(), yc.mean()):.1f} µm", fontsize=9)
        self.ax_cg_barre.clear()
        barre = self.fig_cg.colorbar(trace, cax=self.ax_cg_barre)
        barre.set_label("angle du balancier (°)", fontsize=8)
        barre.ax.tick_params(labelsize=7)

        self.ax_cg_angle.plot(angles, xc - xc.mean(), "o-", ms=3, label="x (écart à la moyenne)")
        self.ax_cg_angle.plot(angles, yc - yc.mean(), "s-", ms=3, label="y (écart à la moyenne)")
        self.ax_cg_angle.axvline(0, color="k", lw=0.8)
        self.ax_cg_angle.set_xlabel("Angle du balancier (°)")
        self.ax_cg_angle.set_ylabel("Déplacement du CG (µm)")
        self.ax_cg_angle.grid(True, linestyle=":", alpha=0.5)
        self.ax_cg_angle.legend(fontsize=8)
        self.ax_cg_angle.set_title(
            f"{len(angles)} positions de déploiement calculées", fontsize=9)
        self.fig_cg.tight_layout()
        self.canvas_cg.draw_idle()

    def _tracer_epaisseur(self):
        self.ax_epaisseur.clear()
        if self.geom is not None:
            if self.geom.epaisseur_oxyde > 0.0:
                self.ax_epaisseur.plot(self.geom.s, self.geom.epaisseur_finie * 1000.0, "-",
                                        label="lame finie (oxyde compris)")
                self.ax_epaisseur.plot(self.geom.s, self.geom.epaisseur * 1000.0, "--",
                                        color="0.5", label="cote de gravure")
                self.ax_epaisseur.legend(fontsize=8)
            else:
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
            bouton.setStyleSheet("background-color: #f5c542; font-weight: bold;")
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
            reglage_thermique=self._reglage_thermique(),
        )
        self._thread.signal_termine.connect(self.on_sondage_termine)
        self._thread.start()

    def on_sondage_termine(self, resultats: dict):
        self._definir_bouton_actif(self.btn_sondage, False)
        self.sondage_resultats = resultats
        n_ok = int(np.sum(resultats["jeu_mm"] >= self.spec.jeu_min_spires))
        self.label_sondage.setText(
            f"{len(resultats['jeu_mm'])} tirages -- {n_ok} atteignent déjà le jeu requis "
            f"({self.spec.jeu_min_spires*1000:.1f} µm), à titre indicatif seulement.\n"
            f"Les {len(resultats['classement'])} meilleurs (classés du moins mauvais) sont "
            f"listés ci-dessous -- sélectionnez-en un comme base, ou passez à l'Optimisation "
            f"(le meilleur est utilisé par défaut)."
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
        geom, diag = construire_geometrie(self.base, return_diagnostics=True,
                                              reglage_thermique=self._reglage_thermique())
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
            pas_balayage_deg=self.sb_pas_balayage.value(),
            reglage_echappement=self._reglage_echappement(),
            reglage_thermique=self._reglage_thermique(),
        )
        self._thread.signal_progres.connect(self.on_optimisation_progres)
        self._thread.signal_termine.connect(self.on_optimisation_terminee)
        self._thread.start()

    def on_optimisation_progres(self, gen: int, dispersion: float, jeu: float):
        self.label_optim.setText(
            f"Génération {gen}/{self.hyperparams.n_gen} -- "
            f"meilleure dispersion pop. : {dispersion:.2f} s/j, jeu : {jeu*1000:.1f} µm"
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
            self.label_optim.setText(f"Terminée -- {len(self.solutions_pareto)} solution(s) non-dominée(s).")

        self.table_pareto.setRowCount(len(self.solutions_pareto))
        for i, sol in enumerate(self.solutions_pareto):
            self.table_pareto.setItem(i, 0, QTableWidgetItem(f"{sol['dispersion']:.2f}"))
            self.table_pareto.setItem(i, 1, QTableWidgetItem(f"{sol['jeu']*1000:.1f}"))
            ecart = sol.get("ecart_compensation")
            self.table_pareto.setItem(
                i, 2, QTableWidgetItem("--" if ecart is None else f"|{ecart:.2f}|"))
        self._tracer_pareto()

    def _tracer_pareto(self):
        self.ax_pareto.clear()
        if self.solutions_pareto:
            disp = [s["dispersion"] for s in self.solutions_pareto]
            jeu_um = [s["jeu"] * 1000 for s in self.solutions_pareto]
            self.ax_pareto.scatter(disp, jeu_um)
            self.ax_pareto.set_xlabel("Dispersion de marche (s/j)")
            self.ax_pareto.set_ylabel("Jeu min inter-spires (µm)")
            self.ax_pareto.axhline(self.spec.jeu_min_spires * 1000, color="r", ls="--", lw=0.8, label="jeu min requis")
            self.ax_pareto.legend()
        self.canvas_pareto.draw_idle()

    def on_selectionner_solution(self):
        lignes = self.table_pareto.selectionModel().selectedRows()
        if not lignes or not self.solutions_pareto:
            return
        idx = lignes[0].row()
        sol = self.solutions_pareto[idx]
        self.base = params_depuis_vecteur(sol["x"], self.base)
        self.geom = construire_geometrie(self.base, reglage_thermique=self._reglage_thermique())
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
        self.ax_marche.plot(amplitudes, resultat["isochronisme"], ":", color="0.35",
                            lw=1.5, label="isochronisme propre")
        self.ax_marche.plot(amplitudes, resultat["echappement"], ":", color="tab:red",
                            lw=1.5, label="échappement")
        if resultat["diag_jeu_goupilles"].get("actif"):
            self.ax_marche.plot(amplitudes, resultat["jeu_goupilles"], ":", color="tab:purple",
                                lw=1.5, label="jeu de goupilles")
        self.ax_marche.axvline(reglage.amplitude_basse_deg, color="0.7", ls="-.", lw=0.8)
        self.ax_marche.set_xlabel("Amplitude (deg)")
        self.ax_marche.set_ylabel("Marche diurne (s/j)")
        self.ax_marche.set_title(
            f"Compensation à {reglage.amplitude_basse_deg:.0f}° : "
            f"spiral {resultat['iso_amplitude_basse']:+.2f} s/j "
            f"contre {-reglage.retard_s_par_j:+.2f} s/j d'échappement "
            f"(écart {resultat['ecart_compensation']:+.2f} s/j)", fontsize=9)
        self.ax_marche.axhline(0, color="k", lw=0.8)
        self.ax_marche.grid(True, linestyle=":", alpha=0.5)
        self.ax_marche.legend(fontsize=8, ncol=2)
        self.fig_marche.tight_layout()
        self.canvas_marche.draw_idle()
        self.tabs.setCurrentWidget(self.canvas_marche.parentWidget())

        self._dernier_dispersion = resultat["dispersion_gravite"]
        self._dernier_iso_basse = resultat["iso_amplitude_basse"]
        self._dernier_ecart_compensation = resultat["ecart_compensation"]
        self._dernier_diag_iso = resultat["diag_isochronisme"]
        self._dernier_jeu_basse = resultat["jeu_amplitude_basse"]
        self._dernier_diag_jeu = resultat["diag_jeu_goupilles"]
        self._tracer_orbite_cg(resultat["diag_orbite"])
        self._dernier_force_pivot = axe_force_pivot(
            self.geom, r_min, self.spec.amplitude_max)
        self._afficher_resultats_numeriques()

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
        self.ax_thermique.axhspan(-0.6 * etendue / 2, 0.6 * etendue / 2, color="0.9",
                                   zorder=0, label="bande COSC (0.6 s/j/°C)")
        self.ax_thermique.axhline(0, color="k", lw=0.8)
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
        self.label_validation.setText("\n".join(lignes))

    def on_exporter_dxf(self):
        if self.geom is None:
            QMessageBox.information(self, "Export", "Générez ou chargez d'abord une géométrie.")
            return
        chemin, _ = QFileDialog.getSaveFileName(self, "Exporter en DXF", "", "DXF (*.dxf)")
        if chemin:
            chemin = _assurer_extension(chemin, ".dxf")
            # Deux fichiers quand la lame est oxydée : la cote de GRAVURE part
            # à la DRIE, la pièce FINIE sert à vérifier et à mesurer. Les
            # confondre donnerait une lame trop épaisse de ~4 µm.
            chemins = dxf_io.exporter_dxf_paire(self.geom, chemin)
            self._dernier_export = chemins[-1]
            QMessageBox.information(self, "Export", "Exporté :\n" + "\n".join(chemins))

    def on_valider_boucle_fermee(self):
        """Réexporte la géométrie active en DXF (déjà fait par "Exporter en
        DXF..."), la réimporte comme un fichier externe quelconque, et
        compare : ça vérifie que le pipeline d'export/import (ray-casting,
        ré-échantillonnage, cf. dxf_io.py) ne dégrade pas la géométrie -- une
        vérification de FIDÉLITÉ NUMÉRIQUE, pas une validation physique
        (la physique elle-même est déjà validée par Sondage/Optimisation)."""
        self._synchroniser_spec()
        chemin = getattr(self, "_dernier_export", None)
        if chemin is None or self.geom is None:
            QMessageBox.information(self, "Validation", "Exportez d'abord un DXF.")
            return
        try:
            geom_reimporte = dxf_io.charger_dxf(chemin, n_elements=self.geom.n, hauteur_h=self.geom.hauteur)
        except Exception as e:
            QMessageBox.critical(self, "Erreur de réimport", str(e))
            return

        r_avant, r_apres = self.geom.rayon(), geom_reimporte.rayon()
        ecart_r_min = abs(r_avant.min() - r_apres.min()) * 1000
        ecart_r_max = abs(r_avant.max() - r_apres.max()) * 1000

        amplitude_max = self.spec.amplitude_max
        pas = self.sb_pas_balayage.value()
        _, _, _, jeu_avant, _ = balayage_angulaire(self.geom, float(r_avant[0]), amplitude_max, pas)
        _, _, _, jeu_apres, _ = balayage_angulaire(geom_reimporte, float(r_apres[0]), amplitude_max, pas)

        self.label_validation.setText(
            f"Écart r_min après réimport : {ecart_r_min:.2f} µm\n"
            f"Écart r_max après réimport : {ecart_r_max:.2f} µm\n"
            f"Jeu min avant export : {jeu_avant.min()*1000:.1f} µm\n"
            f"Jeu min après réimport : {jeu_apres.min()*1000:.1f} µm"
        )

    def _metriques_geometrie(self, geom, r_min: float = None):
        """(freq_hz, jeu_min_um, dispersion_s_j, iso_amplitude_basse_s_j) d'une
        géométrie, en UN seul balayage elastica (marche, jeu et isochronisme
        propre partagent le même).

        r_min : rayon de collet imposé comme condition limite. Passé
        EXPLICITEMENT par l'appelant plutôt que déduit de la géométrie, pour
        que deux géométries comparées entre elles (maillage de travail vs
        maillage fin) utilisent forcément la même condition limite -- sans
        quoi on compare deux calculs qui ne diffèrent pas que par le maillage.
        Par défaut : le rayon effectif du premier point."""
        if r_min is None:
            r_min = float(geom.rayon()[0])
        resultat = analyser_marches(
            geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
            amplitudes_deg=self._amplitudes_analyse(),
            reglage_echappement=self._reglage_echappement(),
            pas_balayage_deg=self.sb_pas_balayage.value(),
            jeu_goupilles_um=self._jeu_goupilles_um(),
        )
        return (frequence_propre(geom, self.spec.inertie_balancier),
                float(resultat["diag_orbite"]["jeu_mm"].min()) * 1000.0,
                resultat["dispersion_gravite"], resultat["iso_amplitude_basse"])

    def on_verifier_haute_resolution(self):
        """Compare les métriques au maillage de travail et à un maillage 4x
        plus fin (plafonné à 8000 éléments), pour quantifier directement le
        biais de discrétisation sur LA géométrie active (cf. convergence.py)."""
        self._synchroniser_spec()
        if self.geom is None or self.base is None or self.geom.n != self.base.n_elements:
            QMessageBox.information(
                self, "Haute résolution",
                "Disponible pour une géométrie générée (panneaux 2 à 4). Pour un DXF "
                "importé, l'import se fait déjà à 3000 éléments.",
            )
            return

        facteur = min(4.0, 8000.0 / self.base.n_elements)
        base_fine = replace(
            self.base,
            n_elements_corps=int(self.base.n_elements_corps * facteur),
            grossmann_n_pas=int(self.base.grossmann_n_pas * facteur),
            phillips_n_pas=int(self.base.phillips_n_pas * facteur),
        )
        self.label_verif_haute_res.setText(f"Calcul à {base_fine.n_elements} éléments...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            f0, j0, d0, i0 = self._metriques_geometrie(self.geom, r_min=self.base.r_min)
            f1, j1, d1, i1 = self._metriques_geometrie(
                construire_geometrie(base_fine, reglage_thermique=self._reglage_thermique()),
                r_min=self.base.r_min)
        finally:
            QApplication.restoreOverrideCursor()

        # Incertitude de discrétisation : l'écart travail -> fin donne un ORDRE
        # DE GRANDEUR de l'erreur restante. S'il est comparable à la valeur
        # elle-même, le chiffre n'est pas résolu et ne doit pas être lu comme
        # un résultat physique (cas constaté : une dispersion annoncée à
        # 0.57 s/j alors que l'incertitude était de l'ordre de 1 s/j).
        reglage = self._reglage_echappement()
        incertitude = abs(d1 - d0)
        incertitude_iso = abs(i1 - i0)
        lignes = [
            f"Maillage de travail ({self.base.n_elements}) -> fin ({base_fine.n_elements}) :",
            f"Fréquence : {f0:.4f} -> {f1:.4f} Hz",
            f"Jeu minimal : {j0:.1f} -> {j1:.1f} µm",
            f"Dispersion : {d0:.2f} -> {d1:.2f} s/j",
            f"Isochronisme propre à {reglage.amplitude_basse_deg:.0f}° : "
            f"{i0:+.2f} -> {i1:+.2f} s/j (cible {reglage.retard_s_par_j:+.2f})",
        ]
        if incertitude > 0.5 * abs(d1):
            lignes.append(
                f"ATTENTION : l'écart de maillage ({incertitude:.2f} s/j) est du même "
                f"ordre que la dispersion elle-même -- ce chiffre n'est pas résolu. "
                f"Raffiner le maillage de travail avant de comparer des solutions."
            )
        # La compensation d'échappement se juge contre sa TOLÉRANCE, pas contre
        # sa propre valeur : viser +/-1 s/j avec un maillage qui bouge de 2 s/j
        # revient à contraindre du bruit de discrétisation.
        if incertitude_iso > reglage.tolerance_s_par_j:
            lignes.append(
                f"ATTENTION : l'isochronisme propre bouge de {incertitude_iso:.2f} s/j "
                f"avec le maillage, soit plus que la tolérance de compensation demandée "
                f"({reglage.tolerance_s_par_j:.2f} s/j) -- élargir la tolérance ou "
                f"raffiner le maillage, sinon l'optimiseur contraint du bruit."
            )
        self.label_verif_haute_res.setText("\n".join(lignes))

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
        self.ax_jeu_min.plot(ext0[:, 0], ext0[:, 1], "-", color="lightgray", lw=0.7, label="repos (0°)")
        self.ax_jeu_min.plot(inte0[:, 0], inte0[:, 1], "-", color="lightgray", lw=0.7)
        self.ax_jeu_min.plot(Xe, Ye, "r-", lw=0.9, label=f"déformé à {angle_pire:+.0f}°")
        self.ax_jeu_min.plot(Xi, Yi, "r-", lw=0.9)
        self.ax_jeu_min.plot(0, 0, "k+", ms=10)
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

    def on_exporter_gif(self):
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "GIF", "Générez ou chargez d'abord une géométrie.")
            return
        chemin, _ = QFileDialog.getSaveFileName(self, "Enregistrer le GIF", "deploiement.gif", "GIF (*.gif)")
        if not chemin:
            self.label_gif.setText("(annulé -- aucun chemin choisi)")
            return
        chemin = _assurer_extension(chemin, ".gif")
        if self._anim_frames is not None:
            self._exporter_gif_depuis_frames(chemin)
        else:
            self._lancer_calcul_animation(chemin_gif_a_exporter=chemin)

    def _lancer_calcul_animation(self, chemin_gif_a_exporter):
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Animation", "Générez ou chargez d'abord une géométrie.")
            return
        self._chemin_gif_a_exporter = chemin_gif_a_exporter
        self._definir_bouton_actif(self.btn_calculer_anim, True, "Calcul en cours...")
        self._definir_bouton_actif(self.btn_exporter_gif, True)
        self.label_gif.setText(f"Calcul des frames en cours (amplitude ±{self.spec.amplitude_max:.0f}°)...")
        self._thread = AnimationThread(self.geom, float(self.geom.rayon()[0]), self.spec.amplitude_max, n_frames=41)
        self._thread.signal_termine.connect(self.on_animation_calculee)
        self._thread.start()

    def on_animation_calculee(self, angles_deg, frames):
        self._definir_bouton_actif(self.btn_calculer_anim, False)
        self._definir_bouton_actif(self.btn_exporter_gif, False)
        self._anim_charger(angles_deg, frames)
        self.label_gif.setText(f"Animation calculée ({len(frames)} frames) -- voir l'onglet \"Animation\".")
        self.tabs.setCurrentWidget(self.canvas_anim.parentWidget())
        if getattr(self, "_chemin_gif_a_exporter", None):
            self._exporter_gif_depuis_frames(self._chemin_gif_a_exporter)
            self._chemin_gif_a_exporter = None

    def _exporter_gif_depuis_frames(self, chemin: str):
        self.label_gif.setText("Encodage du GIF...")
        try:
            sauvegarder_gif_deploiement(self._anim_angles, self._anim_frames, chemin, aller_retour=True)
        except Exception as e:
            self.label_gif.setText("(échec de l'enregistrement)")
            QMessageBox.critical(self, "Erreur d'enregistrement du GIF", str(e))
            return
        self.label_gif.setText(f"Enregistré : {chemin}")

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
        self.ax_anim.clear()
        ext0, inte0 = self.geom.contour_flancs() if self.geom is not None else (None, None)
        if ext0 is not None:
            self.ax_anim.plot(ext0[:, 0], ext0[:, 1], "-", color="lightgray", lw=0.6)
            self.ax_anim.plot(inte0[:, 0], inte0[:, 1], "-", color="lightgray", lw=0.6)
        self.ax_anim.plot(Xe, Ye, "b-", lw=0.9)
        self.ax_anim.plot(Xi, Yi, "b-", lw=0.9)
        self.ax_anim.plot(0, 0, "k+", ms=10)
        self.ax_anim.set_xlim(*self._anim_xlim)
        self.ax_anim.set_ylim(*self._anim_ylim)
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

    def on_sauver_configuration(self):
        # Rafraîchit spec ET base depuis les champs affichés (pas seulement
        # spec) : on veut sauver ce qui est VISIBLE à l'écran, y compris si
        # l'utilisateur a changé un angle Grossmann/Phillips sans avoir
        # ensuite cliqué "Générer géométrie".
        self._rafraichir_base_depuis_spec()
        chemin, _ = QFileDialog.getSaveFileName(self, "Sauver la configuration complète", "", "JSON (*.json)")
        if chemin:
            sauvegarder_configuration(_assurer_extension(chemin, ".json"), self.spec, self.base,
                                       self.bornes, self.hyperparams, self._reglage_echappement())

    def on_charger_configuration(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Charger une configuration complète", "", "JSON (*.json)")
        if not chemin:
            return
        try:
            spec, base, bornes, hyperparams, echappement = charger_configuration(chemin)
        except Exception as e:
            QMessageBox.critical(self, "Erreur de chargement", str(e))
            return
        self.spec, self.base, self.bornes, self.hyperparams = spec, base, bornes, hyperparams

        self.sb_r_min.setValue(spec.r_min); self.sb_r_max.setValue(spec.r_max)
        self.sb_n_tours.setValue(spec.n_tours_reels); self.sb_n_tours_fictifs.setValue(spec.n_tours_fictifs)
        self.sb_hauteur.setValue(spec.hauteur_h); self.sb_freq.setValue(spec.freq_cible)
        self.sb_inertie.setValue(spec.inertie_balancier); self.sb_jeu_min.setValue(spec.jeu_min_spires)
        self.sb_amplitude_max.setValue(spec.amplitude_max); self.sb_amp_max.setValue(spec.amplitude_max)

        self.sb_gross_angle.setValue(base.grossmann_angle_deg); self.sb_gross_n_pas.setValue(base.grossmann_n_pas)
        self.sb_gross_n_points.setValue(base.grossmann_n_points)
        self.sb_phil_angle.setValue(base.phillips_angle_deg); self.sb_phil_n_pas.setValue(base.phillips_n_pas)
        self.sb_phil_n_points.setValue(base.phillips_n_points)
        self.sb_n_corps.setValue(base.n_elements_corps)
        self.chk_raquette.setChecked(base.raquette_active)
        self.sb_raq_plage.setValue(base.raquette_plage_s_par_j)
        self.sb_raq_trans_angle.setValue(base.raquette_transition_angle_deg)
        self.sb_raq_n_points.setValue(base.raquette_n_points)
        self.sb_raq_trans_n_pas.setValue(base.raquette_transition_n_pas)
        self.sb_raq_n_pas_arc.setValue(base.raquette_n_pas_arc)
        self.sb_raq_piton_angle.setValue(base.raquette_piton_angle_deg)
        self.sb_raq_piton_n_pas.setValue(base.raquette_piton_n_pas)
        self.sb_raq_jeu.setValue(base.raquette_jeu_um)
        self._maj_etat_raquette()
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

        geom, diag = construire_geometrie(self.base, return_diagnostics=True,
                                              reglage_thermique=self._reglage_thermique())
        self.geom = geom
        self.label_dxf.setText("(configuration complète chargée)")
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()
