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
    4. Optimisation  : NSGA-II, objectif = vraie dispersion de marche.
    5. Résultats     : marche diurne par position, export DXF, validation
                       par réimport (boucle fermée).

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
    QAbstractItemView, QSlider, QApplication, QHeaderView,
)
from PySide6.QtCore import Qt, QTimer

from .parametrage import ParametresSpiral, construire_geometrie, epaisseur_theorique
from .presets import (
    SpecMouvement, BornesOptimisation, HyperparametresGA,
    sauvegarder_configuration, charger_configuration,
)
from .physique import resoudre_elastica, balayage_angulaire, frequence_propre
from .constantes import RHO_SILICIUM
from .isochronisme import marche_positions_verticales, critere_dispersion_marche
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
        box = QGroupBox("2. Paramètres de forme (Grossmann intérieur / Phillips extérieur)")
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

        self.table_sondage = QTableWidget(0, 2)
        self.table_sondage.setHorizontalHeaderLabels(["Dispersion (s/j)", "Jeu min (µm)"])
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

        self.btn_optimisation = QPushButton("Lancer l'optimisation")
        self.btn_optimisation.setToolTip(
            "Objectif réel minimisé : dispersion de marche entre positions verticales. "
            "Contrainte réelle : jeu >= jeu_min_spires (panneau 1). Démarre depuis le "
            "meilleur point du Sondage si disponible, sinon depuis la géométrie active."
        )
        self.btn_optimisation.clicked.connect(self.on_lancer_optimisation)
        layout.addWidget(self.btn_optimisation)

        self.label_optim = QLabel("(pas encore lancée)")
        layout.addWidget(self.label_optim)

        self.table_pareto = QTableWidget(0, 2)
        self.table_pareto.setToolTip("Solutions non-dominées du front de Pareto. Cliquer une ligne charge cette solution comme géométrie active.")
        self.table_pareto.setHorizontalHeaderLabels(["Dispersion (s/j)", "Jeu min (µm)"])
        self.table_pareto.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_pareto.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_pareto.itemSelectionChanged.connect(self.on_selectionner_solution)
        self.table_pareto.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table_pareto.setMaximumHeight(160)
        layout.addWidget(self.table_pareto)

        return box

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
        )

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
            geom, diag = construire_geometrie(self.base, return_diagnostics=True)
        except ValueError as e:
            QMessageBox.warning(self, "Erreur", str(e))
            return
        self.geom = geom
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()

    def _afficher_diagnostics(self, diag: dict):
        lignes = [
            f"geometrie_valide      : {diag['geometrie_valide']}",
            f"grossmann_monotone    : {diag['grossmann_monotone']}",
            f"phillips_monotone     : {diag['phillips_monotone']}",
            f"phillips_plausibilite : {diag['phillips_ratio_plausibilite']:.2f}  (repère ~1.0, pas une contrainte)",
            f"freq avant scaling    : {diag['freq_avant_scaling_hz']:.3f} Hz  (beta={diag['beta_scaling']:.3f})",
        ]
        self.label_diag.setText("\n".join(lignes))

    def _tracer_geometrie(self):
        # Une nouvelle géométrie invalide le jeu minimal et la dispersion de
        # marche affichés (calculés pour l'ANCIENNE géométrie) -- on les
        # efface plutôt que de laisser afficher une valeur périmée.
        self._dernier_jeu_min_um = None
        self._dernier_dispersion = None

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
                n_g, n_p = self.base.grossmann_n_pas, self.base.phillips_n_pas
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
            self.ax_zoom_phillips.set_title("Zoom piton (Phillips)", fontsize=9)
            self.ax_zoom_phillips.tick_params(labelsize=7)

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
        if self._dernier_jeu_min_um is not None:
            lignes.append(f"Jeu minimal : {self._dernier_jeu_min_um:.1f} µm à {self._dernier_jeu_min_angle:+.0f}°")
        else:
            lignes.append("Jeu minimal : non calculé")
        if self._dernier_dispersion is not None:
            lignes.append(f"Dispersion de marche : {self._dernier_dispersion:.2f} s/j")
        else:
            lignes.append("Dispersion de marche : non calculée")
        self.label_resultats_num.setText("\n".join(lignes))

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
            bouton.setStyleSheet("background-color: #f5c542; font-weight: bold;")
        else:
            texte_repos = bouton.property("_texte_repos")
            if texte_repos:
                bouton.setText(texte_repos)
            bouton.setEnabled(True)
            bouton.setStyleSheet("")

    def on_lancer_sondage(self):
        self._synchroniser_spec()
        self._definir_bouton_actif(self.btn_sondage, True, "Sondage en cours...")
        self.label_sondage.setText("Sondage en cours...")
        self._thread = SondageThread(
            self.base, self.bornes, self._amplitudes_analyse(),
            n_echantillons=self.sb_sondage_n.value(), seed=self.sb_sondage_seed.value(),
            pas_balayage_deg=self.sb_pas_balayage.value(),
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
            self.solutions_pareto = [
                {"x": X[i], "dispersion": float(F[i, 0]), "jeu": float(-F[i, 1])} for i in ordre
            ]
            self.label_optim.setText(f"Terminée -- {len(self.solutions_pareto)} solution(s) non-dominée(s).")

        self.table_pareto.setRowCount(len(self.solutions_pareto))
        for i, sol in enumerate(self.solutions_pareto):
            self.table_pareto.setItem(i, 0, QTableWidgetItem(f"{sol['dispersion']:.2f}"))
            self.table_pareto.setItem(i, 1, QTableWidgetItem(f"{sol['jeu']*1000:.1f}"))
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
        marches = marche_positions_verticales(
            self.geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
            amplitudes_deg=self._amplitudes_analyse(), pas_balayage_deg=self.sb_pas_balayage.value(),
        )
        self.ax_marche.clear()
        for nom, valeurs in marches.items():
            self.ax_marche.plot(self._amplitudes_analyse(), valeurs, "o-", label=nom)
        self.ax_marche.set_xlabel("Amplitude (deg)")
        self.ax_marche.set_ylabel("Marche diurne (s/j)")
        self.ax_marche.axhline(0, color="k", lw=0.8)
        self.ax_marche.legend()
        self.canvas_marche.draw_idle()
        self.tabs.setCurrentWidget(self.canvas_marche)

        self._dernier_dispersion = critere_dispersion_marche(marches)
        self._afficher_resultats_numeriques()

    def on_exporter_dxf(self):
        if self.geom is None:
            QMessageBox.information(self, "Export", "Générez ou chargez d'abord une géométrie.")
            return
        chemin, _ = QFileDialog.getSaveFileName(self, "Exporter en DXF", "", "DXF (*.dxf)")
        if chemin:
            chemin = _assurer_extension(chemin, ".dxf")
            dxf_io.exporter_dxf(self.geom, chemin)
            self._dernier_export = chemin
            QMessageBox.information(self, "Export", f"Exporté : {chemin}")

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
        _, _, _, jeu_avant = balayage_angulaire(self.geom, float(r_avant[0]), amplitude_max, pas)
        _, _, _, jeu_apres = balayage_angulaire(geom_reimporte, float(r_apres[0]), amplitude_max, pas)

        self.label_validation.setText(
            f"Écart r_min après réimport : {ecart_r_min:.2f} µm\n"
            f"Écart r_max après réimport : {ecart_r_max:.2f} µm\n"
            f"Jeu min avant export : {jeu_avant.min()*1000:.1f} µm\n"
            f"Jeu min après réimport : {jeu_apres.min()*1000:.1f} µm"
        )

    def _metriques_geometrie(self, geom, r_min: float = None):
        """(freq_hz, jeu_min_um, dispersion_s_j) d'une géométrie, en UN seul
        balayage elastica (marche + jeu partagent le même).

        r_min : rayon de collet imposé comme condition limite. Passé
        EXPLICITEMENT par l'appelant plutôt que déduit de la géométrie, pour
        que deux géométries comparées entre elles (maillage de travail vs
        maillage fin) utilisent forcément la même condition limite -- sans
        quoi on compare deux calculs qui ne diffèrent pas que par le maillage.
        Par défaut : le rayon effectif du premier point."""
        if r_min is None:
            r_min = float(geom.rayon()[0])
        marches, diag = marche_positions_verticales(
            geom, r_min, self.spec.freq_cible, self.spec.inertie_balancier,
            amplitudes_deg=self._amplitudes_analyse(), pas_balayage_deg=self.sb_pas_balayage.value(),
            return_diagnostics=True,
        )
        return (frequence_propre(geom, self.spec.inertie_balancier),
                float(diag["jeu_mm"].min()) * 1000.0, critere_dispersion_marche(marches))

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
            f0, j0, d0 = self._metriques_geometrie(self.geom, r_min=self.base.r_min)
            f1, j1, d1 = self._metriques_geometrie(construire_geometrie(base_fine), r_min=self.base.r_min)
        finally:
            QApplication.restoreOverrideCursor()

        # Incertitude de discrétisation : l'écart travail -> fin donne un ORDRE
        # DE GRANDEUR de l'erreur restante. S'il est comparable à la valeur
        # elle-même, le chiffre n'est pas résolu et ne doit pas être lu comme
        # un résultat physique (cas constaté : une dispersion annoncée à
        # 0.57 s/j alors que l'incertitude était de l'ordre de 1 s/j).
        incertitude = abs(d1 - d0)
        lignes = [
            f"Maillage de travail ({self.base.n_elements}) -> fin ({base_fine.n_elements}) :",
            f"Fréquence : {f0:.4f} -> {f1:.4f} Hz",
            f"Jeu minimal : {j0:.1f} -> {j1:.1f} µm",
            f"Dispersion : {d0:.2f} -> {d1:.2f} s/j",
        ]
        if incertitude > 0.5 * abs(d1):
            lignes.append(
                f"ATTENTION : l'écart de maillage ({incertitude:.2f} s/j) est du même "
                f"ordre que la dispersion elle-même -- ce chiffre n'est pas résolu. "
                f"Raffiner le maillage de travail avant de comparer des solutions."
            )
        self.label_verif_haute_res.setText("\n".join(lignes))

    def on_afficher_jeu_minimal(self):
        self._synchroniser_spec()
        if self.geom is None:
            QMessageBox.information(self, "Jeu minimal", "Générez ou chargez d'abord une géométrie.")
            return
        r_min = float(self.geom.rayon()[0])
        angles_deg, _, _, jeu_mm = balayage_angulaire(
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
            sauvegarder_configuration(_assurer_extension(chemin, ".json"), self.spec, self.base, self.bornes, self.hyperparams)

    def on_charger_configuration(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Charger une configuration complète", "", "JSON (*.json)")
        if not chemin:
            return
        try:
            spec, base, bornes, hyperparams = charger_configuration(chemin)
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
        self.sb_n_noeuds_ep.setValue(len(base.epaisseur_noeuds)); self.sb_n_noeuds_dr.setValue(len(base.delta_r_noeuds))

        self.sb_pop.setValue(hyperparams.pop_size); self.sb_gen.setValue(hyperparams.n_gen)
        self.sb_ga_seed.setValue(hyperparams.seed)

        geom, diag = construire_geometrie(self.base, return_diagnostics=True)
        self.geom = geom
        self.label_dxf.setText("(configuration complète chargée)")
        self._afficher_diagnostics(diag)
        self._tracer_geometrie()
