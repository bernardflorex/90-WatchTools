"""Fenêtre de l'outil d'échappement à ancre suisse.

À gauche, quatre volets :
  1. Tracé : fichier DXF, paramètres physiques, variables du tracé et bouton
     « Tracer l'échappement » ;
  2. Calculs : calcul statique et dynamique, résultats numériques ;
  3. Optimisation : variables, objectifs, contraintes, réglages NSGA-II ;
  4. Résultats : solutions de l'optimisation, à appliquer ou exporter.

À droite : l'échappement et son animation (100 positions, statique ou
dynamique), le couple sur le balancier, la cinématique, le bilan d'énergie,
le retard selon l'amplitude et le front de Pareto.
"""

import json
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Canvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as Barre
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from PySide6.QtCore import QThread, QTimer, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFileDialog, QFormLayout, QGridLayout, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QMainWindow, QMessageBox, QProgressBar,
                               QPushButton, QScrollArea, QSlider, QSpinBox, QSplitter,
                               QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from commun import navigation, theme
from commun.widgets import aide

from . import dynamique, statique
from .cinematique import Releve
from .dxf_import import lire_trace
from .geometrie import rotation
from .optimisation import (CONTRAINTES, OBJECTIFS, ConfigOpti, VariableOpti, optimiser,
                           variables_par_defaut)
from .parametres import Parametres
from .tables import DEG, Tables
from .trace_parametrique import LIBELLES, VariablesTrace, construire, mesurer

DXF_DEFAUT = Path(__file__).with_name("simulation.dxf")
N_IMAGES = 100

#: paramètres physiques : (attribut, libellé, unité affichée, facteur affichage -> SI,
#: min, max, pas, décimales, aide)
CHAMPS_PHYSIQUES = [
    ("couple_roue", "Couple sur la roue", "µN·m", 1e-6, 0.01, 20, 0.05, 3,
     "Couple moteur transmis par le rouage à la roue d'échappement."),
    ("mu", "Frottement µ", "", 1, 0, 0.6, 0.01, 3,
     "Coefficient de frottement de Coulomb aux contacts dent-palette et ellipse-fourchette."),
    ("e_chute", "Restitution dent-palette", "", 1, 0, 1, 0.05, 2,
     "Coefficient de restitution du choc de la dent sur la palette (chute). "
     "0 par défaut : aucun rebond n'est visible à la caméra rapide."),
    ("e_fourchette", "Restitution ellipse-fourchette", "", 1, 0, 1, 0.05, 2,
     "Coefficient de restitution du choc de l'ellipse dans la fourchette."),
    ("e_butee", "Restitution fourchette-butée", "", 1, 0, 1, 0.05, 2,
     "Coefficient de restitution du choc de la fourchette sur la butée."),
    ("inertie_roue", "Inertie de la roue", "g·mm²", 1e-9, 0, 10, 0.001, 5,
     "Inertie de la roue d'échappement et de son pignon autour de leur axe "
     "(roue acier 0,021 + pignon 0,0002 g·mm²)."),
    ("inertie_rouage", "Inertie du rouage ramenée", "g·mm²", 1e-9, 0, 10, 0.001, 5,
     "Inertie du rouage ramenée à l'axe de la roue d'échappement (0 si inconnue)."),
    ("inertie_ancre", "Inertie de l'ancre (dard compris)", "g·mm²", 1e-9, 0, 10, 0.0001, 5,
     "Inertie de l'ancre complète (palettes, fourchette, dard) autour de son axe."),
    ("inertie_balancier", "Inertie du balancier", "mg·cm²", 1e-10, 0.1, 100, 0.1, 3,
     "Même valeur que dans le programme spiral."),
    ("frequence", "Fréquence", "Hz", 1, 1, 10, 0.5, 3,
     "Fréquence propre du balancier-spiral (3,5 Hz = 25 200 alternances par heure). "
     "Avec l'inertie du balancier, elle fixe la raideur du spiral."),
    ("amplitude_deg", "Amplitude", "°", 1, 100, 360, 5, 1,
     "Amplitude du balancier (demi-oscillation) pour les calculs."),
    ("facteur_qualite", "Facteur de qualité", "", 1, 10, 2000, 10, 0,
     "Facteur de qualité du balancier seul, frottement de l'air compris "
     "(220 en position horizontale entre 310° et 250°)."),
]


def _spin(mini, maxi, val, pas, dec):
    sb = QDoubleSpinBox()
    sb.setRange(mini, maxi)
    sb.setDecimals(dec)
    sb.setSingleStep(pas)
    sb.setValue(val)
    sb.setKeyboardTracking(False)
    return sb


def _ligne_tableau(table, valeurs, gras=False):
    r = table.rowCount()
    table.insertRow(r)
    for c, v in enumerate(valeurs):
        it = QTableWidgetItem(str(v))
        if c > 0:
            it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if gras:
            f = it.font()
            f.setBold(True)
            it.setFont(f)
        table.setItem(r, c, it)


class Travail(QThread):
    """Exécute une fonction hors de l'interface."""
    fini = Signal(object)
    echec = Signal(str)
    avance = Signal(str)

    def __init__(self, fonction):
        super().__init__()
        self.fonction = fonction

    def run(self):
        try:
            self.fini.emit(self.fonction(self.avance.emit))
        except Exception as e:  # remonté à l'interface
            import traceback
            traceback.print_exc()
            self.echec.emit(f"{type(e).__name__} : {e}")


class FenetreEchappement(QMainWindow):

    def __init__(self, chemin_dxf=None):
        super().__init__()
        self.setWindowTitle("Échappement à ancre suisse")
        self.resize(1500, 920)
        self.ref = None              # tracé lu dans le DXF
        self.trace = None            # tracé courant (variables appliquées)
        self.rel = None              # relevé cinématique
        self.res_stat = self.eff_stat = None
        self.sim = self.res_dyn = None
        self.retards_dyn = None
        self.amp_equilibre = None
        self.solutions = []
        self.noms_opti = []
        self._travail = None
        self._arret_opti = False
        self.images_stat = self.images_dyn = None

        cadre = QSplitter(Qt.Orientation.Horizontal)
        cadre.addWidget(self._panneau())
        cadre.addWidget(self._vues())
        cadre.setStretchFactor(1, 1)
        cadre.setSizes([640, 860])
        self.setCentralWidget(cadre)
        self.statusBar().showMessage("Ouvrez un DXF ou tracez l'échappement de référence.")
        self._minuterie = QTimer(self)
        self._minuterie.setInterval(60)
        self._minuterie.timeout.connect(self._image_suivante)
        self.ouvrir_dxf(chemin_dxf or DXF_DEFAUT)

    # ======================================================================
    # Panneau de gauche
    # ======================================================================
    def _panneau(self):
        self.onglets = QTabWidget()
        self.onglets.addTab(self._page_trace(), "1. Tracé")
        self.onglets.addTab(self._page_calculs(), "2. Calculs")
        self.onglets.addTab(self._page_optimisation(), "3. Optimisation")
        self.onglets.addTab(self._page_resultats(), "4. Résultats")
        for i, txt in enumerate(("Tracé : DXF, paramètres physiques et variables du tracé.",
                                 "Calculs : statique, dynamique, retard, amplitude d'équilibre.",
                                 "Optimisation NSGA-II : variables, objectifs, contraintes.",
                                 "Résultats de l'optimisation : choisir, appliquer, exporter.")):
            self.onglets.setTabToolTip(i, txt)
        defile = QScrollArea()
        defile.setWidgetResizable(True)
        defile.setWidget(self.onglets)
        defile.setMinimumWidth(560)
        return defile

    # -- 1. Tracé --------------------------------------------------------------
    def _page_trace(self):
        page = QWidget()
        lay = QVBoxLayout(page)

        g = QGroupBox("Fichiers")
        gl = QGridLayout(g)
        b = QPushButton("Ouvrir un DXF…")
        aide(b, "Lit un tracé d'échappement (roue, ancre, ellipse, butées) dans un DXF. "
             "Il est orienté automatiquement en vue standard (balancier en haut, roue en bas, "
             "palette d'entrée à gauche) et devient le tracé de référence.")
        b.clicked.connect(lambda: self.ouvrir_dxf())
        self.lbl_dxf = QLabel("—")
        self.lbl_dxf.setWordWrap(True)
        gl.addWidget(b, 0, 0)
        gl.addWidget(self.lbl_dxf, 0, 1)
        aide(self.lbl_dxf, "Fichier DXF de référence, nombre de dents de la roue et "
             "transformation appliquée pour obtenir la vue standard.")
        b1 = QPushButton("Ouvrir un projet…")
        aide(b1, "Relit un projet .json : DXF de référence, paramètres physiques et variables "
             "du tracé, puis trace l'échappement.")
        b1.clicked.connect(self.ouvrir_projet)
        b2 = QPushButton("Enregistrer le projet…")
        aide(b2, "Enregistre dans un .json le chemin du DXF, les paramètres physiques et les "
             "variables du tracé.")
        b2.clicked.connect(self.enregistrer_projet)
        gl.addWidget(b1, 1, 0)
        gl.addWidget(b2, 1, 1)
        lay.addWidget(g)

        g = QGroupBox("Paramètres physiques")
        f = QFormLayout(g)
        self.champs = {}
        defaut = Parametres()
        for nom, lib, unite, fac, mi, ma, pas, dec, txt in CHAMPS_PHYSIQUES:
            sb = _spin(mi, ma, getattr(defaut, nom) / fac, pas, dec)
            if txt:
                aide(sb, txt)
            self.champs[nom] = (sb, fac)
            f.addRow(f"{lib}" + (f" ({unite})" if unite else ""), sb)
        lay.addWidget(g)

        g = QGroupBox("Variables du tracé")
        f = QFormLayout(g)
        self.vars = {}
        for nom, (lib, unite) in LIBELLES.items():
            if unite == "mm":
                sb = _spin(0.01, 20, 0, 0.005, 4)
            elif unite == "µm":
                sb = _spin(-200, 200, 0, 1, 1)
            else:
                sb = _spin(-30, 60, 0, 0.1, 3)
            self.vars[nom] = sb
            f.addRow(f"{lib} ({unite})", sb)
        aide(self.vars["entraxe_roue_ancre"], "Distance entre les axes de la roue et de l'ancre. "
             "La roue glisse le long de leur ligne des centres.")
        aide(self.vars["entraxe_ancre_balancier"], "Distance entre les axes de l'ancre et du "
             "balancier. Le balancier glisse le long de leur ligne des centres.")
        aide(self.vars["rayon_ellipse"], "Distance du centre de l'ellipse (cheville) à l'axe du "
             "balancier : règle le rapport de levée ancre/balancier.")
        aide(self.vars["largeur_fourchette"], "Largeur de l'entrée de fourchette : les deux flancs "
             "et les cornes s'écartent symétriquement. Règle l'ébat de l'ellipse.")
        aide(self.vars["course_ancre"], "Course totale de l'ancre entre les deux butées : les "
             "butées tournent symétriquement autour de l'axe de l'ancre.")
        for c in ("entree", "sortie"):
            aide(self.vars[f"position_{c}"], "Rotation de la palette autour de l'axe de l'ancre "
                 "(écart au dessin) : répartit repos et chute.")
            aide(self.vars[f"impulsion_{c}"], "Rotation du plan d'impulsion autour de l'arête "
                 "de repos (écart au dessin) : règle la levée.")
            aide(self.vars[f"longueur_{c}"], "Allongement du plan d'impulsion, arête de sortie "
                 "déplacée le long du plan (écart au dessin).")
            aide(self.vars[f"tirage_{c}"], "Rotation de la face de repos autour de l'arête de "
                 "repos (écart au dessin).")
        bz = QPushButton("Revenir au dessin")
        aide(bz, "Remet toutes les variables du tracé aux valeurs mesurées sur le DXF "
             "(annule les modifications et une solution d'optimisation appliquée). "
             "Ne recalcule rien : cliquez ensuite sur « Tracer l'échappement ».")
        bz.clicked.connect(self.variables_du_dessin)
        f.addRow(bz)
        lay.addWidget(g)

        self.btn_tracer = QPushButton("Tracer l'échappement")
        aide(self.btn_tracer, "Construit l'échappement avec les variables ci-dessus, calcule ses "
             "tables de contact, le relevé cinématique et le bilan statique, et met à jour "
             "l'animation et les courbes. À relancer après chaque modification des variables ou "
             "des paramètres physiques (quelques secondes).")
        self.btn_tracer.setStyleSheet("font-weight: bold; padding: 8px;")
        self.btn_tracer.clicked.connect(self.tracer)
        lay.addWidget(self.btn_tracer)
        lay.addStretch(1)
        return page

    def parametres(self):
        return Parametres(**{nom: sb.value() * fac for nom, (sb, fac) in self.champs.items()})

    def variables(self):
        return VariablesTrace(**{nom: sb.value() for nom, sb in self.vars.items()})

    def variables_du_dessin(self):
        if self.ref is None:
            return
        m = mesurer(self.ref)
        for nom, sb in self.vars.items():
            sb.setValue(m.get(nom, 0.0))

    # -- 2. Calculs ------------------------------------------------------------
    def _page_calculs(self):
        page = QWidget()
        lay = QVBoxLayout(page)
        g = QGroupBox("Calculs")
        gl = QGridLayout(g)
        b = QPushButton("Calcul dynamique (une oscillation)")
        b.clicked.connect(self.calcul_dynamique)
        aide(b, "Simule une oscillation complète du balancier à l'amplitude donnée, départ au "
             "repos sur la palette d'entrée : chocs, frottements, inertie de la roue et de "
             "l'ancre. Environ une seconde.")
        gl.addWidget(b, 0, 0, 1, 3)
        gl.addWidget(QLabel("Retard dynamique de"), 1, 0)
        self.amp_min = _spin(100, 360, 200, 10, 0)
        self.amp_max = _spin(100, 360, 320, 10, 0)
        aide(self.amp_min, "Plus petite amplitude (°) de la courbe de retard dynamique.")
        aide(self.amp_max, "Plus grande amplitude (°) de la courbe de retard dynamique.")
        gl.addWidget(self.amp_min, 1, 1)
        gl.addWidget(self.amp_max, 1, 2)
        self.n_amp = QSpinBox()
        self.n_amp.setRange(2, 20)
        self.n_amp.setValue(5)
        aide(self.n_amp, "Nombre d'amplitudes simulées entre les deux bornes (une oscillation "
             "chacune).")
        gl.addWidget(QLabel("en points"), 2, 0)
        gl.addWidget(self.n_amp, 2, 1)
        b = QPushButton("Retard selon l'amplitude")
        aide(b, "Marche due à l'échappement (s/j) en fonction de l'amplitude, une oscillation "
             "simulée par point ; courbe dans l'onglet « Retard », avec le calcul statique.")
        b.clicked.connect(self.retard_dynamique)
        gl.addWidget(b, 2, 2)
        b = QPushButton("Amplitude d'équilibre")
        aide(b, "Cherche l'amplitude stable du balancier : celle où l'énergie que l'échappement "
             "lui rend sur une oscillation compense exactement ses pertes (facteur de qualité). "
             "Au-dessus l'amplitude baisse d'une oscillation à l'autre, au-dessous elle monte. "
             "Quatre à six oscillations simulées, quelques secondes.")
        b.clicked.connect(self.amplitude_equilibre)
        gl.addWidget(b, 3, 0, 1, 3)
        self.lbl_calcul = QLabel("")
        self.lbl_calcul.setWordWrap(True)
        aide(self.lbl_calcul, "Résumé du dernier calcul. La levée est l'angle parcouru par le "
             "balancier entre l'entrée de l'ellipse dans la fourchette et la fin de "
             "l'impulsion, pour chaque palette.")
        gl.addWidget(self.lbl_calcul, 4, 0, 1, 3)
        lay.addWidget(g)

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Grandeur", "Entrée", "Sortie"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setMinimumHeight(560)
        aide(self.table, "Colonnes : alternance où travaille la palette d'entrée / de sortie. "
             "Chute, repos et chemin perdu sont mesurés sur la palette qui reçoit la dent.")
        lay.addWidget(self.table, 1)
        return page

    # -- 3. Optimisation -------------------------------------------------------
    def _page_optimisation(self):
        page = QWidget()
        lay = QVBoxLayout(page)
        g = QGroupBox("Variables (bornes)")
        self.grille_vars = QGridLayout(g)
        self.opti_vars = {}
        lay.addWidget(g)

        g = QGroupBox("Objectifs")
        v = QVBoxLayout(g)
        self.opti_obj = {}
        for k, (lib, sens) in OBJECTIFS.items():
            cb = QCheckBox(lib + (" — maximiser" if sens > 0 else " — minimiser"))
            cb.setChecked(k in ("rendement", "retard"))
            aide(cb, {"rendement": "Énergie reçue par le balancier / travail du couple moteur "
                                   "sur un pas de roue (statique ou dynamique selon l'évaluation).",
                      "retard": "Valeur absolue de la marche due à l'échappement, à l'amplitude "
                                "visée (paramètres physiques).",
                      "amplitude": "Amplitude stable du balancier (énergie rendue = pertes du "
                                   "facteur de qualité). En dynamique, elle coûte quatre à six "
                                   "oscillations par candidat."}[k])
            self.opti_obj[k] = cb
            v.addWidget(cb)
        lay.addWidget(g)

        g = QGroupBox("Contraintes (sur les deux alternances)")
        gl = QGridLayout(g)
        gl.addWidget(QLabel("min"), 0, 1)
        gl.addWidget(QLabel("max"), 0, 2)
        defauts = {"chute": (1.0, 3.0), "chemin_perdu": (10.0, 40.0),
                   "repos_chute": (20.0, 80.0), "levee": (40.0, 56.0)}
        self.opti_contr = {}
        for i, (k, (lib, _)) in enumerate(CONTRAINTES.items(), start=1):
            cb = QCheckBox(lib)
            cb.setChecked(k in ("chute", "chemin_perdu"))
            aide(cb, {"chute": "Rotation libre de la roue entre l'échappement d'une dent et la "
                               "chute de la suivante sur l'autre palette.",
                      "chemin_perdu": "Course de la dent sur la face de repos entre le premier "
                                      "contact de chute et le repos, ancre sur la butée.",
                      "repos_chute": "Recouvrement de la dent sur la face de repos au moment de "
                                     "la chute (sécurité contre le renversement).",
                      "levee": "Angle du balancier entre l'entrée de l'ellipse dans la "
                               "fourchette et la fin de l'impulsion."}[k]
                 + " Contrainte respectée sur les deux alternances quand la case est cochée.")
            lo, hi = _spin(-1000, 1000, defauts[k][0], 1, 2), _spin(-1000, 1000, defauts[k][1], 1, 2)
            aide(lo, "Valeur minimale admise.")
            aide(hi, "Valeur maximale admise.")
            gl.addWidget(cb, i, 0)
            gl.addWidget(lo, i, 1)
            gl.addWidget(hi, i, 2)
            self.opti_contr[k] = (cb, lo, hi)
        lay.addWidget(g)

        g = QGroupBox("NSGA-II")
        f = QFormLayout(g)
        self.pop = QSpinBox()
        self.pop.setRange(4, 400)
        self.pop.setValue(24)
        self.gen = QSpinBox()
        self.gen.setRange(1, 500)
        self.gen.setValue(12)
        self.proc = QSpinBox()
        self.proc.setRange(1, 128)
        import os
        self.proc.setValue(max(1, (os.cpu_count() or 2) - 1))
        self.graine = QSpinBox()
        self.graine.setRange(0, 10 ** 6)
        self.graine.setValue(1)
        self.mode = QComboBox()
        self.mode.addItems(["statique", "dynamique"])
        aide(self.mode, "Statique : bilan quasi statique, environ deux secondes par candidat. "
             "Dynamique : une oscillation simulée par candidat (rendement et marche "
             "dynamiques), environ une seconde de plus ; quatre à six si l'amplitude "
             "d'équilibre est un objectif.")
        aide(self.pop, "Nombre de candidats par génération de NSGA-II.")
        aide(self.gen, "Nombre de générations ; l'arrêt reste possible à tout moment.")
        aide(self.proc, "Processus de calcul en parallèle (un candidat chacun).")
        aide(self.graine, "Graine du tirage aléatoire : la même graine refait la même "
             "optimisation.")
        f.addRow("Population", self.pop)
        f.addRow("Générations", self.gen)
        f.addRow("Processus", self.proc)
        f.addRow("Graine", self.graine)
        f.addRow("Évaluation", self.mode)
        lay.addWidget(g)

        h = QHBoxLayout()
        self.btn_opti = QPushButton("Lancer l'optimisation")
        aide(self.btn_opti, "Lance NSGA-II sur les variables cochées. Les variables non cochées "
             "gardent leur valeur actuelle (volet Tracé). Le front de Pareto s'affiche dans le "
             "volet Résultats et l'onglet Pareto.")
        self.btn_opti.clicked.connect(self.lancer_optimisation)
        self.btn_stop = QPushButton("Arrêter")
        aide(self.btn_stop, "Arrête l'optimisation à la fin de la génération en cours ; le "
             "front obtenu jusque-là est conservé.")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.arreter_optimisation)
        h.addWidget(self.btn_opti)
        h.addWidget(self.btn_stop)
        lay.addLayout(h)
        self.barre_opti = QProgressBar()
        aide(self.barre_opti, "Générations effectuées.")
        lay.addWidget(self.barre_opti)
        self.lbl_opti = QLabel("")
        self.lbl_opti.setWordWrap(True)
        lay.addWidget(self.lbl_opti)
        lay.addStretch(1)
        return page

    def _remplir_variables_opti(self):
        gl = self.grille_vars
        while gl.count():
            w = gl.takeAt(0).widget()
            if w:
                w.deleteLater()
        gl.addWidget(QLabel("min"), 0, 1)
        gl.addWidget(QLabel("max"), 0, 2)
        self.opti_vars = {}
        for i, v in enumerate(variables_par_defaut(self.ref), start=1):
            lib, unite = LIBELLES[v.nom]
            cb = QCheckBox(f"{lib} ({unite})")
            cb.setChecked(v.nom in ("impulsion_entree", "impulsion_sortie"))
            aide(cb, "Cochée : variable libre de l'optimisation, entre les bornes min et max. "
                 "Décochée : gardée à sa valeur du volet Tracé.")
            dec = 4 if unite == "mm" else 2
            lo, hi = _spin(-1000, 1000, v.mini, 0.01, dec), _spin(-1000, 1000, v.maxi, 0.01, dec)
            aide(lo, "Borne inférieure de la variable.")
            aide(hi, "Borne supérieure de la variable.")
            gl.addWidget(cb, i, 0)
            gl.addWidget(lo, i, 1)
            gl.addWidget(hi, i, 2)
            self.opti_vars[v.nom] = (cb, lo, hi)

    # -- 4. Résultats ------------------------------------------------------------
    def _page_resultats(self):
        page = QWidget()
        lay = QVBoxLayout(page)
        self.table_sol = QTableWidget(0, 0)
        self.table_sol.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_sol.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_sol.verticalHeader().setVisible(False)
        self.table_sol.setMinimumHeight(420)
        self.table_sol.itemSelectionChanged.connect(self._solution_choisie)
        aide(self.table_sol, "Solutions du front de Pareto : objectifs, contraintes (entrée / "
             "sortie) et variables. Une ligne sélectionnée est mise en évidence dans l'onglet "
             "Pareto.")
        lay.addWidget(self.table_sol, 1)
        h = QHBoxLayout()
        b = QPushButton("Appliquer au tracé")
        aide(b, "Reporte les variables de la solution sélectionnée dans le volet Tracé et trace "
             "l'échappement correspondant.")
        b.clicked.connect(self.appliquer_solution)
        h.addWidget(b)
        b = QPushButton("Exporter les solutions (CSV)…")
        aide(b, "Écrit toutes les solutions du front (objectifs, contraintes, variables) dans "
             "un fichier CSV, séparateur point-virgule.")
        b.clicked.connect(self.exporter_solutions)
        h.addWidget(b)
        lay.addLayout(h)
        b = QPushButton("Exporter le tracé courant en DXF…")
        aide(b, "Écrit le tracé courant (variables appliquées) dans un DXF, un calque par pièce.")
        b.clicked.connect(self.exporter_dxf)
        lay.addWidget(b)
        self.lbl_sol = QLabel("Aucune optimisation lancée.")
        self.lbl_sol.setWordWrap(True)
        lay.addWidget(self.lbl_sol)
        return page

    # ======================================================================
    # Vues de droite
    # ======================================================================
    def _canvas(self, figsize=(7, 6)):
        # « constrained » : la mise en page suit la taille du canevas, la
        # figure occupe toute la place disponible
        fig = Figure(figsize=figsize, layout="constrained")
        canvas = Canvas(fig)
        navigation.activer(canvas)
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(Barre(canvas, w))
        v.addWidget(canvas, 1)
        return fig, canvas, w, v

    def _vues(self):
        self.vues = QTabWidget()
        # échappement et animation
        self.fig_ech, self.can_ech, w, v = self._canvas((7, 8))
        self.ax_ech = self.fig_ech.add_subplot(111)
        h = QHBoxLayout()
        self.chk_dyn = QCheckBox("Dynamique")
        aide(self.chk_dyn, "Animation et courbes issues du calcul dynamique plutôt que du "
             "relevé statique.")
        self.chk_dyn.toggled.connect(self._changer_source)
        self.btn_lire = QPushButton("▶")
        aide(self.btn_lire, "Lance ou arrête l'animation (100 positions, en boucle).")
        self.btn_lire.setCheckable(True)
        self.btn_lire.toggled.connect(self._lire)
        self.curseur = QSlider(Qt.Orientation.Horizontal)
        self.curseur.setRange(0, N_IMAGES - 1)
        self.curseur.valueChanged.connect(self._afficher_image)
        aide(self.curseur, "Position dans le cycle : 50 images par alternance, palette "
             "d'entrée d'abord.")
        self.lbl_image = QLabel("")
        self.lbl_image.setMinimumWidth(230)
        aide(self.lbl_image, "Statique : alternance, phase de la fonction et angle d'ancre. "
             "Dynamique : instant depuis le départ et contacts fermés.")
        h.addWidget(self.chk_dyn)
        h.addWidget(self.btn_lire)
        h.addWidget(self.curseur, 1)
        h.addWidget(self.lbl_image)
        v.addLayout(h)
        self.vues.addTab(w, "Échappement")
        # couple au balancier
        self.fig_couple, self.can_couple, w, _ = self._canvas()
        self.vues.addTab(w, "Couple au balancier")
        # cinématique
        self.fig_cin, self.can_cin, w, _ = self._canvas()
        self.vues.addTab(w, "Cinématique")
        # énergie
        self.fig_en, self.can_en, w, _ = self._canvas()
        self.vues.addTab(w, "Énergie")
        # retard
        self.fig_ret, self.can_ret, w, _ = self._canvas()
        self.vues.addTab(w, "Retard")
        # Pareto
        self.fig_par, self.can_par, w, _ = self._canvas()
        self.can_par.mpl_connect("pick_event", self._clic_pareto)
        self.vues.addTab(w, "Pareto")
        for i, txt in enumerate((
                "L'échappement et son animation, en statique (relevé) ou en dynamique.",
                "Couple de la fourchette sur l'axe du balancier : positif quand l'ancre est "
                "motrice, négatif quand le balancier fournit le dégagement, nul hors contact ; "
                "chocs reçus par le balancier.",
                "Positions de la roue et de l'ancre au cours du relevé.",
                "Bilan d'énergie : travail moteur, énergie reçue par le balancier et pertes.",
                "Marche due à l'échappement selon l'amplitude, statique et dynamique.",
                "Front de Pareto de la dernière optimisation (cliquez sur un point).")):
            self.vues.setTabToolTip(i, txt)
        return self.vues

    # ======================================================================
    # Fichiers
    # ======================================================================
    def ouvrir_dxf(self, chemin=None):
        if chemin is None:
            chemin, _ = QFileDialog.getOpenFileName(self, "Ouvrir un tracé", str(DXF_DEFAUT.parent),
                                                    "DXF (*.dxf)")
            if not chemin:
                return
        try:
            self.ref = lire_trace(str(chemin))
        except Exception as e:
            QMessageBox.critical(self, "DXF", f"Lecture impossible : {e}")
            return
        self.chemin_dxf = str(chemin)
        self.lbl_dxf.setText(f"{Path(chemin).name}  ({self.ref.nombre_dents} dents, "
                             f"vue standard : {self.ref.transformation})")
        self.variables_du_dessin()
        self._remplir_variables_opti()
        self.trace = self.ref
        self._dessiner_trace_seul()
        self.tracer()

    def enregistrer_projet(self):
        chemin, _ = QFileDialog.getSaveFileName(self, "Enregistrer le projet", "", "JSON (*.json)")
        if not chemin:
            return
        d = {"dxf": self.chemin_dxf, "parametres": self.parametres().en_dict(),
             "variables": self.variables().en_dict()}
        Path(chemin).write_text(json.dumps(d, indent=2, ensure_ascii=False), encoding="utf-8")
        self.statusBar().showMessage(f"Projet enregistré : {chemin}")

    def ouvrir_projet(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Ouvrir un projet", "", "JSON (*.json)")
        if not chemin:
            return
        d = json.loads(Path(chemin).read_text(encoding="utf-8"))
        if d.get("dxf") and Path(d["dxf"]).exists() and d["dxf"] != getattr(self, "chemin_dxf", None):
            self.ref = lire_trace(d["dxf"])
            self.chemin_dxf = d["dxf"]
            self._remplir_variables_opti()
        for nom, val in d.get("parametres", {}).items():
            if nom in self.champs:
                sb, fac = self.champs[nom]
                sb.setValue(val / fac)
        for nom, val in d.get("variables", {}).items():
            if nom in self.vars and val is not None and np.isfinite(val):
                self.vars[nom].setValue(val)
        self.tracer()

    def exporter_dxf(self):
        if self.trace is None:
            return
        chemin, _ = QFileDialog.getSaveFileName(self, "Exporter le tracé", "trace.dxf", "DXF (*.dxf)")
        if not chemin:
            return
        import ezdxf
        doc = ezdxf.new()
        msp = doc.modelspace()
        for corps, calque in ((self.trace.roue, "ROUE"), (self.trace.ancre, "ANCRE"),
                              (self.trace.balancier, "BALANCIER")):
            doc.layers.add(calque)
            for e in corps.elements:
                msp.add_lwpolyline([tuple(p) for p in e.points], dxfattribs={"layer": calque})
        doc.layers.add("BUTEES")
        for b in self.trace.butees:
            msp.add_circle(tuple(b.axe), b.elements[0].rayon, dxfattribs={"layer": "BUTEES"})
        for corps in (self.trace.roue, self.trace.ancre, self.trace.balancier):
            msp.add_circle(tuple(corps.axe), 0.05, dxfattribs={"layer": "AXES"})
        doc.saveas(chemin)
        self.statusBar().showMessage(f"Tracé exporté : {chemin}")

    # ======================================================================
    # Travaux de fond
    # ======================================================================
    def _lancer(self, fonction, fini, message):
        if self._travail is not None and self._travail.isRunning():
            QMessageBox.information(self, "Calcul en cours", "Un calcul est déjà en cours.")
            return False
        self._travail = Travail(fonction)
        self._travail.fini.connect(fini)
        self._travail.echec.connect(self._echec)
        self._travail.avance.connect(self.statusBar().showMessage)
        self.statusBar().showMessage(message)
        self._travail.start()
        return True

    def _echec(self, msg):
        self.statusBar().showMessage("Échec du calcul")
        self.btn_opti.setEnabled(True)
        self.btn_stop.setEnabled(False)
        QMessageBox.critical(self, "Calcul", msg)

    # -- tracé + statique ------------------------------------------------------
    def tracer(self):
        if self.ref is None:
            return
        ref, var, par = self.ref, self.variables(), self.parametres()

        def calcul(avance):
            avance("Construction du tracé…")
            trace = construire(ref, var)
            avance("Tables de contact…")
            T = Tables(trace)
            avance("Relevé cinématique…")
            rel = Releve(T)
            avance("Efforts statiques…")
            res, eff = statique.bilan(rel, par)
            return trace, rel, res, eff
        self.btn_tracer.setEnabled(False)
        self._lancer(calcul, self._trace_fini, "Tracé de l'échappement…")

    def _trace_fini(self, r):
        self.btn_tracer.setEnabled(True)
        self.trace, self.rel, self.res_stat, self.eff_stat = r
        self.res_dyn = self.sim = None
        self.amp_equilibre = None
        self.retards_dyn = None
        self.images_dyn = None
        self.chk_dyn.setChecked(False)
        self._preparer_images_statiques()
        self._dessiner_trace_seul()
        self._afficher_image(self.curseur.value())
        self._remplir_tableau()
        self._tracer_courbes()
        self.lbl_calcul.setText(
            f"Statique : rendement {100 * self.res_stat['rendement']:.1f} %, retard "
            f"{self.res_stat['retard_s_j']:+.2f} s/j. {self._texte_levees()}")
        self.statusBar().showMessage(
            f"Tracé fait : rendement statique {100 * self.res_stat['rendement']:.1f} %, "
            f"retard {self.res_stat['retard_s_j']:+.2f} s/j à {self.parametres().amplitude_deg:.0f}°")

    # -- dynamique ----------------------------------------------------------------
    def calcul_dynamique(self):
        if self.rel is None:
            return
        rel, par = self.rel, self.parametres()

        def calcul(avance):
            sim = dynamique.Simulateur(rel, par)
            return sim, sim.simuler()
        self._lancer(calcul, self._dynamique_fini, "Simulation dynamique d'une oscillation…")

    def _dynamique_fini(self, r):
        self.sim, self.res_dyn = r
        b = self.res_dyn.bilan
        pas = self.trace.pas / DEG
        avance = np.degrees(b["travail_moteur"] / self.parametres().couple_roue)
        ok = abs(avance - pas) < 0.05 * pas and self.res_dyn.valide
        self.lbl_calcul.setText(
            f"Dynamique : rendement {100 * b['rendement']:.1f} %, marche {b['marche_s_j']:+.2f} s/j, "
            f"amplitude {b['amplitude_depart_deg']:.0f}° → {b['amplitude_fin_deg']:.2f}° "
            f"en une oscillation. {self._texte_levees()}"
            + ("" if ok else f"ATTENTION : la roue a tourné de {avance:.1f}° au lieu de {pas:.0f}° "
               "(échappement défaillant)."))
        self._preparer_images_dynamiques()
        self._remplir_tableau()
        self._tracer_courbes()
        self.chk_dyn.setChecked(True)
        self.statusBar().showMessage("Simulation dynamique terminée")

    def _texte_levees(self):
        lev = {a.palette: a.mesures["levee_deg"] for a in self.rel.alternances}
        return (f"Levée {lev.get('entrée', float('nan')):.2f}° (entrée) / "
                f"{lev.get('sortie', float('nan')):.2f}° (sortie).")

    def amplitude_equilibre(self):
        if self.rel is None:
            return
        rel, par = self.rel, self.parametres()

        def calcul(avance):
            def rappel(A, d):
                avance(f"Amplitude d'équilibre : {A:.2f}° → {d:+.3f}° par oscillation…")
            return dynamique.amplitude_equilibre(rel, par, rappel=rappel)
        self._lancer(calcul, self._equilibre_fini, "Recherche de l'amplitude d'équilibre…")

    def _equilibre_fini(self, r):
        A, res, pts = r
        self.amp_equilibre = (A, res.bilan)
        b = res.bilan
        self.lbl_calcul.setText(
            f"Amplitude d'équilibre : {A:.1f}° (visée {self.parametres().amplitude_deg:.0f}°). "
            f"Rendement {100 * b['rendement']:.1f} %, marche {b['marche_s_j']:+.2f} s/j à cette "
            f"amplitude. Une oscillation : " +
            ", ".join(f"{a:.0f}° {d:+.2f}°" for a, d in pts[:3]) + ". " + self._texte_levees())
        self._remplir_tableau()
        self.statusBar().showMessage(f"Amplitude d'équilibre : {A:.1f}°")

    def retard_dynamique(self):
        if self.rel is None:
            return
        rel, par = self.rel, self.parametres()
        amps = np.linspace(self.amp_min.value(), self.amp_max.value(), self.n_amp.value())

        def calcul(avance):
            out = []
            for k, a in enumerate(amps):
                avance(f"Retard dynamique : {a:.0f}° ({k + 1}/{len(amps)})…")
                r = dynamique.Simulateur(rel, par).simuler(amplitude_deg=a)
                out.append((a, r.bilan["marche_s_j"], r.bilan["rendement"]))
            return np.array(out)
        self._lancer(calcul, self._retard_fini, "Retard dynamique selon l'amplitude…")

    def _retard_fini(self, r):
        self.retards_dyn = r
        self._tracer_courbes()
        self.vues.setCurrentIndex(4)
        self.statusBar().showMessage("Retard dynamique calculé")

    # -- optimisation ---------------------------------------------------------------
    def config_opti(self):
        vs = [VariableOpti(n, cb.isChecked(), lo.value(), hi.value())
              for n, (cb, lo, hi) in self.opti_vars.items()]
        obj = [k for k, cb in self.opti_obj.items() if cb.isChecked()]
        contr = {k: (lo.value(), hi.value()) for k, (cb, lo, hi) in self.opti_contr.items()
                 if cb.isChecked()}
        return ConfigOpti(vs, obj, contr, self.pop.value(), self.gen.value(), self.proc.value(),
                          self.graine.value(), self.mode.currentText())

    def lancer_optimisation(self):
        if self.ref is None:
            return
        cfg = self.config_opti()
        if not any(v.actif for v in cfg.variables) or not cfg.objectifs:
            QMessageBox.information(self, "Optimisation",
                                    "Choisissez au moins une variable et un objectif.")
            return
        ref, par = self.ref, self.parametres()
        # Les variables non optimisées restent à leur valeur courante : elles
        # sont appliquées une fois pour toutes au tracé de départ.
        courant = self.variables()
        fixes = {v.nom: getattr(courant, v.nom) for v in cfg.variables if not v.actif}
        cfg.variables = [v for v in cfg.variables if v.actif]
        ref_fixe = construire(ref, VariablesTrace(**fixes))
        self._arret_opti = False
        self.barre_opti.setRange(0, cfg.generations)
        self.barre_opti.setValue(0)

        def calcul(avance):
            def rappel(g, n, front):
                avance(f"Génération {g}/{cfg.generations}, {n} évaluations, "
                       f"{len(front)} solutions au front")
            return optimiser(ref_fixe, par, cfg, rappel=rappel, arreter=lambda: self._arret_opti)
        self._fixes_opti = fixes
        if self._lancer(calcul, self._opti_finie, "Optimisation NSGA-II…"):
            self.btn_opti.setEnabled(False)
            self.btn_stop.setEnabled(True)
            self._travail.avance.connect(self._avance_opti)

    def _avance_opti(self, msg):
        self.lbl_opti.setText(msg)
        if msg.startswith("Génération"):
            try:
                self.barre_opti.setValue(int(msg.split()[1].split("/")[0]))
            except ValueError:
                pass

    def arreter_optimisation(self):
        self._arret_opti = True
        self.lbl_opti.setText("Arrêt demandé : fin de la génération en cours…")

    def _opti_finie(self, r):
        self.btn_opti.setEnabled(True)
        self.btn_stop.setEnabled(False)
        sol, noms = r
        for s in sol:
            s["variables"].update(self._fixes_opti)
        self.solutions, self.noms_opti = sol, noms
        self._remplir_solutions()
        self._tracer_pareto()
        self.onglets.setCurrentIndex(3)
        self.vues.setCurrentIndex(5)
        self.statusBar().showMessage(f"Optimisation terminée : {len(sol)} solutions au front")

    # ======================================================================
    # Résultats
    # ======================================================================
    def _remplir_tableau(self):
        t = self.table
        t.setRowCount(0)
        if self.rel is None:
            return
        _ligne_tableau(t, ["Cinématique", "", ""], gras=True)
        for nom, e, s in self.rel.tableau():
            _ligne_tableau(t, [nom, e, s])
        _ligne_tableau(t, ["Statique", "", ""], gras=True)
        for nom, e, s in statique.tableau(self.res_stat, self.parametres()):
            _ligne_tableau(t, [nom, e, s])
        if self.res_dyn is not None:
            b = self.res_dyn.bilan
            _ligne_tableau(t, ["Dynamique (une oscillation)", "", ""], gras=True)
            for nom, val in (("Rendement dynamique (%)", f"{100 * b['rendement']:.1f}"),
                             ("Travail moteur (nJ)", f"{1e9 * b['travail_moteur']:.1f}"),
                             ("Énergie reçue par le balancier (nJ)",
                              f"{1e9 * b['energie_recue_balancier']:.1f}"),
                             ("Pertes par frottement (nJ)", f"{1e9 * b['frottements']:.1f}"),
                             ("Pertes aux chocs (nJ)", f"{1e9 * b['chocs']:.1f}"),
                             ("Pertes du balancier, Q (nJ)", f"{1e9 * b['amortissement']:.1f}"),
                             ("Marche (s/j)", f"{b['marche_s_j']:+.2f}"),
                             ("Amplitude après une oscillation (°)", f"{b['amplitude_fin_deg']:.2f}"),
                             ("Amplitude de départ (°)", f"{b['amplitude_depart_deg']:.2f}"),
                             ("Nombre de chocs", f"{len(dynamique.chocs_significatifs(self.res_dyn))}"),
                             ("Écart du bilan d'énergie (nJ)", f"{1e9 * b['ecart_bilan']:.2e}")):
                _ligne_tableau(t, [nom, val, ""])
            # recul de la roue : la première alternance simulée est celle où
            # le balancier tourne en sens négatif
            reculs = b.get("recul_roue_deg", [])
            ordre = sorted(self.rel.alternances, key=lambda a: a.palette != "entrée")
            par_palette = {a.palette: r for a, r in zip(ordre, reculs)}
            if par_palette:
                _ligne_tableau(t, ["Recul de la roue au dégagement (°)",
                                   f"{par_palette.get('entrée', float('nan')):.3f}",
                                   f"{par_palette.get('sortie', float('nan')):.3f}"])
            # décollement de la dent après le dégagement (comparable à la caméra rapide)
            dec = {d["palette"]: d for d in dynamique.decollement(self.sim, self.res_dyn)}
            for nom, cle, f in (("Décollement maximal de la dent (µm)", "jeu_max_um", "{:.1f}"),
                                ("Durée du décollement (ms)", "duree_ms", "{:.3f}"),
                                ("Reprise de contact, % du plan d'impulsion", "reprise_pct", "{:.0f}")):
                _ligne_tableau(t, [nom] + [f.format(dec[c][cle]) if c in dec else "—"
                                           for c in ("entrée", "sortie")])
        if self.amp_equilibre is not None:
            A, b = self.amp_equilibre
            _ligne_tableau(t, ["Équilibre", "", ""], gras=True)
            for nom, val in (("Amplitude d'équilibre (°)", f"{A:.2f}"),
                             ("Rendement à l'équilibre (%)", f"{100 * b['rendement']:.1f}"),
                             ("Marche à l'équilibre (s/j)", f"{b['marche_s_j']:+.2f}")):
                _ligne_tableau(t, [nom, val, ""])
        t.resizeColumnsToContents()
        t.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)

    def _remplir_solutions(self):
        t = self.table_sol
        sol = self.solutions
        cols = ["#", "Rendement", "Retard (s/j)", "Amplitude d'équilibre (°)"] + \
            [f"{CONTRAINTES[k][0]} E/S" for k in CONTRAINTES] + \
            [LIBELLES[n][0] for n in self.noms_opti]
        t.clear()
        t.setColumnCount(len(cols))
        t.setRowCount(0)
        t.setHorizontalHeaderLabels(cols)
        for i, s in enumerate(sol):
            m = s["mesures"]
            ligne = [str(i + 1), f"{100 * m['rendement']:.1f} %", f"{m['retard']:+.2f}",
                     f"{m['amplitude']:.0f}"]
            ligne += [" / ".join(f"{v:.2f}" for v in m[k]) for k in CONTRAINTES]
            ligne += [f"{s['variables'][n]:.4g}" for n in self.noms_opti]
            _ligne_tableau(t, ligne)
        t.resizeColumnsToContents()
        self.lbl_sol.setText(f"{len(sol)} solutions au front de Pareto"
                             + ("" if all(s["faisable"] for s in sol) else
                                " (certaines violent encore une contrainte)") +
                             ". Sélectionnez une ligne puis « Appliquer au tracé ».")

    def _solution_choisie(self):
        rows = {i.row() for i in self.table_sol.selectedItems()}
        self._sol_courante = min(rows) if rows else None
        self._tracer_pareto()

    def appliquer_solution(self):
        k = getattr(self, "_sol_courante", None)
        if k is None or k >= len(self.solutions):
            return
        for nom, val in self.solutions[k]["variables"].items():
            self.vars[nom].setValue(val)
        self.onglets.setCurrentIndex(0)
        self.tracer()

    def exporter_solutions(self):
        if not self.solutions:
            return
        chemin, _ = QFileDialog.getSaveFileName(self, "Exporter les solutions", "solutions.csv",
                                                "CSV (*.csv)")
        if not chemin:
            return
        noms = list(self.solutions[0]["variables"])
        with open(chemin, "w", encoding="utf-8") as f:
            f.write(";".join(["rendement", "retard_s_j", "amplitude_deg"] +
                             [f"{k}_{c}" for k in CONTRAINTES for c in ("entree", "sortie")] +
                             noms) + "\n")
            for s in self.solutions:
                m = s["mesures"]
                f.write(";".join([f"{m['rendement']:.5f}", f"{m['retard']:.4f}",
                                  f"{m['amplitude']:.2f}"] +
                                 [f"{v:.4f}" for k in CONTRAINTES for v in m[k]] +
                                 [f"{s['variables'][n]:.6g}" for n in noms]) + "\n")
        self.statusBar().showMessage(f"Solutions exportées : {chemin}")

    # ======================================================================
    # Dessin de l'échappement et animation
    # ======================================================================
    def _dessiner_trace_seul(self):
        ax = self.ax_ech
        ax.clear()
        t = self.trace
        if t is None:
            return
        self._lc = {}
        for nom, corps, couleur, lw in (("roue", t.roue, "0.75", 0.8),
                                        ("ancre", t.ancre, theme.ACC.name(), 1.2),
                                        ("balancier", t.balancier, theme.ACC4.name(), 1.2)):
            lc = LineCollection(corps.segments, colors=couleur, linewidths=lw)
            ax.add_collection(lc)
            self._lc[nom] = lc
        for b in t.butees:
            ax.add_collection(LineCollection(b.segments, colors=theme.ACC2.name(), linewidths=1.0))
        for corps in (t.roue, t.ancre, t.balancier):
            ax.plot(*corps.axe, "+", color=theme.encre(), ms=8)
        # ligne des centres ancre-balancier
        ax.plot(*np.array([t.ancre.axe, t.balancier.axe]).T, ":", color=theme.gris(0.5), lw=0.8)
        tout = np.concatenate([t.roue.segments.reshape(-1, 2), t.balancier.segments.reshape(-1, 2),
                               t.ancre.segments.reshape(-1, 2)])
        lo, hi = tout.min(0), tout.max(0)
        m = 0.05 * (hi - lo).max()
        ax.set_xlim(lo[0] - m, hi[0] + m)
        ax.set_ylim(lo[1] - m, hi[1] + m)
        # échelle égale obtenue en étendant les données, pas en rétrécissant
        # le cadre : le dessin prend toute la largeur de la vue
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(alpha=0.2)
        self.can_ech.draw_idle()

    def _poser(self, phi, alpha, beta):
        t = self.trace
        self._lc["roue"].set_segments(rotation(t.roue.segments, t.roue.axe, phi))
        self._lc["ancre"].set_segments(rotation(t.ancre.segments, t.ancre.axe, alpha))
        self._lc["balancier"].set_segments(rotation(t.balancier.segments, t.balancier.axe, beta))
        self.can_ech.draw_idle()

    def _preparer_images_statiques(self):
        """100 positions du relevé : 50 par alternance, à pas d'ancre réguliers.
        Le balancier suit la fourchette (ellipse poussante au dégagement, menée
        à l'impulsion), et file librement avant et après."""
        rel, T = self.rel, self.rel.T
        images = []
        # le cycle commence sur la palette d'entrée, comme le calcul dynamique
        for alt in sorted(rel.alternances, key=lambda a: a.palette != "entrée"):
            ev = alt.evenements
            pousse, mene = alt.pousse.beta, alt.mene.beta
            s_a = np.sign(alt.arrivee - alt.depart)
            b_deg = T.interp(pousse, ev["fin_degagement"])
            b_fin = T.interp(mene, ev["chute"])
            pente = (b_fin - b_deg) / (ev["chute"] - ev["reprise_ellipse"])
            for k in np.linspace(0, len(alt.alpha) - 1, N_IMAGES // 2).astype(int):
                a, phi = alt.alpha[k], alt.phi[k]
                if (a - ev["fin_degagement"]) * s_a <= 0:
                    beta, phase = T.interp(pousse, a), "dégagement"
                elif (a - ev["reprise_ellipse"]) * s_a <= 0:
                    beta, phase = b_deg, "ébat de l'ellipse"
                elif (a - ev["chute"]) * s_a <= 0:
                    beta, phase = T.interp(mene, a), "impulsion"
                else:
                    beta, phase = b_fin + pente * (a - ev["chute"]), "chute et chemin perdu"
                if not np.isfinite(beta):
                    beta = b_deg
                images.append((phi, a, beta, f"{alt.palette} : {phase}, ancre "
                                             f"{rel.depuis_mediane(a):+.2f}°"))
        self.images_stat = images

    def _preparer_images_dynamiques(self):
        r = self.res_dyn
        fen = dynamique.fenetres_action(r)
        if not fen:
            self.images_dyn = None
            return
        n = N_IMAGES // len(fen)
        images = []
        for a, b in fen:
            for t in np.linspace(a, b, n):
                q = dynamique.etat_a(r, t)
                act = ", ".join(dynamique.actifs_a(r, t)) or "aucun contact"
                images.append((q[0], q[1], q[2], f"t = {1000 * t:.3f} ms : {act}"))
        self.images_dyn = images[:N_IMAGES]

    def _images(self):
        if self.chk_dyn.isChecked() and self.images_dyn:
            return self.images_dyn
        return self.images_stat

    def _afficher_image(self, k):
        imgs = self._images()
        if not imgs or self.trace is None:
            return
        k = min(k, len(imgs) - 1)
        phi, a, beta, txt = imgs[k]
        self._poser(phi, a, beta)
        self.lbl_image.setText(txt)

    def _image_suivante(self):
        self.curseur.setValue((self.curseur.value() + 1) % N_IMAGES)

    def _lire(self, oui):
        self.btn_lire.setText("⏸" if oui else "▶")
        (self._minuterie.start if oui else self._minuterie.stop)()

    def _changer_source(self, _oui):
        if self.chk_dyn.isChecked() and not self.images_dyn:
            self.statusBar().showMessage("Lancez d'abord le calcul dynamique (volet 2).")
        self._afficher_image(self.curseur.value())
        self._tracer_courbes()

    # ======================================================================
    # Courbes
    # ======================================================================
    def _tracer_courbes(self):
        if self.rel is None:
            return
        self._courbe_couple()
        self._courbe_cinematique()
        self._courbe_energie()
        self._courbe_retard()

    def _courbe_couple(self):
        """Couple sur l'axe du balancier pendant toute la fonction.

        Convention : positif quand l'ancre est motrice (elle donne de
        l'énergie au balancier), négatif quand le balancier est moteur
        (dégagement), nul hors contact. En haut, une alternance par colonne,
        en fonction du temps, avec les chocs reçus par le balancier ; en bas,
        en fonction de l'angle du balancier, statique et dynamique ensemble."""
        fig = self.fig_couple
        fig.clear()
        par_pal = {a.palette: a for a in self.rel.alternances}
        dyn = dynamique.couple_fonction(self.sim, self.res_dyn) if self.res_dyn is not None else []
        ordre = ["entrée", "sortie"]
        gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 0.55, 1.2])
        coul = {"entrée": theme.ACC.name(), "sortie": theme.ACC2.name()}
        for col, (pal, d) in enumerate(zip(ordre, dyn or [None, None])):
            ax = fig.add_subplot(gs[0, col])
            axc = fig.add_subplot(gs[1, col], sharex=ax)
            if d is None:
                ax.text(0.5, 0.5, "Lancez le calcul dynamique", ha="center", va="center",
                        transform=ax.transAxes)
            else:
                t = 1000 * d["t"]
                C = 1e9 * d["couple"]
                ax.fill_between(t, C, 0, where=C > 0, step="post", color=theme.ACC3.name(),
                                alpha=0.35, label="ancre motrice")
                ax.fill_between(t, C, 0, where=C < 0, step="post", color=theme.ACC4.name(),
                                alpha=0.35, label="balancier moteur")
                ax.step(t, C, where="post", color=coul[pal], lw=1.0)
                if d["chocs"]:
                    tc = np.array([1000 * c[0] for c in d["chocs"]])
                    ic = np.array([1e12 * c[1] for c in d["chocs"]])
                    ml, sl, bl = axc.stem(tc, ic)
                    plt_col = theme.encre()
                    sl.set_color(plt_col)
                    ml.set_color(plt_col)
                    ml.set_markersize(3)
                    bl.set_color(theme.gris(0.5))
            ax.axhline(0, color=theme.gris(0.5), lw=0.8)
            ax.set_title(f"Alternance où travaille la palette d'{pal}" if pal == "entrée"
                         else "Alternance où travaille la palette de sortie", fontsize=9)
            ax.grid(alpha=0.3)
            axc.axhline(0, color=theme.gris(0.5), lw=0.8)
            axc.grid(alpha=0.3)
            axc.set_xlabel("temps depuis le début de l'action (ms)")
            if col == 0:
                ax.set_ylabel("couple sur le balancier (nN·m)")
                axc.set_ylabel("chocs (nN·m·ms)")
                if d is not None:
                    ax.legend(fontsize=7, loc="upper left")
        # en fonction de l'angle
        ax = fig.add_subplot(gs[2, :])
        for e, alt in zip(self.eff_stat, self.rel.alternances):
            ax.plot(np.degrees(e["theta"]), 1e9 * e["couple"] * e["s_b"], "-", lw=1.2,
                    color=coul[alt.palette], alpha=0.6,
                    label=f"statique, palette d'{alt.palette}" if alt.palette == "entrée"
                    else "statique, palette de sortie")
        for pal, d in zip(ordre, dyn):
            ax.plot(np.degrees(d["theta"]), 1e9 * d["couple"], ".", ms=2, color=coul[pal],
                    label="dynamique, palette d'entrée" if pal == "entrée"
                    else "dynamique, palette de sortie")
        ax.axvline(0, color=theme.gris(0.5), lw=0.8)
        ax.axhline(0, color=theme.gris(0.5), lw=0.8)
        ax.set_xlabel("angle du balancier depuis la ligne des centres (°)")
        ax.set_ylabel("couple (nN·m)\n> 0 : ancre motrice")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
        self.can_couple.draw_idle()

    def _courbe_cinematique(self):
        fig = self.fig_cin
        fig.clear()
        axs = fig.subplots(3, 1, sharex=True)
        dyn = self.chk_dyn.isChecked() and self.res_dyn is not None
        if dyn:
            r = self.res_dyn
            for a, b in dynamique.fenetres_action(r):
                k = (r.t >= a) & (r.t <= b)
                t = 1000 * r.t[k]
                axs[0].plot(t, np.degrees(r.y[k, 0]), ".-", ms=2)
                axs[1].plot(t, np.degrees(r.y[k, 1] - self.rel.a_mid), ".-", ms=2)
                axs[2].plot(t, np.degrees(r.y[k, 2] - self.rel.T.beta_ldc), ".-", ms=2)
            axs[2].set_xlabel("temps (ms)")
        else:
            for alt in self.rel.alternances:
                x = self.rel.depuis_mediane(alt.alpha)
                axs[0].plot(x, np.degrees(alt.phi - alt.phi[0]), "-", label=alt.palette)
                axs[1].plot(x, x, "-")
                T = self.rel.T
                ok = np.isfinite(alt.pousse.beta)
                axs[2].plot(self.rel.depuis_mediane(T.alpha[ok]),
                            np.degrees(alt.pousse.beta[ok] - T.beta_ldc), "-", lw=0.8)
            axs[0].legend(fontsize=8)
            axs[2].set_xlabel("angle d'ancre depuis la médiane (°)")
        axs[0].set_ylabel("roue (°)")
        axs[1].set_ylabel("ancre (°)")
        axs[2].set_ylabel("balancier (°)")
        for a in axs:
            a.grid(alpha=0.3)
        self.can_cin.draw_idle()

    def _courbe_energie(self):
        fig = self.fig_en
        fig.clear()
        ax = fig.add_subplot(111)
        s = self.res_stat
        W = s["travail_par_oscillation"]
        E = s["energie_par_oscillation"]
        series = [("statique", [W, E, W - E, 0, 0])]
        if self.res_dyn is not None:
            b = self.res_dyn.bilan
            series.append(("dynamique", [b["travail_moteur"], b["energie_recue_balancier"],
                                         b["frottements"], b["chocs"],
                                         b["energie_cinetique_roue_ancre"]]))
        noms = ["travail moteur", "reçu par le balancier", "frottements", "chocs",
                "énergie cinétique\nroue + ancre"]
        x = np.arange(len(noms))
        w = 0.8 / len(series)
        for k, (lib, v) in enumerate(series):
            ax.bar(x + (k - (len(series) - 1) / 2) * w, 1e9 * np.array(v), w, label=lib)
        ax.set_xticks(x, noms, fontsize=8)
        ax.set_ylabel("énergie par oscillation (nJ)")
        ax.set_title("Bilan d'énergie (en statique, les pertes sont toutes au frottement)",
                     fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, axis="y")
        self.can_en.draw_idle()

    def _courbe_retard(self):
        fig = self.fig_ret
        fig.clear()
        ax = fig.add_subplot(111)
        amps = np.linspace(120, 360, 49)
        par = self.parametres()
        ax.plot(amps, statique.retard(self.eff_stat, par, amps), "-", label="statique (Airy)")
        if self.retards_dyn is not None:
            ax.plot(self.retards_dyn[:, 0], self.retards_dyn[:, 1], "o-", label="dynamique")
        if self.res_dyn is not None:
            b = self.res_dyn.bilan
            ax.plot([b["amplitude_depart_deg"]], [b["marche_s_j"]], "s",
                    color=theme.ACC2.name(), label="calcul dynamique courant")
        ax.axhline(0, color=theme.gris(0.5), lw=0.8)
        ax.set_xlabel("amplitude du balancier (°)")
        ax.set_ylabel("marche due à l'échappement (s/j)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        self.can_ret.draw_idle()

    def _tracer_pareto(self):
        fig = self.fig_par
        fig.clear()
        ax = fig.add_subplot(111)
        sol = self.solutions
        if sol:
            obj = [k for k, cb in self.opti_obj.items() if cb.isChecked()] or ["rendement", "retard"]
            kx, ky = (obj + ["retard", "rendement"])[:2] if len(obj) < 2 else obj[:2]
            fx = lambda m: 100 * m[kx] if kx == "rendement" else (abs(m[kx]) if kx == "retard" else m[kx])
            fy = lambda m: 100 * m[ky] if ky == "rendement" else (abs(m[ky]) if ky == "retard" else m[ky])
            X = [fx(s["mesures"]) for s in sol]
            Y = [fy(s["mesures"]) for s in sol]
            col = [theme.ACC.name() if s["faisable"] else theme.gris(0.5) for s in sol]
            ax.scatter(X, Y, c=col, picker=5)
            k = getattr(self, "_sol_courante", None)
            if k is not None and k < len(sol):
                ax.scatter([X[k]], [Y[k]], s=120, facecolors="none", edgecolors=theme.ACC2.name())
            ax.set_xlabel(OBJECTIFS[kx][0] + (" (%)" if kx == "rendement" else ""))
            ax.set_ylabel(OBJECTIFS[ky][0] + (" (%)" if ky == "rendement" else ""))
        else:
            ax.text(0.5, 0.5, "Aucune optimisation", ha="center", transform=ax.transAxes)
        ax.grid(alpha=0.3)
        self.can_par.draw_idle()

    def _clic_pareto(self, ev):
        if not len(ev.ind):
            return
        k = int(ev.ind[0])
        self.table_sol.selectRow(k)

    def closeEvent(self, e):
        self._arret_opti = True
        super().closeEvent(e)
