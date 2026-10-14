"""Interface graphique : deux onglets à gauche, engrenage et courbes à droite.

Onglet « Création » : paramètres des engrenages, arcs de cercle tangents et export.
Onglet « Simulation » : simulation d'un engrènement défini par des arcs de cercle
(issus de l'onglet Création ou d'un fichier .json / .Eng / .Eno) ou par des suites de
points (roues taillées par l'onglet Fraise).

Usage : python main.py [fichier]   (nécessite PySide6, matplotlib, numpy)
"""

import json
import math
import os
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import NamedTuple

import numpy as np

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import QLocale, Qt, QThread, Signal
from PySide6.QtGui import QAction, QDoubleValidator, QIntValidator, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSlider,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from commun import navigation, theme
from . import arcs
from . import contact
from . import evolvente

RACINE = Path(__file__).resolve().parent
FICHIER_DEFAUT = RACINE / "test.Cev"
NI = 150  # points par demi-dent
SIGMA_DEFAUT = 20.0  # angle de croisement par défaut (degrés) : valeur courante en skiving
MSG_SCARVING = (
    "Le fichier .Frs est le format Fortran de la fraise-mère (crémaillère) : il n'a pas de place "
    "pour le nombre de dents d'un outil de scarving. Enregistrez plutôt le profil taillé "
    "(.Eno / .dxf), ou repassez l'outil concerné à Z outil = 0.")
MSG_ARCS_REQUIS = (
    "Le profil courant est une suite de points (roues taillées). Reconvertissez-le en arcs de "
    "cercle dans l'onglet Fraise (groupe « reconvertir en arcs »).")
N_DENTS = 3  # dents dessinées par roue
N_POS = contact.N_POS  # positions de la menante sur une période d'engrènement
ROUES = ("Menante", "Menée")


@dataclass(frozen=True)
class Champ:
    cle: str  # nom de l'attribut de evolvente.Roue (ou de evolvente.Projet)
    nom: str  # nom de la variable
    aide: str
    mini: float
    maxi: float
    entier: bool = False
    exclus: tuple = None  # (bas, haut) : intervalle interdit à l'intérieur de [mini, maxi]


# Une ligne du panneau par champ : pour ajouter une variable de roue, ajouter un Champ ici
# (et l'attribut correspondant dans evolvente.Roue).
CHAMPS = (
    Champ("mo", "MO", "Module (mm)", 0.001, 1000),
    Champ("zz", "ZZ", "Nombre de dents (négatif = couronne intérieure)", -10000, 10000,
          entier=True, exclus=(-2, 2)),
    Champ("nw", "NW", "Angle de pression (°)", 1, 60),
    Champ("qe", "QE", "Rapport creux / plein au diamètre primitif", 0.1, 10),
    Champ("qg", "QG", "Raccord de pied (x module), 0 = aucun", 0, 5),
    Champ("ma", "MA", "Hauteur de développante sous le primitif (x module)", 0, 5),
    Champ("me", "ME", "Hauteur de développante au-dessus du primitif (x module)", 0, 5),
)
# Export en arcs de cercle (.Eno / .json) et matériau (données de Hertz), par roue.
CHAMPS_EXPORT = (
    Champ("narcs", "Na", "Arcs de cercle sur la développante", 1, arcs.NARCS_MAX, entier=True),
    Champ("e_mod", "E", "Module d'élasticité (N/mm²)", 1, 1e7),
    Champ("poisson", "ν", "Coefficient de Poisson", 0, 0.5),
    Champ("epaisseur", "H", "Épaisseur de la roue (mm)", 0.001, 1000),
)
# Valeurs des données de Hertz quand le fichier n'en contient pas (exemples de Files/).
HERTZ_DEFAUT = ((100000.0, 0.3, 0.3), (200000.0, 0.3, 0.28))
# Variables de la simulation (communes aux deux roues, enregistrées dans le projet).
CHAMPS_SIM = (
    Champ("mu", "μ", "Coefficient de frottement", 0, 1),
    Champ("mt", "Cm", "Couple sur la menante (mN·m)", 1e-6, 1e6),
    Champ("da", "ΔA", "Variation de l'entraxe simulée (µm)", -1000, 1000),
    Champ("da_neg", "ΔA min", "Plage d'entraxe étudiée, minimum (µm)", -1000, 0),
    Champ("da_pos", "ΔA max", "Plage d'entraxe étudiée, maximum (µm)", 0, 1000),
)

# Critères d'optimisation : (clé, libellé, nom dans les fichiers de l'ancien programme)
OBJECTIFS = (
    ("couple", "Variation de couple mini", "couple"),
    ("jeu", "Jeu mini", "jeux"),
    ("rendement", "Rendement maxi", "rendement"),
    ("hertz", "Pression de Hertz mini", "PHZ"),
)
CHAMPS_OPTI = (
    Champ("ponderation", "Poids", "Poids du 2e paramètre par rapport au 1er : J = f1/f1₀ + poids·f2/f2₀", 0, 1e6),
    Champ("deplacement_um", "Étendue", "Écart typique de profil exploré autour du départ (µm)", 0.5, 500),
    Champ("jeu_min_admis", "Jeu min", "Jeu minimal de la menée à l'entraxe minimal : évite le blocage (°)", 0, 90),
    Champ("jeu_max_admis", "Jeu max", "Jeu maximal de la menée à l'entraxe maximal (°) ; 0 : jeu du profil de départ", 0, 360),
    Champ("r_tete_min", "R tête min", "Rayon minimal de l'arc de tête (mm), durée de vie des fraises ; 0 : pas de limite", 0, 10),
    Champ("population", "Population", "Taille de la population de CMA-ES", 4, 10000, entier=True),
    Champ("evaluations", "Évaluations", "Évaluations par départ", 100, 10**7, entier=True),
    Champ("departs", "Départs", "Nombre de départs (profil courant, puis profils tirés de plus en plus loin)", 1, 1000, entier=True),
    Champ("graine", "Graine", "Graine aléatoire (résultats reproductibles)", 0, 2**31 - 1, entier=True),
)

# Onglet Fraise : création (0 = automatique, voir fraise.periodes_defaut) et taillage.
CHAMPS_FRAISE_CREATION = (
    Champ("n_periodes", "Périodes", "Pas de denture voisins pris en compte (impair ; 0 = auto. selon Z)", 0, 99, entier=True),
    Champ("n_positions", "Positions", "Positions par période (plus grand = plus précis, plus lent)", 5, 2000, entier=True),
    Champ("pc", "Période", "Période de la fraise / période standard π·module (1 = standard)", 0.1, 10),
    Champ("n_points", "Points", "Nombre de points du profil de fraise", 10, 2000, entier=True),
)
CHAMPS_FRAISE_TAILLAGE = (
    Champ("zz_cible", "Z cible", "Nombre de dents de la roue à tailler", 3, 9999, entier=True),
    Champ("oy_um", "Déport Oy", "Déport de la fraise (µm ; positif = fraise éloignée du centre)", -1e6, 1e6),
    Champ("module_cible", "Module cible", "Module voulu (mm) ; prioritaire sur le déport si différent de celui de la fraise (0 : ignoré)", 0, 1000),
    Champ("n_periodes_t", "Périodes", "Pas de denture voisins (0 = auto. selon Z)", 0, 99, entier=True),
    Champ("n_positions_t", "Positions", "Positions par période", 5, 2000, entier=True),
    Champ("n_points_t", "Points", "Nombre de points du profil taillé", 10, 2000, entier=True),
)

STYLE_ERREUR = f"QLineEdit {{ border: 1px solid {theme.ERREUR}; background: #4a2226; }}"
FILTRE_ARCS = "Arcs de cercle (*.json *.Eng *.eng *.Eno *.eno);;Tous les fichiers (*)"


class ProfilPoints(NamedTuple):
    """Demi-profil donné en suite de points, pas en arcs de cercle.

    C'est ce que produit réellement une fraise de taillage : l'enveloppe d'un balayage, pas une
    chaîne d'arcs. Pour simuler l'engrènement, ces points sont la forme fidèle ; l'approximation
    en arcs n'est nécessaire que pour repartir en optimisation ou pour écrire un `.Eno`.
    """

    nom: str
    mo: float
    zz: int
    points: object  # (n,2) demi-profil, repère de evolvente.evolve
    e_mod: float = None
    poisson: float = None
    epaisseur: float = None


def texte_nombre(v):
    return f"{v:.6f}".rstrip("0").rstrip(".") or "0"


def _texte_ecart_arcs(profil_arcs):
    """Écart d'un ajustement en arcs : sur tout le profil et, quand il a été mesuré, sur le seul
    flanc actif -- c'est celui-là qui compte pour l'engrènement, le congé de pied n'étant pas
    une surface de travail (voir `arcs.approximer`)."""
    t = f"{profil_arcs.ecart_max * 1e3:.2f} µm"
    if profil_arcs.ecart_flanc is not None:
        t += f" ({profil_arcs.ecart_flanc * 1e3:.2f} µm sur le flanc actif)"
    return t


def _champ_saisie(c, info):
    le = QLineEdit()
    if c.entier:
        v = QIntValidator(int(c.mini), int(c.maxi), le)
    else:
        v = QDoubleValidator(c.mini, c.maxi, 6, le)
        v.setNotation(QDoubleValidator.StandardNotation)
    v.setLocale(QLocale.c())  # séparateur décimal : le point
    le.setValidator(v)
    le.setAlignment(Qt.AlignRight)
    le.setToolTip(f"{info} : {c.aide} ({c.mini:g} … {c.maxi:g})")
    return le


def _etiq(texte):
    """Étiquette qui passe à la ligne (sinon un texte long élargit tout le panneau)."""
    e = QLabel(texte)
    e.setWordWrap(True)
    return e


def _defilant(widget):
    zone = QScrollArea()
    zone.setWidgetResizable(True)
    zone.setWidget(widget)
    return zone


def _json_est_arcs(chemin):
    """Un .json d'arcs de cercle (schéma profile.json ou .Eng), pas un projet de création."""
    try:
        doc = json.loads(Path(chemin).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return isinstance(doc, list) or (isinstance(doc, dict) and doc.get("format") != evolvente.FORMAT_JSON
                                     and "m" in doc and "s" in doc)


def _libelle_violation(cle):
    from . import recherche
    return recherche.LIBELLES.get(cle, "profil de départ inutilisable")


class Travail(QThread):
    """Recherche d'optima (`recherche.rechercher`) dans un fil à part ; l'interface reste utilisable."""

    avance = Signal(dict)
    fini = Signal(object, dict)
    echec = Signal(str)

    def __init__(self, prep):
        super().__init__()
        self.prep = prep
        self._stop = False

    def arreter(self):
        self._stop = True

    def run(self):
        try:
            from . import recherche
            roues, mobiles, mu, da_neg, da_pos, objectif, options, libres = self.prep
            # la résolution de "jeu_max=auto" se fait dans rechercher() : elle a besoin du profil de
            # référence, qui peut ne pas être le profil courant si celui-ci n'est pas valable (voir
            # _jeu_max_auto, utilisé lui tel quel par le bouton "Évaluer" qui décrit ce profil précis)
            liste, stats = recherche.rechercher(roues, mobiles, mu, da_neg, da_pos, objectif, options, libres,
                                                rappel=self.avance.emit, arret=lambda: self._stop)
            self.fini.emit(liste, stats)
        except Exception as e:  # noqa: BLE001 - toute erreur doit revenir à l'interface
            self.echec.emit(f"{type(e).__name__} : {e}")


class Fenetre(QMainWindow):
    def __init__(self, chemin=None):
        super().__init__()
        self.chemin = None  # fichier .json de travail (projet de création ; None : jamais enregistré)
        self.source = None  # ancien fichier .Cev importé, s'il y en a un
        self.projet = evolvente.Projet(roues=[evolvente.Roue(0.135, 102, 20, 1.2, 0, 1, 0.5),
                                              evolvente.Roue(0.135, 10, 20, 1.2, 0.2, 1, 1)])
        self.modifie = False
        # onglet Création
        self.roues = self.profils = self.limites_c = None
        self.arcs_roues = self.arcs_echantillons = None
        self._cache_arcs = {}  # profil en arcs par (paramètres de la roue, nombre d'arcs)
        # onglet Simulation
        self.mobiles = None  # (menant, mené) : arcs.Mobile
        # profil courant donné en suite de points au lieu d'arcs de cercle (roues taillées) :
        # (menant, mené) de ProfilPoints, ou None. Les deux sont exclusifs.
        self.mobiles_points = None
        self.extras = None  # paramètres du fichier d'arcs conservés pour la suite
        self.origine = ""
        self.sim = self.sim_roues = self.sim_profils = self.sim_limites = self.superposition = None
        self.sim_profils_dos = None  # demi-profils arrière échantillonnés (dents dissymétriques) ou None
        self.da_mm, self.mt = 0.0, 1.0
        # onglet Fraise
        self.fraise_paire = None  # (outil menante, outil menée) : ProfilFraise (crémaillère) ou ProfilOutil
        self.z_outils = [0, 0]  # nombre de dents de l'outil par roue (0 = fraise-mère)
        self.sigmas = [0.0, 0.0]  # angle de croisement des axes par roue (rad)
        self.fraise_source = None  # (menant, mené) : arcs.Mobile qui ont servi à les calculer, ou None
        self.taillage = [None, None]  # par roue : {"zp", "ecart", "zz", "oy"} ou None
        self._vue_ok = False  # la vue du graphique a été fixée (zoom de l'utilisateur conservé)
        self.resize(1400, 880)

        self._creer_menu()
        gauche = self._creer_gauche()
        droite = self._creer_graphiques()
        cadre = QSplitter(Qt.Horizontal)
        cadre.addWidget(gauche)
        cadre.addWidget(droite)
        cadre.setStretchFactor(1, 1)
        cadre.setSizes([540, 860])
        self.setCentralWidget(cadre)

        self.remplir_panneau()
        chemin = Path(chemin) if chemin else None
        est_arcs = chemin is not None and (
            chemin.suffix.lower() in (".eng", ".eno") or (chemin.suffix.lower() == ".json" and _json_est_arcs(chemin)))
        if chemin is not None and not est_arcs and self.charger(chemin):
            pass
        else:  # onglet Création : fichier par défaut ; puis fichier d'arcs éventuel
            self._charger_defaut()
            if est_arcs:
                self.charger_arcs(chemin)
        self._titre()

    def _charger_defaut(self):
        """Onglet Création construit d'après test.Cev (ou les valeurs par défaut)."""
        if FICHIER_DEFAUT.exists():
            self.projet, self.source = evolvente.lire_projet(FICHIER_DEFAUT), FICHIER_DEFAUT
        self.remplir_panneau()
        self.reconstruire()

    # --- interface -------------------------------------------------------
    def _creer_menu(self):
        menu = self.menuBar().addMenu("&File")
        for texte, raccourci, action in (
            ("&Open File…", QKeySequence.Open, self.ouvrir),
            ("&Save", QKeySequence.Save, self.sauver),
            ("Save &As…", QKeySequence.SaveAs, self.sauver_sous),
        ):
            a = QAction(texte, self)
            a.setShortcut(raccourci)
            a.triggered.connect(action)
            menu.addAction(a)
        menu.addSeparator()
        export = QAction("&Export Arcs (.json / .Eno)…", self)
        export.setShortcut("Ctrl+E")
        export.triggered.connect(self.exporter_arcs)
        menu.addAction(export)
        rapport_action = QAction("Exporter un &rapport (.docx)…", self)
        rapport_action.triggered.connect(self.exporter_rapport)
        menu.addAction(rapport_action)
        dxf_action = QAction("Exporter les profils en &DXF…", self)
        dxf_action.triggered.connect(self.exporter_profils_dxf)
        menu.addAction(dxf_action)
        menu.addSeparator()
        quitter = QAction("&Quit", self)
        quitter.setShortcut(QKeySequence.Quit)
        quitter.triggered.connect(self.close)
        menu.addAction(quitter)
        lancer = QAction("Reconstruire / Simuler", self)
        lancer.setShortcut("F5")
        lancer.triggered.connect(self.actualiser)
        menu_run = self.menuBar().addMenu("&Run")
        menu_run.addAction(lancer)
        self.action_optimiser = QAction("&Optimiser", self)
        self.action_optimiser.setShortcut("F6")
        self.action_optimiser.triggered.connect(self.lancer_optimisation)
        self.action_arret = QAction("&Arrêter l'optimisation", self)
        self.action_arret.setEnabled(False)
        self.action_arret.triggered.connect(self.arreter_optimisation)
        menu_run.addAction(self.action_optimiser)
        menu_run.addAction(self.action_arret)

    def _creer_gauche(self):
        self.edits = ({}, {})
        self.edits_sim = {}
        self.onglets = QTabWidget()
        self.onglets.addTab(self._page_creation(), "Création")
        self.onglets.addTab(self._page_simulation(), "Simulation")
        self.onglets.addTab(self._page_optimisation(), "Optimisation")
        self.onglets.addTab(self._page_fraise(), "Fraise")
        self.onglets.currentChanged.connect(self._changer_onglet)
        self.statut = QLabel()
        self.statut.setWordWrap(True)
        gauche = QWidget()
        lay = QVBoxLayout(gauche)
        lay.addWidget(self.onglets, 1)
        lay.addWidget(self.statut)
        gauche.setMinimumWidth(420)
        return gauche

    def _groupe_roues(self, titre, champs):
        """Groupe de champs : une colonne par roue, une ligne par paramètre."""
        groupe = QGroupBox(titre)
        grille = QGridLayout(groupe)
        grille.setColumnStretch(1, 1)
        grille.setColumnStretch(2, 1)
        for k, nom in enumerate(ROUES):
            en_tete = QLabel(f"<b>{nom}</b>")
            en_tete.setAlignment(Qt.AlignCenter)
            grille.addWidget(en_tete, 0, 1 + k)
        for i, c in enumerate(champs, start=1):
            etiquette = QLabel(f"<b>{c.nom}</b>  {c.aide}")
            etiquette.setWordWrap(True)
            grille.addWidget(etiquette, i, 0)
            for k in range(2):
                le = _champ_saisie(c, f"{c.nom} {ROUES[k].lower()}")
                le.textEdited.connect(lambda _t, le=le: self._edite(le))
                grille.addWidget(le, i, 1 + k)
                self.edits[k][c.cle] = le
        return groupe

    def _page_creation(self):
        contenu = QWidget()
        lay = QVBoxLayout(contenu)
        lay.addWidget(self._groupe_roues("Profils", CHAMPS))
        lay.addWidget(self._groupe_roues("Arcs de cercle et matériau", CHAMPS_EXPORT))
        cadre_arcs = QGroupBox("Arcs de cercle tangents")
        lay_arcs = QVBoxLayout(cadre_arcs)
        self.case_arcs = QCheckBox("Arcs (pointillés) et points de tangence")
        self.case_arcs.setChecked(True)
        self.case_arcs.toggled.connect(lambda _c: self._dessiner_creation())
        self.info_arcs = QLabel("-")
        self.info_arcs.setWordWrap(True)
        self.info_arcs.setTextFormat(Qt.RichText)
        lay_arcs.addWidget(self.case_arcs)
        lay_arcs.addWidget(self.info_arcs)
        lay.addWidget(cadre_arcs)
        lay.addStretch(1)

        self.bouton = QPushButton("Reconstruire l'engrenage  (F5)")
        self.bouton.setMinimumHeight(36)
        self.bouton.clicked.connect(self.reconstruire)
        bouton_export = QPushButton("Exporter les arcs (.json / .Eno)…")
        bouton_export.clicked.connect(self.exporter_arcs)

        page = QWidget()
        lay_page = QVBoxLayout(page)
        lay_page.setContentsMargins(0, 0, 0, 0)
        lay_page.addWidget(_defilant(contenu), 1)
        lay_page.addWidget(self.bouton)
        lay_page.addWidget(bouton_export)
        return page

    def _page_simulation(self):
        contenu = QWidget()
        lay = QVBoxLayout(contenu)

        source = QGroupBox("Profil simulé")
        lay_source = QVBoxLayout(source)
        self.etiq_origine = QLabel("Aucun fichier d'arcs.")
        self.etiq_origine.setWordWrap(True)
        self.etiq_origine.setTextFormat(Qt.RichText)
        lay_source.addWidget(self.etiq_origine)
        b_creation = QPushButton("Utiliser les arcs de l'onglet Création")
        b_creation.clicked.connect(self.utiliser_arcs_creation)
        b_ouvrir = QPushButton("Ouvrir un fichier d'arcs…")
        b_ouvrir.clicked.connect(self.ouvrir_arcs)
        lay_source.addWidget(b_creation)
        lay_source.addWidget(b_ouvrir)
        lay.addWidget(source)

        self.groupe_modif = modif = QGroupBox("Modifier les arcs (sans reprendre la création)")
        lay_modif = QVBoxLayout(modif)
        lay_modif.addWidget(_etiq(
            "Sélectionnez un arc : « Ajouter » le scinde en deux sur le même cercle (la forme ne "
            "change pas, un arc de plus pour l'optimisation) ; « Supprimer » fusionne son point de "
            "départ avec l'arc précédent (la forme change après ce point)."))
        self.liste_arcs_sim = []
        for nom in ROUES:
            lay_modif.addWidget(QLabel(f"<b>{nom}</b>"))
            liste = QListWidget()
            liste.setMaximumHeight(110)
            liste.itemSelectionChanged.connect(lambda l=liste: self._exclure_autres_listes(l))
            self.liste_arcs_sim.append(liste)
            lay_modif.addWidget(liste)
        boutons_arcs = QHBoxLayout()
        b_ajouter = QPushButton("Ajouter un point")
        b_ajouter.clicked.connect(self.ajouter_point_arc)
        b_retirer = QPushButton("Supprimer ce point")
        b_retirer.clicked.connect(self.retirer_point_arc)
        boutons_arcs.addWidget(b_ajouter)
        boutons_arcs.addWidget(b_retirer)
        lay_modif.addLayout(boutons_arcs)
        lay.addWidget(modif)

        sim = QGroupBox("Simulation")
        grille_sim = QGridLayout(sim)
        grille_sim.setColumnStretch(0, 1)
        for i, c in enumerate(CHAMPS_SIM):
            etiquette = QLabel(f"<b>{c.nom}</b>  {c.aide}")
            etiquette.setWordWrap(True)
            grille_sim.addWidget(etiquette, i, 0)
            le = _champ_saisie(c, c.nom)
            le.setMaximumWidth(110)
            le.textEdited.connect(lambda _t, le=le: self._edite(le))
            grille_sim.addWidget(le, i, 1)
            self.edits_sim[c.cle] = le
        lay.addWidget(sim)

        self.bouton_sim = QPushButton("Simuler  (F5)")
        self.bouton_sim.setMinimumHeight(36)
        self.bouton_sim.clicked.connect(self.simuler)
        lay.addWidget(self.bouton_sim)

        # curseur des positions d'engrènement
        self.curseur = QSlider(Qt.Horizontal)
        self.curseur.setRange(0, N_POS - 1)
        self.curseur.setEnabled(False)
        self.curseur.valueChanged.connect(self.afficher_position)
        self.etiq_pos = _etiq("-")
        cadre_pos = QGroupBox("Position d'engrènement")
        lay_pos = QVBoxLayout(cadre_pos)
        lay_pos.addWidget(self.etiq_pos)
        lay_pos.addWidget(self.curseur)
        self.resultats = QLabel()
        self.resultats.setWordWrap(True)
        self.resultats.setTextFormat(Qt.RichText)
        lay_pos.addWidget(self.resultats)
        lay.addWidget(cadre_pos)
        self.stats = QLabel()
        self.stats.setWordWrap(True)
        self.stats.setTextFormat(Qt.RichText)
        cadre_stats = QGroupBox("Variation sur la période (100 % = idéal)")
        QVBoxLayout(cadre_stats).addWidget(self.stats)
        lay.addWidget(cadre_stats)
        lay.addStretch(1)
        return _defilant(contenu)

    def _page_optimisation(self):
        contenu = QWidget()
        lay = QVBoxLayout(contenu)
        self.edits_opti = {}

        obj = QGroupBox("Paramètres à optimiser")
        grille = QGridLayout(obj)
        self.combo_principal, self.combo_secondaire = QComboBox(), QComboBox()
        for combo in (self.combo_principal, self.combo_secondaire):
            combo.setMinimumWidth(190)
        self.combo_secondaire.addItem("Aucun", "")
        for cle, libelle, _ in OBJECTIFS:
            self.combo_principal.addItem(libelle, cle)
            self.combo_secondaire.addItem(libelle, cle)
        for combo in (self.combo_principal, self.combo_secondaire):
            combo.currentIndexChanged.connect(lambda _i: self._marquer_modifie())
        grille.addWidget(_etiq("<b>Principal</b>"), 0, 0)
        grille.addWidget(self.combo_principal, 0, 1)
        grille.addWidget(_etiq("<b>Secondaire</b>"), 1, 0)
        grille.addWidget(self.combo_secondaire, 1, 1)
        grille.setColumnStretch(0, 1)
        self._ajouter_champs(grille, CHAMPS_OPTI[:1], 2)
        # Troisième terme, facultatif, de l'objectif : l'équilibre de la menée de part et d'autre
        # de la ligne des centres. Sans lui l'optimisation va volontiers vers beaucoup de menée
        # après et peu avant, ce qui n'est pas toujours voulu -- d'où la case plutôt qu'un terme
        # imposé.
        self.case_equilibre = QCheckBox("Équilibrer la menée avant / après la ligne des centres")
        self.case_equilibre.setToolTip(
            "Ajoute à l'objectif l'écart entre les deux étendues de la menée de part et d'autre du "
            "point primitif. Coché, l'optimisation cherche à les égaliser ; décoché, elle les "
            "laisse libres et donne souvent beaucoup plus de menée après la ligne des centres.")
        self.case_equilibre.stateChanged.connect(lambda _s: self._marquer_modifie())
        grille.addWidget(self.case_equilibre, 3, 0, 1, 2)
        lay.addWidget(obj)

        contr = QGroupBox("Contraintes (une solution qui les viole n'est pas valable)")
        grille_c = QGridLayout(contr)
        grille_c.setColumnStretch(0, 1)
        self._ajouter_champs(grille_c, CHAMPS_OPTI[2:5], 0)
        grille_c.addWidget(_etiq("Pas d'interférence entre les dents : toujours contrôlé sur les contours "
                                 "complets (les deux flancs, tête et pied) à tous les entraxes."), 3, 0, 1, 2)
        self.case_dissym = QCheckBox("Réparer par le demi-profil arrière (dents dissymétriques)")
        self.case_dissym.setToolTip("Si le profil optimisé viole une contrainte, les flancs arrière des dents des deux roues "
                                    "sont modifiés pour la respecter ; le flanc menant reste inchangé. "
                                    "Export en .json seulement (flag sym: false).")
        self.case_dissym.stateChanged.connect(lambda _s: self._marquer_modifie())
        grille_c.addWidget(self.case_dissym, 4, 0, 1, 2)
        lay.addWidget(contr)

        var = QGroupBox("Variables : arcs de cercle libres")
        lay_var = QVBoxLayout(var)
        lay_var.addWidget(_etiq("Cochez les arcs dont la courbure et la longueur varient. Le premier point "
                                "reste fixe ; le dernier arc est refermé automatiquement sur l'axe de la dent."))
        self.liste_points = []
        for nom in ROUES:
            lay_var.addWidget(QLabel(f"<b>{nom}</b>"))
            liste = QListWidget()
            liste.setMaximumHeight(120)
            liste.itemChanged.connect(lambda _it: self._marquer_modifie())
            self.liste_points.append(liste)
            lay_var.addWidget(liste)
        lay.addWidget(var)

        algo = QGroupBox("Recherche (CMA-ES avec redémarrages)")
        grille_algo = QGridLayout(algo)
        grille_algo.setColumnStretch(0, 1)
        self._ajouter_champs(grille_algo, CHAMPS_OPTI[1:2] + CHAMPS_OPTI[5:], 0)
        lay.addWidget(algo)

        ligne = QHBoxLayout()
        self.bouton_evaluer = QPushButton("Évaluer les critères et le temps")
        self.bouton_evaluer.clicked.connect(self.evaluer_criteres)
        self.bouton_optimiser = QPushButton("Optimiser")
        self.bouton_optimiser.setStyleSheet("font-weight: bold;")
        self.bouton_optimiser.setToolTip("Lance la recherche (menu Run > Optimiser, F6)")
        self.bouton_optimiser.clicked.connect(self.lancer_optimisation)
        self.bouton_arret = QPushButton("Arrêter")
        self.bouton_arret.setEnabled(False)
        self.bouton_arret.clicked.connect(self.arreter_optimisation)
        for b in (self.bouton_evaluer, self.bouton_optimiser, self.bouton_arret):
            ligne.addWidget(b)
        lay.addLayout(ligne)
        self.barre_opti = QProgressBar()
        self.barre_opti.setVisible(False)
        lay.addWidget(self.barre_opti)
        self.info_opti = QLabel("-")
        self.info_opti.setWordWrap(True)
        self.info_opti.setTextFormat(Qt.RichText)
        lay.addWidget(self.info_opti)

        res = QGroupBox("Optima trouvés")
        lay_res = QVBoxLayout(res)
        self.table_opti = QTableWidget(0, 7)
        self.table_opti.setHorizontalHeaderLabels(
            ["J", "Couple %", "Hertz", "Jeu °", "Équilibre", "Dents", "État"])
        self.table_opti.setSelectionBehavior(QTableWidget.SelectRows)
        self.table_opti.setSelectionMode(QTableWidget.SingleSelection)
        self.table_opti.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table_opti.setMinimumHeight(170)
        self.table_opti.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_opti.horizontalHeader().setStretchLastSection(True)
        self.table_opti.itemDoubleClicked.connect(lambda _it: self.afficher_optimum())
        lay_res.addWidget(self.table_opti)
        boutons = QHBoxLayout()
        b_aff = QPushButton("Afficher (onglet Simulation)")
        b_aff.clicked.connect(self.afficher_optimum)
        b_exp = QPushButton("Exporter…")
        b_exp.clicked.connect(self.exporter_optimum)
        boutons.addWidget(b_aff)
        boutons.addWidget(b_exp)
        lay_res.addLayout(boutons)
        lay_res.addWidget(_etiq("Couple : variation aux 3 entraxes ; Hertz : pression maximale (N/mm² pour 1 mN·m) ; "
                                "jeu de la menée de l'entraxe minimal à l'entraxe maximal."))
        lay.addWidget(res)
        lay.addStretch(1)
        self.optima, self.travail = [], None
        return _defilant(contenu)

    def _ajouter_champs(self, grille, champs, ligne0):
        for i, c in enumerate(champs, start=ligne0):
            grille.addWidget(_etiq(f"<b>{c.nom}</b>  {c.aide}"), i, 0)
            self.edits_opti[c.cle] = self._saisie_opti(c)
            grille.addWidget(self.edits_opti[c.cle], i, 1)

    def _saisie_opti(self, c):
        le = _champ_saisie(c, c.nom)
        le.setMaximumWidth(110)
        le.textEdited.connect(lambda _t, le=le: self._edite(le))
        return le

    def _page_fraise(self):
        contenu = QWidget()
        lay = QVBoxLayout(contenu)
        self.edits_fraise = {}

        defauts_fraise = {"n_periodes": "0", "n_positions": "101", "pc": "1", "n_points": "200",
                          "zz_cible": "20", "oy_um": "0", "module_cible": "0",
                          "n_periodes_t": "0", "n_positions_t": "101", "n_points_t": "200"}

        def champs(groupe, table):
            grille = QGridLayout(groupe)
            grille.setColumnStretch(0, 1)
            for i, c in enumerate(table):
                grille.addWidget(_etiq(f"<b>{c.nom}</b>  {c.aide}"), i, 0)
                le = _champ_saisie(c, c.nom)
                le.setMaximumWidth(110)
                le.setText(defauts_fraise[c.cle])
                self.edits_fraise[c.cle] = le
                grille.addWidget(le, i, 1)

        creer = QGroupBox("Créer l'outil (roulement sans glissement roue / outil)")
        lay_creer = QVBoxLayout(creer)
        lay_creer.addWidget(_etiq("Le profil de l'outil ne change pas la forme de la dent : c'est le "
                                  "profil conjugué exact, dans le repère de l'outil."))
        ligne = QHBoxLayout()
        ligne.addWidget(QLabel("Roue affichée :"))
        self.combo_fraise_roue = QComboBox()
        self.combo_fraise_roue.addItems(ROUES)
        self.combo_fraise_roue.currentIndexChanged.connect(self._roue_fraise_changee)
        ligne.addWidget(self.combo_fraise_roue)
        ligne.addStretch(1)
        lay_creer.addLayout(ligne)
        # Outil par roue : une couronne ne peut être taillée que par un outil qui passe dans son
        # ouverture, donc par scarving ; un pignon se taille à la fraise-mère. 0 = crémaillère.
        lay_creer.addWidget(_etiq("Outil par roue : nombre de dents de l'outil de scarving, "
                                  "<b>0 = fraise-mère</b> (crémaillère). Une couronne intérieure "
                                  "exige un outil (Z ≥ 3) qui passe dans son ouverture. Σ est "
                                  "l'angle de croisement des axes : c'est lui qui donne la vitesse "
                                  "de coupe, <b>Σ = 0 ne coupe pas</b> (c'est du mortaisage)."))
        self.edits_z_fraise, self.edits_sigma = [], []
        for champ, liste, defaut, info in (
                ("Z outil", self.edits_z_fraise, "0",
                 "0 : fraise-mère (crémaillère). Sinon nombre de dents de l'outil de scarving."),
                ("Σ (°)", self.edits_sigma, f"{SIGMA_DEFAUT:g}",
                 "Angle de croisement des axes, en degrés. 0 = axes parallèles (mortaisage).")):
            ligne_o = QHBoxLayout()
            ligne_o.addWidget(QLabel(f"{champ} :"))
            for k in range(2):
                ligne_o.addWidget(QLabel(ROUES[k]))
                le = QLineEdit(defaut)
                le.setMaximumWidth(70)
                le.setAlignment(Qt.AlignRight)
                le.setToolTip(info)
                if liste is self.edits_z_fraise:
                    le.setValidator(QIntValidator(0, 9999, le))
                else:
                    le.setValidator(QDoubleValidator(0.0, 60.0, 2, le))
                liste.append(le)
                ligne_o.addWidget(le)
            ligne_o.addStretch(1)
            lay_creer.addLayout(ligne_o)
        sous = QGroupBox()
        champs(sous, CHAMPS_FRAISE_CREATION)
        lay_creer.addWidget(sous)
        b_calc = QPushButton("Calculer les outils")
        b_calc.clicked.connect(self.calculer_fraise)
        lay_creer.addWidget(b_calc)
        self.info_fraise = QLabel("-")
        self.info_fraise.setWordWrap(True)
        self.info_fraise.setTextFormat(Qt.RichText)
        lay_creer.addWidget(self.info_fraise)
        lay.addWidget(creer)

        fichier = QGroupBox("Fichier fraise (.Frs / .dxf)")
        lay_fichier = QHBoxLayout(fichier)
        b_save = QPushButton("Enregistrer .Frs…")
        b_save.clicked.connect(self.sauver_fraise)
        b_open = QPushButton("Ouvrir…")
        b_open.clicked.connect(self.ouvrir_fraise)
        b_dxf = QPushButton("Exporter .dxf…")
        b_dxf.clicked.connect(self.exporter_fraise_dxf)
        lay_fichier.addWidget(b_save)
        lay_fichier.addWidget(b_open)
        lay_fichier.addWidget(b_dxf)
        lay.addWidget(fichier)

        tailler = QGroupBox("Tailler un mobile par génération")
        lay_tailler = QVBoxLayout(tailler)
        sous2 = QGroupBox()
        champs(sous2, CHAMPS_FRAISE_TAILLAGE)
        lay_tailler.addWidget(sous2)
        lay_tailler.addWidget(_etiq("Le module cible n'est pris en compte que s'il diffère du module de "
                                    "la fraise ; sinon c'est le déport qui compte (0 : denture standard) -- "
                                    "utilisés seulement par « Tailler cette roue »."))
        boutons_t = QHBoxLayout()
        b_tailler = QPushButton("Tailler cette roue  (F5)")
        b_tailler.setToolTip("Utilise Z cible / déport ci-dessus, pour la roue choisie plus haut.")
        b_tailler.clicked.connect(self.tailler_mobile)
        b_tailler_2 = QPushButton("Tailler les deux roues (Z et module d'origine)")
        b_tailler_2.setToolTip("Cas le plus courant : reproduit chaque roue avec son propre nombre de "
                               "dents, sans déport -- ignore les champs Z cible / déport ci-dessus.")
        b_tailler_2.clicked.connect(self.tailler_les_deux)
        boutons_t.addWidget(b_tailler)
        boutons_t.addWidget(b_tailler_2)
        lay_tailler.addLayout(boutons_t)
        self.info_taillage = QLabel("-")
        self.info_taillage.setWordWrap(True)
        self.info_taillage.setTextFormat(Qt.RichText)
        lay_tailler.addWidget(self.info_taillage)
        lay.addWidget(tailler)

        points_t = QGroupBox("Roue(s) taillée(s) : utiliser telles quelles (suite de points)")
        lay_points_t = QVBoxLayout(points_t)
        lay_points_t.addWidget(_etiq("La fraise engendre une enveloppe, pas des arcs de cercle : la suite "
                                     "de points est la forme fidèle. C'est ce qu'il faut pour simuler."))
        b_sim_points = QPushButton("Simuler les roues taillées (points)")
        b_sim_points.setToolTip("Charge les profils taillés tels quels comme profil courant et simule "
                                "l'engrènement entre lignes brisées ; la roue non taillée garde le sien.")
        b_sim_points.clicked.connect(self.utiliser_taillage_points)
        lay_points_t.addWidget(b_sim_points)
        b_dxf_points = QPushButton("Exporter les profils en DXF (polyligne)…")
        b_dxf_points.setToolTip("Contour complet de chaque roue taillée, en polyligne (LWPOLYLINE).")
        b_dxf_points.clicked.connect(self.exporter_taillage_points_dxf)
        lay_points_t.addWidget(b_dxf_points)
        lay.addWidget(points_t)

        arcs_t = QGroupBox("Roue(s) taillée(s) : reconvertir en arcs")
        lay_arcs_t = QVBoxLayout(arcs_t)
        lay_arcs_t.addWidget(_etiq("Approxime le(s) profil(s) pratique(s) (taillés) par des arcs de cercle, "
                                   "comme dans l'onglet Création. Nécessaire seulement pour repartir en "
                                   "optimisation ou pour écrire un .Eno."))
        ligne_na = QHBoxLayout()
        ligne_na.addWidget(_etiq("Nombre d'arcs"))
        self.edits_fraise["n_arcs_taillage"] = QLineEdit("4")
        self.edits_fraise["n_arcs_taillage"].setMaximumWidth(60)
        self.edits_fraise["n_arcs_taillage"].setValidator(QIntValidator(1, 9))
        ligne_na.addWidget(self.edits_fraise["n_arcs_taillage"])
        ligne_na.addStretch(1)
        lay_arcs_t.addLayout(ligne_na)
        b_utiliser = QPushButton("Utiliser dans Simulation / Optimisation")
        b_utiliser.setToolTip("Charge les roues taillées (arcs) comme profil courant, sans passer par un "
                              "fichier ; la roue non taillée garde son profil d'origine.")
        b_utiliser.clicked.connect(self.utiliser_taillage_simulation)
        lay_arcs_t.addWidget(b_utiliser)
        b_exp_t = QPushButton("Exporter (.Eno / .dxf)…")
        b_exp_t.clicked.connect(self.exporter_taillage_arcs)
        lay_arcs_t.addWidget(b_exp_t)
        lay.addWidget(arcs_t)
        lay.addStretch(1)
        return _defilant(contenu)

    def _roue_fraise_changee(self, w):
        """Changement de roue (Menante/Menée) : réaffiche ce qui concerne cette roue (sa fraise,
        et son taillage s'il a déjà été fait -- chaque roue garde le sien)."""
        if self.fraise_paire is None:
            return
        self.edits_fraise["zz_cible"].setText(str(self.fraise_paire[w].zz_origine))
        self._proposer_n_arcs(w)
        t = self.taillage[w]
        self.info_taillage.setText(self._texte_taillage(w, t) if t else "-")
        self._dessiner_fraise()

    def _proposer_n_arcs(self, w):
        """Propose, pour la reconversion en arcs de la roue taillée, le nombre d'arcs de la roue
        d'origine : c'est celui qui permet à `arcs.approximer` de repartir de son découpage, et
        donc de retrouver la forme d'origine au lieu d'une chaîne voisine (voir `approximer`)."""
        if self.fraise_source is None:
            return
        self.edits_fraise["n_arcs_taillage"].setText(str(max(1, len(self.fraise_source[w].points) - 3)))

    def _texte_taillage(self, w, t):
        # Avec une fraise-mère, éloigner la crémaillère revient à changer le module engendré
        # (`module + 2*oy/zz`, le choix O/M de `Fraesen.F90`). Avec un outil de scarving le module
        # est fixé par les deux nombres de dents : `oy` n'est plus qu'un déport de denture, à
        # annoncer comme tel.
        if self.z_outils[w] == 0:
            quoi = (f"module effectif <b>"
                    f"{self.fraise_paire[w].module + 2 * t['oy'] / t['zz']:.6f} mm</b>")
        else:
            quoi = (f"module <b>{self.fraise_paire[w].module:g} mm</b>, déport "
                    f"<b>{t['oy'] * 1e3:.1f} µm</b>")
        texte = (f"Roue taillée ({ROUES[w].lower()}) : Z=<b>{t['zz']}</b>, {quoi} "
                f"({len(t['zp'])} points).<br>")
        if t["ecart"]:
            exc, _, pen, _ = t["ecart"]
            texte += (f"Écart à la référence : excès max <b>{exc * 1e3:.2f} µm</b> (★), "
                     f"pénétration max <b>{abs(pen) * 1e3:.2f} µm</b> (▽).")
        else:
            texte += ("<i>Comparaison non pertinente (nombre de dents ou déport différent de la "
                     "roue d'origine de la fraise).</i>")
        return texte

    def _sigmas(self):
        """Angle de croisement par roue, en radians ; ValueError si incohérent.

        C'est l'angle entre l'axe de l'outil et celui de la pièce. Sans lui il n'y a pas de
        vitesse de coupe : un outil à axes parallèles ne fait que du mortaisage. Il n'a aucun sens
        pour une fraise-mère, où il est simplement ignoré."""
        out = []
        for k, le in enumerate(self.edits_sigma):
            texte = (le.text() or "0").replace(",", ".")
            try:
                v = float(texte)
            except ValueError:
                raise ValueError(f"Σ ({ROUES[k].lower()}) : nombre de degrés attendu")
            if not 0.0 <= v <= 60.0:
                raise ValueError(f"Σ ({ROUES[k].lower()}) : entre 0 et 60 degrés")
            out.append(math.radians(v))
        return out

    def _z_outils(self):
        """Nombre de dents de l'outil par roue (0 = fraise-mère) ; ValueError si incohérent."""
        zs = []
        for k, le in enumerate(self.edits_z_fraise):
            if not le.hasAcceptableInput():
                raise ValueError(f"Z outil ({ROUES[k].lower()}) : entier entre 0 et 9999 attendu")
            z = int(le.text() or 0)
            if z != 0 and z < 3:
                raise ValueError(f"Z outil ({ROUES[k].lower()}) : 0 (fraise-mère) ou au moins 3")
            if self.mobiles is not None and self.mobiles[k].zz < 0 and z == 0:
                raise ValueError(f"{ROUES[k]} est une couronne intérieure : une fraise-mère ne peut "
                                 "pas y entrer, donnez un nombre de dents d'outil de scarving")
            if (self.mobiles is not None and self.mobiles[k].zz < 0
                    and z >= abs(self.mobiles[k].zz)):
                raise ValueError(f"Z outil ({ROUES[k].lower()}) : l'outil doit être plus petit que "
                                 f"la couronne (Z < {abs(self.mobiles[k].zz)}) pour y entrer")
            zs.append(z)
        return zs

    def calculer_fraise(self):
        """Calcule l'outil conjugué de chaque roue de l'onglet Simulation : fraise-mère
        (crémaillère, `Z outil = 0`) ou outil de scarving. Les deux roues sont traitées ensemble,
        comme le Fortran (ENGOPT, option "F"), pour pouvoir enregistrer la paire."""
        if self.mobiles is None:
            QMessageBox.information(self, "Outil", "Ouvrez ou simulez d'abord un profil en arcs "
                                    "de cercle (onglet Simulation).")
            return
        try:
            z_outils = self._z_outils()
            sigmas = self._sigmas()
            v, erreur = self._lire(self.edits_fraise, CHAMPS_FRAISE_CREATION, "fraise")
            if erreur:
                raise ValueError(erreur)
        except ValueError as e:
            self.info_fraise.setText(f"<span style='color:{theme.ERREUR}'>{e}</span>")
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from . import fraise as fraise_mod
            paire = tuple(
                fraise_mod.profil_fraise(m, n_periodes=int(v["n_periodes"]) or None,
                                         n_positions=int(v["n_positions"]), pc=v["pc"],
                                         n_points=int(v["n_points"]))
                if z == 0 else
                fraise_mod.profil_outil(m, z, sg, n_points=int(v["n_points"]))
                for m, z, sg in zip(self.mobiles, z_outils, sigmas))
        except (ValueError, ArithmeticError, NotImplementedError) as e:
            self.info_fraise.setText(f"<span style='color:{theme.ERREUR}'>{type(e).__name__} : {e}</span>")
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.fraise_paire, self.fraise_source = paire, self.mobiles
        self.z_outils, self.sigmas = z_outils, sigmas
        self.taillage = [None, None]
        self.info_taillage.setText("-")
        w = self.combo_fraise_roue.currentIndex()
        self.edits_fraise["zz_cible"].setText(str(self.mobiles[w].zz))
        self._proposer_n_arcs(w)
        detail = []
        for m, z, o in zip(self.mobiles, z_outils, paire):
            if z == 0:
                detail.append(f"<b>{m.nom}</b> (Z={m.zz}) : fraise-mère, module {o.module:g} mm")
            else:
                taux = fraise_mod.derive_affutage(m, o.zz, o.sigma)
                detail.append(
                    f"<b>{m.nom}</b> (Z={m.zz}) : outil de scarving Z={o.zz}, Σ="
                    f"{math.degrees(o.sigma):g}°, rayon primitif {o.rayon:.4f} mm, module apparent "
                    f"{o.module_apparent:.6f} mm, hélice {math.degrees(o.sigma):g}°"
                    + (" — <i>dent pointue</i>." if o.tronque else
                       " — pas de tête imposée par la pièce.")
                    + ("" if not o.sigma else
                       f" Dérive au réaffûtage : <b>{taux * 1e3 / 0.05:.1f} µm/mm</b>.")
                    + (f"<br><span style='color:{theme.ERREUR}'>Σ = 0 : axes parallèles, l'outil ne coupe "
                       "pas (mortaisage).</span>" if not o.sigma else ""))
        self.info_fraise.setText("Outils calculés :<br>" + "<br>".join(detail)
                                 + f"<br>{len(paire[0].points)} points chacun.")
        self._dessiner_fraise()
        self._aller_onglet(3)

    def sauver_fraise(self):
        if self.fraise_paire is None:
            QMessageBox.information(self, "Fraise", "Calculez d'abord un outil.")
            return
        if any(z != 0 for z in self.z_outils):
            QMessageBox.information(self, "Fraise", MSG_SCARVING)
            return
        ref = self.chemin or self.source
        proposition = ref.with_name(f"{ref.stem}.Frs") if ref else RACINE / "Files" / "fraise.Frs"
        nom, _ = QFileDialog.getSaveFileName(self, "Enregistrer la fraise", str(proposition), "Fraise (*.Frs)")
        if not nom:
            return
        chemin = Path(nom)
        if chemin.suffix.lower() != ".frs":
            chemin = chemin.with_suffix(".Frs")
        from . import fraise as fraise_mod
        hertz = None
        if self.fraise_source is not None and all(m.e_mod for m in self.fraise_source):
            hertz = [(m.e_mod, m.poisson, m.epaisseur) for m in self.fraise_source]
        try:
            fraise_mod.ecrire_frs(chemin, self.fraise_paire, hertz)
        except OSError as e:
            QMessageBox.warning(self, "Fraise", f"Écriture impossible : {e}")
            return
        self._message(f"Fraises enregistrées (menante + menée) : {chemin}.")

    def ouvrir_fraise(self):
        ref = self.chemin or self.source
        dossier = ref.parent if ref else RACINE / "Files"
        nom, _ = QFileDialog.getOpenFileName(self, "Ouvrir une fraise", str(dossier), "Fraise (*.Frs)")
        if not nom:
            return
        from . import fraise as fraise_mod
        try:
            f_m, f_e, hertz = fraise_mod.lire_frs(nom)
        except (OSError, ValueError, IndexError) as e:
            QMessageBox.critical(self, "Fraise", f"Lecture impossible de {nom} :\n{e}")
            return
        self.fraise_paire, self.fraise_source = (f_m, f_e), None  # pas de référence : pas de comparaison
        self.taillage = [None, None]
        w = self.combo_fraise_roue.currentIndex()
        self.edits_fraise["zz_cible"].setText(str(self.fraise_paire[w].zz_origine))
        self.info_fraise.setText(
            f"Fraises chargées : <b>{f_m.nom}</b> (module {f_m.module:g} mm, {len(f_m.points)} points, "
            f"issue d'une roue à {f_m.zz_origine} dents) et <b>{f_e.nom}</b> (module {f_e.module:g} mm, "
            f"{len(f_e.points)} points, issue d'une roue à {f_e.zz_origine} dents)."
            + (f"<br>Données de Hertz du fichier : {hertz}" if hertz else ""))
        self._dessiner_fraise()
        self._aller_onglet(3)

    def _tailler_roue(self, w, zz, oy, n_periodes, n_positions, n_points):
        """Taille la roue `w`, remplit `self.taillage[w]` ; retourne une erreur (str) ou None."""
        from . import fraise as fraise_mod
        f = self.fraise_paire[w]
        # fraise-mère : roulement contre une droite (`tailler`) ; outil de scarving : roulement
        # contre une roue (`tailler_outil`). `oy` a le même sens dans les deux cas, un déport de
        # denture, et le même signe (outil éloigné de la pièce).
        tailleur = fraise_mod.tailler if self.z_outils[w] == 0 else fraise_mod.tailler_outil
        try:
            zp = tailleur(f, zz, oy=oy, n_periodes=n_periodes, n_positions=n_positions,
                          n_points=n_points)
        except (ValueError, ArithmeticError) as e:
            return f"{type(e).__name__} : {e}"
        ecart = None
        if self.fraise_source is not None and zz == self.fraise_source[w].zz and abs(oy) < 1e-9:
            ecart = fraise_mod.ecart(zp, self.fraise_source[w].echantillonner(2000))
        self.taillage[w] = {"zp": zp, "ecart": ecart, "zz": zz, "oy": oy}
        return None

    def tailler_mobile(self):
        """Taille la roue choisie avec Z cible / déport (cas avancé, un déport de denture ou un
        nombre de dents différent de l'original)."""
        if self.fraise_paire is None:
            QMessageBox.information(self, "Fraise", "Calculez ou ouvrez d'abord une fraise.")
            return
        v, erreur = self._lire(self.edits_fraise, CHAMPS_FRAISE_TAILLAGE, "taillage")
        if erreur:
            self.info_taillage.setText(f"<span style='color:{theme.ERREUR}'>{erreur}</span>")
            return
        w = self.combo_fraise_roue.currentIndex()
        f = self.fraise_paire[w]
        zz = int(v["zz_cible"])
        oy = v["oy_um"] / 1000.0
        if v["module_cible"] and abs(v["module_cible"] - f.module) > 1e-9:
            oy = (v["module_cible"] - f.module) * zz / 2.0
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            erreur = self._tailler_roue(w, zz, oy, int(v["n_periodes_t"]) or None,
                                        int(v["n_positions_t"]), int(v["n_points_t"]))
        finally:
            QApplication.restoreOverrideCursor()
        if erreur:
            self.info_taillage.setText(f"<span style='color:{theme.ERREUR}'>{erreur}</span>")
            return
        self.info_taillage.setText(self._texte_taillage(w, self.taillage[w]))
        self._dessiner_fraise()

    def tailler_les_deux(self):
        """Taille les deux roues avec leur propre nombre de dents et sans déport (le cas
        courant : vérifier comment un couple déjà optimisé survit au taillage) -- ignore les
        champs Z cible / déport, garde la résolution du taillage déjà saisie."""
        if self.fraise_paire is None:
            QMessageBox.information(self, "Fraise", "Calculez ou ouvrez d'abord une fraise.")
            return
        v, erreur = self._lire(self.edits_fraise, CHAMPS_FRAISE_TAILLAGE[3:], "taillage")
        if erreur:
            self.info_taillage.setText(f"<span style='color:{theme.ERREUR}'>{erreur}</span>")
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            erreurs = [self._tailler_roue(w, self.fraise_paire[w].zz_origine, 0.0,
                                          int(v["n_periodes_t"]) or None, int(v["n_positions_t"]),
                                          int(v["n_points_t"])) for w in (0, 1)]
        finally:
            QApplication.restoreOverrideCursor()
        w = self.combo_fraise_roue.currentIndex()
        if any(erreurs):
            self.info_taillage.setText(f"<span style='color:{theme.ERREUR}'>{' | '.join(e for e in erreurs if e)}</span>")
        else:
            self.info_taillage.setText("<b>Les deux roues sont taillées.</b><br>" + self._texte_taillage(w, self.taillage[w]))
        self._dessiner_fraise()

    def utiliser_taillage_simulation(self):
        """Charge les roues taillées (reconverties en arcs) comme profil courant, sans fichier
        intermédiaire ; la roue qui n'a pas été taillée garde son profil d'origine."""
        if self.fraise_paire is None or not any(self.taillage):
            QMessageBox.information(self, "Fraise", "Taillez d'abord au moins une roue.")
            return
        try:
            n_arcs = max(1, min(9, int(self.edits_fraise["n_arcs_taillage"].text() or 4)))
        except ValueError:
            n_arcs = 4
        try:
            mobiles, ecarts = [], []
            for w in (0, 1):
                m, profil_arcs = self._taillage_vers_mobile(w, n_arcs)
                if m is None:
                    m = self.fraise_source[w] if self.fraise_source is not None else self.mobiles[w]
                mobiles.append(m)
                ecarts.append(profil_arcs)
        except (ValueError, ArithmeticError) as e:
            QMessageBox.warning(self, "Fraise", f"Conversion en arcs impossible : {e}")
            return
        self.mobiles, self.mobiles_points = tuple(mobiles), None
        etats = self._charger_taillage()
        details = "; ".join(f"{ROUES[w].lower()} {_texte_ecart_arcs(p)}"
                            for w, p in enumerate(ecarts) if p is not None)
        # Rien ne garantit qu'un profil taillé soit approximable à un nombre d'arcs donné : la
        # chaîne peut revenir sur elle-même et être refusée par SWTE. Sans ce contrôle, l'échec
        # ne se voyait qu'en aval, par un refus de la simulation (mesuré sur une couronne taillée
        # au scarving, invalide à 5 arcs mais correcte à 4 et à 6).
        from . import optim as optim_mod
        mauvaises = [ROUES[w].lower() for w, (m, e) in enumerate(zip(mobiles, ecarts))
                     if e is not None and not optim_mod.profil_valide(m)]
        self._message(f"Profil courant = roues taillées en arcs ({', '.join(etats)})." +
                      (f" Écart d'approximation en arcs : {details}." if details else "") +
                      (f" ATTENTION : l'approximation à {n_arcs} arcs de la {' et de la '.join(mauvaises)}"
                       " se replie sur elle-même (contrôle SWTE) ; essayez un autre nombre d'arcs,"
                       " ou simulez la suite de points." if mauvaises else ""),
                      erreur=bool(mauvaises))

    def utiliser_taillage_points(self):
        """Charge les roues taillées telles quelles, en suite de points : c'est la forme que la
        fraise engendre réellement, donc celle qu'il faut simuler. L'onglet Optimisation et
        l'export .Eno demandent des arcs et renverront vers le groupe de reconversion."""
        if self.fraise_paire is None or not any(self.taillage):
            QMessageBox.information(self, "Fraise", "Taillez d'abord au moins une roue.")
            return
        profils = []
        for w in (0, 1):
            t = self.taillage[w]
            precedent = self._profil_precedent(w)
            if t is None:  # roue non taillée : elle garde son profil, échantillonné
                profils.append(ProfilPoints(
                    precedent.nom, precedent.mo, precedent.zz, precedent.echantillonner(NI).tolist(),
                    precedent.e_mod, precedent.poisson, precedent.epaisseur))
                continue
            module_eff = self.fraise_paire[w].module + 2.0 * t["oy"] / t["zz"]
            profils.append(ProfilPoints(
                (precedent.nom if precedent is not None else self.fraise_paire[w].nom) + "-taillee",
                module_eff, t["zz"], t["zp"].tolist(),
                *((precedent.e_mod, precedent.poisson, precedent.epaisseur)
                  if precedent is not None else (None, None, None))))
        self.mobiles, self.mobiles_points = None, tuple(profils)
        etats = self._charger_taillage()
        self._message(f"Profil courant = roues taillées en suite de points ({', '.join(etats)}) : "
                      f"{len(profils[0].points)} et {len(profils[1].points)} points, sans approximation.")

    def _profil_precedent(self, w):
        """Mobile en arcs d'où vient la fraise de la roue `w` (ou, à défaut, le profil courant)."""
        if self.fraise_source is not None:
            return self.fraise_source[w]
        return self.mobiles[w] if self.mobiles is not None else None

    def _charger_taillage(self):
        """Passe à l'onglet Simulation et simule le profil qui vient d'être chargé ; retourne
        l'état (taillée / d'origine) de chaque roue."""
        etats = ["taillée" if self.taillage[w] else "d'origine (non taillée)" for w in (0, 1)]
        self.origine = f"roues taillées par génération (menante {etats[0]}, menée {etats[1]})"
        self._aller_onglet(1)
        self.simuler()
        return etats

    def exporter_taillage_points_dxf(self):
        """Exporte le contour complet des roues taillées en polyligne (elles ne sont pas des
        chaînes d'arcs : c'est l'enveloppe du balayage de la fraise)."""
        if self.fraise_paire is None or not any(self.taillage):
            QMessageBox.information(self, "DXF", "Taillez d'abord au moins une roue.")
            return
        ref = self.chemin or self.source
        proposition = ref.with_name(f"{ref.stem}_taillee.dxf") if ref else RACINE / "Files" / "taillee.dxf"
        nom, _ = QFileDialog.getSaveFileName(self, "Exporter les roues taillées en DXF",
                                             str(proposition), "DXF (*.dxf)")
        if not nom:
            return
        base = Path(nom)
        if base.suffix.lower() == ".dxf":
            base = base.with_suffix("")
        from . import dxf as dxf_mod
        ecrits = []
        try:
            for w in (0, 1):
                t = self.taillage[w]
                if t is None:
                    continue
                chemin = base.with_name(f"{base.name}_{('menante', 'menee')[w]}.dxf")
                dxf_mod.ecrire_polyligne(chemin, t["zp"], t["zz"])
                ecrits.append(chemin.name)
        except OSError as e:
            QMessageBox.warning(self, "DXF", f"Écriture impossible : {e}")
            return
        self._message(f"Roue(s) taillée(s) exportée(s) en DXF (polyligne) : {', '.join(ecrits)}.")

    def _dessiner_fraise(self):
        self.ax.clear()
        w = self.combo_fraise_roue.currentIndex()
        t = self.taillage[w] if self.fraise_paire is not None else None
        if t is not None:
            zp = t["zp"].tolist()
            gauche = [(-x, y) for x, y in reversed(zp)]
            xs, ys = zip(*(zp + gauche[1:]))
            self.ax.plot(xs, ys, color="tab:blue", lw=1.4, label="profil taillé")
            if t["ecart"]:
                exc, p_exc, pen, p_pen = t["ecart"]
                # la référence est une chaîne d'arcs : on la dessine échantillonnée sur ses arcs,
                # pas en reliant ses points de transition par des droites (c'est aussi sur les
                # arcs que `fraise.ecart` mesure l'écart affiché)
                ref = self.fraise_source[w].echantillonner(600).tolist()
                gauche_r = [(-x, y) for x, y in reversed(ref)]
                xr, yr = zip(*(ref + gauche_r[1:]))
                self.ax.plot(xr, yr, "--", color="tab:red", lw=1.0, label="profil de référence")
                self.ax.plot(*p_exc, "*", color="orange", ms=14, mec=theme.encre(), label=f"excès max {exc * 1e3:.2f} µm")
                self.ax.plot(*p_pen, "v", color="purple", ms=10, mec=theme.encre(),
                            label=f"pénétration max {abs(pen) * 1e3:.2f} µm")
            self.ax.set_title(f"Roue taillée par génération ({ROUES[w].lower()})")
        elif self.fraise_paire is not None:
            f = self.fraise_paire[w]
            pts = np.asarray(f.points, dtype=float)
            if self.z_outils[w] == 0:  # fraise-mère : profil rectiligne, répété le long du pas
                pas = f.pc * math.pi * f.module
                dent = pts.tolist() + [(-x, y) for x, y in reversed(pts.tolist())][1:]
                for k in (-1, 0, 1):
                    self.ax.plot([x + k * pas for x, _ in dent], [y for _, y in dent],
                                 color="tab:green", lw=1.4, label="fraise-mère" if k == 0 else None)
                self.ax.axhline(0.0, color=theme.encre(), lw=0.8, ls="--", label="droite primitive")
                self.ax.set_title(f"Fraise-mère ({f.nom}, module {f.module:g} mm)")
            else:  # outil de scarving : une roue, on dessine quelques flancs conjugués
                pas = 2.0 * math.pi / f.zz
                for k in (-2, -1, 0, 1, 2):
                    for s in (1.0, -1.0):
                        q = np.column_stack((s * pts[:, 0], pts[:, 1]))
                        c, si = math.cos(k * pas), math.sin(k * pas)
                        self.ax.plot(q[:, 0] * c - q[:, 1] * si, q[:, 0] * si + q[:, 1] * c,
                                     color="tab:green", lw=1.4,
                                     label="flanc actif de l'outil" if (k == 0 and s > 0) else None)
                th = np.linspace(0.0, 2.0 * math.pi, 400)
                self.ax.plot(f.rayon * np.cos(th), f.rayon * np.sin(th), "--", color=theme.encre(), lw=0.8,
                             label="primitif de l'outil")
                self.ax.set_title(f"Outil de scarving ({f.nom}, Z={f.zz}, "
                                  f"rayon primitif {f.rayon:g} mm)")
        else:
            self.ax.set_title("Aucune fraise calculée")
        self.ax.set_aspect("equal")
        self.ax.grid(True, alpha=0.3)
        if self.ax.get_legend_handles_labels()[0]:
            self.ax.legend(loc="best", fontsize="small")
        self.canvas.draw_idle()

    def exporter_fraise_dxf(self):
        """Exporte la fraise choisie (5 dents) en DXF (LWPOLYLINE, ce n'est pas une chaîne d'arcs)."""
        if self.fraise_paire is None:
            QMessageBox.information(self, "Fraise", "Calculez ou ouvrez d'abord une fraise.")
            return
        w = self.combo_fraise_roue.currentIndex()
        f = self.fraise_paire[w]
        ref = self.chemin or self.source
        proposition = ref.with_name(f"{ref.stem}_fraise.dxf") if ref else RACINE / "Files" / "fraise.dxf"
        nom, _ = QFileDialog.getSaveFileName(self, "Exporter la fraise en DXF", str(proposition), "DXF (*.dxf)")
        if not nom:
            return
        chemin = Path(nom)
        if chemin.suffix.lower() != ".dxf":
            chemin = chemin.with_suffix(".dxf")
        from . import dxf as dxf_mod
        try:
            dxf_mod.ecrire_courbe(chemin, f.echantillonner(400), n_periodes=5, pas=f.pc * math.pi * f.module)
        except OSError as e:
            QMessageBox.warning(self, "Fraise", f"Écriture impossible : {e}")
            return
        self._message(f"Fraise ({ROUES[w].lower()}) exportée en DXF : {chemin}.")

    def _taillage_vers_mobile(self, w, n_arcs):
        """Mobile en arcs (module effectif, données de Hertz reprises de la roue d'origine si
        connue) à partir du taillage de la roue `w`, et son écart d'approximation ; (None, None)
        si cette roue n'a pas encore été taillée.

        Quand la roue d'origine est connue (fraise calculée ici, pas chargée d'un .Frs), ses
        points de transition sont passés à `approximer` : le congé de pied, que la fraise creuse
        en trochoïde et qui n'est pas une surface de travail, ne fausse alors plus le choix du
        découpage (voir `arcs.approximer`)."""
        t = self.taillage[w]
        if t is None:
            return None, None
        from . import arcs as arcs_mod
        ref = self.fraise_source[w] if self.fraise_source is not None else None
        module_eff = self.fraise_paire[w].module + 2.0 * t["oy"] / t["zz"]
        reference = ref.points if ref is not None and ref.zz == t["zz"] else None
        profil = arcs_mod.approximer(t["zp"].tolist(), t["zz"], 0.0, n_arcs, reference=reference)
        n = profil.n_points
        drapeaux = [0 if i in (0, n - 1) else 1 for i in range(n)]
        hertz = (ref.e_mod, ref.poisson, ref.epaisseur) if ref else (None, None, None)
        nom = (ref.nom if ref else self.fraise_paire[w].nom) + "-taillee"
        mobile = arcs_mod.Mobile(nom, module_eff, t["zz"], profil.points, drapeaux, *hertz)
        return mobile, profil

    def exporter_taillage_arcs(self):
        """Reconvertit la roue taillée en arcs de cercle et l'exporte en .Eno (les deux roues
        taillées nécessaires) ou en DXF (une seule suffit)."""
        if self.fraise_paire is None:
            QMessageBox.information(self, "Fraise", "Taillez d'abord un mobile (bouton « Tailler »).")
            return
        w = self.combo_fraise_roue.currentIndex()
        try:
            n_arcs = max(1, min(9, int(self.edits_fraise["n_arcs_taillage"].text() or 4)))
        except ValueError:
            n_arcs = 4
        try:
            mobile_w, ecart_w = self._taillage_vers_mobile(w, n_arcs)
            mobile_autre, _ = self._taillage_vers_mobile(1 - w, n_arcs)
        except (ValueError, ArithmeticError) as e:
            QMessageBox.warning(self, "Fraise", f"Conversion en arcs impossible : {e}")
            return
        if mobile_w is None:
            QMessageBox.information(self, "Fraise", "Taillez d'abord cette roue (bouton « Tailler »).")
            return
        paire_ok = mobile_autre is not None
        ref = self.chemin or self.source
        proposition = (ref.with_name(f"{ref.stem}_taillee") if ref else RACINE / "Files" / "taillee")
        filtre = "Fichiers Eno (*.Eno);;DXF (*.dxf)" if paire_ok else "DXF (*.dxf)"
        nom, filtre_choisi = QFileDialog.getSaveFileName(self, "Exporter la roue taillée (arcs)",
                                                         str(proposition), filtre)
        if not nom:
            return
        chemin = Path(nom)
        if not chemin.suffix:
            chemin = chemin.with_suffix(".Eno" if "Eno" in filtre_choisi else ".dxf")
        from . import arcs as arcs_mod
        from . import dxf as dxf_mod
        try:
            if chemin.suffix.lower() == ".eno":
                if not paire_ok:
                    QMessageBox.warning(self, "Fraise", "Le format .Eno demande les deux roues taillées : "
                                        "taillez aussi l'autre roue avant d'exporter en .Eno.")
                    return
                paire = (mobile_w, mobile_autre) if w == 0 else (mobile_autre, mobile_w)
                arcs_mod.ecrire_eno(chemin, paire)
            else:
                dxf_mod.ecrire_profil(chemin, mobile_w)
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, "Fraise", f"Écriture impossible : {e}")
            return
        self._message(f"Roue taillée ({ROUES[w].lower()}) reconvertie en {n_arcs} arcs et exportée : "
                      f"{chemin} (écart d'approximation {_texte_ecart_arcs(ecart_w)}).")

    def _marquer_modifie(self):
        if not self.modifie:
            self.modifie = True
            self._titre()

    def _creer_graphiques(self):
        # marges fixes : la mise en page automatique est trop lente pour suivre le curseur
        self.figure = Figure(figsize=(8, 6))
        self.figure.subplots_adjust(left=0.06, right=0.98, bottom=0.06, top=0.91)
        self.ax = self.figure.add_subplot()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.barre = NavigationToolbar2QT(self.canvas, self)
        navigation.activer(self.canvas)  # molette, glisser, double-clic : comme la vue de la came
        cadre = QWidget()
        lay = QVBoxLayout(cadre)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.barre)
        lay.addWidget(self.canvas, 1)
        return cadre

    # --- panneau <-> modèle ---------------------------------------------
    def remplir_panneau(self):
        for k, roue in enumerate(self.projet.roues):
            for c in CHAMPS + CHAMPS_EXPORT:
                v = getattr(roue, c.cle)
                if v is None:  # donnée de Hertz absente du fichier
                    v = HERTZ_DEFAUT[k][("e_mod", "poisson", "epaisseur").index(c.cle)]
                le = self.edits[k][c.cle]
                le.setText(str(int(v)) if c.entier else texte_nombre(v))
                le.setStyleSheet("")
        for c in CHAMPS_SIM:
            le = self.edits_sim[c.cle]
            le.setText(texte_nombre(getattr(self.projet, c.cle)))
            le.setStyleSheet("")
        self.remplir_optimisation()

    def _lire(self, edits, champs, nom_de):
        """Valeurs des champs ; ValueError (après avoir marqué les champs) si l'un est invalide."""
        valeurs, erreur = {}, None
        for c in champs:
            le = edits[c.cle]
            if not le.hasAcceptableInput():
                le.setStyleSheet(STYLE_ERREUR)
                erreur = erreur or f"{c.nom} ({nom_de}) : valeur entre {c.mini:g} et {c.maxi:g} attendue"
                continue
            v = int(le.text()) if c.entier else float(le.text())
            if c.exclus is not None and c.exclus[0] <= v <= c.exclus[1]:
                le.setStyleSheet(STYLE_ERREUR)
                erreur = erreur or (f"{c.nom} ({nom_de}) : valeur entre {c.mini:g} et {c.maxi:g} "
                                    f"attendue, hors de [{c.exclus[0]:g}, {c.exclus[1]:g}]")
                continue
            le.setStyleSheet("")
            valeurs[c.cle] = v
        return valeurs, erreur

    def lire_roues(self):
        """Roues d'après l'onglet Création ; ValueError si une valeur est invalide."""
        roues, erreur = [], None
        for k in range(2):
            valeurs, e = self._lire(self.edits[k], CHAMPS + CHAMPS_EXPORT, ROUES[k].lower())
            erreur = erreur or e
            if not e:
                roues.append(replace(self.projet.roues[k], **valeurs))  # garde les données Hertz
        if erreur:
            raise ValueError(erreur)
        return roues

    def lire_simulation(self):
        """Variables de simulation (onglet Simulation) ; ValueError si l'une est invalide."""
        valeurs, erreur = self._lire(self.edits_sim, CHAMPS_SIM, "simulation")
        if erreur:
            raise ValueError(erreur)
        return valeurs

    def lire_optimisation(self):
        """Paramètres de l'onglet Optimisation ; ValueError si une valeur est invalide."""
        opti, erreur = self._lire(self.edits_opti, CHAMPS_OPTI, "optimisation")
        if erreur:
            raise ValueError(erreur)
        opti["principal"] = self.combo_principal.currentData()
        opti["secondaire"] = self.combo_secondaire.currentData()
        opti["dissymetrie"] = int(self.case_dissym.isChecked())
        opti["equilibre"] = int(self.case_equilibre.isChecked())
        opti["libres"] = {}
        for role, liste in zip("me", self.liste_points):
            coches = [liste.item(i).data(Qt.UserRole) for i in range(liste.count())
                      if liste.item(i).checkState() == Qt.Checked]
            opti["libres"][role] = None if len(coches) == liste.count() else coches
        return opti

    def remplir_optimisation(self):
        opti = self.projet.optimisation
        for combo, cle in ((self.combo_principal, "principal"), (self.combo_secondaire, "secondaire")):
            combo.blockSignals(True)
            combo.setCurrentIndex(max(0, combo.findData(opti[cle])))
            combo.blockSignals(False)
        # `blockSignals` obligatoire : sans lui `stateChanged` marque le projet modifié à chaque
        # remplissage programmé (ouverture, rechargement), et l'application réclame une sauvegarde
        # alors que l'utilisateur n'a rien touché -- le bug déjà rencontré sur les deux combos.
        for case, cle, defaut in ((self.case_dissym, "dissymetrie", 1),
                                  (self.case_equilibre, "equilibre", 0)):
            case.blockSignals(True)
            case.setChecked(bool(opti.get(cle, defaut)))
            case.blockSignals(False)
        for c in CHAMPS_OPTI:
            le = self.edits_opti[c.cle]
            le.setText(str(int(opti[c.cle])) if c.entier else texte_nombre(opti[c.cle]))
            le.setStyleSheet("")

    def _remplir_variables(self):
        """Listes des arcs libres, d'après les arcs simulés (le dernier arc est refermé, pas libre)."""
        libres = self.projet.optimisation["libres"]
        for liste in self.liste_points:  # vidées si le profil courant n'est pas en arcs
            liste.clear()
        for role, liste, mobile in zip("me", self.liste_points, self.mobiles or ()):
            liste.blockSignals(True)
            # convexité physique : elle est inversée sur une couronne intérieure
            _, rayons = mobile.arcs()
            rayons = arcs.convexite(mobile.zz, rayons)
            for i in range(len(rayons) - 1):
                r = rayons[i]
                item = QListWidgetItem(f"Arc {i + 1} : rayon {r:+.4f} mm ({'convexe' if r > 0 else 'concave'})")
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setData(Qt.UserRole, i)
                cochee = libres[role] is None or i in libres[role]
                item.setCheckState(Qt.Checked if cochee else Qt.Unchecked)
                liste.addItem(item)
            liste.blockSignals(False)

    def _exclure_autres_listes(self, active):
        """Sélection exclusive entre les listes d'arcs de la menante et de la menée (Simulation)."""
        for liste in self.liste_arcs_sim:
            if liste is not active and liste.currentRow() >= 0:
                liste.blockSignals(True)
                liste.clearSelection()
                liste.blockSignals(False)

    def _remplir_arcs_sim(self):
        """Listes des arcs du profil actuel (onglet Simulation, tous les arcs y compris le dernier,
        contrairement aux arcs libres de l'onglet Optimisation)."""
        # sans arcs (profil en suite de points), il n'y a rien à modifier : listes vides et
        # groupe désactivé plutôt qu'un groupe actif qui ne répondrait pas
        self.groupe_modif.setEnabled(self.mobiles is not None)
        for liste in self.liste_arcs_sim:
            liste.clear()
        for liste, mobile in zip(self.liste_arcs_sim, self.mobiles or ()):
            liste.blockSignals(True)
            # convexité physique : elle est inversée sur une couronne intérieure
            _, rayons = mobile.arcs()
            rayons = arcs.convexite(mobile.zz, rayons)
            for i, r in enumerate(rayons):
                liste.addItem(f"Arc {i + 1} : rayon {r:+.4f} mm ({'convexe' if r > 0 else 'concave'})")
            liste.blockSignals(False)

    def _arc_selectionne(self):
        """(roue 0 ou 1, indice de l'arc) sélectionné dans l'onglet Simulation, ou None."""
        for w, liste in enumerate(self.liste_arcs_sim):
            if liste.currentRow() >= 0:
                return w, liste.currentRow()
        return None

    def ajouter_point_arc(self):
        """Scinde l'arc sélectionné en deux, sur le même cercle : la forme du profil ne change
        pas, c'est donc toujours possible sans risque."""
        if self.mobiles is None:
            return
        choix = self._arc_selectionne()
        if choix is None:
            QMessageBox.information(self, "Modifier les arcs", "Sélectionnez d'abord un arc.")
            return
        w, i = choix
        try:
            nouveau = arcs.ajouter_point(self.mobiles[w], i)
        except ValueError as e:
            QMessageBox.warning(self, "Modifier les arcs", str(e))
            return
        mobiles = list(self.mobiles)
        mobiles[w] = nouveau
        self.mobiles = tuple(mobiles)
        self.simuler()
        self._message(f"Point ajouté sur l'arc {i + 1} de la {ROUES[w].lower()} : "
                      f"{len(nouveau.points) - 1} arcs maintenant, forme inchangée.")

    def retirer_point_arc(self):
        """Fusionne l'arc sélectionné avec son voisin (change la forme après ce point) ; demande
        confirmation en montrant l'écart et la validité du nouveau profil."""
        if self.mobiles is None:
            return
        choix = self._arc_selectionne()
        if choix is None:
            QMessageBox.information(self, "Modifier les arcs", "Sélectionnez d'abord un arc.")
            return
        w, i = choix
        try:
            nouveau = arcs.retirer_point(self.mobiles[w], i)
        except ValueError as e:
            QMessageBox.warning(self, "Modifier les arcs", str(e))
            return
        from . import optim
        ecart = arcs.ecart_um(self.mobiles[w], nouveau)
        valide = optim.profil_valide(nouveau)
        reponse = QMessageBox.question(
            self, "Supprimer ce point",
            f"La forme du profil changera après ce point : écart maximal ≈ {ecart:.2f} µm.\n"
            + ("Le nouveau profil reste valable." if valide else
               "Le nouveau profil n'est PAS valable pour l'optimisation (rayons de courbure) : "
               "il pourra quand même servir de point de départ, la recherche en trouvera un autre à proximité.")
            + "\n\nConfirmer la suppression ?",
            QMessageBox.Yes | QMessageBox.No)
        if reponse != QMessageBox.Yes:
            return
        mobiles = list(self.mobiles)
        mobiles[w] = nouveau
        self.mobiles = tuple(mobiles)
        self.simuler()
        self._message(f"Point supprimé sur la {ROUES[w].lower()} : {len(nouveau.points) - 1} arcs "
                      f"maintenant, écart ≈ {ecart:.2f} µm.")

    def _extras_courants(self):
        """Paramètres du fichier d'arcs : plage d'entraxe (onglet Simulation) et critères
        d'optimisation (onglet Optimisation), avec les noms de l'ancien programme."""
        extras = json.loads(json.dumps(self.extras or arcs.EXTRAS_DEFAUT))
        try:
            v = self.lire_simulation()
            extras["dA_neg"], extras["dA_pos"] = v["da_neg"], v["da_pos"]
        except ValueError:
            pass
        try:
            opti = self.lire_optimisation()
            noms = {cle: nom for cle, _, nom in OBJECTIFS}
            extras["simulation"].update(opti_1=noms[opti["principal"]], ponderation=opti["ponderation"],
                                        opti_2=noms.get(opti["secondaire"], "aucun"))
        except ValueError:
            pass
        return extras

    def _preparer_optimisation(self):
        """(roues, mobiles, mu, da_neg, da_pos, objectif, options, libres) d'après l'interface, ou
        None (message affiché) si les paramètres ou les arcs sont inutilisables."""
        if self.mobiles_points is not None:
            self.info_opti.setText(f"<span style='color:{theme.ERREUR}'>{MSG_ARCS_REQUIS}</span>")
            return None
        if self.mobiles is None and not self._arcs_de_creation():
            return None
        try:
            v, opti = self.lire_simulation(), self.lire_optimisation()
            from . import optim
            from . import recherche
        except Exception as e:  # noqa: BLE001 - toute erreur doit s'afficher, jamais rester muette
            self.info_opti.setText(f"<span style='color:{theme.ERREUR}'>{type(e).__name__} : {e}</span>")
            return None
        roues = [evolvente.Roue(mob.mo, mob.zz, None, None, 0.0, None, None) for mob in self.mobiles]
        libres = [range(len(mob.points) - 2) if opti["libres"][role] is None else
                  [i for i in opti["libres"][role] if 0 <= i < len(mob.points) - 2]
                  for role, mob in zip("me", self.mobiles)]
        if not sum(len(list(li)) for li in libres):
            self.info_opti.setText(f"<span style='color:{theme.ERREUR}'>Aucun arc libre.</span>")
            return None
        objectif = optim.Objectif(opti["principal"], opti["secondaire"] or None, opti["ponderation"],
                                  equilibre=bool(opti.get("equilibre", 0)),
                                  jeu_min=opti["jeu_min_admis"], jeu_max=opti["jeu_max_admis"] or "auto",
                                  r_tete_min=opti["r_tete_min"])
        options = recherche.Options(departs=opti["departs"], population=opti["population"],
                                    evaluations=opti["evaluations"], etendue_um=opti["deplacement_um"],
                                    graine=opti["graine"], dissymetrie=bool(opti["dissymetrie"]))
        return roues, list(self.mobiles), v["mu"], v["da_neg"], v["da_pos"], objectif, options, libres

    @staticmethod
    def _jeu_max_auto(objectif, roues, mobiles, mu, da_neg, da_pos, libres):
        """Objectif dont le jeu maximal est celui du profil de départ à l'entraxe maximal, si demandé."""
        from . import optim
        if isinstance(objectif.jeu_max, str):
            ref = optim.Probleme(roues, mobiles, mu, da_neg, da_pos, libres, replace(objectif, jeu_max=None)).rapport(None)
            objectif = replace(objectif, jeu_max=ref[2]["jeu_max"] + 0.02)
        return objectif

    def evaluer_criteres(self):
        """Critères du profil courant à -ΔA, 0 et +ΔA, contraintes et durée estimée d'une optimisation."""
        prep = self._preparer_optimisation()
        if prep is None:
            return
        roues, mobiles, mu, da_neg, da_pos, objectif, options, libres = prep
        import numpy as np
        from . import optim
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            objectif = self._jeu_max_auto(objectif, roues, mobiles, mu, da_neg, da_pos, libres)
            pb = optim.Probleme(roues, mobiles, mu, da_neg, da_pos, libres, objectif)
            rapport = pb.rapport(None)
            x = np.random.default_rng(0).normal(0, 0.2, (256, pb.n_var))
            pb(x[:16])
            t0 = time.perf_counter()
            pb(x)
            cadence = len(x) / (time.perf_counter() - t0)
            s, _ = pb.evaluer(np.zeros(pb.n_var))
            viol = {k: float(a[0]) for k, a in pb.violations(s).items()}
        except Exception as e:  # noqa: BLE001 - toute erreur doit s'afficher, jamais rester muette
            self.info_opti.setText(f"<span style='color:{theme.ERREUR}'>{type(e).__name__} : {e}</span>")
            return
        finally:
            QApplication.restoreOverrideCursor()
        lignes = "".join(
            f"<tr><td>ΔA {r['da_um']:+.0f} µm</td><td align=right>{r['couple_var']:.2f} %</td>"
            f"<td align=right>{r['hertz']:.0f}</td><td align=right>{r['rendement']:.2f} %</td>"
            + (f"<td align=right>{r['jeu_min']:.2f}–{r['jeu_max']:.2f}°</td>" if r["jeu_min"] == r["jeu_min"] else "<td></td>")
            + f"<td align=right>{r['interference_um']:.1f} µm</td>"
            f"<td align=right>{r['menee_avant']:.3f} / {r['menee_apres']:.3f}</td>"
            f"<td align=right>{100 * r['symetrie']:.0f} %</td></tr>" for r in rapport)
        libelles = {k: recherche_lib for k, recherche_lib in (
            ("jeu_min", "jeu insuffisant"), ("jeu_max", "jeu excessif"), ("r_tete", "rayon de tête trop petit"),
            ("interference", "interférence"), ("ecart_tete", "têtes différentes"))}
        limites = {"jeu_min": 1e-3, "jeu_max": 1e-3, "r_tete": 1e-3, "interference": optim.EPS_INTERF, "ecart_tete": 1e-4}
        pb_viol = [libelles[k] for k, v in viol.items() if v > limites[k]]
        total = options.departs * options.evaluations * 1.2
        self.info_opti.setText(
            "<table cellspacing=4><tr><th></th><th>Couple</th><th>Hertz</th><th>Rendement</th><th>Jeu</th>"
            "<th>Interférence</th><th>Menée av./ap.</th><th>Équilibre</th></tr>"
            f"{lignes}</table>"
            f"Rayon d'arc de tête minimal : <b>{rapport[0]['r_tete']:.4f} mm</b><br>"
            + (f"Contraintes : <b style='color:{theme.ERREUR}'>" + ", ".join(pb_viol) + "</b><br>" if pb_viol
               else f"Contraintes : <b style='color:{theme.SUCCES}'>respectées</b><br>")
            + ("Équilibre de la menée : <b>pris en compte</b> dans l'objectif<br>"
               if objectif.equilibre else
               "Équilibre de la menée : <i>non pris en compte</i> (case décochée)<br>")
            + f"Variables : {pb.n_var} (courbure et longueur de chaque arc libre)<br>"
            f"Cadence : <b>{cadence:.0f} évaluations/s</b> ({os.cpu_count() or 1} cœurs) ; "
            f"{options.departs} départs de {options.evaluations} évaluations : environ "
            f"<b>{total / cadence / 60:.1f} min</b>")

    # --- optimisation ------------------------------------------------------
    def lancer_optimisation(self):
        if self.travail is not None:
            return
        prep = self._preparer_optimisation()
        if prep is None:
            return
        self.onglets.setCurrentIndex(2)
        self.barre_opti.setRange(0, prep[6].departs * prep[6].evaluations)
        self.barre_opti.setValue(0)
        self.info_opti.setText("Préparation (compilation du noyau au premier lancement)…")
        self._optimisation_en_cours(True)
        self.travail = Travail(prep)
        self.travail.avance.connect(self._avance_optimisation)
        self.travail.fini.connect(self._optimisation_finie)
        self.travail.echec.connect(self._optimisation_echec)
        self.travail.start()

    def arreter_optimisation(self):
        if self.travail is not None:
            self.travail.arreter()
            self.info_opti.setText("Arrêt demandé, fin de la génération en cours…")

    def _optimisation_en_cours(self, oui):
        for b in (self.bouton_evaluer, self.bouton_optimiser):
            b.setEnabled(not oui)
        self.bouton_arret.setEnabled(oui)
        self.action_optimiser.setEnabled(not oui)
        self.action_arret.setEnabled(oui)
        self.barre_opti.setVisible(oui)

    def _avance_optimisation(self, d):
        self.barre_opti.setValue(min(d["evaluations"], self.barre_opti.maximum()))
        self.info_opti.setText(f"{d['phase'].capitalize()} : départ {d['depart']}/{d['departs']}, "
                               f"{d['evaluations']} évaluations, meilleur J = <b>{d['meilleur']:.3f}</b>")

    def _fin_travail(self):
        if self.travail is not None:
            self.travail.wait()
            self.travail = None
        self._optimisation_en_cours(False)

    def _optimisation_echec(self, message):
        self._fin_travail()
        self.info_opti.setText(f"<span style='color:{theme.ERREUR}'>Optimisation impossible : {message}</span>")

    def _optimisation_finie(self, liste, stats):
        self._fin_travail()
        self.optima = liste
        self.table_opti.setRowCount(len(liste))
        for r, o in enumerate(liste):
            rap = o.rapport
            cellules = ["%.3f" % o.j]
            if rap:
                cellules += ["/".join(f"{x['couple_var']:.1f}" for x in rap),
                             f"{max(x['hertz'] for x in rap):.0f}",
                             f"{rap[0]['jeu_min']:.2f}–{rap[2]['jeu_max']:.2f}",
                             # équilibre de la menée : affiché même quand la case est décochée,
                             # c'est une propriété du résultat, pas seulement un terme de l'objectif
                             "/".join(f"{100 * x['symetrie']:.0f}" for x in rap)]
            else:
                cellules += ["-", "-", "-", "-"]
            cellules += ["dissymétriques" if o.dissymetrique else "symétriques",
                         "valable" if o.valable else "non valable : " + ", ".join(
                             _libelle_violation(k) for k in o.violations)]
            for c, texte in enumerate(cellules):
                item = QTableWidgetItem(texte)
                if not o.valable:
                    item.setForeground(Qt.gray)
                self.table_opti.setItem(r, c, item)
        if liste:
            self.table_opti.selectRow(0)
        n_ok = sum(o.valable for o in liste)
        avertissement = ""
        if stats.get("depart_invalide"):
            avertissement = (f"<span style='color:{theme.ERREUR}'>Le profil de départ n'était pas valable "
                             f"({stats['depart_invalide']}) : la recherche est partie d'un profil valide "
                             f"trouvé à proximité à la place.</span><br>")
        self.info_opti.setText(
            avertissement +
            f"<b>{n_ok} optimum(s) valable(s)</b> sur {len(liste)} distincts, "
            f"{stats['evaluations']} évaluations en {stats['temps']:.0f} s "
            f"({stats['evaluations'] / max(stats['temps'], 1e-9):.0f}/s), "
            f"{stats['repares']} réparé(s) par le demi-profil arrière.<br>"
            "Double-clic sur une ligne : affichage dans l'onglet Simulation.")
        self._message(f"Optimisation terminée : {n_ok} solution(s) valable(s).")

    def _optimum_choisi(self):
        r = self.table_opti.currentRow()
        if not 0 <= r < len(self.optima):
            QMessageBox.information(self, "Optimisation", "Choisissez d'abord un optimum dans la liste.")
            return None
        return r, self.optima[r]

    def afficher_optimum(self):
        choix = self._optimum_choisi()
        if choix is None:
            return
        r, o = choix
        try:  # garde les arcs libres cochés (le nombre d'arcs ne change pas)
            self.projet.optimisation["libres"] = self.lire_optimisation()["libres"]
        except ValueError:
            pass
        self.mobiles, self.mobiles_points = tuple(o.mobiles), None
        self.origine = f"optimum {r + 1} de la recherche (J = {o.j:.3f})"
        self._aller_onglet(1)
        self.simuler()

    def exporter_optimum(self):
        choix = self._optimum_choisi()
        if choix is not None:
            self._ecrire_arcs(choix[1].mobiles, f"optimum{choix[0] + 1}")

    def _edite(self, le):
        le.setStyleSheet("")
        if not self.modifie:
            self.modifie = True
            self._titre()

    def _titre(self):
        if self.chemin:
            nom = self.chemin.name
        else:
            nom = f"{self.source.name} (importé)" if self.source else "sans nom"
        self.setWindowTitle(f"Engrenage - {nom}{' *' if self.modifie else ''}")

    def _message(self, texte, erreur=False):
        self.statut.setStyleSheet(f"color: {theme.ERREUR if erreur else theme.SUCCES};")
        self.statut.setText(texte)

    # --- onglets ---------------------------------------------------------
    def _changer_onglet(self, i):
        self._vue_ok = False
        if i == 1 and self.mobiles is None and self.mobiles_points is None:
            self.utiliser_arcs_creation()  # première visite : arcs de la création
        elif i == 2 and self.mobiles is None:
            if self.mobiles_points is not None:
                self.info_opti.setText(f"<span style='color:{theme.ERREUR}'>{MSG_ARCS_REQUIS}</span>")
            elif self._arcs_de_creation():
                self._remplir_variables()
        else:
            self.afficher()
            self.barre.update()

    def actualiser(self):
        """F5 : reconstruit l'onglet Création, relance la simulation, ou retaille (onglet Fraise)."""
        idx = self.onglets.currentIndex()
        if idx == 0:
            self.reconstruire()
        elif idx == 3:
            if self.fraise_paire is not None and self.taillage[self.combo_fraise_roue.currentIndex()]:
                self.tailler_mobile()
            elif self.fraise_paire is not None:
                self._dessiner_fraise()
        else:
            self.simuler()

    def afficher(self):
        idx = self.onglets.currentIndex()
        if idx == 0:
            self._dessiner_creation()
        elif idx == 3:
            self._dessiner_fraise()
        else:
            self.afficher_position(self.curseur.value())

    def _aller_onglet(self, i):
        self.onglets.blockSignals(True)
        self.onglets.setCurrentIndex(i)
        self.onglets.blockSignals(False)
        self._vue_ok = False

    # --- onglet Création -------------------------------------------------
    def _arcs_de(self, roue, profil):
        cle = tuple(getattr(roue, c) for c in ("mo", "zz", "nw", "qe", "qg", "ma", "me", "narcs"))
        if cle not in self._cache_arcs:
            self._cache_arcs[cle] = arcs.approximer(profil, roue.zz, roue.qg, roue.narcs)
        return self._cache_arcs[cle]

    def _calculer_arcs(self, roues, profils):
        """(arcs des deux roues, points échantillonnés pour le dessin, texte d'information)."""
        try:
            arcs_roues = [self._arcs_de(r, p) for r, p in zip(roues, profils)]
        except (ValueError, ArithmeticError, IndexError) as e:
            return None, None, f"<span style='color:{theme.ERREUR}'>Approximation impossible : {e}</span>"
        info = "<br>".join(
            f"{nom} : {a.n_points - 1} arcs, écart max <b>{a.ecart_max * 1000:.2f} µm</b>"
            + (f" <span style='color:{theme.ERREUR}'>(tangente finale à {math.degrees(a.erreur_angle):.2f}°)</span>"
               if abs(a.erreur_angle) > math.radians(0.05) else "")
            for nom, a in zip(ROUES, arcs_roues)
        )
        return arcs_roues, [a.echantillonner().tolist() for a in arcs_roues], info

    def reconstruire(self):
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            roues = self.lire_roues()
            profils = [evolvente.evolve(r, NI) for r in roues]
            arcs_roues, echantillons, info = self._calculer_arcs(roues, profils)
            limites = evolvente.limites_engrenage(roues, profils, N_DENTS, [evolvente.SansContact()])
        except (ValueError, ArithmeticError) as e:
            self._message(f"Construction impossible : {e}", erreur=True)
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.projet.roues = roues
        self.roues, self.profils, self.limites_c = roues, profils, limites
        self.arcs_roues, self.arcs_echantillons = arcs_roues, echantillons
        self.info_arcs.setText(info)
        self._vue_ok = False
        self._dessiner_creation()
        self.barre.update()  # nouvelle vue d'origine pour le zoom
        self._message(f"Engrenage construit ({NI} points par profil).")

    def _dessiner_creation(self):
        if self.roues is None or self.onglets.currentIndex() != 0:
            return
        avec_arcs = self.case_arcs.isChecked() and self.arcs_roues is not None
        limites = (self.ax.get_xlim(), self.ax.get_ylim()) if self._vue_ok else self.limites_c
        self.ax.clear()
        evolvente.dessiner_engrenage(
            self.ax, self.roues, self.profils, N_DENTS, limites=limites,
            profils_arcs=self.arcs_echantillons if avec_arcs else None,
            points_tangence=[(a.points, a.normales()) for a in self.arcs_roues] if avec_arcs else None,
        )
        self._vue_ok = True
        self.canvas.draw_idle()

    # --- onglet Simulation ----------------------------------------------
    def _arcs_de_creation(self):
        """Mobiles (arcs) construits d'après l'onglet Création ; False si l'onglet est invalide."""
        try:
            roues = self.lire_roues()
            profils = [evolvente.evolve(r, NI) for r in roues]
            arcs_roues = [self._arcs_de(r, p) for r, p in zip(roues, profils)]
        except (ValueError, ArithmeticError, IndexError) as e:
            self._message(f"Arcs de la création inutilisables : {e}", erreur=True)
            return False
        self.mobiles = tuple(arcs.vers_mobile(r, a, role) for r, a, role in zip(roues, arcs_roues, "me"))
        self.mobiles_points = None
        self.origine = "arcs de l'onglet Création"
        return True

    def utiliser_arcs_creation(self):
        """Simule les arcs de cercle construits dans l'onglet Création."""
        self._aller_onglet(1)
        if self._arcs_de_creation():
            self.simuler()

    def ouvrir_arcs(self):
        ref = self.chemin or self.source
        dossier = ref.parent if ref else RACINE / "Files"
        nom, _ = QFileDialog.getOpenFileName(self, "Ouvrir un fichier d'arcs", str(dossier), FILTRE_ARCS)
        if nom:
            self.charger_arcs(Path(nom))

    def charger_arcs(self, chemin):
        try:
            mobiles, extras = arcs.lire_arcs(chemin)
        except (OSError, ValueError, KeyError, IndexError, StopIteration, TypeError) as e:
            QMessageBox.critical(self, "Open File", f"Lecture impossible de {chemin.name} :\n{e}")
            return False
        self.mobiles, self.extras, self.origine = mobiles, extras, chemin.name
        self.mobiles_points = None
        manque = False
        for k, m in enumerate(mobiles):  # données de Hertz absentes : valeurs par défaut
            for attr, valeur in zip(("e_mod", "poisson", "epaisseur"), HERTZ_DEFAUT[k]):
                if not getattr(m, attr):
                    setattr(m, attr, valeur)
                    manque = True
        self.origine += " (Hertz par défaut)" if manque else ""
        for cle, v in (("da_neg", extras["dA_neg"]), ("da_pos", extras["dA_pos"])):
            self.edits_sim[cle].setText(texte_nombre(v))
        noms = {nom: cle for cle, _, nom in OBJECTIFS}
        sim_extras = extras["simulation"]
        for combo, cle_fichier in ((self.combo_principal, "opti_1"), (self.combo_secondaire, "opti_2")):
            i = combo.findData(noms.get(sim_extras.get(cle_fichier)))
            if i >= 0:
                combo.setCurrentIndex(i)
        try:
            self.edits_opti["ponderation"].setText(texte_nombre(float(sim_extras["ponderation"])))
        except (KeyError, ValueError, TypeError):
            pass
        self.projet.optimisation["libres"] = {"m": None, "e": None}
        self._aller_onglet(1)
        self.simuler()
        return True

    def _profils_courants(self):
        """(mobiles ou None, roues, profils échantillonnés, demi-profils arrière ou None) du profil
        chargé dans l'onglet Simulation, qu'il soit en arcs de cercle ou en suite de points."""
        if self.mobiles_points is not None:
            p = self.mobiles_points
            roues = [evolvente.Roue(m.mo, m.zz, None, None, 0.0, None, None,
                                    m.e_mod, m.poisson, m.epaisseur) for m in p]
            return None, roues, [list(m.points) for m in p], None
        roues = [evolvente.Roue(m.mo, m.zz, None, None, 0.0, None, None,
                                m.e_mod, m.poisson, m.epaisseur) for m in self.mobiles]
        profils = [m.echantillonner().tolist() for m in self.mobiles]
        dos = (None if all(m.symetrique for m in self.mobiles)
               else [m.echantillonner_dos().tolist() for m in self.mobiles])
        return self.mobiles, roues, profils, dos

    def simuler(self):
        """Simule l'engrènement du profil courant sur une période : contact entre arcs de cercle,
        ou entre lignes brisées si le profil est une suite de points (roues taillées)."""
        if self.mobiles is None and self.mobiles_points is None:
            self._message("Ouvrez un fichier d'arcs ou utilisez les arcs de l'onglet Création.", erreur=True)
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            variables = self.lire_simulation()
            mobiles, roues, profils, dos = self._profils_courants()
            da = variables["da"] / 1000.0  # µm -> mm
            try:
                moteur = contact.Engrenement(roues, profils, da, mobiles=mobiles)
                sim, avertissement = moteur.simuler(variables["mu"], N_POS), None
            except contact.ErreurContact as e:
                sim, avertissement = None, f"Pas de simulation d'engrènement : {e}."
            limites = superposition = None
            positions = [sim.position(i) for i in (0, N_POS // 2, N_POS - 1)] if sim else [evolvente.SansContact()]
            if sim:
                superposition = (sim.contact[:, 0].tolist(), {k: v.tolist() for k, v in sim.pourcent().items()}, 0)
            limites = evolvente.limites_engrenage(
                roues, profils, N_DENTS, positions, evolvente.points_courbes(superposition, roues[0].mo), da, dos)
        except (ValueError, ArithmeticError) as e:
            self._message(f"Simulation impossible : {e}", erreur=True)
            return
        finally:
            QApplication.restoreOverrideCursor()
        for cle in ("mu", "mt", "da", "da_neg", "da_pos"):
            setattr(self.projet, cle, variables[cle])
        self.mt, self.da_mm = variables["mt"], da
        self.sim, self.sim_roues, self.sim_profils, self.sim_profils_dos = sim, roues, profils, dos
        self.sim_limites, self.superposition = limites, superposition
        self._vue_ok = False
        self.curseur.blockSignals(True)
        self.curseur.setValue(0)
        self.curseur.setEnabled(sim is not None)
        self.curseur.blockSignals(False)
        m1, m2 = self.mobiles_points or self.mobiles
        if self.mobiles_points is not None:
            forme = "<small>(contact calculé entre les suites de points, sans approximation)</small>"
            taille = [f"{len(m.points)} points" for m in (m1, m2)]
        else:
            forme = "<small>(contact calculé entre les arcs de cercle, sans points)</small>"
            taille = [f"{len(m.points) - 1} arcs" for m in (m1, m2)]
        self.etiq_origine.setText(
            f"<b>{self.origine}</b> {forme}"
            + ("<br><b>dents dissymétriques</b> <small>(demi-profil arrière propre à chaque dent)</small>" if dos else "")
            + f"<br>Menante {m1.nom} : Z={m1.zz}, m={m1.mo:g}, {taille[0]}<br>"
            f"Menée {m2.nom} : Z={m2.zz}, m={m2.mo:g}, {taille[1]}")
        self._remplir_variables()
        self._remplir_arcs_sim()
        self.afficher_position(0)
        self._afficher_statistiques()
        self.barre.update()  # nouvelle vue d'origine pour le zoom
        if avertissement:
            self._message(avertissement, erreur=True)
        else:
            self._message(f"Simulation faite ({N_POS} positions, ΔA = {variables['da']:+g} µm).")

    def _afficher_statistiques(self):
        """Variation des trois paramètres sur la période, en % de leur valeur idéale."""
        if self.sim is None:
            self.stats.setText("-")
            return
        noms = {"rapport": "Transmission", "couple": "Couple", "rendement": "Rendement"}
        lignes = "".join(
            f"<tr><td>{noms[k]}</td><td align=right>{mn:.2f}</td><td align=right>{moy:.2f}</td>"
            f"<td align=right>{mx:.2f}</td><td align=right><b>{var:.2f}</b></td></tr>"
            for k, (mn, moy, mx, var) in self.sim.statistiques().items()
        )
        self.stats.setText(
            "<table cellspacing=4><tr><th></th><th>min</th><th>moy.</th><th>max</th><th>variation</th></tr>"
            f"{lignes}</table><small>en % de la valeur idéale ; variation = 2(max-min)/(max+min)</small><br>"
            f"Jeu de la menée : <b>{self.sim.jeu_min:.3f}° à {self.sim.jeu_max:.3f}°</b>"
            + (f" <b style='color:{theme.ERREUR}'>!! coincement</b>" if self.sim.jeu_min <= 0 else "")
            + self._texte_hertz()
        )

    def _texte_hertz(self):
        """Pression de Hertz maximale (pour 1 mN·m) et où elle se produit sur les deux profils."""
        pmax, i = self.sim.hertz_max()
        if i < 0:
            return ""
        (x1, y1), (x2, y2) = self.sim.hertz_local
        return (f"<br>Pression de Hertz max (1 mN·m) : <b>{pmax:.1f} N/mm²</b> "
                f"<span style='color:{theme.ERREUR}'>★</span> à la position {i + 1}/{len(self.sim)}<br>"
                f"<small>rayon {math.hypot(x1, y1):.4f} mm sur la menante, {math.hypot(x2, y2):.4f} mm sur la menée</small>")

    def afficher_position(self, i):
        """Dessine les roues et les courbes à la position `i` du curseur, met à jour les valeurs."""
        if self.sim_roues is None or self.onglets.currentIndex() != 1:
            return
        sim = self.sim
        position = sim.position(i) if sim else None
        trajet = sim.contact.tolist() if sim else None
        superposition = (self.superposition[0], self.superposition[1], i) if sim else None
        limites = (self.ax.get_xlim(), self.ax.get_ylim()) if self._vue_ok else self.sim_limites
        self.ax.clear()
        evolvente.dessiner_engrenage(
            self.ax, self.sim_roues, self.sim_profils, N_DENTS, position, trajet, limites,
            superposition, self.da_mm, None,
            # points de tangence des arcs : aucun si le profil courant est une suite de points
            None if self.mobiles is None else
            [(m.points, m.normales()) if m.symetrique else (m.points, m.normales(), m.points_dos, m.normales_dos())
             for m in self.mobiles],
            sim.paire.tolist() if sim else None, sim.hertz_local if sim else None, self.sim_profils_dos)
        self._vue_ok = True
        self.canvas.draw_idle()
        if sim is None:
            self.etiq_pos.setText("-")
            self.resultats.setText("-")
            return
        z1, z2 = self.sim_roues[0].zz, self.sim_roues[1].zz
        self.etiq_pos.setText(
            f"Position {i + 1} / {len(sim)} : menante tournée de "
            f"{math.degrees(sim.alpha1[i] - sim.alpha1[0]):.3f}° sur {360.0 / z1:.3f}°")
        pct = {k: v[i] for k, v in sim.pourcent().items()}
        cx, cy = sim.contact[i]

        def ligne(nom, valeur, pc):
            return f"{nom} : {valeur}<b>{pc:.2f} %</b> <small>(variation {pc - 100:+.2f} %)</small><br>"

        self.resultats.setText(
            ligne("Rapport de transmission ω₂/ω₁", f"{sim.rapport[i]:.4f} = ", pct["rapport"])
            + ligne("Rapport de couple M₂/M₁", f"{sim.couple[i]:.4f} = ", pct["couple"])
            + ligne("Rendement", "", pct["rendement"])
            + f"Couple mené : <b>{abs(sim.couple[i]) * self.mt:.4f} mN·m</b> pour {self.mt:g} mN·m sur la menante<br>"
            f"Force normale : <b>{sim.normale[i] * self.mt:.3f} N</b> <small>(μ = {sim.mu:g})</small><br>"
            + (f"Pression de Hertz : <b>{sim.hertz[i]:.1f} N/mm²</b> <small>(pour 1 mN·m)</small><br>"
               if sim.hertz is not None else "")
            + f"Point d'engrènement : x = {cx:+.4f}, y = {cy:+.4f} mm<br>"
            f"Paires de dents en prise : <b>{'2' if sim.contact2[i, 0] == sim.contact2[i, 0] else '1'}</b>"
        )

    # --- fichiers --------------------------------------------------------
    def charger(self, chemin):
        """Ouvre un projet de création (.json, .Cev) ou un fichier d'arcs (.json, .Eng, .Eno)."""
        ext = chemin.suffix.lower()
        if ext in (".eng", ".eno") or (ext == ".json" and _json_est_arcs(chemin)):
            return self.charger_arcs(chemin)
        try:
            projet = evolvente.lire_projet(chemin)
        except (OSError, ValueError, IndexError) as e:  # JSONDecodeError est une ValueError
            QMessageBox.critical(self, "Open File", f"Lecture impossible de {chemin.name} :\n{e}")
            return False
        est_json = ext == ".json"
        self.projet, self.modifie = projet, False
        self.chemin, self.source = (chemin, None) if est_json else (None, chemin)
        # la simulation repartira des nouveaux arcs
        self.mobiles = self.mobiles_points = self.extras = self.sim = None
        self.remplir_panneau()
        self._aller_onglet(0)
        self.reconstruire()
        self._titre()
        return True

    def _confirmer_perte(self):
        if not self.modifie:
            return True
        r = QMessageBox.question(
            self, "Modifications non enregistrées",
            "Enregistrer les modifications ?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
        )
        if r == QMessageBox.Save:
            return self.sauver()
        return r == QMessageBox.Discard

    def ouvrir(self):
        if not self._confirmer_perte():
            return
        ref = self.chemin or self.source
        dossier = ref.parent if ref else RACINE / "Files"
        nom, _ = QFileDialog.getOpenFileName(
            self, "Open File", str(dossier),
            "Engrenages (*.json *.Cev *.cev *.Eng *.eng *.Eno *.eno);;Tous les fichiers (*)",
        )
        if nom:
            self.charger(Path(nom))

    def sauver(self):
        return self._ecrire(self.chemin) if self.chemin else self.sauver_sous()

    def sauver_sous(self):
        if self.chemin:
            proposition = self.chemin
        elif self.source:
            proposition = self.source.with_suffix(".json")
        else:
            proposition = RACINE / "Files" / "nouveau.json"
        nom, _ = QFileDialog.getSaveFileName(
            self, "Save As", str(proposition), "JSON (*.json);;Tous les fichiers (*)"
        )
        if not nom:
            return False
        chemin = Path(nom)
        if chemin.suffix.lower() != ".json":
            chemin = chemin.with_name(chemin.name + ".json")
        return self._ecrire(chemin)

    def _ecrire(self, chemin):
        """Enregistre le projet de création (paramètres des roues et de la simulation)."""
        try:
            self.projet.roues = self.lire_roues()
            for cle, v in self.lire_simulation().items():
                setattr(self.projet, cle, v)
            self.projet.optimisation = self.lire_optimisation()
            evolvente.ecrire_json(chemin, self.projet)
        except ValueError as e:
            QMessageBox.warning(self, "Save", f"Enregistrement impossible : {e}")
            return False
        except OSError as e:
            QMessageBox.critical(self, "Save", f"Écriture impossible :\n{e}")
            return False
        self.chemin, self.modifie = chemin, False
        self._titre()
        self._message(f"Enregistré : {chemin}")
        return True

    def exporter_arcs(self):
        """Écrit les profils en arcs de cercle (.json au schéma profile.json, ou .Eno)."""
        try:
            if self.onglets.currentIndex() == 0:
                roues = self.lire_roues()
                profils = [evolvente.evolve(r, NI) for r in roues]
                arcs_roues = [self._arcs_de(r, p) for r, p in zip(roues, profils)]
                mobiles = tuple(arcs.vers_mobile(r, a, role) for r, a, role in zip(roues, arcs_roues, "me"))
            elif self.mobiles is not None:
                mobiles = self.mobiles
            elif self.mobiles_points is not None:
                raise ValueError(MSG_ARCS_REQUIS)
            else:
                raise ValueError("aucun profil en arcs de cercle")
        except (ValueError, ArithmeticError, IndexError) as e:
            QMessageBox.warning(self, "Export Arcs", f"Export impossible : {e}")
            return
        self._ecrire_arcs(mobiles)

    def _ecrire_arcs(self, mobiles, suffixe="arcs"):
        """Demande un fichier et y écrit `mobiles` (.json, ou .Eno si les dents sont symétriques)."""
        ref = self.chemin or self.source
        proposition = ref.with_name(f"{ref.stem}_{suffixe}.json") if ref else RACINE / "Files" / f"{suffixe}.json"
        nom, filtre = QFileDialog.getSaveFileName(
            self, "Export Arcs", str(proposition), "JSON (*.json);;Fichiers Eno (*.Eno *.eno)"
        )
        if not nom:
            return
        chemin = Path(nom)
        if not chemin.suffix:
            chemin = chemin.with_suffix(".Eno" if "Eno" in filtre else ".json")
        try:
            if chemin.suffix.lower() == ".eno":
                if any(m.zz < 0 for m in mobiles):  # le format Fortran ne connaît pas les couronnes
                    raise ValueError("une couronne intérieure (Z < 0) ne peut pas s'écrire en .Eno")
                arcs.ecrire_eno(chemin, mobiles, self._extras_courants())
            else:
                arcs.ecrire_json(chemin, mobiles, self._extras_courants())
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, "Export Arcs", f"Écriture impossible : {e}")
            return
        self._message(f"Arcs exportés : {chemin} ({len(mobiles[0].points)} et {len(mobiles[1].points)} points).")

    def exporter_rapport(self):
        """Écrit un rapport de simulation (.docx) pour le profil actuellement chargé (onglet
        Simulation) : conditions, résultats sur la plage d'entraxe par pas de 10 µm, graphiques."""
        if self.mobiles is None:
            QMessageBox.information(self, "Rapport", MSG_ARCS_REQUIS if self.mobiles_points is not None
                                    else "Ouvrez ou simulez d'abord un profil en arcs de cercle.")
            return
        try:
            v = self.lire_simulation()
        except ValueError as e:
            QMessageBox.warning(self, "Rapport", f"Paramètres de simulation invalides : {e}")
            return
        ref = self.chemin or self.source
        proposition = ref.with_name(f"{ref.stem}_rapport.docx") if ref else RACINE / "Result" / "rapport.docx"
        nom, _ = QFileDialog.getSaveFileName(self, "Exporter un rapport", str(proposition), "Word (*.docx)")
        if not nom:
            return
        chemin = Path(nom).with_suffix(".docx")
        roues = [evolvente.Roue(m.mo, m.zz, None, None, 0.0, None, None) for m in self.mobiles]
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from . import rapport
            _, resume = rapport.ecrire_rapport(
                chemin, roues, self.mobiles, v["mu"], v["mt"], v["da_neg"], v["da_pos"],
                optimisation=self.projet.optimisation, origine=self.origine)
        except (ValueError, ArithmeticError, OSError) as e:
            QMessageBox.warning(self, "Rapport", f"Écriture impossible : {e}")
            return
        finally:
            QApplication.restoreOverrideCursor()
        detail = (f" Pire cas : couple {resume['varC_max']:.2f} %, Hertz {resume['hertz_max']:.0f} N/mm², "
                 f"jeu {resume['jeu_min']:.2f}° à {resume['jeu_max']:.2f}°." if resume else "")
        self._message(f"Rapport écrit : {chemin}.{detail}")

    def exporter_profils_dxf(self):
        """Exporte le contour complet (toutes les dents) des deux profils de l'onglet Simulation
        en DXF, un fichier par roue : de vraies entités ARC si le profil est une chaîne d'arcs,
        une polyligne si c'est une suite de points (roue taillée)."""
        if self.mobiles is None and self.mobiles_points is None:
            QMessageBox.information(self, "DXF", "Ouvrez ou simulez d'abord un profil.")
            return
        ref = self.chemin or self.source
        proposition = ref.with_name(f"{ref.stem}.dxf") if ref else RACINE / "Files" / "profil.dxf"
        nom, _ = QFileDialog.getSaveFileName(self, "Exporter les profils en DXF", str(proposition), "DXF (*.dxf)")
        if not nom:
            return
        base = Path(nom)
        if base.suffix.lower() == ".dxf":
            base = base.with_suffix("")
        from . import dxf as dxf_mod
        chemins = [base.with_name(f"{base.name}_{suffixe}.dxf") for suffixe in ("menante", "menee")]
        try:
            if self.mobiles_points is not None:
                for chemin, p in zip(chemins, self.mobiles_points):
                    dxf_mod.ecrire_polyligne(chemin, p.points, p.zz)
            else:
                for chemin, mobile in zip(chemins, self.mobiles):
                    dxf_mod.ecrire_profil(chemin, mobile)
        except OSError as e:
            QMessageBox.warning(self, "DXF", f"Écriture impossible : {e}")
            return
        forme = "polyligne" if self.mobiles_points is not None else "arcs de cercle"
        self._message(f"Profils exportés en DXF ({forme}) : {chemins[0].name} et {chemins[1].name}.")

    def closeEvent(self, event):
        if self._confirmer_perte():
            if self.travail is not None:  # optimisation en cours : on l'interrompt proprement
                self.travail.arreter()
                self.travail.wait()
            event.accept()
        else:
            event.ignore()


def main():
    app = QApplication(sys.argv)
    theme.appliquer(app)
    fenetre = Fenetre(sys.argv[1] if len(sys.argv) > 1 else None)
    fenetre.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
