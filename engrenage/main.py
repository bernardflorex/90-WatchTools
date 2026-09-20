"""Interface graphique : paramètres des roues à gauche, engrenage et courbes à droite.

Usage : python main.py [fichier.json | ancien fichier.Cev]   (nécessite PySide6, matplotlib, numpy)
"""

import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import QLocale, Qt
from PySide6.QtGui import QAction, QDoubleValidator, QIntValidator, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSlider,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

import arcs
import contact
import evolvente

RACINE = Path(__file__).resolve().parent
FICHIER_DEFAUT = RACINE / "test.Cev"
NI = 150  # points par demi-dent
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


# Une ligne du panneau par champ : pour ajouter une variable de roue, ajouter un Champ ici
# (et l'attribut correspondant dans evolvente.Roue).
CHAMPS = (
    Champ("mo", "MO", "Module (mm)", 0.001, 1000),
    Champ("zz", "ZZ", "Nombre de dents", 3, 10000, entier=True),
    Champ("nw", "NW", "Angle de pression (°)", 1, 60),
    Champ("qe", "QE", "Rapport creux / plein au diamètre primitif", 0.1, 10),
    Champ("qg", "QG", "Raccord de pied (x module), 0 = aucun", 0, 5),
    Champ("ma", "MA", "Hauteur de développante sous le primitif (x module)", 0, 5),
    Champ("me", "ME", "Hauteur de développante au-dessus du primitif (x module)", 0, 5),
)
# Export en arcs de cercle (.Eno) et matériau (données de Hertz du fichier .Eno), par roue.
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
    Champ("da", "ΔA", "Variation de l'entraxe (µm)", -1000, 1000),
)

STYLE_ERREUR = "QLineEdit { border: 1px solid #c62828; background: #ffebee; }"


def texte_nombre(v):
    return f"{v:.6f}".rstrip("0").rstrip(".") or "0"


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


class Fenetre(QMainWindow):
    def __init__(self, chemin=None):
        super().__init__()
        self.chemin = None  # fichier .json de travail (None : jamais enregistré)
        self.source = None  # ancien fichier .Cev importé, s'il y en a un
        self.projet = evolvente.Projet(roues=[evolvente.Roue(0.135, 102, 20, 1.2, 0, 1, 0.5),
                                              evolvente.Roue(0.135, 10, 20, 1.2, 0.2, 1, 1)])
        self.modifie = False
        self.roues = self.profils = self.sim = self.limites = self.superposition = None
        self.da_mm, self.arcs_roues, self.arcs_echantillons = 0.0, None, None
        self._cache_arcs = {}  # profil en arcs par (paramètres de la roue, nombre d'arcs)
        self._vue_ok = False  # la vue du graphique a été fixée (zoom de l'utilisateur conservé)
        self.resize(1400, 860)

        self._creer_menu()
        gauche = self._creer_panneau()
        droite = self._creer_graphiques()
        cadre = QSplitter(Qt.Horizontal)
        cadre.addWidget(gauche)
        cadre.addWidget(droite)
        cadre.setStretchFactor(1, 1)
        cadre.setSizes([430, 970])
        self.setCentralWidget(cadre)

        if chemin and not self.charger(Path(chemin)):
            chemin = None
        if not chemin and FICHIER_DEFAUT.exists():
            self.charger(FICHIER_DEFAUT)
        elif not chemin:
            self.remplir_panneau()
            self.reconstruire()
        self._titre()

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
        export = QAction("&Export Arcs (.Eno)…", self)
        export.setShortcut("Ctrl+E")
        export.triggered.connect(self.exporter_eno)
        menu.addAction(export)
        menu.addSeparator()
        quitter = QAction("&Quit", self)
        quitter.setShortcut(QKeySequence.Quit)
        quitter.triggered.connect(self.close)
        menu.addAction(quitter)

    def _creer_panneau(self):
        self.edits = ({}, {})
        groupe = self._groupe_roues("Profils", CHAMPS)
        groupe_export = self._groupe_roues("Arcs de cercle et matériau (export .Eno)", CHAMPS_EXPORT)

        # simulation : variables communes
        sim = QGroupBox("Simulation")
        grille_sim = QGridLayout(sim)
        grille_sim.setColumnStretch(0, 1)
        self.edits_sim = {}
        for i, c in enumerate(CHAMPS_SIM):
            etiquette = QLabel(f"<b>{c.nom}</b>  {c.aide}")
            etiquette.setWordWrap(True)
            grille_sim.addWidget(etiquette, i, 0)
            le = _champ_saisie(c, c.nom)
            le.setMaximumWidth(110)
            le.textEdited.connect(lambda _t, le=le: self._edite(le))
            grille_sim.addWidget(le, i, 1)
            self.edits_sim[c.cle] = le

        # arcs de cercle : affichage et écart d'approximation
        cadre_arcs = QGroupBox("Arcs de cercle tangents")
        lay_arcs = QVBoxLayout(cadre_arcs)
        self.case_arcs = QCheckBox("Afficher sur le dessin (pointillés)")
        self.case_arcs.toggled.connect(lambda _c: self.afficher_position(self.curseur.value()))
        self.info_arcs = QLabel("-")
        self.info_arcs.setWordWrap(True)
        self.info_arcs.setTextFormat(Qt.RichText)
        lay_arcs.addWidget(self.case_arcs)
        lay_arcs.addWidget(self.info_arcs)

        contenu = QWidget()
        lay_contenu = QVBoxLayout(contenu)
        lay_contenu.addWidget(groupe)
        lay_contenu.addWidget(sim)
        lay_contenu.addWidget(groupe_export)
        lay_contenu.addWidget(cadre_arcs)
        lay_contenu.addStretch(1)
        zone = QScrollArea()
        zone.setWidgetResizable(True)
        zone.setWidget(contenu)

        self.bouton = QPushButton("Reconstruire l'engrenage  (F5)")
        self.bouton.setShortcut("F5")
        self.bouton.setMinimumHeight(36)
        self.bouton.clicked.connect(self.reconstruire)

        # curseur des positions d'engrènement
        self.curseur = QSlider(Qt.Horizontal)
        self.curseur.setRange(0, N_POS - 1)
        self.curseur.setEnabled(False)
        self.curseur.valueChanged.connect(self.afficher_position)
        self.etiq_pos = QLabel("-")
        cadre_pos = QGroupBox("Position d'engrènement")
        lay_pos = QVBoxLayout(cadre_pos)
        lay_pos.addWidget(self.etiq_pos)
        lay_pos.addWidget(self.curseur)
        self.resultats = QLabel()
        self.resultats.setWordWrap(True)
        self.resultats.setTextFormat(Qt.RichText)
        lay_pos.addWidget(self.resultats)
        self.stats = QLabel()
        self.stats.setTextFormat(Qt.RichText)
        cadre_stats = QGroupBox("Variation sur la période (100 % = idéal)")
        QVBoxLayout(cadre_stats).addWidget(self.stats)

        self.statut = QLabel()
        self.statut.setWordWrap(True)

        panneau = QWidget()
        lay = QVBoxLayout(panneau)
        lay.addWidget(zone, 1)
        lay.addWidget(self.bouton)
        lay.addWidget(cadre_pos)
        lay.addWidget(cadre_stats)
        lay.addWidget(self.statut)
        panneau.setMinimumWidth(340)
        return panneau

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

    def _creer_graphiques(self):
        # marges fixes : la mise en page automatique est trop lente pour suivre le curseur
        self.figure = Figure(figsize=(8, 6))
        self.figure.subplots_adjust(left=0.06, right=0.98, bottom=0.06, top=0.91)
        self.ax = self.figure.add_subplot()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.barre = NavigationToolbar2QT(self.canvas, self)
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

    def _lire(self, edits, champs, nom_de):
        """Valeurs des champs ; ValueError (après avoir marqué les champs) si l'un est invalide."""
        valeurs, erreur = {}, None
        for c in champs:
            le = edits[c.cle]
            if not le.hasAcceptableInput():
                le.setStyleSheet(STYLE_ERREUR)
                erreur = erreur or f"{c.nom} ({nom_de}) : valeur entre {c.mini:g} et {c.maxi:g} attendue"
                continue
            le.setStyleSheet("")
            valeurs[c.cle] = int(le.text()) if c.entier else float(le.text())
        return valeurs, erreur

    def lire_panneau(self):
        """(roues, variables de simulation) d'après le panneau ; ValueError si une valeur
        est invalide."""
        roues, erreur = [], None
        for k in range(2):
            valeurs, e = self._lire(self.edits[k], CHAMPS + CHAMPS_EXPORT, ROUES[k].lower())
            erreur = erreur or e
            if not e:
                roues.append(replace(self.projet.roues[k], **valeurs))  # garde les données Hertz
        valeurs, e = self._lire(self.edits_sim, CHAMPS_SIM, "simulation")
        erreur = erreur or e
        if erreur:
            raise ValueError(erreur)
        return roues, valeurs

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
        self.statut.setStyleSheet("color: #c62828;" if erreur else "color: #2e7d32;")
        self.statut.setText(texte)

    # --- construction ----------------------------------------------------
    def reconstruire(self):
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            roues, variables = self.lire_panneau()
            profils = [evolvente.evolve(r, NI) for r in roues]
            da = variables["da"] / 1000.0  # µm -> mm
            try:
                sim = contact.Engrenement(roues, profils, da).simuler(variables["mu"], N_POS)
                avertissement = None
            except contact.ErreurContact as e:
                sim, avertissement = None, f"Pas de simulation d'engrènement : {e}."
            limites = superposition = None
            if sim:
                pourcent = sim.pourcent()
                superposition = (sim.contact[:, 0].tolist(), {k: v.tolist() for k, v in pourcent.items()}, 0)
                limites = evolvente.limites_engrenage(
                    roues, profils, N_DENTS, [sim.position(i) for i in (0, N_POS // 2, N_POS - 1)],
                    evolvente.points_courbes(superposition, roues[0].mo), da,
                )
            arcs_roues, echantillons, info = self._calculer_arcs(roues, profils)
        except (ValueError, ArithmeticError) as e:
            self._message(f"Construction impossible : {e}", erreur=True)
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.projet.roues = roues
        self.projet.mu, self.projet.mt, self.projet.da = variables["mu"], variables["mt"], variables["da"]
        self.roues, self.profils, self.sim, self.limites = roues, profils, sim, limites
        self.superposition, self.da_mm = superposition, da
        self.arcs_roues, self.arcs_echantillons = arcs_roues, echantillons
        self.info_arcs.setText(info)
        self._vue_ok = False
        self.curseur.blockSignals(True)
        self.curseur.setValue(0)
        self.curseur.setEnabled(sim is not None)
        self.curseur.blockSignals(False)
        self.afficher_position(0)
        self._afficher_statistiques()
        self.barre.update()  # nouvelle vue d'origine pour le zoom
        if avertissement:
            self._message(avertissement, erreur=True)
        else:
            self._message(f"Engrenage construit ({NI} points par profil, {N_POS} positions).")

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
            return None, None, f"<span style='color:#c62828'>Approximation impossible : {e}</span>"
        info = "<br>".join(
            f"{nom} : {a.n_points - 1} arcs, écart max <b>{a.ecart_max * 1000:.2f} µm</b>"
            + (f" <span style='color:#c62828'>(tangente finale à {math.degrees(a.erreur_angle):.2f}°)</span>"
               if abs(a.erreur_angle) > math.radians(0.05) else "")
            for nom, a in zip(ROUES, arcs_roues)
        )
        echantillons = [arcs.echantillonner(a).tolist() for a in arcs_roues]
        return arcs_roues, echantillons, info

    def exporter_eno(self):
        """Écrit les profils en arcs de cercle tangents au format .Eno (suite du programme)."""
        try:
            roues, _ = self.lire_panneau()
            profils = [evolvente.evolve(r, NI) for r in roues]
            arcs_roues = [self._arcs_de(r, p) for r, p in zip(roues, profils)]
        except (ValueError, ArithmeticError, IndexError) as e:
            QMessageBox.warning(self, "Export Arcs", f"Export impossible : {e}")
            return
        ref = self.chemin or self.source
        proposition = ref.with_suffix(".Eno") if ref else RACINE / "Files" / "nouveau.Eno"
        nom, _ = QFileDialog.getSaveFileName(
            self, "Export Arcs", str(proposition), "Fichiers Eno (*.Eno *.eno);;Tous les fichiers (*)"
        )
        if not nom:
            return
        chemin = Path(nom)
        if not chemin.suffix:
            chemin = chemin.with_suffix(".Eno")
        try:
            arcs.ecrire_eno(chemin, roues, arcs_roues)
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, "Export Arcs", f"Écriture impossible : {e}")
            return
        pire = max(a.ecart_max for a in arcs_roues) * 1000
        self._message(f"Arcs exportés : {chemin} ({arcs_roues[0].n_points} et "
                      f"{arcs_roues[1].n_points} points, écart max {pire:.2f} µm).")

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
            + (" <b style='color:#c62828'>!! coincement</b>" if self.sim.jeu_min <= 0 else "")
        )

    def afficher_position(self, i):
        """Dessine les roues et les courbes à la position `i` du curseur, met à jour les valeurs."""
        if self.roues is None:
            return
        sim = self.sim
        position = sim.position(i) if sim else None
        trajet = sim.contact.tolist() if sim else None
        superposition = (self.superposition[0], self.superposition[1], i) if sim else None
        limites = (self.ax.get_xlim(), self.ax.get_ylim()) if self._vue_ok else self.limites
        self.ax.clear()
        profils_arcs = self.arcs_echantillons if self.case_arcs.isChecked() else None
        evolvente.dessiner_engrenage(self.ax, self.roues, self.profils, N_DENTS, position, trajet,
                                     limites, superposition, self.da_mm, profils_arcs)
        self._vue_ok = True
        self.canvas.draw_idle()
        if sim is None:
            self.etiq_pos.setText("-")
            self.resultats.setText("-")
            return
        z1 = self.roues[0].zz
        pas = 360.0 / z1
        self.etiq_pos.setText(
            f"Position {i + 1} / {len(sim)} : menante tournée de "
            f"{math.degrees(sim.alpha1[i] - sim.alpha1[0]):.3f}° sur {pas:.3f}°"
        )
        pct = {k: v[i] for k, v in sim.pourcent().items()}
        mt = self.projet.mt
        cx, cy = sim.contact[i]
        def ligne(nom, valeur, pc):
            return f"{nom} : {valeur}<b>{pc:.2f} %</b> <small>(variation {pc - 100:+.2f} %)</small><br>"

        self.resultats.setText(
            ligne("Rapport de transmission ω₂/ω₁", f"{sim.rapport[i]:.4f} = ", pct["rapport"])
            + ligne("Rapport de couple M₂/M₁", f"{sim.couple[i]:.4f} = ", pct["couple"])
            + ligne("Rendement", "", pct["rendement"])
            + f"Couple mené : <b>{abs(sim.couple[i]) * mt:.4f} mN·m</b> pour {mt:g} mN·m sur la menante<br>"
            f"Force normale : <b>{sim.normale[i] * mt:.3f} N</b> <small>(μ = {sim.mu:g})</small><br>"
            f"Point d'engrènement : x = {cx:+.4f}, y = {cy:+.4f} mm<br>"
            f"Paires de dents en prise : <b>{'2' if sim.contact2[i, 0] == sim.contact2[i, 0] else '1'}</b>"
        )

    # --- fichiers --------------------------------------------------------
    def charger(self, chemin):
        try:
            projet = evolvente.lire_projet(chemin)
        except (OSError, ValueError, IndexError) as e:  # JSONDecodeError est une ValueError
            QMessageBox.critical(self, "Open File", f"Lecture impossible de {chemin.name} :\n{e}")
            return False
        est_json = chemin.suffix.lower() == ".json"
        self.projet, self.modifie = projet, False
        self.chemin, self.source = (chemin, None) if est_json else (None, chemin)
        self.remplir_panneau()
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
            "Engrenages (*.json *.Cev *.cev);;JSON (*.json);;Ancien format (*.Cev *.cev);;Tous les fichiers (*)",
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
        try:
            self.projet.roues, variables = self.lire_panneau()
            self.projet.mu, self.projet.mt, self.projet.da = variables["mu"], variables["mt"], variables["da"]
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

    def closeEvent(self, event):
        if self._confirmer_perte():
            event.accept()
        else:
            event.ignore()


def main():
    app = QApplication(sys.argv)
    fenetre = Fenetre(sys.argv[1] if len(sys.argv) > 1 else None)
    fenetre.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
