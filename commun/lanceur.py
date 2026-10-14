"""Fenetre d'accueil : une icone par outil, un clic le lance.

Chaque outil est un programme a part, lance dans son propre processus : il
evolue, plante ou se ferme sans toucher aux autres, et l'accueil ne charge
aucun de leurs modules.
"""

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from PySide6.QtCore import QProcess, QSize, Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QHBoxLayout, QLabel, QMainWindow, QToolButton,
                               QVBoxLayout, QWidget)

from . import icones

TITRE = "WatchTools"
RACINE = Path(__file__).resolve().parent.parent


@dataclass(frozen=True)
class Outil:
    """Un outil de l'accueil. `module` se lance par `python -m module`."""
    cle: str
    nom: str
    resume: str
    module: str
    dessin: Callable
    couleur: QColor

    def lancer(self):
        """Lance l'outil dans son propre processus ; False s'il ne part pas."""
        ok, _pid = QProcess.startDetached(sys.executable, ["-m", self.module],
                                          str(RACINE))
        return ok


class Lanceur(QMainWindow):

    COTE = 128

    def __init__(self, outils):
        super().__init__()
        self.setWindowTitle(TITRE)
        centre = QWidget()
        lay = QVBoxLayout(centre)
        lay.setContentsMargins(28, 22, 28, 12)
        titre = QLabel("Outils de calcul horloger")
        police = titre.font()
        police.setPointSizeF(police.pointSizeF() * 1.5)
        titre.setFont(police)
        lay.addWidget(titre)
        rang = QHBoxLayout()
        rang.setSpacing(18)
        for o in outils:
            rang.addWidget(self._bouton(o))
        rang.addStretch(1)
        lay.addSpacing(10)
        lay.addLayout(rang)
        lay.addStretch(1)
        self.setCentralWidget(centre)
        self.statusBar().showMessage("Cliquez sur un outil pour l'ouvrir.")

    def _bouton(self, outil):
        b = QToolButton()
        b.setIcon(icones.icone(outil.dessin, outil.couleur, self.COTE))
        b.setIconSize(QSize(self.COTE, self.COTE))
        b.setText(outil.nom)
        b.setToolTip(outil.resume)
        b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        b.setAutoRaise(True)
        b.setCursor(Qt.CursorShape.PointingHandCursor)
        b.setMinimumWidth(self.COTE + 40)
        b.clicked.connect(lambda _=False, o=outil: self._lancer(o))
        return b

    def _lancer(self, outil):
        if outil.lancer():
            self.statusBar().showMessage("%s : ouverture en cours…" % outil.nom, 8000)
        else:
            self.statusBar().showMessage("%s ne s'est pas lance." % outil.nom)
