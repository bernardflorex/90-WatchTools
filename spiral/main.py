# main.py
"""
Point d'entrée. Import absolu (pas `from .app import ...`) + ajout du
répertoire PARENT de spiral_lab/ à sys.path : permet de lancer ce fichier
directement (VSCode "Run Python File", `python main.py` depuis spiral_lab/,
etc.) sans avoir besoin de `python -m spiral_lab.main` depuis le dossier
parent -- un import relatif échoue ("attempted relative import with no
known parent package") dès que ce fichier est exécuté comme script isolé
plutôt qu'importé comme membre du package.
"""
import os
import sys

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PySide6.QtWidgets import QApplication
from spiral_lab.app import LaboSpiralApp


def main():
    app = QApplication(sys.argv)
    fenetre = LaboSpiralApp()
    fenetre.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
