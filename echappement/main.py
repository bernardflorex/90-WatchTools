"""Outil d'échappement à ancre suisse.

    python -m echappement.main [tracé.dxf]      (depuis la racine WatchTools)
    python echappement/main.py [tracé.dxf]      (lancé directement, VSCode...)

Imports absolus et répertoire PARENT ajouté à sys.path quand le fichier est
exécuté comme script : sans cela, ni `commun` ni l'import relatif de l'app
ne sont trouvés.
"""

import os
import sys

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PySide6.QtWidgets import QApplication

from commun import theme
from echappement.app import FenetreEchappement


def main():
    app = QApplication(sys.argv)
    theme.appliquer(app)
    fenetre = FenetreEchappement(sys.argv[1] if len(sys.argv) > 1 else None)
    fenetre.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
