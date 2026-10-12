"""WatchTools : accueil des outils de calcul horloger.

    python main.py              ouvre l'accueil : une icone par outil
    python main.py engrenage    lance directement un outil

Chaque outil reste un programme a part (`python -m <paquet>.main`). Pour en
ajouter un : une ligne dans `OUTILS` et une fonction de dessin dans
`commun/icones.py`.
"""

import argparse
import runpy
import sys

from PySide6.QtWidgets import QApplication

from commun import icones, theme
from commun.lanceur import Lanceur, Outil

OUTILS = (
    Outil("came", "Came", "Dimensionnement d'un coeur de chronographe et de son marteau",
          "came.main", icones.came, theme.ACC),
    Outil("engrenage", "Engrenage", "Profils de denture, simulation d'engrènement, "
          "optimisation, fraise", "engrenage.main", icones.engrenage, theme.ACC2),
    Outil("spiral", "Spiral", "Spiral : courbes terminales, isochronisme, optimisation",
          "spiral.main", icones.spiral, theme.ACC3),
)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outil", nargs="?", choices=[o.cle for o in OUTILS],
                    help="outil a lancer sans passer par l'accueil")
    args, reste = ap.parse_known_args(argv)
    if args.outil:
        module = next(o.module for o in OUTILS if o.cle == args.outil)
        sys.argv = [module] + reste
        runpy.run_module(module, run_name="__main__")
        return 0
    app = QApplication(sys.argv[:1])
    theme.appliquer(app)
    fenetre = Lanceur(OUTILS)
    fenetre.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
