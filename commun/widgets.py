"""Widgets de panneau partages par tous les outils.

Ils viennent du programme de came, ou ils ont ete mis au point : un panneau
de reglages etroit, defilant, dont les textes se replient sans imposer leur
largeur ni laisser de vide.
"""

from PySide6.QtCore import Qt, QSize, QTimer
from PySide6.QtWidgets import QLabel, QScrollArea, QSizePolicy, QTabWidget


class Etiquette(QLabel):
    """Etiquette repliee qui ne dicte pas la largeur du panneau.

    Avec un texte replie et une politique horizontale `Ignored`, Qt calcule la
    hauteur **mini** pour une largeur minimale : le resume du bas reclamait
    280 px alors qu'il n'en occupe que 110. Cette contradiction faisait
    distribuer a Qt des centaines de pixels orphelins. On annule donc ce
    minimum ; la hauteur utile est celle de `heightForWidth`, que le layout
    sait calculer.
    """

    def minimumSizeHint(self):
        return QSize(0, 0)


def soft(lbl):
    """Empeche un QLabel de dicter la largeur du panneau.

    Avec un texte replie et une politique horizontale `Ignored`, Qt calcule la
    hauteur **mini** du label pour une largeur minimale : le resume du bas
    reclamait 280 px alors qu'il n'en occupe que 110. Cette contradiction (un
    minimum superieur a la place reellement prise) faisait distribuer a Qt des
    centaines de pixels orphelins dans le panneau, d'ou la barre d'onglets qui
    descendait et le resume detache du reste.

    On annule donc ce minimum : la hauteur est celle de `heightForWidth` a la
    largeur reellement attribuee, que le layout sait deja calculer.
    """
    lbl.setWordWrap(True)
    lbl.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Minimum)
    lbl.setMinimumSize(0, 0)
    lbl.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
    if not isinstance(lbl, Etiquette):
        # Poser `minimumSizeHint` sur l'instance suffisait avec PyQt6 ; PySide6
        # n'appelle pas une methode virtuelle ainsi attachee. Il faut donc une
        # sous-classe - et `Etiquette` en est une, deja traitee.
        lbl.setMinimumHeight(0)
    return lbl




def aide(widget, texte):
    """Texte long : en infobulle, pas dans le panneau.

    Les explications detaillees gonflaient la colonne de gauche et Qt leur
    reservait plus de hauteur qu'elles n'en occupaient : le panneau se
    retrouvait avec plus de mille pixels de vide.
    """
    widget.setToolTip("<div style='width:360px'>" + texte + "</div>")
    return widget


class OngletsCompacts(QTabWidget):
    """QTabWidget qui ne reserve que la hauteur de l'onglet affiche.

    Par defaut le QStackedWidget interne annonce le maximum des hauteurs de
    toutes les pages : l'onglet le plus charge imposait sa hauteur au panneau,
    d'ou un grand vide sous les onglets courts.

    Les deux indications de taille renvoient **la meme** hauteur. Les pages
    ont en effet un `minimumSizeHint` plus grand que leur `sizeHint` (661
    contre 474 px pour le chronographe : les textes replies imposent une
    hauteur mini que la largeur ignoree ne reduit pas). Annoncer un maximum
    inferieur au minimum laissait Qt distribuer 505 px orphelins autour du
    bloc d'onglets, ce qui faisait descendre la barre d'onglets de plus de
    200 px au changement d'onglet et detachait d'autant le resume du bas.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Ignored,
                           QSizePolicy.Policy.Fixed)
        self.currentChanged.connect(self._reajuster)

    def _reajuster(self, _i=0):
        # `currentChanged` arrive avant que la nouvelle page soit disposee :
        # a cet instant la hauteur annoncee est encore celle de l'ancienne, et
        # le layout parent la gardait en cache (985 px de l'onglet
        # d'optimisation au lieu des 511 du chronographe), d'ou un bloc
        # d'onglets pousse 237 px plus bas. On invalide donc une seconde fois,
        # apres la disposition.
        self.updateGeometry()
        QTimer.singleShot(0, self._invalider)

    def _invalider(self):
        # Invalider le layout du parent ne suffit pas : QWidget garde **son
        # propre** cache de sizeHint. Sans `updateGeometry` sur l'ancetre, son
        # layout annoncait 1138 px alors que le widget en annoncait encore
        # 1612 (l'onglet precedent), et la zone defilante lui donnait 1612 :
        # les 474 px d'excedent se repartissaient en deux vides de 237 px.
        self.updateGeometry()
        par = self.parentWidget()
        while par is not None:
            lay = par.layout()
            if lay is not None:
                lay.invalidate()
                lay.activate()
            par.updateGeometry()
            if isinstance(par, QScrollArea):
                break
            par = par.parentWidget()

    def _haut(self):
        pg = self.currentWidget()
        if pg is None:
            return super().sizeHint().height()
        h = max(pg.sizeHint().height(), pg.minimumSizeHint().height())
        m = self.contentsMargins()
        return h + self.tabBar().sizeHint().height() + m.top() + m.bottom() + 8

    def sizeHint(self):
        return QSize(super().sizeHint().width(), self._haut())

    def minimumSizeHint(self):
        return QSize(super().minimumSizeHint().width(), self._haut())

    # La zone defilante dimensionne le panneau par `heightForWidth` des qu'un
    # texte y est replie. Or un QTabWidget le calcule sur **toutes** ses pages
    # et renvoyait donc la hauteur de l'onglet d'optimisation : 474 px de trop
    # sur l'onglet du chronographe, repartis en deux vides de 237 px.
    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, _w):
        return self._haut()
