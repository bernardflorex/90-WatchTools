"""Theme commun a tous les outils : celui du programme de came.

Un seul endroit pour les couleurs. Chaque outil appelle `appliquer(app)` une
fois, avant de construire sa fenetre ; les traces matplotlib sont regles sur
les memes teintes, de sorte qu'un graphique matplotlib et un graphique dessine
a la main (`came.ui.Plot`) se lisent de la meme facon.
"""

from contextlib import contextmanager

from PySide6.QtGui import QColor, QPalette

BG = QColor(26, 28, 32)
FG = QColor(226, 230, 236)
GRID = QColor(58, 63, 72)
ACC = QColor(96, 180, 255)
ACC2 = QColor(255, 176, 80)
ACC3 = QColor(120, 220, 160)
ACC4 = QColor(240, 120, 150)

#: fond des champs de saisie et des listes
CHAMP = QColor(34, 38, 44)
#: fond et bord des boutons
BOUTON = QColor(43, 48, 58)
BORD = QColor(69, 75, 86)

#: couleurs d'etat, lisibles sur le fond sombre
ERREUR = "#ff7b7b"
SUCCES = "#78dca0"
#: bouton d'un calcul en cours
ACTIF = "background:#8a6a14;font-weight:bold"

_SAISIE = ("QDoubleSpinBox,QSpinBox,QLineEdit,QComboBox,QPlainTextEdit,"
           "QTextEdit,QAbstractItemView")

FEUILLE = (
    "QWidget{background:#1a1c20;color:#e2e6ec}"
    "QGroupBox{border:1px solid #3a3f48;margin-top:8px;padding-top:10px}"
    "QGroupBox::title{subcontrol-origin:margin;left:8px}"
    "QPushButton{background:#2b303a;border:1px solid #454b56;padding:5px}"
    "QPushButton:hover{background:#353b47}"
    "QPushButton:disabled,QLabel:disabled,QCheckBox:disabled{color:#7a808a}"
    + _SAISIE + "{background:#22262c;border:1px solid #3a3f48}"
    "QDoubleSpinBox:disabled,QSpinBox:disabled,QLineEdit:disabled,"
    "QComboBox:disabled{color:#7a808a}"
    "QAbstractItemView{selection-background-color:#2f5f8a;"
    "selection-color:#ffffff;gridline-color:#3a3f48}"
    "QHeaderView::section{background:#2b303a;border:1px solid #3a3f48;"
    "padding:3px}"
    "QMenuBar::item:selected,QMenu::item:selected{background:#2f5f8a}"
    "QMenu{border:1px solid #454b56}"
    "QToolTip{background:#2b303a;color:#e2e6ec;border:1px solid #454b56}"
    "QToolBar{border:0}"
    "QToolButton:hover,QToolButton:checked{background:#353b47}"
)


def _hex(c):
    return c.name()


def _palette():
    """Palette sombre : la feuille de style ne suffit pas.

    La barre d'outils de matplotlib choisit ses icones claires ou sombres
    d'apres la **palette**, pas d'apres la feuille de style : sans elle, les
    icones restaient noires sur fond noir.
    """
    p = QPalette()
    R = QPalette.ColorRole
    for role, c in ((R.Window, BG), (R.WindowText, FG), (R.Base, CHAMP),
                    (R.AlternateBase, BG), (R.Text, FG), (R.Button, BOUTON),
                    (R.ButtonText, FG), (R.ToolTipBase, BOUTON),
                    (R.ToolTipText, FG), (R.PlaceholderText, QColor(122, 128, 138)),
                    (R.Highlight, QColor(47, 95, 138)),
                    (R.HighlightedText, QColor(255, 255, 255)),
                    (R.Link, ACC), (R.Light, BORD), (R.Midlight, GRID),
                    (R.Mid, GRID), (R.Dark, QColor(16, 17, 20)),
                    (R.Shadow, QColor(0, 0, 0))):
        p.setColor(role, c)
    gris = QColor(122, 128, 138)
    for role in (R.WindowText, R.Text, R.ButtonText):
        p.setColor(QPalette.ColorGroup.Disabled, role, gris)
    return p


def _rc_sombre():
    from cycler import cycler
    fond, trait, grille = _hex(BG), _hex(FG), _hex(GRID)
    return {
        "figure.facecolor": fond, "figure.edgecolor": fond,
        "savefig.facecolor": fond, "savefig.edgecolor": fond,
        "axes.facecolor": fond, "axes.edgecolor": grille,
        "axes.labelcolor": trait, "axes.titlecolor": trait,
        "text.color": trait, "xtick.color": trait, "ytick.color": trait,
        "xtick.labelcolor": trait, "ytick.labelcolor": trait,
        "grid.color": grille,
        "legend.facecolor": _hex(CHAMP), "legend.edgecolor": grille,
        "legend.labelcolor": trait,
        "patch.edgecolor": trait, "hatch.color": trait,
        "boxplot.boxprops.color": trait,
        "axes.prop_cycle": cycler(color=[_hex(c) for c in (ACC, ACC2, ACC3, ACC4)]
                                  + ["#b39ddb", "#80cbc4", "#e6ee9c", "#bcaaa4"]),
    }


def appliquer_matplotlib():
    """Regle matplotlib sur le theme. Sans effet si matplotlib est absent."""
    try:
        import matplotlib
    except ImportError:
        return
    matplotlib.rcParams.update(_rc_sombre())


def appliquer(app):
    """Applique le theme a l'application Qt et a matplotlib."""
    app.setStyle("Fusion")
    app.setPalette(_palette())
    app.setStyleSheet(FEUILLE)
    appliquer_matplotlib()


@contextmanager
def clair():
    """Figures sur fond blanc, le temps d'un export (rapport, GIF).

    Un rapport s'imprime : ses figures gardent le fond blanc de matplotlib,
    quel que soit le theme de l'interface. S'emploie aussi en decorateur.
    """
    import matplotlib
    with matplotlib.rc_context():
        matplotlib.rcParams.update(
            {k: matplotlib.rcParamsDefault[k] for k in _rc_sombre()})
        yield


def encre():
    """Couleur de trait principale d'un trace : le « noir » du theme courant."""
    import matplotlib
    return matplotlib.rcParams["text.color"]


def fond():
    """Couleur de fond d'un trace, pour un cartouche pose sur le dessin."""
    import matplotlib
    return matplotlib.rcParams["axes.facecolor"]


def gris(niveau):
    """Gris matplotlib de niveau 0 (noir) a 1 (blanc), **sur fond blanc**.

    Le niveau est replace entre le fond et l'encre du theme courant, pour
    garder le meme contraste avec le fond : sur le theme sombre un trait
    discret le reste, un trait marque aussi.
    """
    from matplotlib.colors import to_hex, to_rgb
    a, b = to_rgb(fond()), to_rgb(encre())
    return to_hex([x + (1.0 - niveau) * (y - x) for x, y in zip(a, b)])
