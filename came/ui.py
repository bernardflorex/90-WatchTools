"""Interface graphique de dimensionnement du coeur (PySide6)."""

import math
import os

from PySide6.QtCore import (Qt, QPointF, QRectF, QSize, QThread, QTimer,
                            Signal)
from PySide6.QtGui import (QBrush, QColor, QFont, QPainter, QPainterPath, QPen,
                         QPolygonF)
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDoubleSpinBox,
                             QFileDialog, QFormLayout, QFrame, QGridLayout,
                             QGroupBox, QHBoxLayout, QLabel, QListWidget,
                             QMainWindow, QMessageBox, QPlainTextEdit,
                             QPushButton, QScrollArea, QSizePolicy, QSlider,
                             QSpinBox, QSplitter, QTabWidget, QVBoxLayout,
                             QWidget)

from commun import theme
from commun.theme import ACC, ACC2, ACC3, ACC4, BG, FG, GRID
from commun.widgets import Etiquette, OngletsCompacts, aide, soft
from . import arcs, dxfio, projet
from .chrono import ChronoParams, curves, pose_at, synthesize_chrono
from .chrono import _X
from .nsga2 import OBJECTIFS, nsga2, pool_lot, probleme
from .profil import (Profil, ProfilImporte, depuis_points,
                     dimensionner_ressort, evaluer, fenetre_pointe,
                     balayage_marteau, from_chrono, longueur_plat,
                     plan_recommande, secteur_faible, secteur_marche)
from .noyau import DISPONIBLE as ACCEL
from .profil import pose as pose_profil
from .core import Params, norm, pose, synthesize
from .optimize import VARIABLES, Weights, cost, optimize

def wrap(txt, n=46):
    """Coupe un texte en lignes de n caracteres (HTML)."""
    out, line = [], ""
    for word in str(txt).split():
        if len(line) + len(word) + 1 > n:
            out.append(line)
            line = word
        else:
            line = (line + " " + word).strip()
    if line:
        out.append(line)
    return "<br>".join(out)


def _ecart_contour(pts, pr):
    """Ecart maxi entre le contour importe et le profil reconstruit.

    Le marteau a une face plane : il ne voit que l'enveloppe convexe du
    dessin. Un contour non convexe ou tronque s'en ecarte, et c'est cet
    ecart qu'il faut regarder avant de se fier au resultat.
    """
    import math as _m
    ec = 0.0
    n = 720
    for i in range(n):
        b = 2.0 * _m.pi * i / n
        cb, sb = _m.cos(b), _m.sin(b)
        h = max(cb * q[0] + sb * q[1] for q in pts)
        ec = max(ec, abs(h - pr.H(b)))
    return ec


def avec_dxf(path):
    """Ajoute l'extension .dxf si l'utilisateur ne l'a pas saisie."""
    if not path:
        return path
    return path if path.lower().endswith(".dxf") else path + ".dxf"


#: Ce que represente chaque valeur a renseigner, affiche en infobulle quand
#: la souris passe sur le champ ou sur son libelle.
AIDES = {
    "pivot_x":
        "Abscisse de l'axe de rotation du marteau, dans le repere du coeur "
        "dont l'origine est son propre axe de rotation. Le marteau ne fait "
        "que tourner autour de ce point : il ne se translate jamais.",
    "pivot_y":
        "Ordonnee de l'axe de rotation du marteau, meme repere. La distance "
        "O->P fixe le bras de levier du marteau et donc le couple qu'il faut "
        "pour obtenir le couple vise sur le coeur.",
    "r_min":
        "Rayon du coeur au plat de zero, c'est-a-dire la distance de l'axe de "
        "rotation au plat qui vient contre le marteau en position stable. "
        "C'est lui qui fixe le bras de levier au repos.",
    "r_max":
        "Encombrement maximal admissible : rayon de la pointe du coeur. "
        "L'ecart r_max - r_min est la course radiale que le marteau doit "
        "absorber en tournant autour de son pivot.",
    "couple_min":
        "Couple du ressort de marteau au repos, soit son prearmage. Dans le "
        "sens inverse du calcul (onglet d'optimisation) ce n'est plus une "
        "entree : le prearmage devient un **resultat**.",
    "pente":
        "Raideur du ressort de marteau, en couple par radian d'ouverture. Le "
        "couple vaut C0 + k.psi. Comme le prearmage, c'est un resultat dans "
        "l'onglet d'optimisation.",
    "mu":
        "Coefficient de frottement dynamique acier sur acier au contact "
        "marteau / coeur. Il rend le couple different a l'aller et au retour "
        "(frottement rentrant et sortant) et fixe le plancher "
        "d'arc-boutement 2.atan(mu) a la pointe.",
    "angle_plan":
        "Pente du plan de contact du marteau, mesuree **de la direction O->P "
        "a la normale du plan**, au repos. Elle est donc independante de "
        "l'armage du coeur, contrairement a un angle rapporte a "
        "l'horizontale ou a la verticale.",
    "rayon_plan":
        "Rayon du plan de contact s'il est bombe au lieu d'etre droit. "
        "Zero signifie un plan plat, qui est le cas usuel.",
    "longueur_plan":
        "Longueur de la **partie droite** du plan d'impulsion, arrondis de "
        "bout non compris. Zero signifie que le plan est suppose assez long "
        "pour rester tangent au coeur en toute position. Un plan court est "
        "souhaitable : il met l'arrondi en contact a la pointe et divise la "
        "plage morte par plus de cinq.",
    "rayon_coin":
        "Rayon des deux arrondis qui terminent le plan d'impulsion. C'est "
        "cet arrondi qui porte le contact des que le coeur sort du plan ; au "
        "passage de la pointe, il empeche le bras de levier de s'annuler. Un "
        "petit rayon reduit la plage morte mais augmente la pression de "
        "contact.",
    "decalage_plan":
        "Position du milieu du plan le long de lui-meme, comptee depuis le "
        "pied de la normale issue de l'axe du coeur quand le systeme est au "
        "repos (plat du coeur contre plat du marteau). Elle decide lequel des "
        "deux arrondis verra passer la pointe.",
    "sweep":
        "Amplitude de rotation du coeur que la synthese analytique couvre. "
        "Un tour complet n'est pas atteignable d'un seul flanc.",
    "n_lobes":
        "Nombre de lobes du profil. Un coeur de remise a zero simple en a un.",
    "n_steps":
        "Nombre de pas d'integration de la synthese. Plus il est grand, plus "
        "le profil est lisse, et plus le calcul est long.",
    "couple_coeur":
        "Couple de rappel voulu sur l'axe du coeur. C'est la consigne du "
        "calcul inverse : la geometrie et le ressort en decoulent.",
    "corde_pointe":
        "Longueur du petit segment qui coiffe la pointe du coeur. Il reduit "
        "la plage morte au sommet en donnant un bras de levier non nul.",
    "rayon_pointe":
        "Rayon d'usinage de la pointe du coeur. Une pointe vive n'est pas "
        "fabricable ; sur ce rayon le couple ne doit pas s'annuler.",
    "rayon_creux":
        "Rayon du creux entre les deux flancs, au droit de la pointe.",
    "rayon_raccord":
        "Rayon de raccordement entre les flancs et le plat du zero. Dans ces "
        "raccords, la resultante de l'effort ne doit pas passer par l'axe de "
        "rotation, sinon l'aiguille ne revient pas tout a fait a zero.",
    "plat_souhaite":
        "Longueur de plat souhaitee, a titre de controle. La longueur reelle "
        "est imposee par la geometrie : environ 2(r_max - r_min)/course.",
    "marge_blocage":
        "Facteur de securite sur la condition de non-blocage du plat "
        "(longueur superieure a mu.r_min).",
    "cellules": "Nombre d'harmoniques qui decrivent la courbure du profil. "
                "Plus il est grand, plus le profil peut se deformer, et plus "
                "l'optimisation est longue.",
    "population": "Nombre de profils evalues a chaque generation.",
    "generations": "Nombre de generations de l'algorithme NSGA-II.",
    "scan": "Nombre de positions du coeur analysees sur un tour.",
    "procs": "Nombre de processus d'evaluation lances en parallele.",
}


def simuler_plan_fini(params, contour, n=360, echelle=1.0):
    """Rejoue une geometrie synthetisee avec le **plan reel** du marteau.

    Les deux syntheses analytiques cherchent la tangence sans borne
    d'abscisse : elles supposent le plan assez long. Des que l'utilisateur
    renseigne une longueur de plan, un rayon de coin ou un decalage, la
    cinematique qu'elles calculent ne correspond plus a la piece, et le
    marteau dessine quitte le coeur (jusqu'a 0.97 mm d'ecart sur 50 positions
    sur 61 pour un plan de 1.4 mm).

    On refait donc tourner le contour obtenu dans le moteur de contact a plan
    fini, celui qui gere la tangence bornee et les deux arrondis de bout.
    C'est cette simulation qui sert ensuite au trace, au point de contact et
    aux courbes.

    Renvoie `(profil, evaluation)` ou `None` si le plan est suppose infini ou
    si la simulation echoue.
    """
    if params.longueur_plan <= 1e-9 or not contour or len(contour) < 16:
        return None
    try:
        q = params.copy(couple_min=params.couple_min * echelle,
                        pente=params.pente * echelle)
        pr = ProfilImporte([(float(a), float(b)) for a, b in contour])
        e = evaluer(q, pr, n=n, rayon_face=0.0, mu_contact=None)
        if e is not None and e.ok:
            e.marche_deg, e.marche = secteur_marche(pr, e)
    except Exception:
        return None
    return (pr, e) if (e is not None and e.ok) else None


AIDES_OPTIM = {
    "cellules": "Nombre d'harmoniques qui decrivent la courbure du profil. "
                "Plus il est grand, plus le profil peut se deformer, et plus "
                "l'optimisation est longue.",
    "population": "Nombre de profils evalues a chaque generation.",
    "generations": "Nombre de generations de l'algorithme NSGA-II.",
    "scan": "Nombre de positions du coeur analysees sur un tour.",
    "procs": "Nombre de processus d'evaluation lances en parallele.",
}


class _Faux:
    """Case a cocher absente : repond toujours non."""

    def isChecked(self):
        return False

    def setChecked(self, _v=False):
        pass

    def toggled(self, *_a):
        pass


#: phase 1 : ce qui definit le coeur, le pan du marteau etant suppose infini.
#: Le pivot et la pente du pan y entrent : le marteau tourne autour de son
#: pivot, et les bras de levier en dependent.
COEUR = ("mu", "r_min", "r_max", "couple_coeur",
         "pivot_x", "pivot_y", "angle_plan", "n_steps")

#: phase 2 : le pan reel, pose sur le plat obtenu en phase 1
MARTEAU = ("longueur_plan", "decalage_plan", "rayon_coin")

#: champs qui changent la **forme du coeur** : il faut refaire la synthese
FORME_COEUR = ("mu", "r_min", "r_max", "couple_coeur", "pivot_x", "pivot_y",
               "angle_plan", "n_steps", "corde_pointe", "rayon_pointe",
               "rayon_creux", "rayon_raccord", "plat_souhaite",
               "marge_blocage")

#: ancien nom, conserve pour la lecture des etudes
SHARED = COEUR + MARTEAU

#: parametres qui ne decrivent que le plan du marteau : ils ne changent pas la
#: geometrie du coeur, donc la synthese reste valable et seule la simulation a
#: plan fini doit etre refaite
PLAN_MARTEAU = ("longueur_plan", "rayon_coin", "decalage_plan")

#: parametres qui deplacent la normale du plan au repos, gamma0 : le contour
#: importe doit etre repose plat sur plan
GAMMA0 = ("angle_plan", "pivot_x", "pivot_y")

#: tout ce qui change la simulation du profil affiche sur l'onglet
#: d'optimisation. Sans cette liste, changer la longueur du plat ou le
#: frottement ne mettait rien a jour, alors que changer l'angle, lui, relancait
#: le calcul : l'interface reagissait a certains champs et pas a d'autres.
SIMULATION = GAMMA0 + PLAN_MARTEAU + ("mu", "couple_coeur")


# ------------------------------------------------------------------ vue camme

class CamView(QWidget):
    """Vue cinematique : coeur, marteau, contact, efforts."""

    def __init__(self):
        super().__init__()
        self.setMinimumSize(360, 320)
        self.params = ChronoParams()
        self.res = None
        self.cres = None
        self.prof = None
        self.pev = None
        self.importe = []
        self.ccurves = None
        self.mode = "simple"
        self.x_deg = 0.0
        self.index = 0
        self.scale = 12.0
        self.center = (0.0, 0.0)
        self.offset = QPointF(0, 0)
        self._user_zoom = False
        self.show_trace = True
        self.show_forces = True
        self._drag = None
        self._trace = []
        # simulation en plan fini de la geometrie synthetisee, le cas echeant
        self.prof_fini = None
        self.ev_fini = None
        self.arcs = None            # chaine d'arcs a superposer

    def set_data(self, params, res):
        self.mode = "simple"
        self.params, self.res = params, res
        self._trace = []
        self.prof_fini = self.ev_fini = None
        if res is not None and res.ok:
            self._trace = [(q[0], q[1]) for q in res.profile]
            fini = simuler_plan_fini(params, res.closed_profile(params.n_lobes))
            if fini:
                self.prof_fini, self.ev_fini = fini
        self.fit()

    def set_arcs(self, approx):
        self.arcs = approx
        self.update()

    def set_profil(self, params, pr, ev, importe=None):
        self.mode = "profil"
        self.params, self.prof, self.pev = params, pr, ev
        self.importe = list(importe) if importe else []
        self._trace = pr.contour(360) if pr is not None else []
        self.fit()

    def set_chrono(self, params, cres):
        self.mode = "chrono"
        self.params, self.cres = params, cres
        self.ccurves = (curves(params, cres)
                        if cres is not None and cres.ok else None)
        self._trace = list(cres.profile) if cres is not None else []
        self.prof_fini = self.ev_fini = None
        if cres is not None and cres.ok:
            fini = simuler_plan_fini(params, cres.profile,
                                     echelle=cres.echelle)
            if fini:
                self.prof_fini, self.ev_fini = fini
        self.fit()

    def fit(self):
        p = self.params
        m = 0.25 * p.r_max + 1.0
        xs = [-p.r_max - m, p.r_max + m, p.pivot_x - m, p.pivot_x + m]
        ys = [-p.r_max - m, p.r_max + m, p.pivot_y - m, p.pivot_y + m]
        bw = max(max(xs) - min(xs), 1e-6)
        bh = max(max(ys) - min(ys), 1e-6)
        self.center = (0.5 * (max(xs) + min(xs)), 0.5 * (max(ys) + min(ys)))
        self.scale = 0.92 * min(max(self.width(), 50) / bw,
                                max(self.height(), 50) / bh)
        self.offset = QPointF(0, 0)
        self._user_zoom = False
        self.update()

    def resizeEvent(self, e):
        if not self._user_zoom:
            self.fit()
        super().resizeEvent(e)

    # -- interactions
    #: bornes du zoom, en pixels par millimetre
    ZOOM_MINI, ZOOM_MAXI = 0.5, 20000.0

    def wheelEvent(self, e):
        """Zoom centre sur le curseur.

        La transformation `_m` a pour point fixe `center`, qui se projette au
        milieu de la fenetre decale de `offset`. Ne multiplier que `scale`
        prenait donc ce milieu comme origine du zoom - un point qui, apres un
        deplacement de la vue, n'est pas forcement sur la partie affichee : la
        vue zoomee se decalait completement.

        On garde a la place le point du dessin qui se trouve **sous le
        curseur**. En notant `a` la position de la souris comptee depuis le
        milieu de la fenetre et `k` le facteur de zoom, la condition
        `_m(point sous le curseur)` inchangee donne, sur chaque composante :

            offset' = a - (a - offset) . k

        C'est exact, et symetrique en x et y malgre l'axe des ordonnees
        inverse du trace.
        """
        dy = e.angleDelta().y()
        if not dy:
            return                         # molette horizontale : on ignore
        k = 1.0015 ** dy
        neuf = max(self.ZOOM_MINI, min(self.ZOOM_MAXI, self.scale * k))
        k = neuf / self.scale              # facteur reellement applique
        if abs(k - 1.0) < 1e-12:
            return
        m = e.position()
        ax = m.x() - self.width() / 2.0
        ay = m.y() - self.height() / 2.0
        self.offset = QPointF(ax - (ax - self.offset.x()) * k,
                              ay - (ay - self.offset.y()) * k)
        self.scale = neuf
        self._user_zoom = True
        self.update()

    def mousePressEvent(self, e):
        self._drag = e.position()

    def mouseMoveEvent(self, e):
        if self._drag is not None:
            d = e.position() - self._drag
            self.offset += d
            self._user_zoom = True
            self._drag = e.position()
            self.update()

    def mouseReleaseEvent(self, e):
        self._drag = None

    def mouseDoubleClickEvent(self, e):
        self.fit()

    # -- rendu
    def _m(self, p):
        cx = self.width() / 2 + self.offset.x()
        cy = self.height() / 2 + self.offset.y()
        return QPointF(cx + (p[0] - self.center[0]) * self.scale,
                       cy - (p[1] - self.center[1]) * self.scale)

    def paintEvent(self, _e):
        q = QPainter(self)
        q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), BG)
        p = self.params

        # grille
        q.setPen(QPen(GRID, 1, Qt.PenStyle.DotLine))
        step = 5.0
        n = int(max(self.width(), self.height()) / (step * self.scale)) + 2
        for i in range(-n, n + 1):
            a = self._m((i * step, -n * step))
            b = self._m((i * step, n * step))
            q.drawLine(a, b)
            a = self._m((-n * step, i * step))
            b = self._m((n * step, i * step))
            q.drawLine(a, b)

        # axes et cercles de contrainte
        q.setPen(QPen(GRID.lighter(140), 1))
        q.drawLine(self._m((-1e3, 0)), self._m((1e3, 0)))
        q.drawLine(self._m((0, -1e3)), self._m((0, 1e3)))
        q.setPen(QPen(QColor(150, 150, 160), 1, Qt.PenStyle.DashLine))
        for r in (p.r_min, p.r_max):
            d = r * self.scale
            c = self._m((0, 0))
            q.drawEllipse(QRectF(c.x() - d, c.y() - d, 2 * d, 2 * d))

        # pivot du marteau
        pv = self._m(p.pivot)
        q.setPen(QPen(ACC2, 1.5))
        q.drawLine(self._m((0, 0)), pv)
        q.setBrush(QBrush(ACC2))
        q.drawEllipse(pv, 4, 4)
        q.setBrush(Qt.BrushStyle.NoBrush)
        q.setPen(QPen(FG, 1))
        q.drawText(pv + QPointF(8, -6), "pivot marteau")

        if self.mode == "profil":
            if self.prof is None or self.pev is None or not self.pev.ok:
                q.setPen(QPen(FG))
                q.drawText(14, 22, "Pas de profil - lancer l'optimisation")
                return
            ps = pose_profil(p, self.prof, self.pev, self.x_deg)
            if ps is None:
                return
            ps["Cm"] = p.couple_min + p.pente * abs(ps["phi"])
        elif self.prof_fini is not None and self.ev_fini is not None:
            # Une longueur de plan a ete renseignee : la cinematique de la
            # synthese (tangence d'un plan infini) ne vaut plus. On trace la
            # simulation a plan fini, sinon le marteau quitte le coeur.
            ech = (self.cres.echelle
                   if self.mode == "chrono" and self.cres is not None else 1.0)
            ps = pose_profil(p, self.prof_fini, self.ev_fini, self.x_deg)
            if ps is None:
                return
            ps["Cm"] = (p.couple_min + p.pente * abs(ps["phi"])) * ech
        elif False:
            if self.res2 is None or not self.res2.ok:
                q.setPen(QPen(FG))
                q.drawText(14, 22, "Pas de solution - ajuster les parametres")
                return
            ps = pose2(p, self.res2, self.x_deg)
            ps["Cm"] = ((p.couple_min + p.pente * abs(ps["phi"]))
                        * self.res2.echelle)
            ps["face"] = ps["faces"][0]
        elif self.mode == "chrono":
            if self.cres is None or not self.cres.ok:
                q.setPen(QPen(FG))
                q.drawText(14, 22, "Pas de solution - ajuster les parametres")
                return
            ps = pose_at(p, self.cres, self.x_deg)
            ps["N"], ps["Th"] = 0.0, 0.0
            if self.ccurves:
                xs, Th, br, N, _ps = self.ccurves
                i = min(range(len(xs)), key=lambda k: abs(xs[k] - self.x_deg))
                ps["N"], ps["Th"] = N[i], Th[i]
            ps["Cm"] = ((p.couple_min + p.pente * abs(ps["phi"]))
                        * self.cres.echelle)
        else:
            if self.res is None or not self.res.ok:
                q.setPen(QPen(FG))
                q.drawText(14, 22, "Pas de solution - ajuster les parametres")
                return
            ps = pose(p, self.res, self.index)

        # trace du profil au repos (reference)
        if self.show_trace and self._trace:
            q.setPen(QPen(QColor(90, 100, 115), 1, Qt.PenStyle.DotLine))
            q.drawPolyline(QPolygonF([self._m(t) for t in self._trace]))

        # coeur
        poly = QPolygonF([self._m(t) for t in ps["profile"]])
        path = QPainterPath()
        path.addPolygon(poly)
        path.closeSubpath()
        q.setBrush(QBrush(QColor(40, 70, 110)))
        q.setPen(QPen(ACC, 2))
        q.drawPath(path)
        q.setBrush(Qt.BrushStyle.NoBrush)
        # contour lu dans le DXF, tourne comme le coeur : il montre ce qui a
        # ete importe, par-dessus la reconstruction
        if self.importe:
            ca, sa = math.cos(ps["theta"]), math.sin(ps["theta"])
            pp = [self._m((v[0] * ca - v[1] * sa, v[0] * sa + v[1] * ca))
                  for v in self.importe]
            q.setPen(QPen(QColor(255, 120, 150), 2, Qt.PenStyle.DashLine))
            q.drawPolyline(QPolygonF(pp + [pp[0]]))

        # chaine d'arcs superposee : on voit ou l'approximation derive
        if self.arcs is not None and self.arcs.arcs:
            ca, sa = math.cos(ps["theta"]), math.sin(ps["theta"])
            q.setPen(QPen(QColor(120, 190, 255), 1, Qt.PenStyle.DashLine))
            for A in self.arcs.arcs:
                pp = [self._m((v[0] * ca - v[1] * sa, v[0] * sa + v[1] * ca))
                      for v in A.echantillonner(20)]
                q.drawPolyline(QPolygonF(pp))
            q.setBrush(QBrush(QColor(120, 190, 255)))
            for A in self.arcs.arcs:
                v = A.depart
                q.drawEllipse(self._m((v[0] * ca - v[1] * sa,
                                       v[0] * sa + v[1] * ca)), 2, 2)
            q.setBrush(Qt.BrushStyle.NoBrush)

        c0 = self._m((0, 0))
        q.setPen(QPen(FG, 1))
        q.drawEllipse(c0, 3, 3)

        if self.mode == "chrono" and self.cres is not None and self.cres.flat:
            ca, sa = math.cos(ps["theta"]), math.sin(ps["theta"])
            fl = [self._m((v[0] * ca - v[1] * sa, v[0] * sa + v[1] * ca))
                  for v in self.cres.flat]
            q.setPen(QPen(ACC3, 3))
            q.drawLine(fl[0], fl[1])

        # marteau : bras + plan de contact
        ct = self._m(ps["contact"])
        # Sur la marche qui ramene la pointe au plat, le coeur est concave :
        # ni le plan ni son arrondi ne peuvent y toucher, le marteau franchit
        # la corde de l'enveloppe convexe. On le dit, plutot que de laisser
        # croire a un defaut de calcul.
        saut = ps.get("franchit", False)
        col = QColor(255, 170, 80) if saut else ACC2
        q.setPen(QPen(col, 4))
        q.drawPolyline(QPolygonF([self._m(t) for t in ps["face"]]))
        q.setPen(QPen(col, 2))
        q.drawLine(pv, self._m(ps["face"][0]))
        q.drawLine(pv, self._m(ps["face"][-1]))
        if saut:
            q.setPen(QPen(QColor(255, 170, 80), 1))
            q.drawText(14, 22, "sans contact : le marteau franchit la marche "
                               "pointe -> plat (concavite)")

        # contact + efforts
        q.setBrush(QBrush(ACC4))
        q.setPen(QPen(ACC4, 1))
        q.drawEllipse(ct, 4, 4)
        q.setBrush(Qt.BrushStyle.NoBrush)
        if self.show_forces and ps["N"] > 0.0:
            if self.mode == "profil":
                nmax = (self.pev.N_max if self.pev else 1.0) or 1.0
            elif self.mode == "simple":
                nmax = max(self.res.N) or 1.0
            else:
                nmax = (self.cres.N_max if self.cres else 1.0) or 1.0
            L = 0.45 * p.r_max * (ps["N"] / nmax) + 0.15 * p.r_max
            nx, ny = ps["normal"]
            tip = (ps["contact"][0] - L * nx, ps["contact"][1] - L * ny)
            q.setPen(QPen(ACC4, 2))
            q.drawLine(self._m((ps["contact"][0] + L * nx,
                                ps["contact"][1] + L * ny)), self._m(tip))
            q.setPen(QPen(FG, 1))
            q.drawText(self._m(tip) + QPointF(6, -4), "N = %.3f N" % ps["N"])

        q.setPen(QPen(FG, 1))
        f = QFont(); f.setPointSize(9); q.setFont(f)
        txt = ["coeur  %+7.2f deg" % math.degrees(ps["theta"]),
               "marteau %+6.2f deg" % math.degrees(ps["phi"]),
               "H = %.3f mm" % ps["H"],
               "Cm = %.3f N.mm" % ps["Cm"],
               "couple coeur = %+.4f N.mm" % ps["Th"]]
        for i, t in enumerate(txt):
            q.drawText(14, 20 + 15 * i, t)


# ------------------------------------------------------------------ graphiques

class Plot(QWidget):
    """Trace de courbes y(x) multiples, sans dependance externe."""

    pointClicked = Signal(int)

    def __init__(self, title="", xlabel="", ylabel=""):
        super().__init__()
        self.setMinimumHeight(150)
        self.title, self.xlabel, self.ylabel = title, xlabel, ylabel
        self.curves = []        # (nom, xs, ys, couleur)
        self.points = []        # (x, y) d'un nuage cliquable
        self.selected = -1
        self.marker = None
        self._tr = None

    def set_points(self, pts, selected=-1):
        self.points = list(pts)
        self.selected = selected
        self.update()

    def mousePressEvent(self, ev):
        if not self.points or self._tr is None:
            return
        mx, my, L, T, w, h, x0, x1, y0, y1 = self._tr
        px, py = ev.position().x(), ev.position().y()
        best, bd = -1, 1e18
        for i, (x, y) in enumerate(self.points):
            dx = (L + (x - x0) / (x1 - x0) * w) - px
            dy = (T + h - (y - y0) / (y1 - y0) * h) - py
            d = dx * dx + dy * dy
            if d < bd:
                best, bd = i, d
        if best >= 0 and bd < 400:
            self.selected = best
            self.update()
            self.pointClicked.emit(best)

    def set_curves(self, curves):
        self.curves = curves
        self.update()

    def set_marker(self, x):
        self.marker = x
        self.update()

    def paintEvent(self, _e):
        q = QPainter(self)
        q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), BG)
        L, R, T, B = 58, 150, 20, 28
        w = self.width() - L - R
        h = self.height() - T - B
        if w < 20 or h < 20:
            return
        pts = [(x, y) for _n, xs, ys, _c in self.curves
               for x, y in zip(xs, ys) if y == y]
        pts += [q for q in self.points if q[1] == q[1]]
        q.setPen(QPen(GRID, 1))
        q.drawRect(L, T, w, h)
        f = QFont(); f.setPointSize(8); q.setFont(f)
        q.setPen(QPen(FG, 1))
        q.drawText(L, T - 6, self.title)
        if not pts:
            q.drawText(L + 8, T + 18, "pas de donnees")
            return
        x0 = min(p[0] for p in pts); x1 = max(p[0] for p in pts)
        y0 = min(p[1] for p in pts); y1 = max(p[1] for p in pts)
        if x1 - x0 < 1e-12:
            x0, x1 = x0 - 0.5, x1 + 0.5
        elif self.points:
            padx = 0.08 * (x1 - x0)
            x0, x1 = x0 - padx, x1 + padx
        pad = 0.08 * (y1 - y0) if y1 - y0 > 1e-12 else max(1e-9, abs(y1) * 0.1 + 1e-9)
        y0 -= pad; y1 += pad

        def mx(x):
            return L + (x - x0) / (x1 - x0) * w

        def my(y):
            return T + h - (y - y0) / (y1 - y0) * h

        # grille + graduations
        q.setPen(QPen(GRID, 1, Qt.PenStyle.DotLine))
        for i in range(1, 5):
            yy = T + h * i / 5.0
            q.drawLine(L, int(yy), L + w, int(yy))
        for i in range(1, 6):
            xx = L + w * i / 6.0
            q.drawLine(int(xx), T, int(xx), T + h)
        q.setPen(QPen(FG, 1))
        for i in range(6):
            v = y0 + (y1 - y0) * i / 5.0
            q.drawText(2, int(my(v)) + 4, "%8.4g" % v)
        span = x1 - x0
        fmt = "%.0f" if span >= 20 else ("%.1f" if span >= 2 else "%.3g")
        for i in range(7):
            v = x0 + span * i / 6.0
            q.drawText(int(mx(v)) - 16, T + h + 16, fmt % v)
        q.drawText(L + w // 2 - 20, self.height() - 2, self.xlabel)
        if y0 < 0 < y1:
            q.setPen(QPen(GRID.lighter(160), 1))
            q.drawLine(L, int(my(0)), L + w, int(my(0)))

        for i, (name, xs, ys, col) in enumerate(self.curves):
            q.setPen(QPen(col, 1.8))
            poly = QPolygonF([QPointF(mx(x), my(y))
                              for x, y in zip(xs, ys) if y == y])
            q.drawPolyline(poly)
            q.setPen(QPen(col, 1))
            q.drawLine(L + w + 10, T + 12 + 15 * i, L + w + 28, T + 12 + 15 * i)
            q.setPen(QPen(FG, 1))
            q.drawText(L + w + 34, T + 16 + 15 * i, name)

        self._tr = (0, 0, L, T, w, h, x0, x1, y0, y1)
        for i, (x, y) in enumerate(self.points):
            cx, cy = mx(x), my(y)
            if i == self.selected:
                q.setPen(QPen(QColor(255, 255, 255), 2))
                q.setBrush(QBrush(ACC2))
                q.drawEllipse(QPointF(cx, cy), 6, 6)
            else:
                q.setPen(QPen(ACC, 1))
                q.setBrush(QBrush(ACC.darker(130)))
                q.drawEllipse(QPointF(cx, cy), 4, 4)
        q.setBrush(Qt.BrushStyle.NoBrush)

        if self.marker is not None and x0 <= self.marker <= x1:
            q.setPen(QPen(QColor(255, 255, 255, 110), 1, Qt.PenStyle.DashLine))
            q.drawLine(int(mx(self.marker)), T, int(mx(self.marker)), T + h)


# ------------------------------------------------------------------ thread opti

class OptiThread(QThread):
    done = Signal(object, object, float, list)
    step = Signal(int, float)

    def __init__(self, p, w, vars_, mode, bounds=None):
        super().__init__()
        self.p, self.w, self.vars, self.mode = p, w, vars_, mode
        self.bounds = bounds

    def run(self):
        try:
            q, r, c, log = optimize(self.p, self.w, self.vars, mode=self.mode,
                                    restarts=2, bounds=self.bounds,
                                    progress=lambda n, c: self.step.emit(n, c))
        except Exception as exc:          # pragma: no cover
            self.done.emit(None, None, 0.0, ["erreur : %r" % exc])
            return
        self.done.emit(q, r, c, log)


# ------------------------------------------------------------------ fenetre

class ChronoThread(QThread):
    done = Signal(object)
    step = Signal(str)

    def __init__(self, p, equilibrer=True):
        super().__init__()
        self.p, self.equilibrer = p, equilibrer

    def run(self):
        try:
            r = synthesize_chrono(self.p, self.equilibrer,
                                  progress=lambda m: self.step.emit(m))
        except Exception as exc:          # pragma: no cover
            r = None
            self.message = repr(exc)
        self.done.emit(r)


class NsgaThread(QThread):
    done = Signal(object, object)
    step = Signal(int, int, int)

    def __init__(self, p, cellules, n_pop, n_gen, n_scan, objectifs,
                 face_libre=False, procs=None, galet=False, figer=False,
                 marteau=(),
                 graine=None):
        super().__init__()
        self.p, self.cellules = p, cellules
        self.n_pop, self.n_gen, self.n_scan = n_pop, n_gen, n_scan
        self.objectifs, self.face_libre = objectifs, face_libre
        self.procs, self.galet, self.figer = procs, galet, figer
        self.marteau = tuple(marteau)
        self.graine_init = graine
        self.graine_f4 = None
        self.gabarit = (0.0, 0.0)
        self.graine_viol = 0.0
        self.ressort_fige = None
        self.stop = False

    def run(self):
        try:
            from .chrono import synthesize_chrono
            r = synthesize_chrono(self.p)
            pr0 = (self.graine_init if self.graine_init is not None
                   else from_chrono(self.p, r, self.cellules))
            graine = pr0.vector()
            if self.face_libre:
                graine = graine + [0.0]      # on part de la face plane
            if self.marteau:
                # la graine porte les valeurs de depart du marteau en queue
                graine = graine + [0.0 if n == "pivot_t"
                                   else float(getattr(self.p, n))
                                   for n in self.marteau]
            if self.figer:
                # etape 2 : on garde le ressort de la version initiale et on
                # fait varier coeur et marteau a armement identique
                r0 = dimensionner_ressort(self.p, pr0, 0.0, None,
                                          n=self.n_scan)
                self.ressort_fige = (r0.prearmage, r0.raideur)
            # le plat de la graine est la reference : l'optimisation ne doit
            # pas le dissoudre, sans quoi le coeur perd sa position nulle
            plat_ref = 0.9 * longueur_plat(pr0, self.p.gamma0, 0.002,
                                           grille=300)
            self.gabarit = getattr(pr0, "gabarit", (0.0, 0.0))
            ev, bornes, rep, gfaces = probleme(self.p, n_scan=self.n_scan,
                                               face_libre=self.face_libre,
                                               galet=self.galet,
                                               ressort=self.ressort_fige,
                                               plat_mini=plat_ref,
                                               gabarit=self.gabarit,
                                               marteau=self.marteau)
            sel = self.objectifs

            def ev_sel(x):
                f, v = ev(x)
                return tuple(f[i] for i in sel), v

            lo, hi = bornes(graine)
            lot, fermer = (None, lambda: None)
            if self.procs and self.procs > 1:
                try:
                    lot, fermer = pool_lot(
                        self.p, self.n_scan, self.face_libre, self.procs,
                        self.galet, self.ressort_fige, plat_ref,
                        self.gabarit, self.marteau)
                    if lot is not None:
                        base = lot

                        def lot(xs):
                            out = base(xs)
                            return [(x, tuple(f[i] for i in sel), v)
                                    for x, f, v in out]
                except Exception:
                    lot, fermer = None, (lambda: None)
            pop, front = nsga2(
                ev_sel, lo, hi, n_pop=self.n_pop, n_gen=self.n_gen,
                graines=gfaces(graine), reparer_x=rep, evaluer_lot=lot,
                progress=lambda g, n, nf, _p: self.step.emit(g, n, nf))
            fermer()
            # on recalcule les 4 objectifs pour l'affichage
            for q in front:
                q["f4"], q["viol"] = ev(q["x"])
            self.graine_f4, self.graine_viol = ev(graine)
            self.done.emit(front, graine)
        except Exception as exc:          # pragma: no cover
            self.message = repr(exc)
            self.done.emit(None, None)


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Dimensionnement de cames horlogeres - coeur de remise a zero")
        self.params = ChronoParams()
        self.weights = Weights()
        self.res = None
        self.cres = None
        self.res2 = None
        self.front = None
        self.profil_importe = None
        self.prof_dessine = None
        self.points_importes = None
        self.approx = None          # chaine d'arcs du profil affiche
        self._ressort_saisi = None
        self._tmr_coeur = None
        self._vue = "chrono"
        self.ressort = None
        self.prof = None
        self.pev = None
        self.nsga = None
        self.chrono = None
        self.spins = {}
        self.aides_champs = {}
        self.labels = {}
        self.wspins = {}
        self.varboxes = {}
        self.opti = None
        self._build()
        self._update_rules()
        self.recompute()

    # -- construction -------------------------------------------------------
    def _build(self):
        central = QSplitter(Qt.Orientation.Horizontal)
        central.addWidget(self._panel())

        right = QSplitter(Qt.Orientation.Vertical)
        self.view = CamView()
        right.addWidget(self._view_block())
        tabs = QTabWidget()
        self.plot_torque = Plot("Couple sur l'axe du coeur", "angle du coeur [deg]")
        self.plot_eff = Plot("Efforts et couple marteau", "angle du coeur [deg]")
        self.plot_geo = Plot("Geometrie du profil", "angle du coeur [deg]")
        self.plot_pareto = Plot("Front de Pareto", "", "")
        self.plot_pareto.pointClicked.connect(self._on_pareto_click)
        tabs.addTab(self.plot_torque, "Couple")
        tabs.addTab(self.plot_eff, "Efforts")
        tabs.addTab(self.plot_geo, "Geometrie")
        tabs.addTab(self.plot_pareto, "Pareto")
        self.tabs_plot = tabs
        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        self.report.setFont(QFont("monospace", 9))
        tabs.addTab(self.report, "Rapport")
        right.addWidget(tabs)
        right.setSizes([520, 300])
        central.addWidget(right)
        central.setStretchFactor(0, 0)
        central.setStretchFactor(1, 1)
        central.setSizes([360, 980])
        self.setCentralWidget(central)
        self.statusBar().showMessage("pret")
        self.resize(1340, 860)

    def _view_block(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(self.view, 1)
        bar = QHBoxLayout()
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, 100)
        self.slider.valueChanged.connect(self._on_slider)
        self.btn_play = QPushButton("Animer")
        self.btn_play.setCheckable(True)
        self.btn_play.toggled.connect(self._on_play)
        self.lbl_pos = QLabel("0.0 deg")
        self.lbl_pos.setMinimumWidth(70)
        cb1 = QCheckBox("profil au repos"); cb1.setChecked(True)
        cb1.toggled.connect(lambda b: (setattr(self.view, "show_trace", b),
                                       self.view.update()))
        cb2 = QCheckBox("efforts"); cb2.setChecked(True)
        cb2.toggled.connect(lambda b: (setattr(self.view, "show_forces", b),
                                       self.view.update()))
        bar.addWidget(self.btn_play)
        bar.addWidget(self.slider, 1)
        bar.addWidget(self.lbl_pos)
        bar.addWidget(cb1)
        bar.addWidget(cb2)
        v.addLayout(bar)
        self.timer = QTimer(self)
        self.timer.setInterval(30)
        self.timer.timeout.connect(self._tick)
        self._dir = 1
        return w

    def _spin(self, name, label, unit, lo, hi, step):
        sb = QDoubleSpinBox()
        sb.setRange(lo, hi)
        sb.setSingleStep(step)
        sb.setDecimals(0 if name in ("n_lobes", "n_steps") else 3)
        sb.setValue(float(getattr(self.params, name)))
        sb.setKeyboardTracking(False)
        sb.setProperty("pname", name)
        sb.valueChanged.connect(self._on_param)
        # Au survol : ce que la valeur represente, plus sa plage admise.
        txt = AIDES.get(name)
        bulle = "<b>%s</b> [%s]" % (label, unit)
        if txt:
            bulle += "<br>" + txt
        bulle += "<br><i>de %g a %g, pas de %g</i>" % (lo, hi, step)
        sb.setToolTip(bulle)
        self.aides_champs[name] = bulle
        self.spins[name] = sb
        return sb

    def _ligne(self, lay, name, label, unit, lo, hi, step):
        """Ajoute une ligne libelle / champ, les deux porteurs de l'infobulle.

        Le libelle la porte aussi : l'utilisateur survole souvent le texte
        plutot que la case.
        """
        sb = self._spin(name, label, unit, lo, hi, step)
        lbl = QLabel("%s [%s]" % (label, unit))
        lbl.setToolTip(sb.toolTip())
        self.labels[name] = lbl
        lay.addRow(lbl, sb)
        return sb

    def _panel(self):
        host = QWidget()
        lay = QVBoxLayout(host)
        specs = {n: (n, l, u, a, b, st)
                 for n, l, u, a, b, st in Params.SPECS + ChronoParams.EXTRA}

        # Deux phases, dans l'ordre ou l'on dimensionne :
        #   1. le coeur, avec un pan de marteau suppose infini ;
        #   2. le pan reel, pose sur le plat obtenu.
        g1 = QGroupBox("1. Coeur \u2014 pan de marteau infini")
        f1 = QFormLayout(g1)
        for name in COEUR:
            _n, label, unit, lo, hi, step = specs[name]
            self._ligne(f1, name, label, unit, lo, hi, step)
        lay.addWidget(g1)

        g1b = QGroupBox("2. Pan reel du marteau")
        f1b = QFormLayout(g1b)
        for name in MARTEAU:
            _n, label, unit, lo, hi, step = specs[name]
            self._ligne(f1b, name, label, unit, lo, hi, step)
        aide(g1b, "Le couple du marteau vient de la phase 1 : il est calcule "
                  "pour que le coeur rende le couple voulu avec un pan "
                  "infini. Ces champs posent ensuite le pan <b>reel</b>, et "
                  "l'on voit le couple qu'il delivre vraiment.")
        lay.addWidget(g1b)

        # OngletsCompacts fixe lui-meme sa politique de taille
        # Un seul panneau. Les trois onglets - coeur simple, chronographe,
        # optimisation - obligeaient a deviner lequel decrivait la piece, et
        # certains champs n'agissaient que sur l'un d'eux. L'application traite
        # maintenant le coeur de chronographe, d'un bout a l'autre, dans une
        # seule page : on construit le coeur avec un pan de marteau suppose
        # infini, puis on pose le pan reel et on voit ce qu'il delivre.
        host_mode = QWidget()
        lmode = QVBoxLayout(host_mode)
        lmode.setContentsMargins(0, 0, 0, 0)
        lmode.addWidget(self._tab_chrono(specs))
        lmode.addWidget(self._tab_optim(specs))
        self.tabs_mode = host_mode
        lay.addWidget(host_mode)

        g4 = QGroupBox("Fichiers")
        f4 = QGridLayout(g4)
        f4.setContentsMargins(6, 4, 6, 4)
        b_imp = QPushButton("Importer DXF")
        b_imp.setToolTip("Lit le coeur, le pivot et le plan du marteau dans "
                         "un DXF. Le coeur est ramene a sa position de repos, "
                         "plat sur plan.")
        b_imp.clicked.connect(self._on_import)
        b_exp = QPushButton("Exporter DXF")
        b_exp.setToolTip("Exporte ce que montre l'onglet actif.")
        b_exp.clicked.connect(self._on_export)
        b_ouv = QPushButton("Ouvrir etude")
        b_ouv.setToolTip("Relit une etude .json : tous les parametres et le "
                         "contour du coeur d'un seul coup.")
        b_ouv.clicked.connect(self._on_ouvrir)
        b_sav = QPushButton("Enregistrer etude")
        b_sav.setToolTip("Ecrit un .json avec les parametres, le contour du "
                         "coeur, le ressort calcule et les reglages.")
        b_sav.clicked.connect(self._on_sauver)
        self.cb_cine = QCheckBox("exporter les 12 positions du marteau")
        self.cb_cine.setToolTip(
            "Decoche : un seul marteau, a la position de repos, plat sur "
            "plat. Coche : douze positions reparties sur la course, pour "
            "voir la cinematique relative.")
        f4.addWidget(b_imp, 0, 0)
        f4.addWidget(b_exp, 0, 1)
        f4.addWidget(b_ouv, 1, 0)
        f4.addWidget(b_sav, 1, 1)
        f4.addWidget(self.cb_cine, 2, 0, 1, 2)
        g4.setSizePolicy(QSizePolicy.Policy.Ignored,
                         QSizePolicy.Policy.Fixed)

        g5 = QGroupBox("Profil en arcs de cercle")
        f5 = QFormLayout(g5)
        f5.setContentsMargins(6, 4, 6, 4)
        self.sp_arc = QDoubleSpinBox()
        self.sp_arc.setRange(0.0005, 0.5)
        self.sp_arc.setDecimals(4)
        self.sp_arc.setSingleStep(0.005)
        self.sp_arc.setValue(0.010)
        self.sp_arc.setToolTip(
            "<b>Ecart maxi admissible</b> [mm]<br>Ecart entre le profil et "
            "son approximation par arcs tangents, mesure le long de la "
            "normale. C'est vous qui le fixez : le <b>nombre d'arcs</b> en "
            "decoule et s'affiche dans le rapport. Un nombre d'arcs impose ne "
            "dirait rien de la precision obtenue.")
        self.sp_arc.valueChanged.connect(self._maj_arcs)
        lbl_arc = QLabel("Ecart maxi [mm]")
        lbl_arc.setToolTip(self.sp_arc.toolTip())
        f5.addRow(lbl_arc, self.sp_arc)
        self.cb_arcs = QCheckBox("montrer le profil en arcs")
        self.cb_arcs.setToolTip(
            "Superpose la chaine d'arcs au contour, pour voir la derive de "
            "l'approximation. Les arcs sont aussi ecrits dans le DXF.")
        self.cb_arcs.toggled.connect(self._maj_arcs)
        f5.addRow(self.cb_arcs)
        g5.setSizePolicy(QSizePolicy.Policy.Ignored,
                         QSizePolicy.Policy.Fixed)
        lay.addWidget(g5)
        # en tete du panneau : le contenu des onglets depasse la hauteur de
        # l'ecran, et en bas ces deux boutons devenaient inatteignables
        lay.insertWidget(0, g4)

        self.lbl_ressort = soft(Etiquette())
        self.lbl_ressort.setFrameShape(QFrame.Shape.StyledPanel)
        aide(self.lbl_ressort,
             "Le ressort est ce qu'il reste a determiner : il est calcule "
             "pour que le couple de rappel voulu soit atteint. Ce n'est donc "
             "pas une donnee a saisir.")
        lay.addWidget(self.lbl_ressort)

        self.lbl_sum = soft(Etiquette())
        self.lbl_sum.setFrameShape(QFrame.Shape.StyledPanel)
        lay.addWidget(self.lbl_sum)
        lay.addStretch(1)

        host.setSizePolicy(QSizePolicy.Policy.Ignored,
                           QSizePolicy.Policy.Preferred)
        sc = QScrollArea()
        sc.setWidget(host)
        sc.setWidgetResizable(True)
        sc.setMinimumWidth(300)
        sc.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Expanding)
        sc.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.panel = sc
        return sc

    def _tab_optim(self, specs):
        w = QWidget()
        lay = QVBoxLayout(w)
        _role = soft(Etiquette('<b>Simule et optimise</b> un coeur existant, avec le plan reel.'))
        aide(_role, "Cet onglet ne construit pas de coeur : il prend celui qui existe deja - importe d'un DXF, relu d'une etude, ou issu d'une des deux syntheses - et c'est le <b>seul</b> a modeliser le plan du marteau tel qu'il est : longueur finie, arrondis de bout, decalage. C'est donc ici que les chiffres decrivent la piece. NSGA-II deforme ensuite le profil pour chercher mieux que la graine, a couple vise constant.")
        lay.addWidget(_role)

        g0 = QGroupBox("3. Optimisation du profil (second plan)")
        f0 = QVBoxLayout(g0)
        lb = soft(Etiquette(
            "<b>Donnees</b> : encombrement, pivot, pente, couple vise.<br>"
            "<b>Sorties</b> : geometrie du coeur, prearmage et raideur du "
            "ressort.<br>"
            "<span style='color:#ffb050'>C0 et k du panneau du haut sont "
            "grises : ils affichent le resultat.</span>"))
        aide(lb, "Le ressort est calcule pour que le couple vise soit "
                 "atteint : on ajuste Cm = C0 + k.psi sur le couple marteau "
                 "exige en chaque position. Les champs C0 et k du panneau du "
                 "haut ne sont pas lus dans cet onglet, ils affichent le "
                 "ressort de la solution choisie.")
        f0.addWidget(lb)
        lay.addWidget(g0)

        lb2 = soft(Etiquette(
            "<span style='color:#7ed9a0'>Noyau accelere : actif</span>"
            if ACCEL else
            "<span style='color:#ffb050'>Noyau en Python pur</span>"))
        aide(lb2, "numpy et numba compilent le noyau de calcul (facteur 68 "
                  "par individu). Le noyau est mono-thread : repartir la "
                  "population sur plusieurs processus occupe tous les "
                  "coeurs." if ACCEL else
                  "Installer numpy et numba accelere l'optimisation d'un "
                  "facteur 68 par individu.")
        lay.addWidget(lb2)

        g = QGroupBox("Reglages de l'optimisation")
        f = QFormLayout(g)
        self.sp_opt = {}
        for name, label, lo, hi, val in (
                ("cellules", "Harmoniques du profil", 6, 60, 24),
                ("population", "Population", 10, 400, 60),
                ("generations", "Generations", 1, 500, 35),
                ("scan", "Pas d'analyse (/tour)", 60, 720, 360),
                ("procs", "Processus d'evaluation", 1, 64,
                 max(1, (__import__("os").cpu_count() or 2) - 1))):
            sb = QSpinBox()
            sb.setRange(lo, hi)
            sb.setValue(val)
            sb.setToolTip("<b>%s</b><br>%s<br><i>de %g a %g</i>"
                          % (label, AIDES_OPTIM.get(name, ""), lo, hi))
            self.sp_opt[name] = sb
            lbl = QLabel(label)
            lbl.setToolTip(sb.toolTip())
            f.addRow(lbl, sb)
        lay.addWidget(g)

        g2 = QGroupBox("Objectifs de l'optimisation (au moins deux)")
        f2 = QVBoxLayout(g2)
        self.cb_obj = []
        for nom in OBJECTIFS:
            cb = QCheckBox(nom)
            cb.setChecked(True)
            self.cb_obj.append(cb)
            f2.addWidget(cb)
        aide(g2,
             "<b>Ce que cherche l'optimisation.</b><br>"
             "Le couple vise est atteint par construction, pour toutes les "
             "solutions : le ressort est recalcule pour chacune. On ne compare "
             "donc pas le couple, mais <i>ce qu'il coute</i>.<br><br>"
             "<b>Fixe</b> : encombrement, pivot, pente et longueur du plan, "
             "couple vise.<br>"
             "<b>Varie</b> : la courbure du coeur, harmonique par harmonique."
             "<br><b>Sortie</b> : le ressort (prearmage et raideur).<br><br>"
             "<b>tour sous la consigne</b> : part du tour ou le couple reste "
             "sous la cible. Les rampes des deux equilibres y mettent un "
             "plancher incompressible.<br>"
             "<b>plage morte</b> : secteur ou le mecanisme s'arc-boute. Le "
             "profil du coeur n'a presque pas de prise dessus : le levier est "
             "le <i>plan du marteau</i> - un plan court met son arrondi au "
             "contact de la pointe et fait tomber la plage de 14.8 a 1.5 deg - "
             "ou un galet.<br>"
             "<b>ondulation</b> : ecart du couple a sa moyenne sur la partie "
             "utile. C'est la l'essentiel de la prise qu'a le profil.<br>"
             "<b>effort de contact</b> : pression sur la pointe du marteau et "
             "sur le coeur, donc l'usure.<br><br>"
             "La <b>graine analytique</b> (l'onglet du chronographe) est dans "
             "la population : le front ne peut donc pas etre plus mauvais "
             "qu'elle, et elle y figure si rien ne la domine.<br><br>"
             "Contraintes toujours actives : convexite, fermeture du contour, "
             "encombrement (a 5 um pres), rayon mini, deux changements de "
             "signe du couple, rayon de raccord usinable, pas de changement "
             "de cote du marteau.")
        lay.addWidget(g2)

        self.cb_figer = QCheckBox("Figer le ressort sur la version initiale")
        self.cb_figer.setToolTip(
            "Etape 2 : on calcule d'abord le ressort du coeur analytique, "
            "puis on le garde fige et on fait varier coeur et marteau. "
            "Toutes les solutions sont alors comparees a armement identique, "
            "seul moyen de savoir laquelle est meilleure.")
        lay.addWidget(self.cb_figer)
        # Le galet est abandonne : il divisait la plage morte par quinze a
        # trente, mais ce n'est pas la piece etudiee. Un faux objet garde les
        # appels existants sans ajouter de reglage au panneau.
        self.cb_galet = _Faux()
        self.cb_face = QCheckBox("Optimiser aussi la courbure du pan")
        self.cb_face.setToolTip(
            "Ajoute la courbure de la face du marteau aux variables.")
        lay.addWidget(self.cb_face)

        g3b = QGroupBox("Grandeurs du marteau a optimiser")
        f3b = QVBoxLayout(g3b)
        f3b.setContentsMargins(6, 4, 6, 4)
        self.cb_mart = {}
        for nom, lib in (("longueur_plan", "longueur du pan"),
                         ("decalage_plan", "decalage le long du pan"),
                         ("rayon_coin", "rayons de coin"),
                         ("angle_plan", "pente du pan"),
                         ("pivot_t", "pivot, perpendiculairement a O->P"),
                         ("pivot_x", "pivot X (peu utile)"),
                         ("pivot_y", "pivot Y (peu utile)")):
            cb = QCheckBox(lib)
            self.cb_mart[nom] = cb
            f3b.addWidget(cb)
        aide(g3b,
             "<b>Sur quoi l'optimisation a prise.</b><br>"
             "Sans ces cases, elle ne fait varier que la courbure du coeur, "
             "harmonique par harmonique. Or le profil du coeur n'a presque "
             "pas de prise sur la <b>plage morte</b> ni sur l'<b>effort de "
             "contact</b> : les leviers sont la longueur du pan, ses "
             "arrondis, sa pente et la position du pivot. Les laisser fixes "
             "revenait a demander d'ameliorer ce qui ne pouvait pas "
             "bouger.<br><br>"
             "Chaque case ajoute une variable, donc agrandit l'espace de "
             "recherche : a budget egal la convergence est plus lente. "
             "L'ondulation y gagne tout de suite, les autres objectifs "
             "demandent plus de generations.<br><br>"
             "<b>Le pivot : une seule direction porte.</b> Sur +/- 4 mm, le "
             "deplacer <i>radialement</i> fait varier l'ondulation de 1.6 "
             "point (41.3 a 42.9 %) ; <i>perpendiculairement</i> a O-&gt;P, de "
             "68.8 points (13.5 a 82.3 %). La distance a l'axe n'apporte donc "
             "quasiment rien : c'est la position <b>angulaire</b> du pivot qui "
             "compte. D'ou la case « perpendiculairement a O-&gt;P », qui ne "
             "coute qu'une variable au lieu de deux dont une inutile. X et Y "
             "restent disponibles pour qui veut le carre complet.<br><br>"
             "Bornes : pan de 0.15 a 1.4 fois l'encombrement, pente a "
             "+/- 28 deg, pivot a +/- 30 % de sa distance a l'axe. La pente "
             "et le pivot deplacent la normale de repos : le plat du coeur "
             "suit, pour que le coeur reste plat contre plat au depart.")
        lay.addWidget(g3b)

        self.cb_crit = QComboBox()
        for lib, cle in (("ondulation du couple", "ondulation"),
                         ("effort de contact", "effort"),
                         ("plage morte", "plage"),
                         ("tour sous la consigne", "consigne")):
            self.cb_crit.addItem(lib, cle)
        f3c = QFormLayout()
        f3c.addRow("Critere du balayage", self.cb_crit)
        lay.addLayout(f3c)

        self.btn_balay = QPushButton("Balayer pivot x pente")
        self.btn_balay.setToolTip(
            "Balaye les deux grandeurs du marteau qui portent reellement : "
            "la position angulaire du pivot et la pente du pan, a coeur fixe. "
            "Quelques secondes, et plus efficace que NSGA-II sur ces deux "
            "nombres - il repartit son effort sur quatre objectifs et sur "
            "tout l'espace des formes de coeur.")
        self.btn_balay.clicked.connect(self._on_balayage)
        lay.addWidget(self.btn_balay)

        self.btn_nsga = QPushButton("Lancer NSGA-II")
        self.btn_nsga.clicked.connect(self._on_nsga)
        lay.addWidget(self.btn_nsga)

        g3 = QGroupBox("Front de Pareto")
        f3 = QVBoxLayout(g3)
        self.lst_front = QListWidget()
        self.lst_front.setMinimumHeight(170)
        self.lst_front.currentRowChanged.connect(self._on_front_row)
        self.lbl_graine = soft(Etiquette())
        lay_g = self.lbl_graine
        f3.addWidget(lay_g)
        aide(self.lst_front,
             "Cliquez une solution ici ou un point du nuage (onglet Pareto) "
             "pour afficher son profil, ses courbes et son ressort.")
        f3.addWidget(self.lst_front)
        lay.addWidget(g3)

        self.lbl_opt = soft(Etiquette())
        self.lbl_opt.setFrameShape(QFrame.Shape.StyledPanel)
        lay.addWidget(self.lbl_opt)
        lay.addStretch(1)
        return w

    def _tab_chrono(self, specs):
        w = QWidget()
        lay = QVBoxLayout(w)
        _role = soft(Etiquette('<b>Construit</b> le coeur a deux flancs du chronographe.'))
        aide(_role, "Meme synthese inverse que l'onglet precedent, mais pour un coeur a deux flancs, avec le plat du zero, les raccords et la pointe. Elle suppose aussi le plan du marteau assez long. C'est elle qui fournit la <b>graine</b> de l'optimisation.")
        lay.addWidget(_role)
        g = QGroupBox("1b. Forme de la pointe, du plat et des raccords")
        f = QFormLayout(g)
        for name, label, unit, lo, hi, step in ChronoParams.EXTRA:
            self._ligne(f, name, label, unit, lo, hi, step)
        self.btn_chrono = QPushButton("Calculer le coeur")
        aide(self.btn_chrono,
             "Les courses des deux flancs sont calculees : les deux sommets "
             "sur le diametre maxi, et couple identique dans les deux sens. "
             "Elles dependent donc du frottement.")
        self.btn_chrono.clicked.connect(self._on_chrono)
        f.addRow(self.btn_chrono)
        lay.addWidget(g)

        self.lbl_rules = soft(Etiquette())
        self.lbl_rules.setFrameShape(QFrame.Shape.StyledPanel)
        lay.addWidget(self.lbl_rules)
        lay.addStretch(1)
        return w

    # -- evenements ---------------------------------------------------------
    def _on_param(self, _v):
        sb = self.sender()
        name = sb.property("pname")
        val = sb.value()
        if name in ("n_lobes", "n_steps"):
            val = int(round(val))
        setattr(self.params, name, val)
        self._update_rules()
        # Tout champ agit, et l'on voit le resultat evoluer. La phase 1
        # reconstruit le coeur (quelques secondes, en tache de fond), la
        # phase 2 ne refait que la simulation. Une temporisation evite de
        # relancer la synthese a chaque cran de la molette.
        if name in FORME_COEUR:
            self._differer_coeur()
            return
        # gamma0 = angle(O->P) + angle_plan : ces trois champs deplacent la
        # normale du plan du marteau, donc la position de repos du coeur
        # importe. Il faut le reposer plat sur plan.
        if name in GAMMA0:
            self._realigner_profil()
        if self.is_optim:
            if name in SIMULATION and self.prof is not None:
                self._resimuler_profil()
            else:
                self.statusBar().showMessage(
                    "parametres modifies - relancer l'optimisation")
        elif False:
            self.res2 = None
            self.statusBar().showMessage(
                "parametres modifies - relancer le calcul")
        elif self.is_chrono:
            # Les parametres du plan du marteau ne changent pas la geometrie
            # du coeur : la synthese reste valable, seule la re-simulation a
            # plan fini est a refaire. Sans cela l'affichage gardait l'ancien
            # plan et les modifications semblaient ignorees.
            if name in PLAN_MARTEAU and self.cres is not None and self.cres.ok:
                self._show_chrono()
                self.statusBar().showMessage(
                    "plan du marteau modifie - simulation remise a jour")
            else:
                self.cres = None
                self.statusBar().showMessage(
                    "parametres modifies - relancer "
                    "\u00ab Calculer le coeur \u00bb")
        else:
            self.recompute()

    def _on_slider(self, v):
        if self.is_optim:
            if self.prof is None or self.pev is None:
                return
            x = 360.0 * v / 100.0
            self.view.x_deg = x
            self.view.update()
            self.lbl_pos.setText("%.1f deg" % x)
            for pl in (self.plot_torque, self.plot_eff, self.plot_geo):
                pl.set_marker(x)
            return
        if False:
            if self.res2 is None or not self.res2.ok:
                return
            x = -180.0 + 360.0 * v / 100.0
            self.view.x_deg = x
            self.view.update()
            self.lbl_pos.setText("%.1f deg" % x)
            for pl in (self.plot_torque, self.plot_eff, self.plot_geo):
                pl.set_marker(x)
            return
        if self.is_chrono:
            if self.cres is None or not self.cres.ok:
                return
            xs = self.view.ccurves[0] if self.view.ccurves else [0.0]
            x0, x1 = min(xs), max(xs)
            x = x0 + (x1 - x0) * v / 100.0
            self.view.x_deg = x
            self.view.update()
            self.lbl_pos.setText("%.1f deg" % x)
            for pl in (self.plot_torque, self.plot_eff, self.plot_geo):
                pl.set_marker(x)
            return
        if self.res is None or not self.res.ok:
            return
        n = len(self.res.theta) - 1
        i = int(round(v / 100.0 * n))
        self.view.index = i
        self.view.update()
        th = math.degrees(self.res.theta[i])
        self.lbl_pos.setText("%.1f deg" % th)
        self.plot_torque.set_marker(th)
        self.plot_eff.set_marker(th)
        self.plot_geo.set_marker(th)

    def _on_play(self, on):
        self.timer.start() if on else self.timer.stop()

    def _tick(self):
        v = self.slider.value() + 2 * self._dir
        if v >= 100:
            v, self._dir = 100, -1
        elif v <= 0:
            v, self._dir = 0, 1
        self.slider.setValue(v)

    @property
    def mode(self):
        return "retour"

    @property
    def is_chrono(self):
        """Toujours vrai : l'application ne traite que le coeur de chrono."""
        return not self.is_optim

    @property
    def is_optim(self):
        """Vrai quand c'est un contour discretise qui est **affiche**.

        A distinguer de son existence : des la construction du coeur, un
        profil discretise est disponible et simule - sans passer par
        NSGA-II - mais c'est la vue du chronographe qui reste a l'ecran,
        parce qu'elle montre le plat, la corde et le creux de pointe.
        """
        return self._vue == "profil"

    #: le ressort n'est une donnee que pour le coeur simple ; ailleurs il est
    #: calcule a partir du couple vise sur le coeur
    RESSORT = ("couple_min", "pente")

    def _maj_ressort(self, C0=None, k=None, calcule=True):
        """Affiche le ressort du marteau, qui est une **sortie**.

        Le prearmage et la raideur ne sont plus des champs a saisir : ils
        decoulent du couple de rappel voulu. Les ecrire dans les parametres
        serait d'ailleurs faux - la forme du coeur depend de C0, et la
        synthese deriverait a chaque relance.
        """
        if C0 is None or k is None:
            self.lbl_ressort.setText(
                "<b>Ressort du marteau</b><br>"
                "<span style='color:#9aa'>en attente du calcul</span>")
            return
        self.lbl_ressort.setText(
            "<b>Ressort du marteau</b> (sortie du calcul)<br>"
            "prearmage C0 = <b>%.4f</b> N.mm<br>"
            "raideur k = <b>%.4f</b> N.mm/rad (%.4f N.mm/deg)"
            % (C0, k, k * math.pi / 180.0))

    def _differer_coeur(self):
        """Relance la synthese du coeur apres une courte pause.

        Chaque cran de molette ne doit pas lancer un calcul de plusieurs
        secondes : on attend que l'utilisateur ait fini de regler.
        """
        if self._tmr_coeur is None:
            self._tmr_coeur = QTimer(self)
            self._tmr_coeur.setSingleShot(True)
            self._tmr_coeur.timeout.connect(self._refaire_coeur)
        self.statusBar().showMessage("parametre modifie - calcul en attente")
        self._tmr_coeur.start(450)

    def _refaire_coeur(self):
        """Reconstruit le coeur, puis repose le pan reel dessus."""
        self.prof = self.prof_dessine = None
        self.points_importes = None
        self.cres = None
        self._on_chrono()

    def _on_mode(self, _i=0):
        """Affiche l'etat courant : un profil s'il en existe, sinon la synthese."""
        self.btn_play.setChecked(False)
        self._maj_ressort(calcule=True)
        if self.is_optim:
            if self.prof is not None:
                self._show_profil()
            else:
                self.statusBar().showMessage(
                    "lancer NSGA-II pour obtenir un front de solutions")
            return
        if self.is_chrono:
            if self.cres is None:
                self._on_chrono()
            else:
                self._show_chrono()
        else:
            self.recompute()

    def _collect_weights(self):
        for n, sb in self.wspins.items():
            setattr(self.weights, n, sb.value())

    def _on_optimize(self):
        if self.opti is not None and self.opti.isRunning():
            return
        self._collect_weights()
        vars_ = [n for n, cb in self.varboxes.items() if cb.isChecked()]
        if not vars_:
            QMessageBox.information(self, "Optimisation",
                                    "Selectionner au moins une variable.")
            return
        self.btn_opt.setEnabled(False)
        self.btn_opt.setText("Optimisation...")
        p = self.params.copy(n_steps=min(self.params.n_steps, 160))
        bnds = {n: (self.varbounds[n][0].value(), self.varbounds[n][1].value())
                for n in vars_}
        self.opti = OptiThread(p, self.weights, vars_, "retour", bnds)
        self.opti.step.connect(
            lambda n, c: self.statusBar().showMessage(
                "optimisation : %d evaluations, cout = %.5g" % (n, c)))
        self.opti.done.connect(self._on_opti_done)
        self.opti.start()

    def _on_opti_done(self, q, _r, _c, log):
        self.btn_opt.setEnabled(True)
        self.btn_opt.setText("Optimiser")
        if q is None:
            QMessageBox.warning(self, "Optimisation", "\n".join(log))
            return
        for n in VARIABLES:
            setattr(self.params, n, getattr(q, n))
            sb = self.spins.get(n)
            if sb:
                sb.blockSignals(True)
                sb.setValue(getattr(q, n))
                sb.blockSignals(False)
        self.recompute()
        self.report.appendPlainText("\n--- optimisation ---\n" + "\n".join(log))
        self.statusBar().showMessage("optimisation terminee")

    def _orientations(self, psi=0.0):
        """Rappel des orientations absolues du plan, pour la CAO."""
        p = self.params
        eps = 1.0 if p.angle_plan < 0 else -1.0
        return (
            "<b>Pente du marteau</b><br>" + wrap(
                "L'angle saisi est mesure depuis la direction O->P (centre du "
                "coeur vers le pivot) jusqu'a la **normale** du plan, et au "
                "repos. Rapporte a OP il ne depend pas de la position du "
                "pivot autour du coeur, et c'est lui qui fixe le bras de "
                "levier q = |OP|.sin(pente).")
            + "<br>" + wrap(
                "O->P est a %.2f deg de l'horizontale, donc au repos la "
                "normale est a %.2f deg et le plan lui-meme a %.2f deg."
                % (math.degrees(p.alpha_p), math.degrees(p.gamma0),
                   p.orientation_plan()))
            + "<br>" + wrap(
                "En cours de course le marteau tourne de phi : le plan "
                "tourne d'autant. C'est pour cela que la pente est definie "
                "au repos." + ("" if psi <= 0 else
                               " Ici la course est de %.2f deg, le plan va "
                               "donc de %.2f a %.2f deg."
                               % (math.degrees(psi), p.orientation_plan(),
                                  p.orientation_plan(eps * psi)))))

    def _update_rules(self):
        p = self.params
        est = (2.0 * (p.r_max - p.r_min) / math.radians(180.0)
               if p.r_max > p.r_min else 0.0)
        muc = p.mu
        # La ligne de la plage morte annoncait le **plancher theorique**
        # 2.atan(mu), valable pour une pointe vive en appui plan. On lisait
        # donc 13.7 deg alors que la simulation en donnait 1.05 : la valeur
        # ne changeait pas quand on raccourcissait le pan, parce que ce
        # n'etait pas la meme grandeur. Les deux sont desormais distinguees.
        e = self.pev
        mesure = ("<b>%.2f deg</b> (simulee)" % e.plage_morte
                  if e is not None and e.ok else "en attente du calcul")
        # Depart plat sur plan : c'est l'hypothese de tout le calcul, et
        # elle se perd des que la pente ou le pivot bougent sans que le plat
        # du coeur suive. Autant pouvoir le verifier d'un coup d'oeil.
        pr = self.prof
        depart = "en attente du calcul"
        if pr is not None and e is not None and e.ok and e.phi:
            ec = math.degrees(abs(e.phi[0]))
            if hasattr(pr, "beta0") and getattr(pr, "plat", 0.0) > 0.0:
                d = math.degrees(abs((p.gamma0 - pr.beta0 + math.pi)
                                     % (2.0 * math.pi) - math.pi))
                ec = max(ec, d)
            depart = ("<span style='color:#7ed9a0'>plat sur plan "
                      "(%.3f deg)</span>" % ec if ec < 0.2 else
                      "<span style='color:#f08'><b>decale de %.2f deg</b>"
                      "</span>" % ec)
        self.lbl_rules.setText(
            "<b>Reperes</b><br>"
            "depart de la simulation : %s<br>"
            "demi-plat mini (non-blocage) : %.4f mm<br>"
            "plat impose par la geometrie : ~%.2f mm<br>"
            "plage morte : %s<br>"
            "&nbsp;&nbsp;plancher d'une pointe vive en appui plan : "
            "2.atan(mu) = %.1f deg<br>"
            "plan au repos : %.2f deg / horizontale"
            % (depart, p.plat_mini, est, mesure,
               math.degrees(2.0 * math.atan(muc)), p.orientation_plan()))
        aide(self.lbl_rules,
             "Demi-plat mini = mu.r_min : en dessous, la resultante passe "
             "par l'axe au bord du plat et l'aiguille reste coincee.<br><br>"
             "La longueur du plat n'est pas reglable : elle vaut "
             "|s_A(0)| + s_B(0), de l'ordre de 2(r_max - r_min)/course."
             "<br><br>"
             "La plage morte au sommet vaut 2.atan(mu au contact) et ne "
             "depend presque que du frottement ; un galet monte fou la "
             "divise par 15 a 30.<br><br>"
             "La pente saisie est mesuree de la direction O->P a la "
             "<b>normale</b> du plan, et <b>au repos</b>. Rapportee a OP "
             "elle ne depend pas de la position du pivot. En cours de "
             "course le marteau tourne, donc l'orientation absolue du plan "
             "change : c'est pour cela qu'elle est definie au repos.")

    def _on_chrono(self):
        if self.chrono is not None and self.chrono.isRunning():
            return
        self.btn_chrono.setEnabled(False)
        self.btn_chrono.setText("Calcul en cours...")
        self.chrono = ChronoThread(self.params.copy(), True)
        self.chrono.step.connect(lambda m: self.statusBar().showMessage(m))
        self.chrono.done.connect(self._on_chrono_done)
        self.chrono.start()

    def _on_chrono_done(self, r):
        self.btn_chrono.setEnabled(True)
        self.btn_chrono.setText("Calculer le coeur")
        self.cres = r
        if r is None:
            msg = getattr(self.chrono, "message", "")
            self.statusBar().showMessage("calcul impossible : " + msg)
            self.report.setPlainText("Calcul impossible.\n" + msg)
            return
        self._show_chrono()

    def _profil_de_la_synthese(self):
        """Discretise le coeur construit, pour pouvoir le simuler tout de suite.

        Le profil discretise etait auparavant le produit d'un import ou d'une
        optimisation : il fallait donc lancer NSGA-II pour obtenir le couple
        reellement delivre, les contacts sur coin, la plage morte ou
        l'approximation par arcs d'un coeur qu'on venait de construire. Il est
        maintenant deduit de la synthese, avec son plat exact.
        """
        r = self.cres
        if r is None or not r.ok:
            return
        try:
            pr = from_chrono(self.params, r,
                             self.sp_opt["cellules"].value())
            self.ressort = dimensionner_ressort(
                self.params, pr, 0.0, None,
                n=self.sp_opt["scan"].value())
        except Exception as exc:
            self.statusBar().showMessage("simulation du coeur : %r" % exc)
            return
        self.prof = self.prof_dessine = pr
        self.points_importes = None
        self.pev = self.ressort.eval

    def _show_chrono(self):
        r = self.cres
        self._vue = "chrono"
        if r is not None and r.ok:
            self._profil_de_la_synthese()
            self._maj_ressort(self.params.couple_min * r.echelle,
                              self.params.pente * r.echelle, calcule=True)
            if self.prof is not None and self.pev is not None:
                # C'est la **simulation** qu'il faut montrer, pas la synthese.
                # Celle-ci donne un couple constant par construction - une
                # droite a la consigne - et sa plage morte est celle d'un pan
                # infini. Le trace, les courbes et le curseur passent donc sur
                # le profil discretise et son pan reel, ou theta = 0 est la
                # position de repos, plat contre plat, au debut de la course.
                self._vue = "profil"
                self._show_profil()
                self.report.setPlainText(
                    "\n".join(self._rapport_construction()))
                return
        self.view.set_chrono(self.params, r)
        self._chrono_plots()
        self._chrono_report()
        # Le curseur part a zero, ce qui sur cet onglet est le **debut de la
        # course**, soit la pointe du coeur : le marteau y apparaissait
        # toujours sur la pointe. On l'amene a la position de repos, x = 0,
        # ou le plat du coeur est contre le plan du marteau.
        if r is not None and r.ok and self.view.ccurves:
            xs = self.view.ccurves[0]
            x0, x1 = min(xs), max(xs)
            if x1 > x0:
                v0 = int(round(100.0 * (0.0 - x0) / (x1 - x0)))
                self.slider.blockSignals(True)
                self.slider.setValue(max(0, min(100, v0)))
                self.slider.blockSignals(False)
        self._on_slider(self.slider.value())
        if not r.ok:
            self.statusBar().showMessage("echec : " + r.message)
        elif r.valid:
            self.statusBar().showMessage(
                "ok - courses %.2f / %.2f deg, couple coeur %.4f N.mm"
                % (r.course_a, r.course_b, r.couple_coeur))
        else:
            self.statusBar().showMessage(" ; ".join(r.warnings)[:200])

    def _chrono_plots(self):
        r = self.cres
        if r is None or not r.ok:
            for pl in (self.plot_torque, self.plot_eff, self.plot_geo):
                pl.set_curves([])
            return
        x, Th, bras, N, psi = curves(self.params, r)
        mini = self.params.mu * self.params.r_min
        self.plot_torque.set_curves([
            ("|couple coeur| [N.mm]", x, Th, ACC),
        ])
        self.plot_eff.set_curves([
            ("bras resultante [mm]", x, bras, ACC3),
            ("mini mu.r_min", [min(x), max(x)], [mini, mini], ACC4),
        ])
        self.plot_geo.set_curves([
            ("N contact [N]", x, N, ACC4),
            ("ouverture marteau [deg]", x, psi, ACC2),
        ])
        self.plot_torque.title = "Couple sur l'axe du coeur (tour complet)"
        self.plot_eff.title = ("Bras de levier de la resultante "
                               "(0 = blocage)")
        self.plot_geo.title = "Efforts et course du marteau"

    def _rapport_construction(self):
        """Les lignes du rapport de construction, sans les afficher."""
        self._lignes_constr = None
        self._chrono_report(rendre=True)
        return self._lignes_constr or []

    def _chrono_report(self, rendre=False):
        r, p = self.cres, self.params
        sp = {n: (l, u) for n, l, u, *_ in Params.SPECS + ChronoParams.EXTRA}
        L = ["=== Coeur de chronographe ==="]
        for name in SHARED + tuple(n for n, *_ in ChronoParams.EXTRA):
            lb, un = sp[name]
            L.append("  %-26s %10.4f %s" % (lb, getattr(p, name), un))
        L.append("")
        if r is None or not r.ok:
            L.append("=== Echec ===")
            L.append("  " + (r.message if r else "pas de calcul"))
            self._lignes_constr = list(L)
            if not rendre:
                self.report.setPlainText("\n".join(L))
            return
        L.append("=== Courses (calculees) ===")
        L.append("  flanc A (sens direct)        %10.2f deg" % r.course_a)
        L.append("  flanc B (sens inverse)       %10.2f deg" % r.course_b)
        L.append("  somme                        %10.2f deg"
                 % (r.course_a + r.course_b))
        L.append("  balayage des sommets         %10.2f deg"
                 % r.ouverture_pointe)
        L += self._bloc_encombrement(p, r)
        L.append("")
        L.append("=== Couples (le marteau est moteur) ===")
        L.append("  couple sur le coeur          %10.5f N.mm" % r.couple_coeur)
        L.append("  dissymetrie des deux sens    %10.3f %%" % (100.0 * r.dissym))
        L.append("  couple marteau C0 requis     %10.4f N.mm"
                 % (p.couple_min * r.echelle))
        L.append("  pente marteau k requise      %10.4f N.mm/rad"
                 % (p.pente * r.echelle))
        L.append("  couple marteau maxi          %10.4f N.mm" % r.Cm_max)
        L.append("  (facteur d'echelle applique  %10.4f)" % r.echelle)
        L += self._bloc_plan(p, getattr(r, "demi_requis", 0.0))
        L += self._bloc_plan_fini()
        L.append("  course du marteau            %10.3f deg"
                 % math.degrees(r.psi_max))
        L.append("  effort normal maxi           %10.4f N" % r.N_max)
        L.append("")
        L.append("=== Plat du zero ===")
        L.append("  longueur imposee             %10.4f mm" % r.plat)
        L.append("    cote A                     %10.4f mm" % r.plat_a)
        L.append("    cote B                     %10.4f mm" % r.plat_b)
        L.append("  longueur souhaitee           %10.4f mm" % p.plat_souhaite)
        L.append("  demi-plat mini (mu.r_min)    %10.4f mm" % r.plat_mini)
        L.append("  bras de levier au bord       %10.4f mm" % r.bras_plat)
        L.append("  marge / mini                 %10.1f x"
                 % (min(r.plat_a, r.plat_b) / max(r.plat_mini, 1e-12)))
        L.append("  raccord plat / flanc : la tangente est commune, le raccord")
        L.append("    est donc naturel. Rayon de courbure du flanc au depart :")
        L.append("    cote A %.4f mm, cote B %.4f mm (raccord demande %.4f mm)"
                 % (r.A.curv[0], r.B.curv[0], p.rayon_raccord))
        L.append("")
        L.append("=== Pointe ===")
        L.append("  corde entre les deux sommets %10.4f mm" % r.corde)
        L.append("  rayon des sommets            %10.4f mm" % p.rayon_pointe)
        L.append("  rayon du creux de pointe     %10.4f mm" % p.rayon_creux)
        L.append("  plage morte au sommet        %10.2f deg" % r.zone_morte)
        L.append("    (balayage de normale       %10.2f deg)" % r.zone_morte_beta)
        L.append("    plancher theorique 2.atan(mu)%9.2f deg"
                 % math.degrees(2.0 * math.atan(p.mu)))
        L.append("    une pointe vive est optimale ; tout creux ou arrondi")
        L.append("    au sommet elargit la plage. Seul levier reel : mu.")
        if r.horns:
            L.append("  rayon des deux sommets       %10.5f / %.5f mm"
                     % (norm(r.horns[0]), norm(r.horns[1])))
        L.append("")
        L.append("=== Diagnostic ===")
        if r.valid:
            L.append("  profil exploitable")
        for wm in r.warnings:
            L.append("  ! " + wm)
        L += self._bloc_simulation()
        L.append("")
        L.append("=== Profil (repere du coeur) ===")
        n = len(r.profile)
        L.append("        x[mm]      y[mm]     rayon[mm]")
        for i in range(0, n, max(1, n // 60)):
            qq = r.profile[i]
            L.append("  %10.4f %10.4f %10.4f" % (qq[0], qq[1], norm(qq)))
        self._lignes_constr = list(L)
        if not rendre:
            self.report.setPlainText("\n".join(L))
        col = "#7ed9a0" if r.valid else "#ffb050"
        self.lbl_sum.setText(
            "<b>Courses</b><br>%.2f / %.2f deg<br>"
            "<b>Couple coeur</b> %.4f N.mm<br>"
            "<b>Couple coeur</b> %.4f N.mm<br>"
            "<b>Plat</b> %.3f mm<br>"
            "<b>Corde pointe</b> %.3f mm<br>"
            "<b>Plage morte sommet</b><br>%.2f deg (plancher %.2f)<br>"
            "<span style='color:%s'>%s</span>"
            % (r.total_a, r.total_b, r.couple_coeur,
               p.couple_min * r.echelle, r.plat, r.corde,
               r.zone_morte, math.degrees(2.0 * math.atan(p.mu)), col,
               wrap("valide" if r.valid else " ; ".join(r.warnings))))

    def _face_radius(self, x):
        if not self.cb_face.isChecked():
            return 0.0
        # kappa est avant les grandeurs du marteau, qui ferment le vecteur
        nh = len(getattr(self.nsga, "marteau", ()) or ())
        k = x[-nh - 1] if nh else x[-1]
        return 0.0 if abs(k) < 1.0 / (20.0 * self.params.r_max) else 1.0 / k

    def _mu_contact(self, _R=0.0):
        """Frottement au contact : le glissement, le galet etant abandonne."""
        return None

    def _on_balayage(self):
        """Balaye pivot x pente et applique le meilleur point.

        Les deux variables sont celles qui portent : la position angulaire du
        pivot et la pente du pan. Le resultat est applique aux champs, donc
        visible tout de suite dans la simulation.
        """
        pr = self.prof
        if pr is None:
            QMessageBox.information(self, "Balayage",
                                    "Construire le coeur d'abord.")
            return
        crit = self.cb_crit.currentData()
        self.statusBar().showMessage("balayage en cours...")
        QApplication.processEvents()
        try:
            grille, best = balayage_marteau(self.params, pr, critere=crit)
        except Exception as exc:
            QMessageBox.critical(self, "Balayage", "Echec : %r" % exc)
            return
        if not best:
            self.statusBar().showMessage("balayage : aucun point exploitable")
            return
        depart = dict(pivot=(self.params.pivot_x, self.params.pivot_y),
                      angle=self.params.angle_plan)
        for nom, val in (("pivot_x", best["pivot"][0]),
                         ("pivot_y", best["pivot"][1]),
                         ("angle_plan", best["angle"])):
            setattr(self.params, nom, float(val))
            sb = self.spins.get(nom)
            if sb:
                sb.blockSignals(True)
                sb.setValue(float(val))
                sb.blockSignals(False)
        self._realigner_profil()
        self._resimuler_profil()
        self.report.appendPlainText(
            "\n=== Balayage pivot x pente (%d points) ===\n"
            "  critere                      %10s\n"
            "  depart  : pivot (%.2f, %.2f), pente %.2f deg\n"
            "  retenu  : pivot (%.2f, %.2f), pente %.2f deg\n"
            "  ondulation %.1f %%   effort %.4f N   plage morte %.2f deg\n"
            "  residu de l'ajustement du ressort %.1f %%\n"
            "  Le pivot n'a presque pas de prise radialement : sur +/- 4 mm,\n"
            "  l'ondulation varie de 1.6 point radialement et de 68.8 points\n"
            "  perpendiculairement a O->P. Seule sa position angulaire compte."
            % (len(grille), crit, depart["pivot"][0], depart["pivot"][1],
               depart["angle"], best["pivot"][0], best["pivot"][1],
               best["angle"], best["ondulation"], best["N_max"],
               best["plage_morte"], best["residu"]))
        self.statusBar().showMessage(
            "balayage : ondulation %.1f %%, effort %.4f N, plage %.2f deg"
            % (best["ondulation"], best["N_max"], best["plage_morte"]))

    def _on_nsga(self):
        if self.nsga is not None and self.nsga.isRunning():
            return
        sel = [i for i, cb in enumerate(self.cb_obj) if cb.isChecked()]
        if len(sel) < 2:
            QMessageBox.information(self, "Optimisation",
                                    "Cocher au moins deux objectifs.")
            return
        self.btn_nsga.setEnabled(False)
        self.btn_nsga.setText("Optimisation en cours...")
        self.nsga = NsgaThread(self.params.copy(),
                               self.sp_opt["cellules"].value(),
                               self.sp_opt["population"].value(),
                               self.sp_opt["generations"].value(),
                               self.sp_opt["scan"].value(), sel,
                               self.cb_face.isChecked(),
                               self.sp_opt["procs"].value(),
                               self.cb_galet.isChecked(),
                               self.cb_figer.isChecked(),
                               tuple(n for n, cb in self.cb_mart.items()
                                     if cb.isChecked()),
                               self.profil_importe)
        self.nsga.step.connect(
            lambda g, n, nf: self.statusBar().showMessage(
                "generation %d/%d - %d solutions faisables" % (g, n, nf)))
        self.nsga.done.connect(self._on_nsga_done)
        self.nsga.start()

    def _on_nsga_done(self, front, _graine):
        self.btn_nsga.setEnabled(True)
        self.btn_nsga.setText("Lancer NSGA-II")
        if not front:
            self.statusBar().showMessage(
                "optimisation impossible : "
                + getattr(self.nsga, "message", ""))
            return
        front = sorted(front, key=lambda q: (q["f4"][0], q["f4"][1]))
        self.front = front
        self.lst_front.blockSignals(True)
        self.lst_front.clear()
        for q in front:
            f = q["f4"]
            R = self._face_radius(q["x"])
            self.lst_front.addItem(
                "sous T %5.1f%%  plage %5.2f  ondul %5.0f%%  N %6.4f  %s%s"
                % (f[0], f[1], f[2], f[3],
                   "plan" if R == 0 else
                   ("galet %.1f" % R if self.cb_galet.isChecked()
                    else "R %+.1f" % R),
                   "" if q["viol"] <= 1e-9 else "  (!)"))
        self.lst_front.blockSignals(False)
        self.statusBar().showMessage(
            "front de %d solutions - cliquer pour en afficher une" % len(front))
        self._pareto_plot()
        self.lst_front.setCurrentRow(0)
        # apres la selection : afficher une solution reecrit le rapport
        self._bilan_front(front)

    def _bilan_front(self, front):
        """Compare le front a la graine, objectif par objectif.

        C'est la lecture qui manquait : savoir si l'optimisation a trouve
        mieux que la synthese dont elle part, et sur quoi. La graine etant
        dans la population, le front ne peut pas etre globalement plus
        mauvais ; ce tableau dit sur quels objectifs il gagne.
        """
        g = getattr(self.nsga, "graine_f4", None)
        L = ["", "=== Ce que l'optimisation a trouve ==="]
        if g is None:
            L.append("  graine non evaluee")
        else:
            L.append("  objectif                       graine    meilleur"
                     "   gain")
            for i, nom in enumerate(OBJECTIFS):
                best = min(q["f4"][i] for q in front)
                L.append("  %-28s %8.2f  %8.2f  %+7.2f"
                         % (nom, g[i], best, best - g[i]))
            L.append("")
            L.append("  Le couple vise est atteint par toutes les solutions :")
            L.append("  le ressort est recalcule pour chacune. On compare donc")
            L.append("  ce que ce resultat coute, pas le couple lui-meme.")
            if self.params.longueur_plan <= 1e-9:
                L.append("")
                L.append("  ! plan du marteau suppose infini : la plage morte")
                L.append("    reste alors au plancher d'arc-boutement, que le")
                L.append("    profil du coeur ne peut pas franchir. Renseigner")
                L.append("    une longueur de plan la fait tomber de 14.8 a")
                L.append("    1.5 deg - c'est le vrai levier.")
        self.report.appendPlainText("\n".join(L))

    def _pareto_plot(self):
        if not self.front:
            self.plot_pareto.set_points([])
            return
        sel = [i for i, cb in enumerate(self.cb_obj) if cb.isChecked()]
        ix, iy = (sel + [0, 1])[0], (sel + [0, 1])[1]
        self.plot_pareto.title = ("Front de Pareto : %s (x) contre %s (y)"
                                  % (OBJECTIFS[ix], OBJECTIFS[iy]))
        self.plot_pareto.xlabel = OBJECTIFS[ix]
        self.plot_pareto.set_curves([])
        self.plot_pareto.set_points(
            [(q["f4"][ix], q["f4"][iy]) for q in self.front], 0)

    def _on_pareto_click(self, i):
        if self.front and 0 <= i < len(self.front):
            self.lst_front.setCurrentRow(i)

    def _on_front_row(self, i):
        if not self.front or not (0 <= i < len(self.front)):
            return
        q = self.front[i]
        nq = ((3 if self.cb_face.isChecked() else 2)
              + len(getattr(self.nsga, "marteau", ()) or ()))
        # le gabarit du plat : sans lui, la solution affichee serait un
        # autre coeur que celui qui a ete evalue
        self._vue = "profil"
        # Une solution peut proposer un autre marteau : on reprend ses
        # valeurs dans les champs, sinon la simulation affichee ne serait pas
        # celle qui a ete evaluee.
        mart = getattr(self.nsga, "marteau", ()) or ()
        if mart:
            vals = dict(zip(mart, q["x"][len(q["x"]) - len(mart):]))
            t = vals.pop("pivot_t", None)
            if t is not None:
                # `pivot_t` est un deplacement compte depuis le pivot avec
                # lequel l'optimisation a ete lancee - **pas** depuis le pivot
                # courant. Le compter depuis le pivot courant le faisait
                # s'accumuler a chaque selection : de (16.00, 5.00) il passait
                # a (13.02, 10.98) en cinq clics, a distance a l'axe constante
                # donc en tournant autour du coeur, sans jamais revenir.
                base = getattr(self.nsga, "p", None) or self.params
                P0 = (base.pivot_x, base.pivot_y)
                d = norm(P0) or 1.0
                ut = (-P0[1] / d, P0[0] / d)
                vals.setdefault("pivot_x", P0[0] + t * ut[0])
                vals.setdefault("pivot_y", P0[1] + t * ut[1])
            for nom, val in vals.items():
                setattr(self.params, nom, float(val))
                sb = self.spins.get(nom)
                if sb:
                    sb.blockSignals(True)
                    sb.setValue(float(val))
                    sb.blockSignals(False)
            self._update_rules()
        # Le plat du coeur est porte a la normale de repos : si la solution
        # deplace la pente ou le pivot, gamma0 change et le plat doit suivre,
        # sinon le coeur ne part plus plat contre plat.
        gab = getattr(self.nsga, "gabarit", (0.0, 0.0))
        if mart and gab and gab[0] > 0.0:
            gab = (gab[0], self.params.gamma0)
        self.prof = Profil.from_vector(
            list(q["x"][:-nq]) + [q["x"][-nq], q["x"][-nq + 1]], *gab)
        R = self._face_radius(q["x"])
        fige = getattr(self.nsga, "ressort_fige", None)
        if fige:
            from .profil import Ressort
            self.ressort = Ressort()
            self.ressort.prearmage, self.ressort.raideur = fige
            self.ressort.eval = evaluer(
                self.params.copy(couple_min=fige[0], pente=fige[1]),
                self.prof, n=self.sp_opt["scan"].value(), rayon_face=R,
                mu_contact=self._mu_contact(R))
            self.ressort.course = self.ressort.eval.psi_max
            self.ressort.couple_maxi = fige[0] + fige[1] * self.ressort.course
            self.ressort.ok = True
        else:
            self.ressort = dimensionner_ressort(
                self.params, self.prof, R, self._mu_contact(R),
                n=self.sp_opt["scan"].value())
        self.pev = self.ressort.eval
        if self.ressort.ok:
            self._maj_ressort(self.ressort.prearmage, self.ressort.raideur,
                              calcule=True)
        self.plot_pareto.set_points(self.plot_pareto.points, i)
        self._show_profil()

    def _show_profil(self):
        if self.prof is None or self.pev is None:
            return
        self.lbl_graine.setText(wrap(
            "Geometrie simulee : contour importe du DXF, tel quel. "
            "L'optimisation partira de sa reconstruction."
            if getattr(self.prof, "M", 1) == 0 else
            ("Geometrie simulee : solution du front. Graine : profil "
             "importe." if self.profil_importe is not None else
             "Geometrie simulee : solution du front. Graine : synthese "
             "analytique.")))
        self.view.set_profil(self.params, self.prof, self.pev,
                             self.points_importes)
        self._profil_plots()
        self._profil_report()
        # les reperes portent la plage morte simulee : elle vient de changer
        self._update_rules()
        self._on_slider(self.slider.value())

    def _profil_plots(self):
        e = self.pev
        if e is None or not e.ok:
            for pl in (self.plot_torque, self.plot_eff, self.plot_geo):
                pl.set_curves([])
            return
        x = e.theta
        cible = self.params.couple_coeur
        self.plot_torque.set_curves([
            ("couple coeur [N.mm]", x, e.Th, ACC),
            ("+/- couple vise", [min(x), max(x)], [cible, cible], ACC3),
            ("", [min(x), max(x)], [-cible, -cible], ACC3)])
        self.plot_eff.set_curves([
            ("N contact [N]", x, e.N, ACC4),
            ("ouverture marteau [deg]", x,
             [math.degrees(v) for v in e.psi], ACC2)])
        beta = [1.0 * i for i in range(360)]
        self.plot_geo.set_curves([
            ("rayon de courbure [mm]", beta,
             [self.prof.courbure(math.radians(b)) for b in beta], ACC3),
            ("rayon du profil [mm]", x,
             [norm(self.prof.point(math.radians(b))) for b in x], ACC)])
        self.plot_torque.title = "Couple sur l'axe du coeur"
        self.plot_eff.title = "Efforts et course du marteau"
        self.plot_geo.title = "Courbure (contre la normale) et rayon"

    def _profil_report(self):
        p, pr, e, re_ = self.params, self.prof, self.pev, self.ressort
        L = ["=== Donnees ==="]
        L.append("  encombrement du coeur        %10.3f a %.3f mm"
                 % (p.r_min, p.r_max))
        L.append("  pivot du marteau             %10.3f , %.3f mm"
                 % (p.pivot_x, p.pivot_y))
        L.append("  pente du plan du marteau     %10.2f deg" % p.angle_plan)
        L.append("  longueur du plan             %10.4f mm%s"
                 % (p.longueur_plan,
                    "  (0 = suppose infini)" if p.longueur_plan <= 1e-9
                    else ""))
        if p.rayon_coin > 0:
            L.append("  rayon des coins du plan      %10.4f mm" % p.rayon_coin)
        L.append("    mesuree de O->P a la normale, au repos")
        L.append("    direction O->P               %10.2f deg / horizontale"
                 % math.degrees(p.alpha_p))
        L.append("    plan au repos                %10.2f deg / horizontale"
                 % p.orientation_plan())
        if e is not None and e.psi:
            i0 = min(range(len(e.phi)), key=lambda k: abs(e.phi[k]))
            i1 = max(range(len(e.phi)), key=lambda k: abs(e.phi[k]))
            L.append("    plan a pleine ouverture      %10.2f deg / horizontale"
                     % p.orientation_plan(e.phi[i1]))
            L.append("    soit une rotation de         %10.2f deg"
                     % abs(math.degrees(e.phi[i1] - e.phi[i0])))
        L.append("  couple de rappel vise        %10.4f N.mm" % p.couple_coeur)
        L.append("  frottement                   %10.4f" % p.mu)
        L.append("")
        L.append("=== Ressort de marteau (%s) ==="
                 % ("fige sur la version initiale"
                    if self.cb_figer.isChecked() else "sortie"))
        if re_ and re_.ok:
            L.append("  prearmage C0                 %10.4f N.mm"
                     % re_.prearmage)
            L.append("  raideur k                    %10.4f N.mm/rad"
                     % re_.raideur)
            L.append("    soit                       %10.4f N.mm/deg"
                     % (re_.raideur * math.pi / 180.0))
            L.append("  couple a pleine ouverture    %10.4f N.mm"
                     % re_.couple_maxi)
            L.append("  course du marteau            %10.3f deg"
                     % math.degrees(re_.course))
            if not self.cb_figer.isChecked():
                L.append("  residu de l'ajustement affine%10.1f %%"
                         % re_.ecart)
                L.append("    (ce qu'un ressort lineaire ne peut pas suivre)")
        L.append("")
        importe = getattr(pr, "M", 1) == 0
        L.append("=== Profil %s ==="
                 % ("importe (simule tel quel)" if importe else "obtenu"))
        if importe:
            L.append("  contour du DXF               %10d points"
                     % len(self.points_importes or []))
            L.append("  la courbure ci-dessous est indicative : sur un")
            L.append("  contour discretise elle vaut 0 sur les segments et")
            L.append("  diverge aux sommets, elle est donc lissee")
        else:
            L.append("  harmoniques                  %10d" % pr.M)
        L.append("  plat du marteau              %10s"
                 % ("plan" if abs(e.rayon_face) < 1e-9
                    else "R = %+.3f mm" % e.rayon_face))
        L.append("  rayon de courbure            %10.4f a %.4f mm"
                 % (e.rho_min, e.rho_max))
        L.append("  rayon de raccord impose      %10.4f mm" % p.rayon_raccord)
        L.append("  rayon du profil              %10.4f a %.4f mm"
                 % (e.r_min, e.r_max))
        L.append("  centre (cx, cy)              %10.4f , %.4f mm"
                 % (pr.cx, pr.cy))
        # Le plat du zero donne la position nulle : il faut pouvoir le lire.
        try:
            lp = longueur_plat(pr, p.gamma0, 0.002)
            L.append("  plat du zero (a 2 um pres)   %10.4f mm" % lp)
            vise = float(getattr(pr, "plat", 0.0) or 0.0)
            if vise > 1e-9:
                L.append("    vise par la synthese       %10.4f mm" % vise)
                if lp < 0.75 * vise:
                    L.append("    ! une serie a bande limitee ne porte pas un")
                    L.append("      plat exact : il faudrait un Dirac de")
                    L.append("      courbure. L'optimisation ne peut pas le")
                    L.append("      degrader davantage, mais pas le retablir.")
            elif lp < 1e-4:
                L.append("    ! aucun plat : pas de position nulle definie")
        except Exception:
            pass
        if self.points_importes and self.profil_importe is not None:
            ec = _ecart_contour(self.points_importes, self.profil_importe)
            L.append("  ecart dessin / reconstruction%10.4f mm (maxi)" % ec)
            if ec > 0.05 * p.r_max:
                L.append("  ! la reconstruction s'ecarte du dessin : le "
                         "contour est peut-etre non convexe, ou tronque")
        L.append("")
        L.append("=== Couple obtenu ===")
        from .profil import secteur_faible
        L.append("  part du tour sous la consigne%10.1f %%"
                 % secteur_faible(e, p.couple_coeur))
        L.append("  couple mini garanti          %10.5f N.mm (vise %.4f)"
                 % (e.couple_utile, p.couple_coeur))
        L.append("    sur %.0f %% du tour, hors les deux equilibres"
                 % (100.0 * e.part_utile))
        L.append("  couple moyen                 %10.5f N.mm" % e.couple_moyen)
        L.append("  ondulation                   %10.2f %%"
                 % (100.0 * e.ondulation))
        L.append("  effort de contact maxi       %10.4f N" % e.N_max)
        if p.longueur_plan > 1e-9:
            L.append("  contacts sur un coin du plan %10d / %d positions"
                     % (e.contacts_coin, len(e.theta)))
            if e.contacts_coin:
                L.append("    le plan est trop court pour rester tangent :")
                L.append("    c'est son coin qui touche le coeur")
        else:
            L.append("  longueur du plan non renseignee : il est suppose")
            L.append("    toujours assez long (tangence pure)")
        L += self._bloc_pointe(p, pr, e)
        muc = self._mu_contact(e.rayon_face)
        muc = p.mu if muc is None else muc
        L.append("  plage morte                  %10.2f deg" % e.plage_morte)
        L.append("    mesuree sur la marge |s| - mu.H interpolee,")
        L.append("    donc independante du pas de calcul")
        L.append("    pointe vive en appui plan : 2.atan(mu) = %.2f deg ;"
                 % math.degrees(2.0 * math.atan(muc)))
        L.append("    un coin arrondi passe sous cette valeur")
        L.append("")
        L.append("=== Coherence mecanique ===")
        L.append("  plages de signe du couple    %10d (attendu 2)"
                 % e.runs_couple)
        L.append("  changements de cote marteau  %10d (attendu 0)"
                 % e.cote_change)
        L.append("  sauts de contact             %10d (attendu 0)"
                 % e.saut_contact)
        L.append("  positions en interference    %10d (attendu 0)"
                 % e.interference)
        L.append("")
        L.append("=== Profil (repere du coeur) ===")
        pts = pr.contour(360)
        L.append("        x[mm]      y[mm]     rayon[mm]")
        for i in range(0, len(pts), 6):
            L.append("  %10.4f %10.4f %10.4f"
                     % (pts[i][0], pts[i][1], norm(pts[i])))
        self.report.setPlainText("\n".join(L))
        ok = (e.couple_utile >= 0.95 * p.couple_coeur
              and e.rho_min >= p.rayon_raccord - 1e-6 and e.runs_couple == 2
              and e.r_max <= p.r_max * 1.001 and not e.interference
              and not e.cote_change and not e.saut_contact)
        self.lbl_opt.setText(
            "<b>Ressort : C0</b> %.3f N.mm<br>"
            "<b>raideur k</b> %.3f N.mm/rad<br>"
            "<b>Couple coeur</b> %.4f N.mm<br>(vise %.3f)<br>"
            "<b>Plage morte</b> %.2f deg<br>"
            "<b>Ondulation</b> %.1f %%<br>"
            "<b>N contact</b> %.4f N<br>"
            "<span style='color:%s'>%s</span>"
            % ((re_.prearmage if re_ else 0.0), (re_.raideur if re_ else 0.0),
               e.couple_utile, p.couple_coeur, e.plage_morte,
               100.0 * e.ondulation, e.N_max,
               "#7ed9a0" if ok else "#ffb050",
               wrap("contraintes respectees" if ok
                    else "contrainte non respectee")))

    def _on_export_profil(self, pr=None, ev=None, defaut="coeur_optimise.dxf"):
        pr = self.prof if pr is None else pr
        ev = self.pev if ev is None else ev
        if pr is None or ev is None:
            QMessageBox.information(self, "Export", "Aucun profil a exporter.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Exporter en DXF",
                                              defaut, "DXF (*.dxf)")
        if not path:
            return
        path = avec_dxf(path)
        from .dxfio import Dessin, _dans_repere_coeur
        d = Dessin()
        d.poly(pr.contour(720), "COEUR")
        d.circle((0.0, 0.0), self.params.r_min, "CONSTRUCTION")
        d.circle((0.0, 0.0), self.params.r_max, "CONSTRUCTION")

        # le coeur est dans son repere : on y ramene aussi le marteau, sinon
        # le plan de contact traverse le coeur sur le plan
        def trace(x, calque):
            ps = pose_profil(self.params, pr, ev, x)
            if not ps:
                return
            th = ps["theta"]
            f = _dans_repere_coeur(ps["face"], th)
            piv = _dans_repere_coeur([self.params.pivot], th)[0]
            d.poly(f, calque, closed=False)
            d.line(piv, _dans_repere_coeur([ps["contact"]], th)[0], calque)
            d.circle(piv, max(0.2, 0.05 * self.params.r_max), calque)

        # Un seul marteau par defaut, a la position de repos : douze
        # marteaux autour du coeur ne servent qu'a verifier la cinematique.
        trace(0.0, "PLAN_CONTACT")
        if self.cb_cine.isChecked():
            for i in range(1, 12):
                trace(360.0 * i / 12, "MARTEAU_CINEMATIQUE")
        if self.cb_arcs.isChecked() and self.approx is not None:
            for A in self.approx.arcs:
                c, r, a0, a1 = A.pour_dxf()
                d.arc(c, r, a0, a1, "COEUR_ARCS")
        d.enregistrer(path)
        self.statusBar().showMessage("ecrit : " + path)

    def _profil_affiche(self):
        """Le profil dont on parle : importe ou optimise, sinon re-simule."""
        if self.prof is not None:
            return self.prof
        return getattr(self.view, "prof_fini", None)

    def _maj_arcs(self, *_a):
        """Recalcule la chaine d'arcs et rafraichit l'affichage."""
        pr = self._profil_affiche()
        self.approx = None
        if pr is not None and self.cb_arcs.isChecked():
            try:
                self.approx = arcs.approcher(pr, self.sp_arc.value())
            except Exception as exc:
                self.statusBar().showMessage("arcs : %r" % exc)
        self.view.set_arcs(self.approx if self.cb_arcs.isChecked() else None)
        if self.approx is not None:
            self.statusBar().showMessage(
                "%d arcs pour un ecart de %.4f mm (vise %.4f)"
                % (self.approx.n, self.approx.ecart, self.sp_arc.value()))
        if self.is_optim and self.prof is not None:
            self._refresh_report_profil_arcs()

    def _refresh_report_profil_arcs(self):
        """Ajoute le bilan des arcs au rapport, sans tout recalculer."""
        if self.approx is None:
            return
        a = self.approx
        R = a.rayons()
        self.report.appendPlainText(
            "\n=== Profil en arcs de cercle ===\n"
            "  arcs                         %10d\n"
            "  ecart maxi obtenu            %10.4f mm (vise %.4f)\n"
            "  residu de fermeture          %10.2e mm\n"
            "  rayons                       %10.3f a %.3f mm\n"
            "  dont concaves                %10d"
            % (a.n, a.ecart, self.sp_arc.value(), a.fermeture,
               min(R), max(R), sum(1 for x in a.arcs if x.concave)))

    def _on_sauver(self):
        path, _ = QFileDialog.getSaveFileName(self, "Enregistrer l'etude",
                                              "etude.json", "JSON (*.json)")
        if not path:
            return
        # Le sujet de l'etude est le contour importe ou optimise, quel que
        # soit l'onglet affiche. Prendre celui de l'onglet actif enregistrait
        # le coeur **synthetise** quand on revenait sur l'onglet d'une
        # synthese apres un import : l'etude relue ne contenait plus la
        # geometrie de l'utilisateur.
        pr = self.prof
        src = "importe ou optimise"
        if pr is None:
            pr = getattr(self.view, "prof_fini", None)
            src = "synthese, re-simulee avec le plan reel"
        reglages = {k: sb.value() for k, sb in self.sp_opt.items()}
        try:
            f = projet.sauver(path, self.params, pr, self.ressort, reglages,
                              commentaire="contour : " + src)
        except Exception as exc:
            QMessageBox.critical(self, "Enregistrer",
                                 "Ecriture impossible : %r" % exc)
            return
        self.statusBar().showMessage("etude enregistree : " + f)

    def _on_ouvrir(self):
        path, _ = QFileDialog.getOpenFileName(self, "Ouvrir une etude", "",
                                              "JSON (*.json);;Tous (*)")
        if not path:
            return
        try:
            vals, pts, reglages, log = projet.charger(path, self.params)
        except Exception as exc:
            QMessageBox.critical(self, "Ouvrir",
                                 "Lecture impossible : %r" % exc)
            return
        for n, sb in self.spins.items():
            if n in vals:
                sb.blockSignals(True)
                sb.setValue(float(getattr(self.params, n)))
                sb.blockSignals(False)
        for k, sb in self.sp_opt.items():
            if k in reglages:
                sb.setValue(int(reglages[k]))
        self._update_rules()
        if pts:
            # un contour enregistre se simule tel quel, comme un import DXF
            self._charger_profil(list(pts), log)
        else:
            self.cres = None
            self._on_chrono()
        self.report.appendPlainText("\n--- etude %s ---\n%s"
                                    % (os.path.basename(path), "\n".join(log)))
        self.statusBar().showMessage("etude ouverte : " + path)

    def _marquer(self, noms, couleur="#ffb050", motif=""):
        """Met en evidence des champs du panneau, et explique pourquoi."""
        for n, lbl in self.labels.items():
            txt = lbl.text()
            if n in noms:
                lbl.setStyleSheet("color:%s; font-weight:bold" % couleur)
                lbl.setToolTip((motif + "<br><br>" if motif else "")
                               + self.aides_champs.get(n, txt))
            else:
                lbl.setStyleSheet("")

    def _mesurer_depuis_dxf(self, pr, log):
        """Remplit les champs de construction a partir du contour lu.

        Repartir d'un dessin pour definir un nouveau coeur demande de
        retrouver, dans le contour, ce que la synthese prend en entree. Le
        frottement et le couple de rappel voulu n'y sont pas : l'un est une
        donnee physique, l'autre une specification. Ils sont mis en evidence
        dans le panneau pour que vous les renseigniez.
        """
        try:
            mes, manque = arcs.mesurer_coeur(pr, self.params.gamma0)
            mes["plat_souhaite"] = longueur_plat(pr, self.params.gamma0, 0.002)
        except Exception as exc:
            log.append("mesures sur le contour impossibles : %r" % exc)
            return
        for n, v in sorted(mes.items()):
            sb = self.spins.get(n)
            if sb is None:
                continue
            setattr(self.params, n, float(v))
            sb.blockSignals(True)
            sb.setValue(float(v))
            sb.blockSignals(False)
        log.append("mesure sur le contour : " + ", ".join(
            "%s = %.4f" % (n, v) for n, v in sorted(mes.items())))
        log.append("mesures approchees, a verifier : elles viennent de "
                   "l'approximation par arcs du dessin")
        self._marquer(manque, motif=(
            "<b>A renseigner</b> : aucun contour ne porte cette grandeur. "
            "Le frottement est une donnee physique, le couple de rappel une "
            "specification."))
        log.append("a renseigner a la main : " + ", ".join(manque))

    def _on_import(self):
        path, _ = QFileDialog.getOpenFileName(self, "Importer un DXF", "",
                                              "DXF (*.dxf);;Tous (*)")
        if not path:
            return
        try:
            ents = dxfio.read_entities(path)
            vals, log = dxfio.interpret(ents, self.params)
        except Exception as exc:
            QMessageBox.critical(self, "Import DXF", "Lecture impossible : %r" % exc)
            return
        pts = vals.pop("profil", None)
        for n, v in vals.items():
            setattr(self.params, n, v)
            sb = self.spins.get(n)
            if sb:
                sb.blockSignals(True)
                sb.setValue(v)
                sb.blockSignals(False)
        self._update_rules()
        if pts:
            self._charger_profil(pts, log)
            if self.prof is not None:
                self._mesurer_depuis_dxf(self.prof, log)
        elif self.is_chrono:
            self.cres = None
        elif not self.is_optim:
            self.recompute()
        self.report.appendPlainText("\n--- import %s ---\n%s"
                                    % (os.path.basename(path), "\n".join(log)))
        self.statusBar().showMessage("DXF importe : " + ", ".join(sorted(vals)))

    def _resimuler_profil(self):
        """Rejoue la simulation du profil affiche, sans toucher a sa forme.

        La longueur du plat, le rayon des coins, le decalage, le frottement et
        le couple vise ne changent pas la geometrie du coeur : seule la
        simulation est a refaire. Le profil, lui, reste celui qui est affiche.
        """
        if self.prof is None:
            return
        # Le plat du coeur est porte a la normale de repos gamma0. Changer la
        # pente ou le pivot la deplace : le profil doit suivre, sinon le coeur
        # ne part plus plat contre plat et la simulation ne decrit plus rien.
        # (`_realigner_profil` ne s'applique qu'a un contour importe ; un
        # profil issu de la synthese porte son plat dans son gabarit.)
        pr = self.prof
        if (hasattr(pr, "a0") and getattr(pr, "plat", 0.0) > 0.0
                and abs(pr.beta0 - self.params.gamma0) > 1e-9):
            self.prof = pr = Profil(pr.a0, pr.an, pr.bn, pr.cx, pr.cy,
                                    pr.plat, self.params.gamma0)
            self.prof_dessine = pr
        self.ressort = dimensionner_ressort(
            self.params, pr, 0.0, None,
            n=self.sp_opt["scan"].value())
        self.pev = self.ressort.eval
        self._show_profil()
        self._maj_arcs()
        self.statusBar().showMessage(
            "simulation remise a jour : plage morte %.2f deg, effort %.4f N"
            % (self.pev.plage_morte, self.pev.N_max))

    def _realigner_profil(self):
        """Repose le plat du coeur sur le plan apres un changement de gamma0.

        Le recalage fait a l'import vaut pour l'angle de plan et le pivot de
        ce moment-la. Les changer deplace la normale du plan du marteau : sans
        recalage, le marteau recommence sa course deja ouvert et tout le
        parcours se decale.
        """
        pts = self.points_importes
        if not pts:
            return
        npts, ang = dxfio.aligner_contour(pts, self.params.pivot,
                                          self.params.angle_plan)
        if not ang:
            return
        self._charger_profil(npts, [], basculer=False)
        self.statusBar().showMessage(
            "coeur retourne de %.2f deg : plat sur plan au repos" % ang)

    def _charger_profil(self, pts, log, basculer=True):
        """Reprend le profil lu dans le DXF : l'origine en est l'axe.

        On construit deux objets :
          - le contour **tel quel**, qui est celui que l'on simule ;
          - sa reconstruction a bande limitee, qui ne sert que de graine a
            l'optimisation (celle-ci a besoin d'un vecteur de coefficients).
        """
        try:
            self.prof_dessine = ProfilImporte(pts)
            self.profil_importe = depuis_points(
                pts, self.sp_opt["cellules"].value(),
                self.params.rayon_raccord, self.params.r_max)
        except Exception as exc:
            log.append("profil non exploitable : %r" % exc)
            self.profil_importe = self.prof_dessine = None
            return
        self.points_importes = pts
        self.prof = self.prof_dessine          # on simule le dessin
        self.front = None
        self.lst_front.clear()
        self.ressort = dimensionner_ressort(
            self.params, self.prof, 0.0, None,
            n=self.sp_opt["scan"].value())
        self.pev = self.ressort.eval
        self._vue = "profil"
        self._show_profil()
        if basculer:
            log.append("le contour importe est simule tel quel ; sa "
                       "reconstruction sert de graine a l'optimisation")

    def _on_export(self):
        # Le bouton est contextuel : sur l'onglet d'optimisation c'est le
        # profil affiche (importe ou optimise) qu'il faut ecrire, pas le
        # coeur de la synthese analytique.
        if self.is_optim:
            return self._on_export_profil()
        # Sur les onglets analytiques, si un plan fini est renseigne, c'est la
        # re-simulation qui decrit la piece : exporter la pose de la synthese
        # ecrirait un marteau decolle du coeur.
        v = self.view
        if getattr(v, "prof_fini", None) is not None and v.ev_fini is not None:
            return self._on_export_profil(
                v.prof_fini, v.ev_fini,
                "coeur_chrono.dxf" if self.is_chrono else "coeur.dxf")
        if self.is_chrono:
            if self.cres is None or not self.cres.ok:
                QMessageBox.information(self, "Export",
                                        "Aucun profil a exporter.")
                return
            path, _ = QFileDialog.getSaveFileName(
                self, "Exporter en DXF", "coeur_chrono.dxf", "DXF (*.dxf)")
            if not path:
                return
            path = avec_dxf(path)
            dxfio.write_dxf_chrono(
                path, self.params, self.cres, 0.0,
                n_cine=12 if self.cb_cine.isChecked() else 0)
            self.statusBar().showMessage("ecrit : " + path)
            return
        if self.res is None or not self.res.ok:
            QMessageBox.information(self, "Export", "Aucun profil a exporter.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Exporter en DXF", "coeur.dxf",
                                              "DXF (*.dxf)")
        if not path:
            return
        path = avec_dxf(path)
        dxfio.write_dxf(path, self.params, self.res, self.view.index)
        self.statusBar().showMessage("ecrit : " + path)

    # -- calcul -------------------------------------------------------------
    def recompute(self):
        """La came simple n'existe plus : on construit le coeur de chrono."""
        self._on_chrono()

    def _refresh_plots(self):
        r = self.res
        if r is None or not r.ok:
            for p in (self.plot_torque, self.plot_eff, self.plot_geo):
                p.set_curves([])
            return
        x = [math.degrees(t) for t in r.theta]
        self.plot_torque.set_curves([
            ("|C| retour", x, [abs(v) for v in r.Th_retour], ACC),
            ("|C| armement", x, [abs(v) for v in r.Th_armement], ACC2),
            ("|C| sans frott.", x, [abs(v) for v in r.Th_ideal], ACC3),
        ])
        self.plot_eff.set_curves([
            ("Cm marteau [N.mm]", x, r.Cm, ACC2),
            ("N contact [N]", x, r.N, ACC4),
            ("ouverture [deg]", x, [math.degrees(v) for v in r.psi], ACC),
        ])
        self.plot_geo.set_curves([
            ("rayon [mm]", x, [norm(q) for q in r.profile], ACC),
            ("H appui [mm]", x, r.H, ACC3),
            ("courbure H+H'' [mm]", x, r.curv, ACC2),
        ])
        # Si un plan fini est renseigne, la courbe de couple de la synthese ne
        # decrit plus la piece : on superpose celle de la re-simulation, pour
        # que le graphique ne contredise pas le trace.
        e = getattr(self.view, "ev_fini", None)
        if e is not None and e.ok and e.theta:
            self.plot_torque.set_curves([
                ("|C| retour (plan infini)", x,
                 [abs(v) for v in r.Th_retour], ACC3),
                ("|C| plan reel", list(e.theta),
                 [abs(v) for v in e.Th], ACC),
            ])
            self.plot_eff.set_curves([
                ("Cm marteau [N.mm]", x, r.Cm, ACC2),
                ("N contact, plan reel [N]", list(e.theta), e.N, ACC4),
                ("ouverture, plan reel [deg]", list(e.theta),
                 [math.degrees(v) for v in e.phi], ACC),
            ])

    def _bloc_simulation(self):
        """Simulation du coeur construit, sans passer par l'optimisation.

        Le profil discretise est deduit de la synthese des que le coeur est
        construit : le couple reellement delivre, la plage morte, les contacts
        sur coin et le ressort sont donc disponibles tout de suite.
        """
        pr, e = self.prof, self.pev
        if pr is None or e is None or not e.ok:
            return []
        p = self.params
        L = ["", "=== Simulation du coeur construit ==="]
        L.append("  harmoniques du profil        %10d" % getattr(pr, "M", 0))
        L.append("  plat du zero (a 2 um pres)   %10.4f mm"
                 % longueur_plat(pr, p.gamma0, 0.002))
        L.append("  couple moyen sur le coeur    %10.5f N.mm (vise %.4f)"
                 % (e.couple_moyen, p.couple_coeur))
        L.append("  part du tour sous la consigne%10.1f %%"
                 % secteur_faible(e, p.couple_coeur))
        L.append("  ondulation                   %10.2f %%"
                 % (100.0 * e.ondulation))
        L.append("  effort de contact maxi       %10.4f N" % e.N_max)
        L.append("  plage morte                  %10.2f deg" % e.plage_morte)
        if p.longueur_plan > 1e-9:
            L.append("  contacts sur un coin du pan  %10d / %d"
                     % (e.contacts_coin, len(e.theta)))
            L.append("  hors arrondi dans +/-10 deg  %10.1f %% de la fenetre"
                     % fenetre_pointe(e, 10.0))
        else:
            L.append("  pan suppose infini : renseigner sa longueur pour")
            L.append("    voir ce que le marteau reel delivre")
        if self.ressort is not None and self.ressort.ok:
            L.append("  residu de l'ajustement affine%10.1f %%"
                     % self.ressort.ecart)
        return L

    def _bloc_encombrement(self, p, r):
        """Encombrement reellement utilise, et monotonie du rayon.

        Le diametre maxi est une **limite**, pas un but : il est mesure sur
        tout le profil, pas seulement a la pointe. Un coeur non monotone peut
        donc tenir dans moins de place que son rayon de pointe ne le
        suggererait, et c'est cette mesure qui le dit.
        """
        if not r.profile:
            return []
        ray = [math.hypot(q[0], q[1]) for q in r.profile]
        re_ = max(ray)
        i_max = ray.index(re_)
        L = ["", "=== Encombrement ==="]
        L.append("  rayon maxi atteint           %10.4f mm sur %.4f admis"
                 % (re_, p.r_max))
        L.append("  rayon mini (plat)            %10.4f mm" % min(ray))
        # monotonie du rayon le long de chaque flanc
        souple = []
        for nom, res in (("A", r.A), ("B", r.B)):
            if not (res and res.ok):
                continue
            rr = [math.hypot(*_X(res.H[i], res.s[i], res.beta[i]))
                  for i in range(len(res.beta))]
            baisses = sum(1 for i in range(len(rr) - 1) if rr[i + 1] < rr[i])
            souple.append((nom, baisses, len(rr)))
        if souple:
            mono = all(b == 0 for _n, b, _t in souple)
            L.append("  rayon monotone croissant     %10s"
                     % ("oui" if mono else "non"))
            if not mono:
                for nom, b, t in souple:
                    if b:
                        L.append("    flanc %s : %d descentes sur %d points"
                                 % (nom, b, t))
                L.append("    un coeur non monotone peut tenir dans moins de")
                L.append("    place : le rayon maxi ci-dessus est la mesure")
                L.append("    qui compte, pas celui de la pointe.")
        return L

    def _bloc_plan_fini(self):
        """Resultats de la re-simulation avec le plan reel du marteau.

        Les deux syntheses analytiques raisonnent sur un plan infini. Des que
        l'utilisateur renseigne une longueur de plan, un rayon de coin ou un
        decalage, c'est cette re-simulation qui dit ce que la piece fait
        vraiment : c'est elle qui est tracee, et ce sont ses chiffres qui
        comptent.
        """
        v = self.view
        pr, e = getattr(v, "prof_fini", None), getattr(v, "ev_fini", None)
        if pr is None or e is None or not e.ok:
            return []
        p = self.params
        L = ["", "=== Avec le plan reel du marteau (re-simulation) ==="]
        L.append("  les chiffres ci-dessus supposent le plan infini ;")
        L.append("  ceux-ci tiennent compte de sa longueur, de ses")
        L.append("  arrondis et de son decalage")
        L.append("  couple moyen sur le coeur     %9.5f N.mm" % e.couple_moyen)
        L.append("  effort de contact maxi        %9.4f N" % e.N_max)
        L.append("  plage morte                   %9.2f deg" % e.plage_morte)
        L.append("  contacts sur un coin du plan  %9d / %d"
                 % (e.contacts_coin, len(e.theta)))
        hors = fenetre_pointe(e, 10.0)
        L.append("  hors arrondi dans +/-10 deg   %9.1f %% de la fenetre"
                 % hors)
        if hors > 1e-9:
            L.append("  ! l'arete du plan rencontre la pointe :")
            L.append("    la plage morte y remonte au plancher")
        marche = getattr(e, "marche_deg", 0.0)
        if marche > 1e-9:
            L.append("  sans contact (marche pointe -> plat) %4.1f deg"
                     % marche)
            L.append("    le coeur y est concave : ni le plan ni son arrondi")
            L.append("    n'y touchent, le marteau y tombe. Il est trace en")
            L.append("    orange sur ce secteur.")
        L += self._bloc_compromis(p, pr)
        return L

    def _bloc_pointe(self, p, pr, e):
        """Controle du contact sur l'arrondi au passage de la pointe.

        A la pointe, un contact sur l'arete vive du plan annule le bras de
        levier : la plage morte vaut alors 2.atan(mu), pres de 14 deg. Quand
        c'est l'arrondi du marteau qui porte, elle tombe a quelques degres. On
        exige donc l'arrondi sur toute la fenetre de +/- 10 deg autour de la
        pointe, et on annonce la longueur de plan maxi qui le garantit.
        """
        hors = fenetre_pointe(e, 10.0)
        L = ["  hors arrondi dans +/-10 deg   %9.1f %% de la fenetre" % hors]
        if hors <= 1e-9:
            L.append("    l'arrondi porte sur tout le passage de la pointe")
        else:
            L.append("  ! l'arete du plan rencontre la pointe : le bras de")
            L.append("    levier s'y annule et la plage morte remonte")
        L += self._bloc_compromis(p, pr, e.rayon_face)
        return L

    def _bloc_compromis(self, p, pr, rayon_face=0.0):
        """Compromis sur la longueur du plan, balaye et non tranche.

        Deux exigences s'opposent : un plan court met son arrondi a la pointe
        et effondre la plage morte, mais trop court il ne peut plus enjamber
        la marche qui ramene la pointe au plat, et le marteau tombe sans rien
        toucher. La fenetre de pointe n'etant pas monotone en longueur de
        plan, on balaye au lieu de chercher une limite par dichotomie.
        """
        try:
            best, tab, front = plan_recommande(
                p, pr, n_pts=11, n=180, rayon_face=rayon_face,
                mu_contact=self._mu_contact(rayon_face))
        except Exception:
            return []
        if not tab:
            return []
        ids = {id(d) for d in front}
        L = ["", "  compromis sur la longueur du plan :",
             "    L [mm]  hors pointe  plage morte  sans contact   N [N]"]
        for d in tab:
            L.append("  %s %6.2f    %6.1f %%     %6.2f deg    %5.1f deg  %6.4f"
                     % ("->" if id(d) in ids else "  ",
                        d["L"], d["hors_pointe"], d["plage_morte"],
                        d["marche"], d["N_max"]))
        if best:
            L.append("  -> conseille : %.3f mm (plage morte %.2f deg, "
                     "sans contact %.1f deg)"
                     % (best["L"], best["plage_morte"], best["marche"]))
            L.append("     plus petite plage morte sans lacher le coeur")
        else:
            L.append("  aucune longueur ne garde le contact : la marche de ce")
            L.append("  coeur est trop profonde pour qu'un plan l'enjambe.")
            L.append("  Les lignes marquees -> sont les compromis non domines")
            L.append("  entre plage morte et secteur sans contact ; l'arbitrage")
            L.append("  vous revient. Attention : la plage morte perd son sens")
            L.append("  quand une grande part du tour est sans contact.")
        return L

    def _bloc_plan(self, p, demi_requis):
        """Longueur de plan exigee par la synthese, et avertissement.

        La synthese analytique cherche la tangence **sans borne d'abscisse** :
        elle suppose donc le plan assez long. Si le plan reel est plus court,
        le contact passe sur son coin et le couple n'est plus celui calcule
        ici. Seul l'onglet d'optimisation modelise le plan fini.
        """
        if demi_requis <= 0.0:
            return []
        L = ["  longueur de plan exigee      %10.4f mm" % (2.0 * demi_requis)]
        if p.longueur_plan <= 1e-9:
            L.append("    plan suppose assez long (longueur_plan = 0)")
        elif p.longueur_plan < 2.0 * demi_requis - 1e-9:
            L.append("  ! votre plan ne fait que %.4f mm : le contact sort"
                     % p.longueur_plan)
            L.append("    du plan et passe sur son coin. Le couple ci-dessus")
            L.append("    ne vaut plus ; utilisez l'onglet d'optimisation,")
            L.append("    seul a modeliser le plan fini et ses arrondis.")
        else:
            L.append("    votre plan de %.4f mm suffit" % p.longueur_plan)
        return L

    def _refresh_report(self):
        r = self.res
        p = self.params
        L = ["=== Parametres ==="]
        for name, label, unit, _a, _b, _s in Params.SPECS:
            L.append("  %-26s %10.4f %s" % (label, getattr(p, name), unit))
        L.append("")
        if r is None or not r.ok:
            L.append("=== Echec ===")
            L.append("  " + (r.message if r else "pas de calcul"))
            self.report.setPlainText("\n".join(L))
            self.lbl_sum.setText("<b style='color:#f08'>%s</b>"
                                 % (r.message if r else ""))
            return
        L.append("=== Resultats ===")
        L.append("  couple coeur constant        %10.5f N.mm" % abs(r.T0))
        L.append("  rendement (avec/sans frott.) %10.3f" % r.rendement)
        L.append("  couple marteau C0 -> C max   %10.4f -> %.4f N.mm"
                 % (min(r.Cm), max(r.Cm)))
        L.append("  augmentation du couple       %10.4f N.mm (%.2f %%)"
                 % (r.delta_Cm, 100.0 * r.delta_Cm / max(min(r.Cm), 1e-9)))
        L.append("  course du marteau            %10.3f deg" % math.degrees(r.psi_max))
        L.append("  effort normal maxi           %10.4f N" % r.N_max)
        L.append("  bras de levier mini (/pivot) %10.4f mm" % r.jam_min)
        L.append("  rayon de courbure mini       %10.4f mm" % r.curv_min)
        L.append("  rayon du profil              %10.4f -> %.4f mm"
                 % (min(norm(q) for q in r.profile), r.r_eff))
        L.append("  ecart couple rentrant/sortant%10.2f %%" % (100.0 * r.ripple))
        L += self._bloc_plan(p, getattr(r, "demi_requis", 0.0))
        L += self._bloc_plan_fini()
        L.append("")
        L.append("=== Diagnostic ===")
        if r.valid:
            L.append("  profil exploitable (convexe, monotone, sans arc-boutement)")
        for wmsg in r.warnings:
            L.append("  ! " + wmsg)
        if r.message:
            L.append("  > " + r.message)
        L.append("")
        L.append("=== Profil (repere du coeur) ===")
        L.append("   theta[deg]   phi[deg]    H[mm]     s[mm]      x[mm]     y[mm]")
        n = len(r.theta)
        for i in range(0, n, max(1, n // 40)):
            L.append("  %10.2f %10.3f %9.4f %9.4f %10.4f %10.4f"
                     % (math.degrees(r.theta[i]), math.degrees(r.phi[i]),
                        r.H[i], r.s[i], r.profile[i][0], r.profile[i][1]))
        self.report.setPlainText("\n".join(L))

        col = "#7ed9a0" if r.valid else "#ffb050"
        self.lbl_sum.setText(
            "<b>Couple coeur</b> %.5f N.mm<br>"
            "<b>Augmentation Cm</b> %.4f N.mm (%.1f %%)<br>"
            "<b>Course marteau</b> %.2f deg<br>"
            "<b>N maxi</b> %.4f N<br>"
            "<b>R profil</b> %.3f mm<br>"
            "<span style='color:%s'>%s</span>"
            % (abs(r.T0), r.delta_Cm, 100.0 * r.delta_Cm / max(min(r.Cm), 1e-9),
               math.degrees(r.psi_max), r.N_max, r.r_eff, col,
               "valide" if r.valid else " ; ".join(r.warnings) or r.message))


def run():
    app = QApplication([])
    theme.appliquer(app)
    w = MainWindow()
    w.show()
    return app.exec()
