"""Navigation a la souris dans un graphique matplotlib, comme dans une CAO.

Meme logique que la vue du programme de came (`came.ui.CamView`) :

    molette          zoom centre sur le curseur
    glisser          deplace la vue
    double-clic      revient a la vue d'ensemble

La barre d'outils de matplotlib reste utilisable : quand un de ses outils
(loupe, main) est enclenche, c'est lui qui recoit la souris.
"""

#: meme progressivite que la vue de la came : 1.0015 par unite de molette,
#: soit 120 unites par cran
ZOOM_PAR_CRAN = 1.0015 ** 120


class Navigation:

    def __init__(self, canvas):
        self.canvas = canvas
        self._axe = None            # axes en cours de deplacement
        self._vues = {}             # axes -> (vue d'ensemble, derniere vue posee ici)
        self._a_relire = set()      # axes dont la vue posee attend son trace
        for evenement, f in (("scroll_event", self._molette),
                             ("button_press_event", self._appui),
                             ("motion_notify_event", self._glisse),
                             ("button_release_event", self._relache),
                             ("draw_event", self._trace)):
            canvas.mpl_connect(evenement, f)

    # -- vue d'ensemble ------------------------------------------------------
    @staticmethod
    def _limites(ax):
        return ax.get_xlim(), ax.get_ylim()

    def _avant(self, ax):
        """A appeler avant de toucher a la vue : retient la vue d'ensemble.

        Si les limites ne sont plus celles posees ici, c'est que le programme
        a retrace le graphique : sa vue devient la nouvelle vue d'ensemble.
        """
        v = self._vues.get(ax)
        if v is None or v[1] != self._limites(ax):
            self._vues[ax] = (self._limites(ax), None)
            barre = self.canvas.toolbar
            if barre is not None:   # le bouton « origine » de la barre y revient aussi
                barre.update()
                barre.push_current()

    def _apres(self, ax):
        self._vues[ax] = (self._vues[ax][0], self._limites(ax))
        self._a_relire.add(ax)
        self.canvas.draw_idle()

    def _trace(self, _e):
        # Des axes a echelle egale ajustee sur les donnees (`datalim`) voient
        # leurs limites retouchees par matplotlib au trace : on relit ce qu'il
        # en a fait, sinon cette retouche passerait pour un nouveau trace du
        # programme et la vue d'ensemble serait perdue.
        for ax in self._a_relire:
            if ax in self._vues:
                self._vues[ax] = (self._vues[ax][0], self._limites(ax))
        self._a_relire.clear()

    def vue_utilisateur(self, ax):
        """(xlim, ylim) si l'utilisateur a zoome ou deplace ces axes, sinon None.

        Pour un trace refait a chaque image (animation) : il reprend cette vue
        au lieu de ses limites par defaut, et le zoom tient pendant la lecture.
        """
        v = self._vues.get(ax)
        if v is not None and v[1] is not None and v[1] == self._limites(ax):
            return v[1]
        return None

    def ajuster(self, ax):
        """Revient a la vue d'ensemble."""
        v = self._vues.pop(ax, None)
        if v is not None and v[1] == self._limites(ax):
            ax.set_xlim(v[0][0])
            ax.set_ylim(v[0][1])
            self.canvas.draw_idle()

    # -- souris --------------------------------------------------------------
    def _libre(self, e):
        """Axes sous le curseur, si aucun outil de la barre n'est enclenche."""
        barre = self.canvas.toolbar
        if barre is not None and str(barre.mode):
            return None
        ax = e.inaxes
        if ax is None or not ax.get_navigate() or hasattr(ax, "_colorbar"):
            return None
        return ax

    def _molette(self, e):
        ax = self._libre(e)
        if ax is None or not e.step:
            return
        self._avant(ax)
        k = ZOOM_PAR_CRAN ** e.step
        # Le rectangle des axes est retreci autour du curseur, en pixels, puis
        # relu en coordonnees du trace : le point sous le curseur ne bouge pas,
        # quels que soient l'echelle (lineaire ou non) et le sens des axes.
        b = ax.bbox
        coins = [(e.x + (x - e.x) / k, e.y + (y - e.y) / k)
                 for x, y in ((b.x0, b.y0), (b.x1, b.y1))]
        (x0, y0), (x1, y1) = ax.transData.inverted().transform(coins)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        self._apres(ax)

    def _appui(self, e):
        ax = self._libre(e)
        if ax is None:
            return
        if e.dblclick:
            self._axe = None
            self.ajuster(ax)
            return
        self._avant(ax)
        ax.start_pan(e.x, e.y, 1)
        self._axe = ax

    def _glisse(self, e):
        if self._axe is not None:
            self._axe.drag_pan(1, None, e.x, e.y)
            self._apres(self._axe)

    def _relache(self, _e):
        if self._axe is not None:
            self._axe.end_pan()
            self._axe = None


def activer(canvas):
    """Donne a `canvas` la navigation a la souris et la renvoie.

    Elle est accrochee au canvas : matplotlib ne garde qu'une reference faible
    sur ses rappels, sans cela elle disparaitrait aussitot.
    """
    canvas.navigation = Navigation(canvas)
    return canvas.navigation
