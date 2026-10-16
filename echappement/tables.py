"""Tables de contact d'un tracé, en fonction de l'angle d'ancre.

Pour chaque angle d'ancre alpha d'une grille fine couvrant toute la course :

  - pour chaque palette, la position de la roue (poussée dans son sens) où
    une dent touche cette palette : Phi(alpha), avec le point et la normale
    de contact et les éléments touchés. La roue ayant des dents identiques,
    Phi n'est connu qu'à un pas près : c'est la dent qui précède la palette.
    Phi saute d'un coup quand la dent s'échappe (fin d'impulsion) : la
    suivante prend le relais. Hors de la roue, la palette ne touche rien
    (nan) ;
  - pour l'ellipse, les deux positions du balancier où elle touche un flanc
    de la fourchette (bornes basse et haute), nan quand l'ellipse peut
    sortir librement de ce côté (aux butées).

Toute la mécanique (relevé, statique, dynamique) lit ces tables : le calcul
géométrique n'est fait qu'une fois par tracé. Les contacts sont cherchés sur
la géométrie exacte du DXF (segments et vrais arcs, voir
geometrie.premier_contact_exact) : les rapports de transmission qu'on en
dérive n'ont pas les crans d'un contour en facettes.

Angles en radians, comptés depuis la pose dessinée (sens trigonométrique de
la vue standard).
"""

from dataclasses import dataclass

import numpy as np

from .geometrie import distance, premier_contact_exact, primitives, rotation

DEG = np.pi / 180.0
#: recul de la roue avant la recherche de contact (le recul au dégagement
#: reste bien en deçà)
RECUL = 1.0 * DEG
#: course de roue explorée d'abord en suivant une table (rad)
COURSE_COURTE = 4 * RECUL
#: débord de la grille au-delà des butées
DEBORD = 0.5 * DEG


def _sous_ensemble(corps, predicat):
    """Éléments de `corps` dont l'étiquette vérifie `predicat` : leurs
    segments (pose dessinée, pour le test de liberté), leurs primitives
    exactes (pour les contacts) et les éléments eux-mêmes."""
    el = [e for e in corps.elements if predicat(e.etiquette)]
    S = np.concatenate([np.stack([e.points[:-1], e.points[1:]], axis=1) for e in el])
    return S, primitives(el), el


def course_butees(trace):
    """Angles d'ancre (depuis la pose dessinée) des deux butées.

    La recherche part de l'ancre placée au milieu des deux butées (baguette
    dirigée vers le milieu de leurs centres) : la pose dessinée peut, quand
    on resserre les butées, déjà les chevaucher."""
    t = trace
    bag = t.ancre.element("baguette").points
    A0 = t.ancre.axe
    bout = max((bag[0], bag[-1]), key=lambda p: np.hypot(*(p - A0)))
    milieu = np.mean([b.axe for b in t.butees], axis=0)
    u, v = bout - A0, milieu - A0
    a_c = float(np.arctan2(u[0] * v[1] - u[1] * v[0], u @ v))
    P_anc = primitives(t.ancre.elements).tournees(A0, a_c)
    P_but = primitives([e for b in t.butees for e in b.elements])
    a = [a_c + sens * premier_contact_exact(P_anc, A0, sens, P_but, course=60 * DEG)[0]
         for sens in (-1, 1)]
    return min(a), max(a)


@dataclass
class TableRoue:
    """Contact roue-palette pour une palette."""
    palette: str
    phi: np.ndarray          # angle de roue au contact (nan : palette hors de la roue)
    P: np.ndarray            # point de contact (n, 2)
    N: np.ndarray            # normale palette -> roue (n, 2)
    sur_ancre: list          # étiquette de l'élément de palette touché
    sur_roue: list           # étiquette de l'élément de dent touché


@dataclass
class TableFourche:
    """Contact ellipse-fourchette, borne basse ou haute du balancier."""
    beta: np.ndarray         # angle du balancier au contact (nan : sortie libre)
    P: np.ndarray
    N: np.ndarray            # normale fourchette -> ellipse
    sur_ancre: list


class Tables:
    """Tables de contact d'un Trace."""

    def __init__(self, trace, pas=0.02 * DEG):
        self.trace = trace
        t = trace
        self.s = t.sens_roue
        self.pas_roue = t.pas
        # éléments utiles de chaque pièce
        self.seg_palettes = {}
        for cote in ("entrée", "sortie"):
            self.seg_palettes[cote] = _sous_ensemble(t.ancre, lambda e, c=cote: e.endswith(c))
        fourche = ("baguette", "flanc", "corne", "fond")
        self.seg_fourche = _sous_ensemble(t.ancre, lambda e: e.startswith(fourche))
        r_min = min(np.hypot(*(S.reshape(-1, 2) - t.roue.axe).T).min()
                    for S, _, _ in self.seg_palettes.values())
        # de la roue, seuls les éléments qui montent assez haut peuvent
        # atteindre une palette
        garde = {id(e) for e in t.roue.elements
                 if np.hypot(*(e.points - t.roue.axe).T).max() > r_min - 0.3}
        self.seg_roue = _sous_ensemble(_Filtre(t.roue, garde), lambda e: True)
        self.seg_ellipse = _sous_ensemble(t.balancier, lambda e: True)

        self.a_min, self.a_max = self._butees()
        # La grille déborde des butées : le lissage des tables est mal tenu
        # à ses extrémités, et c'est justement aux butées que l'ancre repose.
        n = max(3, int(np.ceil((self.a_max - self.a_min) / pas)) + 1)
        h = (self.a_max - self.a_min) / (n - 1)
        k = int(np.ceil(DEBORD / h))
        self.alpha = self.a_min + h * np.arange(-k, n + k)
        self.i_min, self.i_max = k, k + n - 1
        self.beta_ldc = self._ligne_des_centres()
        self.roue = {c: self._balayer_palette(c) for c in ("entrée", "sortie")}
        self.fourche_bas, self.fourche_haut = self._bornes_fourche()

    def _butees(self):
        return course_butees(self.trace)

    def _ligne_des_centres(self):
        t = self.trace
        c = next(e.centre for e in t.balancier.elements if e.est_arc and e.rayon > 0.1)
        v, w = c - t.balancier.axe, t.ancre.axe - t.balancier.axe
        return float(np.arctan2(v[0] * w[1] - v[1] * w[0], v @ w))

    # -- roue ---------------------------------------------------------------
    def contact_roue(self, cote, alpha, phi_ref, course=None):
        """Premier contact de la roue, partant de phi_ref (position libre),
        avec la palette `cote` posée à alpha. Rend (phi, P, n, élément de
        palette, élément de dent) ou None."""
        t = self.trace
        _, P_pal, el_pal = self.seg_palettes[cote]
        _, P_r, el_r = self.seg_roue
        tt, P, n, jf, jm = premier_contact_exact(P_r.tournees(t.roue.axe, phi_ref), t.roue.axe,
                                                 self.s, P_pal.tournees(t.ancre.axe, alpha),
                                                 course=course or self.pas_roue)
        if not np.isfinite(tt):
            return None
        return phi_ref + self.s * tt, P, n, el_pal[jf].etiquette, el_r[jm].etiquette

    def libre(self, cote, alpha, phi):
        """La roue à phi ne touche pas la palette `cote` à alpha."""
        t = self.trace
        S_pal = rotation(self.seg_palettes[cote][0], t.ancre.axe, alpha)
        S_r = rotation(self.seg_roue[0], t.roue.axe, phi)
        # contour en facettes, inscrit dans le vrai : une marge de 0,1 µm
        # couvre largement la flèche des arcs
        return distance(S_r, S_pal)[0] > 1e-4

    def _balayer_palette(self, cote):
        """Phi(alpha) pour une palette, en suivant la dent de proche en proche
        depuis la pose dessinée."""
        A = self.alpha
        n = len(A)
        phi = np.full(n, np.nan)
        P = np.full((n, 2), np.nan)
        N = np.full((n, 2), np.nan)
        lab_a, lab_r = [""] * n, [""] * n
        i0 = int(np.argmin(np.abs(A)))
        for ordre in (range(i0, n), range(i0 - 1, -1, -1)):
            ref = 0.0
            for i in ordre:
                depart = ref
                # partir d'une position libre : reculer tant qu'il le faut
                for _ in range(int(self.pas_roue / RECUL) + 2):
                    if self.libre(cote, A[i], depart):
                        break
                    depart -= self.s * RECUL
                # le contact est presque toujours à quelques degrés : une
                # course courte réduit d'autant les paires sommet-segment ;
                # le premier contact trouvé sur une course courte l'est aussi
                # sur la longue
                r = self.contact_roue(cote, A[i], depart, course=COURSE_COURTE)
                if r is None:
                    r = self.contact_roue(cote, A[i], depart)
                if r is None:
                    continue
                phi[i], P[i], N[i], lab_a[i], lab_r[i] = r
                ref = phi[i] - self.s * RECUL
        return TableRoue(cote, phi, P, N, lab_a, lab_r)

    # -- fourchette ----------------------------------------------------------
    def centre_fourchette(self, alpha):
        """Angle du balancier qui place le centre de l'ellipse sur l'axe de la
        fourchette posée à alpha."""
        t = self.trace
        bag = t.ancre.element("baguette").points
        u = rotation(bag[-1], t.ancre.axe, alpha) - t.ancre.axe
        u = u / np.hypot(*u)
        c = next(e.centre for e in t.balancier.elements if e.est_arc and e.rayon > 0.1)
        rho = np.hypot(*(c - t.balancier.axe))
        o = t.ancre.axe - t.balancier.axe
        b, cc = o @ u, o @ o - rho ** 2
        tt = -b - np.sqrt(max(b * b - cc, 0.0))
        p = t.ancre.axe + tt * u - t.balancier.axe
        v = c - t.balancier.axe
        return float(np.arctan2(v[0] * p[1] - v[1] * p[0], v @ p))

    def contact_fourche(self, alpha, beta_ref, sens, course=40 * DEG):
        """Rotation du balancier depuis beta_ref (libre) jusqu'au contact de
        l'ellipse avec la fourchette posée à alpha."""
        t = self.trace
        _, P_f, el_f = self.seg_fourche
        P_e = self.seg_ellipse[1].tournees(t.balancier.axe, beta_ref)
        tt, P, n, jf, _ = premier_contact_exact(P_e, t.balancier.axe, sens,
                                                P_f.tournees(t.ancre.axe, alpha), course=course)
        if not np.isfinite(tt):
            return None
        return beta_ref + sens * tt, P, n, el_f[jf].etiquette

    def _bornes_fourche(self):
        A = self.alpha
        n = len(A)
        out = []
        for sens in (-1, 1):
            beta = np.full(n, np.nan)
            P = np.full((n, 2), np.nan)
            N = np.full((n, 2), np.nan)
            lab = [""] * n
            for i, a in enumerate(A):
                r = self.contact_fourche(a, self.centre_fourchette(a), sens)
                if r is not None:
                    beta[i], P[i], N[i], lab[i] = r
            out.append(TableFourche(beta, P, N, lab))
        return out

    # -- lecture -------------------------------------------------------------
    def indice(self, alpha):
        return int(np.clip(np.searchsorted(self.alpha, alpha), 1, len(self.alpha) - 1))

    def interp(self, tableau, alpha):
        """Interpolation linéaire d'un tableau de la grille ; nan si l'un des
        deux voisins l'est."""
        i = self.indice(alpha)
        a0, a1 = self.alpha[i - 1], self.alpha[i]
        w = (alpha - a0) / (a1 - a0)
        return (1 - w) * tableau[i - 1] + w * tableau[i]


class _Filtre:
    """Vue d'un corps réduite à certains de ses éléments."""

    def __init__(self, corps, garde):
        self.elements = [e for e in corps.elements if id(e) in garde]
