"""Simulation dynamique de l'échappement : roue, ancre et balancier.

Trois corps tournent autour d'axes fixes : la roue (phi), l'ancre (alpha) et
le balancier (beta). Ils interagissent par des contacts unilatéraux dont la
géométrie est lue dans les tables (tables.py), exprimés comme des JEUX en
angle :

  - dent-palette, pour chaque palette : g = s.(Phi(alpha) - phi), la rotation
    qui reste à la roue avant de toucher (s : sens de la roue) ;
  - ellipse-fourchette : g = beta - B_bas(alpha) et g = B_haut(alpha) - beta,
    seulement quand l'ellipse est près de la fourchette ;
  - fourchette-butées : g = alpha - alpha_min et g = alpha_max - alpha.

Efforts extérieurs : couple moteur constant sur la roue ; spiral linéaire et
amortissement visqueux sur le balancier, réglé sur le facteur de qualité (il
contient le frottement de l'air). Frottement de Coulomb aux contacts
dent-palette et ellipse-fourchette, régularisé à très basse vitesse de
glissement.

Méthode : par événements. Entre deux événements, les contacts actifs sont
des liaisons (multiplicateurs de Lagrange, effort normal >= 0) ; un contact
qui se ferme avec une vitesse d'approche est un CHOC (loi de Newton,
coefficient de restitution, impulsions simultanées sur les contacts déjà
fermés) ; un contact dont l'effort s'annule s'ouvre ; une dent qui quitte la
palette (fin de table) s'échappe.

Le bilan d'énergie est tenu poste par poste (travail moteur, balancier,
frottements, chocs, amortissement, énergie cinétique de la roue et de
l'ancre) ; sa fermeture sert de contrôle.
"""

import math
from bisect import bisect_right
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from .cinematique import Releve
from .parametres import Parametres
from .geometrie import rotation
from .tables import DEG, RECUL

#: vitesse de glissement (mm/s) en dessous de laquelle le frottement est régularisé
V_REG = 0.01
#: fenêtre (rad) autour de la ligne des centres où l'ellipse peut toucher la fourchette
FENETRE_FOURCHE = 70 * DEG
#: une vitesse de rebond plus faible que cette fraction de la vitesse d'impact
#: (ou que V_COLLE) referme le contact : on évite la cascade de micro-rebonds
V_COLLE = 1e-3          # rad/s en vitesse de jeu
#: rebonds en cascade : au-delà de N_REBONDS chocs sur un même contact, chacun
#: à moins de FENETRE_ZENON du précédent, ou sous V_ZENON (rad/s de vitesse de
#: jeu), le choc est traité comme mou
N_REBONDS = 6
FENETRE_ZENON = 2e-4
V_ZENON = 0.02
#: le balancier part à cet angle avant l'entrée de l'ellipse dans la fourchette
MARGE_ENTREE = 2.0 * DEG
#: tolérance de pénétration (rad) avant qu'un contact ouvert soit déclaré fermé
G_TOL = 1e-9
#: constante de temps de la stabilisation des liaisons (Baumgarte)
TAU_STAB = 2e-4


def _cr(a, b):
    return a[0] * b[1] - a[1] * b[0]


class TableLisse:
    """Table y(alpha) lissée par morceaux (splines cubiques), nan hors des
    morceaux ; points et normales interpolés linéairement."""

    def __init__(self, alpha, y, P=None, N=None, n_min=5):
        self.morceaux = []
        # bords de la grille (au-delà des butées : l'ancre n'y va pas)
        self.a_lo, self.a_hi = alpha[0], alpha[-1]
        ok = np.isfinite(y)
        i = 0
        while i < len(y):
            if not ok[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(y) and ok[j + 1]:
                j += 1
            if j - i + 1 >= n_min:
                a = alpha[i:j + 1]
                # spline d'interpolation : un lissage déformerait les passages
                # d'arête (congés de quelques µm), où le rapport de transmission
                # change vite mais réellement
                sp = CubicSpline(a, y[i:j + 1])
                m = {"lo": a[0], "hi": a[-1], "s": sp, "a": a,
                     # pour l'évaluation rapide, en flottants Python : nœuds
                     # et coefficients (c0.dx³ + c1.dx² + c2.dx + c3)
                     "x": a.tolist(), "c": sp.c.T.tolist()}
                if P is not None:
                    m["P"], m["N"] = P[i:j + 1], N[i:j + 1]
                    m["xy"] = np.c_[P[i:j + 1], N[i:j + 1]].tolist()
                self.morceaux.append(m)
            i = j + 1

    def morceau(self, a):
        for k, m in enumerate(self.morceaux):
            # un morceau qui touche le bout de la grille (une butée) se
            # prolonge au-delà : l'ancre peut y déborder d'un rien
            lo = -np.inf if m["lo"] <= self.a_lo else m["lo"]
            hi = np.inf if m["hi"] >= self.a_hi else m["hi"]
            if lo <= a <= hi:
                return k
        return None

    def marge_interieure(self, a, k):
        """Distance (rad) au bout du morceau k qui n'est pas une butée : le
        contact cesse quand elle s'annule (la dent s'échappe)."""
        m = self.morceaux[k]
        d = []
        if m["lo"] > self.a_lo:
            d.append(a - m["lo"])
        if m["hi"] < self.a_hi:
            d.append(m["hi"] - a)
        return min(d) if d else 1.0

    def proche(self, a):
        """Morceau qui contient a, ou à défaut le plus proche : un contact
        fermé reste évalué (par prolongement de la spline) quand un pas
        d'intégration déborde un peu de la table, avant que l'événement de fin
        de table ne soit localisé."""
        k = self.morceau(a)
        if k is not None or not self.morceaux:
            return k
        d = [min(abs(a - m["lo"]), abs(a - m["hi"])) for m in self.morceaux]
        return int(np.argmin(d))

    def eval(self, a, k=None):
        """(y, y', y'') en a sur le morceau k. Évaluation directe du polynôme
        de la spline (hors du morceau : prolongement du polynôme d'extrémité,
        comme CubicSpline) ; c'est le cœur de la boucle d'intégration."""
        k = self.morceau(a) if k is None else k
        if k is None:
            return None
        m = self.morceaux[k]
        x = m["x"]
        i = min(max(bisect_right(x, a) - 1, 0), len(x) - 2)
        c0, c1, c2, c3 = m["c"][i]
        d = a - x[i]
        return (((c0 * d + c1) * d + c2) * d + c3, (3 * c0 * d + 2 * c1) * d + c2,
                6 * c0 * d + 2 * c1)

    def point_xy(self, a, k):
        """Point et normale de contact interpolés, en flottants :
        (px, py, nx, ny), normale unitaire."""
        m = self.morceaux[k]
        x = m["x"]
        a = min(max(a, x[0]), x[-1])
        i = min(max(bisect_right(x, a) - 1, 0), len(x) - 2)
        u = (a - x[i]) / (x[i + 1] - x[i])
        p0, p1 = m["xy"][i], m["xy"][i + 1]
        px, py, nx, ny = (p0[j] + u * (p1[j] - p0[j]) for j in range(4))
        r = math.hypot(nx, ny)
        return px, py, nx / r, ny / r

    def point(self, a, k):
        px, py, nx, ny = self.point_xy(a, k)
        return np.array([px, py]), np.array([nx, ny])


@dataclass
class Contact:
    nom: str
    genre: str               # "roue", "fourche", "butee"
    e: float
    signe: int = 1           # fourche : +1 borne basse, -1 borne haute ; butée : +1 min, -1 max
    table: object = None
    valeur: float = 0.0      # butée : angle d'ancre de la butée


@dataclass
class Resultat:
    """Une simulation : historique et bilan."""
    t: np.ndarray = None
    y: np.ndarray = None                 # (n, 8) : phi, alpha, beta, vitesses, pertes cumulées
    couple_balancier: np.ndarray = None  # couple de la fourchette sur le balancier (hors chocs)
    actifs: list = field(default_factory=list)
    chocs: list = field(default_factory=list)       # (t, contact, énergie perdue J,
                                                     #  impulsion sur le balancier N.m.s)
    bilan: dict = field(default_factory=dict)
    fenetres: list = field(default_factory=list)     # intervalles de temps où l'échappement agit
    valide: bool = True
    message: str = ""


class Simulateur:
    """Dynamique d'un tracé (tables + relevé) avec des paramètres donnés."""

    def __init__(self, rel: Releve, par: Parametres, trace_ev=False):
        self.rel, self.par = rel, par
        self.trace_ev = trace_ev
        T = rel.T
        self.T = T
        t = T.trace
        self.s = t.sens_roue
        self.p = t.pas
        self.O_r, self.O_a, self.O_b = t.roue.axe, t.ancre.axe, t.balancier.axe
        self.J = np.array([par.inertie_roue + par.inertie_rouage, par.inertie_ancre,
                           par.inertie_balancier])
        self.Minv = 1.0 / self.J
        self.k, self.c = par.raideur, par.amortissement
        self.beta0 = T.beta_ldc
        self.contacts = []
        for cote in ("entrée", "sortie"):
            tab = T.roue[cote]
            self.contacts.append(Contact(f"dent-palette {cote}", "roue", par.e_chute,
                                         table=TableLisse(T.alpha, tab.phi, tab.P, tab.N)))
        for nom, signe, tab in (("ellipse-flanc bas", 1, T.fourche_bas),
                                ("ellipse-flanc haut", -1, T.fourche_haut)):
            self.contacts.append(Contact(nom, "fourche", par.e_fourchette, signe,
                                         table=TableLisse(T.alpha, tab.beta, tab.P, tab.N)))
        self.contacts.append(Contact("butée min", "butee", par.e_butee, 1, valeur=T.a_min))
        self.contacts.append(Contact("butée max", "butee", par.e_butee, -1, valeur=T.a_max))

    # -- géométrie des contacts ---------------------------------------------
    def jeu(self, c, q):
        """Jeu du contact c (rad), ou None s'il ne peut pas agir ici."""
        phi, a, b = q
        if c.genre == "butee":
            return c.signe * (a - c.valeur)
        if c.genre == "roue":
            v = c.table.eval(a)
            if v is None:
                return None
            return np.mod(self.s * (v[0] - phi) + RECUL, self.p) - RECUL
        if abs(b - self.beta0) > FENETRE_FOURCHE:
            return None
        v = c.table.eval(a)
        if v is None:
            return None
        return c.signe * (b - v[0])

    def jacobien(self, c, q, v):
        """(gradient du jeu, terme de courbure h, vecteur de frottement f par
        unité de multiplicateur) au contact c, en tuples de flottants."""
        a = q[1]
        if c.genre == "butee":
            return (0.0, float(c.signe), 0.0), 0.0, (0.0, 0.0, 0.0)
        k = c.table.proche(a)
        y, y1, y2 = c.table.eval(a, k)
        px, py, nx, ny = c.table.point_xy(a, k)
        tx, ty = -ny, nx
        mu = self.par.mu
        ax, ay = px - self.O_a[0], py - self.O_a[1]
        va = v[1]
        if c.genre == "roue":
            s = self.s
            grad = (-s, s * y1, 0.0)
            h = s * y2 * va * va
            # effort normal physique N = -s.lambda / ((P - O_r) x n)
            rx, ry = px - self.O_r[0], py - self.O_r[1]
            N_par_lambda = -s / (rx * ny - ry * nx)
            vt = (-v[0] * ry + va * ay) * tx + (v[0] * rx - va * ax) * ty
            F = -mu * N_par_lambda * math.tanh(vt / V_REG)
            fx, fy = F * tx, F * ty
            return grad, h, (rx * fy - ry * fx, -(ax * fy - ay * fx), 0.0)
        # fourchette
        sg = c.signe
        grad = (0.0, -sg * y1, float(sg))
        h = -sg * y2 * va * va
        bx, by = px - self.O_b[0], py - self.O_b[1]
        K_par_lambda = sg / (bx * ny - by * nx)
        vb = v[2]
        vt = (-vb * by + va * ay) * tx + (vb * bx - va * ax) * ty
        F = -mu * K_par_lambda * math.tanh(vt / V_REG)
        fx, fy = F * tx, F * ty
        return grad, h, (0.0, -(ax * fy - ay * fx), bx * fy - by * fx)

    # -- dynamique liée --------------------------------------------------------
    def efforts_ext(self, q, v):
        return np.array([self.s * self.par.couple_roue, 0.0,
                         -self.k * (q[2] - self.beta0) - self.c * v[2]])

    def resoudre(self, q, v, actifs):
        """Accélérations et multiplicateurs pour un ensemble de contacts
        fermés. Rend (a, lambdas, colonnes de forces généralisées)."""
        a, lam, A, _ = self._resoudre(q, v, actifs)
        return a, lam, A

    def _resoudre(self, q, v, actifs):
        """Comme resoudre, et la puissance dissipée par frottement."""
        Q = self.efforts_ext(q, v)
        Minv = self.Minv
        if not actifs:
            return Q * Minv, np.zeros(0), np.zeros((3, 0)), 0.0
        W, A, r, F = [], [], [], []
        for c in actifs:
            g = self.jeu(c, q)
            grad, h, f = self.jacobien(c, q, v)
            gp = grad[0] * v[0] + grad[1] * v[1] + grad[2] * v[2]
            W.append(grad)
            A.append((grad[0] + f[0], grad[1] + f[1], grad[2] + f[2]))
            F.append(f)
            r.append(-h - 2 * gp / TAU_STAB - (g or 0.0) / TAU_STAB ** 2)
        W, A, r = np.array(W), np.array(A).T, np.array(r)
        Mm = W @ (Minv[:, None] * A)
        b = r - W @ (Minv * Q)
        if len(actifs) == 1:
            lam = b / Mm[0, 0] if Mm[0, 0] != 0 else np.zeros(1)
        else:
            try:
                lam = np.linalg.solve(Mm, b)
            except np.linalg.LinAlgError:
                lam = np.linalg.lstsq(Mm, b, rcond=None)[0]
        acc = Minv * (Q + A @ lam)
        frott = -sum(lam[k] * (f[0] * v[0] + f[1] * v[1] + f[2] * v[2]) for k, f in enumerate(F))
        return acc, lam, A, float(frott)

    def rhs(self, t, y, actifs):
        q, v = y[:3], y[3:6]
        a, _, _, frott = self._resoudre(q, v, actifs)
        return np.concatenate([v, a, [self.c * v[2] ** 2, frott]])

    # -- chocs ----------------------------------------------------------------
    def choc(self, q, v, impact, actifs):
        """Vitesses après le choc du contact `impact`, les contacts `actifs`
        restant fermés s'ils le peuvent. Rend (v+, contacts fermés après)."""
        Minv = 1.0 / self.J
        grads = {}
        for c in [impact] + list(actifs):
            grad, _, f = self.jacobien(c, q, v)
            grad, f = np.array(grad), np.array(f)
            # frottement d'impulsion au sens du glissement d'avant choc
            grads[c.nom] = (grad, grad + f)
        v_imp = grads[impact.nom][0] @ v
        cible_imp = -impact.e * v_imp
        meilleur = None
        n = len(actifs)
        for masque in range(2 ** n - 1, -1, -1):
            garde = [actifs[i] for i in range(n) if masque >> i & 1]
            jeu_c = [impact] + garde
            W = np.array([grads[c.nom][0] for c in jeu_c])
            A = np.array([grads[c.nom][1] for c in jeu_c]).T
            cible = np.array([cible_imp] + [0.0] * len(garde))
            Mm = W @ (Minv[:, None] * A)
            try:
                Pi = np.linalg.solve(Mm, cible - W @ v)
            except np.linalg.LinAlgError:
                continue
            if np.any(Pi < -1e-14):
                continue
            vp = v + Minv * (A @ Pi)
            libres = [c for c in actifs if c not in garde]
            if any(grads[c.nom][0] @ vp < -1e-9 for c in libres):
                continue
            meilleur = (vp, garde)
            break
        if meilleur is None:                  # sans frottement en dernier recours
            W = np.array([grads[impact.nom][0]])
            Pi = (cible_imp - v_imp) / (W @ (Minv * W[0]))
            meilleur = (v + Minv * W[0] * Pi[0], [])
        vp, garde = meilleur
        rebond = grads[impact.nom][0] @ vp
        fermes = list(garde)
        if rebond < max(V_COLLE, 0.01 * abs(v_imp)):
            fermes.append(impact)
        return vp, fermes

    def energie_cinetique(self, v):
        return 0.5 * float(np.sum(self.J * v * v))

    # -- intégration ----------------------------------------------------------
    def etat_initial(self, amplitude_deg):
        """État de départ : le balancier arrive sur la fourchette.

        C'est l'ellipse qui lance le dégagement : elle arrive avec la vitesse
        angulaire que lui donne son amplitude, theta' = omega0.sqrt(Theta0² -
        theta²). Le balancier est donc posé un peu avant l'entrée de l'ellipse
        dans la fourchette (MARGE_ENTREE), avec exactement cette vitesse ;
        l'ancre est sur la butée où l'a laissée l'alternance précédente, la
        roue au repos sur la palette d'entrée."""
        # départ au repos sur la palette d'entrée : la première alternance
        # simulée est celle où travaille l'entrée
        alt = next(a for a in self.rel.alternances if a.palette == "entrée")
        prec = next(a for a in self.rel.alternances if a is not alt)
        phi = prec.phi[-1]
        # la roue est posée exactement sur la table de la palette qui la
        # retient (jeu nul), et non sur la table brute du relevé
        c = next(c for c in self.contacts if c.nom == f"dent-palette {prec.palette_suivante}")
        v = c.table.eval(prec.arrivee)
        if v is not None:
            g = np.mod(self.s * (v[0] - phi) + RECUL, self.p) - RECUL
            phi = phi + self.s * g
        theta_entree = alt.s_b * np.radians(alt.mesures["entree_ellipse_deg"])
        theta0 = theta_entree - alt.s_b * MARGE_ENTREE
        Th = np.radians(amplitude_deg)
        w0 = 2 * np.pi * self.par.frequence
        vitesse = alt.s_b * w0 * np.sqrt(max(Th ** 2 - theta0 ** 2, 0.0))
        self.theta_depart, self.sens_depart = theta0, alt.s_b
        q = np.array([phi, prec.arrivee, self.beta0 + theta0])
        return np.concatenate([q, [0.0, 0.0, vitesse], np.zeros(2)])

    def actifs_initiaux(self, y):
        q = y[:3]
        return [c for c in self.contacts if (g := self.jeu(c, q)) is not None and abs(g) < 1e-9]

    def simuler(self, amplitude_deg=None, n_oscillations=1, rtol=1e-8, atol=1e-11,
                max_evenements=4000):
        """Simule n oscillations complètes, depuis l'arrivée de l'ellipse sur
        la fourchette jusqu'à son retour au même point."""
        par = self.par
        amp = par.amplitude_deg if amplitude_deg is None else amplitude_deg
        y = self.etat_initial(amp)
        actifs = self.actifs_initiaux(y)
        res = Resultat()
        T_libre = 2 * np.pi / (2 * np.pi * par.frequence)
        t, t_fin = 0.0, n_oscillations * T_libre * 1.5
        ts, ys, couples, act_hist = [np.array([0.0])], [y[None]], [], []
        retours = 0
        # le retour au point de départ n'est guetté qu'à partir des trois
        # quarts de période qui suivent le départ ou le retour précédent (le
        # balancier est alors vers l'extrémité opposée : la fonction ne saute
        # pas de signe). Sans ce délai, le segment qui repart exactement sur
        # le point de départ, à l'arrondi près du mauvais côté, compterait le
        # même passage une seconde fois.
        t_garde = 0.75 * T_libre
        E_choc = 0.0
        meme_instant, t_prec = 0, -1.0
        rebonds = {}
        for _ in range(max_evenements):
            fermes = list(actifs)
            ouverts = [c for c in self.contacts if c not in fermes]

            def ev_jeu(c):
                # seuil légèrement négatif : un contact qui vient de s'ouvrir
                # part d'un jeu nul, à l'erreur d'arrondi près ; sans cette
                # marge, une pénétration naissante ne changerait pas de signe
                # et la dent traverserait la palette.
                # Pour la roue, la dent visée est fixée pour toute la phase
                # (pas de réduction modulo le pas) : une roue libre et rapide
                # peut avancer de plus d'un degré entre deux pas d'intégration
                # sans que le contact soit manqué.
                # Le seuil est placé sous le jeu de départ : juste après un
                # rebond, le jeu vaut le seuil précédent à l'arrondi près, et
                # une nouvelle pénétration ne changerait pas de signe.
                if c.genre == "roue":
                    v0 = c.table.eval(y[1])
                    decalage = None if v0 is None else \
                        (np.mod(self.s * (v0[0] - y[0]) + RECUL, self.p) - RECUL) \
                        - self.s * (v0[0] - y[0])
                    g0 = None if v0 is None else self.s * (v0[0] - y[0]) + decalage
                    seuil = -G_TOL if g0 is None else min(g0, 0.0) - G_TOL

                    def f(tt, yy):
                        v = c.table.eval(yy[1])
                        if v is None or decalage is None:
                            return 1.0
                        return self.s * (v[0] - yy[0]) + decalage - seuil
                    f.terminal, f.direction = True, -1
                    return f

                g0 = self.jeu(c, y[:3])
                seuil = -G_TOL if g0 is None else min(g0, 0.0) - G_TOL

                def f(tt, yy):
                    g = self.jeu(c, yy[:3])
                    return 1.0 if g is None else g - seuil
                f.terminal, f.direction = True, -1
                return f

            def ev_effort(k):
                def f(tt, yy):
                    _, lam, _ = self.resoudre(yy[:3], yy[3:6], fermes)
                    return lam[k]
                f.terminal, f.direction = True, -1
                return f

            def ev_domaine(c):
                k0 = c.table.proche(y[1])
                def f(tt, yy):
                    return c.table.marge_interieure(yy[1], k0)
                f.terminal, f.direction = True, -1
                return f

            def ev_extremite(tt, yy, t_garde=t_garde):
                # fin d'une oscillation : le balancier repasse par son angle
                # de départ, dans le même sens
                if tt < t_garde:
                    return -1.0
                return self.sens_depart * (yy[2] - self.beta0 - self.theta_depart)
            ev_extremite.terminal, ev_extremite.direction = True, 1

            evs = [ev_jeu(c) for c in ouverts] + [ev_effort(k) for k in range(len(fermes))]
            evs += [ev_domaine(c) for c in fermes if c.genre != "butee"]
            evs.append(ev_extremite)
            sol = solve_ivp(lambda tt, yy: self.rhs(tt, yy, fermes), (t, t_fin), y,
                            method="DOP853", rtol=rtol, atol=atol, events=evs,
                            max_step=T_libre / 200)
            ts.append(sol.t[1:])
            ys.append(sol.y[:, 1:].T)
            act_hist.append((t, sol.t[-1], [c.nom for c in fermes]))
            t, y = sol.t[-1], sol.y[:, -1].copy()
            if sol.status == 0:
                break
            # quel événement ?
            k_ev = next(i for i, te in enumerate(sol.t_events) if len(te))
            n_o = len(ouverts)
            n_f = len(fermes)
            q, v = y[:3], y[3:6]
            if self.trace_ev:
                noms = [c.nom for c in ouverts] + [f"effort {c.nom}" for c in fermes] + \
                    [f"fin de table {c.nom}" for c in fermes if c.genre != "butee"] + ["extrémité"]
                print(f"  t={1000 * t:9.4f} ms  {noms[k_ev]:30s} fermés={[c.nom for c in fermes]}"
                      f"  phi={np.degrees(y[0]):9.3f} alpha={np.degrees(y[1]):7.3f}"
                      f" theta={np.degrees(y[2] - self.beta0):8.2f} dphi={y[3]:9.2f}")
            meme_instant = meme_instant + 1 if abs(t - t_prec) < 1e-12 else 0
            t_prec = t
            if k_ev < n_o:                                   # un contact se ferme
                c = ouverts[k_ev]
                grad, _, _ = self.jacobien(c, q, v)
                grad = np.array(grad)
                if grad @ v < -1e-12:
                    E0 = self.energie_cinetique(v)
                    e_sauve = c.e
                    # Rebonds en cascade (effet Zénon) : sous un effort qui
                    # reste appliqué, les rebonds se resserrent sans fin. Après
                    # quelques rebonds rapprochés sur le même contact, ou sous
                    # une vitesse d'impact négligeable, le choc devient mou et
                    # le contact se ferme.
                    serie = rebonds.get(c.nom, (0, -1.0))
                    n_reb = serie[0] + 1 if t - serie[1] < FENETRE_ZENON else 1
                    rebonds[c.nom] = (n_reb, t)
                    if meme_instant > 3 or n_reb > N_REBONDS or abs(grad @ v) < V_ZENON:
                        c.e = 0.0
                    vp, nouveaux = self.choc(q, v, c, fermes)
                    c.e = e_sauve
                    E1 = self.energie_cinetique(vp)
                    # impulsion angulaire reçue par le balancier (N.m.s)
                    res.chocs.append((t, c.nom, E0 - E1, self.J[2] * (vp[2] - v[2])))
                    E_choc += E0 - E1
                    y[3:6] = vp
                    actifs = nouveaux
                else:
                    actifs = fermes + [c]
            elif k_ev < n_o + n_f:                           # un effort s'annule
                actifs = [c for i, c in enumerate(fermes) if i != k_ev - n_o]
            elif k_ev < len(evs) - 1:                        # la dent s'échappe
                roues = [c for c in fermes if c.genre != "butee"]
                c = roues[k_ev - n_o - n_f]
                actifs = [x for x in fermes if x is not c]
            else:                                            # retour au point de départ
                if t >= t_garde:
                    retours += 1
                    t_garde = t + 0.75 * T_libre
                    if retours >= n_oscillations:
                        break
                actifs = fermes
        else:
            res.valide, res.message = False, "trop d'événements"
        res.t = np.concatenate(ts)
        res.y = np.concatenate(ys)
        res.actifs = act_hist
        res.bilan = self._bilan(res, E_choc, amp, max(retours, 1))
        # recul de la roue au dégagement, alternance par alternance : la roue
        # repoussée par la palette avant de pouvoir avancer
        reculs = []
        for a, b in fenetres_action(res):
            k = (res.t >= a) & (res.t <= b)
            if k.any():
                phi = res.y[k, 0]
                reculs.append(float(np.degrees(max(0.0, np.max(-self.s * (phi - phi[0]))))))
        res.bilan["recul_roue_deg"] = reculs
        return res

    def _amplitude(self, y):
        """Amplitude (°) d'un état du balancier : celle de l'oscillateur libre
        qui a la même énergie."""
        w0 = 2 * np.pi * self.par.frequence
        th = y[2] - self.beta0
        return float(np.degrees(np.sqrt(th ** 2 + (y[5] / w0) ** 2)))

    def _bilan(self, res, E_choc, amp, n_osc=1):
        par = self.par
        y0, y1 = res.y[0], res.y[-1]
        W = self.s * par.couple_roue * (y1[0] - y0[0])
        Eb = lambda y: 0.5 * self.J[2] * y[5] ** 2 + 0.5 * self.k * (y[2] - self.beta0) ** 2
        Ec = lambda y: 0.5 * self.J[0] * y[3] ** 2 + 0.5 * self.J[1] * y[4] ** 2
        amort, frott = y1[6] - y0[6], y1[7] - y0[7]
        dEb = Eb(y1) - Eb(y0)
        recu = dEb + amort
        b = {
            "travail_moteur": W, "energie_recue_balancier": recu,
            "variation_energie_balancier": dEb, "amortissement": amort,
            "frottements": frott, "chocs": E_choc,
            "energie_cinetique_roue_ancre": Ec(y1) - Ec(y0),
            "rendement": recu / W if W > 0 else np.nan,
            "periode": res.t[-1] - res.t[0],         # durée simulée (n oscillations)
            "oscillations": n_osc,
            "amplitude_depart_deg": self._amplitude(y0),
            "amplitude_fin_deg": self._amplitude(y1),
        }
        b["ecart_bilan"] = W - (dEb + amort + frott + E_choc + b["energie_cinetique_roue_ancre"])
        # Marche : avance de PHASE de l'oscillateur comparée à celle du
        # balancier libre amorti. Comparer des instants de passage à un angle
        # fixe serait faux dès que l'amplitude varie sur l'oscillation (le
        # passage arrive alors plus tôt ou plus tard sans que la marche change :
        # 3,5 s/j d'écart pour 0,8° de perte d'amplitude). La phase de
        # theta = A.exp(-g.t).cos(psi) se lit sans connaître A :
        # psi = atan2(-(theta' + g.theta)/omega_d, theta).
        w0 = 2 * np.pi * par.frequence
        g = w0 / (2 * par.facteur_qualite)
        wd = w0 * np.sqrt(1 - 1 / (4 * par.facteur_qualite ** 2))
        psi = lambda y: np.arctan2(-(y[5] + g * (y[2] - self.beta0)) / wd, y[2] - self.beta0)
        dpsi = np.mod(psi(y1) - psi(y0) + np.pi, 2 * np.pi) - np.pi
        b["marche_s_j"] = 86400.0 * ((2 * np.pi * n_osc + dpsi) / (wd * b["periode"]) - 1.0)
        return b


def simuler(rel, par=None, amplitude_deg=None):
    return Simulateur(rel, par or Parametres()).simuler(amplitude_deg)


def amplitude_equilibre(rel, par, depart_deg=None, tol_deg=0.01, n_max=8, rappel=None):
    """Amplitude stable du balancier : celle où l'énergie que l'échappement
    lui rend sur une oscillation compense exactement ses pertes (facteur de
    qualité). Au-dessus, l'amplitude décroît d'une oscillation à l'autre ;
    au-dessous, elle croît : l'équilibre est stable.

    Méthode de la sécante sur d(A) = A_après une oscillation - A, chaque
    évaluation étant une oscillation simulée. Rend (amplitude °, dernier
    Resultat, liste des (A, d(A)) évalués)."""
    A0 = par.amplitude_deg if depart_deg is None else depart_deg
    points = []

    def d(A):
        r = Simulateur(rel, par).simuler(amplitude_deg=A)
        points.append((A, r.bilan["amplitude_fin_deg"] - A))
        if rappel is not None:
            rappel(A, points[-1][1])
        return points[-1][1], r

    f0, r = d(A0)
    if abs(f0) < tol_deg:
        return A0, r, points
    # second point du côté où l'amplitude va
    A1 = A0 + float(np.clip(10 * f0, -30, 30))
    f1, r = d(A1)
    for _ in range(n_max):
        if abs(f1) < tol_deg or f1 == f0:
            break
        A2 = A1 - f1 * (A1 - A0) / (f1 - f0)
        A2 = float(np.clip(A2, 90.0, 360.0))
        A0, f0 = A1, f1
        A1 = A2
        f1, r = d(A1)
    return A1, r, points


if __name__ == "__main__":
    import sys
    import time
    from pathlib import Path
    from .cinematique import releve
    from .dxf_import import lire_trace
    chemin = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("simulation.dxf"))
    rel = releve(lire_trace(chemin))
    t0 = time.time()
    r = simuler(rel)
    print(f"simulation : {time.time() - t0:.1f} s, {len(chocs_significatifs(r))} chocs, {r.message}")
    for k, v in r.bilan.items():
        print(f"  {k:32s} {v:.6g}")
    for c in chocs_significatifs(r)[:40]:
        print(f"  choc à {1000 * c[0]:8.3f} ms : {c[1]:22s} perte {1e9 * c[2]:.3f} nJ")


# ---------------------------------------------------------------------------
# Exploitation d'une simulation
# ---------------------------------------------------------------------------

def actifs_a(res, t):
    """Noms des contacts fermés à l'instant t."""
    for t0, t1, noms in res.actifs:
        if t0 <= t <= t1:
            return noms
    return []


def fenetres_action(res, marge=3e-4, seuil=1e-3):
    """Intervalles de temps où l'ancre bouge : là où l'échappement agit.
    Un par alternance en principe."""
    bouge = np.abs(res.y[:, 4]) > seuil
    fen = []
    i = 0
    while i < len(bouge):
        if not bouge[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(bouge) and (bouge[j + 1] or res.t[j + 1] - res.t[j] < marge):
            j += 1
        fen.append((max(res.t[0], res.t[i] - marge), min(res.t[-1], res.t[j] + marge)))
        i = j + 1
    # fusionne les fenêtres très proches (rebonds)
    out = []
    for a, b in fen:
        if out and a - out[-1][1] < 2e-3:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def etat_a(res, t):
    """Positions (phi, alpha, beta) interpolées à l'instant t."""
    return np.array([np.interp(t, res.t, res.y[:, k]) for k in range(3)])


def couple_fourchette(sim, res):
    """Couple de la fourchette sur le balancier aux instants calculés, hors
    chocs : (t, theta depuis la ligne des centres, couple N.m)."""
    noms = {c.nom: c for c in sim.contacts}
    T, TH, M = [], [], []
    for t, y in zip(res.t, res.y):
        act = [noms[n] for n in actifs_a(res, t)]
        fourche = [k for k, c in enumerate(act) if c.genre == "fourche"]
        if not fourche:
            continue
        q, v = y[:3], y[3:6]
        try:
            _, lam, A = sim.resoudre(q, v, act)
        except Exception:
            continue
        m = sum(lam[k] * A[2, k] for k in fourche)
        T.append(t)
        TH.append(q[2] - sim.beta0)
        M.append(m)
    return np.array(T), np.array(TH), np.array(M)


def decollement(sim, res):
    """Décollement de la dent après le dégagement, alternance par alternance.

    À la fin du dégagement la palette, menée par le balancier, s'efface plus
    vite que la roue ne peut la suivre : la dent quitte la palette, puis la
    rattrape et retombe sur le plan d'impulsion. On mesure, pour chaque
    alternance, le plus grand écart entre la dent et la palette (µm, mesuré à
    la tête de la roue), la durée du vol libre, et l'endroit où la dent
    retouche le plan d'impulsion, en % de sa longueur depuis l'arête de repos.
    """
    T = sim.T
    trace = T.trace
    Rt = trace.rayon_tete
    out = []
    # une fenêtre d'action par alternance, dans l'ordre simulé : l'entrée
    # d'abord (départ au repos sur l'entrée)
    fenetres = fenetres_action(res)
    for n, alt in enumerate(sorted(sim.rel.alternances, key=lambda a: a.palette != "entrée")):
        if n >= len(fenetres):
            break
        f0, f1 = fenetres[n]
        c = next(x for x in sim.contacts if x.nom == f"dent-palette {alt.palette}")
        imp = trace.ancre.element(f"plan d'impulsion {alt.palette}").points
        arete = trace.ancre.element(f"arête de repos {alt.palette}").points.mean(0)
        L, X = (imp[0], imp[-1]) if np.hypot(*(imp[0] - arete)) < np.hypot(*(imp[-1] - arete)) \
            else (imp[-1], imp[0])
        u = X - L
        a_deg = alt.evenements["fin_degagement"]
        s_a = np.sign(alt.arrivee - alt.depart)
        # instants où la dent est libre alors que l'ancre a dépassé la fin du
        # dégagement et que la dent n'a pas encore échappé
        libre = []
        for t, y in zip(res.t, res.y):
            if t < f0 or t > f1:
                continue
            if (y[1] - a_deg) * s_a < -1e-9 or (y[1] - alt.evenements["chute"]) * s_a > 0:
                continue
            if c.nom in actifs_a(res, t):
                if libre:
                    break
                continue
            g = sim.jeu(c, y[:3])
            if g is not None:
                libre.append((t, g, y))
        if not libre:
            out.append({"palette": alt.palette, "jeu_max_um": 0.0, "duree_ms": 0.0,
                        "reprise_pct": float("nan")})
            continue
        t_fin, _, y_fin = libre[-1]
        k = c.table.proche(y_fin[1])
        P, _ = c.table.point(y_fin[1], k)
        P_loc = rotation(P, trace.ancre.axe, -y_fin[1])
        out.append({"palette": alt.palette,
                    "jeu_max_um": 1000 * Rt * max(g for _, g, _ in libre),
                    "duree_ms": 1000 * (t_fin - libre[0][0]),
                    "reprise_pct": 100 * float((P_loc - L) @ u / (u @ u))})
    return out


#: énergie en deçà de laquelle un choc n'est pas compté (J) : en passant
#: d'une facette à l'autre des arcs échantillonnés, un contact s'ouvre parfois
#: un instant et se referme à vitesse quasi nulle, sans rien dissiper
CHOC_NEGLIGEABLE = 1e-12


def chocs_significatifs(res):
    return [c for c in res.chocs if c[2] > CHOC_NEGLIGEABLE]


def couple_fonction(sim, res):
    """Couple sur l'axe du balancier pendant toute la fonction, alternance par
    alternance, signé par l'ÉCHANGE D'ÉNERGIE : positif quand l'ancre est
    motrice (elle donne de l'énergie au balancier : impulsion), négatif quand
    le balancier est moteur (il fournit le travail du dégagement), nul quand
    l'ellipse ne touche pas la fourchette.

    Rend une liste, une entrée par fenêtre d'action : dict avec t (s, depuis
    le début de la fenêtre), couple (N.m), theta (rad, depuis la ligne des
    centres) et les chocs reçus par le balancier [(t, impulsion N.m.s signée
    de la même façon, contact)]."""
    noms = {c.nom: c for c in sim.contacts}
    out = []
    for a, b in fenetres_action(res):
        k = (res.t >= a) & (res.t <= b)
        T, C, TH = [], [], []
        for t, y in zip(res.t[k], res.y[k]):
            act = [noms[n] for n in actifs_a(res, t)]
            fourche = [i for i, c in enumerate(act) if c.genre == "fourche"]
            m = 0.0
            if fourche:
                try:
                    _, lam, A = sim.resoudre(y[:3], y[3:6], act)
                    m = float(sum(lam[i] * A[2, i] for i in fourche))
                except Exception:
                    m = 0.0
            T.append(t - a)
            C.append(m * np.sign(y[5]))
            TH.append(y[2] - sim.beta0)
        chocs = []
        for c in chocs_significatifs(res):
            if a <= c[0] <= b:
                v = np.interp(c[0], res.t, res.y[:, 5])
                chocs.append((c[0] - a, c[3] * np.sign(v), c[1]))
        out.append({"t": np.array(T), "couple": np.array(C), "theta": np.array(TH),
                    "chocs": chocs, "debut": a})
    return out
