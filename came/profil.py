"""Profil de coeur discretise : espace de conception pour une optimisation.

Pourquoi changer de representation
----------------------------------
La synthese analytique (core.py, chrono.py) **impose** le couple constant et
en deduit le profil par integration : le profil est une sortie, pas une
entree. Il ne reste que 4 ou 5 parametres (pente du plan, pivot, k, T0). On ne
peut donc pas arbitrer entre des objectifs contradictoires (ondulation du
couple, plage morte, effort, encombrement) : on ne peut que rejeter les
solutions qui ne conviennent pas.

Ici le profil devient la variable de conception. La bonne representation n'est
ni le rayon polaire ni une suite de points, mais le **rayon de courbure en
fonction de l'angle de normale** :

    rho(beta) = H(beta) + H''(beta)

parce que les deux proprietes dures du probleme y sont triviales :

  - **convexite** : rho(beta) >= 0, une simple borne sur chaque variable.
    C'est indispensable ici : le marteau a une face plane, il ne touche que
    l'enveloppe convexe du coeur ;
  - **fermeture** du contour : les harmoniques de rang 1 de rho doivent etre
    nulles. Il suffit de ne pas les mettre dans le jeu de variables.

En developpant rho en serie de Fourier

    rho(beta) = a0 + somme_{n>=2} [ an.cos(n.beta) + bn.sin(n.beta) ]

la fonction d'appui se reconstruit analytiquement (pas de derivee numerique) :

    H(beta)  = a0 + somme_{n>=2} [an.cos + bn.sin]/(1-n^2) + cx.cos + cy.sin
    s(beta)  = H'(beta)

Les deux termes cx, cy sont les harmoniques de rang 1 de H : ils ne changent
pas la forme, ils **deplacent le centre de rotation** dans le profil. Ce sont
donc exactement les deux degres de liberte qui comptent pour une came, et ils
sont separes de la forme.

Jeu de variables : a0, a2..aM, b2..bM, cx, cy  (2M variables pour M
harmoniques). Tout jeu respectant rho >= 0 donne un contour ferme convexe.
"""

import math

from . import noyau as _NOYAU
from .core import _tourner_autour as _tourner
from .core import contour_plan, dot, norm, rot, sgn, t_of, u_of

TAU = 2.0 * math.pi


class Profil:
    """Contour convexe ferme decrit par son rayon de courbure rho(beta).

    rho est une **serie de Fourier a bande limitee** :

        rho(beta) = a0 + somme_{n=2..M} [ an.cos(n.beta) + bn.sin(n.beta) ]

    et c'est exactement la courbure de la courbe reconstruite, donc la
    contrainte `rho >= r_outil` porte sur la geometrie que l'on usinera.

    (Une premiere version stockait rho par cellules : les valeurs de cellules
    etaient bien positives, mais la reconstruction a bande limitee sonnait
    autour du Dirac que represente un plat du profil, et la courbure reelle
    descendait a -1.2 mm. Un plat exact n'existe pas dans un espace a bande
    limitee ; c'est de toute facon ce que dit la contrainte d'usinabilite.)

    Le rang 1 est absent du jeu de variables : c'est exactement la condition
    de fermeture du contour. Les deux termes cx, cy sont les harmoniques de
    rang 1 de H : ils ne changent pas la forme, ils deplacent le centre de
    rotation.
    """

    def __init__(self, a0, an, bn, cx=0.0, cy=0.0, plat=0.0, beta0=0.0):
        self.a0 = float(a0)
        self.an = list(an)          # rangs 2..M
        self.bn = list(bn)
        self.cx, self.cy = float(cx), float(cy)
        # `plat` et `beta0` ne sont **pas** portes par le modele : ils sont
        # conserves a titre d'information (longueur du plat visee et normale
        # de repos), et `profil.longueur_plat()` mesure ce que le profil en
        # rend reellement.
        #
        # Un plat exact demanderait un Dirac de courbure a beta0. Il serait
        # admissible en soi - toute fonction d'appui periodique a une courbure
        # de rang 1 nul, donc fermee - mais la decomposition
        # `H = H_plat + H_bande` ne l'est pas : elle exigerait que le coeur
        # contienne la forme du plat comme sommande de Minkowski, ce que rien
        # ne garantit. Essai fait : le profil sortait a 0.24..1.41 mm au lieu
        # de 0.30..3.00, et sa courbure totale devenait negative.
        #
        # Le porter correctement demande de reparametrer l'espace
        # d'optimisation : rho = Dirac + partie lissee, avec le rang 1 de la
        # partie lissee libre pour compenser celui du Dirac - alors que c'est
        # aujourd'hui l'absence de rang 1 qui assure la fermeture.
        self.plat = max(0.0, float(plat))
        self.beta0 = float(beta0)
        self.M = len(self.an) + 1
        self._w = [1.0 / (1.0 - n * n) for n in range(2, self.M + 1)]

    # -- variables d'optimisation -----------------------------------------
    def vector(self):
        return [self.a0] + self.an + self.bn + [self.cx, self.cy]

    @property
    def gabarit(self):
        """(plat, beta0) : ce qui ne fait pas partie des variables."""
        return (self.plat, self.beta0)

    @staticmethod
    def from_vector(v, plat=0.0, beta0=0.0):
        k = (len(v) - 3) // 2
        return Profil(v[0], v[1:1 + k], v[1 + k:1 + 2 * k], v[-2], v[-1],
                      plat, beta0)

    # -- evaluation --------------------------------------------------------
    def courbure(self, beta):
        """Rayon de courbure **total**, Dirac du plat exclu.

        Le Dirac n'y figure pas : un plat est une droite, donc usinable. La
        contrainte de rayon d'outil porte sur ce que rend cette fonction.
        """
        r = self.a0
        for i, n in enumerate(range(2, self.M + 1)):
            r += self.an[i] * math.cos(n * beta) + self.bn[i] * math.sin(n * beta)
        return r + self.rho_plat(beta)

    def serie(self, beta):
        """La serie seule, sans la part lissee du plat."""
        r = self.a0
        for i, n in enumerate(range(2, self.M + 1)):
            r += self.an[i] * math.cos(n * beta) + self.bn[i] * math.sin(n * beta)
        return r

    def plancher(self, mini, grille=720):
        """Minimum impose a la **serie**, echantillon par echantillon.

        La contrainte d'usinabilite porte sur la courbure totale : la serie
        doit donc rester au-dessus de `mini - rho_plat(beta)`.
        """
        return [mini - self.rho_plat(TAU * i / grille) for i in range(grille)]

    def courbure_mini(self, n=720):
        return min(self.courbure(TAU * i / n) for i in range(n))

    # -- le plat du zero, porte explicitement ---------------------------
    def _x(self, beta):
        return (beta - self.beta0) % TAU

    def h_plat(self, beta):
        return (self.plat * math.sin(0.5 * self._x(beta))
                if self.plat > 0.0 else 0.0)

    def s_plat(self, beta):
        return (0.5 * self.plat * math.cos(0.5 * self._x(beta))
                if self.plat > 0.0 else 0.0)

    def rho_plat(self, beta):
        """Part **lissee** que le plat ajoute a la courbure (hors Dirac)."""
        return (0.75 * self.plat * math.sin(0.5 * self._x(beta))
                if self.plat > 0.0 else 0.0)

    def H(self, beta):
        h = self.a0 + self.cx * math.cos(beta) + self.cy * math.sin(beta)
        for i, n in enumerate(range(2, self.M + 1)):
            h += (self.an[i] * math.cos(n * beta)
                  + self.bn[i] * math.sin(n * beta)) * self._w[i]
        return h + self.h_plat(beta)

    def s(self, beta):
        d = -self.cx * math.sin(beta) + self.cy * math.cos(beta)
        for i, n in enumerate(range(2, self.M + 1)):
            d += n * (-self.an[i] * math.sin(n * beta)
                      + self.bn[i] * math.cos(n * beta)) * self._w[i]
        return d + self.s_plat(beta)

    def point(self, beta):
        H, s = self.H(beta), self.s(beta)
        u, t = u_of(beta), t_of(beta)
        return (H * u[0] + s * t[0], H * u[1] + s * t[1])

    def contour(self, n=360):
        return [self.point(TAU * i / n) for i in range(n)]

    def ferme(self):
        """Defaut de fermeture : harmoniques de rang 1 de rho (nulles ici)."""
        return 0.0, 0.0


def _refit(vals, M):
    """Coefficients de rangs 0 et 2..M d'un signal echantillonne (rang 1 jete,
    ce qui garantit la fermeture du contour)."""
    n = len(vals)
    a0 = sum(vals) / n
    an, bn = [], []
    for k in range(2, M + 1):
        c = s_ = 0.0
        for i, v in enumerate(vals):
            x = TAU * i / n
            c += v * math.cos(k * x)
            s_ += v * math.sin(k * x)
        an.append(2.0 * c / n)
        bn.append(2.0 * s_ / n)
    return a0, an, bn


def reparer(vec, mini=0.0, M=None, iters=25, grille=720,
            plat=0.0, beta0=0.0):
    """Ramene un jeu de coefficients dans l'espace admissible.

    Projections alternees entre :
      - l'ensemble convexe { courbure totale >= mini pour tout beta } ;
      - le sous-espace a bande limitee sans harmonique de rang 1
        (fermeture du contour).

    `mini` est le **rayon de courbure usinable** : le profil est convexe donc
    sans angle rentrant, et rho = 0 serait un angle vif. C'est donc exactement
    la contrainte de raccord usinable.

    Avec un plat porte explicitement, la contrainte porte sur la courbure
    **totale** : la serie doit rester au-dessus de `mini - rho_plat(beta)`,
    qui depend de beta. Un plancher scalaire imposait le minimum a la seule
    serie, ce qui gonflait sa courbure la ou le plat en apportait deja, et
    deformait le profil (0.24..1.41 mm au lieu de 0.30..3.00).
    """
    cx, cy = vec[-2], vec[-1]
    k = (len(vec) - 3) // 2
    M = M or (k + 1)
    modele = Profil(vec[0], vec[1:1 + k], vec[1 + k:1 + 2 * k], cx, cy,
                    plat, beta0)
    sol = modele.plancher(mini, grille)
    if _NOYAU.DISPONIBLE:
        np = _NOYAU._np
        a0, an, bn = _NOYAU.reparer_coeffs(
            float(vec[0]),
            np.array(vec[1:1 + k], dtype=np.float64),
            np.array(vec[1 + k:1 + 2 * k], dtype=np.float64),
            float(mini), int(iters), int(grille),
            np.array(sol, dtype=np.float64))
        return [float(a0)] + [float(v) for v in an] + [float(v) for v in bn] \
            + [cx, cy]
    pr = modele
    for _ in range(iters):
        vals = [pr.serie(TAU * i / grille) for i in range(grille)]
        if min(v - sol[i] for i, v in enumerate(vals)) >= -1e-9:
            break
        vals = [max(sol[i], v) for i, v in enumerate(vals)]
        a0, an, bn = _refit(vals, M)
        pr = Profil(a0, an, bn, cx, cy, plat, beta0)
    # garantie finale : on releve le terme constant du deficit residuel
    m = min(pr.serie(TAU * i / grille) - sol[i] for i in range(grille))
    if m < 0.0:
        pr = Profil(pr.a0 - m + 1e-6, pr.an, pr.bn, cx, cy, plat, beta0)
    return pr.vector()


# ------------------------------------------------------- depuis l'analytique

def from_support(betas, Hs, ss, M=24, mini=0.0, n=720, r_max=None,
                 plat=0.0, beta0=0.0):
    """Profil a bande limitee ajuste sur une fonction d'appui echantillonnee.

    On echantillonne la courbure exacte (y compris ses Dirac, integres par
    cellule), on l'ajuste en serie, puis on la repare pour qu'elle respecte
    le rayon de courbure usinable. Un coeur analytique a un **plat** au zero,
    donc une courbure infinie : il n'appartient pas a l'espace a bande
    limitee, et la graine en est l'approximation la plus proche qui soit
    usinable.

    `plat` et `beta0` sont seulement transmis au profil, pour memoire et pour
    la mesure : le plat reel rendu par la serie se lit avec
    `longueur_plat()`. Il ne fait que 0.31 mm d'un plat de 1.72 a douze
    harmoniques, et 0.67 encore a quatre-vingt-seize - voir `Profil` pour
    pourquoi un plat exact demande de reparametrer l'espace.
    """

    b0 = betas[0]
    d = TAU / n

    def interp(x):
        x = b0 + (x - b0) % TAU
        lo, hi = 0, len(betas) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if betas[m] <= x:
                lo = m
            else:
                hi = m
        db = betas[hi] - betas[lo]
        w = 0.0 if db < 1e-15 else (x - betas[lo]) / db
        return (Hs[lo] + w * (Hs[hi] - Hs[lo]),
                ss[lo] + w * (ss[hi] - ss[lo]))

    # --- le plat du zero, sorti de la serie
    #
    # Le plat est un **saut de s** dans la fonction d'appui : on le detecte,
    # on le retire des donnees, et on le porte ensuite explicitement. Il faut
    # le retirer des donnees elles-memes, et non de leur interpolation : la
    # discontinuite y est repartie sur un intervalle d'echantillonnage, et une
    # soustraction analytique laisserait une pointe de L/dbeta, soit plus de
    # deux cents, dans une cellule de courbure.
    if plat <= 0.0:
        j0, saut = -1, 0.0
        for i in range(len(betas) - 1):
            dd = ss[i + 1] - ss[i]
            if abs(dd) > abs(saut):
                j0, saut = i, dd
        if saut > 1e-6:
            plat, beta0 = saut, betas[j0]
    Hs0 = list(Hs)                 # cible d'origine, pour caler le centre
    Hs = list(Hs)
    ss = list(ss)
    if plat > 0.0:
        for i in range(len(betas)):
            x = (betas[i] - beta0) % TAU
            if x < 1e-12 and i <= (j0 if 'j0' in dir() else 0):
                x = TAU                     # limite a gauche du saut
            Hs[i] -= plat * math.sin(0.5 * x)
            ss[i] -= 0.5 * plat * math.cos(0.5 * x)

    vals = []
    for i in range(n):
        b1, b2 = d * i, d * (i + 1)
        H1, s1 = interp(b1)
        H2, s2 = interp(b2)
        vals.append((s2 - s1) / d + 0.5 * (H1 + H2))
    a0, an, bn = _refit(vals, M)
    pr = Profil(a0, an, bn, 0.0, 0.0, plat, beta0)
    # Centre : on cale cx, cy sur la fonction d'appui **d'origine**, plat
    # compris. Le faire contre la cible privee de son plat laissait un rang 1
    # de 0.6 mm dans l'ecart, c'est-a-dire un centre de rotation deplace
    # d'autant, et le profil sortait a 1.42..2.41 mm au lieu de 0.60..3.00.
    def cible(x):
        x = b0 + (x - b0) % TAU
        lo, hi = 0, len(betas) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if betas[m] <= x:
                lo = m
            else:
                hi = m
        db = betas[hi] - betas[lo]
        w = 0.0 if db < 1e-15 else (x - betas[lo]) / db
        return Hs0[lo] + w * (Hs0[hi] - Hs0[lo])

    ex = ey = 0.0
    for i in range(256):
        b = TAU * i / 256
        err = cible(b) - pr.H(b)
        ex += 2.0 * err * math.cos(b) / 256
        ey += 2.0 * err * math.sin(b) / 256
    pr.cx, pr.cy = ex, ey
    pr = Profil.from_vector(
        reparer(pr.vector(), mini, M, plat=plat, beta0=beta0), plat, beta0)
    if r_max:
        # l'ajustement a bande limitee deborde un peu l'encombrement :
        # on remet la graine a l'echelle (homothetie = facteur sur rho)
        re = max(norm(pr.point(TAU * i / 360)) for i in range(360))
        if re > r_max:
            k = r_max / re
            pr = Profil(pr.a0 * k, [v * k for v in pr.an],
                        [v * k for v in pr.bn], pr.cx * k, pr.cy * k,
                        pr.plat * k, beta0)
            pr = Profil.from_vector(
                reparer(pr.vector(), mini, M, plat=pr.plat, beta0=beta0),
                pr.plat, beta0)
    return pr


def from_chrono(p, res, M=24, mini=None, plat=None):
    """Profil a bande limitee issu d'une synthese analytique de chronographe.

    Le plat du zero est repris **explicitement** : c'est lui qui donne la
    position nulle, et une serie a bande limitee ne sait pas le porter.
    """
    from .chrono import support
    bet, Hs, ss = support(p, res)
    if mini is None:
        mini = getattr(p, "rayon_raccord", 0.0)
    if plat is None:
        plat = float(getattr(res, "plat", 0.0) or 0.0)
    return from_support(bet, Hs, ss, M, mini, r_max=p.r_max,
                        plat=plat, beta0=p.gamma0)


# ------------------------------------------------------------- evaluation

class Eval:
    """Resultat d'une evaluation complete d'un profil sur un tour."""

    def __init__(self):
        self.ok = False
        self.theta = []
        self.beta = []        # normale de contact, repere du coeur
        self.Hc = []          # appui au contact
        self.sc = []          # bras de levier de la normale
        self.phi = []         # angle du marteau (signe)
        self.Th = []          # couple signe sur le coeur [N.mm]
        self.N = []           # effort normal [N]
        self.psi = []         # ouverture du marteau [rad]
        self.bloque = []
        self.r_min = 0.0
        self.r_max = 0.0
        self.rho_min = 0.0
        self.plage_morte = 0.0
        self.ondulation = 0.0   # (max-min)/moyenne sur la partie utile
        self.couple_moyen = 0.0
        self.couple_utile = 0.0  # couple mini garanti hors des equilibres
        self.couple_p10 = 0.0    # 10e centile de |couple| sur le tour
        self.part_utile = 0.0    # fraction du tour prise en compte
        self.N_max = 0.0
        self.psi_max = 0.0
        self.transitions = 0     # changements de signe du couple sur un tour
        self.penalite = 0.0
        self.h0 = 0.0
        self.interference = 0
        self.contacts_coin = 0   # positions ou c'est un coin qui touche
        self.coin = []           # drapeau par position : un coin touche
        self.pen_max = 0.0
        self.rho_max = 0.0
        self.rayon_face = 0.0
        self.saut_contact = 0     # discontinuites de la position du marteau
        self.cote_change = 0      # le marteau passe de l'autre cote du coeur
        self.runs_couple = 0      # nombre de plages de signe du couple


def _remplir(e, res, p, pr, n, n_int=0, mu=None):
    """Transfere le resultat du noyau accelere dans un objet Eval."""
    th, be, H, sv, ph, T, N, bl, va = res[:9]
    for i in range(n):
        if not va[i]:
            e.penalite += 10.0
            continue
        e.theta.append(float(th[i]))
        e.beta.append(float(be[i]))
        e.Hc.append(float(H[i]))
        e.sc.append(float(sv[i]))
        e.phi.append(float(ph[i]))
        e.Th.append(float(T[i]))
        e.N.append(float(N[i]))
        e.psi.append(abs(float(ph[i])))
        e.bloque.append(bool(bl[i]))
    e.interference = n_int
    e.penalite += 1.0 * n_int
    _finir(e, p, pr, n, mu)
    return e


def mu_galet(rayon_galet, rayon_axe, mu_palier, f_roulement=0.005):
    """Frottement equivalent au contact d'un **galet** monte fou.

    Un plat glissant transmet mu directement. Un galet, lui, roule : ce qui
    reste au contact est le frottement de son propre axe, ramene par le
    rapport des rayons, plus la resistance au roulement :

        mu_eq = mu_palier . r_axe / R_galet  +  f_roulement / R_galet

    C'est le seul moyen de descendre la plage morte au sommet, qui vaut
    2.atan(mu) et ne depend que du frottement au contact.
    """
    R = max(abs(rayon_galet), 1e-9)
    return mu_palier * abs(rayon_axe) / R + f_roulement / R


def evaluer(p, pr, n=360, couple_cible=None, rayon_face=0.0, accel=True,
            mu_contact=None):
    """Analyse directe : pour chaque position du coeur, trouve le contact.

    `rayon_face` : rayon du plat du marteau.
        0   -> face plane (cas de reference) ;
        > 0 -> face **convexe**, bombee vers le coeur ;
        < 0 -> face **concave**, qui enveloppe le coeur.

    Pour une face courbe on travaille sur la **developpee** du coeur : le
    centre du cercle de la face est en  Y(beta) = (H + R).u + s.t  dans le
    repere du coeur. Comme ce centre est a distance fixe L du pivot, chaque
    position du coeur donne l'equation scalaire |Y(beta) - P_theta| = L, dont
    les racines sont les contacts possibles ; on retient celle qui exige la
    plus grande ouverture du marteau.

    L'interet : avec une face plane la normale de contact est rigidement liee
    a l'angle du marteau, donc dtheta/dbeta vaut environ 1 et la fenetre de
    frottement 2.atan(mu) se transporte telle quelle en angle de coeur. Une
    face courbe decouple les deux et permet de comprimer localement cette
    fenetre.
    """
    e = Eval()
    P = p.pivot
    mu = p.mu if mu_contact is None else mu_contact
    g0 = p.gamma0
    eps = sgn(dot(t_of(g0), P))
    if abs(dot(t_of(g0), P)) < 1e-9:
        return e
    psi_max = abs(math.radians(p.angle_plan)) * 0.999

    # echantillonnage du profil
    nb = 720
    bet = [TAU * i / nb for i in range(nb + 1)]
    # le noyau compile ne demande que les tableaux H et s : il sert aussi
    # pour un profil importe, dont seul l'echantillonnage reste en Python
    rapide = accel and _NOYAU.DISPONIBLE
    if rapide and hasattr(pr, "a0"):
        np = _NOYAU._np
        aH, aS, aR = _NOYAU.echantillonner(
            pr.a0, np.asarray(pr.an, dtype=np.float64),
            np.asarray(pr.bn, dtype=np.float64), pr.cx, pr.cy, nb,
            float(getattr(pr, "plat", 0.0)),
            float(getattr(pr, "beta0", 0.0)))
        Hs, ss = aH, aS
        e.rho_min = float(aR[:nb].min())
        e.rho_max = float(aR[:nb].max())
    elif rapide:
        np = _NOYAU._np
        Hs = np.array([pr.H(b) for b in bet], dtype=np.float64)
        ss = np.array([pr.s(b) for b in bet], dtype=np.float64)
        rc = [pr.courbure(b) for b in bet[:nb]]
        e.rho_min, e.rho_max = min(rc), max(rc)
    else:
        Hs = [pr.H(b) for b in bet]
        ss = [pr.s(b) for b in bet]
        e.rho_min = min(pr.courbure(b) for b in bet[:nb])
        e.rho_max = max(pr.courbure(b) for b in bet[:nb])
    rays = [math.hypot(Hs[i], ss[i]) for i in range(nb)]
    e.r_min, e.r_max = min(rays), max(rays)
    if rapide:
        e.r_min = float(min(rays))
        e.r_max = float(max(rays))
    if min(Hs) <= 0.0 or e.rho_min < -1e-3:
        e.penalite = 1e3 + abs(min(min(Hs), e.rho_min))
        return e

    # cale du marteau : tangent au profil pour l'appui le plus proche
    h0 = float(min(Hs)) - dot(u_of(g0), P)
    e.h0 = h0

    def sup(beta):
        x = beta % TAU
        i = int(x / TAU * nb)
        i = min(i, nb - 1)
        w = (x - bet[i]) / (bet[i + 1] - bet[i])
        return (Hs[i] + w * (Hs[i + 1] - Hs[i]),
                ss[i] + w * (ss[i + 1] - ss[i]))

    C0, k = p.couple_min, p.pente
    e.rayon_face = rayon_face
    if p.longueur_plan > 1e-9 and abs(rayon_face) <= 1e-9:
        if rapide:
            np = _NOYAU._np
            aH = Hs if hasattr(Hs, "dtype") else np.asarray(Hs, np.float64)
            aS = ss if hasattr(ss, "dtype") else np.asarray(ss, np.float64)
            eps = sgn(dot(t_of(g0), P))
            res = _NOYAU.scan_plan_fini(
                aH, aS, bet[0], bet[1] - bet[0], P[0], P[1], g0, h0, eps, mu,
                p.couple_min, p.pente,
                abs(math.radians(p.angle_plan)) * 0.999, n,
                0.5 * p.longueur_plan, p.rayon_coin, max(1, nb // 240),
                float(getattr(p, "decalage_plan", 0.0)))
            r = _remplir(e, res[:9], p, pr, n, mu=mu)
            r.contacts_coin = int(res[9])
            r.coin = [bool(x) for x, v in zip(res[10], res[8]) if v]
            return r
        return _evaluer_plan_fini(p, pr, e, bet, Hs, ss, nb, n, mu,
                                  0.5 * p.longueur_plan, p.rayon_coin,
                                  float(getattr(p, "decalage_plan", 0.0)))
    if rapide and (abs(rayon_face) <= 1e-9 or hasattr(pr, "an")):
        r = _accelere(p, pr, e, bet, Hs, ss, nb, n, rayon_face, h0, mu)
        if r is not None:
            return r
    if abs(rayon_face) > 1e-9:
        return _evaluer_courbe(p, pr, e, bet, Hs, ss, nb, n, rayon_face, mu)

    for i in range(n):
        th = TAU * i / n
        f = lambda psi: dot(u_of(g0 + eps * psi), P) + h0 - sup(g0 + eps * psi - th)[0]
        if f(0.0) >= 0.0:
            psi = 0.0
        elif f(psi_max) < 0.0:
            e.penalite += 10.0
            continue
        else:
            lo, hi = 0.0, psi_max
            for _ in range(50):
                m = 0.5 * (lo + hi)
                if f(m) < 0.0:
                    lo = m
                else:
                    hi = m
            psi = 0.5 * (lo + hi)
        phi = eps * psi
        gamma = g0 + phi
        H, s = sup(gamma - th)
        q = s - dot(t_of(gamma), P)
        Cm = -eps * (C0 + k * psi)

        def torque(sig):
            den = q + sig * mu * h0
            if abs(den) < 1e-12:
                return float("nan"), float("nan")
            return Cm * (s + sig * mu * H) / den, abs(Cm / den)

        TA, NA = torque(+1.0)
        TB, NB = torque(-1.0)
        goA = TA == TA and TA < 0.0           # le coeur peut reculer
        goB = TB == TB and TB > 0.0           # ou avancer
        if goA and not goB:
            T, N = TA, NA
        elif goB and not goA:
            T, N = TB, NB
        elif goA and goB:
            T, N = (TA, NA) if abs(TA) >= abs(TB) else (TB, NB)
        else:
            T, N = 0.0, (NA if NA == NA else NB)
        e.theta.append(math.degrees(th))
        e.beta.append(gamma - th)
        e.Hc.append(H)
        e.sc.append(s)
        e.phi.append(phi)
        e.Th.append(T)
        e.N.append(N if N == N else 0.0)
        e.psi.append(psi)
        e.bloque.append(not (goA or goB))

    if len(e.Th) < n // 2:
        e.penalite += 100.0
        return e
    e.plage_morte = _plage_morte(e, n, mu, pr)

    libres = [abs(t) for t, b in zip(e.Th, e.bloque) if not b and abs(t) > 1e-9]
    if libres:
        e.couple_moyen = sum(libres) / len(libres)
    _utile(e, n, mu)
    e.N_max = max(e.N) if e.N else 0.0
    e.psi_max = max(e.psi) if e.psi else 0.0

    # equilibres : on compte les changements de signe en ignorant les zones
    # bloquees (couple nul). Un coeur correct en a exactement deux : le zero
    # (stable) et le sommet (instable).
    _controles(e)
    e.ok = True
    return e


def pose(p, pr, e, theta_deg):
    """Elements a tracer pour une position du coeur (profil discretise)."""
    if not e.ok or not e.theta:
        return None
    th = theta_deg % 360.0
    i = min(range(len(e.theta)),
            key=lambda k: min(abs(e.theta[k] - th),
                              abs(e.theta[k] - th - 360.0),
                              abs(e.theta[k] - th + 360.0)))
    a = math.radians(e.theta[i])
    ca, sa = math.cos(a), math.sin(a)
    rot = lambda v: (v[0] * ca - v[1] * sa, v[0] * sa + v[1] * ca)
    ub, tb = u_of(e.beta[i]), t_of(e.beta[i])
    contact = rot((e.Hc[i] * ub[0] + e.sc[i] * tb[0],
                   e.Hc[i] * ub[1] + e.sc[i] * tb[1]))
    # Le marteau est un solide : sa face est definie **une fois au repos**
    # puis tournee rigidement de phi autour du pivot. La tracer a partir du
    # contact ou du pied de la perpendiculaire issue de O la ferait glisser,
    # et le marteau semblerait se translater.
    R = e.rayon_face
    phi = e.phi[i]
    L = p.demi_plan()
    if abs(R) > 1e-9:
        j0 = min(range(len(e.phi)), key=lambda k: abs(e.phi[k]))
        u0 = u_of(e.beta[j0])
        t0 = t_of(e.beta[j0])
        c0 = (e.Hc[j0] * u0[0] + e.sc[j0] * t0[0],
              e.Hc[j0] * u0[1] + e.sc[j0] * t0[1])
        Q0 = (c0[0] + R * u0[0], c0[1] + R * u0[1])
        b0 = math.atan2(c0[1] - Q0[1], c0[0] - Q0[0])
        sp = min(0.9, L / max(abs(R), 1e-9))
        face = [(Q0[0] + abs(R) * math.cos(b0 + sp * (2.0 * j / 24 - 1.0)),
                 Q0[1] + abs(R) * math.sin(b0 + sp * (2.0 * j / 24 - 1.0)))
                for j in range(25)]
    else:
        u0, t0 = u_of(p.gamma0), t_of(p.gamma0)
        H0 = dot(u0, p.pivot) + e.h0
        # contour materiel complet : les deux arrondis de bout en font partie,
        # sinon le marteau parait fini par une arete vive alors que c'est
        # l'arrondi qui porte le contact hors du plan
        face = contour_plan(u0, t0, H0, L, getattr(p, "rayon_coin", 0.0),
                            dec=getattr(p, "decalage_plan", 0.0))
    face = [_tourner(q, p.pivot, phi) for q in face]
    drap = getattr(e, "marche", None)
    return {"theta": a, "phi": e.phi[i], "profile": [rot(v) for v in pr.contour(360)],
            "face": face, "contact": contact,
            "normal": u_of(e.beta[i] + a), "H": e.Hc[i],
            "N": e.N[i], "Th": e.Th[i], "bloque": e.bloque[i],
            "franchit": bool(drap[i]) if drap and i < len(drap) else False}


def _degage(pr, R, b, H, sv, prof_pts, tol):
    """La face posee en beta ne traverse pas le profil."""
    u, t = u_of(b), t_of(b)
    Y = ((H + R) * u[0] + sv * t[0], (H + R) * u[1] + sv * t[1])
    if R > 0.0:
        d = min((Y[0] - X[0]) ** 2 + (Y[1] - X[1]) ** 2 for X in prof_pts)
        return d >= (R - tol) ** 2
    d = max((Y[0] - X[0]) ** 2 + (Y[1] - X[1]) ** 2 for X in prof_pts)
    return d <= (abs(R) + tol) ** 2


def _raffiner(pr, R, Pt, L, b, db, iters=4):
    """Secante sur |Y(beta) - Pt|^2 = L^2 avec les valeurs exactes."""
    def g(x):
        u, t = u_of(x), t_of(x)
        H, sv = pr.H(x), pr.s(x)
        Y = ((H + R) * u[0] + sv * t[0], (H + R) * u[1] + sv * t[1])
        return (Y[0] - Pt[0]) ** 2 + (Y[1] - Pt[1]) ** 2 - L * L

    x0, x1 = b - 0.25 * db, b + 0.25 * db
    f0, f1 = g(x0), g(x1)
    for _ in range(iters):
        if abs(f1 - f0) < 1e-18:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) < 1e-12:
            x1 = x2
            break
        x0, f0, x1, f1 = x1, f1, x2, g(x2)
    return x1 if abs(x1 - b) < 2.0 * db else b


def _evaluer_courbe(p, pr, e, bet, Hs, ss, nb, n, R, mu=None):
    """Analyse directe pour une face de marteau circulaire (R signe)."""
    P = p.pivot
    if mu is None:
        mu = p.mu
    g0 = p.gamma0
    C0, k = p.couple_min, p.pente

    # points du profil (repere du coeur), pour le controle d'interference
    step = max(1, nb // 180)
    prof_pts = []
    for i in range(0, nb, step):
        u, t = u_of(bet[i]), t_of(bet[i])
        prof_pts.append((Hs[i] * u[0] + ss[i] * t[0],
                         Hs[i] * u[1] + ss[i] * t[1]))

    # developpee : centre du cercle de la face, dans le repere du coeur
    Y = []
    for i in range(nb + 1):
        b = bet[i]
        u, t = u_of(b), t_of(b)
        Y.append(((Hs[i] + R) * u[0] + ss[i] * t[0],
                  (Hs[i] + R) * u[1] + ss[i] * t[1]))

    # calage : au repos la face touche le profil la ou l'appui est minimal
    i0 = min(range(nb), key=lambda i: Hs[i])
    # on amene ce point de contact sur la normale nominale g0
    th0 = bet[i0] - g0
    Q0 = rot(Y[i0], th0)
    L = norm((Q0[0] - P[0], Q0[1] - P[1]))
    if L < 1e-6:
        e.penalite += 100.0
        return e
    ang_q0 = math.atan2(Q0[1] - P[1], Q0[0] - P[0])
    psi_max = abs(math.radians(p.angle_plan)) * 0.999
    # sens d'ouverture du marteau : fixe une fois pour toutes. Le ressort
    # pousse toujours dans le meme sens, il ne peut pas s'inverser en cours
    # de course.
    eps_c = sgn(dot(t_of(p.gamma0), P))
    # tolerance d'interference : au-dela des micro-concavites laissees par la
    # troncature de la serie (quelques microns), a l'echelle de l'usinage
    tol = max(0.02, 1e-2 * p.r_min)

    # le balayage de racines se fait sur une grille allegee : la racine est
    # ensuite raffinee sur les valeurs exactes, la grille ne sert qu'a la
    # localiser et a separer les branches
    pas = max(1, nb // 180)
    idx = list(range(0, nb + 1, pas))
    if idx[-1] != nb:
        idx.append(nb)

    for i in range(n):
        th = TAU * i / n + th0
        ct, st = math.cos(-th), math.sin(-th)
        Pt = (P[0] * ct - P[1] * st, P[0] * st + P[1] * ct)   # pivot vu du coeur
        f = {j: (Y[j][0] - Pt[0]) ** 2 + (Y[j][1] - Pt[1]) ** 2 - L * L
             for j in idx}
        best = None
        for a, b_ in zip(idx, idx[1:]):
            j, jn = a, b_
            if f[j] == 0.0 or f[j] * f[jn] < 0.0:
                w = 0.0 if f[jn] == f[j] else f[j] / (f[j] - f[jn])
                b = bet[j] + w * (bet[jn] - bet[j])
                # raffinement sur les valeurs exactes : l'interpolation
                # lineaire laisse une erreur de contact de quelques microns,
                # qui ferait croire a une interference
                b = _raffiner(pr, R, Pt, L, b, bet[jn] - bet[j])
                u, t = u_of(b), t_of(b)
                H, sv = pr.H(b), pr.s(b)
                Yb = ((H + R) * u[0] + sv * t[0], (H + R) * u[1] + sv * t[1])
                Qf = rot(Yb, th)
                phi = math.atan2(Qf[1] - P[1], Qf[0] - P[0]) - ang_q0
                phi = (phi + math.pi) % TAU - math.pi
                if abs(phi) > psi_max:
                    continue
                # Le ressort pousse le marteau vers sa position fermee : il
                # s'arrete au **premier** contact rencontre, donc a la plus
                # petite ouverture sans interpenetration. Retenir la plus
                # grande ouverture ferait sauter le contact de l'autre cote
                # du coeur.
                if not _degage(pr, R, b, H, sv, prof_pts, tol):
                    continue
                if best is None or abs(phi) < abs(best[0]):
                    best = (phi, b, H, sv, th)
        if best is None:
            e.penalite += 10.0
            e.interference += 1
            continue
        phi, b, H, sv, _ = best
        gamma = b + th
        u, t = u_of(gamma), t_of(gamma)
        X = (H * u[0] + sv * t[0], H * u[1] + sv * t[1])
        q = sv - dot(t, P)
        ee = H - dot(u, P)
        Cm = -eps_c * (C0 + k * abs(phi))

        def torque(sig):
            den = q + sig * mu * ee
            if abs(den) < 1e-12:
                return float("nan"), float("nan")
            return Cm * (sv + sig * mu * H) / den, abs(Cm / den)

        TA, NA = torque(+1.0)
        TB, NB = torque(-1.0)
        goA = TA == TA and TA < 0.0
        goB = TB == TB and TB > 0.0
        if goA and not goB:
            T, N = TA, NA
        elif goB and not goA:
            T, N = TB, NB
        elif goA and goB:
            T, N = (TA, NA) if abs(TA) >= abs(TB) else (TB, NB)
        else:
            T, N = 0.0, (NA if NA == NA else NB)
        e.theta.append(math.degrees(TAU * i / n))
        e.beta.append(b)
        e.Hc.append(H)
        e.sc.append(sv)
        e.phi.append(phi)
        e.Th.append(T)
        e.N.append(N if N == N else 0.0)
        e.psi.append(abs(phi))
        e.bloque.append(not (goA or goB))

    _finir(e, p, pr, n, mu)
    return e


def _utile(e, n, mu=0.12, marge_deg=None):
    """Restreint les criteres de couple a la partie **utile** du parcours.

    Le couple s'annule forcement aux **deux** equilibres : le zero (stable)
    et le sommet (instable). Les inclure fausse les criteres : la moyenne
    masque un couple qui s'effondre, et l'ondulation rapportee a la moyenne
    vaut mecaniquement 150 a 200 % quel que soit le profil.

    La fenetre ecartee est calee sur ce qui est **physiquement inevitable** :
    le cone de frottement, de demi-largeur atan(mu), plus 3 deg de marge.
    Tout ce qui depasse est compte comme un defaut. On ecarte donc une
    fenetre autour de ces deux positions,
    reperees par les changements de signe du couple **en ignorant les zeros
    pour la detection**. Tout autre echantillon a couple nul reste compte :
    une version precedente excluait tous les zeros, ce qui revenait a cacher
    les plages mortes du critere et poussait l'optimiseur a les allonger.
    """
    if marge_deg is None:
        marge_deg = math.degrees(math.atan(max(mu, 1e-6))) + 3.0
    m = len(e.Th)
    if m < 8:
        return
    sg = [1 if t > 1e-12 else (-1 if t < -1e-12 else 0) for t in e.Th]
    nz = [i for i in range(m) if sg[i] != 0]
    if not nz:
        e.couple_utile = 0.0
        e.ondulation = 10.0
        e.part_utile = 0.0
        return
    bascules = []
    for a, b in zip(nz, nz[1:] + nz[:1]):
        if sg[a] != sg[b]:
            bascules.append((a + ((b - a) % m) // 2) % m)
    larg = max(1, int(marge_deg / 360.0 * m))
    exclu = set()
    for i in bascules:
        for d in range(-larg, larg + 1):
            exclu.add((i + d) % m)
    vals = [abs(e.Th[i]) for i in range(m) if i not in exclu]
    e.part_utile = len(vals) / m
    if len(vals) < m // 4:
        e.couple_utile = 0.0
        e.ondulation = 10.0
        return
    moy = sum(vals) / len(vals)
    e.couple_utile = min(vals)
    e.ondulation = (max(vals) - min(vals)) / max(moy, 1e-12)
    # centile robuste : insensible aux deux rampes inevitables, c'est lui
    # qui sert a caler le ressort
    tri = sorted(abs(t) for t in e.Th)
    e.couple_p10 = tri[max(0, int(0.10 * len(tri)))] if tri else 0.0


def _controles(e, seuil_saut=math.radians(25.0)):
    """Verifie la coherence mecanique du tour.

    - le marteau ne doit pas changer de cote du coeur : phi garde un signe ;
    - sa position ne doit pas sauter d'un pas a l'autre (sinon le contact a
      bascule sur une autre partie du profil) ;
    - le couple ne doit changer de sens que deux fois sur un tour : au zero
      (equilibre stable) et au sommet (instable). Un changement de sens en
      pleine descente de flanc voudrait dire que le coeur repart en arriere,
      ce qui est impossible s'il doit ramener son plat sur celui du marteau.
    """
    m = len(e.phi)
    if m < 3:
        return
    # seuil relatif a la course : au repos phi passe par zero, le bruit
    # numerique ne doit pas compter comme un changement de cote
    amp = max(abs(v) for v in e.phi)
    seuil = max(1e-6, 0.02 * amp)
    pos = sum(1 for v in e.phi if v > seuil)
    neg = sum(1 for v in e.phi if v < -seuil)
    e.cote_change = min(pos, neg)
    e.saut_contact = sum(1 for i in range(m)
                         if abs(e.phi[i] - e.phi[i - 1]) > seuil_saut)
    sg = [1 if t > 1e-12 else (-1 if t < -1e-12 else 0) for t in e.Th]
    nz = [x for x in sg if x != 0]
    if nz:
        e.runs_couple = sum(1 for i in range(len(nz))
                            if nz[i] != nz[(i + 1) % len(nz)])
    e.transitions = e.runs_couple


def _marge_large(m, step):
    """Plus large plage continue ou `m` est negatif, bords interpoles."""
    k = len(m)
    if k == 0:
        return 0.0
    if all(x < 0.0 for x in m):
        return k * step
    if all(x >= 0.0 for x in m):
        return 0.0
    pire = 0.0
    for i in range(k):
        if m[i] >= 0.0 or m[(i - 1) % k] < 0.0:
            continue                       # i n'ouvre pas une plage
        j = i
        while m[(j + 1) % k] < 0.0 and j - i <= k:
            j += 1
        a, b = m[(i - 1) % k], m[i]        # bord gauche : a >= 0 > b
        fl = a / (a - b) if a != b else 0.0
        c, d = m[j % k], m[(j + 1) % k]    # bord droit : c < 0 <= d
        fr = c / (c - d) if c != d else 0.0
        pire = max(pire, (j - i + fl + fr) * step)
    return pire


def _plage_morte(e, n, mu, pr=None, grille=2880):
    """Largeur de la plage morte, **independante du pas d'echantillonnage**.

    Compter les echantillons bloques donnait un multiple de 360/n : un seul
    echantillon isole pesait 3 deg a n=120 et 1 deg a n=360, ce qui faisait
    varier le critere d'optimisation avec la finesse du calcul. On mesure donc
    la marge d'arc-boutement

        marge = |s| - mu*H          (negatif = cale)

    sur la grille propre du profil, qui ne depend pas de n, et on interpole ses
    passages par zero. L'angle de normale beta et l'angle du coeur theta
    varient a l'unite pres (beta = gamma - theta, gamma ne bougeant que de
    quelques degres sur un tour), la largeur en beta vaut donc celle en theta.
    """
    if (mu is not None and pr is not None and hasattr(pr, "H")
            and len(e.beta) >= 8 and len(e.theta) == len(e.beta)):
        # beta(theta) est lisse : on le deroule, on l'interpole sur une grille
        # fine en theta, puis on evalue la marge avec le H et le s **exacts**
        # du profil. Seules les normales reellement visitees comptent ainsi.
        k = len(e.beta)
        bu = [e.beta[0]]
        for i in range(1, k):
            d = e.beta[i] - e.beta[i - 1]
            bu.append(bu[-1] + (d + math.pi) % TAU - math.pi)
        d = e.beta[0] - e.beta[-1]
        bfin = bu[-1] + (d + math.pi) % TAU - math.pi
        pas = 360.0 / grille
        m = []
        for i in range(grille):
            x = i * k / grille
            j = int(x)
            w = x - j
            b0 = bu[j]
            b1 = bu[j + 1] if j + 1 < k else bfin
            b = b0 + w * (b1 - b0)
            m.append(abs(pr.s(b)) - mu * pr.H(b))
        return min(_marge_large(m, pas), 360.0)
    step = 360.0 / n
    if mu is None or not e.Hc:
        runs, cur = [], 0
        for b in e.bloque + e.bloque:
            if b:
                cur += 1
            else:
                runs.append(cur)
                cur = 0
        runs.append(cur)
        return min(max(runs), n) * step

    m = [abs(sv) - mu * H for H, sv in zip(e.Hc, e.sc)]
    return min(_marge_large(m, step), 360.0)


def _finir(e, p, pr, n, mu=None):
    if len(e.Th) < n // 2:
        e.penalite += 100.0
        return
    e.plage_morte = _plage_morte(e, n, mu, pr)
    libres = [abs(t) for t, b in zip(e.Th, e.bloque) if not b and abs(t) > 1e-9]
    if libres:
        e.couple_moyen = sum(libres) / len(libres)
    _utile(e, n, mu)
    e.N_max = max(e.N) if e.N else 0.0
    e.psi_max = max(e.psi) if e.psi else 0.0
    _controles(e)
    e.ok = True


# ------------------------------------------------------- chemin accelere

def _accelere(p, pr, e, bet, Hs, ss, nb, n, R, h0, mu=None):
    """Appelle le noyau compile. Renvoie None si le cas n'est pas gere."""
    np = _NOYAU._np
    P = p.pivot
    if mu is None:
        mu = p.mu
    g0 = p.gamma0
    eps = sgn(dot(t_of(g0), P))
    psi_max = abs(math.radians(p.angle_plan)) * 0.999
    aH = Hs if hasattr(Hs, "dtype") else np.asarray(Hs, dtype=np.float64)
    aS = ss if hasattr(ss, "dtype") else np.asarray(ss, dtype=np.float64)
    dbet = bet[1] - bet[0]
    if abs(R) <= 1e-9:
        res = _NOYAU.scan_plan(aH, aS, bet[0], dbet, P[0], P[1], g0, h0,
                               eps, mu, p.couple_min, p.pente, psi_max, n)
        return _remplir(e, res, p, pr, n, mu=mu)

    # calage de la face circulaire, identique au chemin Python
    i0 = int(min(range(nb), key=lambda i: Hs[i]))
    b0 = bet[i0]
    u0, t0 = u_of(b0), t_of(b0)
    Y0 = ((Hs[i0] + R) * u0[0] + ss[i0] * t0[0],
          (Hs[i0] + R) * u0[1] + ss[i0] * t0[1])
    th0 = b0 - g0
    Q0 = rot(Y0, th0)
    L = norm((Q0[0] - P[0], Q0[1] - P[1]))
    if L < 1e-6:
        return None
    ang_q0 = math.atan2(Q0[1] - P[1], Q0[0] - P[0])
    pas = max(1, nb // 180)
    step = max(1, nb // 180)
    px_, py_ = [], []
    for i in range(0, nb, step):
        u, t = u_of(bet[i]), t_of(bet[i])
        px_.append(Hs[i] * u[0] + ss[i] * t[0])
        py_.append(Hs[i] * u[1] + ss[i] * t[1])
    res = _NOYAU.scan_courbe(
        aH, aS, bet[0], dbet,
        pr.a0, np.asarray(pr.an, dtype=np.float64),
        np.asarray(pr.bn, dtype=np.float64), pr.cx, pr.cy,
        P[0], P[1], R, L, ang_q0, th0, eps, mu, p.couple_min, p.pente,
        psi_max, n, pas, np.asarray(px_, dtype=np.float64),
        np.asarray(py_, dtype=np.float64), max(0.02, 1e-2 * p.r_min))
    return _remplir(e, res[:9], p, pr, n, int(res[9]), mu=mu)


# --------------------------------------------------- dimensionnement du ressort

class Ressort:
    """Ressort de rappel du marteau deduit du couple vise sur le coeur."""

    def __init__(self):
        self.ok = False
        self.message = ""
        self.prearmage = 0.0     # C0, couple a marteau ferme [N.mm]
        self.raideur = 0.0       # k [N.mm/rad]
        self.couple_maxi = 0.0   # couple a pleine ouverture [N.mm]
        self.course = 0.0        # ouverture maxi du marteau [rad]
        self.ecart = 0.0         # ecart type relatif de l'ajustement [%]
        self.eval = None         # Eval obtenue avec ce ressort


def dimensionner_ressort(p, pr, rayon_face=0.0, mu_contact=None, n=360,
                         couple_cible=None):
    """Calcule le ressort du marteau a partir du couple vise sur le coeur.

    C'est le sens naturel du probleme : on connait l'encombrement, le pivot
    du marteau, la pente de son plan et le couple de rappel voulu ; le
    ressort est ce qu'il reste a determiner.

    Le rapport de transmission `T / Cm` ne depend que de la geometrie et du
    frottement, pas du niveau de couple. On le mesure donc une fois avec un
    couple marteau unitaire, puis on en deduit le couple marteau **exige**
    en chaque position pour obtenir le couple vise :

        Cm_exige(theta) = couple_vise / rapport(theta)

    Un ressort reel est affine en l'ouverture du marteau. On ajuste donc
    `Cm = C0 + k.psi` au sens des moindres carres sur la partie utile du
    parcours : **C0 est le prearmage, k la raideur**, et le residu de
    l'ajustement est l'ondulation de couple qu'on ne peut pas eviter avec un
    ressort lineaire.
    """
    r = Ressort()
    cible = couple_cible if couple_cible is not None else p.couple_coeur
    # 1. rapport de transmission, mesure a couple marteau unitaire
    q = p.copy(couple_min=1.0, pente=0.0)
    e1 = evaluer(q, pr, n=n, rayon_face=rayon_face, mu_contact=mu_contact)
    if not e1.ok or not e1.Th:
        r.message = "cinematique non exploitable"
        return r
    m = len(e1.Th)
    sg = [1 if t > 1e-12 else (-1 if t < -1e-12 else 0) for t in e1.Th]
    bascules = [i for i in range(m)
                if sg[i] != 0 and sg[i - 1] != 0 and sg[i] != sg[i - 1]]
    bascules += [i for i in range(m) if sg[i] == 0]
    larg = max(1, int(20.0 / 360.0 * m))
    exclu = set()
    for i in bascules:
        for d in range(-larg, larg + 1):
            exclu.add((i + d) % m)
    idx = [i for i in range(m) if i not in exclu and abs(e1.Th[i]) > 1e-12]
    if len(idx) < m // 8:
        r.message = "partie utile trop courte"
        return r

    # 2. couple marteau exige en chaque position
    psi = [e1.psi[i] for i in idx]
    cm = [cible / abs(e1.Th[i]) for i in idx]

    # 3. ajustement affine Cm = C0 + k.psi, robuste.
    #    Pres des equilibres le rapport de transmission tend vers zero, donc
    #    le couple marteau exige tend vers l'infini : quelques points de
    #    rampe suffisent a tirer un ajustement par moindres carres. On ecarte
    #    donc les valeurs aberrantes avant d'ajuster.
    med = sorted(cm)[len(cm) // 2]
    garde = [i for i in range(len(cm)) if cm[i] <= 2.5 * med]
    if len(garde) >= max(8, len(cm) // 3):
        psi = [psi[i] for i in garde]
        cm = [cm[i] for i in garde]
    npt = len(psi)
    sx = sum(psi); sy = sum(cm)
    sxx = sum(v * v for v in psi); sxy = sum(a * b for a, b in zip(psi, cm))
    det = npt * sxx - sx * sx
    if abs(det) < 1e-18:
        k = 0.0
        C0 = sy / npt
    else:
        k = (npt * sxy - sx * sy) / det
        C0 = (sy - k * sx) / npt
    if C0 <= 0.0:                       # un prearmage negatif n'a pas de sens
        k = max(0.0, sxy / sxx) if sxx > 1e-18 else 0.0
        C0 = max(1e-6, (sy - k * sx) / npt)
    if k < 0.0:                         # ressort qui se detend : impossible
        k = 0.0
        C0 = sy / npt
    r.prearmage, r.raideur = C0, k
    r.course = max(e1.psi)
    r.couple_maxi = C0 + k * r.course
    res = [abs((C0 + k * a) - b) / max(b, 1e-12) for a, b in zip(psi, cm)]
    r.ecart = 100.0 * sum(res) / len(res)

    # 4. verification avec ce ressort. On ne remet **pas** le resultat a
    #    l'echelle : l'ajustement vise deja la consigne en chaque point du
    #    parcours utile, et recaler sur une statistique du couple figerait
    #    mecaniquement le critere qui la mesure. L'ecart qui subsiste est
    #    precisement ce qu'un ressort lineaire ne peut pas suivre, et c'est
    #    ce que l'optimisation doit reduire en changeant le profil.
    e2 = evaluer(p.copy(couple_min=C0, pente=k), pr, n=n,
                 rayon_face=rayon_face, mu_contact=mu_contact)
    r.prearmage, r.raideur = C0, k
    r.couple_maxi = C0 + k * r.course
    r.eval = e2
    r.ok = bool(e2 and e2.ok)
    return r


def fenetre_pointe(e, demi=10.0):
    """Part de la fenetre autour de la pointe ou le **coin** ne touche pas.

    La pointe du coeur est l'endroit ou le contact est au plus grand rayon.
    C'est la qu'un contact sur l'arete vive du plan arc-boute : la plage morte
    y vaut alors 2.atan(mu), soit pres de 14 deg. Quand c'est l'**arrondi** du
    marteau qui porte, le bras de levier ne s'annule plus et la plage morte
    tombe a environ un degre - c'est ce que montre le coeur de
    `Sketches - minute.dxf`, a 1.2 deg.

    On veut donc que l'arrondi soit en contact sur toute la fenetre de +/-
    `demi` degres autour du passage de la pointe. Cette fonction renvoie le
    pourcentage de cette fenetre ou ce n'est **pas** le cas : zero est le but.
    """
    if not e.coin or not e.theta or len(e.coin) != len(e.theta):
        return 100.0
    rayons = [math.hypot(H * math.cos(b) - sv * math.sin(b),
                         H * math.sin(b) + sv * math.cos(b))
              for H, sv, b in zip(e.Hc, e.sc, e.beta)]
    ip = max(range(len(rayons)), key=lambda i: rayons[i])
    th_p = e.theta[ip]
    dans = [i for i in range(len(e.theta))
            if min(abs(e.theta[i] - th_p), 360.0 - abs(e.theta[i] - th_p)) <= demi]
    if not dans:
        return 100.0
    return 100.0 * sum(1 for i in dans if not e.coin[i]) / len(dans)


def secteur_marche(pr, e, tol=None):
    """Positions ou le marteau **ne touche pas** le profil reel.

    Un plan, et meme un coin arrondi, ne peut toucher que l'enveloppe convexe
    du coeur. La marche qui ramene la pointe au plat du zero est une
    concavite : sur ce secteur le marteau franchit la corde de l'enveloppe et
    ne touche rien. Ce n'est pas un defaut de calcul, c'est une chute libre -
    le modele quasi statique n'y a pas de sens.

    Renvoie `(largeur_en_degres, drapeaux)`, un drapeau par position.
    """
    pts = getattr(pr, "pts", None)
    if not pts or not e.theta or not e.Hc:
        return 0.0, []
    if tol is None:
        # trois fois le pas moyen du contour : au-dela, le point de contact
        # n'est plus sur le dessin mais sur une corde de l'enveloppe
        n = len(pts)
        pas = sum(norm((pts[(i + 1) % n][0] - pts[i][0],
                        pts[(i + 1) % n][1] - pts[i][1]))
                  for i in range(n)) / max(n, 1)
        tol = max(3.0 * pas, 1e-4)
    drap = []
    for H, sv, b in zip(e.Hc, e.sc, e.beta):
        ub, tb = u_of(b), t_of(b)
        X = (H * ub[0] + sv * tb[0], H * ub[1] + sv * tb[1])
        d = min(norm((X[0] - q[0], X[1] - q[1])) for q in pts)
        drap.append(d > tol)
    m = len(drap)
    return 360.0 * sum(1 for x in drap if x) / max(m, 1), drap


def balayage_plan(p, pr, n_pts=13, demi=10.0, n=180, rayon_face=0.0,
                  mu_contact=None, lo=None, hi=None):
    """Compromis sur la longueur du plan du marteau, longueur par longueur.

    Deux exigences s'opposent et aucune ne se tranche seule :

    * **court**, le plan met son arrondi en contact au passage de la pointe,
      ou le bras de levier ne s'annule plus : la plage morte s'effondre ;
    * **trop court**, il ne peut plus enjamber la marche qui ramene la pointe
      au plat du zero - le coeur y est concave et le marteau tombe sans rien
      toucher, sur un secteur qui peut depasser 60 deg.

    La fenetre de pointe n'est d'ailleurs **pas monotone** en longueur de plan
    (80 % a 3.0 mm, 40 % a 1.6, 45 % a 1.4 sur le coeur synthetise) : la
    chercher par dichotomie, comme je le faisais, n'a pas de sens. On balaye
    donc, et on rend le tableau.

    Chaque entree : `{"L", "hors_pointe", "plage_morte", "marche", "N_max"}`.
    """
    q = p.copy()
    if lo is None:
        lo = 0.15 * max(p.r_max, 1e-3)
    if hi is None:
        hi = 1.2 * max(p.r_max, 1e-3)
    out = []
    for i in range(n_pts):
        L = lo + (hi - lo) * i / max(n_pts - 1, 1)
        q.longueur_plan = L
        try:
            e = evaluer(q, pr, n=n, rayon_face=rayon_face,
                        mu_contact=mu_contact)
        except Exception:
            continue
        if e is None or not e.ok:
            continue
        marche, _ = secteur_marche(pr, e)
        out.append({"L": L, "hors_pointe": fenetre_pointe(e, demi),
                    "plage_morte": e.plage_morte, "marche": marche,
                    "N_max": e.N_max})
    return out


def plan_recommande(p, pr, marche_maxi=5.0, **kw):
    """Longueur de plan conseillee, et front de compromis a defaut.

    Quand une longueur garde le secteur sans contact sous `marche_maxi`
    degres, on prend parmi celles-la la plus petite plage morte, et a egalite
    le plus faible effort de contact.

    Sinon il n'y a **pas** de bon choix a annoncer : la marche du coeur est
    trop profonde pour qu'un plan l'enjambe. Choisir au plus petit secteur
    sans contact donnerait le plan le plus long, donc la plage morte la pire
    (11.75 deg sur le coeur simple par defaut) ; choisir a la plus petite
    plage morte donnerait l'inverse. On rend donc le front non domine sur
    (plage morte, secteur sans contact) et on laisse l'arbitrage au
    concepteur.

    Renvoie `(entree_ou_None, tableau, front)`.
    """
    tab = balayage_plan(p, pr, **kw)
    if not tab:
        return None, tab, []
    bons = [d for d in tab if d["marche"] <= marche_maxi]
    if bons:
        best = min(bons, key=lambda d: (round(d["plage_morte"], 2),
                                        d["N_max"]))
        return best, tab, [best]
    front = []
    for d in tab:
        domine = any(o is not d
                     and o["plage_morte"] <= d["plage_morte"]
                     and o["marche"] <= d["marche"]
                     and (o["plage_morte"] < d["plage_morte"]
                          or o["marche"] < d["marche"])
                     for o in tab)
        if not domine:
            front.append(d)
    return None, tab, front


def longueur_plat(pr, beta0, tol=0.002, demi_fenetre=1.2, grille=1200):
    """Longueur du plat du coeur autour de la normale `beta0`.

    Le plat du zero est ce qui donne la **position nulle** : c'est lui qui
    vient se poser sur le plan du marteau, et sans lui le coeur n'a pas de
    position de repos definie. Une representation a bande limitee ne peut pas
    porter un vrai plat - il demanderait un Dirac de courbure - donc rien
    n'empeche l'optimisation de le dissoudre si on ne le mesure pas.

    On le mesure sur la geometrie : la tangente a `beta0` est la droite
    `u(beta0).x = H(beta0)` ; le plat est l'etendue, le long de cette droite,
    des points du contour qui s'en ecartent de moins de `tol`.
    """
    u0, t0 = u_of(beta0), t_of(beta0)
    H0 = pr.H(beta0)
    smin = smax = None
    for i in range(grille + 1):
        b = beta0 - demi_fenetre + 2.0 * demi_fenetre * i / grille
        H, sv = pr.H(b), pr.s(b)
        ub, tb = u_of(b), t_of(b)
        X = (H * ub[0] + sv * tb[0], H * ub[1] + sv * tb[1])
        if abs(u0[0] * X[0] + u0[1] * X[1] - H0) > tol:
            continue
        a = t0[0] * X[0] + t0[1] * X[1]
        smin = a if smin is None or a < smin else smin
        smax = a if smax is None or a > smax else smax
    if smin is None:
        return 0.0
    return smax - smin


def balayage_marteau(p, pr, n_t=9, n_a=7, demi_t=0.3, demi_a=20.0,
                     n=180, critere="ondulation"):
    """Balaye les deux grandeurs du marteau qui portent, coeur fixe.

    Ce sont la position **angulaire** du pivot - son deplacement
    perpendiculaire a O->P - et la **pente du pan**. Mesure a l'appui : sur
    +/- 4 mm, deplacer le pivot radialement fait varier l'ondulation de 1.6
    point, perpendiculairement de 68.8 points. La distance a l'axe n'apporte
    donc quasiment rien.

    Un balayage direct de ces deux variables bat ici la recherche
    multi-objectif : NSGA-II plafonne a 23.9 % d'ondulation avec 32 individus
    sur 30 generations, la ou cette grille en trouve 12.5. C'est normal - il
    repartit son effort sur quatre objectifs et sur tout l'espace des formes
    de coeur, alors qu'on cherche ici deux nombres.

    Le ressort est redimensionne a chaque point : toutes les variantes
    delivrent le couple vise, on compare ce que ce resultat coute.

    Renvoie `(grille, meilleur)`. Chaque entree porte `t`, `angle`, `pivot`,
    `ondulation`, `N_max`, `plage_morte`, `sous_consigne` et `residu`.
    """
    P0 = (p.pivot_x, p.pivot_y)
    d0 = norm(P0) or 1.0
    ut = (-P0[1] / d0, P0[0] / d0)
    a0 = p.angle_plan
    plat, b_unused = getattr(pr, "plat", 0.0), None
    grille, meilleur = [], None
    for i in range(n_t):
        t = (-demi_t + 2.0 * demi_t * i / max(n_t - 1, 1)) * d0
        for j in range(n_a):
            ang = a0 - demi_a + 2.0 * demi_a * j / max(n_a - 1, 1)
            if abs(ang) > 88.0:
                continue
            q = p.copy(pivot_x=P0[0] + t * ut[0], pivot_y=P0[1] + t * ut[1],
                       angle_plan=ang)
            # le plat du coeur suit la normale de repos
            pq = (Profil(pr.a0, pr.an, pr.bn, pr.cx, pr.cy, plat, q.gamma0)
                  if hasattr(pr, "a0") else pr)
            try:
                res = dimensionner_ressort(q, pq, 0.0, None, n=n)
            except Exception:
                continue
            e = res.eval
            if e is None or not e.ok or not res.ok:
                continue
            d = {"t": t, "angle": ang, "pivot": (q.pivot_x, q.pivot_y),
                 "ondulation": 100.0 * e.ondulation, "N_max": e.N_max,
                 "plage_morte": e.plage_morte, "residu": res.ecart,
                 "sous_consigne": secteur_faible(e, p.couple_coeur),
                 "C0": res.prearmage, "k": res.raideur}
            grille.append(d)
            cle = {"ondulation": "ondulation", "effort": "N_max",
                   "plage": "plage_morte",
                   "consigne": "sous_consigne"}.get(critere, "ondulation")
            if meilleur is None or d[cle] < meilleur[cle]:
                meilleur = d
    return grille, meilleur


def secteur_faible(e, cible):
    """Part du tour (en %) ou le couple est sous la consigne.

    C'est le critere le plus direct, et le seul sans fenetre arbitraire : une
    longue plage a couple nul le fait monter, point. Les rampes inevitables
    de part et d'autre des deux equilibres y contribuent d'un minimum
    incompressible, identique pour toutes les solutions comparees.
    """
    if not e.Th:
        return 100.0
    return 100.0 * sum(1 for t in e.Th if abs(t) < cible) / len(e.Th)


# -------------------------------------------------- profil importe d'un DXF

def support_de_points(pts, n=720):
    """Fonction d'appui d'un contour donne par ses points.

    L'origine est le **centre de rotation** du coeur. Le marteau a une face
    plane : il ne voit que l'enveloppe convexe du contour, et c'est donc
    elle que decrit la fonction d'appui

        H(beta) = max sur les points de  u(beta).X
        s(beta) = t(beta).X  au point qui realise ce maximum

    Les sauts de s aux sommets sont les segments droits du contour (le plat
    du zero en est un) ; l'integration par cellule de `from_support` les
    traite exactement.
    """
    bet, Hs, ss = [], [], []
    for i in range(n + 1):
        b = TAU * i / n
        cb, sb = math.cos(b), math.sin(b)
        best, bh = None, -1e18
        for q in pts:
            h = cb * q[0] + sb * q[1]
            if h > bh:
                bh, best = h, q
        bet.append(b)
        Hs.append(bh)
        ss.append(-sb * best[0] + cb * best[1])
    return bet, Hs, ss


def depuis_points(pts, M=24, mini=0.0, r_max=None, n=720):
    """Construit un Profil a partir d'un contour importe (origine = axe)."""
    bet, Hs, ss = support_de_points(pts, n)
    return from_support(bet, Hs, ss, M, mini, n, r_max)


class ProfilImporte:
    """Contour lu dans un DXF, simule **tel quel**.

    Meme interface qu'un `Profil`, mais la fonction d'appui n'est pas
    ajustee : elle est calculee directement sur les points,

        H(beta) = max sur les points de  u(beta).X

    c'est-a-dire l'enveloppe convexe du dessin, la seule chose qu'une face
    plane puisse toucher. On simule donc la geometrie dessinee, sans passer
    par la reconstruction a bande limitee qui, elle, ne sert que de graine a
    l'optimisation.

    La courbure n'est qu'indicative : sur un contour discretise elle vaut 0
    le long des segments et diverge aux sommets. Elle est donc lissee, et ne
    doit pas servir de controle d'usinabilite.
    """

    M = 0

    def __init__(self, pts, n=1440):
        self.pts = [(float(q[0]), float(q[1])) for q in pts]
        self.n = n
        self.d = TAU / n
        self._bet, self._H, self._s = support_de_points(self.pts, n)
        self.cx = self.cy = 0.0
        self._rho = self._courbures()

    def _courbures(self, lis=9):
        """rho = H + H'' estime par differences finies, puis lisse."""
        n, d = self.n, self.d
        raw = []
        for i in range(n):
            h0 = self._H[(i - 1) % n]
            h1 = self._H[i]
            h2 = self._H[(i + 1) % n]
            raw.append(h1 + (h2 - 2.0 * h1 + h0) / (d * d))
        out = []
        for i in range(n):
            out.append(sum(raw[(i + k) % n] for k in range(-lis, lis + 1))
                       / (2 * lis + 1))
        return out

    def _i(self, beta):
        return int((beta % TAU) / self.d) % self.n

    def H(self, beta):
        i = self._i(beta)
        w = (beta % TAU) / self.d - i
        return self._H[i] + w * (self._H[i + 1] - self._H[i])

    def s(self, beta):
        i = self._i(beta)
        w = (beta % TAU) / self.d - i
        return self._s[i] + w * (self._s[i + 1] - self._s[i])

    def courbure(self, beta):
        return self._rho[self._i(beta)]

    def courbure_mini(self, n=720):
        return min(self._rho)

    def point(self, beta):
        H, sv = self.H(beta), self.s(beta)
        u, t = u_of(beta), t_of(beta)
        return (H * u[0] + sv * t[0], H * u[1] + sv * t[1])

    def contour(self, n=360):
        """Contour **vu par le contact** : l'enveloppe convexe echantillonnee.

        C'est elle qu'un plan, ou un coin arrondi, peut toucher. Tracer les
        points eux-memes montrerait 0.3 mm d'ecart apparent dans les
        concavites, alors que le calcul y travaille correctement sur la corde
        de l'enveloppe ; la reconstruction, elle, ne s'en ecarte que de 14 um
        au plus (au sommet d'un coeur synthetise).

        `contour_reel()` donne les points tels qu'ils ont ete lus ou
        synthetises, et `secteur_marche()` dit ou les deux divergent.
        """
        return [self.point(TAU * i / n) for i in range(n)]

    def contour_reel(self, n=360):
        """Les points tels qu'ils ont ete lus dans le DXF ou synthetises."""
        return list(self.pts)

    def ferme(self):
        return 0.0, 0.0


# ------------------------------------------- plan de marteau de longueur finie

def _ouverture_nez(Hs, nb, Pt, LA, acA, LB, acB, th, rho, psi_max, eps):
    """Ouverture mini du marteau telle que son nez ne penetre pas le coeur.

    Le nez est le segment des deux centres de coin, dilate du rayon `rho`. Le
    coeur est convexe (son enveloppe, seule chose qu'un plan puisse toucher),
    donc intersection de demi-plans `u_j.x <= H_j`. Un segment est hors d'une
    telle intersection **si et seulement si** un appui j le separe en entier,
    c'est-a-dire met ses **deux** extremites au-dela :

        min(u_j.A, u_j.B) >= H_j + rho

    L'ensemble admissible est alors exactement l'union sur j de ces
    conditions, et la plus petite ouverture admissible est le minimum sur j.

    Traiter les trois contacts separement - plan, coin A, coin B - puis prendre
    la plus grande des trois ouvertures etait faux : l'ensemble admissible d'un
    coin est un **arc**, pas un intervalle croissant. Sur `test.dxf` le coin
    cote -t etait degage a psi = 0, rentrait dans le coeur a psi = 4.6, et la
    regle du maximum retenait 4.6 : le marteau penetrait le coeur de 0.34 mm
    et le couple sautait de 26 deg d'une position a l'autre.

    Renvoie `(psi, j, tendu)` ou `tendu` vaut 0 si c'est l'extremite A qui
    touche, 1 pour B, 2 si les deux (appui parallele au plan), ou
    `(None, -1, -1)` si aucune ouverture ne degage le nez.
    """
    aA = acA - th
    aB = acB - th
    best, jb, tb = None, -1, -1
    for j in range(nb):
        bj = TAU * j / nb
        cb, sb = math.cos(bj), math.sin(bj)
        num = Hs[j] + rho - (cb * Pt[0] + sb * Pt[1])
        cA, cB = num / LA, num / LB
        if cA > 1.0 or cB > 1.0:
            continue                      # cet appui ne peut separer le nez
        bornes = []
        for c, om in ((cA, bj - aA), (cB, bj - aB)):
            if c < -1.0:
                bornes.append(None)       # toujours satisfait
            else:
                w = math.acos(c)
                lo, hi = om - w, om + w
                if eps < 0.0:
                    lo, hi = -hi, -lo
                bornes.append((lo, hi))
        cands = [0.0]
        for b in bornes:
            if b is not None:
                cands.append(b[0] % TAU)

        def dedans(psi):
            for b in bornes:
                if b is not None and (psi - b[0]) % TAU > b[1] - b[0] + 1e-12:
                    return False
            return True

        for psi in sorted(cands):
            if psi > psi_max:
                break
            if dedans(psi):
                if best is None or psi < best:
                    best = psi
                    jb = j
                    # laquelle des deux extremites est tendue ?
                    dA = abs((psi - bornes[0][0]) % TAU) if bornes[0] else 9.0
                    dB = abs((psi - bornes[1][0]) % TAU) if bornes[1] else 9.0
                    tb = 2 if abs(dA - dB) < 1e-9 else (0 if dA < dB else 1)
                break
    if best is None:
        return None, -1, -1
    return best, jb, tb


def _evaluer_plan_fini(p, pr, e, bet, Hs, ss, nb, n, mu, demi, rho,
                       dec=0.0):
    """Contact avec un plan de marteau **de longueur finie**.

    Un plan infini est une commodite de calcul : si le point de tangence sort
    du plan reel, ce n'est plus lui qui touche le coeur mais son **coin**
    (eventuellement arrondi au rayon `rayon_coin`). Le marteau se ferme
    jusqu'au premier obstacle, donc parmi les trois contacts possibles - le
    plan, et chacun des deux coins - c'est celui qui exige la **plus grande
    ouverture** qui s'impose.
    """
    P = p.pivot
    g0 = p.gamma0
    C0, k = p.couple_min, p.pente
    eps = sgn(dot(t_of(g0), P))
    psi_max = abs(math.radians(p.angle_plan)) * 0.999
    u0, t0 = u_of(g0), t_of(g0)
    H0 = float(min(Hs))
    h0 = H0 - dot(u0, P)
    M0 = (H0 * u0[0], H0 * u0[1])          # pied de la normale, au repos
    coins = []
    for sgc in (-1.0, 1.0):
        # centre rentre de rho dans la matiere : tangent au plan, sans
        # depasser devant lui (cf. noyau.scan_plan_fini)
        c0 = (M0[0] + (dec + sgc * demi) * t0[0] + rho * u0[0],
              M0[1] + (dec + sgc * demi) * t0[1] + rho * u0[1])
        d = norm((c0[0] - P[0], c0[1] - P[1]))
        coins.append((d, math.atan2(c0[1] - P[1], c0[0] - P[0]), sgc))

    def interp(beta):
        x = beta % TAU
        i = min(int(x / TAU * nb), nb - 1)
        w = (x - bet[i]) / (bet[i + 1] - bet[i])
        return (Hs[i] + w * (Hs[i + 1] - Hs[i]),
                ss[i] + w * (ss[i + 1] - ss[i]))

    LA, acA, _ = coins[0]
    LB, acB, _ = coins[1]
    for i in range(n):
        th = TAU * i / n
        cands = []
        psi, j, tendu = _ouverture_nez(Hs, nb, (P[0] * math.cos(-th)
                                                - P[1] * math.sin(-th),
                                                P[0] * math.sin(-th)
                                                + P[1] * math.cos(-th)),
                                       LA, acA, LB, acB, th, rho, psi_max, eps)
        if psi is not None:
            phi = eps * psi
            bb = TAU * j / nb
            ct, st = math.cos(-th), math.sin(-th)
            Pt = (P[0] * ct - P[1] * st, P[0] * st + P[1] * ct)
            Ah = (Pt[0] + LA * math.cos(acA - th + eps * psi),
                  Pt[1] + LA * math.sin(acA - th + eps * psi))
            Bh = (Pt[0] + LB * math.cos(acB - th + eps * psi),
                  Pt[1] + LB * math.sin(acB - th + eps * psi))
            ub, tb = u_of(bb), t_of(bb)
            # Point du coeur qui touche : son point d'appui dans la direction
            # de la droite separatrice. S'il se projette **a l'interieur** du
            # segment, c'est le plan qui porte et le bras de levier est celui
            # de ce point ; sinon c'est un coin, et le bras de levier est
            # celui du pied de la perpendiculaire issue de son centre.
            Hj, sj = interp(bb)
            Xj = (Hj * ub[0] + sj * tb[0], Hj * ub[1] + sj * tb[1])
            dxs, dys = Bh[0] - Ah[0], Bh[1] - Ah[1]
            l2 = dxs * dxs + dys * dys
            tproj = (0.5 if l2 < 1e-18 else
                     ((Xj[0] - Ah[0]) * dxs + (Xj[1] - Ah[1]) * dys) / l2)
            if 1e-9 < tproj < 1.0 - 1e-9:
                ee = Hj - dot(u_of(bb + th), P)
                cands.append((psi, phi, bb, Hj, sj, ee, "plan"))
            else:
                Ch = Ah if tproj <= 0.5 else Bh
                Hb = ub[0] * Ch[0] + ub[1] * Ch[1] - rho
                sb = tb[0] * Ch[0] + tb[1] * Ch[1]
                ee = Hb - dot(u_of(bb + th), P)
                cands.append((psi, phi, bb, Hb, sb, ee, "coin"))

        if not cands:
            e.penalite += 10.0
            continue
        # le marteau s'arrete au premier obstacle : la plus grande ouverture
        _a, phi, beta, H, sv, ee, quoi = max(cands, key=lambda c: c[0])
        gam = beta + th
        q = sv - dot(t_of(gam), P)
        Cm = -eps * (C0 + k * abs(phi))

        def torque(sig):
            den = q + sig * mu * ee
            if abs(den) < 1e-12:
                return float("nan"), float("nan")
            return Cm * (sv + sig * mu * H) / den, abs(Cm / den)

        TA, NA = torque(+1.0)
        TB, NB = torque(-1.0)
        goA = TA == TA and TA < 0.0
        goB = TB == TB and TB > 0.0
        if goA and not goB:
            T, N = TA, NA
        elif goB and not goA:
            T, N = TB, NB
        elif goA and goB:
            T, N = (TA, NA) if abs(TA) >= abs(TB) else (TB, NB)
        else:
            T, N = 0.0, (NA if NA == NA else NB)
        e.theta.append(math.degrees(th))
        e.beta.append(beta)
        e.Hc.append(H)
        e.sc.append(sv)
        e.phi.append(phi)
        e.Th.append(T)
        e.N.append(N if N == N else 0.0)
        e.psi.append(abs(phi))
        e.bloque.append(not (goA or goB))
        e.coin.append(quoi == "coin")
        e.contacts_coin += 1 if quoi == "coin" else 0
    e.h0 = h0
    _finir(e, p, pr, n, mu)
    e.ok = True
    return e
