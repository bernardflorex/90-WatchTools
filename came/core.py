"""
Noyau de calcul pour le dimensionnement d'un coeur horloger (came de remise a zero).

Modele
------
Le coeur tourne autour de O = (0,0). Le marteau pivote autour de P. Le marteau
porte un plan (ou un arc) de contact qui appuie sur le coeur ; un couple
Cm = C0 + k*|psi| (psi = ouverture du marteau) le ramene vers le coeur.

Le coeur est decrit par sa FONCTION D'APPUI H(beta) : distance de O a la
tangente de normale sortante u(beta), dans le repere du coeur. Le point de
contact a pour coordonnees (H, s) dans la base (u, t) avec s = dH/dbeta.

Bilan des moments au contact (normale u sortante du coeur, tangente t) :
    - point de contact      C = H*u + s*t
    - force sur le coeur    F = -N*u + S*t,  S = sigma*mu*N
    - couple sur le coeur   Th = N*s + S*H
    - bras de levier /P     q = s - t.P   (normale),  e = H - u.P  (tangente)
    - equilibre du marteau  Cm = N*q + S*e
  d'ou :
    N  = Cm / (q + sigma*mu*e)
    Th = Cm * (s + sigma*mu*H) / (q + sigma*mu*e)

sigma est impose par le sens de glissement relatif :
    v_rel_tangentiel = thetadot*H - phidot*e
    sigma = -sign(v_rel_tangentiel)
ce qui distingue le frottement RENTRANT (retour, thetadot<0) du frottement
SORTANT (armement, thetadot>0).

Synthese : on impose Th = -T0 (constant). L'equation ci-dessus donne s
explicitement, puis la coherence dH/dtheta = s*(dbeta/dtheta) fournit une EDO
du premier ordre sur phi(theta) que l'on integre en RK4. Le profil est ensuite
reconstruit par  X(beta) = H*u(beta) + s*t(beta).

Unites : mm, N, N.mm, rad en interne (degres dans l'interface).
"""

import math

TAU = 2.0 * math.pi


# ---------------------------------------------------------------- petite algebre

def u_of(a):
    return (math.cos(a), math.sin(a))


def t_of(a):
    return (-math.sin(a), math.cos(a))


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def rot(p, a):
    c, s = math.cos(a), math.sin(a)
    return (c * p[0] - s * p[1], s * p[0] + c * p[1])


def norm(p):
    return math.hypot(p[0], p[1])


def sgn(x):
    return 1.0 if x >= 0.0 else -1.0


def _tourner_autour(q, c, a):
    """Rotation rigide du point q de l'angle a autour du centre c."""
    ca, sa = math.cos(a), math.sin(a)
    dx, dy = q[0] - c[0], q[1] - c[1]
    return (c[0] + dx * ca - dy * sa, c[1] + dx * sa + dy * ca)


# ---------------------------------------------------------------- parametres

class Params:
    """Parametres d'entree du dimensionnement."""

    #: (nom, libelle, unite, mini, maxi, pas)
    SPECS = [
        ("pivot_x",    "Pivot marteau X",        "mm",      -60.0, 60.0, 0.5),
        ("pivot_y",    "Pivot marteau Y",        "mm",      -60.0, 60.0, 0.5),
        ("r_min",      "Distance plat - axe",    "mm",       0.05, 30.0, 0.05),
        ("r_max",      "Encombrement maxi (limite)", "mm",     0.5, 60.0, 0.1),
        ("couple_min", "Couple marteau mini C0", "N.mm",     0.01, 500.0, 0.1),
        ("pente",      "Pente marteau k",        "N.mm/rad", 0.0, 500.0, 0.1),
        ("mu",         "Frottement dynamique",   "-",         0.0, 0.6, 0.01),
        ("angle_plan", "Pente du pan du marteau", "deg",    -89.0, 89.0, 0.5),
        ("rayon_plan", "Rayon du plan (0=plat)", "mm",        0.0, 500.0, 0.5),
        ("sweep",      "Rotation synthetisee",   "deg",      10.0, 355.0, 1.0),
        ("n_lobes",    "Nombre de lobes",        "-",         1.0, 4.0, 1.0),
        ("longueur_plan", "Longueur du pan",     "mm",        0.0, 50.0, 0.1),
        ("rayon_coin", "Rayons de coin du pan",  "mm",        0.0, 5.0, 0.01),
        ("decalage_plan", "Decalage le long du pan", "mm",  -30.0, 30.0, 0.05),
        ("n_steps",    "Pas d'integration",      "-",        40.0, 2000.0, 20.0),
    ]

    def __init__(self, **kw):
        self.pivot_x = 16.0
        self.pivot_y = 5.0
        self.r_min = 1.0
        self.r_max = 5.0
        self.couple_min = 4.0
        self.pente = 2.0
        self.mu = 0.12
        self.angle_plan = 60.0
        self.rayon_plan = 0.0
        self.sweep = 350.0
        self.n_lobes = 1
        self.longueur_plan = 0.0     # 0 = plan suppose assez long
        self.rayon_coin = 0.1        # arrondi des extremites du plan
        # Position du milieu du plan le long de lui-meme, comptee depuis le
        # pied de la normale issue de O quand le marteau est au repos (plat du
        # coeur contre plat du marteau). Elle decide quel coin voit la pointe.
        self.decalage_plan = 0.0
        self.n_steps = 360
        for k, v in kw.items():
            setattr(self, k, v)

    # --- grandeurs derivees -------------------------------------------------
    @property
    def pivot(self):
        return (self.pivot_x, self.pivot_y)

    @property
    def alpha_p(self):
        return math.atan2(self.pivot_y, self.pivot_x)

    @property
    def gamma0(self):
        """Angle de la **normale** de contact au repos, repere fixe.

        `angle_plan` est mesure depuis la direction O->P (du centre du coeur
        vers le pivot du marteau) jusqu'a la normale du plan de contact, et
        **au repos**, coeur au zero.

        Il n'est volontairement pas rapporte a l'horizontale : rapporte a OP,
        il ne depend pas de l'endroit ou l'on place le pivot autour du coeur,
        et c'est lui qui fixe le bras de levier au pivot, q = |OP|.sin(pente).

        En cours de course le marteau tourne de phi : la normale devient
        gamma0 + phi et l'orientation absolue du plan change d'autant. C'est
        pourquoi la pente est definie au repos, et non par rapport a la
        verticale ou a l'horizontale.
        """
        return self.alpha_p + math.radians(self.angle_plan)

    def demi_plan(self):
        """Demi-longueur du plan de contact a tracer."""
        return (0.5 * self.longueur_plan if self.longueur_plan > 1e-9
                else 0.7 * self.r_max)

    def orientation_plan(self, phi=0.0):
        """Orientation absolue du plan de contact (deg / horizontale)."""
        return math.degrees(self.gamma0 + phi) + 90.0

    def copy(self, **kw):
        p = type(self)()
        p.__dict__.update(self.__dict__)
        p.__dict__.update(kw)
        return p

    def as_dict(self):
        return {k: getattr(self, k) for k, _l, _u, _a, _b, _s in self.SPECS}


# ---------------------------------------------------------------- face du marteau

def contour_plan(u0, t0, H0, demi, rc, n=10, dec=0.0):
    """Contour **materiel** du plan de contact, arrondis compris.

    `demi` est la demi-longueur de la partie droite ; chaque extremite porte
    un quart de cercle de rayon `rc` tangent au plan, de centre rentre de `rc`
    dans la matiere du marteau (+u0) pour ne pas depasser devant le plan.
    C'est cet arrondi qui porte le contact des que le coeur sort du plan.

    `dec` decale le milieu du plan le long de lui-meme, depuis le pied de la
    normale issue de O au repos : c'est lui qui decide lequel des deux coins
    verra passer la pointe du coeur.
    """
    rc = max(0.0, min(rc, 0.95 * demi)) if demi > 0 else 0.0

    def pt(a, b):
        return (a * u0[0] + b * t0[0], a * u0[1] + b * t0[1])

    if rc <= 1e-9:
        return [pt(H0, dec - demi), pt(H0, dec + demi)]
    pts = []
    for i in range(n, -1, -1):             # coin cote -t0 : flanc -> plan
        a = 0.5 * math.pi * i / n
        pts.append(pt(H0 + rc - rc * math.cos(a),
                      dec - demi - rc * math.sin(a)))
    for i in range(n + 1):                 # coin cote +t0 : plan -> flanc
        a = 0.5 * math.pi * i / n
        pts.append(pt(H0 + rc - rc * math.cos(a),
                      dec + demi + rc * math.sin(a)))
    return pts


class FlatFace:
    """Plan de contact rectiligne solidaire du marteau.

    La normale sortante (du coeur vers le marteau) tourne avec le marteau :
    gamma = gamma0 + phi.  La distance O->plan vaut  H = u(gamma).P + h0.
    """

    kind = "plan"

    def __init__(self, p):
        self.P = p.pivot
        self.gamma0 = p.gamma0
        # h0 cale la face pour que H(phi=0) = r_min
        self.h0 = p.r_min - dot(u_of(self.gamma0), self.P)
        self.rc = max(0.0, getattr(p, "rayon_coin", 0.0))
        self.dec = float(getattr(p, "decalage_plan", 0.0))

    def gamma(self, phi, s_target_fn=None, gamma_hint=None):
        return self.gamma0 + phi

    def H(self, phi, gamma):
        return dot(u_of(gamma), self.P) + self.h0

    def e(self, phi, gamma):
        return self.H(phi, gamma) - dot(u_of(gamma), self.P)   # = h0

    def dphi_dtheta(self, phi, gamma, s, side=1.0):
        """phi' tel que  dH/dtheta = s*(gamma' - side)  avec gamma' = phi'."""
        gp = dot(t_of(gamma), self.P)          # dH/dphi
        den = s - gp
        if abs(den) < 1e-12:
            return None
        return side * s / den

    def segment(self, phi, half=None):
        """Segment **materiel** a tracer, en rotation rigide autour du pivot.

        Les extremites sont definies une fois au repos puis tournees de phi
        autour de P. Les tracer a partir du pied de la perpendiculaire issue
        de O les ferait glisser le long de la face : le marteau semblerait
        se translater alors qu'il ne peut que tourner.
        """
        half = half or 1.0
        u0, t0 = u_of(self.gamma0), t_of(self.gamma0)
        H0 = dot(u0, self.P) + self.h0
        pts = contour_plan(u0, t0, H0, half, self.rc, dec=self.dec)
        return [_tourner_autour(q, self.P, phi) for q in pts]


class ArcFace:
    """Plan de contact en arc de cercle de rayon Rf (convexe vers le coeur).

    Le centre Q tourne avec le marteau. La normale gamma n'est plus liee a phi :
    elle est fixee par la condition geometrique  t(gamma).Q = s  (le point de
    contact appartient au cercle).
    """

    kind = "arc"

    def __init__(self, p):
        self.P = p.pivot
        self.Rf = p.rayon_plan
        self.gamma0 = p.gamma0
        # centre du cercle au repos, place sur la normale au repos
        c = dot(u_of(self.gamma0), self.P) - self.Rf - p.r_min
        self.q0 = (-c * math.cos(self.gamma0), -c * math.sin(self.gamma0))

    def Q(self, phi):
        r = rot(self.q0, phi)
        return (self.P[0] + r[0], self.P[1] + r[1])

    def H(self, phi, gamma):
        return dot(u_of(gamma), self.Q(phi)) - self.Rf

    def e(self, phi, gamma):
        r = rot(self.q0, phi)
        return dot(u_of(gamma), r) - self.Rf

    def gamma(self, phi, s_target_fn, gamma_hint):
        """Resout  t(gamma).Q(phi) - s_target(gamma) = 0  (secante)."""
        Q = self.Q(phi)

        def f(g):
            st = s_target_fn(phi, g)
            if st is None:
                return None
            return dot(t_of(g), Q) - st

        g0 = gamma_hint if gamma_hint is not None else self.gamma0
        g1 = g0 + 1e-3
        f0, f1 = f(g0), f(g1)
        if f0 is None or f1 is None:
            return None
        for _ in range(60):
            if abs(f1 - f0) < 1e-15:
                break
            g2 = g1 - f1 * (g1 - g0) / (f1 - f0)
            if abs(g2 - g1) < 1e-13:
                g1 = g2
                break
            g0, f0 = g1, f1
            g1 = g2
            f1 = f(g1)
            if f1 is None:
                return None
        return g1

    def dphi_dtheta(self, phi, gamma, s, side=1.0):
        """dH/dtheta = s.gamma' + A.phi' et dH/dtheta = s.(gamma' - side)
        donnent  phi' = -side.s / A,  avec A = u(gamma).(z x (Q-P))."""
        r = rot(self.q0, phi)
        A = dot(u_of(gamma), (-r[1], r[0]))
        if abs(A) < 1e-12:
            return None
        return -side * s / A

    def segment(self, phi, half=None):
        """Arc **materiel**, defini au repos puis tourne autour du pivot."""
        Q0 = self.Q(0.0)
        base = math.atan2(-Q0[1], -Q0[0])
        span = min(math.radians(45.0),
                   (half or 1.0) / max(abs(self.Rf), 1e-9))
        n = 24
        pts = []
        for i in range(n + 1):
            a = base - span + 2 * span * i / n
            pts.append((Q0[0] + self.Rf * math.cos(a),
                        Q0[1] + self.Rf * math.sin(a)))
        return [_tourner_autour(q, self.P, phi) for q in pts]


def make_face(p):
    return ArcFace(p) if p.rayon_plan > 1e-9 else FlatFace(p)


# ---------------------------------------------------------------- synthese

class Result:
    """Resultat d'une synthese."""

    def __init__(self):
        self.ok = False
        self.message = ""
        self.T0 = 0.0
        self.side = 1.0
        self.theta = []        # angle du coeur [rad]
        self.phi = []          # angle marteau [rad]
        self.psi = []          # ouverture marteau >=0 [rad]
        self.gamma = []        # normale de contact [rad]
        self.beta = []         # normale dans le repere coeur [rad]
        self.H = []            # fonction d'appui [mm]
        self.s = []            # dH/dbeta [mm]
        self.contact = []      # point de contact, repere coeur [mm]
        self.profile = []      # profil du coeur, repere coeur [mm]
        self.Cm = []           # couple marteau [N.mm]
        self.N = []            # effort normal [N]
        self.sigma = []        # sens du frottement
        self.Th_retour = []    # couple sur l'axe du coeur, retour
        self.Th_armement = []  # idem, sens armement
        self.Th_ideal = []     # couple sans frottement
        self.jam = []          # denominateur q + sigma*mu*e
        self.curv = []         # rayon de courbure H + H''
        # indicateurs
        self.r_eff = 0.0
        self.psi_max = 0.0
        self.delta_Cm = 0.0
        self.N_max = 0.0
        self.jam_min = 0.0
        self.curv_min = 0.0
        self.demi_requis = 0.0  # demi-longueur de plan exigee par le profil
        self.ripple = 0.0      # ondulation du couple (sens oppose)
        self.rendement = 0.0
        self.valid = False     # profil exploitable (geometrie + encombrement)
        self.geom_ok = False   # convexite / monotonie / arc-boutement
        self.warnings = []

    def closed_profile(self, n_lobes=1):
        """Profil ferme (lobes repetes + raccord de la chute)."""
        if not self.profile:
            return []
        base = list(self.profile)
        n = max(1, int(n_lobes))
        if n == 1:
            return base
        out = []
        for i in range(n):
            a = TAU * i / n
            out.extend(rot(q, a) for q in base)
        return out


def _integrate(p, face, T0, eps, mode="retour", side=1.0,
               theta0=0.0, phi0=0.0, theta_end=None, n=None, gamma0=None,
               sigma0=1.0):
    """Integre phi(theta) pour un T0 donne.

    `side` = +1 : le coeur tourne dans le sens direct quand theta croit ;
    `side` = -1 : flanc oppose (le coeur tourne dans l'autre sens), ce qui
    inverse le sens de glissement donc le frottement.
    """
    res = Result()
    res.T0 = T0
    res.side = side
    P = p.pivot
    mu = p.mu
    C0, k = p.couple_min, p.pente
    if theta_end is None:
        theta_end = theta0 + math.radians(p.sweep)
    sweep = theta_end - theta0
    n = max(10, int(n if n else p.n_steps))
    dth = sweep / n
    arc = (face.kind == "arc")

    def Cm_of(phi):
        # couple applique : referme le marteau (sens -eps)
        return -eps * (C0 + k * abs(phi))

    def s_target(phi, gamma, sigma):
        """s impose par Th = -T0."""
        Cm = Cm_of(phi)
        if abs(Cm) < 1e-15:
            return None
        lam = (-side * T0) / Cm
        if abs(1.0 - lam) < 1e-9:
            return None
        H = face.H(phi, gamma)
        e = face.e(phi, gamma)
        tP = dot(t_of(gamma), P)
        return (-lam * tP + sigma * mu * (lam * e - H)) / (1.0 - lam)

    state = {"gamma": face.gamma0 if gamma0 is None else gamma0,
             "sigma": sigma0}

    def solve_point(theta, phi):
        """Renvoie (gamma, H, s, sigma, phip) coherents, ou None."""
        best = None
        for sigma in (state["sigma"], -state["sigma"]):
            if arc:
                g = face.gamma(phi, lambda ph, gg: s_target(ph, gg, sigma), state["gamma"])
                if g is None:
                    continue
            else:
                g = face.gamma(phi)
            s = s_target(phi, g, sigma)
            if s is None:
                continue
            H = face.H(phi, g)
            e = face.e(phi, g)
            phip = face.dphi_dtheta(phi, g, s, side)
            if phip is None:
                continue
            # coherence du sens de glissement
            vrel = side * H - phip * e    # au signe de thetadot pres
            sig_pred = sgn(vrel) if mode == "retour" else -sgn(vrel)
            cand = (g, H, s, sigma, phip, e)
            if best is None:
                best = cand
            if abs(sig_pred - sigma) < 1e-9:
                return cand
        return best

    thetas, phis = [], []
    phi = phi0
    theta = theta0
    for i in range(n + 1):
        sp = solve_point(theta, phi)
        if sp is None:
            res.message = "integration impossible (theta = %.1f deg)" % math.degrees(theta)
            return res
        g, H, s, sigma, phip, e = sp
        state["gamma"], state["sigma"] = g, sigma

        thetas.append(theta)
        phis.append(phi)
        res.gamma.append(g)
        res.H.append(H)
        res.s.append(s)
        res.sigma.append(sigma)
        tP = dot(t_of(g), P)
        q = s - tP
        den = q + sigma * mu * e
        den_op = q - sigma * mu * e
        Cm = Cm_of(phi)
        res.Cm.append(abs(Cm))
        res.jam.append(den)
        if abs(den) < 1e-12:
            res.message = "arc-boutement (theta = %.1f deg)" % math.degrees(theta)
            return res
        N = Cm / den
        res.N.append(abs(N))
        res.Th_retour.append(Cm * (s + sigma * mu * H) / den)
        res.Th_armement.append(Cm * (s - sigma * mu * H) / den_op
                               if abs(den_op) > 1e-12 else float('nan'))
        res.Th_ideal.append(Cm * s / q if abs(q) > 1e-12 else float('nan'))

        if i == n:
            break
        # RK4 sur phi(theta)
        def f(th, ph):
            r = solve_point(th, ph)
            return None if r is None else r[4]
        k1 = phip
        k2 = f(theta + dth / 2, phi + dth * k1 / 2)
        if k2 is None:
            k2 = k1
        k3 = f(theta + dth / 2, phi + dth * k2 / 2)
        if k3 is None:
            k3 = k2
        k4 = f(theta + dth, phi + dth * k3)
        if k4 is None:
            k4 = k3
        phi += dth * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        theta += dth
        if abs(phi) > math.pi:
            res.message = "course du marteau excessive"
            return res

    res.theta = thetas
    res.phi = phis
    res.psi = [abs(x) for x in phis]
    res.beta = [res.gamma[i] - side * thetas[i] for i in range(len(thetas))]
    res.profile = [(res.H[i] * math.cos(res.beta[i]) - res.s[i] * math.sin(res.beta[i]),
                    res.H[i] * math.sin(res.beta[i]) + res.s[i] * math.cos(res.beta[i]))
                   for i in range(len(thetas))]
    res.contact = list(res.profile)
    res.ok = True
    res.demi_requis = _demi_requis(p, res)
    _indicators(res, p)
    return res


def _indicators(res, p):
    n = len(res.theta)
    res.r_eff = max(norm(q) for q in res.profile)
    res.psi_max = max(res.psi)
    res.delta_Cm = max(res.Cm) - min(res.Cm)
    res.N_max = max(res.N)
    res.jam_min = min(abs(x) for x in res.jam)
    # courbure : H + d2H/dbeta2, avec dH/dbeta = s
    curv = []
    for i in range(n):
        i0, i1 = max(0, i - 1), min(n - 1, i + 1)
        db = res.beta[i1] - res.beta[i0]
        ds = res.s[i1] - res.s[i0]
        curv.append(res.H[i] + (ds / db if abs(db) > 1e-12 else 0.0))
    res.curv = curv
    res.curv_min = min(curv)
    opp = [x for x in res.Th_armement if x == x]
    if opp:
        m = sum(abs(x) for x in opp) / len(opp)
        res.ripple = (max(abs(x) for x in opp) - min(abs(x) for x in opp)) / max(m, 1e-9)
    ideal = [abs(x) for x in res.Th_ideal if x == x]
    if ideal:
        res.rendement = abs(res.T0) / (sum(ideal) / len(ideal))
    # validite geometrique (independante de l'encombrement)
    w = []
    if min(res.H) <= 1e-6:
        w.append("fonction d'appui negative : profil non convexe")
    if any(res.H[i + 1] < res.H[i] - 1e-9 for i in range(n - 1)):
        w.append("rayon du coeur non monotone")
    if res.curv_min <= 0.0:
        w.append("rayon de courbure negatif (profil non convexe)")
    sj = sgn(res.jam[0])
    if any(sgn(x) != sj for x in res.jam):
        w.append("changement de signe du bras de levier : arc-boutement")
    res.geom_ok = not w
    we = []
    if res.r_eff > p.r_max * 1.001:
        we.append("encombrement maxi depasse (R = %.3f mm)" % res.r_eff)
    if min(norm(q) for q in res.profile) < p.r_min * 0.999:
        we.append("rayon mini non respecte")
    res.warnings = w + we
    res.valid = not res.warnings


def synthesize(p, mode="retour", T0=None, side=1.0):
    """Synthese complete : trouve T0 tel que le rayon maxi atteigne r_max."""
    face = make_face(p)
    gp0 = dot(t_of(p.gamma0), p.pivot)       # dH/dphi au repos
    if abs(gp0) < 1e-9:
        r = Result()
        r.message = "normale de contact passant par le pivot : arc-boutement"
        return r
    eps = sgn(gp0)

    if T0 is not None:
        return _integrate(p, face, T0, eps, mode, side)

    if face.kind == "plan":
        h_geo = norm(p.pivot) - dot(u_of(p.gamma0), p.pivot) + p.r_min
        if h_geo < p.r_max:
            r = Result()
            r.message = ("pente du plan insuffisante : distance maxi atteignable "
                         "%.2f mm < encombrement %.2f mm (augmenter |pente| ou "
                         "eloigner le pivot)" % (h_geo, p.r_max))
            return r

    cache = {}

    def run(T):
        key = round(T, 12)
        if key not in cache:
            cache[key] = _integrate(p, face, T, eps, mode, side)
        return cache[key]

    def err(T):
        """r_eff - r_max, ou None si non exploitable."""
        r = run(T)
        if not (r.ok and r.geom_ok):
            return None
        return r.r_eff - p.r_max

    C0 = max(p.couple_min, 1e-12)
    # 1) encadrement : r_eff croit avec T0, on part petit et on double
    lo, elo = None, None
    T = 1e-4 * C0
    for _ in range(14):
        e = err(T)
        if e is not None and e < 0.0:
            lo, elo = T, e
            break
        T *= 0.25 if e is None else 0.25
    if lo is None:
        # balayage de secours
        best = None
        for i in range(60):
            T = (1e-6 * C0) * ((10.0 * C0) / (1e-6 * C0)) ** (i / 59.0)
            r = run(T)
            if r.ok and (best is None
                         or abs(r.r_eff - p.r_max) < abs(best.r_eff - p.r_max)):
                best = r
        if best is None:
            r = Result()
            r.message = "aucune solution : verifier pivot, pente et couple"
            return r
        if not best.message and not best.valid:
            best.message = "solution non valide : " + " ; ".join(best.warnings)
        elif not best.message:
            best.message = ("encombrement cible non atteint : R = %.3f mm "
                            "(cible %.3f mm)" % (best.r_eff, p.r_max))
        return best

    hi, ehi = None, None
    T = lo
    for _ in range(80):
        Tn = T * 1.8
        e = err(Tn)
        if e is None:
            # la solution cesse d'etre exploitable : on raffine entre lo et Tn
            a, b = lo, Tn
            for _ in range(30):
                m = 0.5 * (a + b)
                em = err(m)
                if em is None:
                    b = m
                else:
                    a, elo = m, em
                    lo = m
                    if em > 0.0:
                        hi, ehi = m, em
                        break
            break
        T = Tn
        if e > 0.0:
            hi, ehi = T, e
            break
        lo, elo = T, e
    if hi is None:
        best = run(lo)
        best.message = ("encombrement cible non atteignable : R = %.3f mm "
                        "(cible %.3f mm) - revoir pente / pivot / rayons"
                        % (best.r_eff, p.r_max))
        return best

    # 2) secante avec repli bissection
    r = run(lo)
    for _ in range(60):
        if abs(ehi - elo) > 1e-14:
            Tm = lo - elo * (hi - lo) / (ehi - elo)
        else:
            Tm = 0.5 * (lo + hi)
        if not (lo < Tm < hi):
            Tm = 0.5 * (lo + hi)
        em = err(Tm)
        if em is None:
            hi = 0.5 * (lo + hi)
            ehi = err(hi)
            if ehi is None:
                break
            continue
        r = run(Tm)
        if abs(em) < 1e-7 * max(p.r_max, 1.0) or (hi - lo) < 1e-13 * C0:
            break
        if em < 0.0:
            lo, elo = Tm, em
        else:
            hi, ehi = Tm, em
    if r.ok and not r.valid and not r.message:
        r.message = "solution non valide : " + " ; ".join(r.warnings)
    return r


# ---------------------------------------------------------------- cinematique

def _demi_requis(p, res):
    """Demi-longueur de plan que le profil exige reellement.

    La synthese analytique suppose le plan **assez long** : la tangence est
    cherchee sans borne d'abscisse. Il faut donc savoir jusqu'ou le point de
    contact s'eloigne du milieu du plan, d'abord pour tracer une face qui
    touche vraiment le coeur, ensuite pour avertir quand le plan reel est plus
    court que cela.
    """
    face = make_face(p)
    if face.kind != "plan":
        return 0.0
    u0 = u_of(p.gamma0)
    M0 = (p.r_min * u0[0], p.r_min * u0[1])
    pire = 0.0
    for i in range(len(res.theta)):
        phi = res.phi[i]
        gam = res.gamma[i]
        Mf = _tourner_autour(M0, p.pivot, phi)
        X = rot(res.contact[i], res.theta[i])
        tg = t_of(gam)
        pire = max(pire, abs(tg[0] * (X[0] - Mf[0]) + tg[1] * (X[1] - Mf[1])))
    return pire


def pose(p, res, index):
    """Elements a tracer pour l'image `index` de la cinematique."""
    face = make_face(p)
    i = max(0, min(len(res.theta) - 1, index))
    th, phi = res.theta[i], res.phi[i]
    prof = [rot(q, th) for q in res.closed_profile(p.n_lobes)]
    # Plan reel si l'utilisateur l'a donne, sinon assez long pour couvrir la
    # tangence : la tracer plus courte que le contact faisait quitter le coeur
    # au marteau alors que la synthese le suppose toujours tangent.
    half = (p.demi_plan() if p.longueur_plan > 1e-9
            else max(0.7 * p.r_max, res.demi_requis + 0.1))
    return {
        "theta": th,
        "phi": phi,
        "profile": prof,
        "face": face.segment(phi, half),
        "contact": rot(res.contact[i], th),
        "normal": u_of(res.gamma[i]),
        "H": res.H[i],
        "N": res.N[i],
        "Th": res.Th_retour[i],
        "Cm": res.Cm[i],
    }
