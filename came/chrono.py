"""Coeur de remise a zero de compteur de chronographe.

Le **marteau est moteur** : son couple Cm = C0 + k.|psi| actionne le coeur, qui
restitue un couple constant Tc sur l'axe du mobile. Depuis n'importe quelle
position le coeur part d'un cote ou de l'autre : il a deux flancs.

Dissymetrie
-----------
D'un flanc a l'autre le coeur tourne en sens inverse, donc le glissement au
contact s'inverse : frottement **rentrant** d'un cote, **sortant** de l'autre.
Les deux flancs n'ont donc ni le meme profil ni la meme course.

Courses des deux flancs
-----------------------
Elles ne sont pas imposees, elles sont **calculees** par deux conditions :

1. la pointe est commune aux deux flancs et sur le diametre maxi ;
2. les deux sens de remise a zero donnent le meme couple sur le coeur.

C'est de la que vient le fait que la course jusqu'a la pointe ne vaut pas
180 deg et qu'elle depend du coefficient de frottement.

Plat du zero
------------
Au zero le plan du marteau est tangent au coeur. Les deux flancs demarrent tous
deux **sur cette tangente**, aux abscisses -s_A(0) et +s_B(0) : le segment qui
les separe EST le plat du zero. Sa longueur n'est donc pas libre :

    L_plat = |s_A(0)| + s_B(0)    avec  s = dH/dbeta

et comme s vaut en moyenne (r_max - r_min) / course, on a l'ordre de grandeur

    L_plat ~ 2 (r_max - r_min) / course[rad]

Raccourcir le plat a l'usinage ne raccourcit pas le plat vu par le marteau :
le plan du marteau enjambe le creux et porte sur les deux aretes. Le contact
saute alors d'une arete a l'autre, et c'est la que la resultante peut passer
par l'axe.

Non-blocage
-----------
Le bras de levier de la resultante par rapport a l'axe du coeur vaut

    d = |s + sigma.mu.H|

Il s'annule quand s = -sigma.mu.H : la resultante passe par l'axe, le couple de
rappel disparait, l'aiguille reste coincee. Au bord du plat s = +/-L et
H = r_min, d'ou la condition

    L > mu . r_min

Pointe et plage morte
---------------------
Au sommet, la normale de contact balaye la direction radiale : s passe par
zero, donc le couple s'annule. Le frottement elargit ce point en une plage ou
le coeur ne part d'aucun cote.

Cette plage est **inevitable** avec un marteau a face plane sur une came
convexe. Le bras de levier est |s|, le blocage a lieu pour |s| < mu.H, et la
vitesse a laquelle s traverse zero vaut ds/dbeta = rho_c - H, bornee par -H
puisque rho_c >= 0. La largeur minimale est donc

    2 . atan(mu)

atteinte pour rho_c = 0, c'est-a-dire une **pointe vive**. Mesures a
r_min 0.3 / r_max 3 / mu 0.12 (plancher 13.69 deg) :

    corde de pointe   0.05   0.40   0.72   1.20   1.80   2.40  mm
    plage morte      14.10  20.56  27.14  35.60  47.12  62.76  deg

Tout creux ou arrondi au sommet **elargit** la plage. La pente du plan du
marteau n'a pratiquement pas d'effet (35.6 a 36.0 deg entre 55 et 82 deg de
pente). Le seul levier est le coefficient de frottement : la plage vaut
2.atan(mu) et rien d'autre ne la reduit. En pratique c'est l'energie du
marteau qui la franchit.
"""

import math

from .core import (Params, dot, make_face, norm, rot, sgn, synthesize, t_of,
                   u_of)

TAU = 2.0 * math.pi


class ChronoParams(Params):
    """Parametres propres au coeur de chronographe."""

    EXTRA = [
        ("couple_coeur",   "Couple de rappel desire",  "N.mm", 0.001, 100.0, 0.05),
        ("corde_pointe",   "Corde de la pointe",       "mm",   0.02, 20.0, 0.05),
        ("rayon_pointe",   "Rayon des pointes",        "mm",   0.005, 2.0, 0.01),
        ("rayon_creux",    "Rayon du creux de pointe", "mm",   0.05, 50.0, 0.1),
        ("rayon_raccord",  "Rayon de raccord du plat", "mm",   0.0, 2.0, 0.01),
        ("plat_souhaite",  "Plat souhaite (controle)", "mm",   0.0, 20.0, 0.05),
        ("marge_blocage",  "Marge anti-blocage",       "-",    1.0, 20.0, 0.5),
    ]

    def __init__(self, **kw):
        self.couple_coeur = 0.3
        self.corde_pointe = 0.05
        self.rayon_pointe = 0.05
        self.rayon_creux = 1.0
        self.rayon_raccord = 0.08
        self.plat_souhaite = 0.6
        self.marge_blocage = 3.0
        # Le galet est abandonne. Ses deux grandeurs restent definies pour
        # que `profil.mu_galet` garde un sens si le sujet revient, mais elles
        # ne sont plus des champs du panneau.
        self.axe_galet = 0.15
        self.mu_palier = 0.12
        super().__init__(**kw)
        if "r_min" not in kw:
            self.r_min = 0.3
        if "r_max" not in kw:
            self.r_max = 3.0
        if "r_min" not in kw:
            # nominal donne par l'utilisateur : distance du plat a l'axe
            self.r_min = 0.6

    @property
    def pointe_mini(self):
        """Corde mini de la pointe pour qu'il n'y ait pas de plage morte."""
        return 2.0 * self.mu * self.r_max

    @property
    def plat_mini(self):
        """Demi-plat minimal pour que la resultante ne passe pas par l'axe."""
        return self.mu * self.r_min


# ------------------------------------------------------------------ resultat

class ChronoResult:
    def __init__(self):
        self.ok = False
        self.message = ""
        self.A = None              # Result du flanc A (sens direct)
        self.B = None              # Result du flanc B (sens inverse)
        self.course_a = 0.0        # deg
        self.course_b = 0.0        # deg
        self.ouverture_pointe = 0.0
        self.echelle = 1.0         # facteur applique a C0 et k
        self.couple_coeur = 0.0    # couple obtenu sur le coeur [N.mm]
        self.dissym = 0.0
        self.profile = []
        self.flat = []
        self.plat = 0.0            # longueur naturelle du plat [mm]
        self.plat_a = 0.0
        self.plat_b = 0.0
        self.plat_mini = 0.0
        self.corde = 0.0           # corde de la pointe (segment du sommet)
        self.corde_mini = 0.0      # corde mini pour supprimer la plage morte
        self.bras_pointe = 0.0     # bras de levier sur les rayons de pointe
        self.horns = []            # les deux sommets, repere du coeur
        self.notch = []            # arc du creux de pointe
        self.demi_requis = 0.0     # demi-longueur de plan exigee [mm]
        self.theta_saut = 0.0      # angle du coeur ou le contact saute [deg]
        self.total_a = 0.0         # course totale cote A (flanc + sommet)
        self.total_b = 0.0
        self.horn_a = []
        self.horn_b = []
        self.zone_morte = 0.0      # plage morte au sommet [deg]
        self.zone_morte_beta = 0.0 # idem, en balayage de normale [deg]
        self.bras_plat = 0.0       # bras de levier au bord du plat [mm]
        self.psi_max = 0.0
        self.N_max = 0.0
        self.Cm_max = 0.0
        self.warnings = []
        self.valid = False

    def flank(self, side):
        return self.A if side > 0 else self.B


# ------------------------------------------------------------- deux flancs

def _flank(p, side, course_deg, r_target):
    """Un flanc a couple constant, du zero jusqu'a la course donnee."""
    q = p.copy(sweep=course_deg, r_max=r_target)
    return synthesize(q, side=side)


def _end_angle(res):
    """Angle polaire du point d'extremite d'un flanc (repere du coeur)."""
    if not (res and res.ok):
        return None
    return res.beta[-1] + math.atan2(res.s[-1], res.H[-1])


def _secant(fn, x0, x1, lo, hi, tol, nmax=12):
    f0, f1 = fn(x0), fn(x1)
    if f0 is None or f1 is None:
        return x1 if f1 is not None else x0
    for _ in range(nmax):
        if abs(f1) < tol or abs(f1 - f0) < 1e-15:
            break
        x2 = min(hi, max(lo, x1 - f1 * (x1 - x0) / (f1 - f0)))
        if abs(x2 - x1) < 1e-5:
            x1 = x2
            break
        f2 = fn(x2)
        if f2 is None:
            x2 = 0.5 * (x1 + x2)
            f2 = fn(x2)
            if f2 is None:
                break
        x0, f0, x1, f1 = x1, f1, x2, f2
    return x1


def _close(p, course_a, r_target, hold):
    """Course du flanc B qui amene son extremite sur celle du flanc A."""
    A = _flank(p, +1.0, course_a, r_target)
    aA = _end_angle(A)
    if aA is None:
        return A, None, None
    target = aA + 2 * math.pi

    def g(cb):
        B = _flank(p, -1.0, cb, r_target)
        ab = _end_angle(B)
        if ab is None:
            return None
        hold["B"] = B
        hold["cb"] = cb
        return ab - target

    cb = _secant(g, 150.0, 185.0, 80.0, 280.0, 1e-7)
    B = hold.get("B")
    if B is None or abs(hold.get("cb", -1) - cb) > 1e-6:
        B = _flank(p, -1.0, cb, r_target)
    return A, B, cb


def _courses(p, r_target):
    """Resout (course A, course B) : pointe commune + couples egaux."""
    hold = {}
    best = {}

    def d(ca):
        A, B, cb = _close(p, ca, r_target, hold)
        if not (A and A.ok and B and B.ok):
            return None
        best["v"] = (A, B, ca, cb)
        return abs(A.T0) - abs(B.T0)

    ca = _secant(d, 170.0, 188.0, 110.0, 250.0, 1e-8, nmax=9)
    if "v" not in best or abs(best["v"][2] - ca) > 1e-6:
        A, B, cb = _close(p, ca, r_target, hold)
        if A and A.ok and B and B.ok:
            best["v"] = (A, B, ca, cb)
    return best.get("v")


# ---------------------------------------------------------------- la pointe

def _X(H, s, beta):
    u, t = u_of(beta), t_of(beta)
    return (H * u[0] + s * t[0], H * u[1] + s * t[1])


def _end_point(res):
    return _X(res.H[-1], res.s[-1], res.beta[-1])


def _end_polar(res):
    """Angle polaire de l'extremite d'un flanc (repere du coeur)."""
    return res.beta[-1] + math.atan2(res.s[-1], res.H[-1])


def _courses(p, progress=None):
    """Resout les deux courses.

    Conditions :
      - chaque flanc atteint le diametre maxi a son extremite (les deux
        sommets sont donc sur le diametre maxi) ;
      - les deux sommets sont ecartes de l'ouverture du creux de pointe :
        alpha_B = alpha_A + 360 deg - w, ou w est l'angle au centre de la
        corde demandee. C'est cet ecart qui supprime la plage morte : le
        contact saute d'un sommet a l'autre au lieu de passer par la normale
        qui traverse l'axe ;
      - les deux sens de remise a zero donnent le meme couple.
    """
    w = 2.0 * math.asin(min(0.999, 0.5 * p.corde_pointe / max(p.r_max, 1e-9)))
    hold = {}

    def pair(ca, cb):
        key = (round(ca, 7), round(cb, 7))
        if key not in hold:
            hold[key] = (_flank(p, +1.0, ca, p.r_max),
                         _flank(p, -1.0, cb, p.r_max))
        return hold[key]

    def close(ca):
        """Course B qui place le second sommet a l'ecart voulu du premier."""
        def h(cb):
            A, B = pair(ca, cb)
            if not (A.ok and B.ok):
                return None
            return (_end_polar(B) - _end_polar(A)) - (2 * math.pi - w)
        cb = _secant(h, 140.0, 175.0, 60.0, 300.0, 1e-9)
        A, B = pair(ca, cb)
        return A, B, cb

    best = {}

    def d(ca):
        A, B, cb = close(ca)
        if not (A.ok and B.ok):
            return None
        best["v"] = (A, B, ca, cb)
        if progress:
            progress("courses : A = %.2f deg, B = %.2f deg" % (ca, cb))
        return abs(A.T0) - abs(B.T0)

    ca = _secant(d, 168.0, 182.0, 100.0, 250.0, 1e-10, nmax=14)
    if "v" not in best or abs(best["v"][2] - ca) > 1e-7:
        A, B, cb = close(ca)
        if A.ok and B.ok:
            best["v"] = (A, B, ca, cb)
    return best.get("v")


def _horn_samples(p, res, side, beta_c, rho, scale, n=30):
    """Contact sur un sommet : le point de contact y est fige (coin entre le
    flanc et le creux de pointe), adouci au rayon rho.

    Le balayage va de la normale de fin de flanc a celle de la corde. Le bras
    de levier y vaut au moins corde/2 - rho, ce qui garantit un couple non nul
    sur ce rayon tant que corde/2 - rho > mu.r_max.
    """
    face = make_face(p)
    eps = sgn(dot(t_of(p.gamma0), p.pivot))
    mu, P = p.mu, p.pivot
    C0, k = p.couple_min * scale, p.pente * scale
    psi_max = abs(math.radians(p.angle_plan)) * 0.9995
    b0 = res.beta[-1]
    Pt = _end_point(res)
    bis = 0.5 * (b0 + beta_c)
    M = (Pt[0] - rho * math.cos(bis), Pt[1] - rho * math.sin(bis))

    def phi_of_H(H):
        f = lambda psi: dot(u_of(p.gamma0 + eps * psi), P) + face.h0 - H
        if f(0.0) > 0:
            return 0.0
        if f(psi_max) < 0:
            return None
        lo, hi = 0.0, psi_max
        for _ in range(70):
            m = 0.5 * (lo + hi)
            if f(m) < 0:
                lo = m
            else:
                hi = m
        return eps * 0.5 * (lo + hi)

    out = []
    for i in range(n + 1):
        beta = b0 + (beta_c - b0) * i / n
        H = dot(u_of(beta), M) + rho
        sv = dot(t_of(beta), M)
        phi = phi_of_H(H)
        if phi is None:
            continue
        theta = side * (p.gamma0 + phi - beta)
        e = face.h0
        q = sv - dot(t_of(p.gamma0 + phi), P)
        Cm = -eps * (C0 + k * abs(phi))

        def torque(sig):
            den = q + sig * mu * e
            if abs(den) < 1e-12:
                return float("nan"), float("nan")
            return Cm * (sv + sig * mu * H) / den, abs(Cm / den)

        # le coeur peut partir vers A (theta_coeur decroit) ou vers B ;
        # le frottement s'oppose au glissement, donc sigma differe
        ThA, NA = torque(+1.0)
        ThB, NB = torque(-1.0)
        goA = ThA == ThA and ThA < 0.0
        goB = ThB == ThB and ThB > 0.0
        if goA and not goB:
            Th, N = ThA, NA
        elif goB and not goA:
            Th, N = ThB, NB
        elif goA and goB:
            Th, N = (ThA, NA) if abs(ThA) >= abs(ThB) else (ThB, NB)
        else:
            Th = 0.0
            N = NA if NA == NA else NB
        out.append({"beta": beta, "theta": theta, "phi": phi, "H": H, "s": sv,
                    "N": N if N == N else 0.0, "Th": Th,
                    "bras": abs(Th) / N if N and N == N else 0.0,
                    "bloque": not (goA or goB)})
    return out


def horn_arm(p, corde):
    """Bras de levier minimal garanti sur le rayon d'un sommet.

    Le flanc arrive tangent a la corde : le sommet est une pointe que l'on
    adoucit au rayon r_pointe. Le contact y reste a moins de r_pointe du
    sommet, donc |s| >= corde/2 - r_pointe, et le couple s'annulerait pour
    |s| = mu.H = mu.r_max. D'ou le bras garanti :

        corde/2 - r_pointe - mu.r_max
    """
    return 0.5 * corde - p.rayon_pointe - p.mu * p.r_max


def _notch(p, PA, PB, rayon, n=28):
    """Arc concave reliant les deux sommets (creux de pointe)."""
    cx, cy = 0.5 * (PA[0] + PB[0]), 0.5 * (PA[1] + PB[1])
    dx, dy = PB[0] - PA[0], PB[1] - PA[1]
    L = math.hypot(dx, dy)
    if L < 1e-9:
        return []
    R = max(rayon, 0.5 * L * 1.0001)
    h = math.sqrt(max(0.0, R * R - 0.25 * L * L))
    # centre du cote exterieur : l'arc creuse vers l'interieur
    ux, uy = dx / L, dy / L
    nx, ny = -uy, ux
    if nx * cx + ny * cy < 0:
        nx, ny = -nx, -ny
    C = (cx + h * nx, cy + h * ny)
    a0 = math.atan2(PA[1] - C[1], PA[0] - C[0])
    a1 = math.atan2(PB[1] - C[1], PB[0] - C[0])
    while a1 - a0 > math.pi:
        a1 -= 2 * math.pi
    while a1 - a0 < -math.pi:
        a1 += 2 * math.pi
    return [(C[0] + R * math.cos(a0 + (a1 - a0) * i / n),
             C[1] + R * math.sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]


# ------------------------------------------------------------------ profil

def _build_profile(p, r):
    pts = []
    A, B = r.A, r.B
    for i in range(len(A.theta) - 1, -1, -1):
        pts.append(_X(A.H[i], A.s[i], A.beta[i]))
    pts.append(_X(p.r_min, A.s[0], p.gamma0))
    pts.append(_X(p.r_min, B.s[0], p.gamma0))
    r.flat = [pts[-2], pts[-1]]
    for i in range(len(B.theta)):
        pts.append(_X(B.H[i], B.s[i], B.beta[i]))
    if r.notch:
        pts.extend(r.notch)
    return pts


# ------------------------------------------------------------------ synthese

def _pourquoi_pas(p):
    """Dit ce qui bloque, au lieu de « verifier pivot, pente et rayons ».

    La cause de loin la plus fréquente est geometrique et se calcule : avec un
    plan plat, le marteau ne peut eloigner son plan de l'axe du coeur que de
    `|OP|(1 - cos a) + r_min`. Si cela ne couvre pas l'encombrement demande,
    aucun profil n'existe, quelle que soit la forme du coeur. C'est le cas
    apres l'import d'un dessin dont le pivot et la pente ne viennent pas d'une
    synthese.
    """
    from .core import u_of, dot, norm
    h_geo = norm(p.pivot) - dot(u_of(p.gamma0), p.pivot) + p.r_min
    if h_geo < p.r_max:
        return ("pente du plan insuffisante : le marteau ne peut ecarter son "
                "plan que de %.3f mm de l'axe, il faut atteindre %.3f mm. "
                "Augmenter |pente du plan|, eloigner le pivot, ou reduire "
                "l'encombrement." % (h_geo, p.r_max))
    if p.r_max <= p.r_min + 1e-9:
        return "encombrement maxi au plus egal au rayon mini"
    plat = p.plat_mini
    return ("pas de solution : la geometrie tient (ecart maxi atteignable "
            "%.3f mm pour %.3f demandes), mais l'integration n'aboutit pas. "
            "Essayer une autre pente de plan, un autre pivot, ou un plat de "
            "zero plus grand que %.4f mm." % (h_geo, p.r_max, plat))


def synthesize_chrono(p, equilibrer=True, progress=None):
    """Synthese complete du coeur de chronographe."""
    r = ChronoResult()
    r.plat_mini = p.plat_mini
    r.corde_mini = p.pointe_mini

    fast = p.copy(n_steps=max(90, p.n_steps // 3))
    sol = _courses(fast, progress)
    if sol is None:
        r.message = _pourquoi_pas(p)
        return r
    _a, _b, ca, cb = sol
    A, B = _flank(p, +1.0, ca, p.r_max), _flank(p, -1.0, cb, p.r_max)
    if not (A.ok and B.ok):
        A, B = _a, _b
    r.A, r.B = A, B
    r.course_a, r.course_b = ca, cb

    PA, PB = _end_point(A), _end_point(B)
    r.horns = [PA, PB]
    r.corde = norm((PB[0] - PA[0], PB[1] - PA[1]))
    r.ouverture_pointe = math.degrees(2.0 * math.asin(
        min(1.0, 0.5 * r.corde / max(p.r_max, 1e-9))))
    r.bras_pointe = 0.5 * r.corde - p.mu * p.r_max
    r.notch = _notch(p, PB, PA, p.rayon_creux)

    # echelle des couples : le profil ne depend que de T0/C0 et k/C0
    t_moy = 0.5 * (abs(A.T0) + abs(B.T0))
    r.echelle = p.couple_coeur / t_moy if t_moy > 1e-12 else 1.0
    r.couple_coeur = t_moy * r.echelle
    r.dissym = abs(abs(A.T0) - abs(B.T0)) / max(t_moy, 1e-12)

    # plat du zero
    r.plat_a, r.plat_b = abs(A.s[0]), abs(B.s[0])
    r.plat = r.plat_a + r.plat_b
    r.bras_plat = min(r.plat_a, r.plat_b) - p.mu * p.r_min

    # sommets : coins entre les flancs et le creux, adoucis au rayon r_pointe
    if progress:
        progress("sommets et creux de pointe")
    # normale de la corde : les deux sommets encadrent le secteur de pointe,
    # qui est parcouru apres le flanc B et avant le flanc A du tour suivant
    beta_c = 0.5 * (_end_polar(A) + 2 * math.pi + _end_polar(B))
    hA = _horn_samples(p, A, +1.0, beta_c - 2 * math.pi, p.rayon_pointe,
                       r.echelle)
    hB = _horn_samples(p, B, -1.0, beta_c, p.rayon_pointe, r.echelle)
    r.horn_a, r.horn_b = hA, hB
    r.total_a = math.degrees(max(x["theta"] for x in hA)) if hA else ca
    r.total_b = (math.degrees(max(abs(x["theta"]) for x in hB)) if hB else cb)
    r.theta_saut = r.total_a
    r.bras_pointe = (min(x["bras"] for x in hA + hB) if hA + hB
                     else horn_arm(p, r.corde))
    bl = [x for x in hA + hB if x["bloque"]]
    if len(bl) > 1:
        th = [math.degrees(x["theta"]) if x in hA
              else 360.0 - abs(math.degrees(x["theta"])) for x in bl]
        r.zone_morte = max(th) - min(th)
        r.zone_morte_beta = sum(
            math.degrees(max(x["beta"] for x in h) - min(x["beta"] for x in h))
            for h in ([x for x in bl if x in hA], [x for x in bl if x in hB])
            if len(h) > 1)
    else:
        r.zone_morte = r.zone_morte_beta = 0.0

    r.profile = _build_profile(p, r)

    psis = [abs(x) for x in A.phi] + [abs(x) for x in B.phi]
    psis += [abs(x["phi"]) for x in hA + hB]
    Ns = [x * r.echelle for x in A.N + B.N] + [x["N"] for x in hA + hB]
    r.psi_max = max(psis) if psis else 0.0
    r.N_max = max(Ns) if Ns else 0.0
    r.Cm_max = (p.couple_min + p.pente * r.psi_max) * r.echelle

    # diagnostic
    w = []
    if r.bras_plat <= 0.0:
        w.append("plat trop court : demi-plat %.4f mm <= mu.r_min = %.4f mm, "
                 "la resultante passe par l'axe au bord du plat"
                 % (min(r.plat_a, r.plat_b), p.plat_mini))
    elif min(r.plat_a, r.plat_b) < p.marge_blocage * p.plat_mini:
        w.append("marge faible au bord du plat (%.2f x le mini)"
                 % (min(r.plat_a, r.plat_b) / max(p.plat_mini, 1e-12)))
    plancher = math.degrees(2.0 * math.atan(p.mu))
    if r.zone_morte > 1.35 * plancher:
        w.append("plage morte au sommet %.1f deg, loin du plancher %.1f deg "
                 "(2.atan(mu)) : reduire la corde de pointe"
                 % (r.zone_morte, plancher))
    if p.plat_souhaite > 0 and abs(r.plat - p.plat_souhaite) > 0.1 * max(
            p.plat_souhaite, 1e-9):
        w.append("plat impose par la geometrie : %.3f mm, souhaite %.3f mm "
                 "(L_plat ~ 2(r_max - r_min)/course)"
                 % (r.plat, p.plat_souhaite))
    rc = min(A.curv[0], B.curv[0])
    if p.rayon_raccord > rc:
        w.append("rayon de raccord du plat (%.3f mm) superieur au rayon de "
                 "courbure du flanc a son depart (%.3f mm)"
                 % (p.rayon_raccord, rc))
    if abs(r.corde - p.corde_pointe) > 0.05 * max(p.corde_pointe, 1e-9):
        w.append("corde de pointe obtenue %.3f mm au lieu de %.3f mm demandee"
                 % (r.corde, p.corde_pointe))

    for res, nm in ((A, "A"), (B, "B")):
        if not res.geom_ok:
            w.append("flanc %s : %s" % (nm, " ; ".join(res.warnings)))
    if not (144.0 <= r.total_a <= 216.0):
        w.append("course totale cote A hors de 180 deg +/- 20 %% (%.1f deg)"
                 % r.total_a)
    if not (144.0 <= r.total_b <= 216.0):
        w.append("course totale cote B hors de 180 deg +/- 20 %% (%.1f deg)"
                 % r.total_b)
    r.warnings = w
    r.valid = not w
    # demi-longueur de plan exigee par les deux flancs
    from .core import _demi_requis
    r.demi_requis = max(_demi_requis(p, x) for x in (r.A, r.B) if x and x.ok) \
        if (r.A and r.A.ok) or (r.B and r.B.ok) else 0.0
    r.ok = True
    return r


# ------------------------------------------------------- courbes et poses

def curves(p, r):
    """Courbes continues sur 360 deg : (x [deg], couple, bras, N, ouverture).

    x va de -course_B a +course_A, les rayons de sommet se trouvant aux deux
    extremites. Le contact saute d'un sommet a l'autre : le couple change de
    signe sans passer par zero, il n'y a donc pas de plage morte.
    """
    x, Th, bras, N, psi = [], [], [], [], []
    k = r.echelle

    def push(xx, th, br, nn, ph):
        x.append(xx); Th.append(abs(th)); bras.append(br)
        N.append(nn); psi.append(abs(math.degrees(ph)))

    for q in sorted(r.horn_b, key=lambda z: -abs(z["theta"])):
        push(-abs(math.degrees(q["theta"])), q["Th"], q["bras"], q["N"],
             q["phi"])
    B = r.B
    if B and B.ok:
        for i in range(len(B.theta) - 1, -1, -1):
            n = B.N[i] * k
            push(-math.degrees(B.theta[i]), B.Th_retour[i] * k,
                 abs(B.Th_retour[i] * k) / n if n else 0.0, n, B.phi[i])
    push(0.0, 0.0, min(r.plat_a, r.plat_b) - p.mu * p.r_min, 0.0, 0.0)
    A = r.A
    if A and A.ok:
        for i in range(len(A.theta)):
            n = A.N[i] * k
            push(math.degrees(A.theta[i]), A.Th_retour[i] * k,
                 abs(A.Th_retour[i] * k) / n if n else 0.0, n, A.phi[i])
    for q in sorted(r.horn_a, key=lambda z: z["theta"]):
        push(math.degrees(q["theta"]), q["Th"], q["bras"], q["N"], q["phi"])
    return x, Th, bras, N, psi


def pose_at(p, r, x_deg):
    """Pose du mecanisme pour l'abscisse x (deg) des courbes ci-dessus."""
    face = make_face(p)
    side = 1.0 if x_deg >= 0 else -1.0
    a = abs(x_deg)
    phi, H, s, beta = 0.0, p.r_min, 0.0, p.gamma0
    res = r.flank(side)
    cand = []
    if res and res.ok:
        cand = [(math.degrees(res.theta[i]), res.phi[i], res.H[i], res.s[i],
                 res.beta[i]) for i in range(len(res.theta))]
    cand += [(abs(math.degrees(q["theta"])), q["phi"], q["H"], q["s"],
              q["beta"]) for q in (r.horn_a if side > 0 else r.horn_b)]
    if cand:
        t, phi, H, s, beta = min(cand, key=lambda q: abs(q[0] - a))
        a = t
    contact = rot(_X(H, s, beta), side * math.radians(a))
    prof = [rot(q, side * math.radians(a)) for q in r.profile]
    # plan reel si donne, sinon assez long pour couvrir la tangence
    half = (p.demi_plan() if p.longueur_plan > 1e-9
            else max(0.7 * p.r_max, r.demi_requis + 0.1))
    return {"theta": side * math.radians(a), "phi": phi, "profile": prof,
            "face": face.segment(phi, half),
            "contact": contact, "normal": u_of(p.gamma0 + phi), "H": H}


# ------------------------------------------------ fonction d'appui du coeur

def support(p, r):
    """Fonction d'appui du coeur sur un tour : (betas, H, s) croissants.

    Reconstruite a partir des deux flancs, du plat du zero et du coin de
    pointe (ou le point de contact est fige). C'est le point de depart de la
    discretisation du profil (voir profil.from_chrono).
    """
    A, B = r.A, r.B
    bet, Hs, ss = [], [], []

    def add(b, H, s):
        bet.append(b); Hs.append(H); ss.append(s)

    for i in range(len(A.beta) - 1, -1, -1):          # flanc A : beta croissant
        add(A.beta[i], A.H[i], A.s[i])
    for i in range(len(B.beta)):                      # flanc B
        if B.beta[i] <= bet[-1] + 1e-12:
            continue
        add(B.beta[i], B.H[i], B.s[i])
    # coin de pointe : contact fige sur le sommet
    u0, t0 = u_of(A.beta[-1]), t_of(A.beta[-1])
    Xt = (A.H[-1] * u0[0] + A.s[-1] * t0[0],
          A.H[-1] * u0[1] + A.s[-1] * t0[1])
    b0, b1 = bet[-1], bet[0] + TAU
    n = 180
    for i in range(1, n + 1):
        b = b0 + (b1 - b0) * i / n
        add(b, dot(u_of(b), Xt), dot(t_of(b), Xt))
    return bet, Hs, ss
