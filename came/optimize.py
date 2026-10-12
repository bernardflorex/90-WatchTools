"""Optimisation des parametres du coeur (Nelder-Mead, Python pur)."""

import math

from .core import Params, synthesize


class Weights:
    """Ponderations du critere."""

    SPECS = [
        ("w_dcm",    "Augmentation du couple"),
        ("w_couple", "Couple utile du coeur"),
        ("w_ripple", "Ecart rentrant/sortant"),
        ("w_force",  "Effort de contact"),
        ("w_course", "Course du marteau"),
    ]

    def __init__(self, **kw):
        self.w_dcm = 1.0        # minimiser dCm / C0  (objectif principal)
        self.w_couple = 1.0     # maximiser T0 / C0
        self.w_ripple = 0.5     # minimiser l'ecart entre frottement rentrant/sortant
        self.w_force = 0.2      # minimiser l'effort normal maxi
        self.w_course = 0.1     # minimiser la course du marteau
        self.__dict__.update(kw)


#: variables optimisables : nom -> (libelle, mini, maxi)
VARIABLES = {
    "angle_plan": ("Pente du plan [deg]", 20.0, 88.0),
    "pente":      ("Pente marteau k [N.mm/rad]", 0.0, 100.0),
    "rayon_plan": ("Rayon du plan [mm]", 0.0, 200.0),
    "pivot_x":    ("Pivot X [mm]", -60.0, 60.0),
    "pivot_y":    ("Pivot Y [mm]", -60.0, 60.0),
}

PENALTY = 1e6


def cost(p, w, mode="retour"):
    """Critere scalaire a minimiser. Renvoie (cout, resultat)."""
    try:
        r = synthesize(p, mode=mode)
    except Exception:
        return PENALTY * 10, None
    if not r.ok:
        return PENALTY * 10, None
    c = 0.0
    c0 = max(p.couple_min, 1e-9)
    c += w.w_dcm * (r.delta_Cm / c0)
    c += w.w_couple * (1.0 / (1.0 + 20.0 * abs(r.T0) / c0))
    c += w.w_ripple * min(r.ripple, 50.0)
    c += w.w_force * (r.N_max * c0 / max(abs(r.T0), 1e-9)) * 1e-3
    c += w.w_course * (r.psi_max / math.radians(90.0))
    if not r.geom_ok:
        c += PENALTY
    # ecart a l'encombrement cible (le profil doit exploiter la place dispo)
    c += 5.0 * abs(r.r_eff - p.r_max) / max(p.r_max, 1e-9)
    if r.jam_min < 0.5:
        c += 10.0 * (0.5 - r.jam_min)
    return c, r


def _nelder_mead(f, x0, lo, hi, maxiter=300, tol=1e-6):
    n = len(x0)
    def clamp(x):
        return [min(hi[i], max(lo[i], x[i])) for i in range(n)]
    step = [max(1e-6, 0.08 * (hi[i] - lo[i])) for i in range(n)]
    simplex = [clamp(list(x0))]
    for i in range(n):
        y = list(x0)
        y[i] += step[i]
        simplex.append(clamp(y))
    fv = [f(s) for s in simplex]
    for _ in range(maxiter):
        order = sorted(range(n + 1), key=lambda i: fv[i])
        simplex = [simplex[i] for i in order]
        fv = [fv[i] for i in order]
        if abs(fv[-1] - fv[0]) < tol * (abs(fv[0]) + tol):
            break
        centroid = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
        def pt(coef):
            return clamp([centroid[j] + coef * (centroid[j] - simplex[-1][j])
                          for j in range(n)])
        xr = pt(1.0); fr = f(xr)
        if fr < fv[0]:
            xe = pt(2.0); fe = f(xe)
            simplex[-1], fv[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < fv[-2]:
            simplex[-1], fv[-1] = xr, fr
        else:
            xc = pt(-0.5); fc = f(xc)
            if fc < fv[-1]:
                simplex[-1], fv[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    simplex[i] = clamp([simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j])
                                        for j in range(n)])
                    fv[i] = f(simplex[i])
    k = min(range(n + 1), key=lambda i: fv[i])
    return simplex[k], fv[k]


def optimize(p, w, variables, mode="retour", restarts=3, progress=None,
             bounds=None):
    """Optimise `variables` (liste de noms) a partir de `p`.

    `bounds` : dict optionnel nom -> (mini, maxi) remplacant les bornes par defaut.
    Renvoie (params_optimises, resultat, cout, journal).
    """
    variables = [v for v in variables if v in VARIABLES]
    if not variables:
        c, r = cost(p, w, mode)
        return p, r, c, ["aucune variable selectionnee"]
    bounds = bounds or {}
    lo = [bounds.get(v, VARIABLES[v][1:])[0] for v in variables]
    hi = [bounds.get(v, VARIABLES[v][1:])[1] for v in variables]
    for i, v in enumerate(variables):
        if hi[i] <= lo[i]:
            hi[i] = lo[i] + 1e-6
    log = []
    n_eval = [0]

    def f(x):
        q = p.copy(**{variables[i]: x[i] for i in range(len(variables))})
        c, _ = cost(q, w, mode)
        n_eval[0] += 1
        if progress and n_eval[0] % 10 == 0:
            progress(n_eval[0], c)
        return c

    best_x, best_c = None, float("inf")
    starts = [[getattr(p, v) for v in variables]]
    for j in range(1, max(1, restarts)):
        starts.append([lo[i] + (hi[i] - lo[i]) * ((j * 0.37 + 0.2 * i) % 1.0)
                       for i in range(len(variables))])
    for s in starts:
        x, c = _nelder_mead(f, s, lo, hi)
        log.append("depart %s -> cout %.5g" % (
            ", ".join("%.3f" % v for v in s), c))
        if c < best_c:
            best_x, best_c = x, c
    q = p.copy(**{variables[i]: best_x[i] for i in range(len(variables))})
    c, r = cost(q, w, mode)
    log.append("optimum : " + ", ".join("%s = %.4f" % (variables[i], best_x[i])
                                        for i in range(len(variables))))
    log.append("cout final %.5g (%d evaluations)" % (c, n_eval[0]))
    return q, r, c, log
