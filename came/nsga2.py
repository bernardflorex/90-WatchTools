"""NSGA-II (Python pur) et son application au profil de coeur.

Pourquoi un algorithme multi-objectif ici : les criteres du coeur sont
reellement contradictoires et sans taux d'echange evident.

  - plage morte au sommet         <-> ondulation du couple
  - couple utile (rendement)      <-> effort de contact, course du marteau
  - encombrement                  <-> longueur du plat, courbure usinable

Avec la synthese analytique on ne peut pas arbitrer : le couple constant est
impose, le profil en decoule, et on ne peut que rejeter ce qui ne va pas.
En prenant le profil lui-meme comme variable (voir profil.py) on obtient un
front de Pareto et on choisit ensuite.

La population de depart est semee avec la solution analytique, qui est deja
un bon point : l'optimiseur part d'un coeur valide et l'ameliore, au lieu de
chercher la convexite et la fermeture au hasard.
"""

import math
import random

from .core import norm
from .profil import (Profil, dimensionner_ressort, evaluer, longueur_plat,
                     mu_galet, reparer, secteur_faible)


# ------------------------------------------------------------------ moteur

def domine(a, b):
    """a domine b (tous les objectifs <=, au moins un <)."""
    if a["viol"] > 0.0 or b["viol"] > 0.0:
        if a["viol"] < b["viol"]:
            return True
        if a["viol"] > b["viol"]:
            return False
    f, g = a["f"], b["f"]
    return all(x <= y for x, y in zip(f, g)) and any(x < y for x, y in zip(f, g))


def fronts_non_domines(pop):
    S = [[] for _ in pop]
    n = [0] * len(pop)
    fronts = [[]]
    for i, pi in enumerate(pop):
        for j, pj in enumerate(pop):
            if i == j:
                continue
            if domine(pi, pj):
                S[i].append(j)
            elif domine(pj, pi):
                n[i] += 1
        if n[i] == 0:
            pi["rang"] = 0
            fronts[0].append(i)
    k = 0
    while fronts[k]:
        suiv = []
        for i in fronts[k]:
            for j in S[i]:
                n[j] -= 1
                if n[j] == 0:
                    pop[j]["rang"] = k + 1
                    suiv.append(j)
        k += 1
        fronts.append(suiv)
    return fronts[:-1]


def distance_entassement(pop, front):
    for i in front:
        pop[i]["dist"] = 0.0
    if len(front) <= 2:
        for i in front:
            pop[i]["dist"] = float("inf")
        return
    m = len(pop[front[0]]["f"])
    for k in range(m):
        order = sorted(front, key=lambda i: pop[i]["f"][k])
        pop[order[0]]["dist"] = pop[order[-1]]["dist"] = float("inf")
        lo, hi = pop[order[0]]["f"][k], pop[order[-1]]["f"][k]
        if hi - lo < 1e-15:
            continue
        for t in range(1, len(order) - 1):
            pop[order[t]]["dist"] += (pop[order[t + 1]]["f"][k]
                                      - pop[order[t - 1]]["f"][k]) / (hi - lo)


def _tournoi(pop, rng):
    a, b = rng.randrange(len(pop)), rng.randrange(len(pop))
    pa, pb = pop[a], pop[b]
    if pa["rang"] != pb["rang"]:
        return pa if pa["rang"] < pb["rang"] else pb
    return pa if pa["dist"] > pb["dist"] else pb


def _sbx(x1, x2, lo, hi, rng, eta=15.0, pc=0.9):
    y1, y2 = list(x1), list(x2)
    if rng.random() > pc:
        return y1, y2
    for i in range(len(y1)):
        if rng.random() > 0.5:
            continue
        u = rng.random()
        beta = (2 * u) ** (1 / (eta + 1)) if u <= 0.5 \
            else (1 / (2 * (1 - u))) ** (1 / (eta + 1))
        a, b = y1[i], y2[i]
        y1[i] = min(hi[i], max(lo[i], 0.5 * ((1 + beta) * a + (1 - beta) * b)))
        y2[i] = min(hi[i], max(lo[i], 0.5 * ((1 - beta) * a + (1 + beta) * b)))
    return y1, y2


def _mutation(x, lo, hi, rng, eta=20.0, pm=None):
    n = len(x)
    pm = pm if pm is not None else 1.0 / n
    y = list(x)
    for i in range(n):
        if rng.random() > pm:
            continue
        d = hi[i] - lo[i]
        if d <= 0:
            continue
        u = rng.random()
        delta = (2 * u) ** (1 / (eta + 1)) - 1 if u < 0.5 \
            else 1 - (2 * (1 - u)) ** (1 / (eta + 1))
        y[i] = min(hi[i], max(lo[i], y[i] + delta * d))
    return y


def nsga2(evaluer_x, lo, hi, n_pop=60, n_gen=30, graines=None, reparer_x=None,
          rng=None, progress=None, evaluer_lot=None):
    """NSGA-II generique. `evaluer_x(x) -> (objectifs, violation)`.

    `evaluer_lot(liste_de_x) -> liste de (x_repare, objectifs, violation)`
    permet d'evaluer toute une generation d'un coup, par exemple sur
    plusieurs coeurs (voir `pool_lot`). Les individus d'une generation sont
    independants, c'est le seul endroit ou le parallelisme est utile.
    """
    rng = rng or random.Random(12345)
    n = len(lo)

    def lot(xs):
        if evaluer_lot is not None:
            return [{"x": x, "f": f, "viol": v, "rang": 0, "dist": 0.0}
                    for x, f, v in evaluer_lot(xs)]
        out = []
        for x in xs:
            if reparer_x:
                x = reparer_x(x)
            f, v = evaluer_x(x)
            out.append({"x": x, "f": f, "viol": v, "rang": 0, "dist": 0.0})
        return out

    depart = [list(g) for g in (graines or [])]
    while len(depart) < n_pop:
        base = rng.choice(graines) if graines else None
        if base is not None and rng.random() < 0.8:
            x = [min(hi[i], max(lo[i],
                     base[i] + rng.gauss(0, 0.08 * (hi[i] - lo[i]))))
                 for i in range(n)]
        else:
            x = [rng.uniform(lo[i], hi[i]) for i in range(n)]
        depart.append(x)
    pop = lot(depart[:max(n_pop, len(graines or []))])

    for gen in range(n_gen):
        fr = fronts_non_domines(pop)
        for f in fr:
            distance_entassement(pop, f)
        bruts = []
        while len(bruts) < n_pop:
            p1, p2 = _tournoi(pop, rng), _tournoi(pop, rng)
            c1, c2 = _sbx(p1["x"], p2["x"], lo, hi, rng)
            bruts.append(_mutation(c1, lo, hi, rng))
            if len(bruts) < n_pop:
                bruts.append(_mutation(c2, lo, hi, rng))
        melange = pop + lot(bruts)
        fr = fronts_non_domines(melange)
        nouvelle = []
        for f in fr:
            distance_entassement(melange, f)
            if len(nouvelle) + len(f) <= n_pop:
                nouvelle.extend(melange[i] for i in f)
            else:
                reste = sorted(f, key=lambda i: -melange[i]["dist"])
                nouvelle.extend(melange[i] for i in reste[:n_pop - len(nouvelle)])
                break
        pop = nouvelle
        if progress:
            faisables = [q for q in pop if q["viol"] <= 0.0]
            progress(gen + 1, n_gen, len(faisables), pop)
    fr = fronts_non_domines(pop)
    for f in fr:
        distance_entassement(pop, f)
    return pop, [pop[i] for i in fr[0]]


# ------------------------------------------------- le probleme du coeur

#: noms des objectifs, dans l'ordre
OBJECTIFS = ("tour sous la consigne [%]", "plage morte [deg]",
             "ondulation du couple [%]", "effort de contact [N]")

#: Sens de calcul : on connait l'encombrement, le pivot du marteau, la pente
#: de son plan et le **couple de rappel voulu** sur le coeur. Le ressort du
#: marteau (prearmage C0 et raideur k) est une **sortie** : il est calcule
#: pour que le couple vise soit atteint partout. Toutes les solutions
#: delivrent donc le meme couple, et on les compare sur ce que coute ce
#: resultat : ondulation, plage morte, effort, taille du ressort.


#: rayon de face mini en valeur absolue (au-dela on prend la face plane)
R_FACE_MINI = 1.2

#: ecart a la tangente sous lequel on compte un point comme etant sur le plat
TOL_PLAT = 0.002

#: tolerance sur r_min et r_max [mm]. Un depassement de quelques microns est
#: un residu de reconstruction, pas un defaut de conception : le refuser
#: ecartait la graine analytique au profit de solutions bien plus mauvaises.
TOL_ENCOMBREMENT = 0.005



#: grandeurs du marteau que l'optimisation peut faire varier.
#:
#: `pivot_t` n'est pas un champ du panneau mais un **deplacement du pivot
#: perpendiculairement a O->P**. La mesure est sans appel : sur +/- 4 mm,
#: deplacer le pivot radialement fait varier l'ondulation de 1.6 point
#: (41.3 a 42.9 %), tangentiellement de 68.8 points (13.5 a 82.3 %). La
#: distance a l'axe n'apporte quasiment rien ; c'est la position **angulaire**
#: du pivot qui compte. Une seule variable suffit donc la ou `pivot_x` et
#: `pivot_y` en demandaient deux, dont une inutile.
MARTEAU_LIBRE = ("longueur_plan", "decalage_plan", "rayon_coin",
                 "angle_plan", "pivot_t", "pivot_x", "pivot_y")

#: celles que l'on conseille : le pan, ses arrondis, la pente, et le pivot
#: dans la seule direction qui porte
MARTEAU_UTILE = ("longueur_plan", "decalage_plan", "rayon_coin",
                 "angle_plan", "pivot_t")


def probleme(p, n_scan=240, face_libre=True, galet=False, ressort=None,
             plat_mini=0.0, gabarit=(0.0, 0.0), marteau=()):
    """Renvoie (evaluer_x, bornes, reparer_x) pour le profil de coeur.

    `ressort` : (C0, k) a **figer**. C'est la deuxieme etape du travail :
    une fois le coeur initial defini et son ressort calcule, on garde ce
    ressort et on fait varier le coeur et le marteau pour ameliorer la
    version initiale. Toutes les solutions sont alors comparees a armement
    identique, ce qui est le seul moyen de savoir laquelle est meilleure.
    Si `ressort` vaut None, le ressort reste une sortie, recalculee pour
    chaque profil.

    Si `face_libre`, une variable supplementaire est ajoutee en queue : la
    **courbure du plat du marteau** kappa = 1/R. On optimise la courbure et
    non le rayon, parce que la face plane (R infini) est alors simplement
    kappa = 0, au milieu de l'intervalle, et que les deux signes sont traites
    symetriquement : kappa > 0 face convexe (bombee vers le coeur),
    kappa < 0 face concave (qui enveloppe le coeur).
    """
    mini = max(0.0, getattr(p, "rayon_raccord", 0.0))
    # Le plat du zero est porte **explicitement** : il n'est pas une variable
    # d'optimisation mais un gabarit, repris de la synthese. Toutes les
    # solutions gardent donc le meme plat, exact, et le coeur sa position
    # nulle. Sans cela la serie n'en rendait que 0.31 mm sur 1.46.
    plat, beta0 = float(gabarit[0]), float(gabarit[1])

    # Grandeurs du marteau mises en variables. Le profil du coeur n'a presque
    # pas de prise sur la plage morte ni sur l'effort de contact : les leviers
    # sont la **longueur du pan**, ses **arrondis**, sa **pente** et la
    # **position du pivot**. Les laisser fixes revenait a demander a
    # l'optimisation d'ameliorer ce qu'elle ne pouvait pas toucher.
    marteau = tuple(n for n in marteau if n in MARTEAU_LIBRE)
    nh = len(marteau)
    # queue du vecteur : cx, cy [, kappa] [, grandeurs du marteau]
    nq = (3 if face_libre else 2) + nh
    i_kappa = -nh - 1 if face_libre else None

    def mu_du_contact(R):
        """Un galet monte fou remplace le glissement par du roulement :
        il ne reste que le frottement de son axe, ramene par le rapport des
        rayons. C'est le seul levier sur la plage morte, qui vaut
        2.atan(mu_contact)."""
        if not galet or abs(R) < 1e-9:
            return None
        return mu_galet(R, p.axe_galet, p.mu_palier)

    P0 = (p.pivot_x, p.pivot_y)
    d0 = norm(P0) or 1.0
    UT = (-P0[1] / d0, P0[0] / d0)        # perpendiculaire a O->P

    def params_de(x):
        """Parametres de cet individu : ceux du probleme, marteau ajuste."""
        if not nh:
            return p, beta0
        vals = {n: float(v) for n, v in zip(marteau, x[len(x) - nh:])}
        t = vals.pop("pivot_t", None)
        if t is not None:
            # deplacement perpendiculaire a O->P, a partir du pivot nominal
            vals.setdefault("pivot_x", P0[0] + t * UT[0])
            vals.setdefault("pivot_y", P0[1] + t * UT[1])
        q = p.copy(**vals)
        # gamma0 bouge avec la pente et le pivot : le plat du coeur doit
        # suivre, sinon le coeur ne part plus de sa position de repos
        b0 = (q.gamma0
              if any(n in ("angle_plan", "pivot_t", "pivot_x", "pivot_y")
                     for n in marteau) else beta0)
        return q, b0

    def reparer_x(x):
        """Projette le vecteur dans l'espace admissible.

        `profil.reparer` attend le vecteur **complet**, cx et cy compris : il
        les lit en queue et en deduit le nombre d'harmoniques. Lui passer
        `x[:-nq]` les retirait une premiere fois, puis il en retirait deux
        autres : le nombre d'harmoniques tombait de 11 a 10 et les series an
        et bn se trouvaient decalees d'un rang. La graine analytique en
        ressortait brouillee - irrecevable et plus mauvaise sur les quatre
        objectifs - et c'etait elle que l'optimisation finissait par rendre.
        """
        n_prof = len(x) - nh - (1 if face_libre else 0)
        tete = list(x[:n_prof])
        queue = list(x[n_prof:])
        _q, b0 = params_de(x)
        return list(reparer(tete, mini, plat=plat, beta0=b0)) + queue

    def face(x):
        if not face_libre:
            return 0.0
        k = x[i_kappa]
        return 0.0 if abs(k) < 1.0 / (20.0 * max(p.r_max, 1e-9)) else 1.0 / k

    def graines_face(x0):
        """Points de depart semes sur les directions qui portent.

        La mutation seule trouve mal une bande utile etroite : on la seme.
        C'etait vrai de la courbure de face, ca l'est encore plus du
        **pivot**. Sa direction utile est la perpendiculaire a O->P : sur
        +/- 4 mm, l'ondulation y varie de 68.8 points, contre 1.6 point
        radialement. Sans semis, dix generations plafonnaient a 35 %
        d'ondulation la ou un balayage a la main en trouvait 12.5.
        """
        base = [list(x0)]
        for nom, valeurs in (("pivot_t", (-0.25, -0.15, -0.08, 0.08, 0.15)),
                             ("angle_plan", (-12.0, -6.0, 6.0, 12.0))):
            if nom not in marteau:
                continue
            j = len(x0) - nh + marteau.index(nom)
            ech = d0 if nom == "pivot_t" else 1.0
            nouveaux = []
            for g in base:
                for v in valeurs:
                    h = list(g)
                    h[j] = (v * ech if nom == "pivot_t" else h[j] + v)
                    nouveaux.append(h)
            base += nouveaux
        if not face_libre:
            return base
        i = i_kappa
        d = norm(p.pivot) or 10.0
        out = []
        choix = ((0.5, 0.7, 1.0, 1.4, 2.0) if galet else
                 (0.0, -d, -0.75 * d, -0.6 * d, -1.3 * d, 0.5, 1.0, 2.0))
        for g0 in base:
            for f_ in choix:
                R = f_ * p.r_max if galet or f_ > 0 else f_
                g = list(g0)
                g[i] = 0.0 if R == 0.0 else 1.0 / R
                out.append(g)
        return out

    def ev(x):
        q, b0 = params_de(x)
        pr = Profil.from_vector(list(x[:-nq]) + [x[-nq], x[-nq + 1]],
                                plat, b0)
        R = face(x)
        mu_c = mu_du_contact(R)
        if ressort is None:
            res = dimensionner_ressort(q, pr, R, mu_c, n=n_scan)
            e = res.eval
        else:
            res = None
            e = evaluer(q.copy(couple_min=ressort[0], pente=ressort[1]), pr,
                        n=n_scan, rayon_face=R, mu_contact=mu_c)
        if e is None:
            return (1e3, 1e3, 1e3, 1e3), 1e3
        viol = e.penalite
        if not e.ok:
            return (1e3, 1e3, 1e3, 1e3), viol + 1e3
        # Tolerance de fabrication sur l'encombrement. Sans elle, la graine
        # analytique etait **irrecevable** pour 0.9 um de depassement
        # (3.0009 contre 3.0000 mm), simple residu de la reconstruction a
        # bande limitee. En domination contrainte, toute solution recevable
        # battait alors le meilleur coeur, aussi mauvaise soit-elle : c'est ce
        # qui faisait rendre a l'optimisation des profils plus mauvais que la
        # synthese dont elle part.
        if e.r_max > q.r_max + TOL_ENCOMBREMENT:
            viol += 10.0 * (e.r_max - q.r_max) / q.r_max
        if e.r_min < q.r_min - TOL_ENCOMBREMENT:
            viol += 10.0 * (q.r_min - e.r_min) / max(q.r_min, 1e-9)
        if e.couple_utile <= 1e-9 or e.couple_moyen <= 1e-9:
            return (1e3, 1e3, 1e3, 1e3), viol + 1e3
        if res is not None and (res.raideur < 0.0 or res.prearmage <= 0.0):
            viol += 5.0
        # raccords usinables : la contrainte porte sur la courbure **reelle**
        # de la courbe reconstruite, pas sur la valeur de cellule
        if e.rho_min < mini - 1e-9:
            viol += 10.0 * (mini - e.rho_min) / max(mini, 1e-3)
        if e.interference:
            viol += 0.05 * e.interference
        if e.interference:
            viol += 0.05 * e.interference
        # coherence mecanique (voir profil._controles)
        viol += 2.0 * abs(e.runs_couple - 2)
        viol += 0.5 * e.cote_change
        viol += 0.2 * e.saut_contact
        # Le plat du zero est ce qui donne la **position nulle** : sans lui le
        # coeur n'a plus de repos defini. Rien ne l'imposait, et l'optimisation
        # le dissolvait. On exige donc qu'il ne descende pas sous celui de la
        # graine - la reference est passee par l'appelant, car une serie a
        # bande limitee n'en rend qu'une fraction (voir `Profil`).
        if plat_mini > 1e-9:
            lp = longueur_plat(pr, b0, TOL_PLAT, grille=300)
            if lp < plat_mini:
                viol += 5.0 * (plat_mini - lp) / plat_mini
        # Pas de contrainte sur la fenetre de pointe ici. Elle servait de
        # proxy a la plage morte, qui est **deja** un objectif (le second), et
        # elle est non monotone en longueur de plan : sur le coeur synthetise
        # l'imposer menait a un plan de 1.05 mm qui lache le coeur sur 62 deg,
        # la ou 3.0 mm donne la meme plage morte avec 2 deg sans contact. Le
        # profil optimise est de toute facon convexe (rho >= rayon d'outil),
        # donc sans marche a enjamber.
        # le couple marteau est fige : le couple obtenu sur le coeur est une
        # sortie, et on mesure son ecart a la cible
        # part du tour sous la consigne : c'est le critere qui refuse
        # directement les longues plages a couple faible ou nul
        return (secteur_faible(e, q.couple_coeur), e.plage_morte,
                100.0 * e.ondulation, e.N_max), viol

    def bornes(x0):
        k = len(x0) - nq
        rmax = 5.0 * p.r_max
        lo = [mini] + [-rmax] * (k - 1) + [-p.r_max, -p.r_max]
        hi = [rmax] * k + [p.r_max, p.r_max]
        if face_libre:
            km = 1.0 / R_FACE_MINI
            # un galet est forcement convexe
            lo.append(1.0 / (6.0 * p.r_max) if galet else -km)
            hi.append(km)
        # Bornes du marteau : une fenetre autour de la valeur de depart, pour
        # rester dans un dimensionnement voisin plutot que d'explorer des
        # geometries ou la synthese n'a plus de sens.
        d_piv = 0.25 * max(norm(p.pivot), 1e-6)
        for n in marteau:
            v = 0.0 if n == "pivot_t" else float(getattr(p, n))
            if n == "longueur_plan":
                a, b = 0.15 * p.r_max, 1.4 * p.r_max
            elif n == "decalage_plan":
                a, b = -0.8 * p.r_max, 0.8 * p.r_max
            elif n == "rayon_coin":
                a, b = 0.02, max(0.06, 0.6 * p.r_min)
            elif n == "angle_plan":
                sg = 1.0 if v >= 0.0 else -1.0
                a, b = sorted((sg * 12.0, v + sg * 28.0))
                a, b = max(-88.0, a), min(88.0, b)
            elif n == "pivot_t":
                a, b = -0.3 * d0, 0.3 * d0
            else:                           # pivot_x, pivot_y
                a, b = v - d_piv, v + d_piv
            lo.append(a)
            hi.append(b)
        return lo, hi

    return ev, bornes, reparer_x, graines_face


# ------------------------------------------------- evaluation multi-coeurs

_W = {}


def _init_worker(p, n_scan, face_libre, galet=False, ressort=None,
                 plat_mini=0.0, gabarit=(0.0, 0.0), marteau=()):
    """Reconstruit le probleme dans le processus fils.

    `plat_mini` et `gabarit` doivent etre transmis : sans eux les fils
    evaluaient sans le plat du zero, donc un autre coeur que celui du
    processus principal. Toute la generation sortait irrecevable - « 0
    faisables » - et le front ne contenait plus la graine.
    """
    ev, _b, rep, _g = probleme(p, n_scan, face_libre, galet, ressort,
                               plat_mini, gabarit, marteau)
    _W["ev"], _W["rep"] = ev, rep


def _eval_worker(x):
    x = _W["rep"](x)
    f, v = _W["ev"](x)
    return x, f, v


def pool_lot(p, n_scan, face_libre, procs=None, galet=False, ressort=None,
             plat_mini=0.0, gabarit=(0.0, 0.0), marteau=()):
    """Renvoie (evaluer_lot, fermer) pour evaluer une generation en parallele.

    Le noyau compile est mono-thread : un seul coeur travaille, ce qui sur une
    machine a 16 coeurs plafonne a 7 % d'occupation processeur. Les individus
    d'une generation etant independants, c'est la qu'il faut paralleliser.
    """
    import multiprocessing as mp
    procs = procs or max(1, (mp.cpu_count() or 2) - 1)
    if procs <= 1:
        return None, (lambda: None)
    ctx = mp.get_context("fork")
    pool = ctx.Pool(procs, initializer=_init_worker,
                    initargs=(p, n_scan, face_libre, galet, ressort,
                              plat_mini, gabarit, marteau))

    def evaluer_lot(xs):
        return pool.map(_eval_worker, xs, chunksize=max(1, len(xs) // (4 * procs)))

    def fermer():
        pool.close()
        pool.join()

    return evaluer_lot, fermer
