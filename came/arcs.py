"""Approximation du profil par une suite d'arcs de cercle tangents.

Un dessin d'atelier ne se cote pas sur une polyligne de 720 points : il se
cote sur une suite d'arcs, chacun par son centre et son rayon. L'enjeu est de
savoir combien il en faut pour rester sous un ecart donne.

Le calcul se fait dans la **fonction d'appui**, ce qui rend la tangence
automatique. Un cercle de centre C et de rayon R a pour fonction d'appui

    H(beta) = C.u(beta) + R

Deux arcs voisins qui se raccordent a la normale `beta_k` partagent ce point
de contact si

    C_k = C_{k-1} + (R_{k-1} - R_k) u(beta_k)

et ils y partagent alors la **meme normale**, donc la meme tangente : le
raccordement est G1 par construction, sans biarc ni point de jonction libre.
Toute la chaine se deduit de `C_0` et des rayons.

L'ecart entre le profil et son approximation se lit directement :
`H(beta) - H_arcs(beta)` est la distance entre les deux tangentes de meme
normale, donc l'ecart mesure le long de la normale. C'est l'erreur que l'on
borne.

La fermeture du contour impose deux equations lineaires sur les rayons :

    somme_k (R_{k-1} - R_k) u(beta_k) = 0

Elles sont posees comme contraintes de l'ajustement, et non corrigees apres
coup : le contour obtenu se referme.
"""

import math

try:
    import numpy as _np
except ImportError:                        # pragma: no cover
    _np = None

TAU = 2.0 * math.pi


class Arc:
    """Un arc : centre, rayon, et la plage de **normales** qu'il couvre."""

    __slots__ = ("cx", "cy", "r", "b0", "b1")

    def __init__(self, cx, cy, r, b0, b1):
        self.cx, self.cy, self.r = cx, cy, r
        self.b0, self.b1 = b0, b1

    @property
    def centre(self):
        return (self.cx, self.cy)

    def point(self, beta):
        """Point de contact pour la normale `beta`."""
        return (self.cx + self.r * math.cos(beta),
                self.cy + self.r * math.sin(beta))

    @property
    def depart(self):
        return self.point(self.b0)

    @property
    def arrivee(self):
        return self.point(self.b1)

    def echantillonner(self, n=24):
        n = max(2, int(n))
        return [self.point(self.b0 + (self.b1 - self.b0) * i / (n - 1))
                for i in range(n)]

    def pour_dxf(self):
        """(centre, rayon positif, angle debut, angle fin) en degres.

        Un rayon negatif decrit un arc **concave** - les raccords entre les
        flancs et le plat du zero en sont. Le point de contact vaut alors
        `C + |R| u(beta + pi)` : on decale donc les deux angles de 180 deg et
        on rend le rayon positif, ce qu'attend une entite ARC.
        """
        if self.r >= 0.0:
            return (self.cx, self.cy), self.r, \
                math.degrees(self.b0), math.degrees(self.b1)
        return (self.cx, self.cy), -self.r, \
            math.degrees(self.b0 + math.pi), math.degrees(self.b1 + math.pi)

    @property
    def concave(self):
        return self.r < 0.0

    def __repr__(self):                    # pragma: no cover
        return ("Arc(centre=(%.4f, %.4f), r=%.4f, %.2f..%.2f deg)"
                % (self.cx, self.cy, self.r,
                   math.degrees(self.b0), math.degrees(self.b1)))


class Approximation:
    """Resultat : la chaine d'arcs, l'ecart obtenu et la fermeture."""

    def __init__(self, arcs, ecart, fermeture, betas):
        self.arcs = arcs
        self.ecart = ecart              # ecart maxi le long de la normale [mm]
        self.fermeture = fermeture      # residu de fermeture du contour [mm]
        self.betas = betas              # normales de raccordement

    @property
    def n(self):
        return len(self.arcs)

    def contour(self, par_arc=24):
        out = []
        for a in self.arcs:
            out.extend(a.echantillonner(par_arc)[:-1])
        return out

    def rayons(self):
        return [a.r for a in self.arcs]


def _ajuster(Hs, bet, coupes):
    """Ajuste la chaine d'arcs pour un decoupage donne.

    Inconnues : `C_0` (deux composantes) et les rayons. Les centres suivants
    se deduisent par la relation de raccordement, donc `H_arcs(beta)` est
    **lineaire** en ces inconnues : l'ajustement est un moindres carres, avec
    les deux equations de fermeture ajoutees en contrainte ponderee.
    """
    N = len(coupes)
    nb = len(Hs)
    # u(beta_k) aux coupes
    uk = [(math.cos(bet[i]), math.sin(bet[i])) for i in coupes]

    # C_k = C_0 + somme_{m=1..k} (R_{m-1} - R_m) u_m
    # -> coefficients de C_k sur (cx0, cy0, R_0..R_{N-1})
    coefC = []                              # coefC[k] = (ax, ay) listes
    ax = [0.0] * (N + 2)
    ay = [0.0] * (N + 2)
    ax[0] = 1.0
    ay[1] = 1.0
    coefC.append((list(ax), list(ay)))
    for k in range(1, N):
        ux, uy = uk[k]
        ax[2 + k - 1] += ux
        ax[2 + k] -= ux
        ay[2 + k - 1] += uy
        ay[2 + k] -= uy
        coefC.append((list(ax), list(ay)))

    lignes = []
    cibles = []
    for k in range(N):
        i0 = coupes[k]
        i1 = coupes[k + 1] if k + 1 < N else coupes[0] + nb
        axk, ayk = coefC[k]
        for i in range(i0, i1):
            b = bet[i % nb]
            cb, sb = math.cos(b), math.sin(b)
            L = [axk[j] * cb + ayk[j] * sb for j in range(N + 2)]
            L[2 + k] += 1.0                 # le rayon R_k
            lignes.append(L)
            cibles.append(Hs[i % nb])

    # fermeture : somme_k (R_{k-1} - R_k) u_k = 0, indices modulo N
    for comp in (0, 1):
        L = [0.0] * (N + 2)
        for k in range(N):
            u = uk[k][comp]
            L[2 + (k - 1) % N] += u
            L[2 + k] -= u
        lignes.append([v * 1.0e4 for v in L])
        cibles.append(0.0)

    A = _np.array(lignes, dtype=float)
    y = _np.array(cibles, dtype=float)
    x, *_ = _np.linalg.lstsq(A, y, rcond=None)

    # reconstruction des arcs
    R = [float(x[2 + k]) for k in range(N)]
    C = [(float(x[0]), float(x[1]))]
    for k in range(1, N):
        ux, uy = uk[k]
        d = R[k - 1] - R[k]
        C.append((C[-1][0] + d * ux, C[-1][1] + d * uy))
    arcs = []
    for k in range(N):
        b0 = bet[coupes[k]]
        b1 = (bet[coupes[k + 1]] if k + 1 < N else bet[coupes[0]] + TAU)
        arcs.append(Arc(C[k][0], C[k][1], R[k], b0, b1))

    # ecart le long de la normale, et residu de fermeture
    ecart = 0.0
    pires = [0.0] * N
    for k in range(N):
        i0 = coupes[k]
        i1 = coupes[k + 1] if k + 1 < N else coupes[0] + nb
        a = arcs[k]
        for i in range(i0, i1):
            b = bet[i % nb]
            h = a.cx * math.cos(b) + a.cy * math.sin(b) + a.r
            d = abs(Hs[i % nb] - h)
            if d > pires[k]:
                pires[k] = d
        ecart = max(ecart, pires[k])
    ux, uy = uk[0]
    dfx = C[0][0] - (C[-1][0] + (R[-1] - R[0]) * ux)
    dfy = C[0][1] - (C[-1][1] + (R[-1] - R[0]) * uy)
    return Approximation(arcs, ecart, math.hypot(dfx, dfy),
                         [bet[i] for i in coupes]), pires


def approcher(profil, ecart_maxi=0.005, n_grille=720, n_mini=4, n_maxi=80):
    """Approche le profil par des arcs tangents, sous un ecart donne.

    On part de `n_mini` arcs et on coupe a chaque tour l'arc le plus mauvais,
    jusqu'a passer sous `ecart_maxi` ou atteindre `n_maxi` arcs. L'utilisateur
    fixe l'ecart, le nombre d'arcs en decoule - c'est le sens utile : un
    nombre d'arcs impose ne dit rien de la precision obtenue.

    Renvoie une `Approximation`, ou `None` si numpy manque.
    """
    if _np is None:                        # pragma: no cover
        return None
    bet = [TAU * i / n_grille for i in range(n_grille)]
    Hs = [profil.H(b) for b in bet]
    coupes = [int(round(n_grille * k / n_mini)) % n_grille
              for k in range(n_mini)]
    best, pires = _ajuster(Hs, bet, coupes)
    while best.ecart > ecart_maxi and len(coupes) < n_maxi:
        k = max(range(len(coupes)), key=lambda j: pires[j])
        i0 = coupes[k]
        i1 = coupes[k + 1] if k + 1 < len(coupes) else coupes[0] + n_grille
        milieu = (i0 + i1) // 2 % n_grille
        if milieu in coupes or i1 - i0 < 4:
            # cet arc ne peut plus etre coupe : on coupe le suivant en taille
            ordre = sorted(range(len(coupes)), key=lambda j: -pires[j])
            fait = False
            for k2 in ordre[1:]:
                j0 = coupes[k2]
                j1 = (coupes[k2 + 1] if k2 + 1 < len(coupes)
                      else coupes[0] + n_grille)
                m2 = (j0 + j1) // 2 % n_grille
                if m2 not in coupes and j1 - j0 >= 4:
                    coupes = sorted(coupes + [m2])
                    fait = True
                    break
            if not fait:
                break
        else:
            coupes = sorted(coupes + [milieu])
        best, pires = _ajuster(Hs, bet, coupes)
    return best


def mesurer_coeur(pr, beta0, tol_arc=0.004):
    """Deduit du contour les grandeurs de construction du coeur.

    Repartir d'un dessin pour definir un nouveau coeur demande de retrouver,
    dans le contour, ce que la synthese prend en entree. Tout n'y est pas :
    le **frottement** et le **couple de rappel voulu** sont une donnee
    physique et une specification, aucun contour ne les porte. Le reste se
    mesure.

    On s'appuie sur l'approximation par arcs : les rayons qu'elle sort sont
    directement les rayons de la piece - un rayon negatif est un raccord
    concave, un petit rayon positif au plus grand encombrement est la pointe.

    Renvoie un dictionnaire des grandeurs trouvees, et la liste de celles qui
    manquent.
    """
    import math as _m
    out, manque = {}, []
    ap = approcher(pr, tol_arc)
    if ap is None or not ap.arcs:
        return {}, ["approximation par arcs indisponible"]

    pts = [pr.point(TAU * i / 720) for i in range(720)]
    ray = [_m.hypot(q[0], q[1]) for q in pts]
    out["r_max"] = max(ray)
    out["r_min"] = pr.H(beta0)             # distance du plat a l'axe

    # raccords concaves : le plus petit rayon en valeur absolue
    conc = [-a.r for a in ap.arcs if a.r < 0.0]
    if conc:
        out["rayon_raccord"] = min(conc)

    # pointe : parmi les arcs convexes, ceux dont le milieu est au plus grand
    # rayon ; le plus petit de leurs rayons est le rayon de pointe
    bout = []
    for a in ap.arcs:
        if a.r <= 0.0:
            continue
        q = a.point(0.5 * (a.b0 + a.b1))
        if _m.hypot(q[0], q[1]) > 0.9 * out["r_max"]:
            bout.append(a.r)
    if bout:
        out["rayon_pointe"] = min(bout)

    # corde de la pointe : etendue des points au-dela de 99 % du rayon maxi
    haut = [q for q, rr in zip(pts, ray) if rr > 0.99 * out["r_max"]]
    if len(haut) >= 2:
        corde = max(_m.hypot(a[0] - b[0], a[1] - b[1])
                    for a in haut for b in haut)
        if corde > 1e-6:
            out["corde_pointe"] = corde

    for nom in ("mu", "couple_coeur"):
        manque.append(nom)
    return out, manque
