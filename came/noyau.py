"""Noyau numerique accelere (numpy + numba), avec repli en Python pur.

Le calcul est domine par une seule boucle : pour chaque position du coeur,
chercher la tangence de la face du marteau puis en deduire couple et effort.
C'est ce que ce module compile.

`DISPONIBLE` dit si l'acceleration est utilisable. Sans numpy ni numba le
reste du programme fonctionne a l'identique, en Python pur : les deux chemins
sont valides l'un contre l'autre (voir `verifier`).

Installation de l'acceleration :

    python3 -m venv .venv && .venv/bin/pip install numba numpy
"""

import math

TAU = 2.0 * math.pi

try:
    import numpy as _np
    NUMPY = True
except ImportError:                                  # pragma: no cover
    _np = None
    NUMPY = False

try:
    from numba import njit as _njit
    NUMBA = True
except ImportError:                                  # pragma: no cover
    NUMBA = False

    def _njit(*a, **k):
        def deco(f):
            return f
        return deco if not a else a[0]

DISPONIBLE = NUMPY


# ------------------------------------------------------------------ outils

@_njit(cache=True, fastmath=True)
def _interp(bet0, dbet, Hs, ss, beta):
    """Interpolation lineaire 2.pi-periodique de (H, s)."""
    nb = Hs.shape[0] - 1
    x = (beta - bet0) % TAU
    j = int(x / dbet)
    if j >= nb:
        j = nb - 1
    w = (x - j * dbet) / dbet
    return (Hs[j] + w * (Hs[j + 1] - Hs[j]),
            ss[j] + w * (ss[j + 1] - ss[j]))


@_njit(cache=True, fastmath=True)
def _couple(Cm, sv, H, q, e, mu):
    """Couple et effort pour les deux sens de glissement possibles.

    Renvoie (T, N, bloque) : si aucun des deux sens ne donne un couple
    moteur, le coeur est bloque.
    """
    dA = q + mu * e
    dB = q - mu * e
    TA = Cm * (sv + mu * H) / dA if abs(dA) > 1e-12 else 0.0
    NA = abs(Cm / dA) if abs(dA) > 1e-12 else 0.0
    TB = Cm * (sv - mu * H) / dB if abs(dB) > 1e-12 else 0.0
    NB = abs(Cm / dB) if abs(dB) > 1e-12 else 0.0
    goA = abs(dA) > 1e-12 and TA < 0.0
    goB = abs(dB) > 1e-12 and TB > 0.0
    if goA and not goB:
        return TA, NA, 0
    if goB and not goA:
        return TB, NB, 0
    if goA and goB:
        if abs(TA) >= abs(TB):
            return TA, NA, 0
        return TB, NB, 0
    return 0.0, (NA if NA > 0.0 else NB), 1


@_njit(cache=True, fastmath=True)
def echantillonner(a0, an, bn, cx, cy, nb, plat, beta0):
    """H, s et courbure du profil sur une grille de normales.

    `plat` est la longueur du **plat du zero**, porte explicitement : un Dirac
    de courbure de poids L a la normale beta0. La fonction d'appui qui le
    porte, periodique et continue, est  L.sin(x/2)  avec
    x = (beta - beta0) mod 2pi ; `s` y saute de L, ce qui est exactement un
    segment droit de longueur L, et la courbure lissee qu'elle ajoute vaut
    (3L/4).sin(x/2), positive donc convexe. Le Dirac lui-meme n'est pas dans
    `rh` : un plat est usinable, c'est une droite.
    """
    Hs = _np.zeros(nb + 1)
    ss = _np.zeros(nb + 1)
    rh = _np.zeros(nb + 1)
    for i in range(nb + 1):
        b = TAU * i / nb
        h = a0 + cx * math.cos(b) + cy * math.sin(b)
        d = -cx * math.sin(b) + cy * math.cos(b)
        r = a0
        for j in range(an.shape[0]):
            nn = j + 2
            w = 1.0 / (1.0 - nn * nn)
            c = math.cos(nn * b)
            sn = math.sin(nn * b)
            h += (an[j] * c + bn[j] * sn) * w
            d += nn * (-an[j] * sn + bn[j] * c) * w
            r += an[j] * c + bn[j] * sn
        if plat > 0.0:
            x = (b - beta0) % TAU
            sx = math.sin(0.5 * x)
            h += plat * sx
            d += 0.5 * plat * math.cos(0.5 * x)
            r += 0.75 * plat * sx
        Hs[i] = h
        ss[i] = d
        rh[i] = r
    return Hs, ss, rh


@_njit(cache=True, fastmath=True)
def reparer_coeffs(a0, an, bn, mini, iters, grille, plancher):
    """Projections alternees : rho >= mini, et bande limitee sans rang 1.

    C'est la partie la plus couteuse de l'optimisation (appelee pour chaque
    individu a chaque generation), bien avant l'evaluation mecanique. Les
    tables de cosinus et sinus sont calculees une fois pour toutes : sans
    elles les appels trigonometriques dominent.

    `plancher` est le minimum **par echantillon** impose a la serie. Avec un
    plat porte explicitement, la contrainte d'usinabilite porte sur la
    courbure totale : la serie doit donc rester au-dessus de
    `mini - (3L/4).sin(x/2)`, qui depend de beta. Un plancher scalaire
    imposait le minimum a la seule serie et deformait le profil.
    """
    M = an.shape[0]
    C = _np.zeros((M, grille))
    S = _np.zeros((M, grille))
    for j in range(M):
        nn = j + 2
        for i in range(grille):
            x = TAU * i / grille
            C[j, i] = math.cos(nn * x)
            S[j, i] = math.sin(nn * x)
    vals = _np.zeros(grille)
    na = _np.zeros(M)
    nb_ = _np.zeros(M)
    for _ in range(iters):
        m = 1e18
        for i in range(grille):
            r = a0
            for j in range(M):
                r += an[j] * C[j, i] + bn[j] * S[j, i]
            vals[i] = r
            if r - plancher[i] < m:
                m = r - plancher[i]
        if m >= -1e-9:
            break
        for i in range(grille):
            if vals[i] < plancher[i]:
                vals[i] = plancher[i]
        s0 = 0.0
        for i in range(grille):
            s0 += vals[i]
        a0 = s0 / grille
        for j in range(M):
            c = 0.0
            sn = 0.0
            for i in range(grille):
                c += vals[i] * C[j, i]
                sn += vals[i] * S[j, i]
            na[j] = 2.0 * c / grille
            nb_[j] = 2.0 * sn / grille
        for j in range(M):
            an[j] = na[j]
            bn[j] = nb_[j]
    m = 1e18
    for i in range(grille):
        r = a0
        for j in range(M):
            r += an[j] * C[j, i] + bn[j] * S[j, i]
        if r - plancher[i] < m:
            m = r - plancher[i]
    if m < 0.0:
        a0 += (-m) + 1e-6
    return a0, an, bn


# --------------------------------------------------------------- face plane

@_njit(cache=True, fastmath=True)
def scan_plan(Hs, ss, bet0, dbet, px, py, g0, h0, eps, mu, C0, kk,
              psi_max, n):
    """Un tour complet avec une face de marteau plane."""
    th_out = _np.zeros(n)
    be_out = _np.zeros(n)
    H_out = _np.zeros(n)
    s_out = _np.zeros(n)
    ph_out = _np.zeros(n)
    T_out = _np.zeros(n)
    N_out = _np.zeros(n)
    bl_out = _np.zeros(n, dtype=_np.int64)
    valide = _np.zeros(n, dtype=_np.int64)
    for i in range(n):
        th = TAU * i / n
        lo = 0.0
        hi = psi_max
        cg = math.cos(g0)
        sg = math.sin(g0)
        f0 = cg * px + sg * py + h0 - _interp(bet0, dbet, Hs, ss, g0 - th)[0]
        if f0 >= 0.0:
            psi = 0.0
        else:
            gh = g0 + eps * psi_max
            fh = (math.cos(gh) * px + math.sin(gh) * py + h0
                  - _interp(bet0, dbet, Hs, ss, gh - th)[0])
            if fh < 0.0:
                continue
            for _ in range(50):
                m = 0.5 * (lo + hi)
                gm = g0 + eps * m
                fm = (math.cos(gm) * px + math.sin(gm) * py + h0
                      - _interp(bet0, dbet, Hs, ss, gm - th)[0])
                if fm < 0.0:
                    lo = m
                else:
                    hi = m
            psi = 0.5 * (lo + hi)
        phi = eps * psi
        gamma = g0 + phi
        cgm = math.cos(gamma)
        sgm = math.sin(gamma)
        H, sv = _interp(bet0, dbet, Hs, ss, gamma - th)
        q = sv - (-sgm * px + cgm * py)
        Cm = -eps * (C0 + kk * psi)
        T, N, bl = _couple(Cm, sv, H, q, h0, mu)
        th_out[i] = math.degrees(th)
        be_out[i] = gamma - th
        H_out[i] = H
        s_out[i] = sv
        ph_out[i] = phi
        T_out[i] = T
        N_out[i] = N
        bl_out[i] = bl
        valide[i] = 1
    return th_out, be_out, H_out, s_out, ph_out, T_out, N_out, bl_out, valide


# -------------------------------------------------------------- face courbe

@_njit(cache=True, fastmath=True)
def _serie(a0, an, bn, cx, cy, beta, deriv):
    """H(beta) ou s(beta) a partir des coefficients de courbure."""
    if deriv == 0:
        v = a0 + cx * math.cos(beta) + cy * math.sin(beta)
    else:
        v = -cx * math.sin(beta) + cy * math.cos(beta)
    for i in range(an.shape[0]):
        nn = i + 2
        w = 1.0 / (1.0 - nn * nn)
        if deriv == 0:
            v += (an[i] * math.cos(nn * beta) + bn[i] * math.sin(nn * beta)) * w
        else:
            v += nn * (-an[i] * math.sin(nn * beta)
                       + bn[i] * math.cos(nn * beta)) * w
    return v


@_njit(cache=True, fastmath=True)
def scan_courbe(Hs, ss, bet0, dbet, a0, an, bn, cx, cy, px, py, R, L,
                ang_q0, th0, eps, mu, C0, kk, psi_max, n, pas, prof_x,
                prof_y, tol):
    """Un tour complet avec une face de marteau circulaire (R signe)."""
    nb = Hs.shape[0] - 1
    Yx = _np.zeros(nb + 1)
    Yy = _np.zeros(nb + 1)
    for i in range(nb + 1):
        b = bet0 + i * dbet
        cb = math.cos(b)
        sb = math.sin(b)
        Yx[i] = (Hs[i] + R) * cb - ss[i] * sb
        Yy[i] = (Hs[i] + R) * sb + ss[i] * cb

    th_out = _np.zeros(n)
    be_out = _np.zeros(n)
    H_out = _np.zeros(n)
    s_out = _np.zeros(n)
    ph_out = _np.zeros(n)
    T_out = _np.zeros(n)
    N_out = _np.zeros(n)
    bl_out = _np.zeros(n, dtype=_np.int64)
    valide = _np.zeros(n, dtype=_np.int64)
    n_int = 0
    for i in range(n):
        th = TAU * i / n + th0
        ct = math.cos(-th)
        st = math.sin(-th)
        ptx = px * ct - py * st
        pty = px * st + py * ct
        best_phi = 0.0
        best_b = 0.0
        trouve = 0
        j = 0
        fj = (Yx[0] - ptx) ** 2 + (Yy[0] - pty) ** 2 - L * L
        while j + pas <= nb:
            jn = j + pas
            fn = (Yx[jn] - ptx) ** 2 + (Yy[jn] - pty) ** 2 - L * L
            if fj == 0.0 or fj * fn < 0.0:
                w = 0.0 if fn == fj else fj / (fj - fn)
                b = bet0 + (j + w * pas) * dbet
                # raffinement sur les valeurs exactes (secante)
                x0 = b - 0.25 * pas * dbet
                x1 = b + 0.25 * pas * dbet
                g0_ = _gserie(a0, an, bn, cx, cy, R, ptx, pty, L, x0)
                g1_ = _gserie(a0, an, bn, cx, cy, R, ptx, pty, L, x1)
                for _ in range(4):
                    if abs(g1_ - g0_) < 1e-18:
                        break
                    x2 = x1 - g1_ * (x1 - x0) / (g1_ - g0_)
                    if abs(x2 - x1) < 1e-12:
                        x1 = x2
                        break
                    x0 = x1
                    g0_ = g1_
                    x1 = x2
                    g1_ = _gserie(a0, an, bn, cx, cy, R, ptx, pty, L, x1)
                if abs(x1 - b) < 2.0 * pas * dbet:
                    b = x1
                H = _serie(a0, an, bn, cx, cy, b, 0)
                sv = _serie(a0, an, bn, cx, cy, b, 1)
                cb = math.cos(b)
                sb = math.sin(b)
                ybx = (H + R) * cb - sv * sb
                yby = (H + R) * sb + sv * cb
                cth = math.cos(th)
                sth = math.sin(th)
                qfx = ybx * cth - yby * sth
                qfy = ybx * sth + yby * cth
                phi = math.atan2(qfy - py, qfx - px) - ang_q0
                phi = (phi + math.pi) % TAU - math.pi
                if abs(phi) <= psi_max:
                    # plus petite ouverture sans interpenetration : le
                    # ressort arrete le marteau au premier contact
                    ok = 1
                    if R > 0.0:
                        dm = 1e18
                        for t2 in range(prof_x.shape[0]):
                            d2 = ((ybx - prof_x[t2]) ** 2
                                  + (yby - prof_y[t2]) ** 2)
                            if d2 < dm:
                                dm = d2
                        if dm < (R - tol) ** 2:
                            ok = 0
                    else:
                        dm = 0.0
                        for t2 in range(prof_x.shape[0]):
                            d2 = ((ybx - prof_x[t2]) ** 2
                                  + (yby - prof_y[t2]) ** 2)
                            if d2 > dm:
                                dm = d2
                        if dm > (abs(R) + tol) ** 2:
                            ok = 0
                    if ok == 1 and (trouve == 0 or abs(phi) < abs(best_phi)):
                        best_phi = phi
                        best_b = b
                        trouve = 1
            j = jn
            fj = fn
        if trouve == 0:
            n_int += 1
            continue
        b = best_b
        phi = best_phi
        H = _serie(a0, an, bn, cx, cy, b, 0)
        sv = _serie(a0, an, bn, cx, cy, b, 1)
        gamma = b + th
        cgm = math.cos(gamma)
        sgm = math.sin(gamma)
        q = sv - (-sgm * px + cgm * py)
        ee = H - (cgm * px + sgm * py)
        Cm = -eps * (C0 + kk * abs(phi))
        T, N, bl = _couple(Cm, sv, H, q, ee, mu)
        th_out[i] = math.degrees(TAU * i / n)
        be_out[i] = b
        H_out[i] = H
        s_out[i] = sv
        ph_out[i] = phi
        T_out[i] = T
        N_out[i] = N
        bl_out[i] = bl
        valide[i] = 1
    return (th_out, be_out, H_out, s_out, ph_out, T_out, N_out, bl_out,
            valide, n_int)


@_njit(cache=True, fastmath=True)
def _gserie(a0, an, bn, cx, cy, R, ptx, pty, L, x):
    H = _serie(a0, an, bn, cx, cy, x, 0)
    sv = _serie(a0, an, bn, cx, cy, x, 1)
    cb = math.cos(x)
    sb = math.sin(x)
    yx = (H + R) * cb - sv * sb
    yy = (H + R) * sb + sv * cb
    return (yx - ptx) ** 2 + (yy - pty) ** 2 - L * L


# ------------------------------------------- plan de marteau a longueur finie

@_njit(cache=True, fastmath=True)
def scan_plan_fini(Hs, ss, bet0, dbet, px, py, g0, h0, eps, mu, C0, kk,
                   psi_max, n, demi, rho, pas, dec):
    """Un tour avec un plan de longueur finie et des coins arrondis.

    Trois contacts possibles a chaque position : le plan, et chacun de ses
    deux coins. Le marteau se ferme jusqu'au premier obstacle, donc celui qui
    exige la plus grande ouverture s'impose. Sur un coin la normale est celle
    du **coeur** au point touche, pas celle du plan.

    `pas` n'est plus utilise : la separation des coins est resolue en forme
    close, appui par appui, au lieu d'etre cherchee par racine sur la courbe
    offset du coeur.
    """
    nb = Hs.shape[0] - 1
    th_out = _np.zeros(n)
    be_out = _np.zeros(n)
    H_out = _np.zeros(n)
    s_out = _np.zeros(n)
    ph_out = _np.zeros(n)
    T_out = _np.zeros(n)
    N_out = _np.zeros(n)
    bl_out = _np.zeros(n, dtype=_np.int64)
    valide = _np.zeros(n, dtype=_np.int64)
    coin_out = _np.zeros(n, dtype=_np.int64)
    n_coin = 0

    cg0 = math.cos(g0)
    sg0 = math.sin(g0)
    H0 = Hs[0]
    for i in range(nb):
        if Hs[i] < H0:
            H0 = Hs[i]
    m0x = H0 * cg0
    m0y = H0 * sg0
    # les deux coins, dans le repere du marteau au repos
    lc = _np.zeros(2)
    ac = _np.zeros(2)
    # Le centre est rentre de rho dans la matiere (+u0) : l'arrondi est
    # tangent au plan et ne depasse pas devant lui. Place sur le plan, il
    # avancait de rho vers le coeur et le coin touchait a chaque position.
    for c in range(2):
        sgc = -1.0 if c == 0 else 1.0
        cx = m0x + (dec + sgc * demi) * (-sg0) + rho * cg0
        cy = m0y + (dec + sgc * demi) * cg0 + rho * sg0
        lc[c] = math.sqrt((cx - px) ** 2 + (cy - py) ** 2)
        ac[c] = math.atan2(cy - py, cx - px)

    # tables independantes de la position du coeur
    cbt = _np.zeros(nb)
    sbt = _np.zeros(nb)
    hrt = _np.zeros(nb)
    for j in range(nb):
        bj = bet0 + j * dbet
        cbt[j] = math.cos(bj)
        sbt[j] = math.sin(bj)
        hrt[j] = Hs[j] + rho

    for i in range(n):
        th = TAU * i / n
        meilleur = -1.0
        mphi = 0.0
        mbeta = 0.0
        mH = 0.0
        ms = 0.0
        me = 0.0
        mcoin = 0

        # Separation exacte du nez du marteau.
        #
        # Le nez est le segment des deux centres de coin, dilate de rho. Le
        # coeur est convexe, donc intersection des demi-plans u_j.x <= H_j :
        # un segment en est hors **si et seulement si** un appui met ses deux
        # extremites au-dela,  min(u_j.A, u_j.B) >= H_j + rho.  L'ensemble
        # admissible est l'union sur j de ces conditions, et l'ouverture
        # cherchee en est le minimum.
        #
        # Traiter le plan et les deux coins separement puis prendre la plus
        # grande des trois ouvertures etait faux : l'ensemble admissible d'un
        # coin est un **arc**, pas un intervalle croissant. Le coin cote -t
        # etait degage a psi = 0, rentrait dans le coeur a psi = 4.6, et la
        # regle du maximum retenait 4.6 : penetration de 0.34 mm et saut de
        # 26 deg sur le couple.
        ct = math.cos(-th)
        st = math.sin(-th)
        ptx = px * ct - py * st
        pty = px * st + py * ct
        aA = ac[0] - th
        aB = ac[1] - th
        for j in range(nb):
            cj = cbt[j] * ptx + sbt[j] * pty
            num = hrt[j] - cj
            cA = num / lc[0]
            cB = num / lc[1]
            if cA > 1.0 or cB > 1.0:
                continue
            bj = bet0 + j * dbet
            # arc admissible de chaque extremite, en psi
            okA = cA < -1.0
            loA = 0.0
            hiA = 0.0
            if not okA:
                w = math.acos(cA)
                loA = bj - aA - w
                hiA = bj - aA + w
                if eps < 0.0:
                    tmp = loA
                    loA = -hiA
                    hiA = -tmp
            okB = cB < -1.0
            loB = 0.0
            hiB = 0.0
            if not okB:
                w = math.acos(cB)
                loB = bj - aB - w
                hiB = bj - aB + w
                if eps < 0.0:
                    tmp = loB
                    loB = -hiB
                    hiB = -tmp
            for ic in range(3):
                if ic == 0:
                    psi = 0.0
                elif ic == 1:
                    if okA:
                        continue
                    psi = loA % TAU
                else:
                    if okB:
                        continue
                    psi = loB % TAU
                if psi > psi_max:
                    continue
                if not okA and (psi - loA) % TAU > hiA - loA + 1e-12:
                    continue
                if not okB and (psi - loB) % TAU > hiB - loB + 1e-12:
                    continue
                if meilleur < -0.5 or psi < meilleur:
                    meilleur = psi
                    mbeta = bj
                    ahx = ptx + lc[0] * math.cos(ac[0] - th + eps * psi)
                    ahy = pty + lc[0] * math.sin(ac[0] - th + eps * psi)
                    bhx = ptx + lc[1] * math.cos(ac[1] - th + eps * psi)
                    bhy = pty + lc[1] * math.sin(ac[1] - th + eps * psi)
                    cbb = math.cos(bj)
                    sbb = math.sin(bj)
                    Hj, sj = _interp(bet0, dbet, Hs, ss, bj)
                    xjx = Hj * cbb - sj * sbb
                    xjy = Hj * sbb + sj * cbb
                    # Le point du coeur qui touche se projette-t-il **dans**
                    # le segment ? Si oui le plan porte et le bras de levier
                    # est celui de ce point ; sinon c'est un coin, et c'est
                    # le pied de la perpendiculaire issue de son centre.
                    dxs = bhx - ahx
                    dys = bhy - ahy
                    l2 = dxs * dxs + dys * dys
                    if l2 < 1e-18:
                        tp = 0.5
                    else:
                        tp = ((xjx - ahx) * dxs + (xjy - ahy) * dys) / l2
                    if 1e-9 < tp < 1.0 - 1e-9:
                        mcoin = 0
                        Hb = Hj
                        sb2 = sj
                    else:
                        mcoin = 1
                        if tp <= 0.5:
                            chx, chy = ahx, ahy
                        else:
                            chx, chy = bhx, bhy
                        Hb = cbb * chx + sbb * chy - rho
                        sb2 = -sbb * chx + cbb * chy
                    mphi = eps * psi
                    mH = Hb
                    ms = sb2
                    gam = bj + th
                    me = Hb - (math.cos(gam) * px + math.sin(gam) * py)
                break

        # meilleur est desormais l'ouverture retenue, qui peut valoir zero :
        # -1 est la valeur sentinelle d'absence de solution
        if meilleur < -0.5:
            continue
        gam = mbeta + th
        q = ms - (-math.sin(gam) * px + math.cos(gam) * py)
        Cm = -eps * (C0 + kk * abs(mphi))
        T, N, bl = _couple(Cm, ms, mH, q, me, mu)
        th_out[i] = math.degrees(th)
        be_out[i] = mbeta
        H_out[i] = mH
        s_out[i] = ms
        ph_out[i] = mphi
        T_out[i] = T
        N_out[i] = N
        bl_out[i] = bl
        valide[i] = 1
        coin_out[i] = mcoin
        n_coin += mcoin
    return (th_out, be_out, H_out, s_out, ph_out, T_out, N_out, bl_out,
            valide, n_coin, coin_out)
