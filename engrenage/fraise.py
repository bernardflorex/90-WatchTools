"""Fraise de taillage par génération : profil conjugué d'un mobile pour une roue de rayon
infini (crémaillère), et simulation du taillage d'un mobile par cette fraise.

Repris de `Abwfrac.F90`/`Abwfrap.F90` (roue -> fraise), `Radabwp.F90` (fraise -> roue taillée),
`Krvabc.F90`/`Krvabp.F90` (écart maximal) et `Spliph.F90` (K. Trümpy) : le flanc d'un mobile, en
prise avec une crémaillère (roulement sans glissement sur le cercle primitif), balaie une famille
de courbes dans le plan ; l'enveloppe de cette famille est le profil de la fraise, dans le repère
de la crémaillère (origine sur la droite primitive, pas au centre du mobile). Le mécanisme inverse
(le profil de la fraise balaie une famille de courbes dans le plan d'un mobile cible, quelconque)
donne le profil « pratique » réellement taillé, à comparer au profil de référence.

Simplification par rapport au Fortran : `Abwfrac.F90` calcule l'enveloppe roue -> fraise
analytiquement sur les arcs de cercle (intersection cercle/horizontale) ; ici les deux sens
utilisent la même méthode par intersection de segments, sur un échantillonnage fin du profil
(`arcs.echantillonner`), ce qui unifie le code et évite les cas particuliers des arcs (signe,
balayage de plus de 180°) au prix d'une erreur de discrétisation rendue négligeable par un
échantillonnage serré (quelques µm à 2000 points, très en dessous des tolérances de fabrication).
"""

import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

import arcs

N_ECHANT = 2000  # points de l'échantillonnage fin du flanc pour l'enveloppe roue -> fraise


N_PERIODES = 9  # périodes de denture balayées (NG de ABWFRAC, NW de RADABWP)


def periodes_defaut(zz):
    """Nombre de périodes de denture balayées par défaut : 9, comme le Fortran, quel que soit le
    nombre de dents.

    Ce nombre a longtemps été réduit pour les roues à peu de dents, parce qu'au-delà d'un
    balayage d'environ +-90 degrés le résultat devenait aberrant (2.3 mm d'écart sur une roue à
    12 dents). La cause n'était pas le balayage mais `_enveloppe`, qui acceptait les segments
    descendants et allait donc chercher son minimum sur une copie retournée ; une fois ce filtre
    remis (le test de ABWFRAC/RADABWP), 5, 7, 9 et 11 périodes donnent exactement le même profil,
    vérifié de 7 à 150 dents. On garde donc la valeur du Fortran.
    """
    return N_PERIODES


# ---------------------------------------------------------------------------------------
# spline paramétrique de Hermite (SPLIPH) : gère les courbes non monotones en x
# ---------------------------------------------------------------------------------------

def _moments(t, y, d0, dn):
    """Dérivées secondes (moments) d'une spline cubique par les points (t, y) ; `d0`/`dn` :
    dérivée première imposée à chaque extrémité (None : dérivée seconde nulle, spline naturelle)."""
    n = len(t)
    h = np.diff(t)
    a, b, c = np.zeros(n), np.zeros(n), np.zeros(n)  # sous-diag., diag., sur-diag.
    rhs = np.zeros(n)
    if d0 is None:
        b[0] = 1.0
    else:
        b[0], c[0] = 2.0 * h[0], h[0]
        rhs[0] = 6.0 * ((y[1] - y[0]) / h[0] - d0)
    for i in range(1, n - 1):
        a[i], b[i], c[i] = h[i - 1], 2.0 * (h[i - 1] + h[i]), h[i]
        rhs[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    if dn is None:
        b[n - 1] = 1.0
    else:
        a[n - 1], b[n - 1] = h[n - 2], 2.0 * h[n - 2]
        rhs[n - 1] = 6.0 * (dn - (y[n - 1] - y[n - 2]) / h[n - 2])
    # Thomas (tridiagonal)
    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        rhs[i] -= w * rhs[i - 1]
    m = np.zeros(n)
    m[-1] = rhs[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        m[i] = (rhs[i] - c[i] * m[i + 1]) / b[i]
    return m


def _evaluer(t, y, m, ts):
    h = np.diff(t)
    k = np.clip(np.searchsorted(t, ts, side="right") - 1, 0, len(t) - 2)
    hk = h[k]
    u, v = ts - t[k], t[k + 1] - ts
    y0, y1, m0, m1 = y[k], y[k + 1], m[k], m[k + 1]
    val = (m0 * v**3 + m1 * u**3) / (6 * hk) + (y0 / hk - m0 * hk / 6) * v + (y1 / hk - m1 * hk / 6) * u
    d1 = (3 * m1 * u**2 - 3 * m0 * v**2) / (6 * hk) + (y1 - y0) / hk - (m1 - m0) * hk / 6
    return val, d1


def spline_hermite(pts, angle_debut, angle_fin, n):
    """`n` points équidistants (abscisse curviligne) sur la spline de Hermite paramétrique passant
    par `pts`, de pente `tan(angle_debut)`/`tan(angle_fin)` à chaque extrémité (SPLIPH du Fortran).

    Comme le profil peut ne pas être une fonction univoque de x (le flanc de la fraise ou de la
    roue taillée peut « revenir en arrière »), la courbe est d'abord tournée pour rendre la
    pente imposée aussi horizontale que possible (x devient la spline naturelle, y la spline
    de Hermite), puis le résultat est tourné en sens inverse.
    """
    p = np.asarray(pts, dtype=float)
    rt = -(angle_debut + angle_fin) / 2.0
    if abs(angle_debut - angle_fin) > math.pi / 2:
        rt += math.pi / 2
    co, si = math.cos(rt), math.sin(rt)
    x = co * p[:, 0] - si * p[:, 1]
    y = si * p[:, 0] + co * p[:, 1]
    t = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    mx = _moments(t, x, None, None)
    ts = np.linspace(0.0, t[-1], n)
    sx, dxdt = _evaluer(t, x, mx, ts)
    d0 = dxdt[0] * math.tan(angle_debut + rt)
    d1 = dxdt[-1] * math.tan(angle_fin + rt)
    my = _moments(t, y, d0, d1)
    sy, _ = _evaluer(t, y, my, ts)
    return np.column_stack((co * sx + si * sy, -si * sx + co * sy))


# ---------------------------------------------------------------------------------------
# enveloppe (roue -> fraise, fraise -> roue) par balayage numérique
# ---------------------------------------------------------------------------------------

def _enveloppe(cx, cy, hauteurs):
    """Pour chaque valeur de `hauteurs` (m,), l'abscisse minimale parmi les intersections des
    segments d'une famille de courbes (cx, cy : (n_courbes, n_points)) avec l'horizontale à
    cette hauteur -- l'enveloppe inférieure de la famille (EINHÜLLENDE du Fortran).

    Seuls les segments *montants* comptent (`y(i) <= fy <= y(i+1)`, le test de ABWFRAC et de
    RADABWP), pas ceux qui redescendent. Ce n'est pas un détail : les profils sont donnés du
    pied vers la tête, donc un segment qui redescend appartient à une copie tellement tournée
    qu'elle s'est retournée. Accepter aussi ces segments-là fait prendre le minimum sur une
    copie située de l'autre côté de la roue, et le profil engendré part complètement ailleurs
    (roue à 8 dents de `test_jeuxmini.json` : rayon de tête 0.94 mm au lieu de 0.42).
    """
    x0, x1 = cx[:, :-1].ravel(), cx[:, 1:].ravel()
    y0, y1 = cy[:, :-1].ravel(), cy[:, 1:].ravel()
    montant = y0 < y1
    sortie = np.full(len(hauteurs), np.nan)
    for k, fy in enumerate(hauteurs):
        m = montant & (y0 <= fy) & (fy <= y1)
        if m.any():
            t = (fy - y0[m]) / (y1[m] - y0[m])
            sortie[k] = np.min(x0[m] + t * (x1[m] - x0[m]))
    return sortie


def _famille(pts, angle_pas, pas_lineaire, n_periodes, n_positions, inverse=False):
    """Famille de copies de `pts`, tournées de `wi=fg*angle_pas` et translatées de `ps=fg*pas_lineaire`
    (roulement sans glissement), pour `fg` couvrant `n_periodes` périodes consécutives (autour de 0)
    en `n_positions` pas chacune.

    Transformation "miroir" (-cos,-sin ; -sin,cos), pas une rotation standard : c'est la convention
    de ABWFRAC/RADABWP, empiriquement correcte ici aussi (vérifiée par le cas développante théorique
    -> droite à l'angle de pression, résidu 0.0002 mm au cœur du profil). Cette transformation est
    une involution (son carré est l'identité), mais elle n'est pas sa propre inverse *dans cette
    famille* : `profil_fraise` sème `pts` (repère roue, immobile) et translate le long de x (axe
    fixe de la crémaillère) ; `tailler` doit au contraire semer `pts` (repère fraise, immobile) et
    observer où ça tombe dans le repère de la roue *qui tourne* -- une translation d'axe fixe vue
    depuis un repère tournant est elle-même tournée, d'où `inverse=True` : translation
    `(cos(wi)*ps, sin(wi)*ps)` au lieu de `(ps, 0)`.
    """
    moitie = (n_periodes - 1) / 2.0
    l = np.arange(n_periodes) - moitie
    n = np.arange(n_positions) / max(1, n_positions - 1)
    fg = (l[:, None] + n[None, :]).ravel()
    wi, ps = fg * angle_pas, fg * pas_lineaire
    co, si = np.cos(wi), np.sin(wi)
    x, y = pts[:, 0], pts[:, 1]
    tx = ps * co if inverse else ps
    ty = ps * si if inverse else np.zeros_like(ps)
    cx = -co[:, None] * x[None, :] - si[:, None] * y[None, :] + tx[:, None]
    cy = -si[:, None] * x[None, :] + co[:, None] * y[None, :] + ty[:, None]
    return cx, cy


@dataclass
class ProfilFraise:
    """Profil d'une fraise de taillage par génération, dans son propre repère (origine sur la
    droite primitive, x = position le long de la crémaillère, y = hauteur au-dessus/en dessous).

    Comme une vraie fraise, ce profil ne dépend que du module et de l'angle de pression (donc du
    flanc d'origine) -- pas du nombre de dents de la roue qui l'a produit : la même fraise sert à
    tailler des roues de nombre de dents différent (`tailler`)."""

    nom: str
    module: float  # module (mm) : la fraise a une période pi*module*pc le long de la crémaillère
    pc: float  # rapport à la période standard (1 = standard) ; PC de KREBO/Fraesen.F90
    points: np.ndarray  # (n,2), pied à tête (même convention que Mobile.points)
    zz_origine: int  # nombre de dents de la roue qui a servi à calculer cette fraise (information)

    def echantillonner(self, n=400):
        return spline_hermite(self.points, 0.0, 0.0, n)


def profil_fraise(mobile, n_periodes=None, n_positions=61, pc=1.0, n_points=150, dos=False):
    """Fraise conjuguée du demi-profil menant de `mobile` (ou de son dos si `dos`) : porte
    `ABWFRAC`/`ABWFRAP`. `n_periodes` : nombre de pas de denture voisins pris en compte dans le
    balayage (impair ; None : `periodes_defaut(mobile.zz)`, adapté au nombre de dents -- un
    balayage total de plus de 90 degrés environ n'a plus de sens physique, voir `periodes_defaut`)
    ; `n_positions` : positions par pas (plus grand = enveloppe plus précise, plus lent) ; `pc` :
    période de la fraise en fraction de la période standard pi*module (1 = standard) ;
    `n_points` : points du profil de fraise renvoyé.
    """
    zz, mo = mobile.zz, mobile.mo
    n_periodes = periodes_defaut(zz) if n_periodes is None else n_periodes
    pts = mobile.profil_dos() if dos else mobile.points
    pm, ra = mobile.arcs_dos() if dos else mobile.arcs()
    fin = arcs.echantillonner(pts, pm, ra, N_ECHANT)
    ap, pl = 2.0 * math.pi / zz, pc * math.pi * mo
    cx, cy = _famille(fin, ap, pl, n_periodes, n_positions)
    r_pied = math.hypot(*fin[0])
    hauteur_dent = fin[-1, 1] - r_pied
    hauteurs = r_pied + hauteur_dent * np.linspace(0.0, 1.0, n_positions)[1:-1]
    fx = _enveloppe(cx, cy, hauteurs)
    pts_env = np.column_stack((fx, hauteurs))
    # extrémités analytiques (ABWFRAC) : elles se sont révélées correctes telles quelles dans ce
    # repère, vérifiées à moins de 0.3 µm contre une vraie fraise Fortran (Result/Validation_T.Frs)
    pts_env = np.vstack(([-pl / 2.0, r_pied], pts_env, [0.0, r_pied + hauteur_dent]))
    fp = spline_hermite(pts_env, 0.0, 0.0, n_points)
    fp[:, 1] -= zz * mo / 2.0  # repère roue -> repère fraise (origine sur la primitive)
    return ProfilFraise(f"Fraise-{mobile.nom}", mo, pc, fp, zz)


def tailler(fraise, zz, oy=0.0, n_periodes=None, n_positions=61, n_points=150):
    """Profil « pratique » du demi-profil menant d'une roue de `zz` dents taillée par génération
    avec `fraise` (porte `RADABWP`), déport `oy` (mm ; positif = fraise éloignée de la roue, comme
    un déport de denture positif) : la roue effectivement engendrée correspond au module
    `fraise.module + 2*oy/zz`. `oy=0` et `zz=fraise.zz_origine` reproduit (aux erreurs de
    discrétisation près) le profil qui a servi à calculer la fraise -- c'est le cas à comparer
    à la référence avec `ecart`. `n_periodes` : None -> `periodes_defaut(zz)`.
    """
    n_periodes = periodes_defaut(zz) if n_periodes is None else n_periodes
    r_roulement = zz * fraise.module / 2.0 + oy
    fin = fraise.echantillonner(N_ECHANT)
    fin = fin.copy()
    fin[:, 1] += r_roulement  # repère fraise -> repère roue cible
    ap, pl = 2.0 * math.pi / zz, fraise.pc * math.pi * fraise.module
    # dans le repère de la roue : la fraise (fixe) tourne de -wi et se translate de -ps (roulement
    # inverse de celui de profil_fraise, où c'était le flanc qui tournait et la fraise qui restait fixe)
    cx, cy = _famille(fin, ap, pl, n_periodes, n_positions, inverse=True)
    # pied et tête directement sur le profil de la fraise (pas sur le min/max de la famille, qui
    # mélange des copies de dents voisines). `creux` est un *rayon* ; le balayage et le point de
    # pied se font à la *hauteur* correspondante sur l'axe du creux, TY = creux*cos(pi/zz) dans
    # RADABWP (vérifié : son TY, pris dans la famille, vaut ça à moins de 10 nm sur 4 roues)
    creux, sommet = fin[0, 1], fin[-1, 1]
    ty = creux * math.cos(math.pi / zz)
    hauteurs = ty + (sommet - ty) * np.linspace(0.0, 1.0, n_positions)[1:-1]
    fx = _enveloppe(cx, cy, hauteurs)
    pts_env = np.column_stack((fx, hauteurs))
    # extrémités analytiques (RADABWP), par le même raisonnement que profil_fraise
    pts_env = np.vstack(([creux * math.sin(math.pi / zz), ty], pts_env, [0.0, sommet]))
    zp = spline_hermite(pts_env, -math.pi / zz, 0.0, n_points)
    return zp


def ecart(pratique, reference, n=800):
    """Écart signé maximal entre le profil `pratique` (taillé) et le profil `reference`, comme
    `KRVABP` : positif -> `pratique` est en dehors de la matière de `reference` (dent pas assez
    taillée), négatif -> `pratique` pénètre dans `reference` (dent trop taillée). Retourne
    (max positif, point, max négatif -- 0 si aucun --, point) : les deux profils sont
    échantillonnés (abscisse curviligne) et comparés point à point avec le pied de perpendiculaire
    le plus proche sur l'autre courbe.
    """
    def fin_spline(p):
        t = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        ts = np.linspace(0.0, t[-1], n)
        k = np.clip(np.searchsorted(t, ts, side="right") - 1, 0, len(t) - 2)
        u = ((ts - t[k]) / (t[k + 1] - t[k]))[:, None]
        return p[k] + u * (p[k + 1] - p[k])

    ref = fin_spline(np.asarray(reference, dtype=float))
    prat = np.asarray(pratique, dtype=float)
    seg0, seg1 = ref[:-1], ref[1:]
    d = seg1 - seg0
    long2 = np.maximum((d * d).sum(1), 1e-30)
    exces, penetre = 0.0, 0.0
    i_exces = i_penetre = 0
    for i, p in enumerate(prat):
        t = np.clip(((p - seg0) * d).sum(1) / long2, 0.0, 1.0)
        pied = seg0 + t[:, None] * d
        dist = np.hypot(*(p - pied).T)
        j = int(np.argmin(dist))
        # signe : produit vectoriel de la tangente locale de la référence avec (p - pied)
        signe = math.copysign(1.0, d[j, 0] * (p[1] - pied[j, 1]) - d[j, 1] * (p[0] - pied[j, 0]))
        ds = signe * dist[j]
        if ds > exces:
            exces, i_exces = ds, i
        if ds < penetre:
            penetre, i_penetre = ds, i
    return exces, tuple(prat[i_exces]), penetre, tuple(prat[i_penetre])


# ---------------------------------------------------------------------------------------
# fichier .Frs (une paire de fraises menante/menée, format Fortran de ENGOPT)
# ---------------------------------------------------------------------------------------

def ecrire_frs(chemin, fraises, hertz=None):
    """Écrit la paire `fraises` (menante, menée : deux `ProfilFraise`) dans `chemin` (.Frs),
    au format du Fortran (`ENGOPT`, option "F" : label A30, en-tête F10.6/I6/F10.6, nombre de
    points I6, points 2F12.7, CRLF) -- vérifié caractère pour caractère contre un fichier réel
    (`Result/Validation_T.Frs`). `hertz` : [(E, poisson, épaisseur)] des deux mobiles d'origine,
    si connues (facultatif, comme dans le fichier Fortran mais pas nécessaire pour retailler)."""
    lignes = []
    for f in fraises:
        lignes.append(f"{f.nom[:30]:<30s}")
        lignes.append(f"{f.module:10.6f}{f.zz_origine:6d}{f.pc:10.6f}     MO, ND, PC")
        lignes.append(f"{len(f.points):6d}")
        lignes += [f"{x:12.7f}{y:12.7f}" for x, y in f.points]
    if hertz:
        lignes.append(" Pression Hertzienne")
        lignes += [f"{e:10.1f}{poisson:10.3f}{epaisseur:10.3f}" for e, poisson, epaisseur in hertz]
    with open(chemin, "w", encoding="latin-1", errors="replace", newline="") as fh:
        fh.write("\r\n".join(lignes) + "\r\n")


def lire_frs(chemin):
    """Lit une paire de fraises (.Frs) : (menante, menée, hertz ou None). Tolérant aux espaces
    (lecture par `split()`, pas par colonnes fixes), lit aussi bien un fichier écrit par
    `ecrire_frs` qu'un fichier produit par le vrai programme Fortran (ENGOPT, option "F")."""
    lignes = Path(chemin).read_text(encoding="latin-1").splitlines()
    fraises, i = [], 0
    for _ in range(2):
        nom = lignes[i].strip()
        mo, zz, pc = lignes[i + 1].split()[:3]
        n = int(lignes[i + 2].split()[0])
        pts = np.array([[float(v) for v in lignes[i + 3 + k].split()[:2]] for k in range(n)])
        fraises.append(ProfilFraise(nom, float(mo), float(pc), pts, int(float(zz))))
        i += 3 + n
    hertz = None
    if i < len(lignes) and "hertzienne" in lignes[i].strip().lower():
        hertz = [tuple(float(v) for v in lignes[i + 1 + k].split()[:3]) for k in range(2)]
    return fraises[0], fraises[1], hertz
