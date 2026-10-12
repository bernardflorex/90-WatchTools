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

from . import arcs

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


# ---------------------------------------------------------------------------------------
# outil de scarving : roue-outil de rayon fini (la crémaillère ci-dessus en est la limite)
# ---------------------------------------------------------------------------------------
#
# Une couronne intérieure ne peut pas être taillée par une fraise-mère : il faut un outil qui
# passe dans l'ouverture centrale, donc une roue-outil. Le taillage par génération n'en est que le
# cas limite, et c'est `z_outil = 0` qui le désigne dans l'interface -- une crémaillère est une
# roue de rayon *infini*, mais « pas de diamètre » est ce que l'utilisateur tape pour dire
# « fraise », et 0 ne peut pas être un nombre de dents valide par ailleurs.
#
# Roulement sans glissement de la pièce (rayon primitif r1, signé : négatif pour une couronne) et
# de l'outil (r2 > 0) à l'entraxe signé A = r1 + r2 : pour une rotation u de la pièce, l'outil
# tourne de -(r1/r2)*u. Pièce centrée à l'origine, outil en (0, A) ; un point P du profil de la
# pièce s'écrit alors, dans le repère de l'outil,
#     Q(u) = R(u*A/r2) . P  +  (A*sin(u*r1/r2), -A*cos(u*r1/r2))
# Quand r2 tend vers l'infini cette forme dégénère (la translation part à l'infini avec A), d'où
# la crémaillère traitée à part par `_famille` -- celle-là est validée contre une vraie fraise
# Fortran et ne doit pas bouger.


def _rotation(p, a):
    """(n,2) tourné de `a` autour de l'origine."""
    c, s = math.cos(a), math.sin(a)
    return np.column_stack((p[:, 0] * c - p[:, 1] * s, p[:, 0] * s + p[:, 1] * c))


def _famille_outil(pts, r1, r2, pas_u, n_periodes, n_positions):
    """Famille des copies du demi-profil `pts` de la pièce, vue dans le repère de l'outil.

    `pas_u` : rotation de la pièce pour une de ses dents ; le balayage couvre `n_periodes` dents
    autour de 0, en `n_positions` pas chacune. Les dents de la pièce doivent regarder l'outil :
    le profil reçoit un demi-tour quand l'outil est de l'autre côté (couronne, A < 0), le même
    `delta1` que `contact.Engrenement`.
    """
    a = r1 + r2
    p = _rotation(np.asarray(pts, dtype=float), math.pi if a < 0.0 else 0.0)
    moitie = (n_periodes - 1) / 2.0
    l = np.arange(n_periodes) - moitie
    n = np.arange(n_positions) / max(1, n_positions - 1)
    us = ((l[:, None] + n[None, :]).ravel()) * pas_u
    qx = np.empty((len(us), len(p)))
    qy = np.empty((len(us), len(p)))
    for i, u in enumerate(us):
        q = _rotation(p, u * a / r2)
        qx[i] = q[:, 0] + a * math.sin(u * r1 / r2)
        qy[i] = q[:, 1] - a * math.cos(u * r1 / r2)
    return qx, qy


def _enveloppe_polaire(qx, qy, rayons):
    """Pour chaque rayon voulu (depuis le centre de l'outil), l'angle extrême parmi les
    intersections des segments de la famille avec ce cercle : l'enveloppe, en polaire.

    Transposition exacte de `_enveloppe` : l'angle y remplace l'abscisse et le rayon la hauteur,
    et seuls les segments dont le rayon *croît* comptent, pour la même raison (un segment qui
    décroît appartient à une copie assez tournée pour s'être retournée). Les angles sont mesurés
    depuis -pi/2, direction dans laquelle se trouve la zone de contact.
    """
    r = np.hypot(qx, qy)
    al = np.arctan2(qy, qx) + math.pi / 2.0
    al = (al + math.pi) % (2.0 * math.pi) - math.pi
    r0, r1_ = r[:, :-1].ravel(), r[:, 1:].ravel()
    a0, a1_ = al[:, :-1].ravel(), al[:, 1:].ravel()
    da = (a1_ - a0 + math.pi) % (2.0 * math.pi) - math.pi  # pas de saut sur +-pi
    montant = r0 < r1_
    sortie = np.full(len(rayons), np.nan)
    for k, fr in enumerate(rayons):
        m = montant & (r0 <= fr) & (fr <= r1_)
        if m.any():
            t = (fr - r0[m]) / (r1_[m] - r0[m])
            sortie[k] = np.min(a0[m] + t * da[m])
    return sortie


@dataclass
class ProfilOutil:
    """Profil d'un outil de scarving (roue-outil), dans son propre repère, origine au centre.

    `points` : le demi-profil (n,2), du pied à la tête comme pour une roue, l'axe du creux sur +y.
    `zz` : son nombre de dents (toujours > 0 : l'outil est une roue extérieure). `zz_origine`/`mo` :
    la pièce et le module qui l'ont produit, pour savoir quand la comparaison au profil de
    référence a un sens. `tronque` : la tête a été coupée à l'axe de dent, donc la dent de l'outil
    sort pointue et la pièce ne lui laisse pas de tête libre (voir `_rangerter`).
    """

    nom: str
    mo: float
    zz: int
    points: object
    zz_origine: int = 0
    tronque: bool = False
    sigma: float = 0.0  # angle de croisement des axes (rad) ; 0 = mortaisage, axes parallèles

    @property
    def rayon(self):
        return rayon_outil(self.mo, self.zz, self.sigma)

    @property
    def module_apparent(self):
        """Module dans le plan transversal de l'outil : `mo / cos(sigma)`. `mo` est son module
        *normal*, celui de la pièce."""
        return self.mo / math.cos(self.sigma)

    @property
    def module(self):
        """Alias de `mo` (module *normal*), pour que `ProfilOutil` et `ProfilFraise` se lisent
        pareil côté GUI."""
        return self.mo

    def echantillonner(self, n=N_ECHANT):
        # tangentes imposées comme pour une roue : perpendiculaire à l'axe du creux au départ,
        # perpendiculaire à l'axe de la dent à l'arrivée
        return spline_hermite(np.asarray(self.points, dtype=float),
                              -math.pi / self.zz, 0.0, n)


def _flanc_normales(points, centres, rayons, n=N_ECHANT):
    """Échantillonnage fin d'une chaîne d'arcs, avec la normale *exacte* en chaque point.

    La normale d'un arc va du centre au point ; elle sort de la matière du côté de la convexité.
    Les prendre sur les arcs plutôt que par différences finies compte : c'est la normale, pas la
    position, qui détermine le profil conjugué.
    """
    pts = np.asarray(points, dtype=float)
    pm, ra = np.asarray(centres, dtype=float), np.asarray(rayons, dtype=float)
    longueurs = []
    for i in range(len(ra)):
        a0 = math.atan2(pts[i, 1] - pm[i, 1], pts[i, 0] - pm[i, 0])
        a1 = math.atan2(pts[i + 1, 1] - pm[i, 1], pts[i + 1, 0] - pm[i, 0])
        longueurs.append(abs(ra[i] * ((a1 - a0 + math.pi) % (2 * math.pi) - math.pi)))
    total = sum(longueurs) or 1.0
    P, N = [], []
    for i in range(len(ra)):
        c, r = pm[i], abs(ra[i])
        a0 = math.atan2(pts[i, 1] - c[1], pts[i, 0] - c[0])
        a1 = math.atan2(pts[i + 1, 1] - c[1], pts[i + 1, 0] - c[0])
        da = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi
        k = max(8, int(round(n * longueurs[i] / total)))
        th = a0 + np.linspace(0.0, 1.0, k, endpoint=(i == len(ra) - 1)) * da
        u = np.column_stack((np.cos(th), np.sin(th)))
        P.append(c + r * u)
        N.append(u * math.copysign(1.0, ra[i]))
    P, N = np.vstack(P), np.vstack(N)
    # la matière est du côté x < 0 du demi-profil (la dent est centrée sur l'axe x = 0) : la
    # normale sortante pointe donc vers les x croissants
    if N[:, 0].mean() < 0.0:
        N = -N
    return P, N


def rayon_outil(mo, z_outil, sigma=0.0):
    """Rayon primitif de l'outil : `mo z / 2` à axes parallèles, `mo z / (2 cos(sigma))` dès que
    les axes sont croisés de `sigma`.

    C'est la relation classique des roues gauches (hélicoïdales à axes croisés), et elle tombe de
    la condition de roulement sur la perpendiculaire commune (voir `conjuguer`) : pour une pièce
    droite, l'outil est une roue hélicoïdale d'angle d'hélice `sigma`, de module *normal* `mo` et
    donc de module apparent `mo / cos(sigma)`.
    """
    return mo * z_outil / (2.0 * math.cos(sigma))


def _resoudre_contact(cr, nu, p, phi, r1, r2, sigma, u0, h=0.0, iterations=60):
    """Rotation `u` de la pièce qui amène chaque point au contact, axes croisés de `sigma`.

    Équation d'engrènement générale `N . v_rel = 0` (voir `conjuguer`), restreinte au plan
    transversal de l'outil :
    `cr (1 - i cos S) = i [ A cos S cos(nu+u) + tan S sin(nu+u) (sin S p cos(phi+u) - h) ]`
    avec `i = -r1 / (r2 cos S)` le rapport des vitesses, `A = r1 + r2` l'entraxe signé et `h` la
    cote axiale de la section d'outil regardée (`h = 0` : celle qui passe par la perpendiculaire
    commune).

    Pas de forme fermée (il y a du `u` et du `2u`), mais la solution est continue en `sigma` :
    on part de la racine exacte `u0` à axes parallèles et on fait du Newton. C'est aussi ce qui
    choisit la branche, sans vote ni heuristique -- `u0` la porte déjà.
    """
    s, c = math.sin(sigma), math.cos(sigma)
    a = r1 + r2
    i = -r1 / (r2 * c)
    k = s * s / c          # sin(S) tan(S)
    th = s / c             # tan(S)
    gauche = cr * (1.0 - i * c)
    u = np.array(u0, dtype=float)
    for _ in range(iterations):
        f = (gauche - i * (a * c * np.cos(nu + u)
                           + (k * p * np.cos(phi + u) - th * h) * np.sin(nu + u)))
        d = i * (a * c * np.sin(nu + u) + th * h * np.cos(nu + u) - k * p * np.cos(nu + phi + 2.0 * u))
        pas = np.where(np.abs(d) > 1e-12, f / np.where(np.abs(d) > 1e-12, d, 1.0), 0.0)
        u = u - np.clip(pas, -0.2, 0.2)
    f = (gauche - i * (a * c * np.cos(nu + u)
                       + (k * p * np.cos(phi + u) - th * h) * np.sin(nu + u)))
    return u, np.abs(f) < 1e-9 * max(1.0, abs(a))


def conjuguer(points, normales, r1, r2, sigma=0.0, h=0.0):
    """Profil conjugué de (`points`, `normales`) -- repère de la pièce, centre à l'origine -- vu
    dans le plan transversal de l'outil de rayon primitif `r2`, axes croisés de `sigma` (rad).

    `sigma = 0` : axes parallèles, c'est le mortaisage (et, à la limite `r2 -> inf`, la
    fraise-mère). **Le vrai skiving a besoin de `sigma != 0`** : c'est l'angle de croisement qui
    donne la vitesse de coupe, un outil parallèle ne coupant pas. Tout ce qui suit vaut pour une
    pièce à denture *droite*, ce qui est le seul cas du projet.

    Équation d'engrènement. Axe de la pièce = z, axe de l'outil passant par `A = (0, a, 0)` et
    dirigé par `k2 = (sin S, 0, cos S)` ; la perpendiculaire commune est l'axe y, donc `a` est bien
    l'entraxe, signé comme partout ici. Le contact exige `N . v_rel = 0` avec
    `v_rel = w1 (k1 x P) - w2 (k2 x (P - A))`, ce qui donne, la pièce étant droite (sa normale n'a
    pas de composante axiale) :
    `cr (1 - i cos S) = i [ A cos S Nx - sin S Ny Pz ]`, `cr = P x N`, `i = w2/w1`.
    Le rapport vient du roulement sur la perpendiculaire commune : `i = -r1 / (r2 cos S)`, d'où
    `r2 = mo z2 / (2 cos S)` (`rayon_outil`).

    À `S = 0` cela redonne exactement `cos(u + nu) = -cr/r1`, la loi de l'engrènement du cas plan
    (la normale commune passe par le point primitif), résolue en forme fermée -- le chemin
    `sigma = 0` est donc inchangé, au bit près.

    Le flanc de la pièce étant un cylindre, *chaque* point est touché à *une* hauteur axiale `Pz` :
    le contact est ponctuel et se promène le long de l'axe, ce qui est précisément la cinématique
    du skiving (l'outil avance axialement). On prend la section de l'outil par son propre plan
    transversal : la section à la cote axiale `h` impose `Px sin S + Pz cos S = h`, et le point
    conjugué y vaut `R(-i u) . (Px / cos S - h tan S, Py - A)` -- le `1/cos S` est l'allongement
    apparent, le même qui fait passer le module de `mo` à `mo / cos S`. `h = 0` est la section qui
    passe par la perpendiculaire commune ; les autres servent à vérifier si la surface conjuguée
    est bien un hélicoïde (voir CLAUDE.md : elle ne l'est pas tout à fait, et c'est le fait connu
    qui empêche un outil de skiving d'être rigoureusement conjugué).

    Les points sans solution sont rendus en NaN.
    """
    P = np.asarray(points, dtype=float)
    N = np.asarray(normales, dtype=float)
    N = N / np.hypot(*N.T)[:, None]
    a = r1 + r2
    cr = P[:, 0] * N[:, 1] - P[:, 1] * N[:, 0]
    nu = np.arctan2(N[:, 1], N[:, 0])
    c = -cr / r1
    possible = np.abs(c) <= 1.0
    # `cos(u + nu) = c` a deux racines, `u = -nu +- arccos(c)`, qui ne diffèrent que par le signe
    # de la composante y de la normale au contact : `sin(u + nu) = +-sqrt(1-c^2)`. La bonne est
    # toujours celle où cette composante est **positive**, car la normale sortante de la pièce doit
    # pointer vers l'outil, et l'outil est du côté +y du point primitif pour les deux signes de
    # `r1` (roue extérieure : son centre est au-delà ; couronne : la normale sortante pointe vers
    # le centre de la couronne, donc vers l'outil qui est dedans). Donc `arccos`, qui rend [0, pi],
    # et rien d'autre : aucun vote, aucune heuristique. (Un vote sur la direction de la normale se
    # trompait dès que l'outil était petit, un vote sur la plus petite rotation se trompait près du
    # cercle de base d'une petite roue extérieure.)
    u = -nu + np.arccos(np.clip(c, -1.0, 1.0))
    if sigma:
        rho = np.hypot(*P.T)
        phi = np.arctan2(P[:, 1], P[:, 0])
        u, ok = _resoudre_contact(cr, nu, rho, phi, r1, r2, sigma, u, h)
        possible &= ok
    cu, su = np.cos(u), np.sin(u)
    w = np.column_stack((P[:, 0] * cu - P[:, 1] * su, P[:, 0] * su + P[:, 1] * cu))
    j = r1 / (r2 * math.cos(sigma))          # -i : rotation de l'outil pour une rotation u
    cj, sj = np.cos(j * u), np.sin(j * u)
    wx, wy = w[:, 0] / math.cos(sigma) - h * math.tan(sigma), w[:, 1] - a
    q = np.column_stack((wx * cj - wy * sj, wx * sj + wy * cj))
    q[~possible] = np.nan
    return q, u


def _rangerter(q, z_outil, n_points):
    """Nuage conjugué remis dans la convention d'un demi-profil de roue -- pied d'abord, flanc du
    côté x > 0 -- tronqué à la dent de l'outil et rééchantillonné à pas d'abscisse curviligne
    constant. Rend `(points, tronque)`, `tronque` disant si la tête a dû être coupée.

    Repère : `conjuguer` rend le point dans le repère de l'outil, et **+y y est l'axe du creux**
    de l'outil pour les deux signes de la pièce (`u = 0`, c'est la pièce dans son rangement, donc
    une de ses dents sur la ligne des centres, et une dent de la pièce est en face d'un creux de
    l'outil). Après le demi-tour, le flanc part donc de l'axe du creux et monte vers la tête.

    Deux troncatures, et aucune extrémité inventée -- le conjugué est le *flanc actif* de l'outil,
    son pied et sa tête restant des choix de conception :
    - **l'axe de dent**, à un demi-pas (`pi/z`) de l'axe du creux : au-delà, le flanc appartient au
      flanc d'en face de la même dent. La dent sort pointue, et c'est la seule tête que la pièce
      impose. Sans cette coupe, le conjugué du congé de pied de la pièce dépasse et croise le flanc
      du creux voisin -- le crochet visible au sommet de chaque dent dessinée ;
    - **le repli** : le conjugué du congé de pied de la pièce rebrousse (rien n'y touche, mais la
      loi de l'engrènement lui donne quand même une solution). Le flanc actif est monotone en
      rayon, on garde la plage monotone qui part du pied.
    """
    q = np.asarray(q, dtype=float)
    q = q[np.isfinite(q[:, 0]) & np.isfinite(q[:, 1])]
    if len(q) < 8:
        raise ValueError("profil conjugué vide (outil incompatible avec cette pièce)")
    q = -q  # demi-tour : le flanc regarde +y
    if q[:, 0].mean() < 0.0:
        q[:, 0] = -q[:, 0]
    # repli : on garde la plus longue plage monotone en rayon. Pas « celle qui part du bout » :
    # selon la pièce le repli est en tête ou en pied du tableau, et les deux extrémités peuvent
    # être dans la partie repliée (mesuré : pièce Z=30 avec un outil Z=12, où partir du bout donne
    # 0.648..0.763 mm au lieu de 0.769..1.232).
    r = np.hypot(*q.T)
    d = np.diff(r) > 0.0
    debut = fin = meilleur = 0
    i = 0
    while i < len(d):
        j = i
        while j < len(d) and d[j] == d[i]:
            j += 1
        if j - i > meilleur:
            meilleur, debut, fin = j - i, i, j
        i = j
    q = q[debut:fin + 1]
    if len(q) < 8:
        raise ValueError("flanc actif conjugué trop court")
    if np.hypot(*q[0]) > np.hypot(*q[-1]):  # du pied vers la tête
        q = q[::-1]
    # axe de dent : angle polaire depuis +y, coupé à pi/z (interpolation sur le dernier segment)
    demi = math.pi / abs(z_outil)
    a = np.abs(np.arctan2(q[:, 0], q[:, 1]))
    tronque = bool(a[-1] > demi)
    if tronque:
        k = int(np.argmax(a > demi))
        if k < 2:
            raise ValueError("flanc actif conjugué trop court (outil trop petit pour cette pièce)")
        f = (demi - a[k - 1]) / (a[k] - a[k - 1])
        q = np.vstack((q[:k], q[k - 1] + f * (q[k] - q[k - 1])))
    if len(q) < 8:
        raise ValueError("flanc actif conjugué trop court")
    s = np.concatenate(([0.0], np.cumsum(np.hypot(*np.diff(q, axis=0).T))))
    cible = np.linspace(0.0, s[-1], n_points)
    return np.column_stack((np.interp(cible, s, q[:, 0]), np.interp(cible, s, q[:, 1]))), tronque


def _famille_piece(q, r1, r2, us, sigma=0.0):
    """Copies du profil `q` de l'outil (sa section transversale, axe du creux sur +y) vues dans le
    plan transversal de la pièce (centre à l'origine, point primitif en `(0, r1)`), pour les
    rotations `us` de la pièce, axes croisés de `sigma`.

    C'est l'inverse exact de `conjuguer`, pas une cinématique re-dérivée : celle-ci envoie un
    point `P` de la pièce, au contact, sur `R(-i u).(Px/cos S, Py - A)` avec `i = -r1/(r2 cos S)`
    et `A = r1 + r2`, donc `P = R(-u).(cos S (R(i u) q)_x, (R(i u) q)_y + A)`. Les deux signes de
    `r1` passent sans cas particulier, comme dans `conjuguer`.

    **Une seule section d'outil suffit**, et c'est physique, pas une approximation de confort : en
    skiving c'est l'*arête* de coupe qui enlève la matière, pas le flanc, et avec un angle de coupe
    nul cette arête est précisément la section transversale de l'outil dans le plan de coupe. Une
    pièce droite étant invariante le long de son axe, l'avance axiale ne change rien au profil
    transversal engendré, et l'arête reste à la même cote dans le repère de l'outil. Ce qui *ne*
    reste pas constant, c'est le profil de l'arête après réaffûtage (voir `derive_affutage`).
    """
    q = np.asarray(q, dtype=float)
    a, j = r1 + r2, r1 / (r2 * math.cos(sigma))
    cs = math.cos(sigma)
    px = np.empty((len(us), len(q)))
    py = np.empty((len(us), len(q)))
    for i, u in enumerate(us):
        c, s = math.cos(j * u), math.sin(j * u)          # R(i u) = R(-j u)
        x = (q[:, 0] * c + q[:, 1] * s) * cs
        y = -q[:, 0] * s + q[:, 1] * c + a
        c, s = math.cos(u), math.sin(u)                  # R(-u)
        px[i] = x * c + y * s
        py[i] = -x * s + y * c
    return px, py


def _enveloppe_angle(px, py, rayons, axe):
    """Pour chaque rayon voulu, le plus petit angle à l'axe de dent `axe` (vecteur unitaire)
    atteint par la famille : l'outil mange la dent, le flanc taillé est donc au plus près de
    l'axe. Transposition de `_enveloppe`, l'angle remplaçant l'abscisse et le rayon la hauteur,
    avec le même filtre : seuls comptent les segments dont le rayon *croît* le long du profil de
    l'outil rangé pied -> tête de la pièce, un segment décroissant appartenant à une copie assez
    tournée pour s'être retournée (c'est le défaut qui avait fait dériver `_enveloppe`).
    """
    r = np.hypot(px, py)
    be = np.arctan2(px * axe[1] - py * axe[0], px * axe[0] + py * axe[1])
    r0, r1_ = r[:, :-1].ravel(), r[:, 1:].ravel()
    b0, b1_ = be[:, :-1].ravel(), be[:, 1:].ravel()
    db = (b1_ - b0 + math.pi) % (2.0 * math.pi) - math.pi
    montant = r0 < r1_
    sortie = np.full(len(rayons), np.nan)
    for k, fr in enumerate(rayons):
        m = montant & (r0 <= fr) & (fr <= r1_)
        if m.any():
            t = (fr - r0[m]) / (r1_[m] - r0[m])
            sortie[k] = np.min(np.abs(b0[m] + t * db[m]))
    return sortie


def derive_affutage(mobile, z_outil, sigma, h=0.05, n_points=400, dos=False, n=N_ECHANT):
    """Changement de *forme* (mm) entre le profil de l'outil dans son plan de coupe et celui à la
    cote axiale `h`, à rotation rigide près : de combien le profil de l'arête change quand on
    réaffûte l'outil de `h` sur sa face de coupe.

    Le fond du problème : la surface conjuguée exacte d'une pièce droite, à axes croisés, **n'est
    pas tout à fait un hélicoïde**, alors qu'un outil est taillé et rectifié comme une roue
    hélicoïdale. L'outil ne peut donc être rigoureusement conjugué que dans *un* plan ; ailleurs il
    dérive, d'autant plus vite que `sigma` est grand. C'est le défaut connu du skiving, et ce
    nombre est ce qu'il faut regarder pour choisir `sigma` (à `sigma = 0` il est nul : le
    mortaisage n'a pas ce problème, mais il ne coupe pas).

    La rotation rigide entre les deux sections -- à peu près `h tan(sigma) / r2`, le pas de
    l'hélice -- est enlevée : elle ne coûte rien, un décalage angulaire de l'outil se rattrapant
    par la phase de la machine. Elle est *mesurée* à rayon égal, pas posée.

    **Mesure locale** : `h` doit rester petit (jusqu'à environ 0.1 mm). Au-delà les deux sections
    ne se recouvrent plus assez en rayon et `_rangerter` ne les tronque plus au même endroit ; ce
    qui lâche alors est la mesure, pas la géométrie. Le résultat est quasi proportionnel à `h`
    (vérifié de 0.02 à 0.1 mm), donc il s'annonce bien comme un taux, en µm par mm de réaffûtage.
    """
    if not sigma:
        return 0.0
    pts = mobile.profil_dos() if dos else mobile.points
    pm, ra = mobile.arcs_dos() if dos else mobile.arcs()
    P, N = _flanc_normales(pts, pm, arcs.convexite(mobile.zz, ra), n)
    r1, r2 = mobile.mo * mobile.zz / 2.0, rayon_outil(mobile.mo, z_outil, sigma)
    if r1 < 0.0:
        P, N = -P, -N
    # `_rangerter` : sans lui les deux nuages portent encore le repli du congé, le rayon n'y est
    # pas monotone et toute mesure à rayon égal est fausse (c'était une erreur de mesure, pas de
    # géométrie : elle donnait des centaines de microns au lieu de quelques-uns)
    # On ramène la section `h` par la rotation d'hélice **avant** de la ranger : `_rangerter`
    # décide du côté du flanc sur le signe de x moyen, et une section très tournée se ferait
    # retourner, ce qui donnait des centaines de microns de pure orientation.
    qh = conjuguer(P, N, r1, r2, sigma, h)[0]
    beta = h * math.tan(sigma) / r2
    c, s = math.cos(beta), math.sin(beta)
    qh = np.column_stack((qh[:, 0] * c + qh[:, 1] * s, -qh[:, 0] * s + qh[:, 1] * c))
    try:
        a = _rangerter(conjuguer(P, N, r1, r2, sigma, 0.0)[0], z_outil, n_points)[0]
        b = _rangerter(qh, z_outil, n_points)[0]
    except ValueError:
        return float("nan")
    ra_, rb = np.hypot(*a.T), np.hypot(*b.T)
    lo, hi = max(ra_.min(), rb.min()) * 1.001, min(ra_.max(), rb.max()) * 0.999
    if not hi > lo:
        return float("nan")
    bande = np.linspace(lo, hi, 80)

    def angle_a(q, r):
        o = np.argsort(r)
        return np.interp(bande, r[o], np.arctan2(q[:, 0], q[:, 1])[o])

    # Écart **à rayon égal**, dont on retire la moyenne : ce qui reste est le changement de forme,
    # converti en longueur sur le flanc. Mesurer plutôt une distance point-courbe serait faux ici,
    # les deux sections n'ayant ni la même étendue radiale ni la même troncature de tête.
    d = angle_a(b, rb) - angle_a(a, ra_)
    return float(np.ptp(d)) * float(np.mean(bande))


def tailler_outil(outil, zz, oy=0.0, n_periodes=None, n_positions=61, n_points=150):
    """Demi-profil de la roue de `zz` dents (signé : < 0 = couronne) taillée par l'outil de
    scarving `outil`, dans la convention de rangement habituelle (dent sur +y, flanc x > 0).

    Pendant de `tailler` pour un outil de rayon fini : l'outil roule sans glisser sur la pièce et
    son flanc balaie une famille, dont l'enveloppe est le flanc engendré. `oy` (mm) déplace la
    pièce de l'outil, donc augmente `|entraxe|` d'autant -- c'est le déport de denture, l'exact
    équivalent de l'éloignement de la crémaillère dans `tailler`. `oy = 0` et `zz = outil.zz_origine`
    redonne (aux erreurs de discrétisation près) le profil qui avait produit l'outil.

    Deux points à connaître :
    - **la tête de l'outil**. Si `outil.tronque`, la dent de l'outil est pointue et son flanc va
      jusqu'à l'axe de dent : c'est tout l'outil, rien à ajouter. Sinon la pièce ne lui impose pas
      de tête, et il faut bien en choisir une pour tailler : on ferme par un **arc de tête
      circulaire** au rayon du dernier point du flanc, jusqu'à l'axe de dent -- le choix minimal,
      qui n'invente aucune forme. Un vrai outil aurait un congé de tête, qui ne changerait que le
      fond du creux de la pièce.
    - **le pied et la tête de la pièce** viennent du profil de l'outil (`rayon_piece` ci-dessous),
      jamais du min/max de la famille, qui mélange des copies de dents voisines -- la même règle
      que `tailler`, où l'ignorer donnait 400 µm d'erreur sur le pied.
    """
    zz = int(zz)
    n_periodes = periodes_defaut(zz) if n_periodes is None else n_periodes
    sigma = getattr(outil, "sigma", 0.0)
    r2 = outil.rayon                     # tient déjà compte de sigma (`rayon_outil`)
    r1 = outil.mo * zz / 2.0
    a = r1 + r2
    a += math.copysign(oy, a)
    if abs(a) <= 1e-9:
        raise ValueError("entraxe nul : outil et pièce de rayons opposés")
    sg = math.copysign(1.0, r1)
    q = np.asarray(outil.points, dtype=float)
    if not outil.tronque:                       # arc de tête circulaire jusqu'à l'axe de dent
        rt = float(np.hypot(*q[-1]))
        b0 = math.atan2(q[-1, 0], q[-1, 1])
        bs = np.linspace(b0, math.copysign(math.pi / outil.zz, b0), 24)[1:]
        q = np.vstack((q, np.column_stack((rt * np.sin(bs), rt * np.cos(bs)))))
    rq = np.hypot(*q.T)
    q = -q                                      # défaire le demi-tour de rangement de _rangerter
    rayon_piece = lambda r: abs(a) - sg * r     # rayon pièce d'un point d'outil au rayon r
    if rayon_piece(rq[0]) > rayon_piece(rq[-1]):
        q, rq = q[::-1], rq[::-1]               # rangé pied -> tête de la *pièce*
    pas_u = 2.0 * math.pi / abs(zz)             # une dent d'outil = une dent de pièce
    moitie = (n_periodes - 1) / 2.0
    l = np.arange(n_periodes) - moitie
    n = np.arange(n_positions) / max(1, n_positions - 1)
    us = ((l[:, None] + n[None, :]).ravel()) * pas_u
    px, py = _famille_piece(q, r1, r2, us, sigma)
    if r1 < 0.0:                                # repère de rangement : dent sur +y dans les 2 cas
        px, py = -px, -py
    axe = (0.0, 1.0)
    pied, tete = rayon_piece(rq[-1]), rayon_piece(rq[0])
    rayons = pied + (tete - pied) * np.linspace(0.0, 1.0, max(32, n_positions))
    be = _enveloppe_angle(px, py, rayons, axe)
    bon = np.isfinite(be)
    if bon.sum() < 8:
        raise ValueError("enveloppe vide : cet outil ne taille pas cette roue")
    be, rayons = be[bon], rayons[bon]
    pts = np.column_stack((rayons * np.sin(be), rayons * np.cos(be)))
    # Tangentes de bout **mesurées sur l'enveloppe**, pas imposées. `tailler` peut imposer les
    # siennes (-pi/z au pied, 0 à la tête) parce qu'il ajoute d'abord les deux points analytiques
    # de `RADABWP`, sur l'axe du creux et sur l'axe de la dent. Ici l'enveloppe ne va ni jusqu'à
    # l'un ni jusqu'à l'autre -- le fond du creux et la tête de la pièce sont des rayons de brut,
    # que l'outil n'engendre pas -- et les inventer fausse les deux bouts (mesuré : 70 à 180 µm,
    # contre 18 µm pour l'effet de bord de la spline seule, lui-même supprimé par ces tangentes).
    def pente(p0, p1):
        return math.atan2(p1[1] - p0[1], p1[0] - p0[0])
    return spline_hermite(pts, pente(pts[0], pts[1]), pente(pts[-2], pts[-1]), n_points)


def profil_outil(mobile, z_outil, sigma=0.0, n_points=150, dos=False, n=N_ECHANT, **_ignore):
    """Profil de l'outil de scarving (roue de `z_outil` dents, même module) conjugué du demi-profil
    menant de `mobile`, ou de son dos si `dos`.

    `z_outil = 0` : crémaillère, c'est-à-dire la fraise de taillage par génération ; la fonction
    délègue alors à `profil_fraise`, validé contre une vraie fraise Fortran. Sinon le conjugué est
    calculé par `conjuguer` (loi de l'engrènement) -- exact, sans paramètre de résolution, d'où
    l'absence de `n_periodes`/`n_positions` ici."""
    if z_outil == 0:
        return profil_fraise(mobile, None, 61, 1.0, n_points, dos)
    if z_outil < 3:
        raise ValueError("nombre de dents de l'outil : 0 (crémaillère) ou au moins 3")
    if not -1.4 < sigma < 1.4:
        raise ValueError("angle de croisement : moins de 80 degrés en valeur absolue")
    pts = mobile.profil_dos() if dos else mobile.points
    pm, ra = mobile.arcs_dos() if dos else mobile.arcs()
    P, N = _flanc_normales(pts, pm, arcs.convexite(mobile.zz, ra), n)
    r1, r2 = mobile.mo * mobile.zz / 2.0, rayon_outil(mobile.mo, z_outil, sigma)
    if r1 < 0.0:
        # Le demi-profil est toujours rangé dent vers +y, alors que le point primitif est en
        # (0, r1), donc en -y pour une couronne : on lui donne le demi-tour que son signe exige,
        # exactement le `delta1` de `contact.Engrenement`. Sans lui la rotation de contact vaut
        # ~pi, les deux branches sont départagées à l'envers et le profil sort au mauvais rayon.
        P, N = -P, -N
    q, _ = conjuguer(P, N, r1, r2, sigma)
    pts, tronque = _rangerter(q, z_outil, n_points)
    return ProfilOutil(f"Outil-{mobile.nom}", mobile.mo, int(z_outil), pts, mobile.zz, tronque,
                       float(sigma))
