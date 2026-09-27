"""Export DXF des profils (mobiles, fraises, roues taillées reconverties en arcs).

Une chaîne d'arcs de cercle (mobile en arcs, roue taillée reconvertie en arcs) est exportée en
véritables entités ARC/LINE (via `ezdxf`), pas en polyligne échantillonnée : tout ce programme
tourne autour de l'approximation des profils par de vrais arcs de cercle, l'export doit le
refléter (utile en particulier pour l'usinage : mouvements d'arc G02/G03, pas une suite de
petits segments). Le profil de la fraise, qui n'est pas une chaîne d'arcs (c'est l'enveloppe
d'un balayage, une courbe quelconque), est exporté en LWPOLYLINE finement échantillonnée.
"""

import math

import ezdxf
import ezdxf.units

RAYON_DROITE = 1e4  # mm : au-delà, un arc est en pratique une droite (ex. raccord de pied)


def _angle_deg(p, c):
    return math.degrees(math.atan2(p[1] - c[1], p[0] - c[0]))


def _rot(p, a):
    co, si = math.cos(a), math.sin(a)
    return (p[0] * co - p[1] * si, p[0] * si + p[1] * co)


def _nouveau_doc():
    doc = ezdxf.new("R2010")
    doc.units = ezdxf.units.MM
    doc.header["$INSUNITS"] = ezdxf.units.MM
    return doc


def _ajouter_chaine_arcs(msp, points, centres, rayons, decalage_angle=0.0, mirror=False, layer="0"):
    """Ajoute au modelspace `msp` la chaîne d'arcs (points de transition, centres, rayons
    signés, dans le repère d'un mobile), tournée de `decalage_angle` (rad) autour de l'origine
    et mirroée en x si `mirror` (pour le flanc opposé, comme `evolvente.contour_dents`)."""
    s = -1.0 if mirror else 1.0
    for i in range(len(rayons)):
        p = _rot((s * points[i][0], points[i][1]), decalage_angle)
        q = _rot((s * points[i + 1][0], points[i + 1][1]), decalage_angle)
        c = _rot((s * centres[i][0], centres[i][1]), decalage_angle)
        r = abs(rayons[i])
        if r > RAYON_DROITE:  # rayon énorme (quasi rectiligne) : une droite, pas un arc
            msp.add_line(p, q, dxfattribs={"layer": layer})
            continue
        a0, a1 = _angle_deg(p, c), _angle_deg(q, c)
        # DXF trace toujours l'arc de a0 vers a1 dans le sens trigonométrique : si le balayage
        # réel (le plus court, p -> q) va dans l'autre sens, on échange les deux angles
        if (a1 - a0 + 180.0) % 360.0 - 180.0 < 0:
            a0, a1 = a1, a0
        msp.add_arc(center=c, radius=r, start_angle=a0, end_angle=a1, dxfattribs={"layer": layer})


def ecrire_profil(chemin, mobile, n_dents=None, layer="PROFIL"):
    """Exporte le contour complet (`n_dents` dents, toutes par défaut) du mobile en DXF, en
    arcs de cercle réels : le flanc menant (`mobile.points`/`arcs()`) et son dos (identiques si
    la dent est symétrique, sinon `profil_dos()`/`arcs_dos()`), répétés et tournés pour chaque
    dent -- comme `evolvente.contour_dents`, mais avec les vrais arcs plutôt que des points."""
    doc = _nouveau_doc()
    msp = doc.modelspace()
    pm, ra = mobile.arcs()
    pts = mobile.points
    if mobile.symetrique:
        pm_d, ra_d, pts_d = pm, ra, pts
    else:
        pm_d, ra_d = mobile.arcs_dos()
        pts_d = mobile.profil_dos()
    pas = 2.0 * math.pi / mobile.zz
    for i in range(n_dents or mobile.zz):
        a = i * pas
        _ajouter_chaine_arcs(msp, pts, pm, ra, a, mirror=True, layer=layer)
        _ajouter_chaine_arcs(msp, pts_d, pm_d, ra_d, a, mirror=False, layer=layer)
    doc.saveas(chemin)


def ecrire_polyligne(chemin, points, zz, n_dents=None, layer="PROFIL"):
    """Exporte le contour complet d'une roue en LWPOLYLINE, à partir d'un demi-profil donné en
    suite de points (repère de `evolvente.evolve` : départ sur l'axe du creux, fin sur l'axe de
    la dent). Une dent = le demi-profil et son miroir ; la dent est répétée `n_dents` fois
    (toutes par défaut) autour du centre.

    C'est la forme qui convient à un profil qui n'est pas une chaîne d'arcs -- une roue taillée
    par `fraise.tailler` en particulier : la fraise engendre une enveloppe, pas des arcs, et une
    polyligne la décrit telle quelle. Pour une chaîne d'arcs, `ecrire_profil` écrit de vraies
    entités ARC, plus légères et exactes.
    """
    doc = _nouveau_doc()
    msp = doc.modelspace()
    pts = [(float(x), float(y)) for x, y in points]
    # une dent : du creux à la tête, puis le miroir de la tête au creux suivant ; le dernier
    # point d'une dent et le premier de la suivante sont confondus sur l'axe du creux, donc le
    # contour se referme dent après dent
    dent = pts + [(-x, y) for x, y in reversed(pts)][1:]
    pas = 2.0 * math.pi / zz
    for i in range(n_dents or zz):
        msp.add_lwpolyline([_rot(p, i * pas) for p in dent], dxfattribs={"layer": layer})
    doc.saveas(chemin)


def ecrire_courbe(chemin, points, n_periodes=1, pas=None, layer="FRAISE"):
    """Exporte une courbe quelconque déjà échantillonnée (ex. le profil d'une fraise, un demi-
    pas) en LWPOLYLINE : mirroée pour la dent complète, répétée sur `n_periodes` avec un
    décalage `pas` entre chacune (utile pour visualiser plusieurs dents de la fraise)."""
    doc = _nouveau_doc()
    msp = doc.modelspace()
    pts = [(float(x), float(y)) for x, y in points]
    gauche = [(-x, y) for x, y in reversed(pts)]
    dent = pts + gauche[1:]
    for k in range(n_periodes):
        decalage = (k - (n_periodes - 1) / 2.0) * (pas or 0.0)
        msp.add_lwpolyline([(x + decalage, y) for x, y in dent], dxfattribs={"layer": layer})
    doc.saveas(chemin)
