"""Lecture / ecriture DXF minimale (ASCII, sans dependance externe).

Import  : recupere les cercles (rayon mini / encombrement maxi), le pivot du
          marteau et la droite du plan de contact pour en deduire sa pente.
Export  : profil du coeur, cercles de construction, plan de contact et pivot.
"""

import math

from .core import Params, dot, norm, u_of

TAU = 2.0 * math.pi


# ----------------------------------------------------------------- lecture

def read_entities(path):
    """Renvoie la liste des entites utiles : dicts {type, ...}."""
    with open(path, "r", errors="replace") as f:
        raw = [l.strip() for l in f.read().splitlines()]
    pairs = []
    for i in range(0, len(raw) - 1, 2):
        try:
            pairs.append((int(raw[i]), raw[i + 1]))
        except ValueError:
            continue

    ents = []
    cur = None
    in_ent = False
    _ = in_ent
    for code, val in pairs:
        if code == 0:
            if cur:
                ents.append(cur)
            cur = None
            if val == "SECTION" or val == "ENDSEC":
                in_ent = (val == "SECTION")
            if val in ("CIRCLE", "LINE", "POINT", "ARC", "LWPOLYLINE",
                       "POLYLINE", "VERTEX"):
                cur = {"type": val, "x": [], "y": [], "r": None,
                       "layer": "", "a0": None, "a1": None}
            continue
        if cur is None:
            continue
        if code == 8:
            cur["layer"] = val
            continue
        try:
            fv = float(val)
        except ValueError:
            continue
        if code in (10, 11):
            cur["x"].append(fv)
        elif code in (20, 21):
            cur["y"].append(fv)
        elif code == 40:
            cur["r"] = fv
        elif code == 50:
            cur["a0"] = fv
        elif code == 51:
            cur["a1"] = fv
    if cur:
        ents.append(cur)

    out = []
    for e in ents:
        n = min(len(e["x"]), len(e["y"]))
        e["pts"] = [(e["x"][i], e["y"][i]) for i in range(n)]
        e["a0"] = e.get("a0")
        e["a1"] = e.get("a1")
        if e["type"] == "ARC" and e["pts"] and e["r"]:
            # un arc n'est pas un cercle : il a un debut et une fin, et c'est
            # ce qui permet de rechainer un contour dessine en CAO
            out.append({"type": "ARC", "c": e["pts"][0], "r": e["r"],
                        "a0": math.radians(e.get("a0", 0.0)),
                        "a1": math.radians(e.get("a1", 360.0)),
                        "layer": e["layer"]})
        elif e["type"] == "CIRCLE" and e["pts"] and e["r"]:
            out.append({"type": "CIRCLE", "c": e["pts"][0], "r": e["r"],
                        "layer": e["layer"]})
        elif e["type"] == "LINE" and len(e["pts"]) >= 2:
            out.append({"type": "LINE", "a": e["pts"][0], "b": e["pts"][1],
                        "layer": e["layer"]})
        elif e["type"] in ("POINT", "VERTEX") and e["pts"]:
            out.append({"type": "POINT", "c": e["pts"][0], "layer": e["layer"]})
        elif e["type"] in ("LWPOLYLINE", "POLYLINE") and len(e["pts"]) >= 2:
            out.append({"type": "POLY", "pts": e["pts"], "layer": e["layer"]})
    return out


def _extremites(e, n=48):
    """Points discretises d'un arc ou d'un segment, dans l'ordre du dessin."""
    if e["type"] == "LINE":
        return [e["a"], e["b"]]
    a0, a1 = e["a0"], e["a1"]
    while a1 <= a0:
        a1 += TAU
    k = max(2, int(n * (a1 - a0) / TAU) + 2)
    return [(e["c"][0] + e["r"] * math.cos(a0 + (a1 - a0) * i / (k - 1)),
             e["c"][1] + e["r"] * math.sin(a0 + (a1 - a0) * i / (k - 1)))
            for i in range(k)]


def assembler_contours(entities, tol=1e-3):
    """Rechaine les arcs et segments en contours, en gardant les entites.

    Une esquisse CAO n'exporte pas une polyligne : le profil y est une suite
    d'arcs et de segments, chacun une entite independante. On les rechaine
    par proximite de leurs extremites, et on **conserve l'ordre des entites
    sources** : c'est ce qui permet de relire ensuite la longueur du plan du
    marteau et le rayon de ses coins.
    """
    brins = []
    for e in entities:
        if e["type"] in ("ARC", "LINE"):
            pts = _extremites(e)
            if len(pts) >= 2:
                brins.append({"pts": list(pts), "ents": [e]})
    contours = []
    while brins:
        ch = brins.pop(0)
        avance = True
        while avance:
            avance = False
            for i, b in enumerate(brins):
                for sens in (1, -1):
                    bp = b["pts"][::sens]
                    if norm((ch["pts"][-1][0] - bp[0][0],
                             ch["pts"][-1][1] - bp[0][1])) < tol:
                        ch["pts"] += bp[1:]
                        ch["ents"] = ch["ents"] + b["ents"]
                        brins.pop(i)
                        avance = True
                        break
                    if norm((ch["pts"][0][0] - bp[-1][0],
                             ch["pts"][0][1] - bp[-1][1])) < tol:
                        ch["pts"] = bp[:-1] + ch["pts"]
                        ch["ents"] = b["ents"] + ch["ents"]
                        brins.pop(i)
                        avance = True
                        break
                if avance:
                    break
        ch["ferme"] = norm((ch["pts"][0][0] - ch["pts"][-1][0],
                            ch["pts"][0][1] - ch["pts"][-1][1])) < tol
        contours.append(ch)
    return contours


def _aire(pts):
    a = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % len(pts)]
        a += x1 * y2 - x2 * y1
    return 0.5 * abs(a)


def _dedans(pt, poly):
    x, y = pt
    c = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


def _rayons(pts, c=(0.0, 0.0)):
    """Rayons mini et maxi d'un contour autour de c, **segments compris**.

    Un plat n'est exporte que par ses deux extremites : prendre le mini sur
    les seuls sommets surestimerait le rayon interieur, alors que le point le
    plus proche de l'axe est au milieu du segment.
    """
    q0 = [(q[0] - c[0], q[1] - c[1]) for q in pts]
    rmax = max(norm(q) for q in q0)
    rmin = rmax
    n = len(q0)
    for i in range(n):
        a, b = q0[i], q0[(i + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        L2 = dx * dx + dy * dy
        t = 0.0 if L2 < 1e-18 else max(0.0, min(
            1.0, -(a[0] * dx + a[1] * dy) / L2))
        rmin = min(rmin, norm((a[0] + t * dx, a[1] + t * dy)))
    return rmin, rmax


def interpret(entities, base=None):
    """Deduit les parametres d'un dessin DXF.

    Le dessin est lu tel qu'un CAO l'exporte : le coeur et le marteau y sont
    des **suites d'arcs et de segments**, pas des polylignes, et les cercles
    pleins reperent les axes. On rechaine donc les contours, puis :

      - le **coeur** est le contour qui entoure un cercle plein en restant a
        distance de lui de tous les cotes (son axe de rotation) ;
      - le **marteau** est le contour qui contient l'autre cercle plein, le
        pivot etant une de ses pointes ;
      - le **plan de contact** est l'arete du marteau la plus proche du
        coeur, et sa longueur est reprise du dessin ;
      - tout est ensuite ramene dans le repere du coeur, axe a l'origine.

    A defaut d'arcs on retombe sur l'ancienne lecture : deux cercles
    concentriques pour les rayons, un cercle eloigne pour le pivot, une
    droite pour le plan.
    """
    p = base.copy() if base else Params()
    log = []
    out = {}
    cercles = [e for e in entities if e["type"] == "CIRCLE"]
    contours = [c for c in assembler_contours(entities)
                if c["ferme"] and len(c["pts"]) > 3]

    coeur = pivot = axe = None
    for c in contours:
        for ce in cercles:
            if not _dedans(ce["c"], c["pts"]):
                continue
            rmin, rmax = _rayons(c["pts"], ce["c"])
            if rmin > 0.08 * rmax:          # axe bien a l'interieur
                if coeur is None or _aire(c["pts"]) > _aire(coeur["pts"]):
                    coeur, axe = c, ce
    if coeur is None:
        # pas de cercle d'axe dans le dessin : l'origine en tient lieu,
        # c'est la convention du dessin
        for c in contours:
            if not _dedans((0.0, 0.0), c["pts"]):
                continue
            rmin, rmax = _rayons(c["pts"], (0.0, 0.0))
            if rmin > 0.08 * rmax:
                if coeur is None or _aire(c["pts"]) > _aire(coeur["pts"]):
                    coeur = c
        if coeur is not None:
            axe = {"c": (0.0, 0.0), "r": 0.0}
            log.append("aucun cercle d'axe : l'origine est prise comme "
                       "centre de rotation")

    marteau = None
    if coeur is not None:
        reste = [c for c in contours if c is not coeur]
        autres = [ce for ce in cercles if ce is not axe]
        if autres and reste:
            # le pivot est souvent une **pointe** du contour du marteau,
            # donc le test "point interieur" le rate : on prend le contour
            # dont un point est le plus proche de ce cercle
            pivot = max(autres, key=lambda ce: norm(
                (ce["c"][0] - axe["c"][0], ce["c"][1] - axe["c"][1])))
            marteau = min(reste, key=lambda c: min(
                norm((q[0] - pivot["c"][0], q[1] - pivot["c"][1]))
                for q in c["pts"]))
        elif reste:
            # pas de cercle de pivot non plus : c'est le point du marteau le
            # plus eloigne du coeur, la ou aboutissent ses bras
            marteau = max(reste, key=lambda c: max(
                norm((q[0] - axe["c"][0], q[1] - axe["c"][1]))
                for q in c["pts"]))
            pt = max(marteau["pts"], key=lambda q: norm(
                (q[0] - axe["c"][0], q[1] - axe["c"][1])))
            pivot = {"c": pt, "r": 0.0}
            log.append("aucun cercle de pivot : le point du marteau le plus "
                       "eloigne du coeur en tient lieu")

    if coeur is None:
        return _interpret_simple(entities, p)

    o = axe["c"]
    pts = [(q[0] - o[0], q[1] - o[1]) for q in coeur["pts"]]
    out["profil"] = pts
    rmin, rmax = _rayons(pts)
    out["r_min"], out["r_max"] = rmin, rmax
    n_ent = sum(1 for e in entities if e["type"] in ("ARC", "LINE"))
    log.append("coeur : contour ferme de %d points, rechaine depuis les arcs "
               "et segments du dessin (%d au total)" % (len(pts), n_ent))
    log.append("axe de rotation : cercle de %.3f mm en (%.3f, %.3f) ; tout "
               "est ramene la" % (2 * axe["r"], o[0], o[1]))
    log.append("rayons du coeur : %.3f a %.3f mm" % (rmin, rmax))

    if pivot is None or marteau is None:
        log.append("pivot ou contour du marteau non trouve : valeurs "
                   "courantes conservees")
        return out, log

    P = (pivot["c"][0] - o[0], pivot["c"][1] - o[1])
    out["pivot_x"], out["pivot_y"] = P
    log.append("pivot du marteau : (%.3f, %.3f), a %.3f mm de l'axe"
               % (P[0], P[1], norm(P)))

    # Plan de contact : relu sur les **entites** du marteau, pas sur les
    # points discretises (un arc y produirait de minuscules segments). On
    # ecarte les droites qui aboutissent au pivot - ce sont les bras - et on
    # retient la plus proche du coeur. Les arcs qui l'encadrent dans la
    # chaine donnent le rayon de ses coins.
    rp = max(pivot["r"] * 1.5, 0.03 * norm(P))
    ents = marteau["ents"]
    best = None
    for k, ent in enumerate(ents):
        if ent["type"] != "LINE":
            continue
        a = (ent["a"][0] - o[0], ent["a"][1] - o[1])
        b = (ent["b"][0] - o[0], ent["b"][1] - o[1])
        L = norm((b[0] - a[0], b[1] - a[1]))
        if L < 1e-6:
            continue
        if min(norm((a[0] - P[0], a[1] - P[1])),
               norm((b[0] - P[0], b[1] - P[1]))) < rp:
            continue
        m = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))
        d = min(norm((m[0] - q[0], m[1] - q[1])) for q in pts)
        if best is None or d < best[0]:
            best = (d, a, b, L, k)
    if best is None:
        log.append("plan de contact non identifie sur le marteau")
        return out, log
    _d, a, b, L, k = best
    rayons = [ents[jj]["r"] for jj in (k - 1, k + 1)
              if 0 <= jj < len(ents) and ents[jj]["type"] == "ARC"]
    if rayons:
        out["rayon_coin"] = sum(rayons) / len(rayons)
        log.append("coins du plan : %d arc(s), rayon moyen %.3f mm"
                   % (len(rayons), out["rayon_coin"]))
    dx, dy = b[0] - a[0], b[1] - a[1]
    nx, ny = -dy / L, dx / L
    if nx * a[0] + ny * a[1] < 0:
        nx, ny = -nx, -ny
    g = math.atan2(ny, nx)
    ang = math.degrees(g - math.atan2(P[1], P[0]))
    ang = (ang + 180.0) % 360.0 - 180.0
    if ang > 90.0:
        ang -= 180.0
    elif ang < -90.0:
        ang += 180.0
    out["angle_plan"] = ang
    out["longueur_plan"] = L
    log.append("plan de contact : segment de %.3f mm, a %.3f mm de l'axe"
               % (L, abs(nx * a[0] + ny * a[1])))
    log.append("pente du plan : %.2f deg, mesuree de O->P a la normale"
               % ang)
    _aligner_sur_le_plat(out, log)
    return out, log


def _interpret_simple(entities, p):
    """Lecture de secours : cercles concentriques et droites."""
    log = []
    out = {}
    circles = [e for e in entities if e["type"] == "CIRCLE"]
    lines = [e for e in entities if e["type"] == "LINE"]
    points = [e for e in entities if e["type"] == "POINT"]
    polys = [e for e in entities if e["type"] == "POLY" and len(e["pts"]) > 8]
    polys_tous = [e for e in entities if e["type"] == "POLY"]
    if polys:
        prof = max(polys, key=lambda e: len(e["pts"]))
        out["profil"] = prof["pts"]
        out["r_min"], out["r_max"] = _rayons(prof["pts"])
        log.append("profil du coeur : %d points, rayon %.3f a %.3f mm "
                   "(origine = axe de rotation)"
                   % (len(prof["pts"]), out["r_min"], out["r_max"]))
    tol = 1e-3
    if circles:
        tol = max(1e-3, 0.02 * max(c["r"] for c in circles))
    centered = [c for c in circles if norm(c["c"]) <= tol]
    far = [c for c in circles if norm(c["c"]) > tol]
    if len(centered) >= 2:
        rs = sorted(c["r"] for c in centered)
        out["r_min"], out["r_max"] = rs[0], rs[-1]
        log.append("r_min = %.3f et encombrement = %.3f mm (cercles centres)"
                   % (rs[0], rs[-1]))
    piv = None
    if far:
        piv = max(far, key=lambda c: norm(c["c"]))["c"]
    elif points:
        piv = max(points, key=lambda c: norm(c["c"]))["c"]
    if piv:
        out["pivot_x"], out["pivot_y"] = piv
        log.append("pivot marteau = (%.3f, %.3f)" % piv)
    px = out.get("pivot_x", p.pivot_x)
    py = out.get("pivot_y", p.pivot_y)

    def _len(l):
        return norm((l["b"][0] - l["a"][0], l["b"][1] - l["a"][1]))

    def _dist_o(l):
        dx, dy = l["b"][0] - l["a"][0], l["b"][1] - l["a"][1]
        L = math.hypot(dx, dy)
        return abs(dx * l["a"][1] - dy * l["a"][0]) / L if L > 1e-9 else 0.0

    cands = []
    for kw in ("PLAN_CONTACT", "PLAN", "FACE", "CONTACT", "HAMMER"):
        cands = [l for l in lines if kw in l["layer"].upper()]
        if cands:
            break
    if not cands:
        cands = [l for l in lines if _dist_o(l) > 1e-3 * max(1.0, _len(l))]
    if not cands:
        cands = lines
    # Le contour materiel du marteau, s'il est dessine comme une polyligne,
    # est plus fiable que la plus longue droite du calque : celle-ci est le
    # bras du marteau. On le reconnait a sa taille (une vingtaine de points)
    # et, pour nos propres exports, a son calque.
    cont = [e for e in polys_tous
            if 8 <= len(e["pts"]) <= 60
            and (e.get("layer") or "").upper().startswith("PLAN")]
    if not cont:
        cont = [e for e in polys_tous if 8 <= len(e["pts"]) <= 60]
    for e in sorted(cont, key=lambda e: len(e["pts"])):
        got = _plan_depuis_polyligne(e["pts"], (px, py))
        if got:
            out.update(got)
            log.append("plan de contact lu dans le contour du marteau : "
                       "partie droite %.4f mm, pente %.2f deg%s"
                       % (got["longueur_plan"], got["angle_plan"],
                          (", rayon de coin %.4f mm" % got["rayon_coin"])
                          if "rayon_coin" in got else ""))
            _aligner_sur_le_plat(out, log)
            if not log:
                log.append("aucune entite exploitable")
            return out, log
    if cands:
        ln = max(cands, key=_len)
        dx, dy = ln["b"][0] - ln["a"][0], ln["b"][1] - ln["a"][1]
        L = math.hypot(dx, dy)
        if L > 1e-9:
            nx, ny = -dy / L, dx / L
            d = nx * ln["a"][0] + ny * ln["a"][1]
            if d < 0:
                nx, ny, d = -nx, -ny, -d
            a = math.degrees(math.atan2(ny, nx) - math.atan2(py, px))
            a = (a + 180.0) % 360.0 - 180.0
            if a > 90.0:
                a -= 180.0
            elif a < -90.0:
                a += 180.0
            out["angle_plan"] = a
            out["longueur_plan"] = L
            log.append("pente du plan de contact = %.2f deg" % a)
    _aligner_sur_le_plat(out, log)
    if not log:
        log.append("aucune entite exploitable")
    return out, log


# ----------------------------------------------------------------- ecriture

try:
    import ezdxf as _ez
except ImportError:                        # pragma: no cover
    _ez = None

#: Couleurs ACI des calques produits. Un DXF sans table de calques est
#: partiellement relu par beaucoup de logiciels de CAO : les entites y
#: apparaissent sans forme exploitable.
CALQUES = {
    "COEUR": 3,                  # vert
    "PLAN_CONTACT": 2,           # jaune
    "MARTEAU_CINEMATIQUE": 8,    # gris
    "BRAS": 1,                   # rouge
    "PIVOT_MARTEAU": 6,          # magenta
    "CONSTRUCTION": 4,           # cyan
    "PLAT_ZERO": 5,              # bleu
    "POINTE": 1,
    "COEUR_ARCS": 5,       # bleu : le profil approche par des arcs
}


class Dessin:
    """Accumulateur d'entites, ecrit par ezdxf quand il est disponible.

    L'ecriture a la main produisait un DXF a section HEADER vide et **sans
    table de calques**. Les logiciels de CAO le relisent alors mal : il en
    sortait des points, des petits cercles et des droites, mais pas de coeur.
    ezdxf ecrit un R2010 complet - en-tete, tables, poignees - et les contours
    y sont de vraies polylignes fermees.

    Le repli en texte brut reste la, pour ne pas faire d'ezdxf une dependance
    dure du calcul.
    """

    def __init__(self):
        self.ents = []

    def poly(self, pts, layer, closed=True):
        pts = [(float(a), float(b)) for a, b in pts]
        if len(pts) >= 2:
            self.ents.append(("poly", layer, pts, bool(closed)))

    def circle(self, c, r, layer):
        if r > 1e-9:
            self.ents.append(("circle", layer,
                              (float(c[0]), float(c[1])), float(r)))

    def line(self, a, b, layer):
        self.ents.append(("line", layer, (float(a[0]), float(a[1])),
                          (float(b[0]), float(b[1]))))

    def arc(self, centre, r, a0, a1, layer):
        """Vrai arc de cercle. `a0`, `a1` en degres, sens trigonometrique."""
        if r > 1e-9 and abs(a1 - a0) > 1e-9:
            self.ents.append(("arc", layer,
                              (float(centre[0]), float(centre[1])),
                              float(r), float(a0), float(a1)))

    def _calques(self):
        return sorted({e[1] for e in self.ents})

    # -- sortie ----------------------------------------------------------
    def _par_ezdxf(self, path):
        doc = _ez.new("R2010", setup=True)
        doc.header["$INSUNITS"] = 4        # millimetres
        doc.header["$MEASUREMENT"] = 1     # metrique
        msp = doc.modelspace()
        for nom in self._calques():
            if nom not in doc.layers:
                doc.layers.add(nom, color=CALQUES.get(nom, 7))
        for e in self.ents:
            if e[0] == "poly":
                _k, layer, pts, closed = e
                msp.add_lwpolyline(pts, close=closed,
                                   dxfattribs={"layer": layer})
            elif e[0] == "circle":
                _k, layer, c, r = e
                msp.add_circle(c, r, dxfattribs={"layer": layer})
            elif e[0] == "arc":
                _k, layer, c, r, a0, a1 = e
                msp.add_arc(c, r, a0, a1, dxfattribs={"layer": layer})
            else:
                _k, layer, a, b = e
                msp.add_line(a, b, dxfattribs={"layer": layer})
        doc.saveas(path)

    def _par_texte(self, path):
        w = ["0\nSECTION\n2\nHEADER\n0\nENDSEC",
             "0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n70\n%d"
             % max(1, len(self._calques()))]
        for nom in self._calques():
            w.append("0\nLAYER\n2\n%s\n70\n0\n62\n%d\n6\nCONTINUOUS"
                     % (nom, CALQUES.get(nom, 7)))
        w.append("0\nENDTAB\n0\nENDSEC")
        w.append("0\nSECTION\n2\nENTITIES")
        for e in self.ents:
            if e[0] == "poly":
                w.append(_poly(e[2], e[1], e[3]))
            elif e[0] == "circle":
                w.append(_circle(e[2], e[3], e[1]))
            elif e[0] == "arc":
                w.append("0\nARC\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n0.0"
                         "\n40\n%.6f\n50\n%.6f\n51\n%.6f"
                         % (e[1], e[2][0], e[2][1], e[3], e[4], e[5]))
            else:
                w.append(_line(e[2], e[3], e[1]))
        w.append("0\nENDSEC\n0\nEOF")
        with open(path, "w") as fh:
            fh.write("\n".join(w) + "\n")

    def enregistrer(self, path):
        if _ez is not None:
            self._par_ezdxf(path)
        else:                              # pragma: no cover
            self._par_texte(path)
        return path


def _hdr(w):
    w.append("0\nSECTION\n2\nHEADER\n0\nENDSEC")
    w.append("0\nSECTION\n2\nENTITIES")


def _poly(pts, layer, closed=True):
    s = ["0\nLWPOLYLINE", "8\n%s" % layer, "100\nAcDbEntity", "100\nAcDbPolyline",
         "90\n%d" % len(pts), "70\n%d" % (1 if closed else 0)]
    for x, y in pts:
        s.append("10\n%.6f\n20\n%.6f" % (x, y))
    return "\n".join(s)


def _circle(c, r, layer):
    return ("0\nCIRCLE\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n0.0\n40\n%.6f"
            % (layer, c[0], c[1], r))


def _line(a, b, layer):
    return ("0\nLINE\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n0.0\n11\n%.6f\n21\n%.6f\n31\n0.0"
            % (layer, a[0], a[1], b[0], b[1]))


def write_dxf(path, p, res, index=0):
    """Ecrit le profil du coeur et la geometrie de construction."""
    from .core import pose
    d = Dessin()
    prof = res.closed_profile(p.n_lobes)
    if prof:
        d.poly(prof, "COEUR")
    d.circle((0.0, 0.0), p.r_min, "CONSTRUCTION")
    d.circle((0.0, 0.0), p.r_max, "CONSTRUCTION")
    d.circle(p.pivot, max(0.2, 0.05 * p.r_max), "PIVOT_MARTEAU")
    d.line((0.0, 0.0), p.pivot, "CONSTRUCTION")
    if res.ok:
        ps = pose(p, res, index)
        f = ps["face"]
        if len(f) >= 2:
            d.poly(f, "PLAN_CONTACT", closed=False)
        d.line(p.pivot, ps["contact"], "BRAS")
    return d.enregistrer(path)


def _rayon_3pts(a, b, c):
    """Rayon du cercle circonscrit, ou 0 si les points sont alignes."""
    d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1])
               + c[0] * (a[1] - b[1]))
    if abs(d) < 1e-12:
        return 0.0
    ax, ay = a[0] * a[0] + a[1] * a[1], 0.0
    ux = (ax * (b[1] - c[1])
          + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1])
          + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / d
    uy = (ax * (c[0] - b[0])
          + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0])
          + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / d
    return math.hypot(a[0] - ux, a[1] - uy)


def _plan_depuis_polyligne(pts, P):
    """Lit le plan d'impulsion dans le contour materiel du marteau.

    Ce contour - celui qu'ecrit l'export - est une polyligne : flanc, quart de
    cercle, partie droite, quart de cercle, flanc. Son **plus long segment**
    est la partie droite, et les trois points de chaque bout donnent le rayon
    de l'arrondi par cercle circonscrit.

    Sans cela, la lecture de secours prenait la plus longue droite du calque,
    qui est le **bras** du marteau : 6.44 mm au lieu de 1.21 mm sur
    `test.dxf`, et aucun rayon de coin.

    Renvoie `{}` si la polyligne ne ressemble pas a un tel contour.
    """
    n = len(pts)
    if n < 8:
        return {}
    k, L = -1, 0.0
    for i in range(n - 1):
        d = math.hypot(pts[i + 1][0] - pts[i][0], pts[i + 1][1] - pts[i][1])
        if d > L:
            k, L = i, d
    if k < 2 or k > n - 4 or L < 1e-9:
        return {}
    a, b = pts[k], pts[k + 1]
    dx, dy = b[0] - a[0], b[1] - a[1]
    nx, ny = -dy / L, dx / L
    if nx * a[0] + ny * a[1] < 0.0:
        nx, ny = -nx, -ny
    ang = math.degrees(math.atan2(ny, nx) - math.atan2(P[1], P[0]))
    ang = (ang + 180.0) % 360.0 - 180.0
    if ang > 90.0:
        ang -= 180.0
    elif ang < -90.0:
        ang += 180.0
    out = {"angle_plan": ang, "longueur_plan": L}
    rs = [r for r in (_rayon_3pts(pts[k - 2], pts[k - 1], pts[k]),
                      _rayon_3pts(pts[k + 1], pts[k + 2], pts[k + 3]))
          if 1e-4 < r < 5.0]
    if rs:
        out["rayon_coin"] = sum(rs) / len(rs)
    return out


def aligner_contour(pts, pivot, angle_plan, grille=1440):
    """Tourne le contour pour poser son plat sur le plan du marteau au repos.

    Le modele suppose que `theta = 0` **est** la position stable : marteau a
    `phi = 0`, pose sur le plat du zero. La normale de ce plat est la direction
    ou la fonction d'appui est minimale, c'est-a-dire exactement la reference
    que le moteur prend au repos ; on l'amene sur la normale du plan du
    marteau, `gamma0 = angle(O->P) + angle_plan`.

    A appliquer **a chaque fois que gamma0 change** - donc quand l'angle du
    plan ou le pivot bougent - sinon le recalage fait a l'import devient faux
    et le marteau recommence sa course deja ouvert.

    L'operation est idempotente : la reappliquer ne derive pas.

    Renvoie `(points, angle_applique_en_degres)`.
    """
    if not pts or abs(pivot[0]) + abs(pivot[1]) < 1e-9:
        return list(pts), 0.0
    g0 = math.atan2(pivot[1], pivot[0]) + math.radians(angle_plan)
    best, bh = 0.0, None
    for i in range(grille):
        b = TAU * i / grille
        cb, sb = math.cos(b), math.sin(b)
        h = max(cb * q[0] + sb * q[1] for q in pts)
        if bh is None or h < bh:
            bh, best = h, b
    a = (g0 - best + math.pi) % TAU - math.pi
    if abs(a) < 1e-6:
        return list(pts), 0.0
    ca, sa = math.cos(a), math.sin(a)
    return ([(q[0] * ca - q[1] * sa, q[0] * sa + q[1] * ca) for q in pts],
            math.degrees(a))


def _aligner_sur_le_plat(out, log):
    """Applique `aligner_contour` au resultat d'une lecture DXF."""
    pts = out.get("profil")
    if not pts or "angle_plan" not in out:
        return
    P = (out.get("pivot_x", 0.0), out.get("pivot_y", 0.0))
    npts, ang = aligner_contour(pts, P, out["angle_plan"])
    if ang:
        out["profil"] = npts
        log.append("coeur tourne de %.2f deg pour poser son plat sur le plan "
                   "du marteau : theta = 0 est la position de repos" % ang)


def _dans_repere_coeur(pts, theta):
    """Ramene une geometrie du repere fixe vers celui du coeur.

    Le coeur est dessine dans **son** repere (position d'usinage). Le marteau,
    lui, est a la position correspondant a l'angle courant : il faut donc le
    tourner de -theta, sinon le plan de contact traverse le coeur sur le plan.
    Dans ce repere le pivot du marteau decrit un cercle autour de O, ce qui
    est bien le mouvement relatif du marteau vu du coeur.
    """
    c, s_ = math.cos(-theta), math.sin(-theta)
    return [(q[0] * c - q[1] * s_, q[0] * s_ + q[1] * c) for q in pts]


def write_dxf_chrono(path, p, r, x_deg=0.0, n_cine=0):
    """Ecrit le profil du coeur de chronographe et sa geometrie de construction.

    `n_cine` : nombre de positions du marteau a tracer en plus, pour voir la
    cinematique relative sur un tour.
    """
    from .chrono import pose_at
    d = Dessin()
    if r.profile:
        d.poly(r.profile, "COEUR")
    if r.flat:
        d.line(r.flat[0], r.flat[1], "PLAT_ZERO")
    if getattr(r, "horns", None) and len(r.horns) == 2:
        d.line(r.horns[0], r.horns[1], "POINTE")
    if getattr(r, "notch", None):
        d.poly(r.notch, "POINTE", closed=False)
    d.circle((0.0, 0.0), p.r_min, "CONSTRUCTION")
    d.circle((0.0, 0.0), p.r_max, "CONSTRUCTION")
    if r.ok:
        ps0 = pose_at(p, r, x_deg)
        piv0 = _dans_repere_coeur([p.pivot], ps0["theta"])[0]
        d.circle(piv0, max(0.2, 0.05 * p.r_max), "PIVOT_MARTEAU")
        d.line((0.0, 0.0), piv0, "CONSTRUCTION")

        def ecrire(x, calque_plan, calque_bras, calque_pivot):
            ps = pose_at(p, r, x)
            th = ps["theta"]
            f = _dans_repere_coeur(ps["face"], th)
            piv = _dans_repere_coeur([p.pivot], th)[0]
            ct = _dans_repere_coeur([ps["contact"]], th)[0]
            d.poly(f, calque_plan, closed=False)
            d.line(piv, ct, calque_bras)
            if calque_pivot:
                d.circle(piv, max(0.2, 0.05 * p.r_max), calque_pivot)

        ecrire(x_deg, "PLAN_CONTACT", "BRAS", None)
        for i in range(max(0, n_cine)):
            ecrire(-r.course_b + (r.course_a + r.course_b) * i
                   / max(1, n_cine - 1), "MARTEAU_CINEMATIQUE",
                   "MARTEAU_CINEMATIQUE", "MARTEAU_CINEMATIQUE")
    return d.enregistrer(path)
