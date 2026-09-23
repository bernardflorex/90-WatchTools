# dxf_io.py
"""
Import/export DXF, unique pour tout le programme.

Import : reconstruit la fibre neutre + l'épaisseur locale d'un spiral à partir
d'un DXF quelconque (dessin machine avec flancs séparés, ou contour fermé
unique), par ray-casting orthogonal à la normale locale d'un flanc de
référence. Portage de la logique de spiran_v1.py, avec deux corrections :

  1. `ezdxf.path.flattening` n'existe pas en tant que fonction de module —
     c'est une méthode d'instance de `Path` (`chemin.flattening(distance=...)`).
     L'ancien code appelait `path.flattening(p, distance=...)` sans même avoir
     importé `path` depuis `ezdxf` (seul `matplotlib.pyplot` était importé) :
     le bloc `try/except` avalait silencieusement une `NameError`, et les
     entités ARC (non gérées par le fallback manuel) disparaissaient de
     l'extraction sans avertissement. Vérifié sur spimo02.dxf : 2 ARC sur
     9 entités auraient été perdus.
  2. Le fallback manuel pour SPLINE/LWPOLYLINE reste en secours si le
     flattening ezdxf échoue pour une raison quelconque, mais n'est plus le
     chemin normal.

Export : régénère les deux flancs (offset ±epaisseur/2 le long de la normale
locale) en LWPOLYLINE, pour (a) livrer une géométrie optimisée à la DRIE et
(b) boucler la validation (réimporter ce qu'on vient d'exporter et vérifier
qu'on retrouve la même physique).
"""
import numpy as np
import ezdxf
from ezdxf import path as ezpath
from scipy.ndimage import gaussian_filter1d

from .modele import GeometrieSpiral

TYPES_SUPPORTES = ("LWPOLYLINE", "POLYLINE", "SPLINE", "ARC", "LINE", "CIRCLE", "ELLIPSE")


def _recuperer_toutes_entites(layout):
    entites = []
    for entity in layout:
        if entity.dxftype() == "INSERT":
            try:
                entites.extend(_recuperer_toutes_entites(entity.virtual_entities()))
            except Exception:
                pass
        else:
            entites.append(entity)
    return entites


def _flatten_entite(entity, distance=0.0005):
    """Retourne les points (N,2) d'une entité, quelle que soit sa nature
    géométrique (arc, spline, polyligne à bulges, ...). Chemin normal :
    ezdxf.path (gère tout). Fallback : extraction manuelle des points bruts."""
    dxftype = entity.dxftype()
    try:
        chemin = ezpath.make_path(entity)
        pts = [(v.x, v.y) for v in chemin.flattening(distance=distance)]
        if len(pts) > 1:
            return np.array(pts)
    except Exception:
        pass

    try:
        if dxftype in ("LWPOLYLINE", "POLYLINE"):
            pts = [(pt[0], pt[1]) for pt in entity.get_points()]
        elif dxftype == "SPLINE":
            fit = getattr(entity, "fit_points", [])
            ctrl = getattr(entity, "control_points", [])
            src = fit if len(fit) > 0 else ctrl
            pts = [(pt[0], pt[1]) for pt in src]
        elif dxftype == "LINE":
            pts = [
                (entity.dxf.start.x, entity.dxf.start.y),
                (entity.dxf.end.x, entity.dxf.end.y),
            ]
        else:
            pts = []
    except Exception:
        pts = []

    return np.array(pts) if len(pts) > 1 else np.empty((0, 2))


def _chainer_segments(raw_segments, tol_jonction):
    """Chaîne bout-à-bout des segments bruts (dans n'importe quel ordre ou
    sens) en chaînes continues, en recollant les extrémités les plus proches
    en dessous de tol_jonction."""
    pool = list(raw_segments)
    chains = []
    while pool:
        curr = pool.pop(0)
        changed = True
        while changed and pool:
            changed = False
            best_i, best_mode, min_d = -1, None, float("inf")
            for i, cand in enumerate(pool):
                d_es = np.linalg.norm(curr[-1] - cand[0])
                d_ee = np.linalg.norm(curr[-1] - cand[-1])
                d_se = np.linalg.norm(curr[0] - cand[-1])
                d_ss = np.linalg.norm(curr[0] - cand[0])
                d_min = min(d_es, d_ee, d_se, d_ss)
                if d_min < min_d:
                    min_d = d_min
                    best_i = i
                    best_mode = {d_es: "end-start", d_ee: "end-end",
                                 d_se: "start-end", d_ss: "start-start"}[d_min]
            if min_d < tol_jonction and best_i != -1:
                cand = pool.pop(best_i)
                if best_mode == "end-start":
                    curr = np.vstack((curr, cand[1:]))
                elif best_mode == "end-end":
                    curr = np.vstack((curr, cand[:-1][::-1]))
                elif best_mode == "start-end":
                    curr = np.vstack((cand[:-1], curr))
                elif best_mode == "start-start":
                    curr = np.vstack((cand[1:][::-1], curr))
                changed = True
        chains.append(curr)
    return chains


def _separer_flancs(chains):
    """Cas 1 : un contour fermé unique -> détection des deux points de
    pliure (bouchons virole/piton) par EXTREMA DE RAYON. Cas 2 : plusieurs
    chaînes ouvertes -> les deux plus longues sont les deux flancs.

    NOTE (remplace la détection par pic de courbure de spiran_v1.py) : la
    détection par courbure locale se laisse tromper par de petits éléments
    accessoires du dessin (pattes de fixation, encoches) dont la courbure
    ponctuelle dépasse celle, très douce, du vrai pli virole/piton — vérifié
    sur spimo02.dxf, qui donnait un rayon max reconstruit de 1.95mm au lieu
    de 2.67mm. Les deux plis réels sont là où la fibre neutre du ruban
    rebrousse chemin radialement, donc là où r(index) le long du contour
    fermé passe par un minimum (virole) et un maximum (piton) global —
    hypothèse fiable tant que les accessoires du dessin restent dans
    l'enveloppe [r_min, r_max] du spiral (nécessairement vrai : au-delà, ils
    entreraient en collision avec le reste du mouvement)."""
    if len(chains) == 1:
        pts = chains[0]
        dists = np.hypot(np.diff(pts[:, 0], append=pts[0, 0]),
                          np.diff(pts[:, 1], append=pts[0, 1]))
        pts = pts[dists > 1e-7]

        r = np.hypot(pts[:, 0], pts[:, 1])
        r_smooth = gaussian_filter1d(r, sigma=5, mode="wrap")
        i1 = int(np.argmin(r_smooth))
        i2 = int(np.argmax(r_smooth))

        idx_start, idx_end = sorted([i1, i2])
        bord_A = pts[idx_start:idx_end + 1]
        bord_B = np.vstack((pts[idx_end:], pts[:idx_start + 1]))
    else:
        chains_triees = sorted(
            chains, key=lambda c: np.sum(np.hypot(*np.diff(c, axis=0).T)), reverse=True
        )
        bord_A, bord_B = chains_triees[0], chains_triees[1]
    return bord_A, bord_B


def charger_dxf(chemin_dxf: str, n_elements: int = 2000, hauteur_h: float = 0.100,
                 scale: float = 1.0, epaisseur_max_attendue: float = 0.200,
                 epaisseur_fallback: float = 0.035) -> GeometrieSpiral:
    """Charge un DXF de spiral quelconque et retourne une GeometrieSpiral
    (index 0 = virole/centre, index -1 = piton/extérieur), à n_elements
    points uniformément répartis en abscisse curviligne.

    hauteur_h (mm) n'est pas encodée dans un DXF 2D : elle doit être fournie
    (100 µm par défaut, conforme au cahier des charges DRIE du projet).
    """
    doc = ezdxf.readfile(chemin_dxf)
    msp = doc.modelspace()

    raw_segments = []
    for entity in _recuperer_toutes_entites(msp):
        if entity.dxftype() not in TYPES_SUPPORTES:
            continue
        pts = _flatten_entite(entity) * scale
        if len(pts) > 1:
            raw_segments.append(pts)

    if not raw_segments:
        raise ValueError(f"Aucune courbe valide trouvée dans {chemin_dxf}.")

    tol_jonction = 0.1 * scale
    chains = _chainer_segments(raw_segments, tol_jonction)

    longueurs = [np.sum(np.hypot(*np.diff(c, axis=0).T)) for c in chains]
    max_len = max(longueurs)
    chains = [c for c, l in zip(chains, longueurs) if l > 0.20 * max_len]

    bord_A, bord_B = _separer_flancs(chains)

    cx, cy = np.mean(bord_A[:, 0]), np.mean(bord_A[:, 1])
    if np.hypot(bord_A[0, 0] - cx, bord_A[0, 1] - cy) > np.hypot(bord_A[-1, 0] - cx, bord_A[-1, 1] - cy):
        bord_A = bord_A[::-1]
    if np.hypot(bord_B[0, 0] - cx, bord_B[0, 1] - cy) > np.hypot(bord_B[-1, 0] - cx, bord_B[-1, 1] - cy):
        bord_B = bord_B[::-1]

    # --- Correspondance des deux flancs par ANGLE POLAIRE (theta), pas par
    # ré-échantillonnage indépendant en abscisse curviligne ---
    # L'ancienne méthode (ray-casting depuis un flanc de référence
    # ré-échantillonné par fraction d'arc-length) suppose implicitement que
    # "même fraction de longueur" ~= "même position radiale sur le ruban".
    # C'est FAUX dès que le spiral fait beaucoup de tours avec un grand
    # rapport r_max/r_min : la longueur d'arc par radian (ds ~ r.dtheta)
    # diffère entre les deux flancs (rayons légèrement différents), donc
    # leurs paramétrages en fraction d'arc-length dérivent l'un par rapport
    # à l'autre au fil des tours. Vérifié sur spimo02.dxf (~14 tours,
    # r: 0.48->2.67mm) : jusqu'à 1.4 rad (~80°) de décalage angulaire accumulé
    # au piton, faisant passer l'épaisseur "mesurée" de ~35 à ~155 µm de façon
    # purement artefactuelle (diagnostic 2026-09-24). theta est un paramètre
    # commun aux deux flancs PAR CONSTRUCTION (même angle polaire = même
    # rayon radial du ruban) : plus de dérive possible.
    theta_A = np.unwrap(np.arctan2(bord_A[:, 1], bord_A[:, 0]))
    r_A = np.hypot(bord_A[:, 0], bord_A[:, 1])
    theta_B = np.unwrap(np.arctan2(bord_B[:, 1], bord_B[:, 0]))
    r_B = np.hypot(bord_B[:, 0], bord_B[:, 1])

    def _theta_r_triable(theta, r):
        idx = np.argsort(theta, kind="stable")
        theta_s, r_s = theta[idx], r[idx]
        _, u = np.unique(theta_s, return_index=True)
        u = np.sort(u)
        return theta_s[u], r_s[u]

    theta_A, r_A = _theta_r_triable(theta_A, r_A)
    theta_B, r_B = _theta_r_triable(theta_B, r_B)

    theta_min = max(theta_A[0], theta_B[0])
    theta_max = min(theta_A[-1], theta_B[-1])
    if theta_max <= theta_min:
        raise ValueError(
            f"Les deux flancs de {chemin_dxf} ne se recouvrent pas angulairement "
            f"(virole/piton mal détectés ?)."
        )
    theta_final = np.linspace(theta_min, theta_max, n_elements)

    r_A_i = np.interp(theta_final, theta_A, r_A)
    r_B_i = np.interp(theta_final, theta_B, r_B)
    epaisseurs = np.abs(r_B_i - r_A_i)
    r_milieu = (r_A_i + r_B_i) / 2.0
    x_ref = r_milieu * np.cos(theta_final)
    y_ref = r_milieu * np.sin(theta_final)

    epaisseurs = np.nan_to_num(epaisseurs, nan=epaisseur_fallback)

    # --- Troncature des pattes de fixation rigides (collet/piton) ---
    # Une patte de fixation n'est pas une section de lame active : le
    # ray-casting y détecte des "épaisseurs" très supérieures à l'épaisseur
    # nominale de la lame (le contour n'est plus un ruban simple à cet
    # endroit). Les inclure fausserait la masse et la rigidité de
    # l'elastica. Détection : depuis chaque extrémité, on avance tant que
    # l'épaisseur dépasse `facteur_stub` x la médiane du corps central
    # (90% du profil, robuste aux extrémités déjà suspectes), en exigeant
    # `k_consecutifs` valeurs redevenues normales pour confirmer qu'on est
    # sorti de la patte (une patte peut présenter des creux locaux, d'où le
    # besoin d'une fenêtre plutôt qu'un simple seuil ponctuel). Recherche
    # bornée à 40% de la longueur de chaque côté par sécurité : au-delà, on
    # renonce plutôt que de risquer de tronquer une vraie portion de lame.
    ref_e = np.median(epaisseurs[int(0.1 * n_elements):int(0.9 * n_elements)])
    seuil_stub = 2.5 * ref_e
    k_consecutifs = 15
    limite_max = n_elements * 4 // 10

    i_debut = 0
    consecutifs = 0
    for i in range(limite_max):
        if epaisseurs[i] <= seuil_stub:
            consecutifs += 1
            if consecutifs >= k_consecutifs:
                i_debut = i - k_consecutifs + 1
                break
        else:
            consecutifs = 0

    i_fin = n_elements - 1
    consecutifs = 0
    for i in range(n_elements - 1, n_elements - 1 - limite_max, -1):
        if epaisseurs[i] <= seuil_stub:
            consecutifs += 1
            if consecutifs >= k_consecutifs:
                i_fin = i + k_consecutifs - 1
                break
        else:
            consecutifs = 0

    if i_debut > 0 or i_fin < n_elements - 1:
        x_ref = x_ref[i_debut:i_fin + 1]
        y_ref = y_ref[i_debut:i_fin + 1]
        epaisseurs = epaisseurs[i_debut:i_fin + 1]
        n_elements = len(epaisseurs)

    n_edge = max(5, int(n_elements * 0.01))
    epaisseurs[:n_edge] = epaisseurs[n_edge]
    epaisseurs[-n_edge:] = epaisseurs[-n_edge - 1]
    epaisseurs = gaussian_filter1d(epaisseurs, sigma=3)
    epaisseurs = np.clip(epaisseurs, a_min=0.008, a_max=epaisseur_max_attendue)

    Xm, Ym = x_ref, y_ref

    if np.hypot(Xm[0], Ym[0]) > np.hypot(Xm[-1], Ym[-1]):
        Xm, Ym, epaisseurs = Xm[::-1], Ym[::-1], epaisseurs[::-1]

    return GeometrieSpiral(x=Xm, y=Ym, epaisseur=epaisseurs, hauteur=hauteur_h)


def exporter_dxf(geom: GeometrieSpiral, chemin_dxf: str):
    """Exporte les deux flancs (offset ±epaisseur/2) en LWPOLYLINE, pour la
    fabrication DRIE et pour la boucle de validation (charger_dxf(exporter_dxf(...))
    doit redonner la même physique à la précision du maillage près)."""
    ext, inte = geom.contour_flancs()
    doc = ezdxf.new(setup=False)
    doc.header["$INSUNITS"] = 4  # millimètres
    msp = doc.modelspace()
    msp.add_lwpolyline(ext, format="xy", dxfattribs={"layer": "spiral_flanc_ext"})
    msp.add_lwpolyline(inte, format="xy", dxfattribs={"layer": "spiral_flanc_int"})
    doc.saveas(chemin_dxf)
