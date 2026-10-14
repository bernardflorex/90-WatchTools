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
    sens) en chaînes continues, en recollant les extrémités distantes de moins
    de tol_jonction.

    La recherche du voisin le plus proche passe par un ARBRE K-D construit une
    fois sur toutes les extrémités, et non par un balayage du pool restant.
    La version précédente rescannait tout le pool à chaque rattachement, soit
    un coût cubique en nombre de segments : invisible sur un DXF fait de
    quelques LWPOLYLINE (spimo02.dxf : 9 entités, instantané), rédhibitoire
    sur un DXF exporté en segments élémentaires. Mesuré sur un fichier réel de
    6 087 LINE, l'import passait plusieurs minutes dans cette seule fonction —
    l'application paraissait figée (diagnostic 2026-10-07).

    L'arbre ne sait pas supprimer de point : les segments déjà consommés sont
    marqués et simplement ignorés à l'interrogation. C'est sans conséquence
    tant que la tolérance reste petite devant l'espacement des segments, ce
    qui est le cas par construction (tol_jonction est une tolérance de
    recollement, pas un rayon de recherche).
    """
    from scipy.spatial import cKDTree

    segments = [np.asarray(seg, dtype=float) for seg in raw_segments]
    n = len(segments)
    if n == 0:
        return []

    # Index 2i = début du segment i, 2i+1 = sa fin.
    extremites = np.array([[seg[0], seg[-1]] for seg in segments],
                           dtype=float).reshape(2 * n, 2)
    arbre = cKDTree(extremites)
    utilise = np.zeros(n, dtype=bool)
    chaines = []

    for depart in range(n):
        if utilise[depart]:
            continue
        utilise[depart] = True
        courant = segments[depart]
        # On étend d'abord par la fin, puis par le début : une chaîne est un
        # chemin, la parcourir dans les deux sens depuis la graine suffit.
        for par_la_fin in (True, False):
            while True:
                point = courant[-1] if par_la_fin else courant[0]
                # Tangente sortante de la chaîne, pour départager les candidats
                # par CONTINUITÉ DE DIRECTION et non par la seule distance.
                voisin = courant[-2] if par_la_fin else courant[1]
                tangente = point - voisin
                norme = float(np.hypot(*tangente))
                tangente = tangente / norme if norme > 1e-12 else None

                alignes, tous = [], []
                for j in arbre.query_ball_point(point, tol_jonction):
                    i_c = j // 2
                    if utilise[i_c]:
                        continue
                    d = float(np.hypot(*(extremites[j] - point)))
                    tous.append((d, j))
                    if tangente is None:
                        continue
                    # Direction dans laquelle part le segment candidat, depuis
                    # le bout apparié.
                    cand = segments[i_c]
                    sortie = (cand[1] - cand[0]) if (j % 2 == 0) else (cand[-2] - cand[-1])
                    n_s = float(np.hypot(*sortie))
                    if n_s <= 1e-12:
                        continue
                    if float(np.dot(tangente, sortie / n_s)) > 0.3:
                        alignes.append((d, j))
                # Un dessin fragmenté présente à chaque coupure DEUX candidats à
                # la même distance : la vraie suite du flanc, alignée à un ou
                # deux degrés, et le flanc d'en face, à 90 ou 180°. Mesuré sur un
                # export réel (spiral-111.dxf) : trous de 38 µm le long du flanc
                # pour 38 µm d'écart entre flancs — la distance seule ne peut pas
                # trancher, et une tolérance plus large ne fait qu'empirer les
                # choses en facilitant le saut d'un flanc à l'autre.
                # On privilégie donc les candidats alignés ; s'il n'y en a aucun,
                # on retombe sur le plus proche, ce qui laisse passer le demi-tour
                # légitime au bout du ruban (diagnostic 2026-10-08).
                choix = alignes or tous
                if not choix:
                    break
                d_min, meilleur = min(choix)
                i_seg, cote = meilleur // 2, meilleur % 2
                cand = segments[i_seg]
                utilise[i_seg] = True
                if par_la_fin:
                    # Le bout de `cand` qui a été apparié doit devenir la
                    # jonction, donc venir en tête du morceau ajouté.
                    morceau = cand if cote == 0 else cand[::-1]
                    courant = np.vstack((courant, morceau[1:]))
                else:
                    morceau = cand[::-1] if cote == 0 else cand
                    courant = np.vstack((morceau[:-1], courant))
        chaines.append(courant)
    return chaines


def _densifier_flanc(flanc, pas_max):
    """Insère des points LE LONG des segments d'un flanc, sans en changer la
    géométrie : les points ajoutés sont sur la droite qui joint déjà deux
    points voisins.

    Nécessaire avant l'appariement en angle polaire. Celui-ci interpole r(θ)
    linéairement entre deux points consécutifs, or r(θ) n'est PAS linéaire le
    long d'une corde : il passe par un minimum au pied de la perpendiculaire
    issue du centre. Sur une corde longue, l'interpolation surestime donc r, et
    l'épaisseur — différence de deux r interpolés — devient la différence de
    deux erreurs, c'est-à-dire du bruit.

    Constaté sur un DXF réel (Spiral-K2-2025) : la courbe terminale y comporte
    une PORTION RECTILIGNE d'environ 1,5 mm, dessinée en deux cordes de 948 et
    593 µm là où le reste du contour est échantillonné tous les 41 µm.
    L'épaisseur reconstruite y tombait à 10 µm au lieu de 36, sur 14 éléments
    (diagnostic 2026-10-07). Densifier supprime l'artefact sans rien inventer :
    la flèche d'une corde de 50 µm à r = 2,4 mm vaut 0,13 µm, négligeable
    devant l'épaisseur de la lame.
    """
    seg = np.hypot(*np.diff(flanc, axis=0).T)
    if len(seg) == 0 or seg.max() <= pas_max:
        return flanc
    morceaux = [flanc[:1]]
    for i, longueur in enumerate(seg):
        n_sous = max(1, int(np.ceil(longueur / pas_max)))
        # On exclut le point de départ (déjà posé) et on garde le point
        # d'arrivée : les morceaux se recollent sans doublon.
        t = np.linspace(0.0, 1.0, n_sous + 1)[1:, None]
        morceaux.append(flanc[i] + t * (flanc[i + 1] - flanc[i]))
    return np.vstack(morceaux)


def _separer_flancs(chains):
    """Cas 1 : un contour unique, qui parcourt un flanc puis l'autre -> on le
    scinde aux deux PLIS (virole et piton). Cas 2 : plusieurs chaînes ouvertes
    -> les deux plus longues sont les deux flancs.

    Les plis sont cherchés aux extrema de l'ANGLE POLAIRE DÉROULÉ, pas du
    rayon. C'est la correction du 2026-10-07, et elle n'est pas cosmétique :

      - θ croît de façon strictement monotone le long d'un flanc et décroît le
        long de l'autre. Son maximum est donc AU bout du ruban, par
        construction, quelle que soit la forme de la terminaison.
      - r, lui, est maximal sur le flanc EXTÉRIEUR, qui est plus loin du centre
        que l'intérieur. En contournant la pointe, le rayon continue donc de
        croître un peu. Mesuré sur un DXF réel (Spiral-K2-2025) : le maximum de
        rayon tombe 127 points et 102° APRÈS la pointe, déjà sur le flanc de
        retour. La découpe était donc décalée d'autant, un flanc empiétait sur
        l'autre, et l'intersection angulaire des deux (voir charger_dxf)
        amputait le spiral de 102° au piton -- silencieusement.

    Une détection par pic de COURBURE avait déjà été écartée auparavant : elle
    se laissait tromper par les petits éléments accessoires du dessin. Le
    critère angulaire ne souffre ni de l'un ni de l'autre défaut.
    """
    if len(chains) == 1:
        pts = chains[0]
        dists = np.hypot(np.diff(pts[:, 0], append=pts[0, 0]),
                          np.diff(pts[:, 1], append=pts[0, 1]))
        pts = pts[dists > 1e-7]

        # Les points très proches du centre ont un angle polaire dénué de sens
        # (bruit pur) et créeraient de faux extrema. Ce sont des restes de
        # lignes de construction, jamais de la matière : le collet d'un spiral
        # est à un rayon du même ordre que le reste de la pièce.
        r = np.hypot(pts[:, 0], pts[:, 1])
        valides = r > 0.05 * np.median(r)
        if valides.sum() < 8:
            raise ValueError("contour trop pauvre pour être scindé en deux flancs")
        pts, r = pts[valides], r[valides]

        theta = np.unwrap(np.arctan2(pts[:, 1], pts[:, 0]))
        if theta.max() - theta.min() > np.pi:
            i1, i2 = int(np.argmin(theta)), int(np.argmax(theta))
        else:
            # Pas un spiral enroulé (moins d'un demi-tour) : l'angle ne
            # discrimine plus, on retombe sur les extrema de rayon lissés.
            r_lisse = gaussian_filter1d(r, sigma=5, mode="nearest")
            i1, i2 = int(np.argmin(r_lisse)), int(np.argmax(r_lisse))

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
                 epaisseur_fallback: float = 0.035,
                 plan_cristallin: str = "111",
                 orientation_cristal_deg: float = 0.0) -> GeometrieSpiral:
    """Charge un DXF de spiral quelconque et retourne une GeometrieSpiral
    (index 0 = virole/centre, index -1 = piton/extérieur), à n_elements
    points uniformément répartis en abscisse curviligne.

    hauteur_h (mm) n'est pas encodée dans un DXF 2D : elle doit être fournie
    (100 µm par défaut, conforme au cahier des charges DRIE du projet). Le plan
    cristallin du wafer non plus, et il pèse jusqu'à 13 % sur la raideur : voir
    constantes.module_young_silicium.
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

    # --- Tolérance de recollement ADAPTATIVE ---
    # Une tolérance fixe ne convient pas : trop serrée, elle laisse en morceaux
    # un dessin dont les extrémités ne coïncident qu'approximativement ; trop
    # lâche, elle recolle au spiral des entités qui n'en font pas partie.
    # La valeur de 0.1 mm employée jusqu'ici était presque trois fois
    # l'épaisseur de la lame : sur un DXF réel (Spiral-K2-2025, 6 087 LINE),
    # elle absorbait une ligne de construction partant du centre et ramenait le
    # rayon de collet reconstruit à 0.249 mm au lieu de 0.468 -- silencieusement
    # (diagnostic 2026-10-07).
    # On essaie donc du plus serré au plus lâche, et on retient la PREMIÈRE
    # tolérance qui produise deux flancs plausibles. Le critère de plausibilité
    # est que les deux flancs d'un ruban ont forcément des longueurs voisines :
    # ils ne diffèrent que par l'épaisseur intégrée sur la rotation totale,
    # soit quelques pour cent.
    essais = []
    bord_A = bord_B = None
    for facteur in (0.002, 0.01, 0.05, 0.1, 0.25, 0.5):
        tol_jonction = facteur * scale
        chains = _chainer_segments(raw_segments, tol_jonction)
        longueurs = [np.sum(np.hypot(*np.diff(c, axis=0).T)) for c in chains]
        gardees = [c for c, l in zip(chains, longueurs) if l > 0.20 * max(longueurs)]
        try:
            A, B = _separer_flancs(gardees)
        except (ValueError, IndexError):
            essais.append((tol_jonction, len(chains), None, 0.0))
            continue
        LA = float(np.sum(np.hypot(*np.diff(A, axis=0).T)))
        LB = float(np.sum(np.hypot(*np.diff(B, axis=0).T)))
        # COUVERTURE : part de la matière dessinée que représentent les deux
        # flancs RETENUS. Un dessin fragmenté donne des dizaines de morceaux de
        # longueurs voisines ; on en garderait alors deux au hasard, et la
        # géométrie reconstituée ne serait qu'un bout du spiral. Mesurer la part
        # des seules chaînes « gardées » ne suffit pas : quand tous les morceaux
        # se valent, ils passent tous le filtre des 20 % et la couverture paraît
        # parfaite. Sans ce contrôle, un fichier cassé en 151 morceaux
        # s'importait « avec succès » en un spiral de 33 mm au lieu de 120
        # (diagnostic 2026-10-08).
        couverture = float((LA + LB) / max(sum(longueurs), 1e-12))
        essais.append((tol_jonction, len(chains), (LA, LB), couverture))
        if (min(LA, LB) > 0 and max(LA, LB) / min(LA, LB) < 1.5
                and min(len(A), len(B)) >= 8 and couverture > 0.80):
            bord_A, bord_B = A, B
            break
    if bord_A is None:
        detail = " ; ".join(
            f"tol {t*1000:.1f} µm : {n} chaînes, {100*cv:.0f} % de la matière retenue"
            + (f", flancs {f[0]:.1f}/{f[1]:.1f} mm" if f else ", séparation impossible")
            for t, n, f, cv in essais)
        raise ValueError(
            f"Impossible de reconstituer deux flancs de longueurs comparables dans "
            f"{chemin_dxf}. Essais : {detail}. Le dessin contient probablement des "
            f"entités étrangères au spiral, ou les deux flancs n'y sont pas tous les deux."
        )

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
    # Densification avant l'appariement angulaire : le pas visé est la longueur
    # MÉDIANE des segments du dessin, bornée pour ne pas exploser le nombre de
    # points sur un dessin déjà fin ni rester trop grossier sur un dessin
    # pauvre. Seuls les rares segments plus longs que ce pas sont subdivisés.
    pas_cible = float(np.median(np.concatenate([
        np.hypot(*np.diff(bord_A, axis=0).T), np.hypot(*np.diff(bord_B, axis=0).T)])))
    pas_cible = float(np.clip(pas_cible, 0.005 * scale, 0.050 * scale))
    bord_A = _densifier_flanc(bord_A, pas_cible)
    bord_B = _densifier_flanc(bord_B, pas_cible)

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

    return GeometrieSpiral(x=Xm, y=Ym, epaisseur=epaisseurs, hauteur=hauteur_h,
                            plan_cristallin=plan_cristallin,
                            orientation_cristal_deg=orientation_cristal_deg)


def exporter_dxf(geom: GeometrieSpiral, chemin_dxf: str, cote_gravure: bool = False,
                  n_points: int = None):
    """Exporte les deux flancs (offset ±epaisseur/2) en LWPOLYLINE.

    cote_gravure=False (défaut) : la pièce FINIE, oxyde compris -- c'est la
    géométrie dont la physique a été calculée, et celle qu'il faut réimporter
    pour la boucle de validation (charger_dxf(exporter_dxf(...)) doit redonner
    la même physique à la précision du maillage près).

    cote_gravure=True : la cote de GRAVURE, avant oxydation -- c'est le masque
    dont l'atelier a besoin. Sur une lame oxydée les deux diffèrent de
    2 x 0.56 x épaisseur d'oxyde, soit typiquement 4 µm : confondre les deux
    fichiers donnerait une lame trop épaisse d'autant, donc une fréquence et
    un jeu entre spires faux. Le nom du calque porte la cote employée pour
    qu'un fichier isolé reste identifiable.

    n_points : nombre de points par flanc. Par défaut, celui du maillage de
    calcul. La densité qu'exige un fichier de gravure n'a aucune raison d'être
    celle qu'exige le calcul -- le maillage de calcul est un compromis entre
    précision et temps, et il densifie volontairement les zones terminales.
    Ré-échantillonner ne change pas la courbe, seulement sa description (voir
    GeometrieSpiral.reechantillonner)."""
    if n_points is not None and n_points != geom.n:
        geom = geom.reechantillonner(n_points)
    ext, inte = geom.contour_flancs(cote_gravure=cote_gravure)
    suffixe = "gravure" if cote_gravure else "fini"
    doc = ezdxf.new(setup=False)
    doc.header["$INSUNITS"] = 4  # millimètres
    msp = doc.modelspace()
    msp.add_lwpolyline(ext, format="xy", dxfattribs={"layer": f"spiral_flanc_ext_{suffixe}"})
    msp.add_lwpolyline(inte, format="xy", dxfattribs={"layer": f"spiral_flanc_int_{suffixe}"})
    doc.saveas(chemin_dxf)


def exporter_dxf_paire(geom: GeometrieSpiral, chemin_dxf: str, n_points: int = None):
    """Écrit les DEUX fichiers : la cote de gravure (suffixe `_gravure`) et la
    pièce finie (suffixe `_fini`). Retourne les deux chemins.

    Deux fichiers plutôt qu'un seul, parce que les deux servent à des choses
    différentes et que se tromper coûte cher : le masque part à la gravure, la
    pièce finie sert à vérifier et à mesurer. Quand il n'y a pas d'oxyde, les
    deux sont identiques et un seul fichier est écrit.

    n_points : voir exporter_dxf. Le ré-échantillonnage est fait UNE fois ici,
    pour que les deux fichiers décrivent rigoureusement la même fibre neutre --
    les ré-échantillonner séparément donnerait deux découpages identiques, mais
    reposerait sur un calcul refait deux fois pour rien."""
    base = chemin_dxf[:-4] if chemin_dxf.lower().endswith(".dxf") else chemin_dxf
    if n_points is not None and n_points != geom.n:
        geom = geom.reechantillonner(n_points)
    if geom.epaisseur_oxyde <= 0.0:
        exporter_dxf(geom, base + ".dxf", cote_gravure=False)
        return [base + ".dxf"]
    chemins = [base + "_gravure.dxf", base + "_fini.dxf"]
    exporter_dxf(geom, chemins[0], cote_gravure=True)
    exporter_dxf(geom, chemins[1], cote_gravure=False)
    return chemins
