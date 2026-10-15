"""Relevé cinématique du tracé : ce que fait l'échappement, sans les efforts.

L'ancre est menée pas à pas d'une butée à l'autre. À chaque position, la roue
est poussée dans son sens de rotation jusqu'au contact (position d'appui) et
le balancier est placé contre l'un ou l'autre flanc d'entrée de la fourchette.
On en déduit, pour chaque palette :

  - la chute : rotation libre de la roue entre la fin de l'impulsion et le
    premier contact sur la face de repos de l'autre palette ;
  - le chemin perdu : course de la dent sur la face de repos entre ce premier
    contact et la position de repos, ancre sur la butée (en µm) ;
  - le repos : distance, sur la face de repos, de la dent à l'arête de repos,
    ancre sur la butée ;
  - le dégagement, l'impulsion sur le plan de la palette puis sur le plan de
    la dent, en angle d'ancre ;
  - le recul de la roue au dégagement (effet du tirage) ;
  - côté balancier : l'angle d'entrée de l'ellipse dans la fourchette, la fin
    du dégagement, le début et la fin de l'impulsion, comptés depuis la LIGNE
    DES CENTRES ancre-balancier. Leur répartition de part et d'autre de cette
    ligne est ce qui fait l'erreur de marche de l'échappement.

Convention : angles dans le sens trigonométrique de la vue standard ; les
angles d'ancre sont comptés depuis la position médiane entre les butées.
"""

from dataclasses import dataclass, field

import numpy as np

from .geometrie import point_de_contact, rotation, tourner_jusquau_contact

DEG = np.pi / 180.0


def appui_roue(trace, alpha, phi_depart):
    """Position de la roue poussée contre l'ancre posée à `alpha`."""
    return tourner_jusquau_contact(trace.roue, phi_depart, trace.sens_roue,
                                   [(trace.ancre, alpha)], pas=0.02 * DEG,
                                   course=30 * DEG)


def contact_roue(trace, alpha, phi):
    """(point, étiquette sur l'ancre, étiquette sur la roue)."""
    P, e_roue, e_ancre = point_de_contact(trace.roue, phi, trace.ancre, alpha)
    return P, (e_ancre.etiquette if e_ancre else ""), (e_roue.etiquette if e_roue else "")


def course_butees(trace):
    """Angles d'ancre (depuis la pose dessinée) des deux butées."""
    a = [tourner_jusquau_contact(trace.ancre, 0.0, s, [(b, 0.0) for b in trace.butees],
                                 pas=0.05 * DEG, course=30 * DEG) for s in (-1, 1)]
    return min(a), max(a)


def angle_ligne_des_centres(trace):
    """Angle du balancier (depuis la pose dessinée) qui met le centre de
    l'ellipse sur la ligne balancier -> ancre."""
    c = next(e.centre for e in trace.balancier.elements if e.est_arc and e.rayon > 0.1)
    v = c - trace.balancier.axe
    w = trace.ancre.axe - trace.balancier.axe
    return float(np.arctan2(v[0] * w[1] - v[1] * w[0], v @ w))


def bornes_balancier(trace, alpha, beta_depart):
    """Les deux positions du balancier où l'ellipse touche la fourchette
    posée à `alpha` : (beta_bas, beta_haut). L'ellipse y est libre entre les
    deux (ébat)."""
    bas = tourner_jusquau_contact(trace.balancier, beta_depart, -1, [(trace.ancre, alpha)],
                                  pas=0.05 * DEG, course=40 * DEG)
    haut = tourner_jusquau_contact(trace.balancier, beta_depart, +1, [(trace.ancre, alpha)],
                                   pas=0.05 * DEG, course=40 * DEG)
    return bas, haut


def centre_fourchette(trace, alpha):
    """Angle du balancier qui place le centre de l'ellipse sur l'axe de la
    fourchette posée à `alpha` (milieu de l'entrée de fourchette)."""
    bag = trace.ancre.element("baguette").points
    u = rotation(bag[-1], trace.ancre.axe, alpha) - trace.ancre.axe
    u /= np.hypot(*u)
    c = next(e.centre for e in trace.balancier.elements if e.est_arc and e.rayon > 0.1)
    rho = np.hypot(*(c - trace.balancier.axe))
    # intersection de la droite axe_ancre + t.u avec le cercle de l'ellipse
    o = trace.ancre.axe - trace.balancier.axe
    b, cc = o @ u, o @ o - rho ** 2
    t = -b - np.sqrt(max(b * b - cc, 0.0))
    p = trace.ancre.axe + t * u - trace.balancier.axe
    v = c - trace.balancier.axe
    return float(np.arctan2(v[0] * p[1] - v[1] * p[0], v @ p))


@dataclass
class Alternance:
    """Relevé d'une alternance : l'ancre va de `butee_depart` à `butee_arrivee`."""
    palette: str                  # palette qui travaille (repos, dégagement, impulsion)
    alpha: np.ndarray = None      # angles d'ancre parcourus (rad, pose dessinée)
    phi: np.ndarray = None        # position d'appui de la roue
    etiquettes: list = field(default_factory=list)
    evenements: dict = field(default_factory=dict)


def _classe(etiquette_ancre, etiquette_roue):
    if etiquette_ancre.startswith("face de repos"):
        return "repos"
    if etiquette_ancre.startswith("arête de sortie") and etiquette_roue == "plan de la dent":
        return "impulsion dent"
    return "impulsion palette"


def balayer(trace, a_depart, a_arrivee, phi_depart, pas=0.05 * DEG):
    """Mène l'ancre de a_depart à a_arrivee, la roue en appui."""
    n = max(2, int(np.ceil(abs(a_arrivee - a_depart) / pas)) + 1)
    A = np.linspace(a_depart, a_arrivee, n)
    phi, P, lab = [], [], []
    f = phi_depart
    for a in A:
        f = appui_roue(trace, a, f)
        pt, ea, er = contact_roue(trace, a, f)
        phi.append(f)
        P.append(trace.ancre.vers_local(pt, a))
        lab.append((ea, er))
    return A, np.array(phi), np.array(P), lab


def _affiner(trace, a0, a1, f0, critere, n=30):
    """Angle d'ancre où `critere(étiquettes)` bascule, entre a0 (faux) et a1 (vrai)."""
    for _ in range(n):
        m = 0.5 * (a0 + a1)
        fm = appui_roue(trace, m, f0)
        if critere(contact_roue(trace, m, fm)):
            a1 = m
        else:
            a0, f0 = m, fm
    return a0, a1, f0


def releve(trace, pas=0.05 * DEG):
    """Relevé complet du tracé, à partir de la pose dessinée."""
    r = {}
    a_min, a_max = course_butees(trace)
    a_mid = 0.5 * (a_min + a_max)
    r["course_ancre_deg"] = (a_max - a_min) / DEG
    r["ancre_dessinee_deg"] = -a_mid / DEG
    beta_ldc = angle_ligne_des_centres(trace)
    r["balancier_dessine_deg"] = -beta_ldc / DEG

    # Pose dessinée : la roue tombe, puis l'ancre gagne la butée la plus proche.
    phi0 = appui_roue(trace, 0.0, 0.0)
    P0, ea, er = contact_roue(trace, 0.0, phi0)
    r["chute_dessinee_deg"] = abs(phi0) / DEG
    r["chute_dessinee_contact"] = ea
    butee_proche = a_min if abs(a_min) < abs(a_max) else a_max
    butee_loin = a_max if butee_proche == a_min else a_min

    alternances = []
    depart, arrivee, f = butee_proche, butee_loin, appui_roue(trace, butee_proche, phi0)
    for _ in range(2):
        A, F, P, lab = balayer(trace, depart, arrivee, f, pas)
        cote = lab[0][0].split()[-1]                     # palette au repos au départ
        alt = Alternance(palette=cote, alpha=A, phi=F, etiquettes=lab)
        classes = [_classe(*l) if l[0].endswith(cote) else "chute" for l in lab]
        ev = alt.evenements
        # fin du dégagement : la dent quitte la face de repos
        k = next(i for i, c in enumerate(classes) if c != "repos")
        ev["fin_degagement"], _, _ = _affiner(
            trace, A[k - 1], A[k], F[k - 1], lambda c: not c[1].startswith("face de repos"))
        # passage du plan de la palette au plan de la dent
        k2 = next((i for i, c in enumerate(classes) if c == "impulsion dent"), None)
        if k2 is not None:
            ev["debut_impulsion_dent"], _, _ = _affiner(
                trace, A[k2 - 1], A[k2], F[k2 - 1],
                lambda c: c[1].startswith("arête de sortie") and c[2] == "plan de la dent")
        # chute : la dent quitte la palette, la suivante tombe sur l'autre
        k3 = next(i for i, c in enumerate(classes) if c == "chute")
        a_chute, a_apres, f_av = _affiner(
            trace, A[k3 - 1], A[k3], F[k3 - 1], lambda c: not c[1].endswith(cote))
        # juste après la fin de l'impulsion, la roue avance librement jusqu'à
        # l'autre palette
        f_ap = appui_roue(trace, a_apres, f_av)
        ev["chute"] = a_apres
        alt.chute_deg = abs(f_ap - f_av) / DEG
        # chemin perdu et repos sur l'autre palette
        autre = contact_roue(trace, a_apres, f_ap)[1].split()[-1]
        p_chute = trace.ancre.vers_local(contact_roue(trace, a_apres, f_ap)[0], a_apres)
        f_b = appui_roue(trace, arrivee, f_ap)
        p_butee = trace.ancre.vers_local(contact_roue(trace, arrivee, f_b)[0], arrivee)
        face = trace.ancre.element(f"face de repos {autre}").points
        arete = trace.ancre.element(f"arête de repos {autre}").points.mean(0)
        coin = face[0] if np.hypot(*(face[0] - arete)) < np.hypot(*(face[-1] - arete)) else face[-1]
        alt.palette_suivante = autre
        alt.chemin_perdu_um = 1000 * np.hypot(*(p_butee - p_chute))
        alt.repos_chute_um = 1000 * np.hypot(*(p_chute - coin))
        alt.repos_butee_um = 1000 * np.hypot(*(p_butee - coin))
        # recul de la roue pendant le dégagement
        i_deg = [i for i, c in enumerate(classes) if c == "repos"]
        alt.recul_deg = max(0.0, float(np.max(-trace.sens_roue * (F[i_deg] - F[0])))) / DEG
        alt.depart, alt.arrivee = depart, arrivee
        alternances.append(alt)
        depart, arrivee, f = arrivee, depart, f_b

    # Balancier : bornes de l'ellipse dans la fourchette sur toute la course.
    A = np.linspace(a_min, a_max, 145)
    # None : aucun flanc de ce côté, l'ellipse sort librement par l'entrée de
    # fourchette (c'est le cas aux butées, où elle entre et sort)
    bornes = np.array([[np.nan if v is None else v
                        for v in bornes_balancier(trace, a, centre_fourchette(trace, a))]
                       for a in A], dtype=float)
    for alt in alternances:
        s_a = np.sign(alt.arrivee - alt.depart)
        s_b = -s_a                                    # le balancier tourne en sens inverse
        pousse = bornes[:, 1] if s_b > 0 else bornes[:, 0]     # l'ellipse mène
        mene = bornes[:, 0] if s_b > 0 else bornes[:, 1]       # la fourchette mène
        interp = lambda B, a: float(np.interp(a, A, B))
        ev = alt.evenements
        b_entree = interp(pousse, alt.depart)
        b_deg = interp(pousse, ev["fin_degagement"])
        # l'ancre rattrape l'ébat : la fourchette reprend l'ellipse en b_deg
        ok = np.isfinite(mene)
        ordre = np.argsort((mene * s_b)[ok])
        a_reprise = float(np.interp(b_deg * s_b, (mene * s_b)[ok][ordre], A[ok][ordre]))
        b_fin = interp(mene, ev["chute"])
        rel = lambda b: s_b * (b - beta_ldc) / DEG    # >0 : après la ligne des centres
        alt.balancier = {
            "entree_ellipse": rel(b_entree), "fin_degagement": rel(b_deg),
            "reprise_impulsion": rel(b_deg), "fin_impulsion": rel(b_fin),
            "levee_totale": abs(b_fin - b_entree) / DEG,
        }
        alt.a_mid = a_mid
        alt.ancre = {
            "depart": (alt.depart - a_mid) / DEG, "fin_degagement": (ev["fin_degagement"] - a_mid) / DEG,
            "reprise_ellipse": (a_reprise - a_mid) / DEG,
            "debut_impulsion_dent": ((ev["debut_impulsion_dent"] - a_mid) / DEG
                                     if "debut_impulsion_dent" in ev else None),
            "chute": (ev["chute"] - a_mid) / DEG, "arrivee": (alt.arrivee - a_mid) / DEG,
        }
    r["alternances"] = alternances
    return r


def imprimer(r, trace):
    print(f"Roue : {trace.nombre_dents} dents, rayon de tête {trace.rayon_tete:.4f} mm, "
          f"pas à la tête {1000 * trace.pas * trace.rayon_tete:.0f} µm")
    print(f"Ancre : course entre butées {r['course_ancre_deg']:.3f}°, dessinée à "
          f"{r['ancre_dessinee_deg']:+.3f}° de la médiane ; balancier dessiné à "
          f"{r['balancier_dessine_deg']:+.2f}° de la ligne des centres")
    print(f"Pose dessinée : chute de {r['chute_dessinee_deg']:.3f}° de roue avant contact "
          f"sur « {r['chute_dessinee_contact']} »")
    for alt in r["alternances"]:
        a, b = alt.ancre, alt.balancier
        print(f"\nAlternance « {alt.palette} » (ancre de {a['depart']:+.2f}° à {a['arrivee']:+.2f}°)")
        print(f"  ancre  : fin du dégagement {a['fin_degagement']:+.3f}°, ellipse reprise "
              f"{a['reprise_ellipse']:+.3f}°, impulsion sur la dent dès "
              + (f"{a['debut_impulsion_dent']:+.3f}°" if a['debut_impulsion_dent'] is not None else "—")
              + f", chute à {a['chute']:+.3f}°")
        print(f"  roue   : recul au dégagement {alt.recul_deg:.3f}°, chute {alt.chute_deg:.3f}° "
              f"({1000 * alt.chute_deg * DEG * trace.rayon_tete:.1f} µm à la tête)")
        print(f"  palette {alt.palette_suivante} reçoit la dent : repos à la chute "
              f"{alt.repos_chute_um:.1f} µm, chemin perdu {alt.chemin_perdu_um:.1f} µm, "
              f"repos sur butée {alt.repos_butee_um:.1f} µm")
        print(f"  balancier (depuis la ligne des centres, <0 avant) : entrée de l'ellipse "
              f"{b['entree_ellipse']:+.2f}°, fin du dégagement {b['fin_degagement']:+.2f}°, "
              f"fin d'impulsion {b['fin_impulsion']:+.2f}° ; levée {b['levee_totale']:.2f}°")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    from .dxf_import import lire_trace
    chemin = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("simulation.dxf"))
    t = lire_trace(chemin)
    imprimer(releve(t), t)
