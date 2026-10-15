"""Efforts quasi statiques, rendement et retard de l'échappement.

Quasi statique : la roue et l'ancre sont supposées sans inertie. À chaque
position du relevé cinématique, l'équilibre de la roue (couple moteur contre
effort de la palette) donne l'effort au contact dent-palette ; l'équilibre de
l'ancre donne l'effort au contact fourchette-ellipse, donc le couple sur le
balancier. Le frottement de Coulomb (coefficient mu) incline chaque effort de
l'angle de frottement, dans le sens opposé au glissement relatif.

On en tire :
  - le couple sur le balancier M(theta) au dégagement (il freine, avant la
    ligne des centres) et à l'impulsion (il pousse) ;
  - le rendement : énergie nette donnée au balancier / travail du couple
    moteur sur l'alternance (recul, chute et chemin perdu compris) ;
  - le retard dû à l'échappement, par la formule d'Airy. Pour theta =
    Theta0.cos(psi), un couple perturbateur M change la fréquence de

        df/f = -1 / (2.pi.k.Theta0) . integrale sur une période de M.cos(psi) dpsi
             = -1 / (2.pi.k.Theta0²) . somme des alternances de
               integrale de M(theta).theta / sqrt(Theta0² - theta²) dtheta

    (theta croissant dans chaque intégrale). Un couple moteur reçu avant la
    ligne des centres fait avancer, après elle fait retarder ; un couple de
    freinage fait l'inverse.

Ce modèle ne voit ni les chocs, ni l'inertie de la roue et de l'ancre : c'est
la référence statique du rendement, que la simulation dynamique corrigera.
"""

from dataclasses import dataclass

import numpy as np

from .cinematique import (DEG, angle_ligne_des_centres, appui_roue, bornes_balancier,
                          centre_fourchette, releve)
from .geometrie import distance


@dataclass
class Conditions:
    couple_roue: float = 1.3e-6        # N.m, couple moteur sur la roue d'échappement
    mu: float = 0.15                   # frottement silicium/silicium (dent-palette et fourchette-ellipse)
    frequence: float = 3.5             # Hz
    amplitude_deg: float = 310.0       # amplitude du balancier (demi-oscillation)
    inertie_balancier: float = 8e-10   # kg.m² -- PROVISOIRE, en attente des données CAO

    @property
    def raideur(self):
        """Couple de rappel du spiral, N.m/rad."""
        return self.inertie_balancier * (2 * np.pi * self.frequence) ** 2


def _contact(c1, a1, c2, a2):
    """Point de contact et normale unitaire dirigée de c2 vers c1."""
    _d, _i, _j, P1, P2 = distance(c1.segments_a(a1), c2.segments_a(a2))
    n = P1 - P2
    return 0.5 * (P1 + P2), n / np.hypot(*n)


def _vitesse(P, O, omega):
    r = P - O
    return omega * np.array([-r[1], r[0]])


def _cr(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def efforts_alternance(trace, alt, cond, pas=0.1 * DEG):
    """Couple sur le balancier le long d'une alternance.

    Rend un dict de tableaux : angle d'ancre, angle du balancier (depuis la
    ligne des centres), couple sur le balancier (N.m, sens trigonométrique),
    efforts normaux aux deux contacts, phase.
    """
    ev = alt.evenements
    s_a = np.sign(alt.arrivee - alt.depart)
    s_b = -s_a
    reprise = alt.ancre["reprise_ellipse"] * DEG + alt.a_mid
    zones = [("dégagement", alt.depart, ev["fin_degagement"]),
             ("impulsion", reprise, ev["chute"])]
    beta_ldc = angle_ligne_des_centres(trace)
    phi = alt.phi[0]
    lignes = []
    for phase, a0, a1 in zones:
        n = max(3, int(np.ceil(abs(a1 - a0) / pas)) + 1)
        # bornes légèrement rentrées : aux extrémités le contact change de nature
        for a in np.linspace(a0, a1, n)[1:-1]:
            phi = appui_roue(trace, a, phi)
            bas, haut = bornes_balancier(trace, a, centre_fourchette(trace, a))
            pousse = haut if s_b > 0 else bas
            mene = bas if s_b > 0 else haut
            beta = pousse if phase == "dégagement" else mene
            if beta is None:
                continue
            P, n_w = _contact(trace.roue, phi, trace.ancre, a)        # ancre -> roue
            Q, n_b = _contact(trace.balancier, beta, trace.ancre, a)  # fourchette -> ellipse
            lignes.append((phase, a, phi, beta, P, n_w, Q, n_b))
    if len(lignes) < 3:
        raise RuntimeError("alternance trop courte pour le calcul des efforts")

    A = np.array([l[1] for l in lignes])
    F = np.array([l[2] for l in lignes])
    B = np.array([l[3] for l in lignes])
    # rapports de vitesse par phase (dérivées le long de chaque zone)
    dphi, dbeta = np.zeros_like(A), np.zeros_like(A)
    phases = np.array([l[0] for l in lignes])
    for ph in ("dégagement", "impulsion"):
        k = phases == ph
        if k.sum() >= 2:
            dphi[k] = np.gradient(F[k], A[k])
            dbeta[k] = np.gradient(B[k], A[k])

    O_r, O_a, O_b = trace.roue.axe, trace.ancre.axe, trace.balancier.axe
    mu, M = cond.mu, cond.couple_roue
    couple_b, N_w, N_b = [], [], []
    for k, (ph, a, phi_k, beta, P, n_w, Q, n_b) in enumerate(lignes):
        w_a = s_a
        w_r, w_b = dphi[k] * w_a, dbeta[k] * w_a
        # contact roue-ancre : effort sur la roue N.(n - mu.s.t)
        t_w = np.array([-n_w[1], n_w[0]])
        s = np.sign((_vitesse(P, O_r, w_r) - _vitesse(P, O_a, w_a)) @ t_w)
        d_w = n_w - mu * s * t_w
        N = -trace.sens_roue * M / _cr(P - O_r, d_w)
        F_w = N * d_w
        # contact fourchette-ellipse : effort sur le balancier K.(m - mu.s'.t')
        t_b = np.array([-n_b[1], n_b[0]])
        s2 = np.sign((_vitesse(Q, O_b, w_b) - _vitesse(Q, O_a, w_a)) @ t_b)
        d_b = n_b - mu * s2 * t_b
        # équilibre de l'ancre : -P x F_w - Q x F_b = 0
        K = _cr(P - O_a, -F_w) / _cr(Q - O_a, d_b)
        couple_b.append(_cr(Q - O_b, K * d_b))
        N_w.append(N)
        N_b.append(K)
    return {"phase": phases, "alpha": A, "phi": F,
            "theta": B - beta_ldc, "couple": np.array(couple_b),
            "N_roue": np.array(N_w), "N_fourchette": np.array(N_b), "s_b": s_b}


def retard(efforts, cond, amplitudes_deg):
    """Marche due à l'échappement (s/j, négatif = retard) pour chaque amplitude."""
    out = []
    for Th in np.radians(np.atleast_1d(amplitudes_deg)):
        somme = 0.0
        for e in efforts:
            th, M = e["theta"], e["couple"]
            o = np.argsort(th)
            th, M = th[o], M[o]
            somme += np.trapezoid(M * th / np.sqrt(Th ** 2 - th ** 2), th)
        df_f = -somme / (2 * np.pi * cond.raideur * Th ** 2)
        out.append(86400.0 * df_f)
    return np.array(out)


def bilan(trace, cond=None, rel=None):
    """Rendement quasi statique et retard. Rend (resultats, efforts)."""
    cond = cond or Conditions()
    rel = rel or releve(trace)
    efforts, res = [], {"alternances": []}
    for alt in rel["alternances"]:
        e = efforts_alternance(trace, alt, cond)
        efforts.append(e)
        # énergie donnée au balancier : integrale de M dtheta le long du mouvement
        E = {}
        for ph in ("dégagement", "impulsion"):
            k = e["phase"] == ph
            E[ph] = float(np.trapezoid(e["couple"][k], e["theta"][k]))
        travail = trace.sens_roue * cond.couple_roue * (alt.phi[-1] - alt.phi[0])
        recul = cond.couple_roue * alt.recul_deg * DEG
        res["alternances"].append({
            "palette": alt.palette,
            "energie_degagement_J": E["dégagement"], "energie_impulsion_J": E["impulsion"],
            "travail_moteur_J": travail, "rendement": (E["dégagement"] + E["impulsion"]) / travail,
            "avance_roue_deg": abs(alt.phi[-1] - alt.phi[0]) / DEG,
            "couple_max_balancier": float(np.max(np.abs(e["couple"]))),
            "travail_recul_J": recul,
            "effort_negatif": bool(np.any(e["N_roue"] <= 0) or np.any(e["N_fourchette"] <= 0)),
        })
    E_tot = sum(a["energie_degagement_J"] + a["energie_impulsion_J"] for a in res["alternances"])
    W_tot = sum(a["travail_moteur_J"] for a in res["alternances"])
    res["rendement"] = E_tot / W_tot
    res["energie_par_oscillation_J"] = E_tot
    res["energie_balancier_J"] = 0.5 * cond.raideur * np.radians(cond.amplitude_deg) ** 2
    res["retard_s_j"] = float(retard(efforts, cond, cond.amplitude_deg)[0])
    return res, efforts


def imprimer(res, cond):
    print(f"Conditions : couple roue {cond.couple_roue * 1e6:.2f} µN.m, mu = {cond.mu}, "
          f"f = {cond.frequence} Hz, amplitude {cond.amplitude_deg:.0f}°, "
          f"inertie balancier {cond.inertie_balancier:.3g} kg.m² (provisoire)")
    for a in res["alternances"]:
        print(f"  alternance {a['palette']:7s}: roue {a['avance_roue_deg']:.3f}°, travail moteur "
              f"{a['travail_moteur_J'] * 1e9:.2f} nJ ; balancier : dégagement "
              f"{a['energie_degagement_J'] * 1e9:+.2f} nJ, impulsion {a['energie_impulsion_J'] * 1e9:+.2f} nJ"
              f" -> rendement {100 * a['rendement']:.1f} % ; couple max {a['couple_max_balancier'] * 1e9:.0f} nN.m"
              + ("  [effort négatif : contact perdu]" if a["effort_negatif"] else ""))
    print(f"Rendement quasi statique sur une oscillation : {100 * res['rendement']:.1f} %")
    print(f"Énergie reçue par oscillation {res['energie_par_oscillation_J'] * 1e9:.2f} nJ, "
          f"énergie du balancier à {cond.amplitude_deg:.0f}° : {res['energie_balancier_J'] * 1e9:.0f} nJ")
    print(f"Retard dû à l'échappement à {cond.amplitude_deg:.0f}° : {res['retard_s_j']:+.2f} s/j")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    from .dxf_import import lire_trace
    chemin = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("simulation.dxf"))
    t = lire_trace(chemin)
    c = Conditions()
    rel = releve(t)
    res, eff = bilan(t, c, rel)
    imprimer(res, c)
    amps = np.array([150, 180, 200, 230, 260, 290, 310, 330])
    print("Retard (s/j) selon l'amplitude : "
          + ", ".join(f"{a}° {r:+.2f}" for a, r in zip(amps, retard(eff, c, amps))))
