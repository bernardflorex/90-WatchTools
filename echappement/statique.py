"""Efforts quasi statiques, rendement et retard de l'échappement.

Quasi statique : la roue et l'ancre sont supposées sans inertie. À chaque
position du relevé cinématique, l'équilibre de la roue (couple moteur contre
effort de la palette) donne l'effort au contact dent-palette ; l'équilibre de
l'ancre donne l'effort au contact fourchette-ellipse, donc le couple sur le
balancier. Le frottement de Coulomb (mu) incline chaque effort de l'angle de
frottement, contre le glissement relatif.

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
la référence statique, que la simulation dynamique corrige.
"""

import numpy as np

from .cinematique import Releve, releve
from .parametres import Parametres
from .tables import DEG


def _vitesse(P, O, omega):
    r = P - O
    return omega * np.array([-r[1], r[0]])


def _cr(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def efforts_contact(trace, P, n_w, Q, n_b, w_r, w_a, w_b, mu, couple, sens):
    """Efforts quasi statiques pour une position donnée.

    P, n_w : contact dent-palette, normale palette -> roue ;
    Q, n_b : contact fourchette-ellipse, normale fourchette -> ellipse ;
    w_* : vitesses angulaires (à un facteur près), pour le sens du glissement.
    Rend (couple sur le balancier, effort normal dent-palette, effort normal
    fourchette-ellipse)."""
    O_r, O_a, O_b = trace.roue.axe, trace.ancre.axe, trace.balancier.axe
    t_w = np.array([-n_w[1], n_w[0]])
    s = np.sign((_vitesse(P, O_r, w_r) - _vitesse(P, O_a, w_a)) @ t_w)
    d_w = n_w - mu * s * t_w
    N = -sens * couple / _cr(P - O_r, d_w)
    F_w = N * d_w
    t_b = np.array([-n_b[1], n_b[0]])
    s2 = np.sign((_vitesse(Q, O_b, w_b) - _vitesse(Q, O_a, w_a)) @ t_b)
    d_b = n_b - mu * s2 * t_b
    K = _cr(P - O_a, -F_w) / _cr(Q - O_a, d_b)
    return _cr(Q - O_b, K * d_b), N, K


def efforts_alternance(rel: Releve, alt, par: Parametres):
    """Couple sur le balancier le long d'une alternance, aux points de grille.

    Rend un dict de tableaux : alpha, theta (balancier depuis la ligne des
    centres, rad), couple (N.m, sens trigonométrique), efforts normaux, phase.
    """
    T = rel.T
    trace = T.trace
    ev = alt.evenements
    s_a = np.sign(alt.arrivee - alt.depart)
    tab = T.roue[alt.palette]
    dphi = np.gradient(tab.phi, T.alpha)
    lignes = []
    for phase, a0, a1, borne in (("dégagement", alt.depart, ev["fin_degagement"], alt.pousse),
                                 ("impulsion", ev["reprise_ellipse"], ev["chute"], alt.mene)):
        dbeta = np.gradient(borne.beta, T.alpha)
        lo, hi = sorted((a0, a1))
        for i in np.nonzero((T.alpha > lo) & (T.alpha < hi))[0]:
            if not (np.isfinite(tab.phi[i]) and np.isfinite(borne.beta[i])
                    and np.isfinite(dphi[i]) and np.isfinite(dbeta[i])):
                continue
            M, N, K = efforts_contact(trace, tab.P[i], tab.N[i], borne.P[i], borne.N[i],
                                      dphi[i] * s_a, s_a, dbeta[i] * s_a,
                                      par.mu, par.couple_roue, trace.sens_roue)
            lignes.append((phase, T.alpha[i], borne.beta[i] - T.beta_ldc, M, N, K))
    lignes.sort(key=lambda l: l[1] * s_a)
    col = lambda k: np.array([l[k] for l in lignes])
    return {"phase": np.array([l[0] for l in lignes]), "alpha": col(1), "theta": col(2),
            "couple": col(3), "N_roue": col(4), "N_fourchette": col(5), "s_b": alt.s_b}


def retard(efforts, par: Parametres, amplitudes_deg):
    """Marche due à l'échappement (s/j, négatif = retard) pour chaque amplitude."""
    out = []
    for Th in np.radians(np.atleast_1d(amplitudes_deg)):
        somme = 0.0
        for e in efforts:
            th, M = e["theta"], e["couple"]
            o = np.argsort(th)
            th, M = th[o], M[o]
            somme += np.trapezoid(M * th / np.sqrt(Th ** 2 - th ** 2), th)
        out.append(86400.0 * (-somme / (2 * np.pi * par.raideur * Th ** 2)))
    return np.array(out)


def bilan(rel: Releve, par: Parametres = None):
    """Rendement quasi statique et retard. Rend (résultats, efforts)."""
    par = par or Parametres()
    efforts, res = [], {"alternances": {}}
    for alt in rel.alternances:
        e = efforts_alternance(rel, alt, par)
        efforts.append(e)
        E = {}
        for ph in ("dégagement", "impulsion"):
            k = e["phase"] == ph
            E[ph] = float(np.trapezoid(e["couple"][k], e["theta"][k])) if k.sum() > 1 else 0.0
        travail = par.couple_roue * alt.mesures["avance_roue_deg"] * DEG
        res["alternances"][alt.palette] = {
            "energie_degagement": E["dégagement"], "energie_impulsion": E["impulsion"],
            "travail_moteur": travail, "rendement": (E["dégagement"] + E["impulsion"]) / travail,
            "couple_max": float(np.max(np.abs(e["couple"]))),
            "contact_perdu": bool(np.any(e["N_roue"] <= 0) or np.any(e["N_fourchette"] <= 0)),
        }
    alts = res["alternances"].values()
    E_tot = sum(a["energie_degagement"] + a["energie_impulsion"] for a in alts)
    W_tot = sum(a["travail_moteur"] for a in alts)
    res["rendement"] = E_tot / W_tot
    res["energie_par_oscillation"] = E_tot
    res["travail_par_oscillation"] = W_tot
    res["retard_s_j"] = float(retard(efforts, par, par.amplitude_deg)[0])
    # amplitude d'équilibre : énergie reçue = pertes du balancier (Q), qui
    # croissent comme le carré de l'amplitude
    res["amplitude_equilibre_deg"] = float(np.degrees(np.sqrt(
        E_tot * par.facteur_qualite / (np.pi * par.raideur))))
    return res, efforts


def tableau(res, par: Parametres):
    a = res["alternances"]
    l = [("Rendement statique (%)", f"{100 * a['entrée']['rendement']:.1f}",
          f"{100 * a['sortie']['rendement']:.1f}"),
         ("Travail moteur (nJ)", f"{1e9 * a['entrée']['travail_moteur']:.1f}",
          f"{1e9 * a['sortie']['travail_moteur']:.1f}"),
         ("Énergie au dégagement (nJ)", f"{1e9 * a['entrée']['energie_degagement']:+.1f}",
          f"{1e9 * a['sortie']['energie_degagement']:+.1f}"),
         ("Énergie à l'impulsion (nJ)", f"{1e9 * a['entrée']['energie_impulsion']:+.1f}",
          f"{1e9 * a['sortie']['energie_impulsion']:+.1f}"),
         ("Couple maximal au balancier (nN.m)", f"{1e9 * a['entrée']['couple_max']:.0f}",
          f"{1e9 * a['sortie']['couple_max']:.0f}"),
         ("Rendement statique sur une oscillation (%)", f"{100 * res['rendement']:.1f}", ""),
         (f"Retard à {par.amplitude_deg:.0f}° (s/j)", f"{res['retard_s_j']:+.2f}", ""),
         ("Amplitude d'équilibre, statique (°)", f"{res['amplitude_equilibre_deg']:.0f}", "")]
    return l


if __name__ == "__main__":
    import sys
    from pathlib import Path
    from .dxf_import import lire_trace
    chemin = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("simulation.dxf"))
    rel = releve(lire_trace(chemin))
    p = Parametres()
    res, eff = bilan(rel, p)
    for nom, e, s in tableau(res, p):
        print(f"{nom:46s} {e:>10s} {s:>10s}")
    amps = np.array([150, 200, 250, 310, 330])
    print("Retard (s/j) : " + ", ".join(f"{a}° {r:+.2f}" for a, r in zip(amps, retard(eff, p, amps))))
