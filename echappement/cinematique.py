"""Relevé cinématique du tracé : ce que fait l'échappement, sans les efforts.

L'ancre est menée d'une butée à l'autre ; à chaque position la roue est
poussée dans son sens jusqu'au contact (position d'appui), et le balancier est
placé contre l'un ou l'autre flanc de la fourchette. Tout est lu dans les
tables de contact (tables.py). On en déduit, pour chaque alternance :

  - le dégagement, l'impulsion sur le plan de la palette puis sur le plan de
    la dent, en angle d'ancre, et le recul de la roue au dégagement ;
  - la chute : rotation libre de la roue entre l'échappement de la dent et
    le premier contact sur la face de repos de l'autre palette ;
  - le chemin perdu : course de la dent sur la face de repos entre ce premier
    contact et la position de repos, ancre sur la butée (µm) ;
  - le repos à la chute et sur la butée (µm, depuis l'arête de repos) ;
  - côté balancier : entrée de l'ellipse dans la fourchette, fin du
    dégagement, fin de l'impulsion, comptés depuis la LIGNE DES CENTRES
    ancre-balancier, positifs après elle dans le sens du mouvement.

Les angles d'ancre affichés sont comptés depuis la médiane des butées.
"""

from dataclasses import dataclass, field

import numpy as np

from .geometrie import rotation
from .tables import DEG, RECUL, Tables

#: au-delà de ce saut de la roue entre deux pas, il y a eu chute
SAUT_CHUTE = 0.2 * DEG


def _wrap(x, p, marge):
    """Ramène x dans [-marge, p - marge)."""
    return np.mod(x + marge, p) - marge


@dataclass
class Alternance:
    palette: str                      # palette qui travaille
    palette_suivante: str = ""        # palette qui reçoit la dent à la chute
    depart: float = 0.0               # angle d'ancre (rad, pose dessinée)
    arrivee: float = 0.0
    s_b: int = 0                      # sens de rotation du balancier
    alpha: np.ndarray = None          # grille parcourue
    phi: np.ndarray = None            # position d'appui de la roue
    evenements: dict = field(default_factory=dict)   # angles d'ancre (rad, pose dessinée)
    mesures: dict = field(default_factory=dict)


class Releve:
    """Relevé cinématique complet d'un tracé."""

    def __init__(self, T: Tables):
        self.T = T
        t = T.trace
        self.s = t.sens_roue
        self.p = t.pas
        self.a_mid = 0.5 * (T.a_min + T.a_max)
        # pose dessinée : la roue finit sa chute, l'ancre gagne la butée proche
        # pose dessinée, ramenée entre les butées si celles-ci ont été resserrées
        i0 = int(np.clip(np.argmin(np.abs(T.alpha)), T.i_min, T.i_max))
        phi0 = self.appui(i0, 0.0)[0]
        self.chute_dessinee = abs(phi0)
        proche = T.i_min if abs(T.a_min) < abs(T.a_max) else T.i_max
        phi = phi0
        for i in (range(i0, proche - 1, -1) if proche < i0 else range(i0, proche + 1)):
            phi = self.appui(i, phi)[0]
        self.alternances = []
        for _ in range(2):
            alt = self._alternance(proche, phi)
            self.alternances.append(alt)
            proche = T.i_max if proche == T.i_min else T.i_min
            phi = alt.phi[-1]

    # -- roue en appui -------------------------------------------------------
    def appui(self, i, phi):
        """Position d'appui de la roue à l'angle d'ancre d'indice i, en partant
        de phi : (phi, palette touchée)."""
        best, cote = None, None
        for c, tab in self.T.roue.items():
            if not np.isfinite(tab.phi[i]):
                continue
            d = _wrap(self.s * (tab.phi[i] - phi), self.p, RECUL)
            if best is None or d < best:
                best, cote = d, c
        if best is None:
            return phi, None
        return phi + self.s * best, cote

    def _etiquette(self, cote, i):
        tab = self.T.roue[cote]
        return tab.sur_ancre[i], tab.sur_roue[i]

    # -- une alternance ------------------------------------------------------
    def _alternance(self, i_dep, phi):
        T = self.T
        i_arr = T.i_max if i_dep == T.i_min else T.i_min
        pas = 1 if i_arr > i_dep else -1
        idx = list(range(i_dep, i_arr + pas, pas))
        Phi, cotes = [], []
        for i in idx:
            phi, c = self.appui(i, phi)
            Phi.append(phi)
            cotes.append(c)
        Phi = np.array(Phi)
        travail = cotes[0]
        alt = Alternance(palette=travail, depart=T.alpha[i_dep], arrivee=T.alpha[i_arr],
                         alpha=T.alpha[idx], phi=Phi)
        s_a = np.sign(alt.arrivee - alt.depart)
        alt.s_b = int(-s_a)
        alt.idx = idx
        ev, me = alt.evenements, alt.mesures

        def premier(cond):
            return next(k for k in range(1, len(idx)) if cond(k))

        # fin du dégagement : la dent quitte la face de repos
        k = premier(lambda k: cotes[k] == travail
                    and not self._etiquette(travail, idx[k])[0].startswith("face de repos"))
        ev["fin_degagement"] = self._affiner(travail, idx[k - 1], idx[k],
                                             lambda ea, er: not ea.startswith("face de repos"))
        # impulsion sur le plan de la dent
        k2 = next((k for k in range(1, len(idx)) if cotes[k] == travail
                   and self._etiquette(travail, idx[k])[1] == "plan de la dent"
                   and self._etiquette(travail, idx[k])[0].startswith("arête de sortie")), None)
        if k2 is not None:
            ev["debut_impulsion_dent"] = self._affiner(
                travail, idx[k2 - 1], idx[k2],
                lambda ea, er: er == "plan de la dent" and ea.startswith("arête de sortie"))
        # chute : la palette de travail lâche la dent
        k3 = premier(lambda k: cotes[k] != travail or abs(Phi[k] - Phi[k - 1]) > SAUT_CHUTE)
        a_e, phi_e = self._echappement(travail, idx[k3 - 1], idx[k3], Phi[k3 - 1])
        ev["chute"] = a_e
        autre = "sortie" if travail == "entrée" else "entrée"
        alt.palette_suivante = autre
        atterr = T.contact_roue(autre, a_e, phi_e)
        me["chute_deg"] = abs(atterr[0] - phi_e) / DEG
        alt.phi_echappement, alt.phi_atterrissage = phi_e, atterr[0]
        p_chute = rotation(atterr[1], T.trace.ancre.axe, -a_e)
        # repos sur butée : contact de l'autre palette à l'arrivée
        tab = T.roue[autre]
        p_butee = rotation(tab.P[i_arr], T.trace.ancre.axe, -T.alpha[i_arr])
        face = T.trace.ancre.element(f"face de repos {autre}").points
        arete = T.trace.ancre.element(f"arête de repos {autre}").points.mean(0)
        coin = min((face[0], face[-1]), key=lambda q: np.hypot(*(q - arete)))
        me["chemin_perdu_um"] = 1000 * np.hypot(*(p_butee - p_chute))
        me["repos_chute_um"] = 1000 * np.hypot(*(p_chute - coin))
        me["repos_butee_um"] = 1000 * np.hypot(*(p_butee - coin))
        deg = [k for k in range(len(idx)) if T.alpha[idx[k]] * s_a <= ev["fin_degagement"] * s_a]
        me["recul_deg"] = max(0.0, float(np.max(-self.s * (Phi[deg] - Phi[0])))) / DEG
        me["avance_roue_deg"] = abs(Phi[-1] - Phi[0]) / DEG

        # balancier
        bas, haut = T.fourche_bas.beta, T.fourche_haut.beta
        alt.pousse, alt.mene = (T.fourche_haut, T.fourche_bas) if alt.s_b > 0 \
            else (T.fourche_bas, T.fourche_haut)
        pousse, mene = alt.pousse.beta, alt.mene.beta
        b_entree = T.interp(pousse, alt.depart)
        b_deg = T.interp(pousse, ev["fin_degagement"])
        ok = np.isfinite(mene)
        o = np.argsort((mene * alt.s_b)[ok])
        ev["reprise_ellipse"] = float(np.interp(b_deg * alt.s_b, (mene * alt.s_b)[ok][o],
                                                T.alpha[ok][o]))
        b_fin = T.interp(mene, a_e)
        rel = lambda b: alt.s_b * (b - T.beta_ldc) / DEG
        me["entree_ellipse_deg"] = rel(b_entree)
        me["fin_degagement_balancier_deg"] = rel(b_deg)
        me["fin_impulsion_balancier_deg"] = rel(b_fin)
        me["levee_deg"] = abs(b_fin - b_entree) / DEG
        return alt

    def _affiner(self, cote, i0, i1, critere, n=30):
        """Angle d'ancre où le contact sur `cote` bascule (critere sur les
        étiquettes devient vrai), entre les indices i0 et i1."""
        T = self.T
        a0, a1 = T.alpha[i0], T.alpha[i1]
        phi0 = T.roue[cote].phi[i0]
        for _ in range(n):
            m = 0.5 * (a0 + a1)
            r = T.contact_roue(cote, m, phi0 - self.s * RECUL)
            if r is not None and critere(r[3], r[4]):
                a1 = m
            else:
                a0 = m
                if r is not None:
                    phi0 = r[0]
        return 0.5 * (a0 + a1)

    def _echappement(self, cote, i0, i1, phi0, n=30):
        """Dernier angle d'ancre où la palette `cote` retient la dent, et la
        position de la roue à cet instant."""
        T = self.T
        a0, a1 = T.alpha[i0], T.alpha[i1]
        for _ in range(n):
            m = 0.5 * (a0 + a1)
            r = T.contact_roue(cote, m, phi0 - self.s * RECUL)
            if r is not None and abs(r[0] - phi0) < SAUT_CHUTE:
                a0, phi0 = m, r[0]
            else:
                a1 = m
        return a0, phi0

    # -- présentation --------------------------------------------------------
    def depuis_mediane(self, a):
        return (a - self.a_mid) / DEG

    def tableau(self):
        """Lignes (grandeur, alternance d'entrée, alternance de sortie).
        Chaque alternance est désignée par la palette qui y travaille ; la
        chute, le repos et le chemin perdu sont mesurés sur l'autre palette,
        celle qui reçoit la dent."""
        par = {a.palette: a for a in self.alternances}
        T = self.T
        lignes = [("Course de l'ancre entre butées (°)", f"{(T.a_max - T.a_min) / DEG:.3f}", ""),
                  ("Chute du dessin (°)", f"{self.chute_dessinee / DEG:.3f}", "")]

        def l(nom, f):
            lignes.append((nom, f(par["entrée"]), f(par["sortie"])))
        d = self.depuis_mediane
        l("Fin du dégagement, ancre (°)", lambda a: f"{d(a.evenements['fin_degagement']):+.3f}")
        l("Impulsion sur la dent dès, ancre (°)",
          lambda a: f"{d(a.evenements['debut_impulsion_dent']):+.3f}"
          if "debut_impulsion_dent" in a.evenements else "—")
        l("Échappement de la dent, ancre (°)", lambda a: f"{d(a.evenements['chute']):+.3f}")
        l("Recul de la roue (°)", lambda a: f"{a.mesures['recul_deg']:.3f}")
        l("Chute (°)", lambda a: f"{a.mesures['chute_deg']:.3f}")
        l("Chute à la tête (µm)",
          lambda a: f"{1000 * a.mesures['chute_deg'] * DEG * T.trace.rayon_tete:.1f}")
        l("Repos à la chute (µm)", lambda a: f"{a.mesures['repos_chute_um']:.1f}")
        l("Chemin perdu (µm)", lambda a: f"{a.mesures['chemin_perdu_um']:.1f}")
        l("Repos sur butée (µm)", lambda a: f"{a.mesures['repos_butee_um']:.1f}")
        l("Entrée de l'ellipse / ligne des centres (°)",
          lambda a: f"{a.mesures['entree_ellipse_deg']:+.2f}")
        l("Fin du dégagement, balancier (°)",
          lambda a: f"{a.mesures['fin_degagement_balancier_deg']:+.2f}")
        l("Fin de l'impulsion, balancier (°)",
          lambda a: f"{a.mesures['fin_impulsion_balancier_deg']:+.2f}")
        l("Levée du balancier (°)", lambda a: f"{a.mesures['levee_deg']:.2f}")
        return lignes


def releve(trace_ou_tables):
    T = trace_ou_tables if isinstance(trace_ou_tables, Tables) else Tables(trace_ou_tables)
    return Releve(T)


def imprimer(r):
    print("Colonnes : alternance où travaille la palette d'entrée / de sortie ; chute, repos"
          " et chemin perdu sur la palette qui reçoit la dent.")
    print(f"{'':46s} {'entrée':>10s} {'sortie':>10s}")
    for nom, e, s in r.tableau():
        print(f"{nom:46s} {e:>10s} {s:>10s}")


if __name__ == "__main__":
    import sys
    import time
    from pathlib import Path
    from .dxf_import import lire_trace
    chemin = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("simulation.dxf"))
    t0 = time.time()
    r = releve(lire_trace(chemin))
    imprimer(r)
    print(f"({time.time() - t0:.1f} s)")
