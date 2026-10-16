"""Optimisation NSGA-II du tracé.

Les variables sont celles de trace_parametrique (chacune activable, avec ses
bornes). Chaque candidat est construit, ses tables de contact calculées,
puis évalué soit en statique (environ deux secondes), soit en dynamique (une
oscillation simulée, environ une seconde de plus ; l'amplitude d'équilibre,
si elle est un objectif, en demande quatre à six).

Objectifs proposés (au choix) :
  - rendement (à maximiser) ;
  - retard dû à l'échappement à l'amplitude visée, en valeur absolue
    (à minimiser) ;
  - amplitude d'équilibre (à maximiser) : celle où l'énergie reçue par
    oscillation compense les pertes du balancier (facteur de qualité).

Contraintes proposées (chacune activable, bornes min et max) : chute, chemin
perdu, repos à la chute, levée du balancier. Elles portent sur les deux
alternances. Un candidat dont le tracé ne fonctionne pas (géométrie
impossible, dent qui ne se pose pas) est déclaré infaisable.
"""

import os
from dataclasses import dataclass, field

import numpy as np
from pymoo.core.problem import ElementwiseProblem

from .cinematique import Releve
from .parametres import Parametres
from .statique import bilan
from .tables import DEG, Tables
from .trace_parametrique import LIBELLES, VariablesTrace, construire, mesurer

OBJECTIFS = {
    "rendement": ("Rendement", +1),            # +1 : à maximiser
    "retard": ("|Retard| à l'amplitude visée (s/j)", -1),
    "amplitude": ("Amplitude d'équilibre (°)", +1),
}

CONTRAINTES = {
    "chute": ("Chute (°)", "chute_deg"),
    "chemin_perdu": ("Chemin perdu (µm)", "chemin_perdu_um"),
    "repos_chute": ("Repos à la chute (µm)", "repos_chute_um"),
    "levee": ("Levée du balancier (°)", "levee_deg"),
}

#: pas de la grille des tables pendant l'optimisation (plus grossier qu'à
#: l'affichage : le relevé des événements est affiné à part)
PAS_OPTI = 0.05 * DEG


@dataclass
class VariableOpti:
    nom: str
    actif: bool
    mini: float
    maxi: float


@dataclass
class ConfigOpti:
    variables: list
    objectifs: list = field(default_factory=lambda: ["rendement", "retard"])
    contraintes: dict = field(default_factory=dict)    # nom -> (mini ou None, maxi ou None)
    population: int = 24
    generations: int = 12
    processus: int = max(1, (os.cpu_count() or 2) - 1)
    graine: int = 1
    mode: str = "statique"


def variables_par_defaut(ref):
    """Bornes raisonnables autour du tracé de référence ; rien d'actif."""
    m = mesurer(ref)
    out = []
    for nom in LIBELLES:
        if nom in m:
            v = m[nom]
            ecart = {"course_ancre": 1.5, "entraxe_roue_ancre": 0.05,
                     "entraxe_ancre_balancier": 0.05, "rayon_ellipse": 0.05,
                     "largeur_fourchette": 0.02}[nom]
            out.append(VariableOpti(nom, False, v - ecart, v + ecart))
        elif nom.startswith("longueur"):
            out.append(VariableOpti(nom, False, -30.0, 30.0))
        else:
            out.append(VariableOpti(nom, False, -2.0, 2.0))
    return out


def evaluer(ref, par, x, noms, mode="statique", objectifs=()):
    """Évalue un candidat. Rend un dict de mesures, ou None si le tracé ne
    fonctionne pas."""
    var = VariablesTrace(**dict(zip(noms, x)))
    try:
        trace = construire(ref, var)
        rel = Releve(Tables(trace, pas=PAS_OPTI))
        res, _ = bilan(rel, par)
    except Exception:
        return None
    m = {"rendement": res["rendement"], "retard": res["retard_s_j"],
         "amplitude": res["amplitude_equilibre_deg"]}
    # palette d'entrée d'abord : c'est l'ordre des colonnes « E/S » des résultats
    for alt in sorted(rel.alternances, key=lambda a: a.palette != "entrée"):
        for nom, (_, cle) in CONTRAINTES.items():
            m.setdefault(nom, []).append(alt.mesures[cle])
    if mode == "dynamique":
        from .dynamique import Simulateur, amplitude_equilibre
        try:
            r = Simulateur(rel, par).simuler()
            if not r.valide or abs(r.bilan["travail_moteur"] - par.couple_roue * trace.pas) \
                    > 0.05 * par.couple_roue * trace.pas:
                return None                      # la roue n'a pas avancé d'une dent
            m["rendement"] = r.bilan["rendement"]
            m["retard"] = r.bilan["marche_s_j"]
            # l'amplitude d'équilibre dynamique coûte quelques oscillations :
            # seulement si elle est demandée
            m["amplitude"] = amplitude_equilibre(rel, par)[0] if "amplitude" in objectifs \
                else r.bilan["amplitude_fin_deg"]
        except Exception:
            return None
    return m


class ProblemeEchappement(ElementwiseProblem):
    """Problème NSGA-II ; défini au niveau du module pour être transmis aux
    processus de calcul."""

    def __init__(self, ref, par, config, **kw):
        self.ref, self.par, self.config = ref, par, config
        actives = [v for v in config.variables if v.actif]
        self.noms = [v.nom for v in actives]
        self.contraintes = [(k, b) for k, b in config.contraintes.items() if b is not None]
        self.n_g = max(1, sum(2 * ((b[0] is not None) + (b[1] is not None))
                              for _, b in self.contraintes))
        super().__init__(n_var=len(actives), n_obj=len(config.objectifs),
                         n_ieq_constr=self.n_g,
                         xl=np.array([v.mini for v in actives]),
                         xu=np.array([v.maxi for v in actives]), **kw)

    def _evaluate(self, x, out, *args, **kwargs):
        cfg = self.config
        m = evaluer(self.ref, self.par, x, self.noms, cfg.mode, cfg.objectifs)
        if m is None:
            out["F"] = np.full(len(cfg.objectifs), 1e6)
            out["G"] = np.full(self.n_g, 1e3)
            out["mesures"] = None
            return
        F = []
        for o in cfg.objectifs:
            v = abs(m[o]) if o == "retard" else m[o]
            F.append(-v if OBJECTIFS[o][1] > 0 else v)
        G = []
        for k, (lo, hi) in self.contraintes:
            for v in m[k]:                        # les deux alternances
                if lo is not None:
                    G.append(lo - v)
                if hi is not None:
                    G.append(v - hi)
        G += [-1.0] * (self.n_g - len(G))
        out["F"] = np.array(F)
        out["G"] = np.array(G)
        out["mesures"] = m


class Arret(Exception):
    pass


def optimiser(ref, par, config, rappel=None, arreter=lambda: False):
    """Lance NSGA-II. `rappel(génération, n_évaluations, front)` est appelé
    à chaque génération ; `arreter()` interrompt proprement.
    Rend une liste de solutions (dicts : variables, objectifs, mesures)."""
    import multiprocessing as mp
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.callback import Callback
    from pymoo.core.problem import StarmapParallelization
    from pymoo.optimize import minimize

    noms = [v.nom for v in config.variables if v.actif]
    if not noms:
        raise ValueError("aucune variable active")

    class Suivi(Callback):
        def notify(self, algorithm):
            if arreter():
                algorithm.termination.force_termination = True
            if rappel is not None:
                rappel(algorithm.n_gen, algorithm.evaluator.n_eval, _front(algorithm.opt, noms))

    # « fork » : les processus héritent du tracé de référence sans le
    # recharger (Linux ; sous Windows il faudrait « spawn »)
    pool = mp.get_context("fork").Pool(config.processus) if config.processus > 1 else None
    try:
        kw = {"elementwise_runner": StarmapParallelization(pool.starmap)} if pool else {}
        probleme = ProblemeEchappement(ref, par, config, **kw)
        res = minimize(probleme, NSGA2(pop_size=config.population),
                       ("n_gen", config.generations), seed=config.graine,
                       callback=Suivi(), verbose=False)
    finally:
        if pool:
            pool.close()
            pool.join()
    return _front(res.opt if res.opt is not None else [], noms), noms


def _front(pop, noms):
    sol = []
    if pop is None:
        return sol
    for ind in pop:
        m = ind.get("mesures")
        if m is None or not isinstance(m, dict):
            continue
        sol.append({"variables": dict(zip(noms, map(float, ind.X))),
                    "mesures": m, "faisable": bool(np.all(ind.G <= 0))})
    return sol
