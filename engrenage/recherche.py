"""Recherche d'optima de profils en arcs de cercle : CMA-ES avec redémarrages.

Pilote de `optim.Probleme` (voir optim.py pour la fonction objectif et le choix des variables) :

1. plusieurs départs : le profil courant, puis des profils tirés au hasard de plus en plus loin
   (on n'a pas de raison de croire que l'optimum est voisin du profil de départ : d'autres
   solutions stables existent plus loin), chacun raffiné par un CMA-ES ; la population double un
   départ sur deux ;
2. les profils qui violent encore une contrainte (interférence entre dents, jeu, rayon de tête)
   sont réparés en libérant le demi-profil arrière des deux roues : le demi-profil menant, optimisé,
   reste tel quel, on modifie l'autre flanc des dents (dents dissymétriques) ;
3. la liste finale ne garde que des optima distincts (écart de profil d'au moins `separation_um`),
   les plus faibles J d'abord ; ceux qui restent non valables sont signalés.
"""

import math
import time
from dataclasses import dataclass, field, replace

import cma
import numpy as np

import optim

TOLERANCES = {"jeu_min": 1e-3, "jeu_max": 1e-3, "r_tete": 1e-3, "interference": optim.EPS_INTERF,
              "ecart_tete": 1e-4}  # dépassements admis (°, °, fraction, mm, mm)
LIBELLES = {"jeu_min": "jeu insuffisant (blocage)", "jeu_max": "jeu excessif", "r_tete": "rayon de tête trop petit",
            "interference": "interférence entre dents", "ecart_tete": "têtes avant et arrière différentes"}


@dataclass
class Options:
    departs: int = 6  # nombre de départs (le premier est le profil courant)
    population: int = 24  # taille de la population de CMA-ES du premier départ
    evaluations: int = 12000  # évaluations par départ
    etendue_um: float = 10.0  # écart typique de profil exploré autour du départ (µm)
    graine: int = 1
    dissymetrie: bool = True  # réparer par le demi-profil arrière (dents dissymétriques)
    separation_um: float = 3.0  # écart minimal entre deux optima de la liste
    conserver: int = 10  # taille maximale de la liste


@dataclass
class Optimum:
    j: float
    x: np.ndarray
    mobiles: list  # deux arcs.Mobile (points refermés ; points_dos si dissymétrique)
    rapport: list  # optim.Probleme.rapport : un dictionnaire par entraxe
    violations: dict  # dépassements restants (clé -> valeur > tolérance)
    depart: int
    dissymetrique: bool
    criteres: dict = field(default_factory=dict)  # valeur des quatre critères / valeur de départ

    @property
    def valable(self):
        return not self.violations


def _violations(pb, x):
    s, codes = pb.evaluer(x)
    if codes[0]:
        return {"invalide": float(codes[0])}, s
    v = pb.violations(s)
    return {k: float(a[0]) for k, a in v.items() if a[0] > TOLERANCES[k]}, s


def _cma(pb, x0, sigma0, popsize, budget, graine, rappel, arret, fonction=None, arret_si=None):
    """Un CMA-ES par lots ; retourne (meilleur x, meilleur J, évaluations)."""
    fonction = fonction or pb
    es = cma.CMAEvolutionStrategy(np.asarray(x0, dtype=float), sigma0,
                                  {"popsize": popsize, "seed": graine, "verbose": -9, "maxfevals": budget})
    n, meilleur, xm = 0, math.inf, np.asarray(x0, dtype=float)
    while not es.stop() and n < budget and not (arret and arret()):
        x = np.array(es.ask())
        f = fonction(x)
        es.tell(list(x), list(f))
        n += len(x)
        i = int(np.argmin(f))
        if f[i] < meilleur:
            meilleur, xm = float(f[i]), x[i].copy()
        if rappel:
            rappel(n, meilleur)
        if arret_si and arret_si(meilleur, xm):
            break
    return xm, meilleur, n


def _tirage(pb, rng, amplitude_um, echelle, essais=40):
    """Vecteur de variables valide dont le profil s'écarte d'environ `amplitude_um` du départ."""
    for _ in range(essais):
        x = rng.normal(size=(32, pb.n_var)) * (amplitude_um / echelle)
        _, codes = pb.evaluer(x, avec_jeu=False)
        ok = np.flatnonzero(codes == 0)
        if len(ok):
            return x[ok[0]]
    return None


def _reparer(pb_dos_modele, opt_x, options, rappel, arret, graine):
    """Répare un optimum en libérant les demi-profils arrière ; retourne le nouveau x complet."""
    pb = pb_dos_modele
    e = pb.echelle_um()
    reg = 1e-4  # préfère la dent la plus proche de la symétrie

    def fonction(x):
        return pb(x) + reg * np.mean(np.asarray(x) ** 2, axis=1)

    def bon(_j, x):
        v, _ = _violations(pb, x)
        return not v

    x, _, n = _cma(pb, np.zeros(pb.n_var), max(0.3, options.etendue_um / max(e, 1e-9)), options.population,
                   options.evaluations // 2, graine, rappel, arret, fonction, arret_si=bon)
    return x, n


def rechercher(roues, mobiles, mu, da_neg, da_pos, objectif, options=None, libres=None, rappel=None, arret=None):
    """Cherche des optima. Retourne (liste d'`Optimum` valables d'abord, statistiques).

    `rappel(dict)` reçoit l'avancement (départ, évaluations, meilleur J, phase) ; `arret()` renvoie
    True pour interrompre (les résultats déjà obtenus sont conservés).
    """
    options = options or Options()
    t0 = time.perf_counter()
    # calibrage différé (calibrer=False) : le profil de départ (x=0) peut être invalide (SWTE,
    # entraxe, période) sans faire échouer toute la recherche, contrairement à un usage direct de
    # Probleme -- on cherche alors un profil valide proche pour servir de départ n°1 et de référence
    pb = optim.Probleme(roues, mobiles, mu, da_neg, da_pos, libres, objectif, calibrer=False)
    echelle = max(pb.echelle_um(), 1e-6)
    rng = np.random.default_rng(options.graine)
    trouves, total = [], 0
    stats = {"evaluations": 0, "departs": 0, "repares": 0, "temps": 0.0, "depart_invalide": None}

    x_ref = np.zeros(pb.n_var)
    try:
        pb.f0 = pb._criteres_ref()
    except ValueError as e:
        for et in (options.etendue_um, 10 * options.etendue_um, 100 * options.etendue_um):
            x_ref = _tirage(pb, rng, et, echelle)
            if x_ref is not None:
                break
        if x_ref is None:  # aucun profil valide trouvé même très loin : rien à optimiser
            raise ValueError(str(e)) from e
        pb.f0 = pb._criteres_ref(x_ref)
        stats["depart_invalide"] = str(e)

    if isinstance(objectif.jeu_max, str):  # "auto" : jeu maximal = celui du profil de référence + marge
        objectif = replace(objectif, jeu_max=pb.rapport(x_ref)[2]["jeu_max"] + 0.02)
        pb.objectif = objectif

    def info(**kw):
        if rappel:
            rappel(kw)

    for d in range(options.departs):
        if arret and arret():
            break
        if d == 0:
            x0 = x_ref  # le profil de départ, ou le profil valide trouvé à proximité s'il ne l'était pas
        else:  # départs de plus en plus éloignés (jusqu'à 8 fois l'étendue)
            x0 = _tirage(pb, rng, options.etendue_um * min(8.0, 1.0 + d), echelle)
            if x0 is None:
                continue
        pop = options.population * (2 if d % 2 == 1 else 1)
        base = total
        x, j, n = _cma(pb, x0, options.etendue_um / echelle, pop, options.evaluations, options.graine + d,
                       lambda n, m, d=d, base=base: info(phase="recherche", depart=d + 1, departs=options.departs,
                                                         evaluations=base + n, meilleur=m), arret)
        total += n
        stats["departs"] += 1
        v, s = _violations(pb, x)
        if options.dissymetrie and v and "invalide" not in v and not (arret and arret()):
            m0 = pb.mobiles(x)  # profil menant optimisé, sans arrière propre
            pb_dos = optim.Probleme(roues, m0, mu, da_neg, da_pos, libres, objectif, dos_libre=True,
                                    avant_libre=False)
            pb_dos.f0 = pb.f0  # même normalisation que la recherche
            xd, n2 = _reparer(pb_dos, x, options, lambda n, m, d=d, base=total: info(
                phase="réparation", depart=d + 1, departs=options.departs, evaluations=base + n, meilleur=m),
                              arret, options.graine + 100 + d)
            total += n2
            v2, s2 = _violations(pb_dos, xd)
            if not v2 or sum(v2.values()) < sum(v.values()):
                pb_use, x_use, v, stats["repares"] = pb_dos, xd, v2, stats["repares"] + int(not v2)
            else:
                pb_use, x_use = pb, x
        else:
            pb_use, x_use = pb, x
        rapport = pb_use.rapport(x_use) if "invalide" not in v else []
        mob = pb_use.mobiles(x_use)
        crit = {}
        if rapport:
            s_use, _ = pb_use.evaluer(x_use)
            crit = {c: float(pb_use.critere(s_use, c)[0] / pb.f0[c]) for c in optim.CRITERES}
            j = float(pb_use.valeur(s_use, np.zeros(1, dtype=np.int64))[0])
        trouves.append(Optimum(j, x_use, mob, rapport, v, d + 1, not all(m.symetrique for m in mob), crit))

    # liste finale : valables d'abord, optima distincts
    trouves.sort(key=lambda o: (not o.valable, o.j))
    liste = []
    for o in trouves:
        if len(liste) >= options.conserver:
            break
        if all(_distance_um(o, p, roues) >= options.separation_um for p in liste):
            liste.append(o)
    stats["evaluations"], stats["temps"] = total, time.perf_counter() - t0
    return liste, stats


def _distance_um(a, b, roues):
    """Écart (µm) entre deux optima : plus grand écart entre demi-profils menants (échantillonnés en arcs)."""
    pire = 0.0
    for ma, mb in zip(a.mobiles, b.mobiles):
        pa, pb_ = ma.echantillonner(160), mb.echantillonner(160)
        d = np.hypot(pa[:, None, 0] - pb_[None, :, 0], pa[:, None, 1] - pb_[None, :, 1])
        pire = max(pire, d.min(axis=1).max(), d.min(axis=0).max())
        if not (ma.symetrique and mb.symetrique):
            pa, pb_ = ma.echantillonner_dos(160), mb.echantillonner_dos(160)
            d = np.hypot(pa[:, None, 0] - pb_[None, :, 0], pa[:, None, 1] - pb_[None, :, 1])
            pire = max(pire, d.min(axis=1).max(), d.min(axis=0).max())
    return 1e3 * pire
