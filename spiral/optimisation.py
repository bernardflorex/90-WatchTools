# optimisation.py
"""
Branchement NSGA-II (pymoo) sur le vecteur de conception de parametrage.py,
avec comme objectif PRINCIPAL le vrai critère de marche diurne
(isochronisme.critere_dispersion_marche), pas un proxy géométrique -- c'est
le changement d'architecture demandé par rapport à spiropt_v1/optimization.py
(qui optimisait symétrie/décentrage, isochronisme.py n'y servait qu'à
l'affichage a posteriori).

Vecteur de conception (ordre fixe, voir vecteur_depuis_params) :
    [epaisseur_noeuds (Me), delta_r_noeuds (Mr),
     grossmann_angle_deg, grossmann_offsets (Ng),
     phillips_angle_deg,  phillips_offsets (Np)]
Me, Mr, Ng, Np sont fixés par la structure de `base` (ParametresSpiral) au
moment de la construction du problème -- ce sont des choix de discrétisation
de la CONCEPTION (nombre de degrés de liberté de forme), pas des variables
d'optimisation elles-mêmes. n_pas (densité de maillage physique) et n_points
(degré des Béziers Grossmann/Phillips) ne sont pas non plus des variables
libres, pour la même raison : pymoo a besoin d'un vecteur de longueur fixe
sur toute la durée d'une même optimisation.
"""
from dataclasses import replace
import numpy as np
from joblib import Parallel, delayed
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize as pymoo_minimize

from .echappement import ReglageEchappement
from .thermique import ReglageThermique
from .parametrage import ParametresSpiral, construire_geometrie, rayons_jonction
from .isochronisme import analyser_marches
from .presets import BornesOptimisation, HyperparametresGA

# Pénalités : très défavorables au tri de dominance NSGA-II, sans être des
# NaN/inf qui casseraient pymoo. Appliquées quand la géométrie elle-même est
# invalide (repli, non-monotonie) -- cf. parametrage.construire_geometrie --
# plutôt que de laisser une exception remonter et interrompre l'optimisation.
PENALITE_MARCHE_S_PAR_J = 1.0e5
PENALITE_JEU_MM = -1.0
PENALITE_ECART_S_PAR_J = 1.0e5


def offsets_exterieurs(p: ParametresSpiral) -> np.ndarray:
    """Variables de forme de la zone terminale EXTÉRIEURE, quel que soit le
    mode : offsets de la courbe de Phillips sur un spiral libre, offsets de la
    courbe en S sur un spiral à raquette. Le mode étant figé pour toute la
    durée d'une optimisation, le vecteur de conception garde bien une longueur
    constante, comme l'exige pymoo."""
    return np.asarray(p.raquette_offsets if p.raquette_active else p.phillips_offsets,
                      dtype=float)


def _remplacer_exterieur(base: ParametresSpiral, angle_deg: float, offsets: np.ndarray):
    if base.raquette_active:
        return {"raquette_transition_angle_deg": angle_deg, "raquette_offsets": offsets}
    return {"phillips_angle_deg": angle_deg, "phillips_offsets": offsets}


def vecteur_depuis_params(p: ParametresSpiral) -> np.ndarray:
    return np.concatenate([
        np.asarray(p.epaisseur_noeuds, dtype=float),
        np.asarray(p.delta_r_noeuds, dtype=float),
        [float(p.grossmann_angle_deg)],
        np.asarray(p.grossmann_offsets, dtype=float),
        [float(p.angle_exterieur_deg)],
        offsets_exterieurs(p),
    ])


def params_depuis_vecteur(x: np.ndarray, base: ParametresSpiral) -> ParametresSpiral:
    n_e, n_r = len(base.epaisseur_noeuds), len(base.delta_r_noeuds)
    n_g, n_p = len(base.grossmann_offsets), len(offsets_exterieurs(base))
    x = np.asarray(x, dtype=float)
    i = 0
    epaisseur_noeuds, i = x[i:i + n_e], i + n_e
    delta_r_noeuds, i = x[i:i + n_r], i + n_r
    grossmann_angle_deg, i = float(x[i]), i + 1
    grossmann_offsets, i = x[i:i + n_g], i + n_g
    angle_exterieur_deg, i = float(x[i]), i + 1
    offsets_ext, i = x[i:i + n_p], i + n_p

    return replace(
        base,
        epaisseur_noeuds=epaisseur_noeuds, delta_r_noeuds=delta_r_noeuds,
        grossmann_angle_deg=grossmann_angle_deg, grossmann_offsets=grossmann_offsets,
        **_remplacer_exterieur(base, angle_exterieur_deg, offsets_ext),
    )


def construire_bornes(base: ParametresSpiral, bornes: BornesOptimisation):
    """Bornes NSGA-II autour du point de référence `base`, dans l'ordre de
    vecteur_depuis_params. Les épaisseurs sont bornées en relatif (spread_
    epaisseur_relatif) ; les offsets Grossmann/Phillips en fraction de la
    largeur radiale RÉELLE de leur zone respective (voir
    parametrage.rayons_jonction -- une borne fixe en mm ne s'adapte pas au
    calibre, cf. presets.BornesOptimisation pour la justification)."""
    ep = np.asarray(base.epaisseur_noeuds, dtype=float)
    ep_low = ep * (1.0 - bornes.spread_epaisseur_relatif)
    ep_high = ep * (1.0 + bornes.spread_epaisseur_relatif)

    dr = np.asarray(base.delta_r_noeuds, dtype=float)
    dr_low = dr - bornes.spread_delta_r_mm
    dr_high = dr + bornes.spread_delta_r_mm

    _, r_body_start, r_body_end = rayons_jonction(base)
    largeur_gross = max(r_body_start - base.r_min, 1e-6)
    largeur_phil = max(base.r_max - r_body_end, 1e-6)

    g_off = np.asarray(base.grossmann_offsets, dtype=float)
    g_off_low = g_off - bornes.spread_grossmann_relatif * largeur_gross
    g_off_high = g_off + bornes.spread_grossmann_relatif * largeur_gross

    # Zone terminale extérieure : courbe de Phillips (spiral libre) ou courbe
    # en S de la raquette. Mêmes bornes angulaires et même largeur radiale de
    # référence dans les deux cas -- les deux courbes partent du corps et
    # arrivent à r_max, seule leur condition de tangence à l'arrivée diffère.
    p_off = offsets_exterieurs(base)
    p_off_low = p_off - bornes.spread_phillips_relatif * largeur_phil
    p_off_high = p_off + bornes.spread_phillips_relatif * largeur_phil

    xl = np.concatenate([ep_low, dr_low, [bornes.angle_grossmann_min_deg], g_off_low,
                          [bornes.angle_phillips_min_deg], p_off_low])
    xu = np.concatenate([ep_high, dr_high, [bornes.angle_grossmann_max_deg], g_off_high,
                          [bornes.angle_phillips_max_deg], p_off_high])
    return xl, xu


PAS_DISPERSION_DEG_DEFAUT = 5.0


def evaluer_candidat(x_libre: np.ndarray, base: ParametresSpiral,
                      amplitudes_analyse_deg: np.ndarray, pas_balayage_deg: float = 20.0,
                      pas_dispersion_deg: float = PAS_DISPERSION_DEG_DEFAUT,
                      reglage_echappement: ReglageEchappement = None,
                      reglage_thermique: ReglageThermique = None):
    """Évalue un vecteur de conception : retourne (dispersion_marche_s_par_j,
    jeu_min_mm, decentrage_mm, ecart_compensation_s_par_j, epaisseur_oxyde_um,
    dvh_s_par_j). C'est la fonction appelée en parallèle (joblib) pour chaque
    individu de la population à chaque génération.

    dvh_s_par_j : critère COSC, différence de marche entre la position 6H et la
    position horizontale, prise AU PIRE CAS sur la plage d'amplitude (comme la
    dispersion). Comme la marche en position verticale suit exactement
    une loi marche(θ) = A·cos(θ − φ) (vérifié à 7e-16 près sur 72 directions),
    le DVH vaut A·cos(90° − φ) : il dépend donc AUSSI de la phase φ, que la
    dispersion seule ignore. C'est pourquoi il est retourné à part.

    ecart_compensation_s_par_j : écart signé entre l'avance que le spiral
    produit de lui-même à l'amplitude basse (isochronisme propre) et celle
    qu'il faudrait pour annuler le retard d'échappement. Zéro = compensation
    exacte. Cette grandeur ne peut PAS être lue sur la dispersion entre
    positions verticales : elle agit identiquement dans toutes les positions,
    horizontales comprises, donc elle ne crée aucune dispersion entre
    positions -- c'est un objectif indépendant, pas un sous-produit du
    précédent (cf. echappement.py).

    Deux points importants (diagnostic 2026-09-23, voir aussi
    physique.balayage_angulaire) :
      - Le jeu inter-spires est le MINIMUM sur TOUT le balayage angulaire
        (angles_deg du diag_orbite, résolution pas_balayage_deg), pas
        seulement à 0/+max/-max comme avant : le jeu minimal peut survenir à
        un angle intermédiaire, une vérification aux seuls extrêmes pouvait
        laisser passer une vraie collision.
      - La dispersion de marche est évaluée sur une grille d'amplitudes
        RE-DENSIFIÉE à pas_dispersion_deg (5° par défaut), PAS sur
        amplitudes_analyse_deg tel quel : ce dernier ne sert qu'à fixer la
        plage [min, max] (piloté par l'utilisateur pour l'affichage, où un
        pas grossier est confortable) ; l'utiliser tel quel comme grille de
        l'OBJECTIF ferait courir le risque de manquer le pire cas entre deux
        points grossièrement échantillonnés -- exactement l'écart observé
        entre la dispersion rapportée par l'optimisation et celle relue sur
        le graphe à une résolution plus fine.
    """
    p = params_depuis_vecteur(x_libre, base)
    geom, diag = construire_geometrie(p, return_diagnostics=True,
                                       reglage_thermique=reglage_thermique)

    if not (diag["geometrie_valide"] and diag["grossmann_monotone"] and diag["phillips_monotone"]):
        return (PENALITE_MARCHE_S_PAR_J, PENALITE_JEU_MM, 0.0, PENALITE_ECART_S_PAR_J,
                0.0, PENALITE_MARCHE_S_PAR_J)

    amp_min, amp_max = float(np.min(amplitudes_analyse_deg)), float(np.max(amplitudes_analyse_deg))
    amplitudes_denses = np.arange(amp_min, amp_max + 1e-6, pas_dispersion_deg)

    resultat = analyser_marches(
        geom, p.r_min, p.freq_cible, p.inertie_balancier,
        amplitudes_deg=amplitudes_denses, reglage_echappement=reglage_echappement,
        pas_balayage_deg=pas_balayage_deg, jeu_goupilles_um=p.raquette_jeu_um,
    )
    diag_orbite = resultat["diag_orbite"]
    dispersion = resultat["dispersion_gravite"]
    ecart = resultat["ecart_compensation"]
    dvh = resultat["dvh_max"]
    decentrage = float(np.hypot(np.mean(diag_orbite["x_cg_mm"]), np.mean(diag_orbite["y_cg_mm"])))
    jeu = diag_orbite["jeu_min_mm"]

    if not np.isfinite(dispersion) or dispersion > PENALITE_MARCHE_S_PAR_J:
        dispersion = PENALITE_MARCHE_S_PAR_J
    if not np.isfinite(jeu):
        jeu = PENALITE_JEU_MM
    if not np.isfinite(ecart):
        ecart = PENALITE_ECART_S_PAR_J
    if not np.isfinite(dvh) or dvh > PENALITE_MARCHE_S_PAR_J:
        dvh = PENALITE_MARCHE_S_PAR_J

    return (float(dispersion), float(jeu), decentrage, float(ecart),
            float(geom.epaisseur_oxyde * 1000.0), float(dvh))


def sonder_faisabilite(base: ParametresSpiral, bornes: BornesOptimisation,
                        amplitudes_analyse_deg: np.ndarray,
                        n_echantillons: int = 200, seed: int = 0,
                        pas_balayage_deg: float = 20.0, n_jobs: int = -1,
                        n_classement: int = 10,
                        reglage_echappement: ReglageEchappement = None,
                        reglage_thermique: ReglageThermique = None) -> dict:
    """Tirage aléatoire uniforme dans les bornes, pour trouver un point de
    départ "le moins mauvais possible" AVANT de lancer NSGA-II -- même
    logique que spiropt_v1/sonder_espace.py. Nécessaire en pratique : partir
    d'un vecteur arbitraire (ex. épaisseur uniforme, aucun offset) est
    presque toujours largement infaisable (spires qui se chevauchent à
    grande amplitude), et NSGA-II part alors d'une population quasi
    entièrement infaisable, ce qui retarde ou empêche la convergence.

    AUCUN seuil de jeu minimal n'est appliqué ici (volontairement -- voir
    échange du 2026-09-22) : le sondage n'a pas vocation à trouver une
    solution faisable, seulement à classer les tirages du moins mauvais au
    plus mauvais pour fournir un bon point de départ à l'étape suivante
    (Optimisation NSGA-II, qui elle applique la vraie contrainte). Filtrer
    ici sur jeu_min_requis n'aurait fait que masquer des points de départ
    par ailleurs utiles (dispersion de marche déjà faible, par exemple)
    simplement parce qu'ils n'atteignaient pas encore le jeu requis.
    """
    xl, xu = construire_bornes(base, bornes)
    rng = np.random.default_rng(seed)
    X = rng.uniform(xl, xu, size=(n_echantillons, len(xl)))

    # La géométrie de DÉPART elle-même est évaluée, en plus des tirages
    # (correction 2026-09-24). Les bornes étant centrées sur elle, un tirage
    # uniforme dans la boîte s'en écarte forcément dans toutes les directions
    # à la fois : il est courant qu'AUCUN des tirages ne fasse mieux que le
    # point central, et l'optimisation démarrait alors d'un warm start pire
    # que la géométrie que l'utilisateur avait sous les yeux. Constaté
    # concrètement : partie du meilleur de 120 tirages, la population restait
    # infaisable sur le jeu après 60 générations, alors que partie de la
    # géométrie neutre elle-même, elle devenait entièrement faisable
    # (jeu ET compensation d'échappement) dès la 20e.
    X = np.vstack([np.clip(vecteur_depuis_params(base), xl, xu), X])

    parallel = Parallel(n_jobs=n_jobs, backend="loky")
    resultats = parallel(
        delayed(evaluer_candidat)(x, base, amplitudes_analyse_deg, pas_balayage_deg,
                                   PAS_DISPERSION_DEG_DEFAUT, reglage_echappement)
        for x in X
    )
    dispersion = np.array([r[0] for r in resultats])
    jeu = np.array([r[1] for r in resultats])
    ecart = np.array([r[3] for r in resultats])
    oxyde = np.array([r[4] for r in resultats])
    dvh = np.array([r[5] for r in resultats])
    ordre = np.argsort(jeu)[::-1]  # du meilleur (jeu le plus grand) au pire
    classement = ordre[:min(n_classement, len(ordre))]

    return {
        "X": X, "dispersion": dispersion, "jeu_mm": jeu, "ecart_compensation": ecart,
        "oxyde_um": oxyde, "dvh": dvh,
        "ordre": ordre,
        "classement": [
            {"x": X[i], "dispersion": float(dispersion[i]), "jeu_mm": float(jeu[i]),
             "ecart_compensation": float(ecart[i]), "oxyde_um": float(oxyde[i]),
             "dvh": float(dvh[i])}
            for i in classement
        ],
        "x_meilleur": X[ordre[0]], "jeu_meilleur_mm": float(jeu[ordre[0]]),
        "dispersion_meilleur": float(dispersion[ordre[0]]),
    }


class ProblemeSpiral(Problem):
    """Problème pymoo à 2 objectifs :
        F[:,0] = critère chronométrique (s/j), à MINIMISER : dispersion de
                 marche entre positions verticales PLUS |DVH|, le critère de
                 contrôle COSC (marche en position 6H moins marche à plat).
        F[:,1] = -jeu inter-spires minimal (mm), plafonné à jeu_max_cap pour
                 éviter une incitation excessive au-delà d'une marge de
                 sécurité déjà suffisante (même logique que spiropt_v1).
    et 2 contraintes :
        G[:,0] = jeu_min_requis - jeu_mesuré <= 0
        G[:,1] = |écart de compensation d'échappement| - tolérance <= 0

    Pourquoi la compensation d'échappement est une CONTRAINTE et pas un 3e
    objectif : ce n'est pas une grandeur à minimiser mais une CIBLE à
    atteindre (+retard_s_par_j à l'amplitude basse). En faire un objectif
    aurait fait passer le front de Pareto de 2 à 3 dimensions, le diluant
    fortement à population constante, pour un arbitrage qui n'a pas lieu
    d'être : au-delà de la tolérance, un spiral qui sur-compense n'est pas
    « meilleur », il est faux dans l'autre sens. NSGA-II minimisant d'abord
    la violation de contrainte, la population CONVERGE VERS la cible avant
    d'arbitrer dispersion contre jeu -- exactement le comportement demandé.
    Si la cible se révèle inatteignable, l'algorithme dégrade proprement en
    s'en approchant au mieux au lieu de tout rejeter.

    Le décentrage moyen (diagnostic COSC) est calculé mais volontairement
    PAS ajouté comme objectif : sur un spiral correctement Grossmann/
    Phillips, il est très fortement corrélé à la dispersion de marche
    (les deux viennent de la même orbite de CG) -- l'ajouter dilue le front
    de Pareto sans apporter d'arbitrage réellement indépendant. Reste
    accessible via evaluer_candidat pour analyse a posteriori.

    Pourquoi le DVH s'AJOUTE à la dispersion au lieu de former un 3e objectif
    (choix de l'utilisateur, 2026-09-29) : les deux se lisent en s/j et vont
    dans le même sens -- aucun arbitrage à exposer, donc rien à gagner à
    diluer le front de Pareto d'une dimension de plus. Ils ne sont pourtant
    pas redondants. La marche en position verticale suit exactement
    marche(θ) = A·cos(θ − φ) ; la dispersion ne voit que le MODULE A, alors
    que le DVH = A·cos(90° − φ) voit aussi la PHASE. À dispersion égale, la
    somme préfère donc les géométries dont la phase reporte l'effet résiduel
    de gravité sur 3H et 9H plutôt que sur la position 6H, celle que le COSC
    note. C'est précisément le degré de liberté que la dispersion seule
    ignorait. Le poids reste exposé (`poids_dvh`) pour pouvoir revenir à
    l'ancien comportement en le mettant à zéro.
    """

    def __init__(self, base: ParametresSpiral, bornes: BornesOptimisation,
                 amplitudes_analyse_deg: np.ndarray, jeu_min_requis: float,
                 jeu_max_cap: float = 0.080, pas_balayage_deg: float = 20.0, n_jobs: int = -1,
                 reglage_echappement: ReglageEchappement = None,
                 reglage_thermique: ReglageThermique = None,
                 poids_dvh: float = 1.0):
        self.base = base
        self.poids_dvh = float(poids_dvh)
        self.amplitudes_analyse_deg = np.asarray(amplitudes_analyse_deg, dtype=float)
        self.jeu_min_requis = jeu_min_requis
        self.jeu_max_cap = jeu_max_cap
        self.pas_balayage_deg = pas_balayage_deg
        self.reglage_echappement = (reglage_echappement if reglage_echappement is not None
                                     else ReglageEchappement())
        self.reglage_thermique = reglage_thermique
        self.parallel = Parallel(n_jobs=n_jobs, backend="loky")

        xl, xu = construire_bornes(base, bornes)
        n_contraintes = 2 if self.reglage_echappement.actif else 1
        self.contrainte_oxyde = (reglage_thermique is not None
                                  and reglage_thermique.contrainte_oxyde_active())
        if self.contrainte_oxyde:
            n_contraintes += 1
        super().__init__(n_var=len(xl), n_obj=2, n_ieq_constr=n_contraintes, xl=xl, xu=xu)

    def _evaluate(self, X, out, *args, **kwargs):
        resultats = self.parallel(
            delayed(evaluer_candidat)(x, self.base, self.amplitudes_analyse_deg,
                                       self.pas_balayage_deg, PAS_DISPERSION_DEG_DEFAUT,
                                       self.reglage_echappement, self.reglage_thermique)
            for x in X
        )
        dispersion = np.array([r[0] for r in resultats])
        jeu = np.array([r[1] for r in resultats])
        ecart = np.array([r[3] for r in resultats])
        oxyde = np.array([r[4] for r in resultats])
        dvh = np.array([r[5] for r in resultats])

        critere = dispersion + self.poids_dvh * np.abs(dvh)
        out["F"] = np.column_stack([critere, -np.minimum(jeu, self.jeu_max_cap)])
        contraintes = [self.jeu_min_requis - jeu]
        if self.reglage_echappement.actif:
            contraintes.append(np.abs(ecart) - self.reglage_echappement.tolerance_s_par_j)
        if self.contrainte_oxyde:
            # Une seule oxydation doit servir toute une série : l'épaisseur
            # requise par chaque candidat doit tomber dans la fenêtre que le
            # procédé sait tenir. Écart au CENTRE de la fenêtre, comparé à sa
            # demi-largeur -- formulation symétrique, qui ne privilégie aucun
            # des deux bords.
            r = self.reglage_thermique
            centre = 0.5 * (r.oxyde_min_um + r.oxyde_max_um)
            demi = 0.5 * (r.oxyde_max_um - r.oxyde_min_um)
            contraintes.append(np.abs(oxyde - centre) - demi)
        out["G"] = np.column_stack(contraintes)


def construire_population_initiale(x0, xl, xu, pop_size, fraction_locale=0.7,
                                    rayon_relatif=0.15, seed=42):
    """Warm-start : une fraction de la population tirée localement autour de
    x0 (déjà connu "pas trop mauvais"), le reste uniformément sur toute la
    boîte pour garder de la diversité génétique -- repris tel quel de
    spiropt_v1/optimization.py (logique déjà validée)."""
    x0 = np.clip(np.asarray(x0, dtype=float), xl, xu)
    largeur = xu - xl
    bas_local = np.clip(x0 - rayon_relatif * largeur, xl, xu)
    haut_local = np.clip(x0 + rayon_relatif * largeur, xl, xu)

    rng = np.random.default_rng(seed)
    n_local = int(round(pop_size * fraction_locale))
    n_global = pop_size - n_local

    pop_locale = rng.uniform(bas_local, haut_local, size=(n_local, len(xl))) if n_local > 0 else np.empty((0, len(xl)))
    pop_globale = rng.uniform(xl, xu, size=(n_global, len(xl))) if n_global > 0 else np.empty((0, len(xl)))
    population = np.vstack([pop_locale, pop_globale])
    if len(population) > 0:
        population[0] = x0
    return population


def lancer_optimisation(base: ParametresSpiral, bornes: BornesOptimisation,
                         hyperparams: HyperparametresGA, jeu_min_requis: float,
                         amplitudes_analyse_deg: np.ndarray = None,
                         warm_start: bool = True, x0: np.ndarray = None, callback=None,
                         pas_balayage_deg: float = 20.0,
                         reglage_echappement: ReglageEchappement = None,
                         reglage_thermique: ReglageThermique = None,
                         poids_dvh: float = 1.0):
    """Lance NSGA-II et retourne le résultat pymoo (res.X, res.F : ensemble
    de solutions non-dominées). Fonction pure, sans dépendance Qt -- la
    future GUI (Phase 3) l'appellera depuis un QThread, un script ou un
    notebook peuvent l'appeler directement.

    x0 : point de départ du warm-start (typiquement `sonder_faisabilite(...)
    ["x_meilleur"]`). Si None, utilise `base` lui-même (raisonnable
    seulement si `base` est déjà proche de la faisabilité)."""
    if amplitudes_analyse_deg is None:
        amplitudes_analyse_deg = np.arange(200.0, 341.0, 20.0)  # cf. spiropt_v1/isochronisme.AMPLITUDES_ANALYSE_DEFAUT

    # pas_balayage_deg est transmis explicitement : il était auparavant figé à
    # la valeur par défaut de ProblemeSpiral, alors que la vérification et les
    # tracés côté interface utilisaient le réglage de l'utilisateur -- deux
    # réglages différents pour la même grandeur, donc deux chiffres différents
    # sans raison (diagnostic 2026-09-26).
    probleme = ProblemeSpiral(base, bornes, amplitudes_analyse_deg, jeu_min_requis,
                               pas_balayage_deg=pas_balayage_deg,
                               reglage_echappement=reglage_echappement,
                               reglage_thermique=reglage_thermique,
                               poids_dvh=poids_dvh)

    if warm_start:
        if x0 is None:
            x0 = vecteur_depuis_params(base)
        sampling = construire_population_initiale(
            x0, probleme.xl, probleme.xu, hyperparams.pop_size,
            fraction_locale=hyperparams.fraction_locale_warmstart,
            rayon_relatif=hyperparams.rayon_relatif_warmstart, seed=hyperparams.seed,
        )
        algorithme = NSGA2(pop_size=hyperparams.pop_size, sampling=sampling)
    else:
        algorithme = NSGA2(pop_size=hyperparams.pop_size)

    kwargs = {"callback": callback} if callback is not None else {}
    return pymoo_minimize(
        probleme, algorithme, ("n_gen", hyperparams.n_gen),
        seed=hyperparams.seed, verbose=False, **kwargs,
    )
