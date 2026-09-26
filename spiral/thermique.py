# thermique.py
"""
Compensation thermique d'un spiral silicium par couche d'oxyde.

Le problème. Le silicium s'assouplit fortement en chauffant : son coefficient
thermoélastique (1/E)·dE/dT vaut environ -60 ppm/K. Un spiral silicium nu
retarde donc d'environ 3.8 s/j/°C -- plus de six fois le critère COSC (0.6
s/j/°C), et 115 s/j d'écart entre 8 et 38 °C. L'oxyde thermique de silicium a
le comportement inverse (+215 ppm/K) : une couche de quelques microns, crûe sur
toutes les faces, ramène le coefficient du composite à zéro.

Pourquoi si peu d'oxyde suffit : en flexion, la contribution d'une couche à la
rigidité varie comme le cube de sa distance à la fibre neutre. L'oxyde étant
justement à la périphérie de la section, il pèse très au-dessus de sa part
volumique -- 3.7 µm d'oxyde sur une lame de 31 µm portent déjà un tiers de la
rigidité de flexion.

Ce que ce module calcule :
  - la rigidité de flexion du composite silicium/oxyde à une température
    donnée, section par section ;
  - l'épaisseur d'oxyde qui annule la variation de marche sur la plage de
    température visée ;
  - la marche résiduelle en fonction de la température ;
  - la sensibilité du résultat à chaque coefficient matériau, puisque ceux-ci
    sont des valeurs de littérature et non des mesures.

Couplage avec le reste de la conception. L'épaisseur d'oxyde requise dépend du
PROFIL d'épaisseur de la lame, que l'optimiseur fait justement varier : la part
de rigidité portée par l'oxyde passe de 48% sur une section de 20 µm à 22% sur
une section de 50 µm. Les deux ne se règlent donc pas séparément, et l'épaisseur
d'oxyde est recalculée pour chaque candidat.

Une limite à connaître : modeler le profil d'épaisseur n'offre AUCUN levier
supplémentaire sur la courbure thermique. Le coefficient global reste
<beta>·alpha_oxyde + (1-<beta>)·alpha_silicium, une famille à un seul paramètre
-- exactement ce que règle déjà l'épaisseur d'oxyde. On ne peut pas annuler à la
fois la pente et la courbure par la géométrie seule.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq

from . import constantes as cst


@dataclass
class ReglageThermique:
    """Plage de température visée et coefficients matériau.

    Les coefficients sont exposés ici plutôt que figés dans constantes.py parce
    qu'ils sont des valeurs de littérature : ils doivent pouvoir être remplacés
    par des mesures, et leur incertitude doit pouvoir être propagée
    (voir `etude_sensibilite`)."""
    temperature_min: float = 8.0
    temperature_max: float = 38.0
    temperature_reference: float = 23.0
    cte_balancier: float = 16.7e-6          # 1/K, balancier CuBe2 (Glucydur)
    tolerance_s_par_j_par_degc: float = 0.1
    actif: bool = True
    # Fenêtre d'épaisseur d'oxyde que le procédé sait tenir. 0/0 = pas de
    # contrainte : l'optimiseur est alors libre, et l'épaisseur requise se lit
    # simplement sur le front de Pareto. La contrainte n'a de sens que si une
    # même oxydation doit servir à toute une série.
    oxyde_min_um: float = 0.0
    oxyde_max_um: float = 0.0

    def contrainte_oxyde_active(self) -> bool:
        return self.oxyde_max_um > self.oxyde_min_um >= 0.0 and self.oxyde_max_um > 0.0

    e_silicium: float = cst.E_SILICIUM
    tce_silicium: float = cst.TCE_SILICIUM
    tce2_silicium: float = cst.TCE2_SILICIUM
    cte_silicium: float = cst.CTE_SILICIUM
    e_oxyde: float = cst.E_OXYDE
    tce_oxyde: float = cst.TCE_OXYDE
    tce2_oxyde: float = cst.TCE2_OXYDE
    cte_oxyde: float = cst.CTE_OXYDE


def rigidite_composite(e_m, h_m, t_m, temperature: float, reglage: ReglageThermique):
    """Rigidité de flexion E·I (N.m²) dans le plan du spiral, d'une lame de
    silicium de section e_m x h_m oxydée sur ses quatre faces d'une épaisseur
    d'oxyde t_m, à la température donnée.

    La section oxydée est le rectangle extérieur MOINS le coeur de silicium,
    donc E·I = E_oxyde·(I_ext - I_coeur) + E_silicium·I_coeur -- exact pour une
    couche d'épaisseur uniforme sur un rectangle.

    e_m est l'épaisseur du silicium AVANT oxydation (la cote de gravure) :
    l'oxydation consomme ensuite 0.44·t de chaque côté et déborde de 0.56·t.
    """
    dT = float(temperature) - reglage.temperature_reference
    dil_si = 1.0 + reglage.cte_silicium * dT
    dil_ox = 1.0 + reglage.cte_oxyde * dT
    conso, croiss = cst.CONSOMMATION_OXYDATION, cst.CROISSANCE_OXYDATION

    e_coeur = np.maximum(e_m - 2.0 * conso * t_m, 1e-12) * dil_si
    h_coeur = np.maximum(h_m - 2.0 * conso * t_m, 1e-12) * dil_si
    e_ext = (e_m + 2.0 * croiss * t_m) * dil_ox
    h_ext = (h_m + 2.0 * croiss * t_m) * dil_ox

    i_coeur = h_coeur * e_coeur ** 3 / 12.0
    i_ext = h_ext * e_ext ** 3 / 12.0

    module_si = reglage.e_silicium * (1.0 + reglage.tce_silicium * dT
                                      + 0.5 * reglage.tce2_silicium * dT ** 2)
    module_ox = reglage.e_oxyde * (1.0 + reglage.tce_oxyde * dT
                                   + 0.5 * reglage.tce2_oxyde * dT ** 2)
    return module_ox * (i_ext - i_coeur) + module_si * i_coeur


def part_oxyde(e_m, h_m, t_m, reglage: ReglageThermique):
    """Fraction de la rigidité de flexion portée par l'oxyde, à la température
    de référence. C'est elle qui fixe le coefficient thermoélastique local :
    alpha = beta·alpha_oxyde + (1-beta)·alpha_silicium."""
    conso, croiss = cst.CONSOMMATION_OXYDATION, cst.CROISSANCE_OXYDATION
    e_c = np.maximum(e_m - 2.0 * conso * t_m, 1e-12)
    h_c = np.maximum(h_m - 2.0 * conso * t_m, 1e-12)
    i_c = h_c * e_c ** 3 / 12.0
    i_e = (h_m + 2.0 * croiss * t_m) * (e_m + 2.0 * croiss * t_m) ** 3 / 12.0
    ox = reglage.e_oxyde * (i_e - i_c)
    return ox / (ox + reglage.e_silicium * i_c)


def _frequence(s_m, e_m, h_m, t_m, temperature, inertie_balancier, reglage):
    """Fréquence propre (Hz) à une température donnée, oxyde compris."""
    dT = float(temperature) - reglage.temperature_reference
    rigidite = rigidite_composite(e_m, h_m, t_m, temperature, reglage)
    longueur = s_m * (1.0 + reglage.cte_silicium * dT)
    couple = 1.0 / trapezoid(1.0 / rigidite, x=longueur)
    inertie = inertie_balancier * (1.0 + 2.0 * reglage.cte_balancier * dT) ** 1.0
    return np.sqrt(couple / inertie) / (2.0 * np.pi)


def marche_temperature(geom, t_oxyde_um: float, temperatures, inertie_balancier: float,
                        reglage: ReglageThermique = None) -> np.ndarray:
    """Marche (s/j) en fonction de la température, rapportée à la température
    de référence. `geom.epaisseur` est la cote de GRAVURE (silicium nu)."""
    reglage = reglage or ReglageThermique()
    actif = geom.partie_active()
    s_m, e_m, h_m = actif.s * 1e-3, actif.epaisseur * 1e-3, actif.hauteur * 1e-3
    t_m = t_oxyde_um * 1e-6
    f_ref = _frequence(s_m, e_m, h_m, t_m, reglage.temperature_reference,
                       inertie_balancier, reglage)
    temperatures = np.atleast_1d(np.asarray(temperatures, dtype=float))
    return np.array([86400.0 * (_frequence(s_m, e_m, h_m, t_m, T, inertie_balancier,
                                            reglage) / f_ref - 1.0)
                     for T in temperatures])


def pente_thermique(geom, t_oxyde_um: float, inertie_balancier: float,
                     reglage: ReglageThermique = None) -> float:
    """Pente moyenne de marche (s/j/°C) entre les deux bornes de température."""
    reglage = reglage or ReglageThermique()
    m = marche_temperature(geom, t_oxyde_um,
                           [reglage.temperature_min, reglage.temperature_max],
                           inertie_balancier, reglage)
    return float((m[1] - m[0]) / (reglage.temperature_max - reglage.temperature_min))


def epaisseur_oxyde_optimale(geom, inertie_balancier: float,
                              reglage: ReglageThermique = None,
                              t_max_um: float = None) -> float:
    """Épaisseur d'oxyde (µm) annulant la pente de marche sur la plage visée.

    Bornée par l'épaisseur qui consommerait entièrement la section la plus
    mince : au-delà il n'y aurait plus de silicium, ce qui n'a pas de sens.
    Retourne NaN si aucune épaisseur de la plage n'annule la pente (cas d'un
    balancier dont la dilatation domine, par exemple)."""
    reglage = reglage or ReglageThermique()
    e_min_um = float(np.min(geom.partie_active().epaisseur)) * 1000.0
    plafond = 0.95 * e_min_um / (2.0 * cst.CONSOMMATION_OXYDATION)
    if t_max_um is not None:
        plafond = min(plafond, t_max_um)
    if plafond <= 1e-3:
        return float("nan")
    # Pour la RECHERCHE de racine, la fréquence de référence se simplifie :
    # la pente est du signe de f(T_max) - f(T_min). On économise ainsi un
    # tiers des évaluations, ce qui compte puisque cette résolution tourne
    # pour chaque candidat de l'optimisation.
    actif = geom.partie_active()
    s_m, e_m, h_m = actif.s * 1e-3, actif.epaisseur * 1e-3, actif.hauteur * 1e-3

    def f(t_um):
        t_m = t_um * 1e-6
        return (_frequence(s_m, e_m, h_m, t_m, reglage.temperature_max,
                           inertie_balancier, reglage)
                - _frequence(s_m, e_m, h_m, t_m, reglage.temperature_min,
                             inertie_balancier, reglage))

    bas, haut = f(0.0), f(plafond)
    if bas * haut > 0.0:
        return float("nan")
    return float(brentq(f, 0.0, plafond, xtol=1e-6, rtol=8.9e-16))


def analyser_thermique(geom, inertie_balancier: float, reglage: ReglageThermique = None,
                        t_oxyde_um: float = None, n_temperatures: int = 13) -> dict:
    """Analyse thermique complète d'une géométrie.

    t_oxyde_um : épaisseur imposée. Si None, l'épaisseur optimale est calculée
    pour CE profil d'épaisseur -- l'optimiseur faisant varier le profil, elle
    change d'un candidat à l'autre."""
    reglage = reglage or ReglageThermique()
    if t_oxyde_um is None:
        t_oxyde_um = epaisseur_oxyde_optimale(geom, inertie_balancier, reglage)
    if not np.isfinite(t_oxyde_um):
        return {"valide": False, "raison": "aucune épaisseur d'oxyde n'annule la pente",
                "epaisseur_oxyde_um": float("nan"), "pente_s_par_j_par_degc": float("nan"),
                "amplitude_s_par_j": float("nan")}

    temperatures = np.linspace(reglage.temperature_min, reglage.temperature_max, n_temperatures)
    marches = marche_temperature(geom, t_oxyde_um, temperatures, inertie_balancier, reglage)
    pente = float((marches[-1] - marches[0])
                  / (reglage.temperature_max - reglage.temperature_min))

    actif = geom.partie_active()
    beta = part_oxyde(actif.epaisseur * 1e-3, actif.hauteur * 1e-3, t_oxyde_um * 1e-6, reglage)
    return {
        "valide": True,
        "epaisseur_oxyde_um": float(t_oxyde_um),
        "croissance_lame_um": float(2.0 * cst.CROISSANCE_OXYDATION * t_oxyde_um),
        "silicium_consomme_um": float(2.0 * cst.CONSOMMATION_OXYDATION * t_oxyde_um),
        "temperatures": temperatures,
        "marches": marches,
        "pente_s_par_j_par_degc": pente,
        "amplitude_s_par_j": float(marches.max() - marches.min()),
        "part_oxyde_min": float(np.min(beta)),
        "part_oxyde_max": float(np.max(beta)),
    }


def etude_sensibilite(geom, inertie_balancier: float, reglage: ReglageThermique = None,
                       incertitudes: dict = None) -> list:
    """De combien l'épaisseur d'oxyde requise et la marche résiduelle bougent
    quand chaque coefficient varie dans son incertitude.

    Les coefficients matériau sont des valeurs de littérature, pas des mesures :
    sans cette étude, le chiffre de résidu affiché donnerait une confiance
    qu'il ne mérite pas. Chaque coefficient est perturbé seul, l'épaisseur
    d'oxyde est RECALCULÉE (on suppose le procédé recalé), et on regarde ce
    qu'il reste -- puis on regarde aussi ce que coûte de NE PAS recaler."""
    reglage = reglage or ReglageThermique()
    if incertitudes is None:
        incertitudes = {
            "tce_silicium": 5e-6,        # -60 +/- 5 ppm/K
            "tce_oxyde": 20e-6,          # +215 +/- 20 ppm/K
            "e_oxyde": 5e9,              # 70 +/- 5 GPa
            "cte_balancier": 1.0e-6,     # selon l'alliage du balancier
            "cte_silicium": 0.1e-6,
        }
    from dataclasses import replace as _replace

    t_ref = epaisseur_oxyde_optimale(geom, inertie_balancier, reglage)
    lignes = []
    for nom, delta in incertitudes.items():
        ligne = {"coefficient": nom, "incertitude": delta}
        for signe in (-1.0, +1.0):
            perturbe = _replace(reglage, **{nom: getattr(reglage, nom) + signe * delta})
            t_perturbe = epaisseur_oxyde_optimale(geom, inertie_balancier, perturbe)
            # pente qui subsisterait si l'on gardait l'épaisseur de référence
            pente_non_recalee = pente_thermique(geom, t_ref, inertie_balancier, perturbe)
            cle = "moins" if signe < 0 else "plus"
            ligne[f"t_{cle}_um"] = float(t_perturbe)
            ligne[f"pente_non_recalee_{cle}"] = float(pente_non_recalee)
        lignes.append(ligne)
    return lignes
