# presets.py
"""
Remplace spiropt_v1/config.py, qui mélangeait dans un seul fichier plat de
constantes Python quatre choses de nature différente :
  1. des constantes matériau (fixes)       -> restent dans constantes.py
  2. la spec d'UN mouvement/calibre donné  -> SpecMouvement, ci-dessous
  3. les facteurs de construction des bornes de recherche NSGA-II
                                            -> BornesOptimisation, ci-dessous
  4. les hyperparamètres de l'algorithme génétique lui-même
                                            -> HyperparametresGA, ci-dessous

Problème avec l'approche "un config.py de constantes" : changer de calibre,
ou comparer deux jeux de bornes, oblige à éditer du code Python (donc à
committer un fichier source pour changer une donnée). Ici, chaque dataclass
peut être sauvegardée/chargée en JSON (`sauvegarder`/`charger`, même
mécanique que spiral_io.Experience) : plusieurs calibres = plusieurs
fichiers dans presets/, pas plusieurs branches de config.py. Les valeurs
ci-dessous ne sont que des DÉFAUTS de dataclass (repris de
spiropt_v1/config.py pour démarrer sur une base connue), pas des constantes
figées.
"""
from dataclasses import asdict, dataclass, field
import json
import numpy as np


def _sauvegarder(obj, chemin: str):
    with open(chemin, "w", encoding="utf-8") as f:
        json.dump(asdict(obj), f, indent=2, ensure_ascii=False)


def _charger(cls, chemin: str):
    with open(chemin, "r", encoding="utf-8") as f:
        return cls(**json.load(f))


@dataclass
class SpecMouvement:
    """Spec d'UN mouvement/calibre : encombrement, cible chronométrique,
    sécurités. Correspond à ce qui était éparpillé en *_DEFAULT dans
    spiropt_v1/config.py."""
    r_min: float = 0.47
    r_max: float = 2.67
    n_tours_reels: float = 12.0
    n_tours_fictifs: float = 14.0
    hauteur_h: float = 0.100          # mm
    freq_cible: float = 3.5           # Hz
    inertie_balancier: float = 8e-10  # kg.m^2
    jeu_min_spires: float = 0.012     # mm
    amplitude_max: float = 340.0      # deg

    def sauvegarder(self, chemin: str):
        _sauvegarder(self, chemin)

    @staticmethod
    def charger(chemin: str) -> "SpecMouvement":
        return _charger(SpecMouvement, chemin)


@dataclass
class BornesOptimisation:
    """Facteurs de construction des bornes NSGA-II autour d'un point de
    référence -- séparés de SpecMouvement (qui décrit le mouvement CIBLE,
    pas la stratégie d'exploration) et de HyperparametresGA (qui décrit
    l'algorithme, pas l'espace de recherche)."""
    spread_epaisseur_relatif: float = 0.35     # +/- 35% autour de l'épaisseur théorique
    spread_delta_r_mm: float = 0.010           # +/- 10 um sur les offsets radiaux libres
    # Bornes des offsets Grossmann/Phillips en fraction de la largeur RADIALE
    # réelle de chaque zone (calculée par parametrage.rayons_jonction), pas en
    # mm fixe -- une valeur fixe (0.03/0.06 mm) s'est révélée beaucoup trop
    # large pour un calibre à pas serré : ~55% des tirages devenaient
    # non-monotones (diagnostic 2026-09-22, r_min=0.47/r_max=2.67/12 tours,
    # zone Grossmann de ~39 um de large avec un offset autorisé de 30 um).
    # AJUSTÉ (diagnostic 2026-09-24) : à 0.15, l'optimiseur poussait contre
    # cette borne (offset trouvé à ~79% du maximum autorisé côté Grossmann)
    # -- signe que la contrainte était active, pas que l'optimum était
    # atteint. Remonté à 0.25 : compromis mesuré (balayage 0.15/0.20/0.25/0.30
    # sur 400 tirages) entre marge de manoeuvre et taux de tirages non-
    # monotones (13% à 0.15, 34% à 0.25, 47% à 0.30 -- monotonie toujours
    # vérifiée par construire_geometrie, donc pas de risque de repli non
    # détecté, juste un budget de tirages "gâchés" à ne pas trop augmenter).
    spread_grossmann_relatif: float = 0.25
    spread_phillips_relatif: float = 0.25
    angle_grossmann_min_deg: float = 30.0
    angle_grossmann_max_deg: float = 150.0
    angle_phillips_min_deg: float = 60.0
    angle_phillips_max_deg: float = 220.0

    def sauvegarder(self, chemin: str):
        _sauvegarder(self, chemin)

    @staticmethod
    def charger(chemin: str) -> "BornesOptimisation":
        return _charger(BornesOptimisation, chemin)


@dataclass
class HyperparametresGA:
    """Réglages de l'algorithme NSGA-II lui-même (pymoo) -- indépendants du
    mouvement visé et des bornes explorées."""
    pop_size: int = 60
    n_gen: int = 40
    seed: int = 42
    sondage_n_echantillons: int = 40
    fraction_locale_warmstart: float = 0.7
    rayon_relatif_warmstart: float = 0.15

    def sauvegarder(self, chemin: str):
        _sauvegarder(self, chemin)

    @staticmethod
    def charger(chemin: str) -> "HyperparametresGA":
        return _charger(HyperparametresGA, chemin)


_CLES_ARRAYS_PARAMETRES = ("epaisseur_noeuds", "delta_r_noeuds", "grossmann_offsets", "phillips_offsets")


def sauvegarder_configuration(chemin: str, spec: SpecMouvement, parametres,
                               bornes: "BornesOptimisation", hyperparams: "HyperparametresGA"):
    """Sauvegarde EN UN SEUL FICHIER toute la configuration nécessaire pour
    reproduire exactement un spiral et la façon dont il a été réglé/optimisé :
    la spec du mouvement, le vecteur de conception complet (ParametresSpiral,
    y compris la structure des courbes de Grossmann/Phillips), les bornes de
    recherche et les hyperparamètres GA. `parametres` est un
    parametrage.ParametresSpiral (non typé ici pour éviter un import
    circulaire présets<->parametrage -- parametrage.py n'importe pas
    présets.py, mais l'inverse aurait pu créer un cycle si parametrage
    importait un jour presets)."""
    p_dict = asdict(parametres)
    for cle in _CLES_ARRAYS_PARAMETRES:
        p_dict[cle] = list(np.asarray(p_dict[cle], dtype=float))
    data = {
        "spec": asdict(spec), "parametres": p_dict,
        "bornes": asdict(bornes), "hyperparametres_ga": asdict(hyperparams),
    }
    with open(chemin, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def charger_configuration(chemin: str):
    """Inverse de sauvegarder_configuration. Retourne (spec, parametres,
    bornes, hyperparams). Import local de ParametresSpiral pour la même
    raison (éviter tout risque de cycle presets<->parametrage)."""
    from .parametrage import ParametresSpiral

    with open(chemin, "r", encoding="utf-8") as f:
        data = json.load(f)

    spec = SpecMouvement(**data["spec"])
    p_dict = dict(data["parametres"])
    for cle in _CLES_ARRAYS_PARAMETRES:
        p_dict[cle] = np.array(p_dict[cle], dtype=float)
    parametres = ParametresSpiral(**p_dict)
    bornes = BornesOptimisation(**data["bornes"])
    hyperparams = HyperparametresGA(**data["hyperparametres_ga"])
    return spec, parametres, bornes, hyperparams
