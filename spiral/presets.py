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

from .echappement import ReglageEchappement


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
    # Côté depuis lequel le spiral est vu et depuis lequel les positions
    # verticales sont nommées. Le spiral étant sous le coq, la vue naturelle
    # est celle du rouage ; l'autre choix échange 3H et 9H.
    vue_cote_rouage: bool = True
    # Repère du MOUVEMENT, vu côté rouage et dessiné 6H EN HAUT -- l'attitude
    # dans laquelle la pièce est examinée, et celle de la position COSC 6H. La
    # tige de remontoir, qui sort à 3H, est alors À DROITE du dessin.
    # La LIGNE DES CENTRES -- celle qui joint l'axe du balancier à
    # celui de l'ancre, dite aussi ligne droite de l'échappement -- fait cet
    # angle avec l'axe de tige, compté dans le sens trigonométrique TEL QU'ON
    # LE VOIT côté fond. C'est sur elle que se fixe le spiral au collet, et, le
    # nombre de tours étant entier, c'est sur elle que retombe le point de
    # comptage, à l'opposé.
    angle_echappement_deg: float = -6.2
    # Sens de développement du spiral depuis le collet, vu côté fond.
    sens_trigonometrique: bool = True
    # LEQUEL DES DEUX BOUTS de la ligne des centres porte le collet. Une droite
    # n'est définie qu'à 180° près par son angle : sans cette information, le
    # spiral peut se retrouver dessiné en demi-tour par rapport au mouvement
    # réel -- ce qui n'altère pas la dispersion mais inverse le signe du DVH.
    # Vrai (défaut) : le bout qui pointe vers l'ancre, à l'opposé de la tige.
    collet_cote_ancre: bool = True
    # Plan cristallin du wafer. Le silicium n'a un module de Young isotrope
    # DANS LE PLAN que sur un wafer (111) : sur un (001) il varie de 130 à
    # 169 GPa selon la direction, avec une période de 90°, et la lame d'un
    # spiral balaie forcément toutes ces directions. Le défaut "111" reproduit
    # le comportement antérieur (169 GPa partout).
    plan_cristallin: str = "111"
    # Angle, dans le repère du dessin, de l'axe cristallin de référence du plan
    # (<100> pour un (001)). Sans effet sur un (111).
    orientation_cristal_deg: float = 0.0

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


_CLES_ARRAYS_PARAMETRES = ("epaisseur_noeuds", "delta_r_noeuds", "grossmann_offsets",
                            "phillips_offsets", "raquette_offsets")


def sauvegarder_configuration(chemin: str, spec: SpecMouvement, parametres,
                               bornes: "BornesOptimisation", hyperparams: "HyperparametresGA",
                               reglage_echappement=None):
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
        if cle in p_dict:
            p_dict[cle] = list(np.asarray(p_dict[cle], dtype=float))
    if reglage_echappement is None:
        reglage_echappement = ReglageEchappement()
    data = {
        "spec": asdict(spec), "parametres": p_dict,
        "bornes": asdict(bornes), "hyperparametres_ga": asdict(hyperparams),
        "echappement": asdict(reglage_echappement),
    }
    with open(chemin, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def charger_configuration(chemin: str):
    """Inverse de sauvegarder_configuration. Retourne (spec, parametres,
    bornes, hyperparams, reglage_echappement). Import local de
    ParametresSpiral pour la même raison (éviter tout risque de cycle
    presets<->parametrage). Un fichier écrit avant l'introduction du réglage
    d'échappement se recharge sans erreur, avec les valeurs par défaut."""
    from .parametrage import ParametresSpiral

    with open(chemin, "r", encoding="utf-8") as f:
        data = json.load(f)

    spec = SpecMouvement(**data["spec"])
    p_dict = dict(data["parametres"])
    for cle in _CLES_ARRAYS_PARAMETRES:
        # `if cle in` : un fichier écrit avant l'introduction de la terminaison
        # à raquette n'a pas la clé, la dataclass reprend alors son défaut
        # (spiral libre) au lieu de lever une KeyError.
        if cle in p_dict:
            p_dict[cle] = np.array(p_dict[cle], dtype=float)
    # Les clés que la dataclass ne connaît plus sont écartées : un fichier
    # écrit avant un remaniement de la paramétrisation reste rechargeable, au
    # lieu de lever un TypeError sur un champ supprimé. Les clés absentes
    # prennent simplement leur valeur par défaut.
    from dataclasses import fields as _fields
    connus = {f.name for f in _fields(ParametresSpiral)}
    ignores = sorted(set(p_dict) - connus)
    parametres = ParametresSpiral(**{k: v for k, v in p_dict.items() if k in connus})
    if ignores:
        print(f"charger_configuration : champs ignorés (paramétrisation modifiée) : "
              f"{', '.join(ignores)}")
    bornes = BornesOptimisation(**data["bornes"])
    hyperparams = HyperparametresGA(**data["hyperparametres_ga"])
    reglage_echappement = ReglageEchappement(**data.get("echappement", {}))
    return spec, parametres, bornes, hyperparams, reglage_echappement
