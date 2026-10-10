# sauvegarde.py
"""
Nommage et écriture des fichiers produits par une session de travail.

Le problème que ce module résout n'est pas d'écrire des fichiers — presets.py
et dxf_io.py savent déjà le faire — mais de les NOMMER de façon à les retrouver
des semaines plus tard, sans les ouvrir. Trois principes :

  1. Un DOSSIER PAR SESSION, préfixé du nom de projet puis horodaté :
     `K2_20261010-1432`. Le nom de projet groupe tout le travail d'un calibre,
     l'horodatage trie chronologiquement à l'intérieur. L'ordre inverse
     (date d'abord) trierait toutes les études mêlées.

  2. Un PRÉFIXE COMMUN à tous les fichiers d'une session, identique au nom du
     dossier. Un fichier retrouvé seul, détaché de son dossier, reste
     identifiable — c'est le cas fréquent d'un DXF envoyé par courriel.

  3. Les SOLUTIONS d'optimisation portent un nom court et prévisible
     (`sol-03_fini.dxf`), et un `index.csv` les décrit toutes. Mettre les
     métriques dans le nom donnerait des noms longs et fragiles ; l'index se
     trie dans un tableur et se lit en une ligne de commande.

Rien ici ne dépend de Qt : l'interface appelle ces fonctions, elle ne les
réimplémente pas.
"""
from dataclasses import asdict
from datetime import datetime
import csv
import json
import os
import re

import numpy as np

from . import dxf_io
from .presets import sauvegarder_configuration

# Caractères refusés par au moins un système de fichiers courant. On les retire
# plutôt que de les remplacer : un nom de projet n'a aucune raison d'en porter.
_INTERDITS = re.compile(r'[\\/:*?"<>|]')


def nettoyer_nom(nom: str) -> str:
    """Nom de projet utilisable comme nom de dossier. Les espaces deviennent
    des tirets bas : un nom avec espaces complique toute manipulation en ligne
    de commande, et ces fichiers sont faits pour être transmis."""
    nom = _INTERDITS.sub("", str(nom)).strip()
    nom = re.sub(r"\s+", "_", nom)
    return nom


def horodatage(instant: datetime = None) -> str:
    """AAAAMMJJ-HHMM : trie chronologiquement en ordre alphabétique, ce qui
    n'est vrai d'aucun format local."""
    return (instant or datetime.now()).strftime("%Y%m%d-%H%M")


def dossier_session(racine: str, projet: str, instant: datetime = None,
                     creer: bool = True) -> str:
    """`racine/<projet>_<horodatage>`, suffixé `_2`, `_3`... si ce dossier
    existe déjà.

    La résolution à la minute suffit en pratique et garde les noms courts ;
    le suffixe couvre le cas de deux enregistrements rapprochés. Réutiliser le
    dossier existant serait plus simple mais écraserait un travail précédent
    sans le dire.
    """
    projet = nettoyer_nom(projet)
    if not projet:
        raise ValueError("Nom de projet vide ou composé uniquement de "
                         "caractères interdits.")
    base = f"{projet}_{horodatage(instant)}"
    chemin = os.path.join(racine, base)
    n = 1
    while os.path.isdir(chemin) and os.listdir(chemin):
        n += 1
        chemin = os.path.join(racine, f"{base}_{n}")
    if creer:
        os.makedirs(chemin, exist_ok=True)
    return chemin


def prefixe(dossier: str) -> str:
    """Préfixe commun des fichiers d'une session : le nom du dossier lui-même."""
    return os.path.basename(os.path.normpath(dossier))


# ------------------------------------------------------------------ JSON
def _convertir(valeur):
    """numpy -> types JSON. Les tableaux deviennent des listes, les scalaires
    des nombres Python : `json` ne sait sérialiser ni les uns ni les autres."""
    if isinstance(valeur, np.ndarray):
        return [_convertir(v) for v in valeur.tolist()]
    if isinstance(valeur, (np.floating, np.integer)):
        return valeur.item()
    if isinstance(valeur, dict):
        return {k: _convertir(v) for k, v in valeur.items()}
    if isinstance(valeur, (list, tuple)):
        return [_convertir(v) for v in valeur]
    if isinstance(valeur, float) and not np.isfinite(valeur):
        # NaN et inf ne sont pas du JSON valide : un lecteur strict les refuse.
        return None
    return valeur


def ecrire_json(chemin: str, donnees: dict):
    with open(chemin, "w", encoding="utf-8") as f:
        json.dump(_convertir(donnees), f, indent=2, ensure_ascii=False)
    return chemin


# ------------------------------------------------------- écritures typées
def enregistrer_configuration(dossier: str, spec, parametres, bornes=None,
                               hyperparametres=None, reglage_echappement=None) -> str:
    """`<stem>_config.json` : de quoi reconstruire exactement la forme."""
    from .presets import BornesOptimisation, HyperparametresGA
    chemin = os.path.join(dossier, f"{prefixe(dossier)}_config.json")
    sauvegarder_configuration(chemin, spec, parametres,
                              bornes or BornesOptimisation(),
                              hyperparametres or HyperparametresGA(),
                              reglage_echappement)
    return chemin


def enregistrer_analyse(dossier: str, resultats: dict, nom: str = "analyse") -> str:
    """`<stem>_analyse.json` : les grandeurs calculées, pas la façon de les
    recalculer. Séparé de la configuration pour que l'un se relise sans
    l'autre — on compare souvent des analyses entre elles."""
    chemin = os.path.join(dossier, f"{prefixe(dossier)}_{nom}.json")
    return ecrire_json(chemin, resultats)


def enregistrer_geometrie(dossier: str, geom, nom: str = None,
                           n_points_dxf: int = None) -> list:
    """DXF de la géométrie : deux fichiers si la lame est oxydée (cote de
    gravure et pièce finie, qui diffèrent d'environ 4 µm), un seul sinon."""
    base = os.path.join(dossier, f"{nom or prefixe(dossier)}.dxf")
    return dxf_io.exporter_dxf_paire(geom, base, n_points=n_points_dxf)


# ------------------------------------------------- solutions d'optimisation
# Colonnes de l'index, dans l'ordre. Le rang d'abord : c'est la clé qui relie
# une ligne à ses fichiers.
COLONNES_INDEX = ("rang", "dispersion_s_j", "dvh_s_j", "jeu_um",
                  "compensation_s_j", "oxyde_um", "frequence_hz", "fichiers")


def nom_solution(rang: int) -> str:
    """`sol-03` : deux chiffres, donc un tri alphabétique correct jusqu'à 99
    solutions, ce qui dépasse largement la taille d'un front de Pareto."""
    return f"sol-{int(rang):02d}"


def enregistrer_solutions(dossier: str, solutions: list,
                           n_points_dxf: int = None) -> dict:
    """Écrit les solutions retenues dans `<dossier>/solutions/` et retourne les
    chemins écrits.

    `solutions` : liste de dicts portant au moins `rang` et `geometrie`, plus
    les métriques à indexer. Seules les solutions fournies sont écrites — c'est
    à l'appelant de n'envoyer que celles que l'utilisateur a cochées.

    L'index est un CSV et non un JSON : il est fait pour être ouvert dans un
    tableur, trié par dispersion ou par DVH, et comparé d'une session à
    l'autre. Le séparateur est le point-virgule, qu'Excel francophone lit sans
    boîte de dialogue d'import.
    """
    if not solutions:
        raise ValueError("Aucune solution sélectionnée à exporter.")
    sous = os.path.join(dossier, "solutions")
    os.makedirs(sous, exist_ok=True)

    lignes, ecrits = [], {}
    for sol in solutions:
        nom = nom_solution(sol["rang"])
        chemins = dxf_io.exporter_dxf_paire(sol["geometrie"],
                                            os.path.join(sous, f"{nom}.dxf"),
                                            n_points=n_points_dxf)
        detail = {k: v for k, v in sol.items() if k != "geometrie"}
        detail["fichiers"] = [os.path.basename(c) for c in chemins]
        ecrire_json(os.path.join(sous, f"{nom}.json"), detail)
        ecrits[nom] = chemins + [os.path.join(sous, f"{nom}.json")]
        lignes.append({
            "rang": sol["rang"],
            "dispersion_s_j": sol.get("dispersion"),
            "dvh_s_j": sol.get("dvh"),
            "jeu_um": sol.get("jeu_um"),
            "compensation_s_j": sol.get("ecart_compensation"),
            "oxyde_um": sol.get("oxyde_um"),
            "frequence_hz": sol.get("frequence_hz"),
            "fichiers": " ".join(os.path.basename(c) for c in chemins),
        })

    chemin_index = os.path.join(sous, "index.csv")
    with open(chemin_index, "w", encoding="utf-8-sig", newline="") as f:
        # utf-8-sig : sans la marque d'ordre, Excel lit les accents de travers.
        ecrivain = csv.DictWriter(f, fieldnames=COLONNES_INDEX, delimiter=";")
        ecrivain.writeheader()
        for ligne in sorted(lignes, key=lambda l: l["rang"]):
            ecrivain.writerow({k: _format_csv(ligne[k]) for k in COLONNES_INDEX})
    ecrits["index"] = [chemin_index]
    return ecrits


def _format_csv(valeur):
    """Décimales en virgule : c'est ce qu'attend un tableur francophone, et
    c'est cohérent avec le point-virgule comme séparateur de colonnes."""
    if valeur is None:
        return ""
    if isinstance(valeur, (float, np.floating)):
        if not np.isfinite(valeur):
            return ""
        return f"{float(valeur):.4f}".replace(".", ",")
    return str(valeur)
