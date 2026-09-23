# convergence.py
"""
Étude de convergence du maillage physique (nombre d'éléments discrets).

Contexte : une note historique du projet (spiropt_v1/config.py) indiquait que
n_elements=600 biaisait l'orbite du CG d'environ +40% par rapport à la valeur
convergée (atteinte vers 2000-3000 éléments), sur la base d'une comparaison
ANSYS ponctuelle. Cette étude généralise et automatise la vérification pour
l'architecture spiral_lab actuelle (theta-interpolation DXF, Bézier
Grossmann/Phillips), sur les deux voies de génération de géométrie (import
DXF réel, génération paramétrique), avec les trois métriques qui comptent
réellement pour l'utilisateur : fréquence propre, jeu inter-spires minimal,
dispersion de marche.

Ne dépend pas de Qt : utilisable en script ou depuis l'appli.
"""
import time
import numpy as np

from .physique import frequence_propre
from .isochronisme import marche_positions_verticales, critere_dispersion_marche


def etude_convergence(constructeur_geom, valeurs_n, r_min: float, freq_cible: float,
                       inertie_balancier: float, amplitudes_deg: np.ndarray,
                       pas_balayage_deg: float = 15.0, verbose: bool = True):
    """constructeur_geom(n) -> GeometrieSpiral pour chaque n de valeurs_n
    (typiquement une fermeture sur dxf_io.charger_dxf ou
    parametrage.construire_geometrie avec n_elements variable). Retourne une
    liste de dicts {n, temps_s, freq_hz, jeu_min_um, dispersion_s_j}, dans
    l'ordre de valeurs_n.

    Une seule résolution elastica par angle est faite par valeur de n (via
    marche_positions_verticales(return_diagnostics=True), qui donne à la fois
    la marche ET le jeu balayé en un seul passage) -- pas de calcul redondant."""
    resultats = []
    for n in valeurs_n:
        t0 = time.time()
        geom = constructeur_geom(n)
        freq = frequence_propre(geom, inertie_balancier)
        marches, diag = marche_positions_verticales(
            geom, r_min, freq_cible, inertie_balancier,
            amplitudes_deg=amplitudes_deg, pas_balayage_deg=pas_balayage_deg,
            return_diagnostics=True,
        )
        dispersion = critere_dispersion_marche(marches)
        jeu_min_um = float(diag["jeu_mm"].min()) * 1000.0
        dt = time.time() - t0
        resultats.append({
            "n": geom.n, "temps_s": dt, "freq_hz": float(freq),
            "jeu_min_um": jeu_min_um, "dispersion_s_j": float(dispersion),
        })
        if verbose:
            print(f"  n={geom.n:5d}  freq={freq:.5f} Hz  jeu_min={jeu_min_um:8.2f} um  "
                  f"dispersion={dispersion:7.3f} s/j  ({dt:.1f}s)")
    return resultats


def rapport_convergence(resultats: list, seuil_pct: float = 1.0) -> str:
    """Compare chaque ligne à la plus fine (dernière) et signale, pour
    chaque métrique, le plus petit n à partir duquel l'écart relatif à la
    référence la plus fine reste sous seuil_pct (%) pour TOUTES les valeurs
    de n suivantes (convergence "stable", pas juste un croisement ponctuel)."""
    ref = resultats[-1]
    lignes = [f"{'n':>6} | {'freq (Hz)':>10} | {'jeu_min (um)':>13} | {'dispersion (s/j)':>17} | {'temps (s)':>9}"]
    lignes.append("-" * len(lignes[0]))
    for r in resultats:
        lignes.append(f"{r['n']:6d} | {r['freq_hz']:10.5f} | {r['jeu_min_um']:13.2f} | "
                       f"{r['dispersion_s_j']:17.3f} | {r['temps_s']:9.1f}")
    lignes.append("")

    for cle, nom in [("freq_hz", "fréquence propre"), ("jeu_min_um", "jeu minimal"),
                      ("dispersion_s_j", "dispersion de marche")]:
        val_ref = ref[cle]
        n_convergence = None
        for i in range(len(resultats)):
            ecarts_suivants = [
                abs(r[cle] - val_ref) / max(abs(val_ref), 1e-12) * 100.0
                for r in resultats[i:]
            ]
            if all(e <= seuil_pct for e in ecarts_suivants):
                n_convergence = resultats[i]["n"]
                break
        lignes.append(
            f"{nom:22s} : convergé (<{seuil_pct}% de la référence n={ref['n']}) "
            f"à partir de n={n_convergence}" if n_convergence else
            f"{nom:22s} : PAS convergé à {seuil_pct}% sur la plage testée"
        )
    return "\n".join(lignes)
