# threads.py
"""
QThread pour les opérations longues (sondage, optimisation NSGA-II), pour ne
pas geler l'interface -- même principe que spiropt_v1/optimization.py,
adapté aux fonctions pures de optimisation.py (qui ne connaissent rien de Qt).
"""
import numpy as np
from PySide6.QtCore import QThread, Signal
from pymoo.core.callback import Callback

from .optimisation import sonder_faisabilite, lancer_optimisation
from .rendu import calculer_frames_deploiement


class SondageThread(QThread):
    signal_termine = Signal(object)

    def __init__(self, base, bornes, amplitudes_analyse_deg,
                 n_echantillons, seed, pas_balayage_deg,
                 reglage_echappement=None, reglage_thermique=None, parent=None):
        super().__init__(parent)
        self.base = base
        self.bornes = bornes
        self.amplitudes_analyse_deg = amplitudes_analyse_deg
        self.n_echantillons = n_echantillons
        self.seed = seed
        self.pas_balayage_deg = pas_balayage_deg
        self.reglage_echappement = reglage_echappement
        self.reglage_thermique = reglage_thermique

    def run(self):
        resultats = sonder_faisabilite(
            self.base, self.bornes, self.amplitudes_analyse_deg,
            n_echantillons=self.n_echantillons, seed=self.seed, pas_balayage_deg=self.pas_balayage_deg,
            reglage_echappement=self.reglage_echappement,
            reglage_thermique=self.reglage_thermique,
        )
        self.signal_termine.emit(resultats)


class _CallbackProgres(Callback):
    def __init__(self, signal_progres):
        super().__init__()
        self.signal_progres = signal_progres

    def notify(self, algorithm):
        gen = algorithm.n_gen
        pop_F = algorithm.pop.get("F")
        if pop_F is not None and len(pop_F) > 0:
            meilleure_dispersion = float(np.min(pop_F[:, 0]))
            meilleur_jeu = float(-np.min(pop_F[:, 1]))
        else:
            meilleure_dispersion, meilleur_jeu = float("nan"), float("nan")
        self.signal_progres.emit(gen, meilleure_dispersion, meilleur_jeu)


class OptimisationThread(QThread):
    signal_progres = Signal(int, float, float)
    signal_termine = Signal(object)

    def __init__(self, base, bornes, hyperparams, jeu_min_requis,
                 amplitudes_analyse_deg, x0=None, pas_balayage_deg=20.0,
                 reglage_echappement=None, reglage_thermique=None, parent=None):
        super().__init__(parent)
        self.base = base
        self.bornes = bornes
        self.hyperparams = hyperparams
        self.jeu_min_requis = jeu_min_requis
        self.amplitudes_analyse_deg = amplitudes_analyse_deg
        self.x0 = x0
        self.pas_balayage_deg = pas_balayage_deg
        self.reglage_echappement = reglage_echappement
        self.reglage_thermique = reglage_thermique

    def run(self):
        callback = _CallbackProgres(self.signal_progres)
        res = lancer_optimisation(
            self.base, self.bornes, self.hyperparams, self.jeu_min_requis,
            amplitudes_analyse_deg=self.amplitudes_analyse_deg,
            warm_start=True, x0=self.x0, callback=callback,
            pas_balayage_deg=self.pas_balayage_deg,
            reglage_echappement=self.reglage_echappement,
            reglage_thermique=self.reglage_thermique,
        )
        self.signal_termine.emit(res)


class AnimationThread(QThread):
    """Calcule les frames (résolutions elastica, la partie lente) en
    arrière-plan ; l'encodage GIF lui-même (rapide) reste sur le thread
    principal (rendu.sauvegarder_gif_deploiement), matplotlib/Qt n'étant pas
    à l'aise pour dessiner depuis un thread secondaire."""
    signal_termine = Signal(object, object)

    def __init__(self, geom, r_min, amplitude_max_deg, n_frames=41, parent=None):
        super().__init__(parent)
        self.geom = geom
        self.r_min = r_min
        self.amplitude_max_deg = amplitude_max_deg
        self.n_frames = n_frames

    def run(self):
        angles_deg, frames = calculer_frames_deploiement(
            self.geom, self.r_min, self.amplitude_max_deg, n_frames=self.n_frames
        )
        self.signal_termine.emit(angles_deg, frames)
