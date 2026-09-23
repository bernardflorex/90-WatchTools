# rendu.py
"""
Utilitaires de rendu graphique (animation GIF du déploiement du spiral).
Séparé de physique.py, qui doit rester utilisable sans dépendance graphique
(appelé en parallèle/headless par l'optimiseur)."""
import numpy as np
import matplotlib.pyplot as plt
# Pas de matplotlib.use() ici : ce module est importé dans le même process que
# l'appli Qt (backend QtAgg déjà actif) -- changer le backend globalement ici
# casserait les FigureCanvasQTAgg déjà créés dans app.py. L'encodage GIF via
# PillowWriter fonctionne indépendamment du backend interactif choisi.
from matplotlib.animation import FuncAnimation, PillowWriter

from .physique import resoudre_elastica
from .modele import GeometrieSpiral


def calculer_frames_deploiement(geom: GeometrieSpiral, r_min: float, amplitude_max_deg: float,
                                 n_frames: int = 41):
    """Résout l'elastica à n_frames angles régulièrement répartis entre
    -amplitude_max_deg et +amplitude_max_deg. Retourne (angles_deg, frames)
    où frames est une liste de (Xe, Ye, Xi, Yi). Séparé de la génération du
    GIF elle-même pour pouvoir l'exécuter dans un QThread (le calcul
    elastica est la partie lente, l'encodage GIF est rapide une fois les
    frames en main)."""
    angles_deg = np.linspace(-amplitude_max_deg, amplitude_max_deg, n_frames)
    frames = []
    for ang in angles_deg:
        Xe, Ye, Xi, Yi, _, _ = resoudre_elastica(geom, float(ang), r_min)
        frames.append((Xe, Ye, Xi, Yi))
    return angles_deg, frames


def indices_aller_retour(n_frames: int):
    """Indices formant un cycle complet aller-retour (0 -> n-1 -> 0), sans
    dupliquer les deux extrémités. Un balancier réel oscille en va-et-vient :
    une animation qui se contente d'aller de -amplitude à +amplitude puis
    boucle instantanément à son point de départ produit un "saut" visuel à
    chaque boucle, différent de ce qu'on voit dans un post-traitement
    ANSYS typique (animation lisse aller-retour). Utilisé à la fois pour
    l'export GIF et pour la lecture interactive (app.py), pour un
    comportement identique aux deux endroits."""
    if n_frames <= 1:
        return list(range(n_frames))
    return list(range(n_frames)) + list(range(n_frames - 2, 0, -1))


def sauvegarder_gif_deploiement(angles_deg, frames, chemin_gif: str, fps: int = 15,
                                 aller_retour: bool = True):
    """Encode les frames déjà calculées (calculer_frames_deploiement) en GIF
    animé, en va-et-vient par défaut (voir indices_aller_retour). Appelé sur
    le thread principal (matplotlib/Qt) après réception des frames calculées
    en arrière-plan."""
    toutes_x = np.concatenate([np.concatenate([f[0], f[2]]) for f in frames])
    toutes_y = np.concatenate([np.concatenate([f[1], f[3]]) for f in frames])
    marge = 0.1 * max(toutes_x.max() - toutes_x.min(), toutes_y.max() - toutes_y.min(), 1e-6)

    fig, ax = plt.subplots(figsize=(6, 6))
    ligne_ext, = ax.plot([], [], "k-", lw=0.9)
    ligne_int, = ax.plot([], [], "k-", lw=0.9)
    ax.plot(0, 0, "k+", ms=10)
    ax.set_xlim(toutes_x.min() - marge, toutes_x.max() + marge)
    ax.set_ylim(toutes_y.min() - marge, toutes_y.max() + marge)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    titre = ax.set_title("")

    sequence = indices_aller_retour(len(frames)) if aller_retour else list(range(len(frames)))

    def update(k):
        i = sequence[k]
        Xe, Ye, Xi, Yi = frames[i]
        ligne_ext.set_data(Xe, Ye)
        ligne_int.set_data(Xi, Yi)
        titre.set_text(f"angle balancier = {angles_deg[i]:+.0f}°")
        return ligne_ext, ligne_int, titre

    anim = FuncAnimation(fig, update, frames=len(sequence), blit=False)
    anim.save(chemin_gif, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return chemin_gif
