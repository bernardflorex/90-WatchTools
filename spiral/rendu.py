# rendu.py
"""
Utilitaires de rendu graphique (animation GIF du déploiement du spiral).
Séparé de physique.py, qui doit rester utilisable sans dépendance graphique
(appelé en parallèle/headless par l'optimiseur)."""
import logging

import numpy as np
import matplotlib.pyplot as plt
# Pas de matplotlib.use() ici : ce module est importé dans le même process que
# l'appli Qt (backend QtAgg déjà actif) -- changer le backend globalement ici
# casserait les FigureCanvasQTAgg déjà créés dans app.py. L'encodage GIF via
# PillowWriter fonctionne indépendamment du backend interactif choisi.
from matplotlib.animation import FuncAnimation, PillowWriter

from commun import theme
from .physique import resoudre_elastica
from .modele import GeometrieSpiral

# L'orbite du centre de gravité est tracée à échelle isotrope avec des limites
# imposées : matplotlib avertit alors, à CHAQUE redessin, qu'il étend l'un des
# deux axes pour respecter l'aspect. C'est le comportement voulu et le message
# n'apprend rien -- il ne ferait que noyer ceux qui comptent. Filtré message par
# message, pour ne pas masquer les autres avertissements du même module. Posé
# ici, dans le module graphique, plutôt que dans l'interface : le rapport de
# simulation trace les mêmes figures sans passer par elle.
logging.getLogger("matplotlib.axes._base").addFilter(
    lambda enregistrement: not enregistrement.getMessage().startswith("Ignoring fixed")
)


def tracer_repere_mouvement(ax, rayon_max: float, cote_rouage: bool = True,
                             angle_collet_deg: float = None, rayon_min: float = None,
                             angle_ancre_deg: float = None):
    """Repère du MOUVEMENT sur un tracé de spiral : repères horaires, axe de
    tige, ligne d'échappement, et mention du côté depuis lequel on regarde.

    Les étiquettes marquent la position des REPÈRES HORAIRES du cadran, pas la
    direction de la gravité. Les confondre décale l'affichage de 180°, car la
    position COSC « 6H » désigne le cas où le 6 est EN HAUT, donc où g pointe
    vers le 12 (erreur corrigée le 2026-09-29).

    Orientation : vue côté rouage, 6H EN HAUT -- l'attitude dans laquelle la
    pièce est examinée, et celle de la position COSC 6H. La tige, qui sort à
    3H, est alors à DROITE. La chiralité reste celle d'une vue de dos :
    12-3-6-9 y tourne dans le sens antihoraire.

    Toutes les décorations restent DANS le rayon du spiral : les limites d'axes
    sont partagées entre les quatre vues d'états et n'ont pas de marge pour
    déborder -- un premier jet les plaçait à 1.16 et 1.35 rayon, et les
    étiquettes sortaient du cadre.
    """
    from .isochronisme import nom_vue, reperes_horaires

    reperes = reperes_horaires(cote_rouage)
    for nom, angle in reperes.items():
        ca, sa = np.cos(angle), np.sin(angle)
        # 1.10 rayon : les repères sortent du spiral, sinon le repère 3H
        # tombe exactement sur le collet (à 6.2° de là, au même rayon) et les
        # deux étiquettes se recouvrent.
        ax.annotate(nom, (1.10 * rayon_max * ca, 1.10 * rayon_max * sa), color=theme.gris(0.30),
                    fontsize=9, ha="center", va="center", zorder=5,
                    bbox=dict(facecolor=theme.fond(), edgecolor="none", alpha=0.8, pad=1.0))

    # Axe de tige : il sort au repère 3H, donc à droite. Tracé d'un seul côté
    # (c'est un demi-axe, la tige ne traverse pas le mouvement) et étiqueté à
    # mi-longueur
    # pour ne pas venir sur le repère horaire qui occupe déjà son extrémité.
    a_tige = reperes["3H"]
    tx, ty = np.cos(a_tige), np.sin(a_tige)
    ax.plot([0.0, 0.88 * rayon_max * tx], [0.0, 0.88 * rayon_max * ty],
            "-", color=theme.gris(0.55), lw=1.2, zorder=4)
    ax.annotate("axe de tige", (0.74 * rayon_max * tx, 0.74 * rayon_max * ty),
                color=theme.gris(0.40), fontsize=7.5, ha="center", va="bottom", zorder=5,
                xytext=(0, 3), textcoords="offset points",
                bbox=dict(facecolor=theme.fond(), edgecolor="none", alpha=0.8, pad=1.0))

    if angle_ancre_deg is not None:
        # Ligne des centres : un DEMI-AXE, tracé depuis l'axe du balancier vers
        # le seul côté où se trouvent l'ancre et la roue d'échappement. Tracée
        # en diamètre entier, elle laissait croire à une symétrie qui n'existe
        # pas : de l'autre côté du balancier il n'y a rien.
        a = np.radians(angle_ancre_deg)
        ex, ey = np.cos(a), np.sin(a)
        ax.plot([0.0, 1.06 * rayon_max * ex], [0.0, 1.06 * rayon_max * ey],
                "--", color="tab:green", lw=1.3, zorder=4)
        # Étiquette ancrée à la POINTE et développée vers l'intérieur : centrée
        # sur la pointe, elle débordait du cadre du côté de l'ancre.
        ax.annotate("ligne des centres\n(vers ancre et roue d'échappement)",
                    (1.06 * rayon_max * ex, 1.06 * rayon_max * ey),
                    color="tab:green", fontsize=7, ha="right", va="top", zorder=5,
                    xytext=(0, -5), textcoords="offset points",
                    bbox=dict(facecolor=theme.fond(), edgecolor="none", alpha=0.85, pad=1.0))

    if angle_collet_deg is not None and rayon_min is not None:
        # Le collet, marqué là où il est réellement : contre la virole. Le point
        # de comptage n'est PAS tracé -- le nombre de tours n'est pas forcément
        # entier, sa position angulaire n'a donc rien d'imposé et l'afficher
        # comme un repère de construction serait trompeur.
        a = np.radians(angle_collet_deg)
        ax.plot([rayon_min * np.cos(a)], [rayon_min * np.sin(a)], "o",
                color="tab:green", ms=5, zorder=6)
        ax.annotate("collet", (rayon_min * np.cos(a), rayon_min * np.sin(a)),
                    color="tab:green", fontsize=7.5, ha="left", va="bottom", zorder=6,
                    xytext=(4, 4), textcoords="offset points",
                    bbox=dict(facecolor=theme.fond(), edgecolor="none", alpha=0.85, pad=1.0))

    ax.annotate(nom_vue(cote_rouage), (0.5, -0.19), xycoords="axes fraction",
                ha="center", va="top", fontsize=7.5, color=theme.gris(0.35))


def tracer_vecteur_g(ax, rayon_max: float, nom_position: str):
    """Vecteur g, tracé vers le BAS dans la marge gauche de la figure.

    La figure qui porte ce vecteur est dessinée dans l'ATTITUDE de la position
    considérée, c'est-à-dire tournée de sorte que la gravité aille vers le bas.
    Un vecteur g pointant vers le haut n'a pas de sens physique et ne se dessine
    jamais ainsi ; c'est le mouvement qui tourne, pas la pesanteur.

    Tracé dans la MARGE et non au centre : au centre, la flèche se superposait
    aux spires et on ne savait plus si elle appartenait au dessin.
    """
    from .constantes import G_TERRE

    x = -1.34 * rayon_max
    haut, bas = 0.62 * rayon_max, -0.22 * rayon_max
    ax.annotate("", xy=(x, bas), xytext=(x, haut),
                arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=2.2,
                                shrinkA=0, shrinkB=0))
    # Étiquette SOUS la pointe, à l'horizontale : posée le long de la flèche,
    # elle se superposait au trait.
    ax.annotate(f"g = {G_TERRE:.2f} m/s²", (x, bas), color="tab:red", fontsize=7.5,
                ha="center", va="top", xytext=(0, -4), textcoords="offset points")
    ax.annotate(f"Position COSC {nom_position} — repère {nom_position} en haut",
                (0.5, -0.06), xycoords="axes fraction", ha="center", va="top",
                fontsize=8.5, color="tab:red")


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


@theme.clair()
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
    ligne_ext, = ax.plot([], [], "-", color=theme.encre(), lw=0.9)
    ligne_int, = ax.plot([], [], "-", color=theme.encre(), lw=0.9)
    ax.plot(0, 0, "+", color=theme.encre(), ms=10)
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
