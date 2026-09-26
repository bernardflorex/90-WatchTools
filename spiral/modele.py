# modele.py
"""
Représentation canonique unique d'un spiral, utilisée par TOUT le reste du
programme (import DXF, génération paramétrique, physique, optimisation,
export DXF). Un seul type d'objet, quelle que soit l'origine de la géométrie.
"""
from dataclasses import dataclass, field
import numpy as np


@dataclass
class GeometrieSpiral:
    x: np.ndarray           # mm, fibre neutre
    y: np.ndarray           # mm, fibre neutre
    epaisseur: np.ndarray   # mm, par élément (peut varier librement)
    hauteur: float          # mm, hors-plan, constante sur tout le spiral (DRIE)
    # Nombre d'éléments ÉLASTIQUEMENT ACTIFS, en partant du collet. Vaut n (tout
    # le spiral) sauf sur un spiral à raquette : la lame y est saisie entre les
    # goupilles au POINT DE COMPTAGE, et la portion qui va de là jusqu'au piton
    # ne participe plus au rappel -- ses deux extrémités étant fixes, elle ne se
    # déforme pas. Elle reste présente dans la géométrie, car elle occupe bien
    # de la place (contrôle du jeu inter-spires) et doit être usinée (export DXF).
    n_actif: int = None
    # Épaisseur de la couche d'OXYDE thermique (mm), nulle par défaut. Quand
    # elle n'est pas nulle, `epaisseur` reste la cote de GRAVURE (silicium nu,
    # ce que doit produire la DRIE) et la lame finie est plus large de
    # 2 x 0.56 x t, tandis que l'oxydation a consommé 2 x 0.44 x t de silicium.
    # Toute la physique utilise alors la rigidité du COMPOSITE et l'encombrement
    # de la lame FINIE : c'est la pièce oxydée qui oscille et qui doit passer
    # entre les spires.
    epaisseur_oxyde: float = 0.0
    # Convention d'unités du module : TOUT est en mm (x, y, epaisseur, hauteur,
    # s), comme un DXF. La conversion en mètres (SI) se fait localement, au
    # point d'usage, dans les formules physiques (cf. physique.py,
    # parametrage.py) -- jamais en stockant hauteur/epaisseur en mètres ici.
    s: np.ndarray = field(init=False)      # abscisse curviligne cumulée (mm)
    theta: np.ndarray = field(init=False)  # angle polaire déroulé (rad, np.unwrap)

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        self.epaisseur = np.asarray(self.epaisseur, dtype=float)
        n = len(self.x)
        if len(self.y) != n or len(self.epaisseur) != n:
            raise ValueError(
                f"x, y, epaisseur doivent avoir la même longueur "
                f"(reçu {n}, {len(self.y)}, {len(self.epaisseur)})"
            )
        ds = np.hypot(np.diff(self.x, prepend=self.x[0]), np.diff(self.y, prepend=self.y[0]))
        ds[0] = 0.0
        self.s = np.cumsum(ds)
        self.theta = np.unwrap(np.arctan2(self.y, self.x))
        if self.n_actif is None or self.n_actif > n:
            self.n_actif = n
        elif self.n_actif < 4:
            raise ValueError(f"n_actif={self.n_actif} : partie active trop courte")

    @property
    def n(self) -> int:
        return len(self.x)

    @property
    def longueur(self) -> float:
        return float(self.s[-1])

    @property
    def a_partie_morte(self) -> bool:
        return self.n_actif < self.n

    def partie_active(self) -> "GeometrieSpiral":
        """Sous-géométrie réellement sollicitée en flexion, du collet au point
        de comptage. C'est elle, et elle seule, que voit le solveur elastica :
        la portion morte ne porte aucun moment et ne se déforme pas.

        Exclure la portion morte du centre de gravité est EXACT, pas approché,
        pour le calcul de marche : la formule d'Airy n'utilise que la dérivée
        du moment statique (masse x position) par rapport à l'angle du
        balancier, et celui de la portion morte est constant."""
        if not self.a_partie_morte:
            return self
        return GeometrieSpiral(
            x=self.x[:self.n_actif], y=self.y[:self.n_actif],
            epaisseur=self.epaisseur[:self.n_actif], hauteur=self.hauteur,
            epaisseur_oxyde=self.epaisseur_oxyde,
        )

    @property
    def epaisseur_finie(self) -> np.ndarray:
        """Épaisseur physique de la lame après oxydation (mm) : c'est elle qui
        occupe l'espace entre spires et qui porte le poids."""
        from .constantes import CROISSANCE_OXYDATION
        return self.epaisseur + 2.0 * CROISSANCE_OXYDATION * self.epaisseur_oxyde

    def rigidite_flexion(self) -> np.ndarray:
        """Rigidité de flexion E·I (N.m²) dans le plan, à la température de
        référence. Réduit EXACTEMENT à E_silicium·I quand l'épaisseur d'oxyde
        est nulle : le spiral nu n'est pas un cas particulier traité à part."""
        from .thermique import ReglageThermique, rigidite_composite
        reglage = ReglageThermique()
        return rigidite_composite(self.epaisseur * 1e-3, self.hauteur * 1e-3,
                                   self.epaisseur_oxyde * 1e-3,
                                   reglage.temperature_reference, reglage)

    def souplesse(self) -> float:
        """Souplesse intégrée de la PARTIE ACTIVE, integrale(ds / E·I) en
        1/(N.m). Son inverse est le couple de rappel par radian."""
        from scipy.integrate import trapezoid
        a = self.partie_active()
        return float(trapezoid(1.0 / a.rigidite_flexion(), x=a.s * 1e-3))

    def rayon(self) -> np.ndarray:
        return np.hypot(self.x, self.y)

    def tangente_normale(self):
        """Tangente/normale unitaires locales par différences finies (np.gradient),
        cohérent avec la convention déjà utilisée par le solveur elastica
        (voir physique.py) — la normale pointe du côté du flanc extérieur
        (grand rayon)."""
        dx, dy = np.gradient(self.x), np.gradient(self.y)
        norme = np.hypot(dx, dy)
        norme[norme == 0] = 1e-9
        tx, ty = dx / norme, dy / norme
        nx, ny = -ty, tx
        # Orientation : la normale doit pointer vers l'extérieur (grand rayon)
        r = self.rayon()
        r[r == 0] = 1e-9
        signe = np.sign(nx * self.x / r + ny * self.y / r)
        signe[signe == 0] = 1.0
        return tx, ty, nx * signe, ny * signe

    def contour_flancs(self, cote_gravure: bool = False):
        """Retourne (flanc_ext, flanc_int), chacun un array (N,2), obtenus en
        décalant la fibre neutre de ±epaisseur/2 le long de la normale locale.

        Par défaut : la lame FINIE (oxyde compris), c'est-à-dire la pièce telle
        qu'elle oscille -- c'est cette cote qui décide du jeu entre spires et
        de ce qu'on voit. Avec cote_gravure=True : la cote de gravure DRIE,
        avant oxydation, qui est ce dont l'atelier a besoin pour le masque."""
        _, _, nx, ny = self.tangente_normale()
        demi_e = (self.epaisseur if cote_gravure else self.epaisseur_finie) / 2.0
        ext = np.column_stack((self.x + demi_e * nx, self.y + demi_e * ny))
        inte = np.column_stack((self.x - demi_e * nx, self.y - demi_e * ny))
        return ext, inte

    def resampler(self, n_elements: int) -> "GeometrieSpiral":
        """Ré-échantillonne à n_elements points uniformément répartis en
        abscisse curviligne. Utilisé pour uniformiser la densité de
        discrétisation physique quelle que soit l'origine (DXF ou génération
        paramétrique)."""
        s_cible = np.linspace(0.0, self.s[-1], n_elements)
        x_new = np.interp(s_cible, self.s, self.x)
        y_new = np.interp(s_cible, self.s, self.y)
        e_new = np.interp(s_cible, self.s, self.epaisseur)
        # Le point de comptage est repéré par son abscisse curviligne, pas par
        # son indice : un ré-échantillonnage ne doit pas le déplacer le long de
        # la lame (ce serait un changement de longueur active, donc de fréquence).
        n_actif = None
        if self.a_partie_morte:
            n_actif = int(np.searchsorted(s_cible, self.s[self.n_actif - 1]) + 1)
            n_actif = int(np.clip(n_actif, 4, n_elements))
        return GeometrieSpiral(x=x_new, y=y_new, epaisseur=e_new, hauteur=self.hauteur,
                                n_actif=n_actif, epaisseur_oxyde=self.epaisseur_oxyde)

    def masse_elementaire(self, rho: float) -> np.ndarray:
        """Masse (kg) de chaque élément, longueurs/épaisseurs en mm -> m.

        `rho` est la masse volumique du silicium ; celle de l'oxyde est prise
        dans les constantes. L'oxyde étant un peu moins dense que le silicium
        qu'il remplace (2200 contre 2330 kg/m³) et la section grossissant,
        l'effet net sur la masse n'est pas négligeable quand la couche fait
        plusieurs microns."""
        from .constantes import CONSOMMATION_OXYDATION, RHO_OXYDE
        ds = np.gradient(self.s) * 1e-3
        h = self.hauteur * 1e-3
        t = self.epaisseur_oxyde * 1e-3
        if t <= 0.0:
            return rho * h * (self.epaisseur * 1e-3) * ds
        e_coeur = np.maximum(self.epaisseur * 1e-3 - 2.0 * CONSOMMATION_OXYDATION * t, 0.0)
        h_coeur = max(h - 2.0 * CONSOMMATION_OXYDATION * t, 0.0)
        aire_coeur = e_coeur * h_coeur
        aire_totale = (self.epaisseur_finie * 1e-3) * (h + 2.0 * (1.0 - CONSOMMATION_OXYDATION) * t)
        return (rho * aire_coeur + RHO_OXYDE * np.maximum(aire_totale - aire_coeur, 0.0)) * ds

    def masse_totale(self, rho: float) -> float:
        return float(np.sum(self.masse_elementaire(rho)))
