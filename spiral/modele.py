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

    @property
    def n(self) -> int:
        return len(self.x)

    @property
    def longueur(self) -> float:
        return float(self.s[-1])

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

    def contour_flancs(self):
        """Retourne (flanc_ext, flanc_int), chacun un array (N,2), obtenus en
        décalant la fibre neutre de ±epaisseur/2 le long de la normale locale.
        Sert à l'export DXF et à l'affichage."""
        _, _, nx, ny = self.tangente_normale()
        demi_e = self.epaisseur / 2.0
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
        return GeometrieSpiral(x=x_new, y=y_new, epaisseur=e_new, hauteur=self.hauteur)

    def masse_elementaire(self, rho: float) -> np.ndarray:
        """Masse (kg) de chaque élément, longueurs/épaisseurs en mm -> m."""
        ds = np.gradient(self.s) * 1e-3
        return rho * self.hauteur * 1e-3 * (self.epaisseur * 1e-3) * ds

    def masse_totale(self, rho: float) -> float:
        return float(np.sum(self.masse_elementaire(rho)))
