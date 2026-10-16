"""Paramètres physiques de l'échappement (ce qui n'est pas dans le tracé).

Unités SI à l'intérieur ; l'interface affiche les unités horlogères usuelles
(µN.m, g.mm², mg.cm²). Les valeurs par défaut sont celles de l'échappement en
silicium étudié : couple 1,13 µN.m sur la roue, frottement silicium/silicium
0,15, balancier de 8 mg.cm² à 3,5 Hz, facteur de qualité 220 en position
horizontale (frottement de l'air compris).

Les coefficients de restitution valent 0 par défaut (chocs mous) : aucun
rebond n'est visible à la caméra rapide sur l'échappement réel, alors qu'avec
0,5 l'ellipse rebondissait une dizaine de fois dans la fourchette. Aucune
mesure publiée n'a été trouvée pour un choc silicium sur silicium dans un
échappement ; le matériau seul rendrait presque tout, ce sont les vibrations
des pièces élancées et le couplage au rouage qui dissipent.
"""

from dataclasses import dataclass, asdict

import numpy as np


@dataclass
class Parametres:
    couple_roue: float = 1.13e-6         # N.m
    mu: float = 0.15                     # frottement dent-palette et ellipse-fourchette
    e_chute: float = 0.0                 # restitution dent-palette
    e_fourchette: float = 0.0            # restitution ellipse-fourchette
    e_butee: float = 0.0                 # restitution fourchette-butée
    inertie_roue: float = 2.12e-11       # kg.m² (roue acier 0,021 + pignon 0,0002 g.mm²)
    inertie_rouage: float = 7e-13        # kg.m² (0,0007 g.mm²), rouage ramené à la roue d'échappement
    inertie_ancre: float = 8.14e-12      # kg.m² (0,00814 g.mm², dard compris)
    inertie_balancier: float = 8e-10     # kg.m² (8 mg.cm², comme le programme spiral)
    frequence: float = 3.5               # Hz
    amplitude_deg: float = 310.0         # amplitude visée (demi-oscillation)
    facteur_qualite: float = 220.0       # balancier seul, air compris (horizontal, 250-310°)

    @property
    def raideur(self):
        """Couple de rappel du spiral, N.m/rad."""
        return self.inertie_balancier * (2 * np.pi * self.frequence) ** 2

    @property
    def amortissement(self):
        """Coefficient d'amortissement visqueux du balancier, N.m.s/rad,
        tel que le facteur de qualité soit celui donné."""
        return self.inertie_balancier * 2 * np.pi * self.frequence / self.facteur_qualite

    def energie_balancier(self, amplitude_deg=None):
        a = np.radians(self.amplitude_deg if amplitude_deg is None else amplitude_deg)
        return 0.5 * self.raideur * a * a

    def perte_par_oscillation(self, amplitude_deg=None):
        """Énergie dissipée par le balancier seul en une période (Q)."""
        return 2 * np.pi * self.energie_balancier(amplitude_deg) / self.facteur_qualite

    def en_dict(self):
        return asdict(self)
