# ancien_format.py
"""
Reprise des « expériences » JSON des anciens programmes (historique/claude_v4
à Claude v7, été 2026) : `{"schema_version": 1, "p": {...}, "mode": ..., "x": [...]}`.

Ces fichiers ne décrivent pas un spiral dans la paramétrisation actuelle
(Grossmann intérieure + corps + Phillips extérieure, profils PCHIP) : ils
donnent un vecteur `x` interprété par l'ancien générateur
`generer_spiral_neutre_flexible`, porté tel quel ci-dessous. Le spiral est
donc RECONSTRUIT point par point, exactement comme l'ancien programme le
dessinait, puis repris comme une géométrie importée (au même titre qu'un DXF) :
il s'analyse (marche, jeu, rapport, export) mais ne se ré-optimise pas avec
les panneaux de forme.

L'ancien format ignorait plusieurs réglages apparus depuis (orientation dans
le mouvement, wafer, oxydation, raquette...). `reprendre_experience` rend,
avec la géométrie, la liste de TOUS les champs et leur statut, pour que
l'interface mette en évidence ceux qu'il faut encore valider.

Ancien générateur, une seule version (v5, v6 et v7 ont un physics.py
identique ; v4 n'en diffère que par l'absence du point milieu) :
  - corps : spirale d'Archimède de r_min, pas (r_max - r_min) / (2.pi.n_tours_fictifs),
    sur n_tours_corps tours, épaisseur interpolée linéairement sur 5 nœuds
    équirépartis en angle ;
  - courbe terminale (côté PITON, bien que l'ancien code la nomme
    « Grossmann ») : raccord de Hermite cubique jusqu'à r_max, pente d'entrée
    égale à celle du corps (bridée à 2.95), pente d'arrivée = cambrure ;
    épaisseur linéaire par morceaux sur (début, milieu, fin) ;
  - rééchantillonnage à n_elements points équidistants en abscisse curviligne ;
  - mise à l'échelle de TOUTE l'épaisseur pour tomber exactement sur
    freq_cible, avec E = 169 GPa (silicium isotrope) -- les épaisseurs
    stockées dans x sont donc celles AVANT cette mise à l'échelle.
"""
import json
from dataclasses import dataclass

import numpy as np

from .modele import GeometrieSpiral

#: module de Young de l'ancien programme (historique/*/config.py)
E_ANCIEN = 169e9

#: clés de `p` dans l'ancien format
CLES_P = ('r_min', 'r_max', 'n_tours_reels', 'n_tours_fictifs', 'jeu_min_spires',
          'amplitude_max', 'freq_cible', 'hauteur_h', 'inertie_balancier', 'n_elements')

# Statuts d'un champ, du plus sûr au moins sûr.
REPRIS = "repris"            # même grandeur, même unité
CONVERTI = "converti"        # même grandeur, unité changée
A_VALIDER = "à valider"      # absent de l'ancien format : valeur supposée
IGNORE = "ignoré"            # présent dans l'ancien fichier, sans usage ici


@dataclass
class Champ:
    cle: str
    libelle: str
    valeur: object
    statut: str
    note: str = ""


def est_experience(donnees) -> bool:
    """Vrai si `donnees` (JSON déjà lu) est une ancienne expérience."""
    return (isinstance(donnees, dict) and donnees.get("schema_version") == 1
            and "x" in donnees and "p" in donnees and "spec" not in donnees)


def _decoder_x(mode, x):
    """(5 épaisseurs du corps, n_tours_corps, (début, milieu, fin), cambrure)."""
    x = [float(v) for v in x]
    if mode == "constante":
        if len(x) != 5:
            raise ValueError(f"mode 'constante' : 5 variables attendues, {len(x)} reçues")
        debut, fin = x[2], x[3]
        return [x[0]] * 5, x[1], (debut, (debut + fin) / 2.0, fin), x[4]
    if len(x) == 10:
        return x[0:5], x[5], (x[6], x[7], x[8]), x[9]
    if len(x) == 9:
        # Format v4, sans point milieu : la moyenne redonne exactement
        # l'interpolation linéaire de l'époque.
        return x[0:5], x[5], (x[6], (x[6] + x[7]) / 2.0, x[7]), x[8]
    raise ValueError(f"mode '{mode}' : 9 ou 10 variables attendues, {len(x)} reçues")


def generer_ancienne_geometrie(p, mode, x):
    """Port fidèle de `generer_spiral_neutre_flexible` (historique/Claude v7).

    `p` est le dictionnaire de l'ancien fichier (hauteur_h en MÈTRES).
    Retourne (x, y, epaisseur_mm, diagnostics) en mm : collet sur l'axe +x,
    développement en sens trigonométrique, comme l'ancien programme.
    """
    ep_corps, n_tours_corps, ep_term, cambrure = _decoder_x(mode, x)
    theta_total = p['n_tours_reels'] * 2 * np.pi
    pas = (p['r_max'] - p['r_min']) / (p['n_tours_fictifs'] * 2 * np.pi)
    theta_corps = n_tours_corps * 2 * np.pi

    r_fin_corps = p['r_min'] + pas * theta_corps
    larg_theta = max(theta_total - theta_corps, 1e-9)
    larg_r = max(p['r_max'] - r_fin_corps, 1e-9)
    S_MAX = 2.95
    s0_brut = pas * larg_theta / larg_r
    s0 = min(max(s0_brut, 0.0), S_MAX)
    s1 = min(max(cambrure, 0.0), S_MAX)

    t = np.linspace(0, theta_total, 5000)
    corps = t <= theta_corps
    u = np.where(corps, 0.0, (t - theta_corps) / larg_theta)
    f = (u**3 - 2*u**2 + u) * s0 + (-2*u**3 + 3*u**2) + (u**3 - u**2) * s1
    r_dense = np.where(corps, p['r_min'] + pas * t, r_fin_corps + larg_r * f)
    ep_dense = np.where(corps,
                        np.interp(t, np.linspace(0, theta_corps, 5), ep_corps),
                        np.interp(u, [0.0, 0.5, 1.0], ep_term))

    x_d, y_d = r_dense * np.cos(t), r_dense * np.sin(t)
    s_dense = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(x_d), np.diff(y_d)))))
    s_u = np.linspace(0, s_dense[-1], int(p['n_elements']))
    xf, yf = np.interp(s_u, s_dense, x_d), np.interp(s_u, s_dense, y_d)
    ep = np.interp(s_u, s_dense, ep_dense)

    # Mise à l'échelle de l'épaisseur sur la fréquence cible, à l'identique
    ds = np.sqrt(np.gradient(xf)**2 + np.gradient(yf)**2) * 1e-3
    souplesse = np.sum(ds / (np.clip(ep * 1e-3, 1e-6, None)**3))
    couple = (E_ANCIEN * p['hauteur_h']) / (12.0 * souplesse)
    freq = (1.0 / (2.0 * np.pi)) * np.sqrt(couple / p['inertie_balancier'])
    beta = (p['freq_cible'] / freq) ** (2.0 / 3.0)

    diagnostics = {
        "s0_bride": bool(s0_brut > S_MAX),
        "cambrure_bridee": bool(cambrure != s1),
        "facteur_epaisseur": float(beta),
        "n_tours_corps": float(n_tours_corps),
    }
    return xf, yf, ep * beta, diagnostics


def _angle_collet(spec):
    """Même règle que ParametresSpiral.angle_collet_deg."""
    angle = spec.angle_echappement_deg + (0.0 if spec.collet_cote_ancre else 180.0)
    return (angle + 180.0) % 360.0 - 180.0


def reprendre_experience(chemin, spec):
    """Lit une ancienne expérience et la transpose.

    `spec` (presets.SpecMouvement) porte les réglages que l'ancien format ne
    connaissait pas, tels qu'ils sont dans l'interface : orientation dans le
    mouvement, wafer. Retourne (spec repris, géométrie, champs, diagnostics).
    """
    from dataclasses import replace
    with open(chemin, encoding="utf-8") as f:
        d = json.load(f)
    if not est_experience(d):
        raise ValueError("ce fichier n'est pas une ancienne expérience (schema_version 1)")
    p = dict(d["p"])
    manquantes = [c for c in CLES_P if c not in p]
    if manquantes:
        raise ValueError(f"clés manquantes dans 'p' : {', '.join(manquantes)}")

    # L'interface v4-v7 saisissait la hauteur en mm et la stockait en m
    # (gui.py : sp_h.value()*1e-3). Le garde-fou couvre un fichier écrit à la main.
    h_m = float(p['hauteur_h'])
    h_mm = h_m * 1e3 if h_m < 0.01 else h_m
    p['hauteur_h'] = h_mm * 1e-3

    x0, y0, ep, diag = generer_ancienne_geometrie(p, d.get("mode", "variable"), d["x"])

    reprise = replace(spec, r_min=float(p['r_min']), r_max=float(p['r_max']),
                      n_tours_reels=float(p['n_tours_reels']),
                      n_tours_fictifs=float(p['n_tours_fictifs']),
                      hauteur_h=h_mm, freq_cible=float(p['freq_cible']),
                      inertie_balancier=float(p['inertie_balancier']),
                      jeu_min_spires=float(p['jeu_min_spires']),
                      amplitude_max=float(p['amplitude_max']))

    # Orientation : l'ancien spiral partait de l'axe +x en sens trigonométrique,
    # sans référence au mouvement. Il est posé ici comme le serait un spiral
    # généré avec les réglages courants.
    # Sens horaire : symétrie par rapport à l'axe du collet, puis rotation.
    if not reprise.sens_trigonometrique:
        y0 = -y0
    a = np.radians(_angle_collet(reprise))
    geom = GeometrieSpiral(x=x0 * np.cos(a) - y0 * np.sin(a),
                           y=x0 * np.sin(a) + y0 * np.cos(a), epaisseur=ep,
                           hauteur=h_mm, plan_cristallin=reprise.plan_cristallin,
                           orientation_cristal_deg=reprise.orientation_cristal_deg)

    champs = [
        Champ("r_min", "Rayon virole r_min (mm)", reprise.r_min, REPRIS),
        Champ("r_max", "Rayon piton r_max (mm)", reprise.r_max, REPRIS),
        Champ("n_tours_reels", "Tours réels", reprise.n_tours_reels, REPRIS),
        Champ("n_tours_fictifs", "Tours fictifs", reprise.n_tours_fictifs, REPRIS),
        Champ("freq_cible", "Fréquence cible (Hz)", reprise.freq_cible, REPRIS),
        Champ("inertie_balancier", "Inertie du balancier (kg.m²)",
              reprise.inertie_balancier, REPRIS),
        Champ("jeu_min_spires", "Jeu minimal entre spires (mm)", reprise.jeu_min_spires, REPRIS),
        Champ("amplitude_max", "Amplitude maximale (°)", reprise.amplitude_max, REPRIS),
        Champ("hauteur_h", "Hauteur de lame (mm)", h_mm, CONVERTI,
              f"stockée en mètres dans l'ancien fichier ({h_m:g})"),
        Champ("x", "Géométrie (vecteur x)", len(d["x"]), CONVERTI,
              "spiral reconstruit point par point par l'ancien générateur, "
              f"{int(p['n_elements'])} éléments ; analysable, pas ré-optimisable "
              "avec les panneaux de forme"),
        Champ("angle_echappement_deg", "Angle de la ligne d'échappement (°)",
              reprise.angle_echappement_deg, A_VALIDER,
              "l'ancien spiral n'était pas placé dans le mouvement : valeur de l'interface ; "
              "le signe du DVH en dépend"),
        Champ("collet_cote_ancre", "Collet du côté de l'ancre",
              reprise.collet_cote_ancre, A_VALIDER, "idem"),
        Champ("sens_trigonometrique", "Développement trigonométrique",
              reprise.sens_trigonometrique, A_VALIDER,
              "l'ancien spiral tournait en sens trigonométrique, sans dire de quel côté il était vu"),
        Champ("vue_cote_rouage", "Vue côté rouage", reprise.vue_cote_rouage, A_VALIDER,
              "absent de l'ancien format"),
        Champ("plan_cristallin", "Plan cristallin du wafer", reprise.plan_cristallin, A_VALIDER,
              "l'ancien calcul supposait E = 169 GPa partout, ce que seul un wafer (111) donne"),
        Champ("orientation_cristal_deg", "Orientation du cristal (°)",
              reprise.orientation_cristal_deg, A_VALIDER, "sans effet sur un wafer (111)"),
        Champ("oxydation", "Oxydation / compensation thermique", "aucune", A_VALIDER,
              "l'ancien spiral était en silicium nu"),
        Champ("raquette", "Terminaison à raquette", "aucune", A_VALIDER,
              "l'ancien spiral se terminait librement au piton"),
        Champ("resultats", "Résultats enregistrés", d.get("resultats"), IGNORE,
              "calculés par l'ancien modèle : à refaire au panneau Résultats"),
        Champ("grossmann_fixe", "Paramètres terminaux figés", d.get("grossmann_fixe"), IGNORE,
              "point de départ de l'ancienne optimisation ; c'est x qui fait foi"),
    ]
    if diag["s0_bride"]:
        champs.append(Champ("s0", "Raccord corps / courbe terminale", "bridé", A_VALIDER,
                            "l'ancien générateur n'a pas pu tenir la continuité de pente"))
    return reprise, geom, champs, diag


if __name__ == "__main__":
    # Contrôle en ligne de commande : python -m spiral.ancien_format fichier.json...
    import sys
    from .physique import frequence_propre
    from .presets import SpecMouvement
    for chemin in sys.argv[1:]:
        spec, geom, champs, diag = reprendre_experience(chemin, SpecMouvement())
        f = frequence_propre(geom, spec.inertie_balancier)
        print(f"{chemin}: {len(geom.x)} points, épaisseur {geom.epaisseur.min()*1e3:.1f}-"
              f"{geom.epaisseur.max()*1e3:.1f} µm, rayon {np.hypot(geom.x, geom.y).max():.4f} mm, "
              f"fréquence recalculée {f:.5f} Hz (cible {spec.freq_cible})")
        for c in champs:
            if c.statut == A_VALIDER:
                print(f"   à valider : {c.libelle} = {c.valeur}")
