#!/usr/bin/env python3
"""Suite de non-régression de spiral_lab.

Lancer :  python tests/test_spiral_lab.py
Sortie : une ligne par vérification, code de retour non nul si l'une échoue.

Les vérifications sont regroupées par sujet et ordonnées du plus élémentaire au
plus intégré : une référence analytique d'abord (ressort linéaire, formule de
Duffing), puis la physique du spiral, puis les chaînes complètes. Un échec en
tête rend les suivants ininterprétables.

Plusieurs de ces contrôles ont attrapé de vrais défauts : le jeu inter-spires
calculé sur la cote de gravure au lieu de la lame oxydée, la souplesse totale
employée là où seule la partie active compte, la fréquence non rejouée après le
calcul d'épaisseur d'oxyde.
"""
import os
import sys
import tempfile
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spiral_lab import constantes as cst
from spiral_lab import dxf_io
from spiral_lab.echappement import ReglageEchappement, marche_echappement
from spiral_lab.isochronisme import (POSITIONS_HORIZONTALES, analyser_marches,
                                      marche_isochronisme_intrinseque)
from spiral_lab.optimisation import (construire_bornes, evaluer_candidat,
                                      lancer_optimisation, params_depuis_vecteur,
                                      sonder_faisabilite, vecteur_depuis_params)
from spiral_lab.modele import GeometrieSpiral
from spiral_lab.parametrage import ParametresSpiral, construire_geometrie
from spiral_lab.physique import (angle_point_comptage, axe_force_pivot,
                                  balayage_angulaire, frequence_propre,
                                  marche_raquette, resoudre_elastica)
from spiral_lab.presets import (BornesOptimisation, HyperparametresGA, SpecMouvement,
                                 charger_configuration, sauvegarder_configuration)
from spiral_lab.raquette import (excursion_laterale_um_par_rad, facteur_zone_morte)
from spiral_lab.thermique import (ReglageThermique, analyser_thermique, epaisseur_oxyde_optimale,
                                   etude_sensibilite, part_oxyde, pente_thermique,
                                   rigidite_composite)

ECHECS = []


def check(nom, condition, info=""):
    etat = "OK " if condition else "ECHEC"
    if not condition:
        ECHECS.append(nom)
    print(f"  [{etat}] {nom} {info}")


def _leve(exception, fonction, *args, **kwargs):
    """Vrai si l'appel lève bien l'exception attendue. Un refus explicite vaut
    mieux qu'un résultat silencieusement faux : on le vérifie donc comme une
    fonctionnalité, pas comme un accident."""
    try:
        fonction(*args, **kwargs)
    except exception:
        return True
    except Exception:
        return False
    return False


def titre(t):
    print(f"\n== {t} ==")


SPEC = SpecMouvement()
J = SPEC.inertie_balancier
AMPS = np.arange(200.0, 341.0, 20.0)
BASE = dict(r_min=SPEC.r_min, r_max=SPEC.r_max, n_tours_reels=SPEC.n_tours_reels,
            n_tours_fictifs=SPEC.n_tours_fictifs, hauteur_h=SPEC.hauteur_h,
            freq_cible=SPEC.freq_cible, inertie_balancier=SPEC.inertie_balancier)
REG_E = ReglageEchappement()
REG_T = ReglageThermique()


def params(**extra):
    return ParametresSpiral(**BASE, **extra)


# ---------------------------------------------------------------- A
titre("A. Isochronisme : références analytiques")
k = 4 * np.pi ** 2 * 3.5 ** 2 * J
ang = np.linspace(-340.0, 340.0, 69)
m, d = marche_isochronisme_intrinseque(ang, 0.5 * k * np.radians(ang) ** 2, J, AMPS, 300.0)
check("ressort linéaire -> isochrone exact", np.max(np.abs(m)) < 1e-6,
      f"({np.max(np.abs(m)):.2e} s/j)")
check("fréquence des petites oscillations exacte",
      abs(d["frequence_petites_oscillations_hz"] - 3.5) < 1e-9)
b4 = 1e-4 * k
m2, _ = marche_isochronisme_intrinseque(
    ang, 0.5 * k * np.radians(ang) ** 2 + 0.25 * b4 * np.radians(ang) ** 4, J, AMPS, 300.0)
attendu = 86400 * (3 * b4 / (8 * k)) * (np.radians(200.) ** 2 - np.radians(300.) ** 2)
check("ressort durcissant : accord Duffing 1er ordre < 1%",
      abs(m2[0] - attendu) / abs(attendu) < 0.01, f"({m2[0]:+.3f} vs {attendu:+.3f} s/j)")

# ---------------------------------------------------------------- B
titre("B. Échappement")
me = marche_echappement(AMPS, REG_E)
check("nul à l'amplitude de référence", abs(np.interp(300., AMPS, me)) < 1e-12)
check("-retard à l'amplitude basse", abs(np.interp(200., AMPS, me) + 8.0) < 1e-12)
check("d'autant plus fort que l'amplitude est faible", bool(np.all(np.diff(me) > 0)))

# ---------------------------------------------------------------- C
titre("C. Composition des marches")
g = construire_geometrie(params())
r = analyser_marches(g, SPEC.r_min, 3.5, J, AMPS, reglage_echappement=REG_E, pas_balayage_deg=30.)
check("les deux horizontales présentes et identiques",
      all(p in r["totales"] for p in POSITIONS_HORIZONTALES)
      and np.allclose(r["totales"]["CH"], r["totales"]["CB"]))
check("horizontale = isochronisme + échappement (gravité nulle à plat)",
      np.allclose(r["totales"]["CH"], r["isochronisme"] + r["echappement"]))
check("verticale = gravité + isochronisme + échappement",
      np.allclose(r["totales"]["3H"], r["gravite"]["3H"] + r["isochronisme"] + r["echappement"]))
r30 = analyser_marches(g, SPEC.r_min, 3.5, J, AMPS,
                       reglage_echappement=ReglageEchappement(retard_s_par_j=30.),
                       pas_balayage_deg=30.)
check("l'échappement ne crée aucune dispersion entre positions",
      abs(r["dispersion_gravite"] - r30["dispersion_gravite"]) < 1e-9)

# ---------------------------------------------------------------- C bis
titre("C bis. Convention de vue et positions verticales")
from spiral_lab.isochronisme import (POSITIONS_VERTICALES_COTE_CADRAN,
                                      POSITIONS_VERTICALES_COTE_ROUAGE,
                                      critere_dispersion_marche,
                                      marche_positions_verticales, nom_vue,
                                      positions_verticales)


def _sens(positions):
    """Sens de rotation de 12-3-6-9, en tenant compte du passage par 0."""
    ang = [np.degrees(positions[n]) % 360 for n in ("12H", "3H", "6H", "9H")]
    pas = [(ang[i + 1] - ang[i]) % 360 for i in range(3)]
    if all(abs(x - 90) < 1e-9 for x in pas):
        return "antihoraire"
    if all(abs(x - 270) < 1e-9 for x in pas):
        return "horaire"
    return "?"


check("côté rouage : 12-3-6-9 antihoraire (cadran vu de dos)",
      _sens(POSITIONS_VERTICALES_COTE_ROUAGE) == "antihoraire")
check("côté cadran : 12-3-6-9 horaire", _sens(POSITIONS_VERTICALES_COTE_CADRAN) == "horaire")
check("le défaut est le côté rouage (le spiral est sous le coq)",
      positions_verticales() == POSITIONS_VERTICALES_COTE_ROUAGE)
check("changer de côté échange 3H et 9H",
      POSITIONS_VERTICALES_COTE_ROUAGE["3H"] == POSITIONS_VERTICALES_COTE_CADRAN["9H"]
      and POSITIONS_VERTICALES_COTE_ROUAGE["9H"] == POSITIONS_VERTICALES_COTE_CADRAN["3H"])
check("6H et 12H ne bougent pas",
      POSITIONS_VERTICALES_COTE_ROUAGE["6H"] == POSITIONS_VERTICALES_COTE_CADRAN["6H"]
      and POSITIONS_VERTICALES_COTE_ROUAGE["12H"] == POSITIONS_VERTICALES_COTE_CADRAN["12H"])
m_r = marche_positions_verticales(g, SPEC.r_min, 3.5, J, AMPS,
                                   positions=POSITIONS_VERTICALES_COTE_ROUAGE,
                                   pas_balayage_deg=40.)
m_c = marche_positions_verticales(g, SPEC.r_min, 3.5, J, AMPS,
                                   positions=POSITIONS_VERTICALES_COTE_CADRAN,
                                   pas_balayage_deg=40.)
check("la dispersion ne dépend PAS du côté (c'est une permutation)",
      abs(critere_dispersion_marche(m_r) - critere_dispersion_marche(m_c)) < 1e-12,
      f"({critere_dispersion_marche(m_r):.6f} s/j)")
check("seules 3H et 9H échangent leurs valeurs",
      np.allclose(m_r["3H"], m_c["9H"]) and np.allclose(m_r["9H"], m_c["3H"])
      and np.allclose(m_r["6H"], m_c["6H"]) and np.allclose(m_r["12H"], m_c["12H"]))
check("la convention est nommée sur les graphes",
      "rouage" in nom_vue(True) and "cadran" in nom_vue(False))
check("la convention est portée par la spec et sauvegardée",
      SpecMouvement().vue_cote_rouage is True)

# --- repère du mouvement et DVH ---------------------------------------
# Le piège corrigé le 2026-09-29 : le nom d'une position COSC désigne le repère
# horaire qui est EN HAUT, donc la direction OPPOSÉE à la gravité. Confondre les
# deux décale tout l'étiquetage de 180° -- une erreur muette, qui ne change ni
# la dispersion ni l'optimisation, mais attribue les valeurs aux mauvaises
# positions dès la comparaison à un relevé.
from spiral_lab.isochronisme import reperes_horaires
for _cote in (True, False):
    _rh, _pv = reperes_horaires(_cote), positions_verticales(_cote)
    check(f"repère horaire et g diamétralement opposés (côté rouage={_cote})",
          all(abs(((_pv[n] - _rh[n]) % (2 * np.pi)) - np.pi) < 1e-12
              for n in ("12H", "3H", "6H", "9H")))
# Orientation du tracé, telle que l'utilisateur regarde ses montres : côté
# rouage, 6H en haut, tige à droite. Un premier jet dessinait 12H en haut, donc
# tige à gauche -- chiralité correcte mais repère tourné d'un demi-tour, ce qui
# ne se voit que sur le signe du DVH (corrigé le 2026-09-29).
check("le repère 6H est en haut du tracé",
      abs(reperes_horaires(True)["6H"] - np.pi / 2) < 1e-12)
check("la tige sort au repère 3H, à droite vu côté rouage",
      abs(reperes_horaires(True)["3H"]) < 1e-12)
check("12H en bas et 9H à gauche",
      abs(reperes_horaires(True)["12H"] - 3 * np.pi / 2) < 1e-12
      and abs(reperes_horaires(True)["9H"] - np.pi) < 1e-12)
# La chiralité, elle, ne change pas : vue de dos, 12-3-6-9 tourne dans le sens
# antihoraire. Testée sur les angles DÉROULÉS, une comparaison brute échouant
# au passage par 360°.
_suite = np.unwrap([reperes_horaires(True)[n] for n in ("12H", "3H", "6H", "9H")])
check("12-3-6-9 tourne en sens antihoraire côté rouage",
      bool(np.all(np.diff(_suite) > 0)), f"({np.degrees(_suite).round(0)})")
check("côté cadran, le sens s'inverse",
      bool(np.all(np.diff(np.unwrap(
          [reperes_horaires(False)[n] for n in ("12H", "3H", "6H", "9H")])) < 0)))
# En position COSC 6H, le 6 est en haut : le tracé par défaut EST déjà dans
# cette attitude, donc g y pointe vers le bas sans rotation.
check("le tracé par défaut est déjà dans l'attitude COSC 6H",
      abs(positions_verticales(True)["6H"] - 3 * np.pi / 2) < 1e-12)

check("la spec porte l'angle d'échappement, le sens et le bout de la ligne",
      abs(SpecMouvement().angle_echappement_deg + 6.2) < 1e-12
      and SpecMouvement().sens_trigonometrique is True
      and SpecMouvement().collet_cote_ancre is True)
# L'angle de la ligne des centres ne la définit qu'à 180° près : le collet part
# du bout qui pointe vers l'ancre, donc à l'angle saisi lui-même, et non à son
# opposé. C'est exactement le demi-tour corrigé le 2026-09-29.
_p_ref = params()
check("le collet part à l'angle de la ligne des centres (-6.2°, pas 173.8°)",
      abs(_p_ref.angle_collet_deg - _p_ref.angle_echappement_deg) < 1e-12,
      f"({_p_ref.angle_collet_deg:+.2f}°)")
check("l'autre bout de la ligne met le collet à 180° de là",
      abs(params(collet_cote_ancre=False).angle_collet_deg
          - (180.0 + _p_ref.angle_echappement_deg)) < 1e-12)
# angle_point_comptage est un angle d'AXE, rendu modulo 180° pour être comparé
# à celui de la force au pivot : 173.8° y désigne la même droite que -6.2°. Un
# premier jet comparait ce modulo 180 à un angle modulo 360 et concluait à un
# demi-tour inexistant -- c'est ce qui faisait dessiner le repère du point de
# comptage à l'opposé de la courbe de Phillips (corrigé le 2026-09-29).
_g_ref = construire_geometrie(_p_ref)
_ecart = (angle_point_comptage(_g_ref) - _p_ref.angle_collet_deg) % 180.0
check("le point de comptage est sur la ligne des centres (à 180° près)",
      min(_ecart, 180.0 - _ecart) < 1e-6,
      f"(comptage {angle_point_comptage(_g_ref):+.2f}° mod 180, "
      f"collet {_p_ref.angle_collet_deg:+.2f}°)")
# La direction de l'ancre ne dépend pas du bout qui porte le collet : le
# demi-axe dessiné doit rester le même quand on retourne le spiral.
check("la direction de l'ancre est indépendante du bout qui porte le collet",
      abs(params(collet_cote_ancre=False).angle_ancre_deg
          - _p_ref.angle_ancre_deg) < 1e-12,
      f"({_p_ref.angle_ancre_deg:+.2f}°)")
check("l'ancre est à l'opposé de la tige, sur la ligne des centres",
      abs(_p_ref.angle_ancre_deg - _p_ref.angle_echappement_deg) < 1e-12)
# Le demi-tour ne change pas la forme : le module de sensibilité est identique,
# seul le signe du DVH bascule. C'est ce qui rendait l'erreur muette.
_r_a = analyser_marches(_g_ref, SPEC.r_min, 3.5, J, AMPS, pas_balayage_deg=40.)
_r_b = analyser_marches(construire_geometrie(params(collet_cote_ancre=False)),
                        SPEC.r_min, 3.5, J, AMPS, pas_balayage_deg=40.)
check("le demi-tour laisse le module de sensibilité inchangé",
      abs(_r_a["sensibilite_module_reference"]
          - _r_b["sensibilite_module_reference"]) < 1e-9)
check("le demi-tour inverse le signe du DVH",
      abs(_r_a["dvh_reference"] + _r_b["dvh_reference"]) < 1e-9,
      f"({_r_a['dvh_reference']:+.4f} contre {_r_b['dvh_reference']:+.4f} s/j)")

# Une rotation du mouvement ne change pas la FORME du spiral : le module du
# vecteur de sensibilité à la gravité doit être rigoureusement invariant, seule
# sa phase tourne. C'est ce qui rend le DVH réglable par orientation -- mais
# seulement en le reportant sur les autres positions, jamais en le supprimant.
_mods, _dvhs = [], []
for _a in (-6.2, 0.0, 45.0, 90.0):
    _r = analyser_marches(construire_geometrie(params(angle_echappement_deg=_a)),
                          SPEC.r_min, 3.5, J, AMPS, pas_balayage_deg=40.)
    _mods.append(_r["sensibilite_module_reference"])
    _dvhs.append(_r["dvh_reference"])
    check(f"DVH = A.cos(90° - phi) à {_a:+.1f}°",
          abs(_r["dvh_reference"] - _r["sensibilite_module_reference"]
              * np.cos(np.radians(90.0 - _r["sensibilite_phase_deg"]))) < 1e-12)
check("le module de sensibilité est invariant par rotation",
      (max(_mods) - min(_mods)) / max(_mods) < 1e-6,
      f"({min(_mods):.9f} à {max(_mods):.9f} s/j)")
check("le DVH, lui, dépend de l'orientation",
      max(_dvhs) - min(_dvhs) > 0.1, f"({min(_dvhs):+.4f} à {max(_dvhs):+.4f} s/j)")

# La vue d'attitude tourne le mouvement pour que g aille vers le BAS : on ne
# dessine jamais la pesanteur vers le haut, c'est la montre qu'on incline.
# Vérifié ici sur la rotation elle-même, pas sur l'image : le repère qui donne
# son nom à la position doit finir en haut de la figure.
for _nom in ("6H", "9H", "12H", "3H"):
    _rot = np.radians(270.0) - positions_verticales(True)[_nom]
    _a_repere = reperes_horaires(True)[_nom] + _rot
    check(f"position COSC {_nom} : le repère {_nom} est en haut, g vers le bas",
          np.sin(_a_repere) > 0.999
          and np.sin(positions_verticales(True)[_nom] + _rot) < -0.999)

# ---------------------------------------------------------------- D
titre("D. Solveur elastica")
*_, diag = resoudre_elastica(g, 280.0, SPEC.r_min, return_diagnostics=True)
# Seuils relâchés d'une décade (2026-09-29) : le spiral est désormais construit
# TOURNÉ dans le repère du mouvement (collet sur la ligne d'échappement, à
# ~174°) au lieu de partir de 0°. Les coordonnées passent donc par des sinus et
# cosinus non triviaux, et le résidu du tir hérite de quelques ULP de plus.
# 1e-5 µm, c'est 10 picomètres sur une pièce de 2.7 mm : le seuil précédent
# mesurait la virgule flottante, pas la convergence.
check("résidu de position négligeable", diag["erreur_position_um"] < 1e-5,
      f"({diag['erreur_position_um']:.2e} µm)")
check("résidu d'angle négligeable", diag["erreur_angle_deg"] < 1e-7,
      f"({diag['erreur_angle_deg']:.2e} deg)")
a_, x_, y_, j_, U_ = balayage_angulaire(g, SPEC.r_min, 340.0, 40.0)
# La fréquence propre n'était vérifiée que par AUTO-COHÉRENCE : la mise à
# l'échelle analytique de l'épaisseur (§4.5) la ramène toujours sur la cible,
# donc elle tombe sur 3,5 Hz quoi qu'il arrive. Cela ne dit rien de la valeur
# ABSOLUE. Le contrôle ci-dessous la confronte à la formule fermée d'une poutre
# droite uniforme, C = E·h·e³/(12·L), seul cas où la raideur se calcule sans
# intégrale. Ajouté le 2026-10-08, après un doute sur un spiral importé dont la
# fréquence sortait 19 % au-dessus de la valeur attendue : le contrôle montre
# que la formule est exacte et que l'écart vient des données d'entrée.
for _L_mm, _e_um, _h_um, _J_mg in ((100.0, 36.0, 100.0, 8.0), (50.0, 20.0, 150.0, 5.0),
                                    (119.707, 35.99, 100.0, 8.0)):
    _n = 4000
    _g_dr = GeometrieSpiral(x=np.linspace(1.0, 1.0 + _L_mm, _n), y=np.zeros(_n),
                            epaisseur=np.full(_n, _e_um * 1e-3), hauteur=_h_um * 1e-3)
    _C = cst.E_SILICIUM * (_h_um * 1e-6) * (_e_um * 1e-6) ** 3 / (12.0 * _L_mm * 1e-3)
    _f_th = np.sqrt(_C / (_J_mg * 1e-10)) / (2 * np.pi)
    _f_pr = frequence_propre(_g_dr, _J_mg * 1e-10)
    check(f"fréquence absolue exacte sur poutre droite "
          f"(L={_L_mm:g} mm, e={_e_um:g} µm, h={_h_um:g} µm)",
          abs(_f_pr - _f_th) / _f_th < 1e-12,
          f"({_f_pr:.5f} Hz pour {_f_th:.5f} attendus)")

check("énergie strictement croissante depuis l'équilibre",
      bool(np.all(np.diff(U_[len(U_) // 2:]) > 0)))
check("fréquence propre sur la cible", abs(frequence_propre(g, J) - 3.5) < 1e-6,
      f"({frequence_propre(g, J):.7f} Hz)")
module, axe = axe_force_pivot(g, SPEC.r_min, 340.0)
poids = g.masse_totale(cst.RHO_SILICIUM) * 9.81 * 1e6
check("force sur le pivot très supérieure au poids du spiral", module > 5 * poids,
      f"({module:.1f} µN pour {poids:.1f} µN de poids)")
check("son axe suit la droite collet -> point de comptage",
      abs((axe - angle_point_comptage(g) + 90.) % 180. - 90.) < 10.0,
      f"(axe {axe:.1f} deg, comptage {angle_point_comptage(g):.1f} deg)")

# ---------------------------------------------------------------- E
titre("E. Terminaison à raquette")
b_raq = params(raquette_active=True)
g_raq, d_raq = construire_geometrie(b_raq, return_diagnostics=True)
n_avant, n_pit = b_raq.n_pas_arc_avant, b_raq.n_pas_piton
check("nombre d'éléments conforme", g_raq.n == b_raq.n_elements,
      f"({g_raq.n} vs {b_raq.n_elements})")
check("arc à rayon constant = r_max",
      np.allclose(g_raq.rayon()[-(n_avant + n_pit):], b_raq.r_max, atol=1e-9))
th = np.degrees(g_raq.theta)
check("nombre ENTIER de tours du collet au point de comptage",
      abs(th[g_raq.n_actif - 1] - th[0] - b_raq.n_tours_reels * 360.) < 1e-6,
      f"({th[g_raq.n_actif-1]-th[0]:.3f} deg)")
check("point de comptage au milieu de l'arc quand les deux angles sont égaux",
      abs(b_raq.raquette_arc_avant_deg - b_raq.raquette_piton_angle_deg) < 1e-9)
check("matière au-delà du point de comptage", g_raq.a_partie_morte,
      f"({g_raq.n - g_raq.n_actif} éléments)")
check("course des goupilles contenue dans l'arc",
      d_raq["raquette_parcours_tient_dans_arc"],
      f"(demi-course {d_raq['raquette_demi_parcours_deg']:.2f} deg)")
_, d_serre = construire_geometrie(params(raquette_active=True, raquette_arc_avant_deg=2.0,
                                          raquette_plage_s_par_j=150.), return_diagnostics=True)
check("arc trop court -> signalé", not d_serre["raquette_parcours_tient_dans_arc"])

Xe0, Ye0, *_ = resoudre_elastica(g_raq, 0.0, b_raq.r_min)
Xe2, Ye2, *_ = resoudre_elastica(g_raq, 300.0, b_raq.r_min)
check("bords rendus sur la pièce entière", len(Xe2) == g_raq.n)
dep = float(np.hypot(Xe2[g_raq.n_actif:] - Xe0[g_raq.n_actif:],
                     Ye2[g_raq.n_actif:] - Ye0[g_raq.n_actif:]).max())
check("portion morte strictement immobile", dep < 1e-12, f"({dep*1e9:.1e} nm)")

for plage in (30.0, 60.0, 150.0):
    gp = construire_geometrie(params(raquette_active=True, raquette_plage_s_par_j=plage))
    etendue = marche_raquette(gp, J)[2]
    check(f"plage ±{plage:.0f} s/j obtenue à mieux que 0.1%",
          abs(etendue - 2 * plage) / (2 * plage) < 1e-3, f"(étendue {etendue:.3f})")
etendues = [marche_raquette(construire_geometrie(
    params(raquette_active=True, raquette_arc_avant_n_pas=n)), J)[2] for n in (20, 100, 400)]
check("plage insensible au maillage de l'arc", max(etendues) - min(etendues) < 1e-6,
      f"({[round(e, 4) for e in etendues]})")

g_libre = construire_geometrie(params())
check("spiral libre : aucune partie morte", not g_libre.a_partie_morte)
check("marche_raquette neutre sur un spiral libre",
      marche_raquette(g_libre, J) == (0.0, 0.0, 0.0))

# ---------------------------------------------------------------- F
titre("F. Jeu entre les goupilles")
check("fonction descriptive N(0)=1", facteur_zone_morte(0.0) == 1.0)
check("fonction descriptive N(1)=0", abs(facteur_zone_morte(1.0)) < 1e-15)
check("fonction descriptive monotone décroissante",
      bool(np.all(np.diff(facteur_zone_morte(np.linspace(0, 1, 201))) <= 1e-15)))


def essai_jeu(piton, jeu):
    bb = params(raquette_active=True, raquette_piton_angle_deg=piton, raquette_jeu_um=jeu)
    return analyser_marches(construire_geometrie(bb), bb.r_min, 3.5, J, AMPS,
                            reglage_echappement=REG_E, pas_balayage_deg=40., jeu_goupilles_um=jeu)


r0 = essai_jeu(30., 0.)
check("jeu nul -> marche de jeu identiquement nulle (= piton fixe)",
      float(np.max(np.abs(r0["jeu_goupilles"]))) == 0.0)
vals = [essai_jeu(30., x)["jeu_amplitude_basse"] for x in (5., 10., 20., 40.)]
check("plus de jeu -> plus de retard", all(vals[i] > vals[i + 1] for i in range(3)),
      f"({[round(v, 1) for v in vals]})")
check("le jeu produit un RETARD", all(v < 0 for v in vals))
exc = [analyser_marches(construire_geometrie(
    params(raquette_active=True, raquette_piton_angle_deg=p_, raquette_jeu_um=20.)),
    SPEC.r_min, 3.5, J, AMPS, reglage_echappement=REG_E, pas_balayage_deg=40.,
    jeu_goupilles_um=20.)["diag_jeu_goupilles"]["excursion_um_par_rad"] for p_ in (15., 30., 60.)]
check("excursion croissante avec la longueur après goupilles", exc[0] < exc[1] < exc[2],
      f"({[round(e, 2) for e in exc]} µm/rad)")
check("excursion superlinéaire (flexion en L²)", exc[1] / exc[0] > 2.0)
r20 = essai_jeu(30., 20.)
check("le jeu ne crée aucune dispersion entre positions verticales",
      abs(r20["dispersion_gravite"] - r0["dispersion_gravite"]) < 1e-9)
for variante, fraction in (("piton", 0.0), ("milieu", 0.5), ("raquette", 1.0)):
    reg = ReglageEchappement(variante_visee=variante)
    check(f"montage visé « {variante} » -> fraction {fraction}",
          reg.fraction_jeu() == fraction)
check("excursion nulle sur un spiral libre",
      excursion_laterale_um_par_rad(g_libre, SPEC.r_min) == 0.0)

# ---------------------------------------------------------------- G
titre("G. Compensation thermique")
e_m, h_m = 30.8e-6, 100e-6
EI0 = rigidite_composite(e_m, h_m, 0.0, REG_T.temperature_reference, REG_T)
check("oxyde nul -> E_Si·I exactement",
      abs(EI0 / (cst.E_SILICIUM * h_m * e_m ** 3 / 12) - 1) < 1e-14)
check("part d'oxyde nulle quand t=0", part_oxyde(e_m, h_m, 0.0, REG_T) == 0.0)
p_nu = pente_thermique(g_libre, 0.0, J, REG_T)
check("silicium nu fortement non compensé", -5.0 < p_nu < -2.5, f"({p_nu:+.3f} s/j/°C)")
a_th = analyser_thermique(g_libre, J, REG_T)
check("épaisseur d'oxyde dans la plage attendue", 2.0 < a_th["epaisseur_oxyde_um"] < 6.0,
      f"({a_th['epaisseur_oxyde_um']:.3f} µm)")
check("pente annulée", abs(a_th["pente_s_par_j_par_degc"]) < 1e-6)
check("croissance = 2·0.56·t",
      abs(a_th["croissance_lame_um"] - 2 * cst.CROISSANCE_OXYDATION * a_th["epaisseur_oxyde_um"]) < 1e-12)
for nom, extra in (("nu", {}), ("oxyde auto", dict(oxydation_active=True)),
                   ("oxyde imposé 3.0", dict(oxydation_active=True, oxyde_impose_um=3.0)),
                   ("oxyde imposé 5.0", dict(oxydation_active=True, oxyde_impose_um=5.0)),
                   ("oxyde + raquette", dict(oxydation_active=True, raquette_active=True))):
    f = frequence_propre(construire_geometrie(params(**extra)), J)
    check(f"{nom:18s} fréquence sur la cible", abs(f - 3.5) < 1e-6, f"({f:.7f} Hz)")
g_ox = construire_geometrie(params(oxydation_active=True))
check("épaisseur finie = gravure + 2·0.56·t",
      np.allclose(g_ox.epaisseur_finie,
                  g_ox.epaisseur + 2 * cst.CROISSANCE_OXYDATION * g_ox.epaisseur_oxyde))
jeu_nu = analyser_marches(g_libre, SPEC.r_min, 3.5, J, AMPS,
                          pas_balayage_deg=60.)["diag_orbite"]["jeu_min_mm"]
jeu_ox = analyser_marches(g_ox, SPEC.r_min, 3.5, J, AMPS,
                          pas_balayage_deg=60.)["diag_orbite"]["jeu_min_mm"]
check("l'oxyde réduit le jeu inter-spires de sa croissance",
      abs((jeu_nu - jeu_ox) * 1000 - 2 * cst.CROISSANCE_OXYDATION * g_ox.epaisseur_oxyde * 1000) < 1.5,
      f"({jeu_nu*1000:.1f} -> {jeu_ox*1000:.1f} µm)")
sens = etude_sensibilite(g_libre, J, REG_T)
check("cinq coefficients étudiés", len(sens) == 5)
check("le coefficient de l'oxyde domine l'incertitude",
      max(sens, key=lambda l: abs(l["pente_non_recalee_plus"]))["coefficient"] == "tce_oxyde")

# ---------------------------------------------------------------- H
titre("H. Export / import DXF")
with tempfile.TemporaryDirectory() as td:
    ch = dxf_io.exporter_dxf_paire(g_ox, os.path.join(td, "s.dxf"))
    check("deux fichiers quand la lame est oxydée", len(ch) == 2,
          f"({[os.path.basename(c) for c in ch]})")
    gg = dxf_io.charger_dxf(ch[0], n_elements=1500, hauteur_h=g_ox.hauteur)
    gf = dxf_io.charger_dxf(ch[1], n_elements=1500, hauteur_h=g_ox.hauteur)
    check("le fichier _gravure relit la cote de gravure",
          abs(gg.epaisseur.mean() - g_ox.epaisseur.mean()) * 1000 < 0.2)
    check("le fichier _fini relit la cote finie",
          abs(gf.epaisseur.mean() - g_ox.epaisseur_finie.mean()) * 1000 < 0.2)
    check("un seul fichier sans oxyde",
          len(dxf_io.exporter_dxf_paire(g_libre, os.path.join(td, "nu.dxf"))) == 1)
    dxf_io.exporter_dxf(g_libre, os.path.join(td, "r.dxf"))
    g2 = dxf_io.charger_dxf(os.path.join(td, "r.dxf"), n_elements=1500, hauteur_h=g_libre.hauteur)
    r1, r2 = g_libre.rayon(), g2.rayon()
    # Seuil ramené de 15 µm à 0,5 µm (2026-10-07). Les ~12 µm d'erreur au collet
    # étaient attribués à une « limite connue » du pipeline d'export/import :
    # c'était en réalité la tolérance de recollement de 0,1 mm du chaînage, qui
    # soudait les deux flancs puis les rescindait par extrema de rayon. Avec une
    # tolérance adaptative, l'aller-retour tombe à 0,04 µm.
    check("r_min conservé à mieux que 0,5 µm",
          abs(r1.min() - r2.min()) * 1000 < 0.5, f"({abs(r1.min()-r2.min())*1000:.3f} µm)")
    check("r_max conservé à mieux que 0,5 µm",
          abs(r1.max() - r2.max()) * 1000 < 0.5, f"({abs(r1.max()-r2.max())*1000:.3f} µm)")
    check("épaisseur conservée à mieux que 3 µm",
          abs(g_libre.epaisseur.mean() - g2.epaisseur.mean()) * 1000 < 3)

    # --- ré-échantillonnage : même COURBE, autre nombre de points ---
    # La densité qu'exige un fichier de gravure n'a aucune raison d'être celle
    # du maillage de calcul. Ce qui doit être préservé, c'est la géométrie --
    # pas le découpage.
    for g_src, nom in ((g_libre, "spiral libre"), (g_raq, "spiral à raquette")):
        for n_cible in (600, 5000):
            r = g_src.reechantillonner(n_cible)
            check(f"ré-échantillonnage à {n_cible} points ({nom})", r.n == n_cible,
                  f"({g_src.n} -> {r.n})")
            # La longueur développée diminue légèrement quand on enlève des
            # points : une polyligne à cordes coupe les virages. L'écart doit
            # rester sous le pour mille, et DÉCROÎTRE quand on raffine.
            ecart = abs(r.longueur - g_src.longueur) / g_src.longueur
            check(f"longueur conservée à 0,1 % près ({n_cible}, {nom})", ecart < 1e-3,
                  f"({100*ecart:.4f} %)")
            check(f"épaisseur dans la plage d'origine ({n_cible}, {nom})",
                  r.epaisseur.min() >= g_src.epaisseur.min() - 1e-12
                  and r.epaisseur.max() <= g_src.epaisseur.max() + 1e-12)
            check(f"rayons extrêmes conservés à 1 µm ({n_cible}, {nom})",
                  abs(r.rayon().min() - g_src.rayon().min()) < 1e-3
                  and abs(r.rayon().max() - g_src.rayon().max()) < 1e-3)
    # Le point de comptage est une frontière PHYSIQUE (les goupilles y saisissent
    # la lame) : il doit rester un point du maillage, jamais tomber entre deux.
    r = g_raq.reechantillonner(5000)
    i0, i1 = g_raq.n_actif - 1, r.n_actif - 1
    check("le point de comptage reste exactement sur le maillage",
          abs(r.x[i1] - g_raq.x[i0]) < 1e-12 and abs(r.y[i1] - g_raq.y[i0]) < 1e-12)
    check("le ré-échantillonnage conserve la partie morte",
          r.a_partie_morte and 0 < r.n_actif < r.n, f"(n_actif {r.n_actif}/{r.n})")
    check("ré-échantillonner au même nombre de points ne change rien",
          g_libre.reechantillonner(g_libre.n) is g_libre)

    # --- export DXF à densité imposée ---
    import ezdxf as _ezdxf
    ch_dense = os.path.join(td, "dense.dxf")
    dxf_io.exporter_dxf(g_libre, ch_dense, n_points=900)
    pts = [len(e) for e in _ezdxf.readfile(ch_dense).modelspace()
           if e.dxftype() == "LWPOLYLINE"]
    check("l'export DXF respecte la densité demandée", pts == [900, 900], f"({pts})")
    g3 = dxf_io.charger_dxf(ch_dense, n_elements=1500, hauteur_h=g_libre.hauteur)
    check("un DXF ré-échantillonné décrit toujours la même pièce",
          abs(g3.rayon().max() - g_libre.rayon().max()) < 5e-3
          and abs(g3.epaisseur.mean() - g_libre.epaisseur.mean()) * 1000 < 3,
          f"(r_max {g3.rayon().max():.4f} vs {g_libre.rayon().max():.4f} mm)")
    # --- DXF exporté en SEGMENTS ÉLÉMENTAIRES, avec géométrie accessoire ---
    # Cas d'un vrai export CAO (Spiral-K2-2025.dxf, 6 087 LINE) : aucune
    # polyligne, et le dessin porte en plus le cercle de virole et une ligne de
    # construction partant du centre. Deux pièges s'y cumulaient : un chaînage
    # en coût cubique (plusieurs minutes, l'application paraissait figée) et une
    # tolérance de recollement assez lâche pour souder la ligne de construction
    # au spiral, ce qui faussait le rayon de collet sans rien signaler.
    import ezdxf as _ez, time as _time
    ext_s, int_s = g_libre.contour_flancs()
    doc_s = _ez.new(setup=False)
    doc_s.header["$INSUNITS"] = 4
    msp_s = doc_s.modelspace()
    for contour in (ext_s, int_s):
        for a, b in zip(contour[:-1], contour[1:]):
            msp_s.add_line(tuple(a), tuple(b), dxfattribs={"layer": "ModelSketch_Visible"})
    # Accessoires : virole et ligne de construction vers le centre, à moins de
    # 0,1 mm du collet — c'est exactement ce qui piégeait l'ancienne tolérance.
    msp_s.add_circle((0.0, 0.0), 0.15, dxfattribs={"layer": "ModelSketch_Visible"})
    msp_s.add_line((0.0, 0.0), tuple(int_s[0] * 0.6), dxfattribs={"layer": "ModelSketch_Visible"})
    ch_seg = os.path.join(td, "segments.dxf")
    doc_s.saveas(ch_seg)

    _t0 = _time.time()
    g_seg = dxf_io.charger_dxf(ch_seg, n_elements=1500, hauteur_h=g_libre.hauteur)
    _dt = _time.time() - _t0
    check("un DXF en segments élémentaires s'importe en moins de 30 s",
          _dt < 30.0, f"({_dt:.2f} s pour {2*(len(ext_s)-1)} segments)")
    r_seg = g_seg.rayon()
    check("la ligne de construction vers le centre n'est pas prise pour le spiral",
          abs(r_seg.min() - g_libre.rayon().min()) * 1000 < 20,
          f"(r_min {r_seg.min():.4f} contre {g_libre.rayon().min():.4f} mm)")
    check("le cercle de virole n'est pas pris pour un flanc",
          abs(r_seg.max() - g_libre.rayon().max()) * 1000 < 20,
          f"(r_max {r_seg.max():.4f} contre {g_libre.rayon().max():.4f} mm)")
    check("épaisseur correcte depuis un DXF en segments",
          abs(g_seg.epaisseur.mean() - g_libre.epaisseur.mean()) * 1000 < 3,
          f"({g_seg.epaisseur.mean()*1000:.2f} contre {g_libre.epaisseur.mean()*1000:.2f} µm)")

    # --- DXF FRAGMENTÉ : trous le long des flancs, de la taille de la lame ---
    # Cas de spiral-111.dxf : le contour y est complet mais cassé en 151
    # morceaux, avec des trous de ~38 µm le long de chaque flanc — exactement
    # l'écart entre les deux flancs. À chaque coupure, deux candidats se
    # présentent donc à la MÊME distance : la vraie suite du flanc et le flanc
    # d'en face. La distance seule ne peut pas trancher, et élargir la
    # tolérance ne fait qu'aggraver le saut d'un flanc à l'autre. Seule la
    # continuité de DIRECTION les sépare (1 à 2° contre 90 ou 180°).
    # Le contour est émis comme UN SEUL parcours — flanc extérieur, bouchon au
    # piton, flanc intérieur en sens inverse — puis fragmenté : c'est la forme
    # qu'a le fichier réel. Les segments font la taille de la lame, donc chaque
    # trou vaut deux fois l'écart entre flancs : le mauvais candidat est plus
    # PROCHE que le bon, et seule la direction permet de choisir.
    _g_fin = g_libre.reechantillonner(int(g_libre.longueur / g_libre.epaisseur.mean()))
    _ext_f, _int_f = _g_fin.contour_flancs()
    _contour_f = np.vstack([_ext_f, _int_f[::-1]])
    _doc_f = _ez.new(setup=False)
    _doc_f.header["$INSUNITS"] = 4
    _msp_f = _doc_f.modelspace()
    _trous = 0
    for _k in range(len(_contour_f) - 1):
        if _k % 40 == 39:              # un segment sur quarante est omis
            _trous += 1
            continue
        _msp_f.add_line(tuple(_contour_f[_k]), tuple(_contour_f[_k + 1]))
    _ch_f = os.path.join(td, "fragmente.dxf")
    _doc_f.saveas(_ch_f)
    check("le DXF d'essai est bien fragmenté",
          _trous > 100, f"({_trous} coupures, lame de {g_libre.epaisseur.mean()*1000:.1f} µm)")

    _g_f = dxf_io.charger_dxf(_ch_f, n_elements=2000, hauteur_h=g_libre.hauteur)
    _t_f = np.unwrap(np.arctan2(_g_f.y, _g_f.x))
    _t_src = np.unwrap(np.arctan2(g_libre.y, g_libre.x))
    check("un contour fragmenté se recolle sans perdre de matière",
          abs(_g_f.longueur - g_libre.longueur) / g_libre.longueur < 0.02,
          f"({_g_f.longueur:.2f} mm pour {g_libre.longueur:.2f} dessinés)")
    check("le nombre de tours est conservé malgré les coupures",
          abs(abs(_t_f[-1] - _t_f[0]) - abs(_t_src[-1] - _t_src[0])) < np.radians(20.0),
          f"({abs(_t_f[-1]-_t_f[0])/(2*np.pi):.3f} contre "
          f"{abs(_t_src[-1]-_t_src[0])/(2*np.pi):.3f} tours)")
    check("le chaînage n'a pas sauté d'un flanc à l'autre",
          abs(np.median(_g_f.epaisseur) - np.median(g_libre.epaisseur)) * 1000 < 3.0,
          f"({np.median(_g_f.epaisseur)*1000:.2f} contre "
          f"{np.median(g_libre.epaisseur)*1000:.2f} µm)")

    # --- Deux pièges d'un DXF de CAO réel, reproduits analytiquement ---
    # 1. La courbe terminale RENTRE sur les derniers degrés : le rayon maximal
    #    n'est donc PAS au bout du ruban mais avant. Chercher les plis aux
    #    extrema de rayon décale la découpe des deux flancs, et l'intersection
    #    angulaire ampute le spiral de toute la portion au-delà. Mesuré sur
    #    Spiral-K2-2025.dxf : 102° perdus au piton, silencieusement.
    # 2. Une PORTION RECTILIGNE de la courbe terminale est dessinée en quelques
    #    longues cordes là où le reste est fin. L'interpolation de r(θ) y
    #    devient fausse et l'épaisseur s'effondre.
    _th = np.linspace(0.0, 16.0 * np.pi, 4000)
    _rc = 0.50 + 0.26 * _th / (2 * np.pi)
    _hook = _th > _th[-1] - np.radians(120.0)      # la terminaison rentre
    _rc[_hook] -= 0.10 * ((_th[_hook] - _th[_hook][0]) / np.radians(120.0)) ** 2
    _ep = 0.035
    _int = np.column_stack([(_rc - _ep / 2) * np.cos(_th), (_rc - _ep / 2) * np.sin(_th)])
    _ext = np.column_stack([(_rc + _ep / 2) * np.cos(_th), (_rc + _ep / 2) * np.sin(_th)])
    check("le ruban d'essai a bien son rayon max AVANT le bout",
          int(np.argmax(np.hypot(*_ext.T))) < len(_th) - 10,
          f"(argmax r à {int(np.argmax(np.hypot(*_ext.T)))}/{len(_th)})")

    _doc_r = _ez.new(setup=False)
    _doc_r.header["$INSUNITS"] = 4
    _msp_r = _doc_r.modelspace()
    for _contour in (_int, _ext):
        _i = 0
        while _i < len(_contour) - 1:
            # Portion rectiligne grossière entre 60 % et 64 % du parcours :
            # une seule corde là où le reste est échantillonné finement.
            _pas = 40 if 0.60 < _i / len(_contour) < 0.64 else 1
            _j = min(_i + _pas, len(_contour) - 1)
            _msp_r.add_line(tuple(_contour[_i]), tuple(_contour[_j]))
            _i = _j
    # Bouchon au bout du ruban, pour que le contour soit un ruban fermé au piton
    _msp_r.add_line(tuple(_int[-1]), tuple(_ext[-1]))
    _ch_r = os.path.join(td, "terminaison_rentrante.dxf")
    _doc_r.saveas(_ch_r)

    _g_r = dxf_io.charger_dxf(_ch_r, n_elements=2000, hauteur_h=0.100)
    _t_imp = np.unwrap(np.arctan2(_g_r.y, _g_r.x))
    _etendue = abs(np.degrees(_t_imp[-1] - _t_imp[0]))
    _attendu = abs(np.degrees(_th[-1] - _th[0]))
    check("l'étendue angulaire est conservée (pas de troncature au piton)",
          _attendu - _etendue < 15.0,
          f"({_etendue:.0f}° importés pour {_attendu:.0f}° dessinés, "
          f"manque {_attendu - _etendue:.0f}°)")
    check("le rayon maximal du ruban est atteint",
          abs(_g_r.rayon().max() - _rc.max()) * 1000 < 25,
          f"({_g_r.rayon().max():.4f} contre {_rc.max():.4f} mm)")
    # L'épaisseur doit rester plate : c'est la portion rectiligne grossière qui
    # la faisait s'effondrer.
    _e_r = _g_r.epaisseur[int(0.05 * _g_r.n):int(0.95 * _g_r.n)] * 1000
    check("pas d'effondrement d'épaisseur sur la portion rectiligne grossière",
          _e_r.min() > 0.80 * _ep * 1000,
          f"({_e_r.min():.2f} µm minimum pour {_ep*1000:.1f} µm dessinés)")

    check("la paire de DXF accepte aussi une densité",
          all(len(e) == 700 for c in dxf_io.exporter_dxf_paire(
                  g_ox, os.path.join(td, "p.dxf"), n_points=700)
              for e in _ezdxf.readfile(c).modelspace() if e.dxftype() == "LWPOLYLINE"))

# --- anisotropie élastique du silicium -------------------------------
# E_SILICIUM = 169 GPa n'est pas une propriété du silicium : c'est sa valeur
# dans le plan d'un wafer (111), le SEUL plan isotrope dans le plan. Prendre
# cette valeur sur un (001) surestime la raideur de 13 % (remonté par
# l'utilisateur, 2026-10-08).
from spiral_lab.constantes import PLANS_CRISTALLINS, module_young_silicium

# Valeurs de référence du silicium monocristallin, indépendantes de ce projet.
_E100 = module_young_silicium("001", 0.0)
_E110 = module_young_silicium("001", np.pi / 4)
check("E⟨100⟩ = 130,2 GPa", abs(_E100 / 1e9 - 130.2) < 0.15, f"({_E100/1e9:.1f} GPa)")
check("E⟨110⟩ = 168,9 GPa", abs(_E110 / 1e9 - 168.9) < 0.15, f"({_E110/1e9:.1f} GPa)")
# ⟨111⟩ appartient au plan (110) — il y fait 54,74° avec ⟨001⟩, pas 45° — et
# c'est la direction la plus raide du cristal, donc le maximum sur ce plan.
_E111 = module_young_silicium("110", np.arctan(np.sqrt(2.0)))
check("E⟨111⟩ = 187,5 GPa", abs(_E111 / 1e9 - 187.5) < 0.2, f"({_E111/1e9:.1f} GPa)")
check("⟨111⟩ est bien la direction la plus raide du plan (110)",
      abs(_E111 - module_young_silicium("110", np.linspace(0, np.pi, 3601)).max())
      / _E111 < 1e-6)

_angles = np.linspace(0.0, 2 * np.pi, 721)
check("le plan (111) est le seul isotrope dans le plan",
      np.ptp(module_young_silicium("111", _angles)) < 1.0
      and np.ptp(module_young_silicium("001", _angles)) > 3e10
      and np.ptp(module_young_silicium("110", _angles)) > 3e10)
check("le (111) isotrope vaut bien la constante historique du projet",
      abs(module_young_silicium("111", 0.0) - cst.E_SILICIUM) / cst.E_SILICIUM < 1e-3,
      f"({module_young_silicium('111', 0.0)/1e9:.2f} contre {cst.E_SILICIUM/1e9:.0f} GPa)")
check("le (001) a bien une période de 90°",
      abs(module_young_silicium("001", 0.3) - module_young_silicium("001", 0.3 + np.pi / 2))
      / _E100 < 1e-12)
check("un plan inconnu est refusé explicitement",
      _leve(ValueError, module_young_silicium, "100", 0.0))

# Le défaut "111" doit reproduire EXACTEMENT le comportement antérieur : sans
# cela, tous les résultats antérieurs du projet deviendraient incomparables.
_g_iso = construire_geometrie(params())
check("le défaut reste le (111), comportement antérieur inchangé",
      SPEC.plan_cristallin == "111" and params().plan_cristallin == "111"
      and np.ptp(_g_iso.module_young()) < 1.0)
# Sur un wafer anisotrope la lame est plus souple : à fréquence imposée, la
# mise à l'échelle doit l'épaissir. Et elle doit tomber sur la cible malgré
# l'absence de formule fermée pour la souplesse.
for _plan in PLANS_CRISTALLINS:
    _g_p = construire_geometrie(params(plan_cristallin=_plan))
    check(f"fréquence cible atteinte sur un wafer ({_plan})",
          abs(frequence_propre(_g_p, J) - 3.5) < 1e-6,
          f"({frequence_propre(_g_p, J):.7f} Hz, lame {_g_p.epaisseur.mean()*1000:.3f} µm)")
_e_iso = construire_geometrie(params(plan_cristallin="111")).epaisseur.mean()
_e_ani = construire_geometrie(params(plan_cristallin="001")).epaisseur.mean()
check("un wafer (001), plus souple, exige une lame plus épaisse",
      _e_ani > _e_iso, f"({_e_ani*1000:.3f} contre {_e_iso*1000:.3f} µm)")
# L'orientation du cristal ne doit rien changer sur un plan isotrope, et
# quelque chose sur un plan qui ne l'est pas.
check("l'orientation du cristal est sans effet sur un (111)",
      abs(construire_geometrie(params(plan_cristallin="111",
                                      orientation_cristal_deg=37.0)).epaisseur.mean()
          - _e_iso) < 1e-12)
check("l'orientation du cristal agit sur un (001)",
      abs(construire_geometrie(params(plan_cristallin="001",
                                      orientation_cristal_deg=20.0)).epaisseur.mean()
          - _e_ani) > 1e-9)
# Le module local suit la TANGENTE : sur une poutre droite en (001) orientée
# selon <100>, il doit valoir exactement E<100>.
_n_dr = 500
_g_dr100 = GeometrieSpiral(x=np.linspace(1.0, 11.0, _n_dr), y=np.zeros(_n_dr),
                           epaisseur=np.full(_n_dr, 0.030), hauteur=0.100,
                           plan_cristallin="001")
check("sur une poutre droite selon ⟨100⟩, le module local est E⟨100⟩",
      abs(np.median(_g_dr100.module_young()) - _E100) / _E100 < 1e-9,
      f"({np.median(_g_dr100.module_young())/1e9:.2f} GPa)")
check("la même poutre, cristal tourné de 45°, voit E⟨110⟩",
      abs(np.median(GeometrieSpiral(
          x=np.linspace(1.0, 11.0, _n_dr), y=np.zeros(_n_dr),
          epaisseur=np.full(_n_dr, 0.030), hauteur=0.100,
          plan_cristallin="001", orientation_cristal_deg=45.0).module_young())
          - _E110) / _E110 < 1e-9)

# Le plan cristallin doit aussi atteindre le calcul THERMIQUE : un silicium
# plus souple laisse à l'oxyde une part plus grande de la rigidité, donc il en
# faut MOINS pour annuler la dérive. Un premier jet n'avait câblé l'anisotropie
# que dans la raideur mécanique, et la compensation thermique restait
# insensible au plan — même épaisseur d'oxyde requise sur (111) et sur (001),
# ce qui ne peut pas être (corrigé le 2026-10-08).
_t_ox_iso = epaisseur_oxyde_optimale(construire_geometrie(params(plan_cristallin="111")),
                                      J, REG_T)
_t_ox_ani = epaisseur_oxyde_optimale(construire_geometrie(params(plan_cristallin="001")),
                                      J, REG_T)
check("la compensation thermique dépend du plan cristallin",
      abs(_t_ox_ani - _t_ox_iso) > 0.05,
      f"({_t_ox_ani:.3f} µm sur (001) contre {_t_ox_iso:.3f} sur (111))")
check("un wafer plus souple demande MOINS d'oxyde",
      _t_ox_ani < _t_ox_iso)

# --- unité d'affichage de l'inertie ---------------------------------
# L'inertie vaut ~8e-10 kg.m^2 : saisie telle quelle, elle est illisible et on
# se trompe d'un facteur dix sans le voir. Elle est donc AFFICHÉE en mg.cm^2,
# mais tout le reste du programme reste en SI. Ce qu'on vérifie ici est la
# propriété qui compte : la conversion ne vit qu'aux deux frontières de
# l'interface et l'aller-retour ne dérive pas.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from spiral_lab.app import (_J_KGM2_VERS_MGCM2, _inertie_affichee, _inertie_interne)

check("1 kg·m² vaut 1e10 mg·cm²", _J_KGM2_VERS_MGCM2 == 1.0e10)
check("l'inertie par défaut s'affiche 8 mg·cm²",
      abs(_inertie_affichee(SpecMouvement().inertie_balancier) - 8.0) < 1e-12,
      f"({_inertie_affichee(SpecMouvement().inertie_balancier):.3f} mg·cm²)")
# Aller-retour EN PASSANT PAR L'ARRONDI DU CHAMP (6 décimales) : c'est le vrai
# trajet d'une valeur saisie, et c'est là qu'une dérive se logerait.
_pires = []
for _j in (1e-11, 8e-10, 1.234e-9, 5.5e-10, 1e-7, 9.999e-8):
    _retour = _inertie_interne(round(_inertie_affichee(_j), 3))
    _pires.append(abs(_retour - _j) / _j)
check("aller-retour kg·m² -> mg·cm² -> kg·m² exact à l'arrondi du champ",
      max(_pires) < 1e-15, f"(écart relatif max {max(_pires):.1e})")
check("la résolution du champ reste celle d'avant (1e-13 kg·m²)",
      abs(_inertie_interne(1e-3) - 1e-13) / 1e-13 < 1e-12,
      f"({_inertie_interne(1e-3):.3e} kg·m²)")
# La spec sauvegardée doit rester en SI : un fichier relu par un autre outil,
# ou par une version antérieure, ne doit pas se retrouver avec un facteur 1e7.
with tempfile.TemporaryDirectory() as _td:
    _f = os.path.join(_td, "spec.json")
    SpecMouvement().sauvegarder(_f)
    import json as _json
    check("le JSON de spec reste en kg·m² (SI)",
          abs(_json.load(open(_f))["inertie_balancier"] - 8e-10) < 1e-20)

# ---------------------------------------------------------------- I
titre("I. Configuration complète")
with tempfile.TemporaryDirectory() as td:
    f = os.path.join(td, "c.json")
    b_tout = params(raquette_active=True, raquette_plage_s_par_j=45., oxydation_active=True)
    sauvegarder_configuration(f, SPEC, b_tout, BornesOptimisation(), HyperparametresGA(),
                              ReglageEchappement(retard_s_par_j=11., variante_visee="piton"))
    sp2, b2, bo2, hp2, ech2 = charger_configuration(f)
    check("cinq objets relus", all(x is not None for x in (sp2, b2, bo2, hp2, ech2)))
    check("raquette, oxydation et montage conservés",
          b2.raquette_active and b2.oxydation_active and b2.raquette_plage_s_par_j == 45.
          and ech2.variante_visee == "piton")
    import json
    data = json.load(open(f))
    data["parametres"]["champ_disparu_d_une_ancienne_version"] = 42
    f2 = os.path.join(td, "vieux.json")
    json.dump(data, open(f2, "w"))
    _, b3, *_ = charger_configuration(f2)
    check("un champ inconnu est ignoré au lieu de faire échouer le chargement",
          b3.raquette_active)

# ---------------------------------------------------------------- J
titre("J. Optimisation")
for nom, extra in (("libre", {}), ("raquette", dict(raquette_active=True)),
                   ("raquette + oxyde", dict(raquette_active=True, oxydation_active=True))):
    b = params(**extra)
    res = evaluer_candidat(vecteur_depuis_params(b), b, AMPS, 40., 10., REG_E, REG_T)
    check(f"evaluer_candidat ({nom}) : 6 grandeurs finies",
          len(res) == 6 and all(np.isfinite(v) for v in res) and res[0] < 1e4,
          f"({[round(v, 3) for v in res]})")
    # Le DVH n'est pas indépendant de la dispersion : les deux se lisent sur la
    # même sinusoïde marche(θ) = A·cos(θ − φ). Le 6H et le 12H valent ±A·sin φ
    # et sont tous deux dans l'éventail, donc l'écart max entre positions
    # verticales vaut AU MOINS deux fois le DVH. C'est une inégalité de
    # structure : si elle tombe, c'est que les positions ne sont plus
    # diamétralement opposées deux à deux.
    check(f"dispersion >= 2 x |DVH| ({nom})", res[0] >= 2.0 * abs(res[5]) - 1e-9,
          f"(dispersion {res[0]:.4f} pour DVH {res[5]:.4f})")
    v = vecteur_depuis_params(b)
    xl, xu = construire_bornes(b, BornesOptimisation())
    check(f"vecteur et bornes cohérents ({nom})", len(v) == len(xl) == len(xu))
    check(f"aller-retour du vecteur exact ({nom})",
          np.allclose(v, vecteur_depuis_params(params_depuis_vecteur(v, b))))

b_court = params(n_elements_corps=250, grossmann_n_pas=30, phillips_n_pas=40)
so = sonder_faisabilite(b_court, BornesOptimisation(), AMPS, n_echantillons=8, seed=2,
                        pas_balayage_deg=60., reglage_echappement=REG_E)
check("sondage : géométrie de départ incluse et en tête",
      len(so["jeu_mm"]) == 9 and so["ordre"][0] == 0)
res_opt = lancer_optimisation(b_court, BornesOptimisation(),
                              HyperparametresGA(pop_size=8, n_gen=2, seed=1),
                              SPEC.jeu_min_spires, amplitudes_analyse_deg=AMPS,
                              pas_balayage_deg=60., reglage_echappement=REG_E,
                              reglage_thermique=replace(REG_T, oxyde_min_um=3., oxyde_max_um=4.))
check("optimisation avec les trois contraintes aboutit", res_opt is not None)

# ---------------------------------------------------------------- K
titre("K. Rapport et enregistrement de projet")
from spiral_lab.rapport import construire_rapport, enregistrer_projet
from docx import Document as _Doc

b_rap = params(n_elements_corps=250, grossmann_n_pas=30, phillips_n_pas=40,
               raquette_active=True, oxydation_active=True, raquette_jeu_um=5.0,
               raquette_transition_n_pas=30, raquette_arc_avant_n_pas=20,
               raquette_piton_n_pas=20)
g_rap = construire_geometrie(b_rap)

with tempfile.TemporaryDirectory() as td:
    etapes = []
    fichiers = enregistrer_projet(td, "essai", g_rap, SPEC, b_rap,
                                   bornes=BornesOptimisation(),
                                   hyperparametres=HyperparametresGA(),
                                   reglage_echappement=REG_E, reglage_thermique=REG_T,
                                   amplitudes_deg=AMPS, pas_balayage_deg=60.0,
                                   jeu_goupilles_um=b_rap.raquette_jeu_um,
                                   n_frames_gif=7, origine="suite de tests",
                                   progression=etapes.append)
    dossier = fichiers["dossier"]
    presents = sorted(os.listdir(dossier))
    check("dossier cree au nom du projet", os.path.basename(dossier) == "essai")
    check("cinq fichiers ecrits (lame oxydee : deux DXF)", len(presents) == 5, f"({presents})")
    check("rapport Word present", presents[0].endswith(".docx") and
          os.path.getsize(fichiers["rapport"]) > 200_000,
          f"({os.path.getsize(fichiers['rapport'])/1e6:.2f} Mo)")
    check("animation GIF autonome", os.path.exists(fichiers["gif"])
          and os.path.getsize(fichiers["gif"]) > 10_000)
    check("deux DXF quand la lame est oxydee", len(fichiers["dxf"]) == 2,
          f"({[os.path.basename(c) for c in fichiers['dxf']]})")
    check("configuration JSON presente", "configuration" in fichiers
          and os.path.exists(fichiers["configuration"]))
    check("progression rapportee", len(etapes) >= 6, f"({len(etapes)} etapes)")

    doc = _Doc(fichiers["rapport"])
    titres = [p.text for p in doc.paragraphs if p.style.name.startswith("Heading 1")]
    check("huit sections dans le rapport", len(titres) == 8, f"({len(titres)})")
    texte = "\n".join(p.text for p in doc.paragraphs)
    check("les quatre etats sont legendes",
          all(m in texte for m in ("Au repos", "Fermeture complète", "Ouverture complète",
                                    "Jeu inter-spires minimal")))
    check("le rapport renvoie vers le GIF anime", os.path.basename(fichiers["gif"]) in texte)
    check("images et animation embarquees",
          len(doc.inline_shapes) >= 9, f"({len(doc.inline_shapes)} images)")
    check("les quatre attitudes COSC sont presentes",
          "Attitudes COSC" in texte
          and all(f"Position COSC {n}" in texte for n in ("6H", "9H", "12H", "3H")))
    # Les resultats numeriques sont dans des TABLEAUX : doc.paragraphs ne les
    # voit pas, il faut parcourir les cellules.
    texte_tab = "\n".join(c.text for t in doc.tables for r in t.rows for c in r.cells)
    check("le DVH et le vecteur de sensibilite sont rapportes",
          "DVH (critère COSC)" in texte_tab and "Vecteur de sensibilité" in texte_tab)
    check("la ligne d'echappement est expliquee", "ligne droite de l'échappement" in texte)
    # Largeur utile lue sur le document lui-meme et non codee en dur : elle
    # depend des marges, qui different d'un document a l'autre.
    from docx.oxml.ns import qn as _qn
    sec = doc.sections[0]
    utile = round((sec.page_width.cm - sec.left_margin.cm - sec.right_margin.cm) / 2.54 * 1440)
    sommes = [sum(int(c.get(_qn("w:w")))
                  for c in t._tbl.find(_qn("w:tblGrid")).findall(_qn("w:gridCol")))
              for t in doc.tables]
    check("tableaux dans la largeur utile", all(x <= utile for x in sommes),
          f"(utile {utile}, sommes {set(sommes)})")

    # dossier deja rempli : refus sans autorisation explicite
    refuse = False
    try:
        enregistrer_projet(td, "essai", g_rap, SPEC, b_rap, amplitudes_deg=AMPS,
                           pas_balayage_deg=60.0, n_frames_gif=5)
    except FileExistsError:
        refuse = True
    check("refus d ecraser un projet existant", refuse)

    # spiral libre : un seul DXF
    f2 = enregistrer_projet(td, "libre", g_libre, SPEC, params(), amplitudes_deg=AMPS,
                            pas_balayage_deg=60.0, n_frames_gif=5)
    check("un seul DXF sans oxyde", len(f2["dxf"]) == 1,
          f"({[os.path.basename(c) for c in f2['dxf']]})")
    check("sept sections sans oxydation ni raquette",
          len([p for p in _Doc(f2["rapport"]).paragraphs
               if p.style.name.startswith("Heading 1")]) == 8)

    # nom comportant des caracteres interdits
    f3 = enregistrer_projet(td, 'proto/v2:essai?', g_libre, SPEC, params(),
                            amplitudes_deg=AMPS, pas_balayage_deg=60.0, n_frames_gif=5)
    check("caracteres interdits retires du nom",
          os.path.basename(f3["dossier"]) == "protov2essai",
          f"({os.path.basename(f3['dossier'])})")
    vide = False
    try:
        enregistrer_projet(td, "  ///  ", g_libre, SPEC, params())
    except ValueError:
        vide = True
    check("nom vide refuse", vide)


print()
print("RESULTAT : " + ("TOUT PASSE" if not ECHECS
                        else f"{len(ECHECS)} ECHEC(S) -- " + " | ".join(ECHECS)))
sys.exit(1 if ECHECS else 0)
