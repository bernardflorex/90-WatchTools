"""Sauvegarde et relecture d'une etude : parametres et geometrie.

Reprendre une etude demandait de re-saisir une vingtaine de champs puis de
re-importer un DXF. Un fichier JSON porte tout d'un coup : les parametres, le
contour du coeur, le ressort calcule et les reglages de l'optimisation.

Le format est volontairement plat et lisible, pour rester modifiable a la main
et diffable. Les angles sont en degres et les longueurs en millimetres, comme
dans l'interface.
"""

import json
import os

VERSION = 1

#: marqueur remplace par la liste de points compacte a l'ecriture
_POINTS = "@@points@@"


def _noms(params):
    """Noms de tous les parametres portes par cet objet."""
    specs = list(getattr(type(params), "SPECS", []))
    specs += list(getattr(type(params), "EXTRA", []))
    vus, out = set(), []
    for spec in specs:
        n = spec[0]
        if n not in vus and hasattr(params, n):
            vus.add(n)
            out.append(n)
    return out


def avec_json(path):
    """Ajoute l'extension si l'utilisateur ne l'a pas tapee."""
    return path if os.path.splitext(path)[1] else path + ".json"


def _points(profil, n=720):
    """Contour, d'ou qu'il vienne.

    Un profil importe porte ses points dans `pts` ; un profil optimise est une
    serie de coefficients et n'en a pas - il faut l'echantillonner. Sans ce
    detour, enregistrer une solution de l'optimisation levait
    `TypeError: 'Profil' object is not iterable`.
    """
    if profil is None:
        return None
    pts = getattr(profil, "pts", None)
    if pts:
        return list(pts)
    if hasattr(profil, "contour"):
        return list(profil.contour(n))
    try:
        return [(float(a), float(b)) for a, b in profil]
    except TypeError:
        return None


def sauver(path, params, profil=None, ressort=None, reglages=None,
           commentaire=""):
    """Ecrit l'etude. `profil` est un objet a `pts` ou une liste de points."""
    d = {"version": VERSION, "commentaire": commentaire,
         "parametres": {n: float(getattr(params, n)) for n in _noms(params)}}
    pts = _points(profil)
    if pts:
        d["coeur"] = {
            "repere": "axe de rotation du coeur a l'origine ; theta = 0 est "
                      "la position de repos, plat du coeur sur plan du "
                      "marteau",
            "points": _POINTS}
    if ressort is not None:
        d["ressort_marteau"] = {
            "prearmage_Nmm": float(ressort.prearmage),
            "raideur_Nmm_par_rad": float(ressort.raideur),
            "course_deg": float(ressort.course) * 180.0 / 3.141592653589793,
            "residu_pourcent": float(ressort.ecart)}
    if reglages:
        d["optimisation"] = {k: (int(v) if isinstance(v, bool) else v)
                             for k, v in reglages.items()}
    txt = json.dumps(d, indent=2, ensure_ascii=True)
    if pts:
        # Un tableau de 720 couples indente fait 2880 lignes : on le remet sur
        # une ligne par point, pour que le fichier reste relisable a l'oeil.
        compact = ",\n      ".join("[%.6f, %.6f]" % (a, b) for a, b in pts)
        txt = txt.replace('"%s"' % _POINTS,
                          "[\n      %s\n    ]" % compact)
    path = avec_json(path)
    with open(path, "w") as fh:
        fh.write(txt + "\n")
    return path


def charger(path, params=None):
    """Relit l'etude. Renvoie `(parametres, points, reglages, journal)`.

    `params`, s'il est fourni, est mis a jour en place : seuls les champs
    qu'il porte reellement sont repris, les autres sont signales au journal.
    """
    with open(path) as fh:
        d = json.load(fh)
    log = []
    if not isinstance(d, dict):
        raise ValueError("ce fichier ne contient pas une etude")
    v = d.get("version")
    if v is not None and v > VERSION:
        log.append("fichier en version %s, lu par une version %d : "
                   "des champs peuvent etre ignores" % (v, VERSION))
    vals = {}
    inconnus = []
    for n, x in (d.get("parametres") or {}).items():
        if params is None or hasattr(params, n):
            vals[n] = float(x)
        else:
            inconnus.append(n)
    if params is not None:
        for n, x in vals.items():
            setattr(params, n, int(round(x))
                    if n in ("n_lobes", "n_steps") else x)
    if inconnus:
        log.append("parametres ignores (inconnus de cette version) : "
                   + ", ".join(sorted(inconnus)))
    pts = None
    coeur = d.get("coeur") or {}
    if coeur.get("points"):
        pts = [(float(a), float(b)) for a, b in coeur["points"]]
        log.append("contour du coeur : %d points" % len(pts))
    if d.get("commentaire"):
        log.append("commentaire : " + str(d["commentaire"]))
    log.append("%d parametres repris" % len(vals))
    return vals, pts, (d.get("optimisation") or {}), log
