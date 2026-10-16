#!/usr/bin/env python3
"""Dimensionnement de cames horlogeres.

Le Python du systeme fait tourner le programme en Python pur. Pour le noyau
compile (numba) et l'interface, utiliser l'environnement du projet :

    /home/florian/DEV/virtual_env/.watchtools/bin/python main.py

depuis ce dossier, ou `python -m came.main` depuis la racine WatchTools.

Usage :
    python3 main.py                 # interface graphique
    python3 main.py --cli           # synthese en console
    python3 main.py --cli --export coeur.dxf
"""

import argparse
import math
import os
import sys

if __package__ in (None, ""):
    # lancé comme script (python main.py depuis ce dossier, bouton « Run » de
    # VSCode) : le répertoire parent, racine de WatchTools, porte `commun` et
    # ce paquet lui-même
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def cli_chrono(args):
    """Synthese du coeur de chronographe, et evaluation d'un marteau a 2 pans."""
    from came.chrono import ChronoParams, synthesize_chrono
    from came.dxfio import write_dxf_chrono
    p = ChronoParams()
    for name, _l, _u, _a, _b, _s in ChronoParams.SPECS + ChronoParams.EXTRA:
        v = getattr(args, name, None)
        if v is not None:
            setattr(p, name, v)
    r = synthesize_chrono(p)
    if not r.ok:
        print("echec :", r.message)
        return 1
    print("courses totales      : %.2f / %.2f deg (somme %.2f)"
          % (r.total_a, r.total_b, r.total_a + r.total_b))
    print("couple sur le coeur  : %.5f N.mm" % r.couple_coeur)
    print("couple marteau C0    : %.4f N.mm" % (p.couple_min * r.echelle))
    print("effort normal maxi   : %.4f N" % r.N_max)
    print("course du marteau    : %.3f deg" % math.degrees(r.psi_max))
    print("plat du zero         : %.4f mm (mini mu.r_min = %.4f)"
          % (r.plat, r.plat_mini))
    print("plage morte au sommet: %.2f deg (plancher 2.atan(mu) = %.2f)"
          % (r.zone_morte, math.degrees(2.0 * math.atan(p.mu))))
    for w in r.warnings:
        print("  !", w)
    if args.nsga:
        import time
        from came.nsga2 import nsga2, pool_lot, probleme
        from came.profil import (evaluer, from_chrono, longueur_plat,
                                 Profil)
        print()
        print("optimisation NSGA-II du profil discretise")
        _pr0 = from_chrono(p, r, args.cellules)
        graine = _pr0.vector()
        # le plat de la graine est la reference : l'optimisation ne doit pas
        # le dissoudre, sous peine de perdre la position nulle du coeur
        plat_ref = 0.9 * longueur_plat(_pr0, p.gamma0, 0.002, grille=300)
        if args.galet:
            args.face = True
        if args.face:
            graine = graine + [0.0]
            print("  plat du marteau libre (courbure optimisee)")
        from came.nsga2 import MARTEAU_LIBRE, MARTEAU_UTILE
        choix = args.marteau.strip()
        mart = (MARTEAU_UTILE if choix in ("tout", "utile")
                else tuple(n.strip() for n in choix.split(",")
                           if n.strip() in MARTEAU_LIBRE))
        if mart:
            # la graine porte les valeurs de depart du marteau en queue ;
            # pivot_t est un deplacement, donc nul au depart
            graine = graine + [0.0 if n == "pivot_t"
                               else float(getattr(p, n)) for n in mart]
            print("  grandeurs du marteau optimisees : " + ", ".join(mart))
        ressort = None
        if args.figer:
            from came.profil import dimensionner_ressort
            r0 = dimensionner_ressort(p, from_chrono(p, r, args.cellules),
                                      0.0, None, n=args.scan)
            ressort = (r0.prearmage, r0.raideur)
            print("  ressort fige sur la version initiale : C0 = %.4f N.mm, "
                  "k = %.4f N.mm/rad" % ressort)
        ev, bornes, rep, gfaces = probleme(p, n_scan=args.scan,
                                           face_libre=args.face,
                                           galet=args.galet, ressort=ressort,
                                           plat_mini=plat_ref,
                                           gabarit=_pr0.gabarit,
                                           marteau=mart)
        graines = gfaces(graine)
        lo, hi = bornes(graine)
        f0, v0 = ev(graine)
        print("  rayon de raccord usinable impose : %.3f mm"
              % p.rayon_raccord)
        print("  donnees : encombrement %.2f a %.2f mm, pivot (%.2f, %.2f),"
              % (p.r_min, p.r_max, p.pivot_x, p.pivot_y))
        print("            pente du plan %.1f deg, couple vise sur le coeur "
              "%.4f N.mm" % (p.angle_plan, p.couple_coeur))
        print("  sortie  : geometrie du coeur + prearmage et raideur du "
              "ressort")
        print("  graine analytique : %.1f %% du tour sous la consigne, "
              "plage %.2f deg, ondulation %.1f %%, N %.3f N" % f0)
        lot, fermer = pool_lot(p, args.scan, args.face, args.procs,
                               args.galet, ressort, plat_ref,
                               _pr0.gabarit, mart)
        if lot:
            import multiprocessing as _mp
            print("  evaluation repartie sur %d processus"
                  % (args.procs or max(1, (_mp.cpu_count() or 2) - 1)))
        t0 = time.time()
        _pop, front = nsga2(
            ev, lo, hi, n_pop=args.population, n_gen=args.generations,
            graines=graines, reparer_x=rep, evaluer_lot=lot,
            progress=lambda g, n, nf, _p: print(
                "  generation %d/%d : %d faisables" % (g, n, nf))
            if g % 5 == 0 else None)
        fermer()
        fais = [q for q in front if q["viol"] <= 1e-9] or front
        print("  %.0f s, front de %d solutions (%d faisables)"
              % (time.time() - t0, len(front), len(fais)))
        def rface(x):
            if not args.face:
                return 0.0
            k = x[-1]
            return 0.0 if abs(k) < 1.0 / (20.0 * p.r_max) else 1.0 / k

        print("  %10s %10s %9s %9s %10s %10s %9s"
              % ("sous T [%]", "plage [deg]", "ondul[%]", "N [N]",
                 "C0 [N.mm]", "k [N.mm/r]", "R face"))
        vus = set()
        from came.profil import Profil, dimensionner_ressort, mu_galet
        nq = (3 if args.face else 2) + len(mart)

        def ressort_de(x):
            if ressort is not None:
                cl = dimensionner_ressort.__class__
                class _R:
                    prearmage, raideur = ressort
                return _R()
            # le gabarit : sans lui le ressort serait dimensionne pour un
            # coeur sans plat, donc pour une autre piece
            pr = Profil.from_vector(list(x[:-nq]) + [x[-nq], x[-nq + 1]],
                                    *_pr0.gabarit)
            R = rface(x)
            mu = (mu_galet(R, p.axe_galet, p.mu_palier)
                  if args.galet and abs(R) > 1e-9 else None)
            return dimensionner_ressort(p, pr, R, mu, n=args.scan)

        for q in sorted(fais, key=lambda z: (z["f"][0], z["f"][1])):
            cle = (round(q["f"][0], 1), round(q["f"][1], 2),
                   round(q["f"][2], 3), round(q["f"][3], 2))
            if cle in vus:
                continue
            vus.add(cle)
            R = rface(q["x"])
            f = q["f"]
            re_ = ressort_de(q["x"])
            print("  %10.1f %10.2f %9.1f %9.4f %10.3f %10.3f %9s"
                  % (f[0], f[1], f[2], f[3], re_.prearmage, re_.raideur,
                     "plane" if R == 0 else "%.2f" % R))
        return 0
    if args.export:
        write_dxf_chrono(args.export, p, r)
        print("DXF ecrit :", args.export)
    return 0


def cli(args):
    from came.core import Params, norm, synthesize
    from came.dxfio import write_dxf
    p = Params()
    for name, _l, _u, _a, _b, _s in Params.SPECS:
        v = getattr(args, name, None)
        if v is not None:
            setattr(p, name, v)
    r = synthesize(p, mode=args.mode)
    if not r.ok:
        print("echec :", r.message)
        return 1
    print("couple coeur constant     : %.5f N.mm" % abs(r.T0))
    print("couple marteau            : %.4f -> %.4f N.mm" % (min(r.Cm), max(r.Cm)))
    print("augmentation du couple    : %.4f N.mm" % r.delta_Cm)
    print("course du marteau         : %.3f deg" % math.degrees(r.psi_max))
    print("effort normal maxi        : %.4f N" % r.N_max)
    print("rayon du profil           : %.4f -> %.4f mm"
          % (min(norm(q) for q in r.profile), r.r_eff))
    print("rendement du contact      : %.3f" % r.rendement)
    print("ecart rentrant/sortant    : %.1f %%" % (100.0 * r.ripple))
    print("valide                    :", r.valid)
    for w in r.warnings:
        print("  !", w)
    if r.message:
        print("  >", r.message)
    if args.export:
        write_dxf(args.export, p, r)
        print("DXF ecrit :", args.export)
    return 0


def main(argv=None):
    from came.core import Params
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cli", action="store_true", help="calcul en console")
    ap.add_argument("--chrono", action="store_true",
                    help="coeur de chronographe (implique --cli)")
    ap.add_argument("--nsga", action="store_true",
                    help="optimise le profil discretise par NSGA-II "
                         "(implique --chrono)")
    ap.add_argument("--cellules", type=int, default=24,
                    help="nombre d'harmoniques du profil")
    ap.add_argument("--face", action="store_true",
                    help="optimise aussi la courbure du plat du marteau")
    ap.add_argument("--marteau", default="",
                    help="grandeurs du marteau a optimiser, separees par des "
                         "virgules : longueur_plan, decalage_plan, "
                         "rayon_coin, angle_plan, pivot_x, pivot_y. "
                         "'tout' les prend toutes.")
    ap.add_argument("--figer", action="store_true",
                    help="fige C0 et k sur la version initiale et optimise "
                         "coeur et marteau a armement identique")
    ap.add_argument("--galet", action="store_true",
                    help="galet monte fou au bout du marteau (implique --face)")
    ap.add_argument("--population", type=int, default=60)
    ap.add_argument("--generations", type=int, default=40)
    ap.add_argument("--scan", type=int, default=240,
                    help="pas d'analyse sur un tour")
    ap.add_argument("--procs", type=int, default=None,
                    help="nombre de processus d'evaluation (defaut : coeurs-1)")
    ap.add_argument("--mode", default="retour", choices=["retour", "armement"])
    ap.add_argument("--export", help="fichier DXF de sortie")
    from came.chrono import ChronoParams
    for name, label, unit, _a, _b, _s in Params.SPECS + ChronoParams.EXTRA:
        ap.add_argument("--" + name.replace("_", "-"), dest=name, type=float,
                        help="%s [%s]" % (label, unit))
    args = ap.parse_args(argv)
    # --nsga n'est implemente que dans la synthese de chronographe : sans
    # cela il tombait sur le lancement de l'interface, sans rien dire.
    # La came simple est remise a plus tard : --cli vaut --chrono.
    if args.chrono or args.nsga or args.cli:
        return cli_chrono(args)
    from came.ui import run
    return run()


if __name__ == "__main__":
    sys.exit(main())
