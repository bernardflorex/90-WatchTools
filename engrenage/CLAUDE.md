# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Python port, in progress, of ENGCRE, a Fortran program by K. Trümpy that generates gear-tooth profiles (involute or cycloid) and analyses the meshing of a gear pair (kinematics, statics, Hertz pressure). The Fortran version works and is the reference for all the maths. Users are French-speaking; Fortran comments are in German.

Ported so far: [evolvente.py](evolvente.py) (port of `Fortran/Egc/Evolve.F90`, file reading/writing, drawing of the gears), [contact.py](contact.py) (numpy port of the contact search of `Kontap.F90`/`Enganap.F90` and of the transmission/torque/efficiency formulas of `Enganap.F90`/`Draweng.F90`) [arcs.py](arcs.py) (port of `Krebax.F90`/`Lagrmo.F90`: profile approximated by tangent circular arcs; reading and writing of arcs files in `.json`, `.Eng` and `.Eno`) [optim.py](optim.py) (numba port of `KRBGS`/`KORZF`/`SWTE`, `Kontac`, `Enganac` and `Resulto` of `EngOpt/`, objective function and constraints), [recherche.py](recherche.py) (CMA-ES driver with restarts, repair by the back half-profile, list of distinct optima), [fraise.py](fraise.py) (port of `Abwfrac.F90`/`Abwfrap.F90`/`Radabwp.F90`/`Krvabc.F90`/`Krvabp.F90`: hobbing-cutter profile conjugate to a wheel, simulated hobbing of a target wheel, deviation from the reference profile) and [main.py](main.py), a PySide6 window: four tabs on the left ("Création": parameters of the involute gears, arcs and export; "Simulation arcs": meshing simulation of circular-arc profiles, from the Création tab or from an opened file; "Optimisation": parameters, constraints, run and list of optima; "Fraise": hobbing-cutter profile, `.Frs` file, simulated cutting), the matplotlib drawing on the right. Everything else (`Cyclo`, the circular-arc contact analysis `Enganac`, Hertz pressure…) is still Fortran-only. `rapport.py` writes a `.docx` simulation report (File > "Exporter un rapport…" in the GUI) from `contact.py`'s results, independent of the optimiser/hobbing modules. `dxf.py` (using `ezdxf`) exports mobile and hobbed-wheel profiles as real arc/line DXF entities and the hob profile as a DXF polyline.

## Commands

```
python main.py                           # from the WatchTools root: the launcher (one icon per tool, each started as its own process)
python main.py engrenage                 # from the WatchTools root: this tool, without the launcher
python -m engrenage.main [file.json|file.Cev]   # from the WatchTools root: this tool alone (default file: test.Cev, imported)
python main.py [file]                    # from THIS folder (or VSCode "Run" / launch config "Engrenage"): this tool alone
python -m engrenage.optim test.Cev 4     # any other module: from the root, or from this folder in the VSCode terminal (PYTHONPATH=..)
python engrenage/evolvente.py [file] [--ni 150] [--plot] [--png out.png] [--dents 3]   # CLI, also writes Result/<name>.Egr
```

- Since the merge into WatchTools this directory is the package `engrenage`: modules import each other with `from . import arcs` etc., so they run as `python -m engrenage.<module>` with the WatchTools root on the path; only `evolvente.py` (no internal import) also runs as a plain script. `main.py` is the exception: it imports `from engrenage import ...` (absolute) and, when run as a script, puts the WatchTools root on `sys.path` itself, so `python main.py` works from this folder. The user opens one VSCode window (and one Claude session) per tool folder to work on several tools in parallel: see the root `CLAUDE.md` for the shared rules (one git working tree for all tools). The Qt theme and matplotlib colours come from `commun/theme.py` (dark, the cam tool's): draw with `theme.encre()` / `theme.gris(n)` instead of `"k"` or a fixed grey, and report figures stay white through `@theme.clair()`. The canvas has the cam view's mouse navigation (`commun/navigation.py`: wheel = zoom on the cursor, drag = pan, double-click = back to the overall view); it steps aside while a toolbar tool is active.

- Interpreter: the system `python3` has neither pip nor matplotlib/PySide6. Use the venv `~/DEV/virtual_env/.watchtools/bin/python` (matplotlib, PySide6, numpy, scipy, numba, cma…; no tkinter, so the matplotlib GUI backend is Qt). Headless GUI test: `QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg`, then `window.grab().save("x.png")`.
- Careful: the CLI overwrites `Result/<name>.Egr`, and `Result/f.Egr` is the tracked **Fortran** reference (restore it with `git checkout -- Result/f.Egr`). Run tests from a scratch directory or pass a file whose name is not `f`.
- `evolvente.py` itself is pure standard library, matplotlib is imported only for drawing; `contact.py` needs numpy. `dessiner_engrenage(ax, …)` draws `--dents` teeth of each wheel almost in mesh (driven wheel turned by 180° minus half a pitch, layout as in `Drawkpo.F90`); the profiles are half-teeth, so an odd tooth count is not symmetric for the driven wheel.
- The Qt font is the shared theme's default (the former +2 pt `agrandir_police` was dropped at the merge); long unwrapped labels or checkbox texts widen the left panel and cause horizontal scrolling, so keep texts short or `_etiq` (word-wrapped). `main.py` state is split per tab: Création (`roues`, `profils`, `arcs_roues`) and Simulation (`mobiles` = two `arcs.Mobile`, `sim`, `sim_roues`, `sim_profils`), both drawing on the same axes (`_vue_ok` keeps the user's zoom between redraws; it is reset on tab change). F5 is one menu action (`actualiser`) that rebuilds or simulates depending on the tab, so the buttons carry no shortcut (two identical shortcuts would be ambiguous). The Simulation tab runs `contact.Engrenement(..., mobiles=self.mobiles)`, the analytic arc-to-arc engine (no profile points involved; the sampled arcs are only passed for the centre-distance test and for drawing). `File > Open` dispatches by extension (and by content for `.json`: `_json_est_arcs`). `main.py` builds the left panel from the `CHAMPS` table (one row per `.Cev` parameter; a new user variable = a new `Champ` plus an attribute on `evolvente.Roue`). Edits do not recompute until the "Reconstruire" button (F5), which runs the whole simulation (about 0.5 s). A slider walks through the 150 positions of one tooth-pair passage; the drawing is redrawn per position (fixed margins, `constrained` layout is too slow). The three curves (W transmission, C torque, N efficiency, as in the Fortran layout `layout de sortie des simulations.png`) are drawn over the gears at the x of the contact point, in % of the ideal value, 100 % at the pitch point (`evolvente.y_pourcent`: 5 % = 0.6 module); the panel also shows min/mean/max/variation per parameter. Files are saved as `.json` only; the old `.Cev` format is import-only (`lire_cev`), so the Fortran program cannot read what the GUI saves. The Hertz data (E-modulus, Poisson, thickness) are kept in the JSON (`hertz` key per wheel) but not editable yet.
- Bug found by the user: closing the app (or opening/reloading a file) always asked to save, even with no edits. Cause: `remplir_optimisation()` set the two objective combo boxes with `setCurrentIndex(...)` without `blockSignals`, so their `currentIndexChanged` (wired to `_marquer_modifie`) fired on every programmatic fill (startup, every load) and set `self.modifie = True` before the user touched anything. `case_dissym` already did this correctly; the fix wraps the two combo updates the same way. `QLineEdit` fields were never affected (`textEdited` only fires on user input, unlike `setText`).
- There is no test suite or build system. The git remote is `https://github.com/bernardflorex/eng_cre.git`; `Fortran/` is in `.gitignore`. To check a change to `evolve()`, compare with the Fortran output `Result/f.Egr`, which was generated from `Files/f.Cev` with **NI=550**: `evolve(roue, 550)` must match it to about 1e-7.
- The Fortran cannot be built or run here: it needs Windows and Intel Fortran (`ifort`, QuickWin graphics, `winapi_c.cpp`). `Fortran/*.Exe` and `libiomp5md.dll` are Windows binaries.
- `.vscode/` (this folder opened on its own in VSCode): `settings.json` sets the `.watchtools` interpreter, `python.analysis.extraPaths: [".."]` (so Pylance resolves `commun`) and `PYTHONPATH=..` for the integrated terminal; `launch.json` has "Engrenage" (`-m engrenage.main`) and "module courant" (`-m engrenage.<current file>`).

## Fortran architecture (the spec to port)

Source is in `Fortran/Egc/`. `Fortran/Conc.F90` is a generated concatenation of all sources (`Conc.Bat`), so ignore it. The `.Bat` files are the Windows build scripts. Files are Latin-1 (cp1252): read them with `iconv -f latin1 -t utf8 <file>`.

Pipeline of the `ENGCRE` program (`Engcre.F90`, an interactive console dialog in D/F/E languages):

1. `Files\<name>.Cev` (involute) or `.Cyl` (cycloid) gives both wheels. The suffix `,I` after the file name swaps driver and driven.
2. `EVOLVE` / `CYCLO` produce a point sequence `ZP(2,NI,wheel)` per wheel (NI=150 by default).
3. Option `Z`: `KREBAX` approximates `ZP` by circular arcs (uses `LAGRMO`, a least-squares fit with Lagrange multipliers), giving `XY`/`PM`/`RA`.
4. Meshing analysis over `IK` driver positions, filling the array `KW(18,IK)` (angles, contact point, sliding direction, transmission ratio, torque, forces, Hertz pressure):
   - `ENGANAP` for point-sequence profiles, using `KONTAP` (contact of polylines) and `SPLIPH` (Hermite spline);
   - `ENGANAC` for circular-arc profiles, using `KONTA`.
5. Option `P`: writes the unformatted binary `TRANSF.SCR` and runs `PLOEGC.Exe` (`Ploegc`, `Winpl`, `Draweng`/`Drawkpo`/`Drawpro`, `Pfeil`).
6. Option `Q`: writes `Result\<name>.Egr` (profiles + Hertz block) and, if `Z` was used, `.Eno`.

`Consini`, `Eapos`, `Datwr` are Windows-console helpers. `Prolic` is a MAC-address licence check, disabled in `Engcre.F90` (`LC='Y'`), so do not port it. `EGCStart.Bat` holds start parameters read by position: `PZ` (tooth pairs checked), `NJ` (which sets the sizes below) and the language.

`NI` and `IK` come from `NJ`, as `IK=101+NJ*100` and `NI=150+NJ*100`. `EGCStart.Bat` has `NJ=0`, hence 150 profile points.

## Data formats

`.json` (current format, written by the GUI): `{"format": "engcre", "version": 1, "roues": {"menante": {...}, "menee": {...}}}`, each wheel holding `mo zz nw qe qg ma me` and optionally `"hertz": {"e_mod", "poisson", "epaisseur"}`, plus `"narcs"` per wheel (arcs on the involute for the export) and `"simulation": {"mu": 0.2, "mt": 1.0, "da": 0}` (friction coefficient, driver torque in mN·m = N·mm, centre-distance variation in µm).

`.Cev` (legacy Fortran format, import only; involute, whitespace separated, first line is a header, CRLF possible):
```
line 2: MO ZZ NW QE QG MA ME      # driver wheel
line 3: MO ZZ NW QE QG MA ME      # driven wheel
then (optional for profile generation): E-mod / Poisson / thickness header + 2 lines for Hertz
```
`MO` module, `ZZ` tooth count, `NW` pressure angle in degrees, `QE` gap/tooth-width ratio at the pitch circle (≥ 1), `QG` length of the straight root transition in module units (0 = none), `MA`/`ME` extent of the involute below/above the pitch circle in module units (some files call them M1/M2). `.Cyl` replaces `NW` with `R1 R2` (rolling-circle radii).

Arcs files (read by the Simulation tab, written by `Export Arcs`): (a) `.json` with the schema of `profile.json`: a list of `{"id","type","details"}` for `Engrenage` (`dA_pos dA_neg h_contact`), `Mobiles` (`mobile_vec`: per wheel `is_master nom z m sym profile[[x,y,flag]] arc_rayon arc_centre epaisseur E poisson ty`) and `Simulation` (`opti_1 opti_2 ponderation centrage`, kept for later use); (b) `.Eng` / older `.json`: `{"m": {...}, "s": {...}, "sim": {...}}` with the same wheel keys; (c) `.Eno`: the Fortran text file. The arcs are always rebuilt from the transition points and the start slope `-pi/Z` (`krbgs`), never from the stored `arc_rayon`: for `Bar-GDM.Eng` the stored radii and centres agree with `krbgs` to 7e-5, but the radii of the example `profile.json` do not match its own points. Asymmetric teeth: `sym: false` in a wheel of the arcs `.json` plus `profile_arriere` (+ `arc_rayon_arriere`, `arc_centre_arriere`), the second half-profile in the same convention as `profile` (`arcs.CLE_DOS`; `arcs.Mobile.points_dos`, None = symmetric); `profile` is the driving half-profile (mirrored, x<0, in the drawing), `profile_arriere` the other flank (x>0). `.Eno` cannot hold it (`ecrire_eno` refuses). In the JSON the profile flags are 0 at the ends and 1 inside (Fortran `.Eno` convention); `arc_rayon` sign: positive when the centre is to the left of the chord P(i)->P(i+1).

`.Egr`: per wheel, 3 header lines (label and parameters; `MO ZZ 1 0`; `NI`) followed by `NI` lines `x y`, then a Hertz block. `evolvente.py` writes the same layout without the Hertz block and with a simpler label.

## Contact and meshing simulation (`contact.py`)

Drawing frame: pitch point at (0,0), driver (wheel 1) centre at (0,-rt1) turning counter-clockwise, driven (wheel 2) centre at (0,rt2), angles counter-clockwise from +y. The driving flank of both wheels is the mirrored profile (`x<0`).

- `premier_contact` (KONTAP): rotates the driven flank from a far reference position until it touches the driver flank; returns the rotation, the contact point and the common tangent. It works on polylines, so the contact point is only accurate to about the square root of the chord sagitta: `Engrenement._contact` first searches with `NA=100` points, then repeats the search on `2*FEN` fine points (`NF=3000`) around the contact. Without this refinement the transmission ratio is noisy at the 1e-3 level.
- Two engines share the same search (`debut`, `simuler`, `jeu`, `efforts`): the polyline one (profile points, above) and the arc one (`Engrenement(..., mobiles=...)`, port of `Konta.F90`/`KONTW`), which is what the Simulation tab uses. `premier_contact_arcs` takes two flanks made of circular arcs (`flanc_arcs`: centres, radii, convexity = positive signed radius, start angle and sweep) and looks for the smallest rotation of the driven wheel that makes an arc of it tangent to an arc of the driver flank: centre distance `r1+r2` (both convex) or `r_concave-r_convex` (one convex inside a concave one), never two concave arcs; the centre of the driven arc runs on a circle around its axis, so the law of cosines gives the rotation; the tangent point lies on the line of centres and must belong to both arcs; the common tangent is perpendicular to it. Cross-check: same driver angles on the arcs sampled with 1500 points and the polyline engine give the same driven angle to 5e-8 rad, contact points to 0.05 µm and tangents to 0.005°. In arc mode the ties are exact (`tol=1e-9`) and the binding pair is the strict minimum, as in the Fortran (`seuil_paire=0`): with few arcs the two pairs alternate (the arcs are not exactly conjugate) and the line of contact and the curves show the jumps (`evolvente._tracer_par_paires` draws them as thin dotted connectors, like `Draweng.F90`). Arc results tend to the 150-point ones when the number of arcs grows (`test.Cev`, 9 arcs: backlash 2.668-2.672 vs 2.666-2.671 degrees, mean efficiency 96.23 vs 96.22 %); with 3 arcs the transmission ratio spans 98.4-101.4 % instead of 99.8-100.8 %. About 0.07 s per simulation.
- Hertz pressure (arc mode only): the contact returns the signed radii of curvature of the two contacting arcs (convex positive); `Engrenement._hertz` reproduces `KW(17)` of `Enganac.F90`: `p = sqrt(F*E*/(pi*H*(1-nu^2)*R))`, `E*=E1E2/(E1+E2)`, mean `nu`, `H` = smaller thickness, `R=|R1R2/(R1+R2)|`, with `F` the normal force for a driver torque of 1 mN.m (always 1, independent of the `Cm` field). `Simulation.hertz_local` is the maximum-pressure contact point in each wheel's own frame (drawn as a star on every tooth). Files without Hertz data get `HERTZ_DEFAUT` from `main.py`.
- Backlash (`Engrenement.jeu`): `theta_max(phi)+theta_max(-phi)` folded modulo one driven pitch; in arc mode it is evaluated at all 150 positions (`simuler(..., n_jeu=...)` lowers that for optimisation loops), with the polyline engine at 10 positions only.
- `Simulation.criteres()` gives the four optimisation criteria: torque variation (%, minimise), maximum backlash (deg, minimise, with a minimum admissible backlash as constraint), mean efficiency (%, maximise), maximum Hertz pressure (N/mm^2, minimise). The Optimisation tab only holds the parameters (`Projet.optimisation`, also written to `Simulation` of the arcs JSON as `opti_1`/`opti_2`/`ponderation`) and an evaluation/timing button; it drives `optim.py` / `recherche.py` (see the section below). A full simulation costs about 0.06 s (0.19 s with the backlash at 150 positions).
- Sensitivity warning (measured on `test.Cev`, 4 arcs): with the chained tangent construction, moving ONE transition point (x, y) by 1 micron rotates the following tangents by 3 to 12 degrees, so with (x, y) variables the acceptable neighbourhood is only 1-3 microns (`optim.py` uses other variables, see below).
- `Engrenement.debut` finds the driver angle where the next tooth pair takes over from the previous one (FW in `Enganap.F90`), i.e. where two pairs touch at once. Pair `k` (driver tooth `k`, driven tooth `-k`) puts the driven axis at `psi(alpha1 + k*A1) + k*A2`; the driven wheel sits on the smallest value. Both pairs must really touch, otherwise `ErreurContact` (contact ratio < 1, e.g. `Files/Evolve3.Cev` with MA=0).
- `Engrenement.simuler` walks `N_POS=150` driver positions from that angle over one pitch `A1`, and returns arrays for the contact point (the line of contact), the second contact during the two-pair overlap, and the three quantities of `Draweng.F90`: transmission ratio `w2/w1` (negative: opposite senses), torque ratio `M2/M1` and efficiency, their product. Forces come from Coulomb friction `mu` opposing the sliding velocity, the driver torque being 1 N·mm (`efforts`); `Simulation.normale` is the normal force for that torque, scaled by the user's `mt` in the GUI. Pair 0 is kept in the overlap unless pair 1 is lower by more than `SEUIL_PAIRE`.
- `evolvente._placement`'s wheel centres use the working (`rayons_fonctionnement`) pitch radii, which always sum to the real centre distance by construction and would draw the dashed "cercle primitif" overlay as forever tangent regardless of `da` (a real bug, reported by the user: a squeeze should make the pitch circles overlap, an opening should separate them). Fixed by drawing that overlay with the fixed nominal radii `m*z/2` instead (unaffected by `da`) at the same, correctly placed, centres: their sum is the nominal centre distance, so a squeeze (negative `da`) makes them overlap by `|da|` and an opening separates them by `da`, checked exactly on `test.Cev` (30 µm each way). The working radii stay unchanged for placing the axes and for `contact.py`'s kinematics, which do need them always tangent.
- Centre distance: `Engrenement(roues, profils, da)` (`da` in mm) places the wheels at `rt1 + rt2 + da` with working pitch circles proportional to the tooth counts, tangent at (0,0); `_tester_entraxe` ports the FI=1/FI=2 tests of `Enganap.F90` (tip against root / tip circles no longer crossing). `Engrenement.jeu` gives the driven-wheel backlash: by symmetry about the line of centres, the lowest driven position with the driver at `phi` is the mirror of the highest one at `-phi`, so backlash = `theta_max(phi) + theta_max(-phi)` folded modulo one driven pitch. It agrees with a brute-force outline-crossing search to about 5e-3 degrees, and changes by the expected `2*dA*tan(alpha)/r2` (30 µm gives +1.57° on `test.Cev`).
- Checks done by hand: with `mu=0` the efficiency is exactly 1; the line of contact is straight at the pressure angle; contact point and tangent are consistent with both base circles to about 1e-8 (driver) and 1e-5 (driven); efficiency drops with the distance to the pitch point as the first-order formula `1 - mu*(1+Z1/Z2)*|x|/rb1` predicts (within 1%); on 6 files of `Files/` the full tooth outlines do not cross at the start position. Not compared with the Fortran numbers here (the Fortran cannot be run in this environment) -- but the **user has compared the Python simulation with the Fortran ENGCRE program himself** and judges it "tout à fait cohérente"; that is the only independent confirmation of `contact.Engrenement` (no Fortran figures were passed on).

## Optimisation (`optim.py`, `recherche.py`, `EngOpt/`)

`EngOpt/` holds the Fortran sources of the optimiser ENGOPT (`Optimal`, `Bidiopt`, `Variopt`, `Resulto`, `Enganac`, `Kontac`, `Krbgs`, `Manipu`; not run here). What it does: the objective `Resulto` is the square root of the sum, over the centre distances -dA, 0, +dA (line `-dA +dA GW ZV GL` of the `.Eno` file), of `(100 varC)^2` (ZV=1) or `(100 (1-mean efficiency))^2` (ZV=2), plus `(Hertz/GW)^2` and `((1-SymL) GL)^2`; ZV=3 minimises the backlash with minimum `|GW|`. A profile is invalid (1D33) when `SWTE` fails (exactly one sign change of the arc radii, tip radius >= `RK` x module), when the centre distance is inconsistent or when there is no mesh period. The variables are the interior transition points (x, y) with flag `IL=1`; the tip point is recomputed by `KORZF` (last arc centred on the tooth axis). `BIDIOPT`: per point a (2MN+1)^2 grid, weighted sum of the best moves, step reduced by 0.618, multi-start only through different initial steps; `VARIOPT`: full factorial (2MN+1)^(NP-2) per axis and wheel.

`optim.py` is the fast Python counterpart, checked against `contact.Engrenement` (`python optim.py [file] [narcs]`: same torque variation, Hertz pressure and backlash to 3e-5; 1000-3000 evaluations/s on 16 cores). Numba compiles on first use (a few seconds, cached in `__pycache__`); it and `cma` are in the `.eng` venv.
- `Probleme(roues, mobiles, mu, da_neg, da_pos, libres, objectif, dos_libre, avant_libre)`: decodes a vector of variables into arcs (`decoder`, `points`, `mobiles`), evaluates N vectors at once (`evaluer` -> array (N, 3 centre distances, `NCOL` columns `VARC…ECART_TETE`), `valeur`, `Probleme(x)`), reports one profile in detail (`rapport`). Variables = curvature and length of each free arc of the first n-2 arcs (`ECHELLE`; first point fixed, last arc closed like `KORZF`); `dos_libre` adds the same variables for the back half-profile (relative to the driving one), `avant_libre=False` keeps the driving profile and moves only the back one (repair).
- `Objectif`: `J = f1/f1₀ + poids f2/f2₀` (`couple`: torque variation, `jeu`: max backlash at ±dA, `rendement`: efficiency loss, `hertz`: max Hertz; quadratic mean over the 3 centre distances; `f₀` = start profile, so J = 1 at the start) + penalties: backlash >= `jeu_min` at -dA (anti-blocking) and <= `jeu_max` at +dA, tip-arc radius >= `r_tete_min` (mm, hob life; the Fortran `RK` is a relative test on x of the tip-arc start, this one is the true radius of the last arc, both flanks), no interference, same tip height on both flanks (<= 1 micron). Invalid profiles get `PENALITE`.
- Interference (`_interference`): arc-to-arc crossings between the complete outlines (both flanks, tip and root) of the teeth of both wheels near the line of centres, at 6 positions of the period and every centre distance, driven wheel placed by the working contact. Tangency (legitimate contact) is excluded (`EPS_INTERF` = 1 micron of overlap). Checked: 99-100 % of the profiles with negative backlash are flagged, 0 % of those with backlash > 0.3 degree, 0 on the reference profile. It costs about 1 ms per evaluation. Not checked: positions between the 6 sampled ones, tool (hob) interference / undercut, tooth thickness.
- Backlash uses the back flanks: `theta_max_avant(phi) + theta_max_arrière(-phi)` folded modulo one pitch (`contact.Engrenement.jeu` does the same with `dos=True`), identical to the numba kernel to 1e-3 degree with asymmetric teeth.
- `recherche.rechercher(...)`: several starts (the current profile, then random valid profiles ever further, up to 8 x `etendue_um`), a CMA-ES per start (population doubled every second start), profiles that still violate a constraint are repaired by freeing the back half-profiles of both wheels (driving half-profile unchanged, small regulariser towards symmetry, stops as soon as no violation is left), then the list of distinct optima (`separation_um` = 3 microns, valid ones first). `Optimum`: J, x, `mobiles` (with `points_dos` if asymmetric), per-centre-distance report, remaining violations. In practice the front search alone satisfied every constraint tried (backlash, tip radius): the repair only triggers when it cannot, and was tested on its own (raising the backlash from 0.93 to 1.5 degree with a 7 micron change of the back flank, front unchanged) and with a starved front budget.
- GUI (tab "Optimisation", menu Run > Optimiser F6 / Arrêter, buttons next to "Évaluer les critères et le temps"): search in a `QThread` (`Travail`), progress bar, table of optima (double-click = show in the Simulation tab, "Exporter…" = `.json`, `.Eno` refused for asymmetric teeth). Free arcs (not points) are chosen per wheel. The evaluation button shows the criteria at -ΔA/0/+ΔA, the constraints, the evaluation rate and the estimated duration. Project keys (`Projet.optimisation`): see `evolvente.OPTIMISATION_DEFAUT`; `jeu_max_admis` 0 = backlash of the start profile at +ΔA (+0.02).

**Choices settled by the user** (do not reopen them unless asked): `cma` accepted as a dependency; free *arcs* (curvature, length) rather than free points; minimum tip-arc radius as a user parameter (hob life); an interference makes a solution invalid; repair through the back half-profiles; save both half-teeth. The key name `profile_arriere` was chosen by Claude; the user may rename it (`arcs.CLE_DOS`). Always check final solutions with `contact.Engrenement`.

Design decisions, from measurements on `test.Cev` (89/12 teeth, 4 arcs, sample Hertz data):
- Variables = curvature and length of each arc, NOT the (x, y) of the points: the valid neighbourhood is 10-30 microns instead of 1-3 (moving ONE point by 1 micron rotates the following tangents by 3-12 degrees), and CMA-ES cannot even start on (x, y).
- Modify both wheels together: at equal budget (12000 evaluations) joint 5.6, driven only 7.1, driver only 8.7-9.3, alternating driver/driven 10.5-12.2 (the valley is diagonal because conjugacy is relative).
- CMA-ES reaches its floor in 10-20 k evaluations; differential evolution gives twice the objective at the same budget; an emulated BIDIOPT is far behind (4.8 vs 0.8 after 61 k vs 10 k evaluations, varC only). Several seeds end within 3 %. With torque variation alone the backlash doubled (2.7 -> 5.3 degrees): the backlash bounds are needed. Typical result: torque variation 6-13 % -> 1-2 %, Hertz 646 -> about 400 N/mm^2 at +dA.
- The optimum found on -30/0/+30 microns stays good at every intermediate centre distance and degrades outside the range (backlash negative at -45 microns); ENGOPT prints statistics on 9 centre distances (-40..+40 step 10), check final solutions the same way. Hertz pressure jumps with the centre distance (arc change).
- Other basins exist (starting from involutes with NW 14.5-25 degrees all reach about the same objective except 14.5, better by 10 %); a start that is invalid stops CMA-ES, hence the random valid starts. More arcs (6, 8) do not lower the floor and need more evaluations.
- Bug found by the user: `arcs.approximer` (KREBAX/LAGRMO) occasionally produced, for a specific number of arcs, a transition-point sequence that folds back on itself (non-monotonic), which fails `SWTE` (code 5, message ends with "essayez un autre nombre d'arcs dans l'onglet Création", which reliably worked around it: e.g. `Files/Evolve109-11.Cev` failed only at narcs=4, not at 3 or 5). **The cause was later found and fixed in `approximer` itself** (index-based breakpoints, see the arcs section); that file now gives 0.37 µm at narcs=4 instead of 92.6 µm, so the workaround should no longer be needed. The robustness added below is still worth keeping, since nothing guarantees an arbitrary profile is approximable at an arbitrary arc count. Before that, `Probleme.__init__` unconditionally required this reference profile valid (`calibrer=True`) to compute `f0`, so the whole multi-start search aborted on construction, trying none of the configured départs; `main.Fenetre._jeu_max_auto` (used whenever "jeu max" is left at its "auto" default) had the same unconditional requirement, called *before* `recherche.rechercher` even started, so it failed the search even for objectives that do not need `f0`. Fixed: `Probleme(..., calibrer=False)` defers `f0`; `rechercher` tries the profile as given first, and if invalid searches nearby (`_tirage`, escalating the search radius 1x/10x/100x `etendue_um`) for a valid substitute to use as `f0` reference and as départ 1, surfacing the substitution (`stats["depart_invalide"]`, shown as a warning in the GUI) instead of raising; the jeu_max "auto" resolution moved from `Travail.run()` into `rechercher`, after that same reference is established. `evaluer_criteres` (button "Évaluer les critères et le temps") is intentionally left strict: it must describe the exact loaded profile, not a substitute.
**Balance of the approach and recess action** (`SymL`, asked for by the user: an optimisation left
to itself gives much more contact after the line of centres than before, which is not always what
is wanted). The contact point travels along the line of action through the pitch point (`x = 0` in
the drawing frame), so the extent on one side is the approach and on the other the recess;
`L_AVANT`/`L_APRES` hold the two (mm) and `SYML` their ratio, 1 = balanced, 0 = all on one side.
The quantity was already in the kernel (it is Resulto's `SymL`, weighted by `GL` in the `.Eno`
line) but was computed and used nowhere. Two changes:
- it is now taken on the **extremes of the path**, not on its two end points as Resulto does. Index
  0 is the pair-handover position, where the contact can still belong to the outgoing pair, so the
  endpoint ratio jumps there: on `test.Cev` it gave 0.70 instead of 0.46 at the nominal centre
  distance, and `x(0) = x(ik-1)` at two centre distances out of three;
- `Objectif.equilibre` (a checkbox, not a third criterion, since the user wants it optional) adds
  `p_equilibre·(1 - SymL)²` to J, as a quadratic mean over the three centre distances -- the
  `((1-SymL)·GL)²` of Resulto. It is **not** normalised by the start profile, unlike the criteria:
  an already balanced start would give a zero reference. `p_equilibre = 1.0` means a mesh entirely
  on one side costs as much as the whole principal criterion does at the start; on `validation.json`
  the term is worth 0.15 at the start (J 1.50 -> 1.65) and 0.49 after an unconstrained optimisation.

How unbalanced it really is, measured on `validation.json` (already optimised, 89/12): approach
0.062 mm against recess 0.299 at the nominal centre distance, i.e. a balance of 21 %, and 2 % at
-30 µm, 98 % at +30. So the user's observation is exactly right, and it is worst at the tight
centre distance.

What ticking it buys, same search twice on `validation.json` (6 starts x 12000 evaluations, 69 s
each, principal `couple`, secondary `hertz` at 0.5):

| | balance -30 / 0 / +30 | approach (mm) | recess (mm) | varC at 0 | Hertz at 0 |
|---|---|---|---|---|---|
| box off | 6 / 19 / 51 % | 0.020-0.112 | 0.220-0.322 | 0.66 % | 354 |
| box on | 28 / 40 / 83 % | 0.073-0.146 | 0.177-0.263 | 0.72 % | 380 |

So it roughly doubles the approach action, and the price is 0.06 points of torque variation and
26 N/mm² of Hertz. Both results satisfy every constraint. Note the J values are not comparable
between the two rows (the second objective has an extra term).

The GUI shows the two extents and the balance per centre distance in the evaluation table, says
whether the term is in the objective, and the optima table has an `Équilibre` column (filled
whether or not the box is ticked -- it is a property of the result, not only a term of the
objective). The project key is `optimisation["equilibre"]` (0/1); old files without it read 0.

## Reverse direction, asymmetric teeth (`sens_inverse`) -- engine done, GUI pending

Asked for by the user: until now nothing said what happens **when the driven wheel becomes the
driver**. The load reverses, so the contact moves to the **back flanks** (drive side / coast side).
Two consequences the user stated and that drive the whole design: the back flank is then a working
surface, so it can no longer be used as a free repair variable -- if it jams, the solution is
rejected -- and one must be able to ask either for a normal reverse behaviour or for a deliberately
bad one (partially or fully blocking gears).

**The reverse case reduces to the existing engine**, which is what makes this affordable:
`contact.inverser(roues, profils, mobiles)` swaps the two wheels and swaps each tooth's two flanks
(the back becomes the active one, the front becomes the dos -- collapsing the dos instead would make
`jeu` measure a flank against itself). The result is the mirror image of the real configuration, and
a reflection changes neither speeds, efficiency, Hertz pressure nor backlash. Verified: at mu = 0 both
directions give exactly 100 %, the transmission ratio inverts exactly (-7.415725 against -0.134767),
at mu = 0.2 the reverse gives 95.37 % against 95.80 % forward (approach and recess swapping roles),
and the minimum backlash agrees between the two directions once converted by Z2/Z1 (3.8389 against
3.8474 degrees at +30 µm). The backlash *maximum* does not -- the modulo-one-pitch folding is on the
driven pitch, which differs -- so the backlash must be read from the forward run.

In the **numba kernel** it costs no new code: `Probleme.evaluer(..., inverse=True)` calls the same
`evaluer_lot` with permuted arguments (back flanks active, front flanks as dos, wheels swapped, and
with them `p0`, `th0`, `z`, `mo` and the Hertz data). Measured cost: exactly x2 (0.007 -> 0.013 s),
and it agrees with `contact.py` to 7e-4, the same spread as the forward direction. `bloc_ok` lets a
jammed position be *counted* rather than invalidate the profile (code 4), which is what the forward
direction must keep doing. New columns: `REND_MAX`, `FRAC_BLOQ`, `BLOC_DEBUT`, `REND_DEBUT`.

`Objectif.sens_inverse` has three modes, all checked at the three centre distances, all keeping the
forward backlash constraint (thickening the back flank eats the backlash, and a negative backlash
would jam the forward direction too -- the user's own condition):
`"reversible"` (no blocked position), `"couple"` (transmitted torque below `couple_inverse` over the
whole mesh) and `"blocage"` (the first `fraction_bloquee` of the period blocked).

**The criterion is the transmitted torque, not the efficiency** (corrected by the user), and that is
not a detail: the efficiency **cannot** serve here. The period passes through the pitch point, where
the sliding velocity is nil, so the efficiency is 1 there whatever the profile -- measured, its
maximum stayed at 0.9995 without moving by 1e-4 over 6000 evaluations at four different targets, and
it is reached at the position whose contact point is 0.001 mm from the pitch point. What is pinned
is the *product*: `rendement = rapport x couple`, so at the pitch point the torque equals `1/rapport`
and stays free (96.4 % forward and 103.6 % reverse on `validation.json`, for an efficiency of 99.98 %
in both cases). Lowering the transmitted torque therefore means raising the transmission ratio at the
same place, i.e. moving away from conjugacy -- which is exactly what a deliberately blocking gear is.
`COUPLE_MAX` and `COUPLE_MIN` (signed, relative to the ideal `Z2/Z1`) agree with `contact.py` to
0.2 %.

**The blocked fraction is a step function, hence gradientless.** With the penalty on `BLOC_DEBUT`,
CMA-ES gave up in 4 s with nothing blocked, at every target. `COUPLE_DEBUT` -- the transmitted torque
over the first `fraction_bloquee` of the period, counted negative where it blocks -- is its
continuous counterpart: it falls smoothly and crosses 0 exactly when that part of the cycle is
blocked.

**Mode `"couple"` works, and it is a band, not a ceiling** (the user: the torque must not exceed the
target but must *stay close to it* over the whole mesh -- a torque limiter has a flat
characteristic). The penalty is therefore two-sided, `max(0, COUPLE_MAX - cible) +
max(0, cible - COUPLE_MIN)`. Phase 2 on `validation.json` from the symmetric start, back flanks only,
20000 evaluations, starting band 83.19-104.45 %:

| target | resulting band | width |
|---|---|---|
| 95 % | 95.13-97.70 % | 2.6 points |
| 90 % | 92.51-95.31 % | 2.8 |
| 80 % | 80.04-95.78 % | 15.7 (not reached) |
| 60 % | 77.81-97.78 % | 20.0 (not reached) |

So the flattening itself works -- a 21-point band collapses to 2.6 and sits on the target -- and the
reachable **level** from a symmetric start with only the back flanks free bottoms out around 92-95 %.
Going lower is a design question, not a search-budget one (see below).

**Correction, after the user pushed back ("je ne vois pas ce qui t'empêche de diminuer encore
plus").** The 92-95 % floor reported here first was a floor of the *search radius*, not of the
physics, and the user was right. Measured: with `ECHELLE = 0.03` and CMA-ES at sigma = 1.5, a
random variable vector moves the back flank by **24 µm median** on a tooth 160 µm tall -- a 15 %
wiggle. Raising sigma does not help either (94.2 % at best at sigma = 3, and the displacement stays
19-77 µm, because CMA-ES shrinks its own step): the shape wanted is not in the basin of the
symmetric start, it has to be **constructed**, not searched. That is exactly why the sweep below is
the right instrument.

**The sweep (`balayer_blocage`), the user's own framing**: tooth count, module and tooth height are
imposed by the driving flank and are *not* variables; what is free on the back flank is its **slope**
and the **backlash**. `flanc_arriere(roue, pente, jeu_qe)` builds a back flank as an involute of
pressure angle `pente` with tooth thickness set by `jeu_qe`, keeping `mo`, `zz`, `ma`, `me` from the
front flank; `balayer_blocage` crosses the two ranges and reports, per combination, validity, the
forward backlash (negative = jams, rejected), whether it blocks at the **three** centre distances,
and the reverse torque band.

**The sweep is blocked by an `approximer` failure on steep flanks**, and this is a bug, not a
geometric wall. Measured on z = 89, `qe = 1.6`, 6 arcs:

| back-flank pressure angle | approximation error | radius sign inversions | `SWTE` |
|---|---|---|---|
| 20 deg | 0.00 µm | 1 | ok |
| 35 deg | 0.00 µm | 1 | ok |
| 50 deg | **184.83 µm** | 5 | rejected |
| 60 deg | 493.43 µm | 5 | rejected |

More arcs do not rescue it (114 µm at 8 arcs), and the error sits at index 33/150, i.e. at the
**root-arc to involute junction** (`qg = 0` halves it, to 86 µm). So this is the same family of
defect already fixed twice in `approximer` (index-based breakpoints, the hobbed-profile layout),
hitting a third sampling regime. Until it is fixed the sweep can only measure up to about 35 degrees
of back-flank slope, which is far short of the blocking regime. Trying several arc counts
(`NARCS_ESSAIS`, first chain `SWTE` accepts) is already in -- it is not enough here.

**Blocking: where the wall is.** The user's limiting intuition is right -- a tooth with nothing on
its back half (a radial back face, a ratchet) locks. The quantity that matters is the back flank's
inclination to the radial direction, and it has a closed form here: the half-pitch arc is always
`pi*m/2` whatever the tooth count, so the inclination is `atan((pi*m/2) / h)` with `h` the tooth
height. **It does not depend on Z**, only on the tooth height in modules, and approaching radial
needs `h >> pi*m/2 = 1.57*m`. Measured (30/90, m = 0.135, qe = 1.8, mu = 0.2, front flank 20 deg):

| lever | reachable before the geometry fails | what fails first |
|---|---|---|
| back-flank pressure angle | blocking appears at 58-70 deg, but **only at isolated centre distances** (2-4 % of the cycle at one of the three) | forward backlash goes negative from 58 deg, `SWTE` rejects the back flank from 50 deg |
| tooth height | inclination 52.6 deg at `h = 1.2 m`, 38.1 deg at 2 m (reverse torque minimum down to 75.1 %, nothing blocked) | tip against root at `h = 3 m`, before the inclination gets small enough |

**No profile valid at the three centre distances blocked.** The walls are tip-to-root interference,
the forward backlash and `SWTE` -- not the optimiser, and not the evaluation budget. So looking for
"the limit plus a safety margin" is a **continuation on the design parameters** (tooth height,
back-flank steepness, backlash allowance, module and tooth count), with a validity gate at the three
centre distances -- not a CMA-ES search on the arc variables. The test geometry above may simply be
unfavourable; the right tool to build is the sweep itself, so the limit can be located on the user's
own proportions.

**What gives the lever on the torque level: the back flank's starting steepness, not the search.** Phase 2 from
the *symmetric* start (front flanks frozen, back flanks free, 20000 evaluations, mu = 0.2) only takes
the maximum reverse torque from 104.45 % to 95.42 % and stops there -- far from a 60 % target.
Starting the back flank from a **steep involute** changes the regime completely (30/90 teeth,
qe = 1.8, mu = 0.2, forward flank 20 degrees):

| back flank | forward torque | reverse torque | share of the cycle blocked |
|---|---|---|---|
| 20 deg | 96.69-99.99 % | 97.08-100.00 % | 0 % |
| 50 deg | 96.69-99.99 % | 80.63-104.68 % | 0 % |
| 60 deg | 96.69-99.99 % | **-320.95-139.35 %** | 4.0 % |
| 65 deg | 96.69-99.99 % | -467.45-161.87 % | 2.7 % |

Two things to read there. The **forward direction is rigorously unaffected** by the back flank, as
it must be -- the two flanks are independent, which is the whole point of an asymmetric tooth. And
the blocking geometry makes the torque *swing*: it collapses and reverses over part of the cycle
while exceeding 100 % elsewhere. So "transmitted torque below 60 % over the whole mesh" (a bound on
`COUPLE_MAX`, i.e. a torque limiter) and "partially blocking" (`COUPLE_MIN` going to zero or below)
are **two different devices**, and the geometry that gives the second pushes the first the wrong way.
Both measures are available; which one the mode should bind is the user's call.

Friction matters too, self-locking being a friction phenomenon: from the symmetric start, mu = 0.2
blocks nothing while mu = 0.5 blocks 19-25 % of the cycle. A starting sigma of 8 is too bold --
every profile comes out invalid.

Still to do: the GUI tab (between Création and Simulation) holding the back-flank geometry and the
reverse-direction requirement, the project keys, and `recherche`'s two-phase driver (forward, then
back flanks only from the symmetric start). Note the **back-flank repair must be disabled** whenever
a reverse requirement is set: the back flank has a job then, and `recherche`'s repair would undo it.

- Not done: tool (hob) interference, `PZ` > 2 tooth pairs, the `GB` blocking option (replaced by `jeu_min`), identical driver/driven (`LI`), reading `.Eno` optimisation parameters beyond `-dA +dA GW` (the `GL` symmetry weight is now `Objectif.p_equilibre`, fixed in code rather than read from the file).

## Circular arcs and `.Eno` export (`arcs.py`)

`approximer(profil, zz, qg, narcs)` = `KREBAX` + `LAGRMO` + `KRBGS`: transition points start at the root-arc end, the straight root part (if `qg>0`), `narcs` evenly spread involute points and the tip-arc end; `lagrmo` moves the interior points (Newton with a diagonal Hessian and a Lagrange multiplier for the final slope 0) to minimise the squared distances of the 150 profile points; the first tangent angle is `-pi/Z`. The arcs are fully described by the transition points, which is why the `.Eno` file only holds those. Checked: arcs pass through their end points, tangential junctions to 1e-6 degrees, first and last points are the profile ends, error at most about 1 µm with 3 arcs and 0.1 µm with 6 (`test.Cev`).

**`approximer` only worked on `evolvente.evolve`'s own point layout** (user-reported, via the hobbed profiles): KREBAX picks its initial breakpoints by *index* -- the first `NI//10` points are assumed to be the root arc, the last `NI//10` the tip arc, the involute in between -- which is true of `evolve` and of nothing else. On a profile sampled by arc length (a wheel hobbed by `fraise.tailler`, an arc chain re-sampled through `echantillonner`) those indices land anywhere, and `lagrmo` only nudges points by small Newton steps, so the fit stalls: 3.9-4.8 µm instead of sub-micron on `validation.json`'s hobbed wheels, with the arc radii changing sign at nearly every arc (up to 7 inversions), which `SWTE` rejects -- the "n'importe quoi comme profil" seen at 5-6 arcs. The same defect, not the arc count, caused the documented `Files/Evolve109-11.Cev` failure at narcs=4 (92.6 µm, 2 sign inversions). Fixed by making the choice of breakpoints sampling-independent: two extra starting layouts (equal arc length, equal cumulative tangent rotation -- `_mesures`/`_bornes`) beside KREBAX's, and from each a chain of re-cuts that equalise the max error per arc (`_redistribuer`, equidistribution in `error^(1/3)`), every layout being fitted by `lagrmo` (`_ajuster`) and the best kept. "Best" = fewest radius sign inversions beyond one first, then smallest max error, so a fit that `SWTE` would reject never wins over a valid one. Measured (narcs 2..8): hobbed wheels 3.9-4.8 µm -> 0.2-1.8 µm, always 1 sign inversion; `Evolve109-11` at narcs=4 92.6 -> 0.37 µm; on every `evolve` profile tried (`test.Cev`, `f.Cev`, `Evolve109-11`, `Evolve102-12`, `Evolve8-96`, `Evolve16-07-0100`, both wheels) identical to the old result or better, the only differences being at the 10-picometre level. End to end, re-using the hobbed pair in the Simulation tab: torque variation 6.4-22.2 % (and a *negative* backlash at 2 arcs, i.e. blocking) -> 2.2-4.9 % with a backlash of 2.35-2.59° against the original profile's 2.38°.

Cost: `approximer` now fits 6-18 layouts instead of one, so `lagrmo` was made much cheaper -- its finite differences are evaluated as one batch (`_krbgs_lot`, `_sqab_lot`, `_affectation`, all vectorised over the batch; checked equal to the one-at-a-time version to 4e-16, and `sqab` equal to a literal point-by-point SQAB to 1e-15 on 480 chains including deliberately folded ones), and it stops when the relative gain in max error falls under `GAIN_MINI` (0.2 %), the Fortran's own stopping rule (error rises) being reached only after ~100 iterations that gain fractions of a nanometre. Net: 0.2-1.4 s per call, the same range as before (the old code hit 1.5 s on the same cases). The Création tab caches per parameter set (`_arcs_de`), so the cost is paid once per edit. `Mobile` (name, module, Z, transition points, flags, Hertz data) is the common representation for reading and writing. `ecrire_eno` reproduces the Fortran layout (`F10.6,2I10` header, `2F12.7,I6` points with flag 0 at both ends, Hertz lines, the fixed `-dA +dA GW ZV` line, CRLF) and was read back with the reader of `~/DEV/Python/001-Eng_open/fileIO.py` (`read_eno`). Its Hertz lines are assigned to master/slave in the opposite order to the Fortran (first line = driver), which is that program's convention, not ours. Files without Hertz data get the sample values `HERTZ_DEFAUT` of `main.py` in the panel.

Editing the arc count of an already-optimised profile (`arcs.ajouter_point`/`retirer_point`, Simulation tab group "Modifier les arcs"), so a fix or a refinement does not require re-running the Création + Optimisation pipeline from the raw involute: `ajouter_point(mobile, i)` inserts a point at the angular midpoint of arc `i` (same formula as `echantillonner`); since a circle is determined by a point and the tangent there, `krbgs` reconstructs that same arc as two halves of the *same* circle and leaves every other arc untouched -- checked both on the arc centres/radii (equal to 1e-12 mm, all arcs, including the tip arc closed by `korzf`) and end-to-end through `contact.Engrenement` (criteria identical to 1e-9 relative). So splitting is always safe (never changes the shape, never turns a valid profile invalid or vice versa) and only gives the optimiser one more free arc; the GUI applies it immediately, no confirmation needed. `retirer_point(mobile, i)` removes the point that starts arc `i` (or that ends arc 0, which has no previous arc), merging two arcs into one that generally does not pass through the removed point any more -- the GUI shows the resulting deviation (`arcs.ecart_um`) and validity (`optim.profil_valide`) and asks for confirmation before applying. Both keep `points_dos` (asymmetric teeth) the same length as `points` by mirroring the operation on it, which `optim.Probleme` requires (`decoder` indexes the back arrays with the front arc's `libres` indices).

## Simulation report (`rapport.py`)

`ecrire_rapport(chemin, roues, mobiles, mu, mt, da_neg, da_pos, optimisation, origine, pas_um=10.0)` writes a
`.docx` report (File > "Exporter un rapport (.docx)…" in the GUI, using the profile currently loaded in the
Simulation tab): profile and simulation conditions, the optimisation conditions if `evolvente.Projet.optimisation`
is passed, a worst-case summary over the whole `da_neg..da_pos` range, a detailed table every ~10 µm, three
engagement plots (nominal, minimal, maximal centre distance, same drawing as the Simulation tab) plus a trend
chart of torque variation / Hertz / backlash vs centre distance, a geometry appendix (transition points and arc
radii), and a closing list of things not yet in the report (Fortran-numbers comparison, interference column, hob
tip-radius check on this specific profile) offered as proposals rather than assumed. Built with `python-docx`
(installed in the `.eng` venv) rather than the docx-js/Node tool this environment's docx-authoring skill defaults
to: a report generated from a button in the shipped PySide6 app must not depend on a Node.js runtime being present
on the user's machine. All figures/numbers come from `contact.Engrenement` (the engine validated against the
Fortran program), never from `optim.py`. Rendering the generated file to PDF for a visual check
(`scripts/office/soffice.py` from the docx skill) failed in this sandbox because `libreoffice-writer` is not
installed and installing it needs `sudo` (not available non-interactively) -- verified instead by reading the
four embedded PNGs directly and by dumping the document structure (headings/paragraph text/table cells) with
`python-docx`; a real visual check on the user's machine (which has proper Word/LibreOffice) is still worth doing
once.

## Hobbing cutter (`fraise.py`)

Repo of `Abwfrac.F90`/`Abwfrap.F90` (wheel -> hob), `Radabwp.F90` (hob -> hobbed wheel), `Krvabc.F90`/`Krvabp.F90`
(max deviation) and `Spliph.F90` (Hermite spline), for a rack (a wheel of infinite radius): the given flank
rolls without slipping against a straight line, sweeping a family of curves; the envelope of that family is the
conjugate rack/hob tooth, in the rack's own frame (origin on the pitch line, not the wheel centre -- a point
users of this feature must get right, flagged explicitly by the user). Rolling the (now fixed) hob back against
a target wheel blank (any tooth count, any profile shift) gives the "practical" cut profile, compared to the
reference with a signed max-deviation function (excess = material left outside the reference, i.e. under-cut;
penetration = hob went inside it, i.e. over-cut).

Deliberate simplification vs the Fortran: `Abwfrac.F90` computes the wheel -> hob envelope analytically on
circular arcs (circle/horizontal-line intersection); here both directions use the same segment-intersection
envelope on a finely sampled polyline (`arcs.echantillonner`, 2000 points), unifying the code at the cost of a
discretization error made negligible by the fine sampling (well under manufacturing tolerances; verified below).

Ported: `_famille` (the rotate+translate family), `_enveloppe` (vectorised lower envelope, replaces the
Fortran's nested loops), `spline_hermite` (`SPLIPH`: parametric cubic spline by arc length, X natural, Y
Hermite with prescribed end slope -- the curve is first rotated to make the prescribed slope near-horizontal,
needed because a hob or hobbed-wheel profile is not a single-valued function of x), `profil_fraise` (`ABWFRAC`/
`ABWFRAP`), `tailler` (`RADABWP`, oy: profile shift in mm, effective module = `fraise.module + 2*oy/zz`, matches
`Fraesen.F90`'s O/M choice), `ecart` (`KRVABP`-style signed perpendicular distance + location).

`.Frs`: a wheel pair per file (menante, menée), matching the real Fortran format exactly -- checked
byte-for-byte against a genuine ENGOPT-written file the user supplied (`Result/Validation_T.Frs`):
label `A30`, header `F10.6,I6,F10.6` + the literal text `     MO, ND, PC`, point count `I6`, points `2F12.7`,
then an optional ` Pression Hertzienne` line and two `F10.1,2F10.3` lines, CRLF, latin-1. `lire_frs` returns
`(menante, menée, hertz_ou_None)`, `ecrire_frs(chemin, (f_m, f_e), hertz=None)` writes it; round-trips a real
file exactly (`ecrire_frs(lire_frs(...)) == fichier d'origine`, bit for bit). One hob is still reusable to cut
any tooth count/`oy` (`tailler`'s `zz` and `oy` are independent of `zz_origine`), but the *file* holds a pair
because that is what the real tool writes and the user needs to reuse hobs already produced that way -- reading
is not limited to files this module wrote.

Sign/frame pitfalls found while porting, first "fixed" wrongly by reasoning about the frame abstractly (which
kept giving wrong answers) and only actually fixed once a real Fortran-computed hob became available to compare
against point-by-point instead of only against this project's own theoretical/round-trip checks:
- `_famille`'s rotate+translate is `(-cos,-sin ; -sin,cos)`, not a standard rotation matrix -- `ABWFRAC`'s own
  convention, correct as-is in this project's frame too despite looking like it should need re-deriving for a
  different handedness.
- The two analytic endpoint formulas of `ABWFRAC` (root at `x=-pl/2`, tip at `x=0`, `pl` the hob's own pitch)
  transfer *exactly* as well, unchanged -- an earlier attempt to "fix" a large-looking mismatch by replacing
  them with a local linear extrapolation of the envelope scan (`_extrapoler`, now removed) was solving the wrong
  problem: the mismatch was in the *validation method* (fitting a line through noisy edge points of a coarse
  scan), not in the formulas, which match a real Fortran hob to within 0.06-0.28 microns once trusted outright.
- `RADABWP`'s own transform for `tailler` (`FX=-FR(1,I)+PS`, then a *standard* rotation, not the mirror one) was
  ported literally too, and is algebraically the same map as the inverse this project first derived by hand
  (`_famille(..., inverse=True)`: rotating the translation, `(cos(wi)*ps, sin(wi)*ps)` instead of `(ps, 0)`,
  since a fixed rack-direction translation seen from the wheel's rotating frame is not fixed either) -- both
  give identical numbers, so either derivation is fine, but trusting the literal port needed less second-guessing.
- **`_enveloppe` must only accept *ascending* segments** (`y(i) <= fy <= y(i+1)`), the test both `ABWFRAC` and
  `RADABWP` use, not `min(y) <= fy <= max(y)`. This was the port's one real bug, found by the user on
  `test_jeuxmini.json` ("le profil menant dérive complètement"): profiles are given root to tip, so a segment
  that descends belongs to a copy rotated far enough to have turned over, and taking the global minimum x over
  those picks a point on the far side of the wheel. On that file's 8-tooth driver the generated tip radius came
  out at 0.94 mm instead of 0.42 (891/820 microns of deviation); with the filter, 0.03/0.85 microns. It showed
  up there and not on `validation.json` because the damage needs a copy rotated far enough to reach the height
  band of the tooth being generated, which happens when the hob tooth is angularly wide relative to the wheel
  (half-pitch / rolling radius 0.39 on the 8-tooth wheel against 0.26 on the 12-tooth one).
- Two consequences of that fix, both verified rather than assumed:
  - **`periodes_defaut` no longer scales with the tooth count**: it returns the Fortran's 9. The +-90 degree cap
    existed only because of the envelope bug (the 2.3 mm error on a 12-tooth wheel at 9 periods). Now 5, 7, 9
    and 11 periods give exactly the same profile on every wheel tried (7 to 150 teeth); 3 differs because it
    under-sweeps and misses part of the neighbouring teeth's undercut. Cost at 9 periods: 0.27 s per wheel
    against 0.11 s at 5.
  - **`RADABWP`'s root endpoint is exact after all**, and the note that used to sit here was wrong. Its `TY` is
    the *height* of the lowest tooth point, read from the family, not a radius: `TY = creux*cos(pi/zz)` (checked
    against the family member `RADABWP` actually uses, to under 10 nm on four wheels), so its
    `(TY*tan(pi/zz), TY)` is exactly the polar point `(creux*sin(pi/zz), creux*cos(pi/zz))` this port uses, and
    both land on the reference root point to 7 decimals. What was genuinely wrong earlier was substituting this
    port's radius-valued `creux` into the Fortran's height formula.
- `tailler`'s root/tip must come directly from the hob's own profile (`fin[0,1]`/`fin[-1,1]`, mirroring how
  `profil_fraise` reads the wheel's own root/tip), not from the min/max of the swept family (which mixes in
  unrelated neighbour-tooth copies and gave a badly wrong, 400+ micron root position). The *height* scan then
  runs from `TY = creux*cos(pi/zz)` to `sommet`, as in the Fortran -- scanning from `creux` instead skipped the
  bottom of the tooth (1.7 % of the height at 89 teeth, 8 % at 8 teeth).

Validated: `profil_fraise` against the real `Result/Validation_T.Frs` (both wheels of `validation.json`, the
arc profile that produced it, `n_positions=151`) to 0.04-0.28 microns, and reproduces that file bit-for-bit
through `ecrire_frs`/`lire_frs`. Full round trip (wheel -> hob -> hobbed wheel, same Z, `oy=0`): 0.12/0.31 and
1.27/2.74 microns (excess/penetration) on `validation.json`'s 89 and 12 teeth, 0.03/0.85 and 0.11/0.51 on
`test_jeuxmini.json`'s 8 and 150; on the involute files (`test.Cev`, `f.Cev`, `Evolve16-07-0100.Cev`, tooth
counts 7 to 89) 0.02-0.78 excess and 0.01-2.6 penetration, the penetration always in the root fillet. On a
theoretical involute the hob is a straight line at the pressure angle: 20.000 degrees with a 0.000 micron
residual on the 89-tooth wheel (20.203 and 1.44 microns on the 12-tooth one, that residual being the 8-arc
approximation of the involute, not this module). Arc-approximated profiles
give the same few-micron round-trip figure regardless of arc count, because the arc count only affects the fine
sampling step, not this module's own math (so this specific check does not, by itself, validate the arc
approximation -- that is covered separately by `arcs.approximer`'s own convergence, documented above).
Profile-shift (`oy`) and different-Z cutting with the same hob both behave as expected (effective module shifts
linearly with `oy`, tooth height stays constant across `zz` for a fixed hob) but have no reference profile to
compare against, so `ecart` is only computed when `zz` and `oy` match the profile that produced the hob (mirrors
`Fraesen.F90` zeroing `US` in that case) -- loading a hob from a `.Frs` file also disables the comparison, since
the originating profile is then unknown.

GUI: tab "Fraise" (`main.py`), independent of the project's saved state (`.json`) -- none of its fields mark
the project modified. "Créer une fraise" computes *both* wheels' hobs at once from the Simulation tab's
`self.mobiles` (`self.fraise_paire`), matching `Fraesen.F90`'s option "F" (`DO K=1,2`); the "Roue" combo picks
which one to view/cut with and switching it drops any previous taillage (it was for the other wheel).
"Fichier fraise" saves/loads the `.Frs` pair (+ Hertz data from `self.mobiles` when saving, if known). "Tailler
un mobile" (target Z, `oy` in µm or a target module -- module wins if it differs from the hob's own) draws the
cut profile, overlaid with the reference and the excess/penetration markers (star/triangle) when the comparison
is relevant, or a plain note when it is not (different Z/`oy`, or a hob loaded from file).

## DXF export (`dxf.py`)

Uses `ezdxf` (already in the `.eng` venv, not a new dependency). `ecrire_profil(chemin, mobile, n_dents=None)`
exports a mobile's full contour as real ARC/LINE entities (not a sampled polyline), rotating and mirroring the
arc geometry itself (`_ajouter_chaine_arcs`) for each of `zz` teeth, both flanks (front + back, or front mirrored
twice if symmetric) -- checked entity count (`2*(n_arcs)*zz`) and by rendering back with
`ezdxf.addons.drawing`/matplotlib into a closed gear outline. An arc whose radius exceeds `RAYON_DROITE` (1e4 mm,
e.g. a near-straight root transition) is written as a LINE instead, since DXF ARC needs a finite radius and
a circle that large is a straight line in practice. `ecrire_courbe(chemin, points, n_periodes, pas)` exports a
generic curve (the hob profile, which is not an arc chain, it is the envelope of a sweep) as a mirrored,
repeated LWPOLYLINE. DXF ARC angles are always CCW start->end; since the arc chain's own short-way sweep
(matching `arcs.echantillonner`'s convention) can go either way, the two endpoint angles are swapped when needed
so the entity draws the same (short) arc regardless of the chain's original direction.

Rings: both writers step by `2*pi/abs(zz)` over `abs(zz)` teeth, the same rule as
`evolvente.contour_dents` (the signed pitch made the contour jump from one tooth to the next). The
arc geometry itself needed nothing. Checked on zz = 30, -40, -90: the entity count is still
`2*n_arcs*|zz|`, consecutive polyline teeth meet to 1e-15 mm, and the exported contour's radii match
the profile's exactly.

GUI: File > "Exporter les profils en DXF…" writes the Simulation tab's two mobiles to `<base>_menante.dxf` and
`<base>_menee.dxf` (real arcs, all teeth). Fraise tab: "Exporter .dxf…" (fichier fraise group) writes the
selected hob's profile (5 periods, LWPOLYLINE); "Exporter la roue taillée (.Eno / .dxf)…" (arc-conversion group,
"Nombre d'arcs" field) re-fits the hobbed profile with `arcs.approximer(..., qg=0, ...)` (no straight root
segment -- the root shape already comes from the hob's own envelope) into a `Mobile` with the effective module
(`fraise.module + 2*oy/zz`) and the origin wheel's Hertz data if known, then offers `.Eno` (needs both wheels
hobbed, since `.Eno`/`lire_arcs` is always a pair like every other arcs format here) or `.dxf` (one wheel is
enough). Checked: `.Eno` round-trip through `arcs.lire_arcs` (module, Z, points, Hertz all preserved).

`self.taillage` (`main.py`) is now a `[roue 0, roue 1]` list, one slot per wheel (`{"zp", "ecart", "zz", "oy"}`
or `None`), not a single "last hobbed wheel" -- switching the "Roue" combo used to discard the other wheel's
result, which made it impossible to hob both wheels before exporting the `.Eno` pair; now each wheel keeps its
own result until a new hob pair is computed or opened (`self.taillage = [None, None]`).

Closing the loop back into Simulation/Optimisation (the user found the wheel->hob->hobbed-wheel pipeline
tedious and asked for both a smoother workflow and a way to reuse the hobbed result without a file round trip):
`tailler_les_deux` hobs both wheels in one click with each one's own `zz_origine` and `oy=0` (the common case --
checking how an already-optimised pair survives hobbing -- needs neither a different tooth count nor a profile
shift), reusing whatever resolution is in the "Périodes/Positions/Points" fields and ignoring "Z cible"/"Déport"
(which stay for the one-wheel-at-a-time advanced case, `tailler_mobile`, unchanged); factored the common part
into `_tailler_roue(w, zz, oy, ...)`. `utiliser_taillage_simulation` re-fits each hobbed wheel to arcs
(`_taillage_vers_mobile`, same "Nombre d'arcs" field as the export) and sets `self.mobiles` to the result
directly, no file involved; a wheel that was not hobbed keeps its previous profile (from `self.fraise_source`,
the pair that produced the current hobs, or `self.mobiles` if that is unknown too, e.g. hobs loaded from a
`.Frs` file) so hobbing only one side of the pair to check it is a normal, supported use, not treated as
incomplete. Checked: hobbing both then reusing runs a normal simulation on the result; hobbing only one wheel
keeps the other's `Mobile` object identical (`is`, not just equal) and still simulates.

Two things the user found on that loop. (1) `_dessiner_fraise` drew the red dashed reference profile by joining
the `Mobile`'s ~7 transition points with straight lines, while the deviation shown beside it was measured on
`echantillonner(2000)`, i.e. on the real arcs -- a visibly polygonal curve that contradicted its own numbers.
It now samples the arcs too (`echantillonner(600)`); that was the only place in the GUI drawing an arc chain
as a raw point list (the hob's own profile legitimately is a polyline). (2) Reusing the hobbed pair gave
"catastrophic" simulation results (3 % torque variation against 1 % before hobbing): first the `approximer`
defect above, then -- once that was fixed and it was still 3 % -- the re-fit itself, see below.

**The hobbing costs nothing; converting the hobbed profile back to arcs is the whole loss.** Measured by
running the *same* (polyline) engine on both sides, `validation.json`: reference resampled 1.3534 % torque
variation, hobbed 1.3733 % at `n_positions=61` and 1.3454 % at 241. So the wheel -> hob -> wheel round trip is
worth 0.02 points. The 2.7 µm penetration on the 12-tooth wheel sits at 16 % of tooth height, below the end of
its concave root arc (24 %) and well below the contact zone (43-65 % on that wheel, 52-89 % on the 89-tooth
one): it is root undercut, on a surface that touches nothing, and it does not shrink with `n_positions`
(2.761 -> 2.728 µm from 61 to 481) although the *excess* does (0.87 -> 0.10 µm). Everything above came from
`approximer` landing on a *different* 6-arc chain than the original: in the contact zone that chain was 1.43 µm
off against 0.16 µm for one keeping the reference's arc structure, and at this sensitivity that is the
difference between 3.2 % and 1.2 %.

Fixed with `approximer(..., reference=<points de transition d'origine>)`, passed by `_taillage_vers_mobile`
whenever the origin wheel is known (hobs computed here, not loaded from a `.Frs`) and has the same tooth count.
When the reference describes the same number of arcs, `lagrmo` starts from **its interior transition points
themselves** (with the hobbed profile's own endpoints) and that fit is kept outright, without competing on a
scalar. Two things had to be got right, both found by measurement: snapping those points to the nearest sampled
index instead fails, because a short arc (a tip, a fillet) then falls between two adjacent indices, is
constrained by a single point and comes out with an arbitrary curvature (radius +3.3 -> -0.16 mm, `SWTE`
rejected); and ranking the reference fit against the others *on any scalar deviation* fails too -- tried on the
max deviation and on the active-flank deviation, both regularly prefer a neighbouring chain that is just as
"accurate" and doubles or triples the torque variation. Result at the reference's own arc count: 1.357 / 1.365 /
1.358 % for `n_positions` 61 / 101 / 241, against 1.346 % for the reference -- stable, and matching the
polyline measurement. It is a real re-fit, not a copy: the 89-tooth wheel's points do not move (0.00 µm, it is
barely undercut) while the 12-tooth wheel's move up to 13 µm, its backlash drops 2.382 -> 2.362° and its Hertz
pressure rises 353 -> 372 N/mm². At any *other* arc count the reference layout does not apply and the generic
search runs, giving 1.7-3.2 % -- a genuinely different shape, so the Fraise tab now proposes the origin wheel's
own arc count in the "Nombre d'arcs" field (`_proposer_n_arcs`, set when the hobs are computed and when the
"Roue" combo changes). `ProfilArcs.ecart_flanc` (deviation above `arcs.rayon_flanc`, the end of the last
concave arc) is reported beside the total one in the GUI messages (`_texte_ecart_arcs`), since the total is
dominated by the root fillet and gives a false idea of the fit.

**Simulating a hobbed wheel as a point sequence, not as arcs** (asked for by the user: what a hob produces is an
envelope, so the point sequence is the faithful shape; the arc fit is only needed to restart an optimisation or
to write a `.Eno`). `main.Fenetre` now holds either `self.mobiles` (a pair of `arcs.Mobile`) **or**
`self.mobiles_points` (a pair of `ProfilPoints`: name, module, Z, the half-profile points, Hertz data), never
both -- `_profils_courants` returns whichever is set and `simuler` passes `mobiles=None` for the point one, so
`contact.Engrenement` uses its polyline engine. Fraise tab, group "utiliser telles quelles (suite de points)":
one button simulates the hobbed pair directly, another exports it to DXF as a polyline
(`dxf.ecrire_polyligne`, one LWPOLYLINE per tooth, consecutive teeth meeting exactly on the gap axis -- checked
to 4e-16 mm and by rendering back). A wheel that was not hobbed keeps its own profile, resampled. In point mode
the Optimisation tab, the `.docx` report and the arcs export (`.json`/`.Eno`) refuse with `MSG_ARCS_REQUIS`
pointing at the arc-conversion group, the "Modifier les arcs" group is disabled and its lists cleared, and
`File > Exporter les profils en DXF` writes polylines instead of ARC entities. The Simulation tab is no longer
called "Simulation arcs" and its banner says which engine ran.

For that mode to be complete, `contact.py` computes the **Hertz pressure with polyline profiles** too:
`rayon_courbure(ligne, i)` fits a least-squares circle (Kåsa) over +-`DEMI_COURBURE` (25) fine points and
returns the signed radius with `flanc_arcs`' convention (convex positive), and `_contact` stores the two radii
in `self._courbures` exactly as the arc engine does; `Engrenement` now takes its Hertz data from `mobiles` when
it has them and from `roues` otherwise (a `.Cev` Hertz block, a `.json` "hertz" key), so the Création-tab
engine could use it as well. Checked against the arc engine on the same geometry (`validation.json` sampled to
150/400/1500 points): median deviation 0.00-0.05 %, max Hertz to 0.6 %; only 12 of 150 positions differ by more
than 2 %, in two short runs where the contact crosses an arc junction -- a curvature discontinuity that an
arc chain has and a hobbed profile does not.

Result, and the real answer to "does my optimum survive hobbing" (`validation.json`, reference pair at 1.3458 %
torque variation / 2.3818° backlash / 95.805 % efficiency / 352.7 N/mm²): hobbed then simulated as points,
1.3664 / 1.3594 / 1.3456 % torque variation for `n_positions` 61 / 101 / 241, backlash 2.3812-2.3816°,
efficiency 95.810 %, Hertz 363.2 / 353.7 / 352.4. It converges to the reference as the hobbing resolution
rises, and Hertz is the criterion most sensitive to it (it depends on curvature), so a Hertz figure on a hobbed
wheel deserves `n_positions` well above the default 101. A point simulation costs 0.38 s against 0.20 s in arc
mode. Note that the backlash is still evaluated at 10 driver positions only with the polyline engine (150 with
arcs), which is enough here (2.3815 against 2.3818) but is worth remembering.

## Internal ring gears (`zz < 0`) -- in progress

Asked for by the user, to be built in stages (geometry, then drawing, then contact, then
optimisation, then a skiving tool). Done so far: `evolvente.evolve` and the drawing, `arcs.py`,
both engines of `contact.py`, `optim.py`'s numba kernel, and the active flank of the skiving tool
in `fraise.py` -- each stage with its own section below. `dxf.py` too. What stays external-only
is the rack half of `fraise.py` (`profil_fraise`, `tailler`, the `.Frs` file): a rack cannot cut a
ring at all, and the GUI refuses that combination rather than computing it. The ring may be either wheel (the user asked for both
roles), and the skiving tool is to be done parallel-axis first, crossed axes afterwards.

Convention: `zz < 0` is an internal ring. The tooth points at the wheel centre, so the root is at
the *large* radius (on the gap axis) and the tip at the *small* one (on the tooth axis); `ma` stays
the involute extent on the root side and `me` on the tip side. Everything reduces to one sign
`sg = +1/-1` in `evolve`:
- involute start radius `rt - sg*ma*mo`, end radius `rt + sg*me*mo` (both clamped to the base
  circle), so for a ring it runs from the outside inwards and the involute parameter *decreases*;
- the involute is taken on its mirror branch (`sg*x`), because for a ring the angle to the tooth
  axis must decrease as the radius decreases, the opposite of an external tooth;
- the tip arc closes on the tooth axis at `w -> 0` (external, the far intersection with the axis)
  or `w -> pi` (ring, the near one);
- the root arc is centred on the gap axis and meets it at `-sg*pi/2 + pi/zz`, i.e. on the wheel-centre
  side for an external tooth and the far side for a ring.
- `contour_dents`/`points_dents`/`traits_dents` step by `2*pi/abs(zz)`, never the signed pitch: a
  tooth is drawn from one gap axis to the other, so teeth only join end to end in that direction
  (with the signed pitch a ring's contour jumped 0.9 mm between consecutive teeth).
- `entraxe_signe(roues, da)` is the new single place that knows `rt1 + rt2` is *negative* for an
  internal pair (both centres on the same side of the pitch point, and the wheels turn the same
  way): `da` moves it away from zero, it does not add algebraically. `rayons_fonctionnement` and
  `_placement` go through it; `contact.Engrenement` does **not** yet.
- **the drawing gives each wheel the half-turn its own sign needs** (`_placement`): the teeth must
  point at the pitch point, which is above wheel 1 when `z1 > 0` and below it when `z1 < 0`, and
  the other way round for wheel 2 -- so `d1 = pi if z1 < 0`, `d2 = pi if z2 > 0`, wheel 1 drawn at
  `angle1 + d1` and wheel 2's default at `d2 - pi/|z2|` (the half-pitch offset that puts a tooth of
  one in a gap of the other). `d1` was missing at first, reported by the user: a *driving* ring had
  its three drawn teeth at the top of the rim instead of the bottom, so it meshed with nothing.
  It is the same `delta1` that `contact.Engrenement` applies to the driver.

Checked on `zz` = -12, -40, -89, -150 and, unchanged, on the external wheels of `test.Cev` and
`f.Cev`: root point exactly on the gap axis, tip exactly on the tooth axis, involute start/end radii
as designed, radius monotonic in the right direction, tooth half-thickness at the pitch circle equal
to `pi/(|zz|*(1+qe))`, and the involute exact -- the polar angle follows `inv(alpha)` to 0.00 nm
(the earlier "distance from the centre to the tangent equals `rb`" test was simply wrong: it is the
*normal* that is at `rb`, the tangent is at `sqrt(r^2-rb^2)`). Tangency at the three junctions is
only checkable on a refined profile: the chord bias falls as 1/ni (8.5 degrees at ni=150, 1.9 at
600, 0.46 at 2400) and is identical for external and internal teeth, so the test runs at ni=2400
with a 0.6 degree tolerance. Non-regression: `evolve` still matches the Fortran `Result/f.Egr`
(NI=550) to 7e-8 mm, and the whole external pipeline (Création, simulation, optimisation, hobbing,
DXF) is unchanged.

**Arc conversion (stage 3a), done.** Three things in `arcs.py`, none of which changes an external
result:
- `angle_depart(zz)` = `-pi/|zz|` -- the *backward* start tangent that `krbgs` expects, for both
  signs (measured: the traversal tangent is `pi - pi/|zz|` on external and internal profiles
  alike). `angle_fin(zz)` = 0 or `pi`: a ring's profile reaches the tooth axis going inwards, so
  the target end tangent is reversed, and `lagrmo`'s angle constraint now wraps the difference
  into (-pi, pi] (`_ecart_angle`) because `krbgs` returns +pi or -pi indifferently.
- `convexite(zz, ra)` = `ra` or `-ra`. `krbgs`' sign is relative to the direction of travel; for a
  ring, whose material is outside, the physical convexity is the opposite. With it the measured
  pattern is root concave, **active flank concave**, tip convex -- which is what a ring is. Since
  multiplying by a constant does not change the *number* of sign changes, `SWTE`'s
  one-sign-change rule needs no special case, and a negative-radius active flank is legitimate
  (the user asked explicitly for this not to be forbidden). `contact.flanc_arcs` goes through it.
- `inversions(kr, ra)` ignores arcs whose sagitta is under `SAGITTA_DROITE` (5e-5 mm). A straight
  root transition (`qg > 0`) comes out of `krbgs` with a huge radius whose sign is pure noise;
  counting it invents two phantom inversions. On a `zz=-150` ring that made `approximer` throw
  away its best fit (0.0009 µm) in favour of a 1.6 µm one. `optim._swte` filters the same way, so
  the two agree; `optim.profil_valide` also uses `arcs.angle_depart` and is correct for a ring
  (the numba kernel `evaluer` around it is not).

Measured after that: `approximer` on rings (zz = -12, -40, -89, -150; narcs 2..7) gives 0.001 to
0.22 µm with the end tangent exact and always exactly one sign inversion, and `optim.profil_valide`
accepts every one of them. External profiles are byte-identical to before on the 12 wheels of the
regression set, and `optim.py test.Cev 4` still matches `contact.py` to 3e-5 with the same 57 % of
valid neighbours.

GUI: the `ZZ` field accepts -10000..10000 with `Champ.exclus=(-2, 2)` (a new field on `Champ`,
checked in `_lire`, since a `QIntValidator` cannot express a hole). The Création tab builds, draws
and converts a ring pair to arcs normally, and exports it to `.json`; `.Eno` is refused for a ring
(the Fortran format has no sign convention for Z). The Simulation and Optimisation tabs accept a
ring since stages 3b-4 below, so `MSG_COURONNE` is gone; what a ring is still refused is the rack
(`MSG_SCARVING` covers the skiving tool's own unfinished half).

**Contact, polyline engine (stage 3b), done.** Five things in `contact.Engrenement`, none of which
changes an external result:
- **signed centre distance**: `entraxe = (r1+r2) + copysign(1, r1+r2)*da`. Note `copysign(da, a)`
  would be wrong -- it takes `|da|` -- and briefly made `da = -30 µm` behave like `+30`.
- **`_tester_entraxe(profils, aa, interne)`**: for an internal pair the root is at the *large*
  radius, so the three tests become pinion tip inside ring root, ring tip clearing pinion root,
  and tip circles still crossing.
- **search sense `self.sens`**. `psi` turns the driven wheel *against its own motion* until it
  meets the driver. An external driven turns the other way from the driver, so that search is
  trigonometric (+1, as before); an internal one turns the same way, so it is clockwise (-1).
  With +1 on a ring the search finds the ring tip against the pinion root fillet -- a real
  contact, but on the wrong side. Driven positions are therefore counted in the search sense
  (`a2u`) and converted back by `self.sens` wherever a rotation is actually applied.
- **positive angular pitches**: `a1 = 2pi/|z1|`, `a2u = 2pi/|z2|`. The driver always advances
  trigonometrically whatever the sign of its tooth count, and in search-sense units the driven
  advances by `2pi/|z2|` per driver tooth for both signs.
- **which half-profile carries the active flank**. The driver pushes towards `-x` at the pitch
  point when its centre is below (`z > 0`) and towards `+x` when it is above (a ring), and a
  flank pushes away from its material. With the driver's own frame given a half-turn when
  `z1 < 0` (`self.delta1`, so its teeth point at the pitch point), this gives `miroir1 = -1`
  always and `miroir2 = -sign(z1*z2)`: on an internal mesh the driven wheel's active flank is
  the **un-mirrored** half-profile. `flanc_arcs(mobile, dos, miroir)` takes the side, and flips
  the radius sign when there is no mirror (the direction of travel is not reversed either).

Acceptance test, on involute pairs built by `evolve` and run through the polyline engine: **where
both flanks are on their involutes, the transmission ratio is the ideal one and the line of action
is exactly the pressure angle**. On 30/-90 teeth it is +0.333333 against an ideal +0.333333 with
the contact at -20.00 +- 0.02 degrees over 90 of the 150 positions, and it is positive (both wheels
turn the same way) as an internal mesh requires. A ring *driver* works too: -90/30 gives
+2.9960..+3.0129 for an ideal +3, an error of 4.3e-3 against the external reference's 1.9e-3.
The deviations that remain are geometric, not numerical, and show up as contacts off the line of
action: the ring's tip arc against the pinion flank at the start of the mesh (positions 0..52 on
30/-90: radius 5.961 on the ring, below its involute start at 5.967), and the last position of the
period, where pair 0 has disengaged and pair 1 has not taken over -- `contact2` is finite at only
1 of the 150 positions, i.e. the contact ratio of these test pairs is about 1. Both need looking at
before the numbers can be trusted end to end on a ring.

**Contact, arc engine (stage 3c), done -- the Simulation tab now runs on a ring.** Two changes,
neither of which touches an external result:
- `premier_contact_arcs(c1, c2, o2, sens)`: the rotation is now counted in the search sense, so
  `d = (sens*(aw - a02 + signe*ac)) % 2pi` and the rotation actually applied to `c2` is `sens*d`
  (used both to place the driven arc centres and to bring `th2` back before rotation). Its
  convex/concave tangency logic needed nothing: a ring's concave active flank simply lands in the
  "convex inside concave" case, which was already there.
- `flanc_arcs`: a mirror is a *reflection*, so it does not change convexity. An earlier attempt to
  flip the radius sign when there is no mirror (reasoning about the direction of travel) was wrong
  -- `arcs.convexite` has already removed the traversal dependence -- and it was what kept the arc
  engine from finding any contact on a ring.

Cross-checked against the polyline engine on the same geometry (30/-90 and -90/30, `me=1.2`,
8 arcs, mu=0.2): mean transmission ratio +0.33335 for an ideal +0.33333 (5.5e-5) and +2.99999 for
an ideal +3 (**2.4e-6**, against 9.9e-6 for the external 30/90 reference), and the **backlash is
identical between the two engines to the fourth decimal** (0.1908 and 0.5716 degrees) -- which also
exercises the `jeu` path and its back flanks. Hertz differs a lot between the two engines on a ring
(447 against 1203 on 30/-90) and jumps with the arc count, as it already does on external pairs
when the contact changes arc; an internal mesh is genuinely conformal (convex in concave), so its
Hertz pressure is far lower than the external pair of the same size (447 against 1283), but the
figure is not yet trustworthy enough to report.

**The driven position has to be unwrapped around a fixed reference** -- this was a real bug, and
the conclusion that used to sit here ("the contact ratio is the thing to watch, not the code") was
wrong. Reported by the user on `ko.json` (pinion z=12 driving a ring z=-76): "deux paires de dents
ne sont jamais en prise simultanément" whatever the tooth geometry. `premier_contact` returns a
rotation in `[0, 2pi[`; when the driven wheel's resting position sits near 0 it flips between ~0
and ~2pi from one driver position to the next, so `debut` sees a 2pi step instead of the crossing
of the two pairs and never finds the period. An internal mesh lands exactly there, because the
search runs the other way and the contact falls just the other side of zero. `Engrenement._reference`
now takes one reference value `u0` (the contact at phi=0) and `_derouler` maps every later value
into the 2pi-wide interval centred on it; the physical gap between neighbouring pairs is one driven
pitch, far smaller than pi, so the ordering between pairs is untouched. The numba kernel carries the
same `u0` through `_contact`/`_c`/`_ecart`/`_debut`/`_analyser` (plus a second one, `u0b`, for the
back flanks used by the backlash). One trap inside the fix: in the polyline engine the refined
search returns a *raw* angle while the coarse one is already unwrapped, so the `abs(d - d0) > 1e-3`
test rejected every refinement until `d` was unwrapped too.

What it was costing: `ko.json` went from refusing outright to torque variation 4.60 %, backlash
0.4326 degrees, Hertz 1132, with 63 of 150 positions carrying two pairs. Every ring pair that had
failed during development now runs, and the arc engine gives the ideal ratio to five decimals with
a healthy overlap: 30/-90 varC 2.20 % and 150/150 positions, 20/-60 4.62 % and 150/150, 12/-40
12.07 % and 123/150, -90/30 2.25 % and 134/150. The kernel agrees with `contact.Engrenement` on
`ko.json` at the three centre distances (Hertz to 4e-4, backlash to 3e-4).

Known limit, measured rather than guessed: the **polyline** engine is at its resolution limit on a
ring with a high contact ratio. The two pairs then agree to about 1e-5 over the whole period, which
is the size of `TOL` (3e-5), so `debut`'s bisection can stop anywhere in that band and lands on a
shifted period (30/-90: it starts half a pitch early and reports 1 overlapping position instead of
150). The arc engine, with `tol = 1e-9`, is unaffected, and it is the one the Simulation tab uses.
The external cases are unaffected (30/90: 103 of 150 positions).

**`da` is an opening in the *backlash* sense, not a centre-distance variation** (user-reported:
"+40 resserre les profils maintenant"). `evolvente.entraxe_signe` and `contact.Engrenement` add it
algebraically to the signed `rt1 + rt2`, so `da > 0` always loosens. On an external pair that is
also a larger centre distance; on an internal one it is a *smaller* one, because moving the pinion
away from the ring's centre drives it deeper into the teeth. Checked: on 30/90, da -20..+40 gives
backlash 0.055 -> 0.470 degrees with the centre distance 8.08 -> 8.14; on 30/-90, da 0 -> +20 gives
0.191 -> 0.327 degrees with the centre distance 4.05 -> 4.03. This is the convention the whole
optimiser rests on (-dA is the tight case, `jeu_min`; +dA the loose one, `jeu_max`). An earlier
`copysign(1, a)*da` "fix" made it the geometric centre distance and inverted the internal case.

**Optimisation kernel (stage 4), done.** The numba kernel got the same five conventions, threaded
as extra arguments so nothing is read from globals: `korzf(..., sg)` (the tip is the near
intersection with the tooth axis on a ring), `_flanc(..., cvs)` (physical convexity),
`_contact(..., sens)` (search direction, rotation applied is `sens*d`), `_c`/`_ecart`/`_debut`/
`_analyser`/`_interference` carrying `sens`, the driver's half-turn `d1` and the driven's tooth-axis
reference `axe2` (pi when its centre is on the far side of the pitch point, 0 when on the same
side), positive angular pitches, and the internal centre-distance tests. `Probleme._vers_kl` and
`_fermer` go through `arcs.angle_depart`.

One trap worth remembering: `korzf` used to keep the sign KRBGS had computed for the tip arc, but
`_construire` leaves the last point at (0,0) until KORZF fixes it, so that sign came from a
*fictitious* arc. It happened to be right on an external tooth and wrong on a ring, which made
`SWTE` reject every ring profile (code 5) while the same points passed when closed from a complete
chain. The tip arc is convex by construction, so its sign now comes from `sg`.

Validated against `contact.Engrenement` on the same geometry (6 arcs, mu=0.2, dA = -10/0/+10 µm):
ring driven 30/-90, Hertz 1420.6/1424.8/1426.4 against 1420.0/1424.1/1426.1 and backlash identical
to the fourth decimal (0.1217, 0.2589); ring driver -90/30, Hertz 184.5/184.6 identical and
backlash 0.3636/0.7757 against 0.3637/0.7758. The torque variation differs by 1-5 %, the same
spread as on the external 30/90 reference (7.372 against 7.120), because the kernel closes the
profile with KORZF and this ad-hoc comparison does not -- `python optim.py test.Cev 4`, which
closes both sides, still agrees to 3.0e-5. A real CMA-ES run on 30/-90 (2 starts, 3000 evaluations,
33 s) took the torque variation from 3.583 % to 2.452 % and the Hertz pressure from 1424 to 845
with no remaining constraint violation, the result confirmed independently by `contact.Engrenement`.

**Skiving tool (stage 5), done, crossed axes included.** Tool profile, tool tip, cutting with it
and the DXF export all work for both signs of the tooth count and for any crossing angle.

A ring can only be cut by a tool that passes through its central opening, so
the user asked for skiving to be generalised to the whole application, with **`z_outil = 0` meaning
the rack**, i.e. today's hobbing cutter (a rack is a wheel of *infinite* radius, but "no diameter"
is what one types to mean "fraise", and 0 cannot be a valid tooth count anyway).
`profil_outil(mobile, z_outil, ...)` delegates to `profil_fraise` when `z_outil == 0` -- checked
identical, point for point -- and otherwise returns a `ProfilOutil` (name, module, Z, the
half-profile root to tip, `zz_origine`; `.rayon` = `mo*zz/2`, `.echantillonner` as for a wheel).
Nothing in the hob path was touched: `fraise.py`'s diff against the last commit is additions only,
and the hob still matches the Fortran `Result/Validation_T.Frs` to 0.071/0.282 µm with the `.Frs`
round trip still bit-for-bit.

`_famille_outil` holds the parallel-axis kinematics (workpiece and tool rolling without slipping at
the signed centre distance `A = r1 + r2`, the tool turning by `-(r1/r2)u`). It is correct but no
longer on the path: the crossed-axis derivation below subsumes it, `S = 0` included.

**Why the envelope route was abandoned, measured.** Carrying the rack's recipe over -- sweep a
family, take the envelope (`_enveloppe_polaire`, the polar transpose of `_enveloppe`) -- gives
plausible-looking tool profiles that are wrong by 40 to 120 µm. The two obvious repairs are both
dead ends: a best-fit tangential shift only brings the deviation down from 165 to 87 µm, and the
conjugate of an involute comes out **not an involute at all** (fitted base radius 3 to 5 % off,
37-119 µm from the best-fitting involute, against 0.0000 µm for the workpiece itself). The shape is
wrong, not its phase, so neither hunting the offset nor re-anchoring on the pitch point could fix it.

**The method that works: the law of gearing** (`conjuguer`), point by point, no family and no
envelope. The general form of it, since the axes may be crossed, is the meshing equation
`N . v_rel = 0` with `v_rel = w1 (k1 x P) - w2 (k2 x (P - A))`. Workpiece axis `k1 = z`, tool axis
`k2 = (sin S, 0, cos S)` through `A = (0, a, 0)`, so the common perpendicular is the y axis and `a`
is the signed centre distance as everywhere else here. A **spur** workpiece -- the only kind in this
project -- has a cylindrical flank, so its normal has no axial component, and the equation collapses
to

    cr (1 - i cos S) = i [ a cos S Nx - sin S Ny Pz ],   cr = P x N,   i = w2/w1

The ratio comes from rolling on the common perpendicular, `i = -r1 / (r2 cos S)`, which is also what
sets the tool's pitch radius: `r2 = mo z2 / (2 cos S)` (`rayon_outil`), the classic crossed-helical
relation -- the tool is a helical gear of helix angle `S`, normal module `mo`, transverse module
`mo / cos S`.

At `S = 0` this is exactly `cos(u + nu) = -cr/r1`, the plane law of gearing (common normal through
the pitch point), kept in closed form, so the parallel-axis path is bit-for-bit unchanged. For
`S != 0` there is no closed form (the unknown appears as both `u` and `2u`), but the solution is
continuous in `S`: `_resoudre_contact` starts Newton from the exact `S = 0` root, which also settles
the branch without any vote.

**Why the crossing angle is not optional** (the user's point, and the reason the first version was
only half the job): with parallel axes there is no cutting speed at all -- that is gear *shaping*
(mortaisage), not skiving. `S` is what produces it.

What crossed axes change geometrically: the contact is a **point that travels along the workpiece
axis**. Every transverse profile point is contacted, each at its own axial height `Pz`, which is
precisely why the tool feeds axially. `conjuguer` therefore takes the section of the conjugate
surface by a plane of the tool, `h` being that plane's axial position (`h = 0`: through the common
perpendicular), and the conjugate point is `R(-i u) . (Px/cos S - h tan S, Py - a)` -- the `1/cos S`
being the same apparent stretch that turns `mo` into `mo / cos S`.

**One section is enough, and that is physics, not convenience**: in skiving the *cutting edge*
removes the material, not the flank, and with zero rake that edge is exactly the tool's transverse
section in the rake plane. A spur workpiece is invariant along its axis, so the axial feed changes
nothing in the transverse profile and the edge stays at a fixed height in the tool frame.

**What does change is the edge after resharpening**, and that is the one real limitation of skiving,
now measured rather than assumed. The exact conjugate surface of a spur workpiece on crossed axes is
**not quite a helicoid**, while a tool is ground as a helical gear: it can only be rigorously
conjugate in *one* plane. Measured on the exact-involute case, the section-to-section rotation does
follow the helix, `h tan S / r2`, to better than a percent (0.1291 against 0.1286 rad at S=10 deg,
h=1 mm) -- what drifts is the *shape*. `derive_affutage` reports it, rigid rotation removed (a
phase shift is free, the machine absorbs it), as a rate in µm per mm of regrind. It is close to
proportional to `h` (checked from 0.02 to 0.1 mm) and grows roughly as `S^2`:

| workpiece / tool | S=5 deg | 10 | 15 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|
| 89 / 30 | 0.2 | 0.5 | 0.7 | 1.7 | 3.3 | 5.8 |
| 30 / 30 | 0.4 | 0.7 | 1.4 | 3.3 | 6.5 | 11.7 |
| -90 / 30 | 0.3 | 1.1 | 3.1 | 6.0 | 10.3 | 15.9 |
| -60 / 20 | 0.5 | 2.0 | 5.5 | 10.6 | 17.5 | 26.4 |

(µm per mm of regrind). That is the number to pick `S` on: it is nil at `S = 0`, where the tool does
not cut. The measure itself is local -- beyond ~0.1 mm the two sections no longer overlap enough in
radius and `_rangerter` truncates them differently, so what breaks down is the measurement, not the
geometry.

Normals come from the arcs analytically (`_flanc_normales`: centre to point, signed by
`arcs.convexite`), because it is the normal, not the position, that sets the conjugate -- finite
differences are not good enough.

Two things had to be got right, and only measurement settled either:
- **which of the two roots**. `cos(u + nu) = c` gives `u = -nu +- arccos(c)`, and the two roots
  differ only by the sign of the contact normal's y component in the fixed frame,
  `sin(u + nu) = +-sqrt(1-c^2)`. The right one is always the **positive** one, because the
  workpiece's outward normal has to point at the tool, and the tool is on the +y side of the pitch
  point for both signs of `r1` (external wheel: the tool's centre is beyond the pitch point; ring:
  the outward normal points at the ring's own centre, which is where the tool is). So the branch is
  plain `arccos`, which returns [0, pi] -- no vote and no heuristic. Two heuristics were tried
  first and both are wrong: scoring the normal's direction towards the tool centre fails as soon as
  the tool is small (the centre is then close to the contact, so the direction swings; a Z=-90 ring
  with a Z=20 tool came out at 3.4-11.0 mm instead of 1.15-1.55), and taking the smaller rotation
  fails near the base circle of a small external wheel, where the physical rotation is the larger
  of the two (Z=30: 90-183 µm instead of 6-13);
- **the workpiece needs the half-turn its sign asks for**, in `profil_outil`, before conjugating:
  a half-profile is always stored tooth-on-+y, while the pitch point is at `(0, r1)`, hence at -y
  for a ring. It is the same `delta1` as `contact.Engrenement`. Without it the contact rotation
  comes out near pi, the two roots are told apart backwards and the profile lands at the wrong
  radius -- which is exactly what the ring was doing.

The conjugate of the workpiece's **root fillet** is meaningless (nothing touches there) yet the
equation still returns a solution, which overshoots the flank; `_rangerter` keeps only the longest
run that is monotonic in radius. The first version of all this also forced the conjugate onto a gap
axis and a tooth axis copied from `profil_fraise`; those anchors are rack-specific, the conjugate
is the tool's *active flank* and nothing else, and inventing its root and tip was worth 20 µm.

Validated, with the two points above in place: on an **exact** involute (normals taken
analytically) the conjugate is the theoretical one to **0.0000 µm**, i.e. machine precision, for
every combination tried -- workpiece Z = 7, 12, 30, 89, 150, -12, -40, -90, -150 and tool Z = 12,
30, 200. The reference radius it is compared against is the closed form: with `s1 = -r1 sin(alpha)`
and `s2 = r2 sin(alpha)` the signed distances from the pitch point to the two base-circle tangency
points, a workpiece radius `rho` meets the tool at radius `hypot(rb2, t - s2)` with
`t = s1 + sign(r1) sqrt(rho^2 - rb1^2)`; the two signs of `r1` are what make an internal pair put
both tangency points on the same side. The *outward* normal has the opposite sign on a ring, which
is precisely what the half-turn supplies.

On a 9-arc approximation of that same involute the conjugate is 0.03 to 13.5 µm off (0.03-0.8 on
the rings, 0.3-13.5 on the external wheels, worst on Z=30), and that deviation falls with the arc
count (Z=30: 23.2 / 16.8 / 13.5 µm at narcs 4 / 6 / 9). It is the arc chain's *normals*, not the
method: conjugation is driven by the normal, so a chain that is 0.1 µm off in position can be ten
times that in conjugate position. Worth remembering when reading a tool profile computed from a
coarse arc chain -- and a reason to keep using analytic arc normals rather than finite differences.

**GUI (asked for by the user): one tool per wheel, in the Fraise tab.** `Z outil` is a pair of
fields (menante / menée) right under the "Roue affichée" combo, **0 meaning the hobbing cutter**
(rack) and anything >= 3 a skiving tool, so a pair can mix the two -- pinion hobbed, ring skived,
which is the normal case since a hob cannot enter a ring. `_z_outils` refuses 1-2, refuses 0 on a
ring ("une fraise-mère ne peut pas y entrer") and refuses a tool at least as big as the ring it has
to fit inside. `calculer_fraise` then calls `profil_fraise` or `profil_outil` per wheel and reports
each tool's kind and whether its tooth came out pointed; `self.z_outils`/`self.sigmas` keep the choice, and
`_tailler_roue` dispatches to `tailler` or `tailler_outil` from it. A second row of fields holds
**Σ**, the crossing angle in degrees per wheel, defaulting to `SIGMA_DEFAUT` = 20 (a usual skiving
value) and accepted from 0 to 60; the panel reports each skiving tool's pitch radius, transverse
module, helix angle and **resharpening drift in µm/mm**, and says in red when Σ = 0 that the tool
does not cut. Σ is ignored for a hob. `_dessiner_fraise` draws a rack
as before and a skiving tool as a wheel. `_texte_taillage` says "module effectif" for a rack
(`module + 2 oy/zz`, `Fraesen.F90`'s O/M choice) and "module + déport" for a skiving tool, where the
module is fixed by the two tooth counts and `oy` is only a profile shift. `ProfilOutil.module` is an
alias of `mo` so both tool classes read alike. The one thing still refused for a skiving tool is the
**`.Frs` file** (`MSG_SCARVING`): that Fortran format has no field for a tool tooth count.

End to end in the GUI, hobbing or skiving each wheel as its sign requires, then re-simulating the
cut pair as point sequences: `validation.json` (89/12, both hobbed) 1.3458 -> 1.3581 % torque
variation, backlash 2.3818 -> 2.3815 degrees; `couronne.json` (ring -60 skived with a Z=20 tool,
pinion 12 hobbed) 4.581 -> 5.056 %, backlash 5.4294 -> 5.4284; `ko.json` (pinion 12 hobbed, ring
-76 skived with Z=16) backlash 0.4326 -> 0.4688. The torque figures of the two ring pairs are not
to be read too closely: in point mode they run on the polyline engine, which is documented above as
being at its resolution limit on a ring with a high contact ratio. The trustworthy number there is
the deviation of the cut profile from the reference, and it behaves as the tool size says it should
(on `ko.json`'s ring: 2.6 µm with a Z=40 tool, 7.2 with Z=16, 18.4 with Z=12, the two small ones
being pointed).

**Converting a cut wheel back to arcs can fail at a given arc count, and now says so.** Nothing
guarantees a swept profile is approximable by `n` tangent arcs: on `couronne.json`'s skived ring the
fit is 2.47 µm and valid at 4 arcs, 2.69 µm but SWTE-invalid at 6, and folds back on itself at 5
(44.70 µm) -- which is exactly the count `_proposer_n_arcs` suggests, since it copies the origin
wheel's. Until now that only showed up downstream, as the simulation refusing the profile
("deux paires de dents ne sont jamais en prise simultanément"), which is unreadable.
`utiliser_taillage_simulation` now runs `optim.profil_valide` on each re-fitted wheel and says in
the status bar which one folded and at what arc count, pointing at another arc count or at the
point-sequence path. The `reference=` layout does not help here: it only applies when the arc count
matches the origin wheel's, and at that count this profile is precisely the one that fails.

The tool frame that comes out of `conjuguer` has its **+y on the tool's gap axis**, for both signs:
`u = 0` is the workpiece in its stored orientation, i.e. a tooth of it on the line of centres, and
a workpiece tooth sits in a tool gap. The two flanks of a tool tooth are therefore mirror images
about that axis, which is what `_dessiner_fraise` draws, and the tool's own tooth axis is half a
pitch away.

**The tool's own tip, `_rangerter`.** The conjugate of the whole given flank also covers the
workpiece's root fillet, and that part runs past the neighbouring flank -- the hook that used to sit
on top of every drawn tool tooth. It is cut in two places, both physical rather than heuristic:
- **the tooth axis**, half a pitch (`pi/z`) from the gap axis: beyond it the curve belongs to the
  opposite flank of the same tooth. The tooth then comes out **pointed**, and that is the only tip
  the workpiece imposes; `ProfilOutil.tronque` says whether it happened, and the GUI reports it.
  When the conjugate stops short of that axis the tool has a free tip, which this module does not
  invent (`tailler_outil` closes it with a plain circular tip arc, and says so);
- **the fold**: the fillet's conjugate cusps back. The *longest* run that is monotone in radius is
  kept -- not the run starting at either end, since depending on the pair the fold sits at the head
  or the foot of the array and both ends can be inside it (measured: workpiece Z=30 with a Z=12
  tool gives 0.648..0.763 mm instead of 0.769..1.232 if one starts from an end).

A pointed tool is a real statement, not a defect: with `ma = me = 1`, `qe = 1.1`, a Z=30 external
workpiece truncates **every** tool from Z=20 to Z=200, i.e. no skiving tool of that family can
generate its whole flank.

**The reverse direction, `tailler_outil` (tool -> cut workpiece), done.** Same role as `tailler`
for the rack: the tool rolls on the workpiece, its flank sweeps a family (`_famille_piece`) and the
envelope (`_enveloppe_angle`, the polar transpose of `_enveloppe`, the angle replacing the abscissa
and the radius the height, same ascending-segments filter) is the generated flank. `_famille_piece`
is not a re-derived kinematic: it is `conjuguer`'s own map inverted algebraically,
`P = R(-u).(R(-j u) q + (0, A))`, checked against it to 9e-16 mm. `oy` is the profile shift, in mm,
exactly as in `tailler`.

Three things it must not do, each found by measurement:
- **do not force the end points**. `tailler` can impose its tangents (`-pi/z` at the root, 0 at the
  tip) because it first adds `RADABWP`'s two analytic points, on the gap axis and on the tooth axis.
  Here the envelope reaches neither -- the bottom of the gap and the workpiece's tip radius are
  blank radii that the tool does not generate -- and adding them cost 70 to 180 µm. The end tangents
  are measured on the envelope instead;
- **undo `_rangerter`'s half-turn** before building the family, and give the ring the half-turn back
  afterwards so the result is in the usual storage frame (tooth on +y). Forgetting the first put the
  whole family at radius 5.5-6.3 mm instead of 1.8-2.2;
- **the workpiece's root and tip radii come from the tool's own profile** (`|A| - sign(r1) R`),
  never from the min/max of the family -- the same rule as `tailler`, where ignoring it was worth
  400 µm.

Validated by the **round trip** workpiece -> tool -> cut workpiece, on the involute part (the root
fillet is a given arc that no tool reproduces: the tool cuts its own trochoid there, and the tip
radius is not generated at all): **0.035 to 0.17 µm** for workpiece Z = 12, 30, 89, -40, -60, -90,
-150 against tool Z = 20 to 200, pointed tools and free-tip tools alike, **and for crossing angles
0, 10, 20 and 30 degrees** (0.055-0.17 µm across that whole grid), flat in `n_positions`
(61 / 151 / 301) -- so what is left is the 9-arc approximation of the reference, not the sweep
resolution. Outside that band the deviation is confined to the fillet and to one or two end points.

The crossed-axis tool is also checked against theory directly: the conjugate of an exact involute
is an involute of base radius `r2 cos(alpha_t)` with `tan(alpha_t) = tan(alpha)/cos(S)`, the
crossed-helical relation. Fitted against theoretical base radius: 1.28582/1.28582 at S=10 deg,
1.33995/1.33966 at 20, 1.43873/1.43709 at 30 (workpiece Z=30, tool Z=20), with an involute residual
of 0.004-0.37 µm up to S=20 deg. The residual grows with S (1-2 µm at 30 deg, 2-9 at 45) and that
growth is the non-helicoid effect above, not a numerical one.

**The rack and the large-tool limit do agree** -- the 166 µm that used to sit here was a wrong frame
mapping, nothing else. The right one, found by trying the four sign pairings against a fitted shift:
`X_rack = -r2*beta`, `Y_rack = -(rho - r2)`, with `beta` the polar angle from the tool's gap axis and
**a fitted offset of 0.00000 mm** -- the two frames already coincide. The two profiles then agree to
**0.01-0.13 µm** over their common range on 89 and 12 teeth, the only larger value being the very
last point, where the tool's conjugate and the rack stop at different places. The three wrong
pairings give 127-141 µm, which is where the old figure came from. It does not move with the rack's
`n_positions` (61 to 1001), i.e. it is not a sweep-resolution effect either way.

Left to do on rings: nothing structural. If it is ever wanted, a file format for a skiving tool
(the `.Frs` layout has no field for a tool tooth count or a crossing angle), and a helical
*workpiece* -- every derivation above assumes a spur one, which is all this project has. Création,
Simulation, Optimisation, the Fraise tab and the DXF export all accept a ring;
`_remplir_arcs_sim`/`_remplir_variables` label the arcs through `arcs.convexite`, so a ring's active
flank reads "concave" and not "convexe".

The static "almost in mesh" drawing is now checked rather than assumed: on 12/40, 12/-40, -40/12,
30/-90 and -90/30, the two full outlines never cross near the pitch point and keep a minimum
clearance of 5.3 to 9.5 microns, the external pair's being 9.5 -- so the half-pitch offset places
the teeth correctly in all four internal configurations. The 5.3 micron case (-40/12) is the ring
tip against the pinion root fillet, a real and tight clearance on that geometry, not a drawing
artefact.

## Profile conventions (easy to get wrong when porting)

- A profile is a half pitch, with the origin at the wheel centre: `NK = NI//10` tip-arc points, the involute `NE = NI-2*NK`, and the root arc plus the 2-point straight transition (only if `QG>0`) filling the start. The last point lies on the tooth axis (`x=0`) and the first point is on the gap axis at `180°/ZZ`. `x` is negated at the end, so the flank lies on the `x>0` side; the other flank is its mirror.
- The involute points are clustered around the pitch circle: index `N0` is where it crosses the pitch circle, with different step sizes below and above it.
- Porting pitfalls: Fortran assignment to an `INTEGER` truncates (`N0`), arrays are 1-based, and in a one-line `IF(...) A; B` (Fortran, `;`) `B` always executes.
