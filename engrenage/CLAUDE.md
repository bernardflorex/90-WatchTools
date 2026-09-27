# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Python port, in progress, of ENGCRE, a Fortran program by K. Trümpy that generates gear-tooth profiles (involute or cycloid) and analyses the meshing of a gear pair (kinematics, statics, Hertz pressure). The Fortran version works and is the reference for all the maths. Users are French-speaking; Fortran comments are in German.

Ported so far: [evolvente.py](evolvente.py) (port of `Fortran/Egc/Evolve.F90`, file reading/writing, drawing of the gears), [contact.py](contact.py) (numpy port of the contact search of `Kontap.F90`/`Enganap.F90` and of the transmission/torque/efficiency formulas of `Enganap.F90`/`Draweng.F90`) [arcs.py](arcs.py) (port of `Krebax.F90`/`Lagrmo.F90`: profile approximated by tangent circular arcs; reading and writing of arcs files in `.json`, `.Eng` and `.Eno`) [optim.py](optim.py) (numba port of `KRBGS`/`KORZF`/`SWTE`, `Kontac`, `Enganac` and `Resulto` of `EngOpt/`, objective function and constraints), [recherche.py](recherche.py) (CMA-ES driver with restarts, repair by the back half-profile, list of distinct optima), [fraise.py](fraise.py) (port of `Abwfrac.F90`/`Abwfrap.F90`/`Radabwp.F90`/`Krvabc.F90`/`Krvabp.F90`: hobbing-cutter profile conjugate to a wheel, simulated hobbing of a target wheel, deviation from the reference profile) and [main.py](main.py), a PySide6 window: four tabs on the left ("Création": parameters of the involute gears, arcs and export; "Simulation arcs": meshing simulation of circular-arc profiles, from the Création tab or from an opened file; "Optimisation": parameters, constraints, run and list of optima; "Fraise": hobbing-cutter profile, `.Frs` file, simulated cutting), the matplotlib drawing on the right. Everything else (`Cyclo`, the circular-arc contact analysis `Enganac`, Hertz pressure…) is still Fortran-only. `rapport.py` writes a `.docx` simulation report (File > "Exporter un rapport…" in the GUI) from `contact.py`'s results, independent of the optimiser/hobbing modules. `dxf.py` (using `ezdxf`) exports mobile and hobbed-wheel profiles as real arc/line DXF entities and the hob profile as a DXF polyline.

## Commands

```
python main.py [file.json|file.Cev]      # GUI (default file: test.Cev, imported)
python evolvente.py [file] [--ni 150] [--plot] [--png out.png] [--dents 3]   # CLI, also writes Result/<name>.Egr
```

- Interpreter: the system `python3` has neither pip nor matplotlib/PySide6. Use the venv `~/DEV/virtual_env/.eng/bin/python` (matplotlib, PySide6, numpy, scipy, numba, cma…; no tkinter, so the matplotlib GUI backend is Qt). Headless GUI test: `QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg`, then `window.grab().save("x.png")`.
- Careful: the CLI overwrites `Result/<name>.Egr`, and `Result/f.Egr` is the tracked **Fortran** reference (restore it with `git checkout -- Result/f.Egr`). Run tests from a scratch directory or pass a file whose name is not `f`.
- `evolvente.py` itself is pure standard library, matplotlib is imported only for drawing; `contact.py` needs numpy. `dessiner_engrenage(ax, …)` draws `--dents` teeth of each wheel almost in mesh (driven wheel turned by 180° minus half a pitch, layout as in `Drawkpo.F90`); the profiles are half-teeth, so an odd tooth count is not symmetric for the driven wheel.
- The Qt font is enlarged in `main()` by `AUGMENTATION_POLICE` points (`agrandir_police`, currently +2); long unwrapped labels or checkbox texts widen the left panel and cause horizontal scrolling, so keep texts short or `_etiq` (word-wrapped). `main.py` state is split per tab: Création (`roues`, `profils`, `arcs_roues`) and Simulation (`mobiles` = two `arcs.Mobile`, `sim`, `sim_roues`, `sim_profils`), both drawing on the same axes (`_vue_ok` keeps the user's zoom between redraws; it is reset on tab change). F5 is one menu action (`actualiser`) that rebuilds or simulates depending on the tab, so the buttons carry no shortcut (two identical shortcuts would be ambiguous). The Simulation tab runs `contact.Engrenement(..., mobiles=self.mobiles)`, the analytic arc-to-arc engine (no profile points involved; the sampled arcs are only passed for the centre-distance test and for drawing). `File > Open` dispatches by extension (and by content for `.json`: `_json_est_arcs`). `main.py` builds the left panel from the `CHAMPS` table (one row per `.Cev` parameter; a new user variable = a new `Champ` plus an attribute on `evolvente.Roue`). Edits do not recompute until the "Reconstruire" button (F5), which runs the whole simulation (about 0.5 s). A slider walks through the 150 positions of one tooth-pair passage; the drawing is redrawn per position (fixed margins, `constrained` layout is too slow). The three curves (W transmission, C torque, N efficiency, as in the Fortran layout `layout de sortie des simulations.png`) are drawn over the gears at the x of the contact point, in % of the ideal value, 100 % at the pitch point (`evolvente.y_pourcent`: 5 % = 0.6 module); the panel also shows min/mean/max/variation per parameter. Files are saved as `.json` only; the old `.Cev` format is import-only (`lire_cev`), so the Fortran program cannot read what the GUI saves. The Hertz data (E-modulus, Poisson, thickness) are kept in the JSON (`hertz` key per wheel) but not editable yet.
- Bug found by the user: closing the app (or opening/reloading a file) always asked to save, even with no edits. Cause: `remplir_optimisation()` set the two objective combo boxes with `setCurrentIndex(...)` without `blockSignals`, so their `currentIndexChanged` (wired to `_marquer_modifie`) fired on every programmatic fill (startup, every load) and set `self.modifie = True` before the user touched anything. `case_dissym` already did this correctly; the fix wraps the two combo updates the same way. `QLineEdit` fields were never affected (`textEdited` only fires on user input, unlike `setText`).
- There is no test suite or build system. The git remote is `https://github.com/bernardflorex/eng_cre.git`; `Fortran/` is in `.gitignore`. To check a change to `evolve()`, compare with the Fortran output `Result/f.Egr`, which was generated from `Files/f.Cev` with **NI=550**: `evolve(roue, 550)` must match it to about 1e-7.
- The Fortran cannot be built or run here: it needs Windows and Intel Fortran (`ifort`, QuickWin graphics, `winapi_c.cpp`). `Fortran/*.Exe` and `libiomp5md.dll` are Windows binaries.
- `.vscode/launch.json` runs `main.py` with the interpreter selected in VSCode (it must be the `.eng` venv).

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
- Checks done by hand: with `mu=0` the efficiency is exactly 1; the line of contact is straight at the pressure angle; contact point and tangent are consistent with both base circles to about 1e-8 (driver) and 1e-5 (driven); efficiency drops with the distance to the pitch point as the first-order formula `1 - mu*(1+Z1/Z2)*|x|/rb1` predicts (within 1%); on 6 files of `Files/` the full tooth outlines do not cross at the start position. Not compared with the Fortran numbers (the Fortran cannot be run here).

## Optimisation (`optim.py`, `recherche.py`, `EngOpt/`)

`EngOpt/` holds the Fortran sources of the optimiser ENGOPT (`Optimal`, `Bidiopt`, `Variopt`, `Resulto`, `Enganac`, `Kontac`, `Krbgs`, `Manipu`; not run here). What it does: the objective `Resulto` is the square root of the sum, over the centre distances -dA, 0, +dA (line `-dA +dA GW ZV GL` of the `.Eno` file), of `(100 varC)^2` (ZV=1) or `(100 (1-mean efficiency))^2` (ZV=2), plus `(Hertz/GW)^2` and `((1-SymL) GL)^2`; ZV=3 minimises the backlash with minimum `|GW|`. A profile is invalid (1D33) when `SWTE` fails (exactly one sign change of the arc radii, tip radius >= `RK` x module), when the centre distance is inconsistent or when there is no mesh period. The variables are the interior transition points (x, y) with flag `IL=1`; the tip point is recomputed by `KORZF` (last arc centred on the tooth axis). `BIDIOPT`: per point a (2MN+1)^2 grid, weighted sum of the best moves, step reduced by 0.618, multi-start only through different initial steps; `VARIOPT`: full factorial (2MN+1)^(NP-2) per axis and wheel.

`optim.py` is the fast Python counterpart, checked against `contact.Engrenement` (`python optim.py [file] [narcs]`: same torque variation, Hertz pressure and backlash to 3e-5; 1000-3000 evaluations/s on 16 cores). Numba compiles on first use (a few seconds, cached in `__pycache__`); it and `cma` are in the `.eng` venv.
- `Probleme(roues, mobiles, mu, da_neg, da_pos, libres, objectif, dos_libre, avant_libre)`: decodes a vector of variables into arcs (`decoder`, `points`, `mobiles`), evaluates N vectors at once (`evaluer` -> array (N, 3 centre distances, `NCOL` columns `VARC…ECART_TETE`), `valeur`, `Probleme(x)`), reports one profile in detail (`rapport`). Variables = curvature and length of each free arc of the first n-2 arcs (`ECHELLE`; first point fixed, last arc closed like `KORZF`); `dos_libre` adds the same variables for the back half-profile (relative to the driving one), `avant_libre=False` keeps the driving profile and moves only the back one (repair).
- `Objectif`: `J = f1/f1₀ + poids f2/f2₀` (`couple`: torque variation, `jeu`: max backlash at ±dA, `rendement`: efficiency loss, `hertz`: max Hertz; quadratic mean over the 3 centre distances; `f₀` = start profile, so J = 1 at the start) + penalties: backlash >= `jeu_min` at -dA (anti-blocking) and <= `jeu_max` at +dA, tip-arc radius >= `r_tete_min` (mm, hob life; the Fortran `RK` is a relative test on x of the tip-arc start, this one is the true radius of the last arc, both flanks), no interference, same tip height on both flanks (<= 1 micron). Invalid profiles get `PENALITE`.
- Interference (`_interference`): arc-to-arc crossings between the complete outlines (both flanks, tip and root) of the teeth of both wheels near the line of centres, at 6 positions of the period and every centre distance, driven wheel placed by the working contact. Tangency (legitimate contact) is excluded (`EPS_INTERF` = 1 micron of overlap). Checked: 99-100 % of the profiles with negative backlash are flagged, 0 % of those with backlash > 0.3 degree, 0 on the reference profile. It costs about 1 ms per evaluation. Not checked: positions between the 6 sampled ones, tool (hob) interference / undercut, tooth thickness.
- Backlash uses the back flanks: `theta_max_avant(phi) + theta_max_arrière(-phi)` folded modulo one pitch (`contact.Engrenement.jeu` does the same with `dos=True`), identical to the numba kernel to 1e-3 degree with asymmetric teeth.
- `recherche.rechercher(...)`: several starts (the current profile, then random valid profiles ever further, up to 8 x `etendue_um`), a CMA-ES per start (population doubled every second start), profiles that still violate a constraint are repaired by freeing the back half-profiles of both wheels (driving half-profile unchanged, small regulariser towards symmetry, stops as soon as no violation is left), then the list of distinct optima (`separation_um` = 3 microns, valid ones first). `Optimum`: J, x, `mobiles` (with `points_dos` if asymmetric), per-centre-distance report, remaining violations. In practice the front search alone satisfied every constraint tried (backlash, tip radius): the repair only triggers when it cannot, and was tested on its own (raising the backlash from 0.93 to 1.5 degree with a 7 micron change of the back flank, front unchanged) and with a starved front budget.
- GUI (tab "Optimisation", menu Run > Optimiser F6 / Arrêter, buttons next to "Évaluer les critères et le temps"): search in a `QThread` (`Travail`), progress bar, table of optima (double-click = show in the Simulation tab, "Exporter…" = `.json`, `.Eno` refused for asymmetric teeth). Free arcs (not points) are chosen per wheel. The evaluation button shows the criteria at -ΔA/0/+ΔA, the constraints, the evaluation rate and the estimated duration. Project keys (`Projet.optimisation`): see `evolvente.OPTIMISATION_DEFAUT`; `jeu_max_admis` 0 = backlash of the start profile at +ΔA (+0.02).

Design decisions, from measurements on `test.Cev` (89/12 teeth, 4 arcs, sample Hertz data):
- Variables = curvature and length of each arc, NOT the (x, y) of the points: the valid neighbourhood is 10-30 microns instead of 1-3 (moving ONE point by 1 micron rotates the following tangents by 3-12 degrees), and CMA-ES cannot even start on (x, y).
- Modify both wheels together: at equal budget (12000 evaluations) joint 5.6, driven only 7.1, driver only 8.7-9.3, alternating driver/driven 10.5-12.2 (the valley is diagonal because conjugacy is relative).
- CMA-ES reaches its floor in 10-20 k evaluations; differential evolution gives twice the objective at the same budget; an emulated BIDIOPT is far behind (4.8 vs 0.8 after 61 k vs 10 k evaluations, varC only). Several seeds end within 3 %. With torque variation alone the backlash doubled (2.7 -> 5.3 degrees): the backlash bounds are needed. Typical result: torque variation 6-13 % -> 1-2 %, Hertz 646 -> about 400 N/mm^2 at +dA.
- The optimum found on -30/0/+30 microns stays good at every intermediate centre distance and degrades outside the range (backlash negative at -45 microns); ENGOPT prints statistics on 9 centre distances (-40..+40 step 10), check final solutions the same way. Hertz pressure jumps with the centre distance (arc change).
- Other basins exist (starting from involutes with NW 14.5-25 degrees all reach about the same objective except 14.5, better by 10 %); a start that is invalid stops CMA-ES, hence the random valid starts. More arcs (6, 8) do not lower the floor and need more evaluations.
- Bug found by the user: `arcs.approximer` (KREBAX/LAGRMO) occasionally produced, for a specific number of arcs, a transition-point sequence that folds back on itself (non-monotonic), which fails `SWTE` (code 5, message ends with "essayez un autre nombre d'arcs dans l'onglet Création", which reliably worked around it: e.g. `Files/Evolve109-11.Cev` failed only at narcs=4, not at 3 or 5). **The cause was later found and fixed in `approximer` itself** (index-based breakpoints, see the arcs section); that file now gives 0.37 µm at narcs=4 instead of 92.6 µm, so the workaround should no longer be needed. The robustness added below is still worth keeping, since nothing guarantees an arbitrary profile is approximable at an arbitrary arc count. Before that, `Probleme.__init__` unconditionally required this reference profile valid (`calibrer=True`) to compute `f0`, so the whole multi-start search aborted on construction, trying none of the configured départs; `main.Fenetre._jeu_max_auto` (used whenever "jeu max" is left at its "auto" default) had the same unconditional requirement, called *before* `recherche.rechercher` even started, so it failed the search even for objectives that do not need `f0`. Fixed: `Probleme(..., calibrer=False)` defers `f0`; `rechercher` tries the profile as given first, and if invalid searches nearby (`_tirage`, escalating the search radius 1x/10x/100x `etendue_um`) for a valid substitute to use as `f0` reference and as départ 1, surfacing the substitution (`stats["depart_invalide"]`, shown as a warning in the GUI) instead of raising; the jeu_max "auto" resolution moved from `Travail.run()` into `rechercher`, after that same reference is established. `evaluer_criteres` (button "Évaluer les critères et le temps") is intentionally left strict: it must describe the exact loaded profile, not a substitute.
- Not done: tool (hob) interference, `PZ` > 2 tooth pairs, the `GB` blocking option (replaced by `jeu_min`), identical driver/driven (`LI`), reading `.Eno` optimisation parameters beyond `-dA +dA GW`.

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
optimisation, then a skiving tool). **Only the first two stages exist so far**: `evolvente.evolve`
and the drawing. `contact.py`, `optim.py`, `arcs.py`, `fraise.py`, `dxf.py` and the GUI have *not*
been adapted and will give wrong answers on a negative tooth count. The ring may be either wheel
(the user asked for both roles), and the skiving tool is to be done parallel-axis first, crossed
axes afterwards.

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
- the drawing places the second wheel's teeth at 180 degrees when it is external and at 0 when it is
  a ring (its centre is then on the same side), minus half a pitch in both cases.

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

Not done yet, and known: the static "almost in mesh" drawing of an internal pair shows the two
outlines crossing slightly; that offset is a half-pitch guess, not a solved contact, so it needs
checking once `contact.py` handles rings.

## Profile conventions (easy to get wrong when porting)

- A profile is a half pitch, with the origin at the wheel centre: `NK = NI//10` tip-arc points, the involute `NE = NI-2*NK`, and the root arc plus the 2-point straight transition (only if `QG>0`) filling the start. The last point lies on the tooth axis (`x=0`) and the first point is on the gap axis at `180°/ZZ`. `x` is negated at the end, so the flank lies on the `x>0` side; the other flank is its mirror.
- The involute points are clustered around the pitch circle: index `N0` is where it crosses the pitch circle, with different step sizes below and above it.
- Porting pitfalls: Fortran assignment to an `INTEGER` truncates (`N0`), arrays are 1-based, and in a one-line `IF(...) A; B` (Fortran, `;`) `B` always executes.
