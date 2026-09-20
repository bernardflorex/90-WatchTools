# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Python port, in progress, of ENGCRE, a Fortran program by K. Trümpy that generates gear-tooth profiles (involute or cycloid) and analyses the meshing of a gear pair (kinematics, statics, Hertz pressure). The Fortran version works and is the reference for all the maths. Users are French-speaking; Fortran comments are in German.

Ported so far: [evolvente.py](evolvente.py) (port of `Fortran/Egc/Evolve.F90`, file reading/writing, drawing of the gears), [contact.py](contact.py) (numpy port of the contact search of `Kontap.F90`/`Enganap.F90` and of the transmission/torque/efficiency formulas of `Enganap.F90`/`Draweng.F90`) [arcs.py](arcs.py) (port of `Krebax.F90`/`Lagrmo.F90`: profile approximated by tangent circular arcs, and the `.Eno` export) and [main.py](main.py), a PySide6 window with the parameters on the left and the matplotlib drawing on the right. Everything else (`Cyclo`, the circular-arc contact analysis `Enganac`, Hertz pressure…) is still Fortran-only.

## Commands

```
python main.py [file.json|file.Cev]      # GUI (default file: test.Cev, imported)
python evolvente.py [file] [--ni 150] [--plot] [--png out.png] [--dents 3]   # CLI, also writes Result/<name>.Egr
```

- Interpreter: the system `python3` has neither pip nor matplotlib/PySide6. Use the venv `~/DEV/virtual_env/.eng/bin/python` (matplotlib, PySide6, numpy, scipy…; no tkinter, so the matplotlib GUI backend is Qt). Headless GUI test: `QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg`, then `window.grab().save("x.png")`.
- Careful: the CLI overwrites `Result/<name>.Egr`, and `Result/f.Egr` is the tracked **Fortran** reference (restore it with `git checkout -- Result/f.Egr`). Run tests from a scratch directory or pass a file whose name is not `f`.
- `evolvente.py` itself is pure standard library, matplotlib is imported only for drawing; `contact.py` needs numpy. `dessiner_engrenage(ax, …)` draws `--dents` teeth of each wheel almost in mesh (driven wheel turned by 180° minus half a pitch, layout as in `Drawkpo.F90`); the profiles are half-teeth, so an odd tooth count is not symmetric for the driven wheel.
- `main.py` builds the left panel from the `CHAMPS` table (one row per `.Cev` parameter; a new user variable = a new `Champ` plus an attribute on `evolvente.Roue`). Edits do not recompute until the "Reconstruire" button (F5), which runs the whole simulation (about 0.5 s). A slider walks through the 150 positions of one tooth-pair passage; the drawing is redrawn per position (fixed margins, `constrained` layout is too slow). The three curves (W transmission, C torque, N efficiency, as in the Fortran layout `layout de sortie des simulations.png`) are drawn over the gears at the x of the contact point, in % of the ideal value, 100 % at the pitch point (`evolvente.y_pourcent`: 5 % = 0.6 module); the panel also shows min/mean/max/variation per parameter. Files are saved as `.json` only; the old `.Cev` format is import-only (`lire_cev`), so the Fortran program cannot read what the GUI saves. The Hertz data (E-modulus, Poisson, thickness) are kept in the JSON (`hertz` key per wheel) but not editable yet.
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

`.Egr`: per wheel, 3 header lines (label and parameters; `MO ZZ 1 0`; `NI`) followed by `NI` lines `x y`, then a Hertz block. `evolvente.py` writes the same layout without the Hertz block and with a simpler label.

## Contact and meshing simulation (`contact.py`)

Drawing frame: pitch point at (0,0), driver (wheel 1) centre at (0,-rt1) turning counter-clockwise, driven (wheel 2) centre at (0,rt2), angles counter-clockwise from +y. The driving flank of both wheels is the mirrored profile (`x<0`).

- `premier_contact` (KONTAP): rotates the driven flank from a far reference position until it touches the driver flank; returns the rotation, the contact point and the common tangent. It works on polylines, so the contact point is only accurate to about the square root of the chord sagitta: `Engrenement._contact` first searches with `NA=100` points, then repeats the search on `2*FEN` fine points (`NF=3000`) around the contact. Without this refinement the transmission ratio is noisy at the 1e-3 level.
- `Engrenement.debut` finds the driver angle where the next tooth pair takes over from the previous one (FW in `Enganap.F90`), i.e. where two pairs touch at once. Pair `k` (driver tooth `k`, driven tooth `-k`) puts the driven axis at `psi(alpha1 + k*A1) + k*A2`; the driven wheel sits on the smallest value. Both pairs must really touch, otherwise `ErreurContact` (contact ratio < 1, e.g. `Files/Evolve3.Cev` with MA=0).
- `Engrenement.simuler` walks `N_POS=150` driver positions from that angle over one pitch `A1`, and returns arrays for the contact point (the line of contact), the second contact during the two-pair overlap, and the three quantities of `Draweng.F90`: transmission ratio `w2/w1` (negative: opposite senses), torque ratio `M2/M1` and efficiency, their product. Forces come from Coulomb friction `mu` opposing the sliding velocity, the driver torque being 1 N·mm (`efforts`); `Simulation.normale` is the normal force for that torque, scaled by the user's `mt` in the GUI. Pair 0 is kept in the overlap unless pair 1 is lower by more than `SEUIL_PAIRE`.
- Centre distance: `Engrenement(roues, profils, da)` (`da` in mm) places the wheels at `rt1 + rt2 + da` with working pitch circles proportional to the tooth counts, tangent at (0,0); `_tester_entraxe` ports the FI=1/FI=2 tests of `Enganap.F90` (tip against root / tip circles no longer crossing). `Engrenement.jeu` gives the driven-wheel backlash: by symmetry about the line of centres, the lowest driven position with the driver at `phi` is the mirror of the highest one at `-phi`, so backlash = `theta_max(phi) + theta_max(-phi)` folded modulo one driven pitch. It agrees with a brute-force outline-crossing search to about 5e-3 degrees, and changes by the expected `2*dA*tan(alpha)/r2` (30 µm gives +1.57° on `test.Cev`).
- Checks done by hand: with `mu=0` the efficiency is exactly 1; the line of contact is straight at the pressure angle; contact point and tangent are consistent with both base circles to about 1e-8 (driver) and 1e-5 (driven); efficiency drops with the distance to the pitch point as the first-order formula `1 - mu*(1+Z1/Z2)*|x|/rb1` predicts (within 1%); on 6 files of `Files/` the full tooth outlines do not cross at the start position. Not compared with the Fortran numbers (the Fortran cannot be run here).

## Circular arcs and `.Eno` export (`arcs.py`)

`approximer(profil, zz, qg, narcs)` = `KREBAX` + `LAGRMO` + `KRBGS`: transition points start at the root-arc end, the straight root part (if `qg>0`), `narcs` evenly spread involute points and the tip-arc end; `lagrmo` moves the interior points (Newton with a diagonal Hessian and a Lagrange multiplier for the final slope 0) to minimise the squared distances of the 150 profile points; the first tangent angle is `-pi/Z`. The arcs are fully described by the transition points, which is why the `.Eno` file only holds those. Checked: arcs pass through their end points, tangential junctions to 1e-6 degrees, first and last points are the profile ends, error at most about 1 µm with 3 arcs and 0.1 µm with 6 (`test.Cev`). `ecrire_eno` reproduces the Fortran layout (`F10.6,2I10` header, `2F12.7,I6` points with flag 0 at both ends, Hertz lines, the fixed `-dA +dA GW ZV` line, CRLF) and was read back with the reader of `~/DEV/Python/001-Eng_open/fileIO.py` (`read_eno`). Its Hertz lines are assigned to master/slave in the opposite order to the Fortran (first line = driver), which is that program's convention, not ours. Files without Hertz data get the sample values `HERTZ_DEFAUT` of `main.py` in the panel.

## Profile conventions (easy to get wrong when porting)

- A profile is a half pitch, with the origin at the wheel centre: `NK = NI//10` tip-arc points, the involute `NE = NI-2*NK`, and the root arc plus the 2-point straight transition (only if `QG>0`) filling the start. The last point lies on the tooth axis (`x=0`) and the first point is on the gap axis at `180°/ZZ`. `x` is negated at the end, so the flank lies on the `x>0` side; the other flank is its mirror.
- The involute points are clustered around the pitch circle: index `N0` is where it crosses the pitch circle, with different step sizes below and above it.
- Porting pitfalls: Fortran assignment to an `INTEGER` truncates (`N0`), arrays are 1-based, and in a one-line `IF(...) A; B` (Fortran, `;`) `B` always executes.
