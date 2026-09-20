# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Python port, in progress, of ENGCRE, a Fortran program by K. Trümpy that generates gear-tooth profiles (involute or cycloid) and analyses the meshing of a gear pair (kinematics, statics, Hertz pressure). The Fortran version works and is the reference for all the maths. Users are French-speaking; Fortran comments are in German.

Only the profile generation is ported so far: [evolvente.py](evolvente.py) is a port of `Fortran/Egc/Evolve.F90`. Everything else (`Cyclo`, `Krebax`, `Enganap`, plotting…) is still Fortran-only.

## Commands

```
python evolvente.py [file.Cev] [--ni 150] [--plot]   # default file: test.Cev -> Result/<name>.Egr
```

- Pure standard library (no numpy). `--plot` lazily imports matplotlib, which is not installed here and is untested.
- There is no test suite or build system. The git remote is `https://github.com/bernardflorex/eng_cre.git`; `Fortran/` is in `.gitignore`. To check a change to `evolve()`, compare with the Fortran output `Result/f.Egr`, which was generated from `Files/f.Cev` with **NI=550**: `evolve(roue, 550)` must match it to about 1e-7.
- The Fortran cannot be built or run here: it needs Windows and Intel Fortran (`ifort`, QuickWin graphics, `winapi_c.cpp`). `Fortran/*.Exe` and `libiomp5md.dll` are Windows binaries.
- `.vscode/launch.json` points to a `main.py` that does not exist; it is stale.

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

`.Cev` (involute, whitespace separated, first line is a header, CRLF possible):
```
line 2: MO ZZ NW QE QG MA ME      # driver wheel
line 3: MO ZZ NW QE QG MA ME      # driven wheel
then (optional for profile generation): E-mod / Poisson / thickness header + 2 lines for Hertz
```
`MO` module, `ZZ` tooth count, `NW` pressure angle in degrees, `QE` gap/tooth-width ratio at the pitch circle (≥ 1), `QG` length of the straight root transition in module units (0 = none), `MA`/`ME` extent of the involute below/above the pitch circle in module units (some files call them M1/M2). `.Cyl` replaces `NW` with `R1 R2` (rolling-circle radii).

`.Egr`: per wheel, 3 header lines (label and parameters; `MO ZZ 1 0`; `NI`) followed by `NI` lines `x y`, then a Hertz block. `evolvente.py` writes the same layout without the Hertz block and with a simpler label.

## Profile conventions (easy to get wrong when porting)

- A profile is a half pitch, with the origin at the wheel centre: `NK = NI//10` tip-arc points, the involute `NE = NI-2*NK`, and the root arc plus the 2-point straight transition (only if `QG>0`) filling the start. The last point lies on the tooth axis (`x=0`) and the first point is on the gap axis at `180°/ZZ`. `x` is negated at the end, so the tooth is on the `x<0` side.
- The involute points are clustered around the pitch circle: index `N0` is where it crosses the pitch circle, with different step sizes below and above it.
- Porting pitfalls: Fortran assignment to an `INTEGER` truncates (`N0`), arrays are 1-based, and in a one-line `IF(...) A; B` (Fortran, `;`) `B` always executes.
