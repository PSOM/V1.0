# Changes on branch `fix-fortran-bugs` (October 2026)

This branch fixes bugs in the Fortran model that were found while validating the Julia port of PSOM
(`https://github.com/PSOM/Julia_Port`, branch `fix-and-validate`). The Julia port was checked line by
line, and routine by routine, against this Fortran. Wherever the two disagreed and the Fortran turned out
to be wrong, the Fortran was fixed here, and the same fix was applied to the Fortran copies inside
`Julia_Port/code/`. One commit per fix; the commit messages give the details.

Paths are relative to `code/`. Every experiment copy of an affected routine was fixed. `doc/BENCHMARKS/`
(a historical snapshot) and `trash/` were left alone.

**Summary of the effect on model output**

| Change | Serial runs | OpenMP runs |
|--------|-------------|-------------|
| F1 `linerelax` wrap | Changes runs where `NI/2` > the coarsest multigrid level's `NI` (e.g. `waves`), at solver-tolerance level | same |
| F2 OpenMP reductions | none | fixes a race in the convergence test |
| F3, F4, F6, F7 | none (out-of-bounds or uninitialised reads that didn't affect results in the tested cases; F7 is a build fix) | none |
| **F5 `staticsigma`** | **Non-hydrostatic runs: the initial state changes slightly** (wiggle: KE at step 0 by 4.5e-7 relative, fields by up to 3e-3) | Same fix; it removes a dependence on uninitialised stack memory |
| D1 waves sponge | **waves: different near the north/south sponges** (now the intended smooth taper) | same |
| D2 `mldepth` namelist | none for the prognostic fields; waves N² diagnostics are now finite | same |
| **D3 multi-colour solvers** | **All runs change at solver-tolerance level** (wiggle KE after 5 steps: 1.7e-9 relative) | **Now bit-identical to the serial run for any thread count** |

---

## 1. Bug fixes

| # | Commit | Fix | Copies fixed | Effect on results |
|---|--------|-----|--------------|-------------------|
| F1 | `c32c53d` | **`linerelax_periodicew.f90`: periodic wrap of the line start, `i = mod(i-1,nxm)+1`.** The old `if (i.gt.nxm) i = i - nxm` only wraps once. On coarse multigrid levels `istart = NI/2` can exceed that level's `nxm`, so `i` stayed out of range and the routine read and wrote outside that level's part of the 1-D pressure array `p`. | model, model_build, wiggle, waves, NA, NP_summer (and `NP_summer/src/with_OMP`), NP_winter, shelfbreak, exp_IRENE_withML_1wiggle_cooled, freewaves (11 files) | None when `istart ≤ nxm` (e.g. all random-start variants such as wiggle). Changes `waves` at the solver-tolerance level. |
| F2 | `8806901` | **OpenMP: `REDUCTION(max:maxres)`** in `resid.f90` and in the residual loops of `hsolve.f90` and `mprove.f90`. `maxres` was `SHARED` and updated by every thread (a data race), so the convergence test could see a wrong residual. | model, model_build, NP_summer/with_OMP, NP_winter | None in serial builds (`parallel` off). |
| F3 | `f7f59ec` | **`rpevalgrad_Song.f90`: k=0 entries at the right index.** `Jx2(j,0)` → `Jx2(i,0)` (out of bounds when NJ > NI). In the k=0 branch, `dzsig(j,1)` and `dzsig_x(i,1)` → `(.,0)`; the k=0 entries were otherwise read uninitialised. | model, leewaves, freewaves | None (the k=0 values are not used downstream). |
| F4 | `f8ff579` | **`diag_n2budget.f90`: `avezf(0:NK)`**, with `avezf(0)` = mean bottom-face depth. `avezf(k-1)` is read at k=1, which was out of bounds. | model, model_build, waves | None in the tested cases. |
| F5 | `025e688` | **`staticsigma.f90`: the `gi3`/`gj3` loops stop at k=NK-1.** The local `g13`/`g23` are only set for k=0..NK-1 in this routine, but the loops ran to NK. That overwrote the top-layer values that `sigma` (called just before, in `main`) had computed, with 0 in a normal build or with garbage in a recursive/OpenMP build. | model, model_build, leewaves, freewaves | **Changes the initial state of non-hydrostatic runs slightly** (wiggle: KE(0) by 4.5e-7 rel., fields up to 3e-3). Hydrostatic runs and waves unchanged. **This is also why OpenMP builds disagreed with serial builds by ~3e-3 even on 1 thread.** |
| F6 | `2a25625` | **waves `diag_n2budget.f90`: periodic halo uses `i2=NI`** (was `NK`, an out-of-bounds read of `pv1`). | waves | Diagnostic halo only. |
| F7 | `81cb49c` | **waves `linerelax_periodicew.f90`: `#ifdef parallel` around `OMP_GET_NUM_THREADS`**, as in the model version. Without it, waves did not link without OpenMP. | waves | None (build fix). |

How they were found:
- F1, F3, F4, F6, F7: by comparing the Fortran with the Julia port routine by routine (see section 4).
- F5: an OpenMP build run on **one** thread differed from the serial build by 2.7e-3. `-fopenmp` implies
  `-frecursive`, which puts local arrays on the stack instead of in static memory. A build with
  `-frecursive -finit-real=snan -ffpe-trap=invalid` then stopped on the uninitialised read in `staticsigma`.
- F2: by reading the OpenMP directives.

## 2. Model changes decided by the PSOM maintainers

| # | Commit | Change | Effect on results |
|---|--------|--------|-------------------|
| D1 | `528f4b8` | **waves `sponge.f90`: smooth taper.** The taper `cos((j-1)/spl*PI/2)` used integer division, so it was a step: full strength for j=0..15 and zero beyond. `(j-1)/spl` and `(j-(NJ-spl))/spl` now use real division (`dble`), giving the intended cosine/sine taper. | waves results change near the south and north sponges. |
| D2 | `6406683` | **`mldepth` from the namelist.** New optional namelist group `&diag mldepth=… /`, read after `&user` in `model/inc/ini_param.h`, default **10 m**. `mldepth` is the mixed-layer depth used by the N² budget diagnostics (`diag_n2budget`). It was never set for waves, so those diagnostics divided by zero. The waves namelists now set it. Experiments whose `ini_st` builds the initial state from a mixed-layer depth still set it there, which overrides the namelist (wiggle 50 m, NA 250 m, NP_winter 90 m). | Prognostic fields unchanged. waves N² diagnostics are now finite. Namelists without a `&diag` group work as before. |
| D3 | `10f2919` | **Multi-colour Gauss–Seidel / SOR in `linerelax_periodicew.f90`, `hsolve.f90`, `mprove.f90`.** See below. | All runs change at the solver-tolerance level (wiggle KE after 5 steps: 1.7e-9 rel.). **OpenMP runs are now bit-identical to serial runs** for 1, 3 and 4 threads (wiggle 10 steps, waves 3 steps), and repeated runs are identical. |

**D3 in more detail.** These three solvers update the solution in place (Gauss–Seidel / SOR), and the
OpenMP versions split the loops over j across threads. Each update reads the neighbours at i±1 and j±1,
some of which another thread may be updating at the same moment. The result therefore depended on thread
timing, and repeated OpenMP runs gave different answers. In `linerelax` the threads also shared the state
of the random number generator `ran3`, used to choose the line start.

The fix splits the points (columns, in `linerelax`) into colours by the parity of i and j. Points of one
colour don't neighbour each other, so they can be updated in any order, in parallel, with the same
result. The colours are swept one after the other. x is periodic, so when the number of points in x is odd,
the last point gets a third colour. All copies of these routines got the same coloured version (they
differed only in the line start: `NI/2`, random, or 0). This also gave `NP_summer/src/hsolve.f90` the
max-reduction it lacked. The Julia port uses the same ordering, which is also what makes it run in
parallel on a GPU.

## 3. Repository and build

- **`code/model_build/` is no longer tracked** (`256aa69`). `tools/compile.sh` regenerates it on every
  build: it copies `model/src` and `model/inc` into it, copies the experiment's `src` and `inc` on top,
  writes `model_build/optfile`, generates `mkfile/Makefile` with the compiler options from `optfile`, and
  compiles from `model_build`. The tracked copy was left over from an old ifort wiggle build and contained
  older versions of `main.f90`, `ini_uv.f90` and `rpevalgrad_Song.f90` than `model/`, so it could only
  mislead. `compile.sh` now creates `model_build/src` and `inc` if they are missing (fresh clone), and
  `clean.sh` no longer fails when `model_build` doesn't exist. No effect on any build.
- **Build products and editor backups removed** (`b2ae763`): the object files, `.mod` files, generated
  `Makefile` and Linux executable in `mkfile/` (the directory stays, with a `.gitignore`), the compiled
  `tools/preproc.exe` and `tools/preprocmine` (compile `tools/preproc.f90` yourself, as GET_STARTED_3 says),
  `*~` backups, `*.bak`, and stray `.mod` files. `.gitignore` now covers `.DS_Store`, `*.o`, `*.mod`, `*~`,
  `*.bak`. A fresh clone was checked to build as before (`compile.sh` stages `model_build` and writes the
  Makefile).
- Documentation: this file; GET_STARTED_3 CASE 7 (OpenMP and stack size); a pointer in `README`.

## 4. Stack size and OpenMP (documented, not changed)

**What happens.** With OpenMP on (`-qopenmp` for ifort, `-fopenmp` for gfortran), the compiler stores every
local variable of every subroutine on the stack, so that each thread gets its own copy. (In a serial build,
large local arrays are stored statically instead.) Several routines have full-grid local arrays (for example
`g13`/`g23` in `sigma_toplayer` and `staticsigma`, `ch` in `hsolve`). These are called from the main thread,
so they all go on the **main thread's stack**, whose size is limited by `ulimit -s` (often 8 MB by default).
If they don't fit, the run dies with a segmentation fault at the first such call. The worker threads only
hold small per-column arrays, so their stack (`OMP_STACKSIZE` / `KMP_STACKSIZE`) is fine at the default.

**How much is needed.** Measured with gfortran `-fopenmp`, waves on a 16×192×96 grid (295k cells), 4
threads: it crashes with a main-thread stack of 8, 16 and 32 MB and runs with 64 MB. The need grows with the
grid, so the full waves grid (16×768×192) needs about 0.5 GB.

**Why it didn't show up before.** The OpenMP builds of this repository have always been made with ifort
(the gfortran options in `tools/genmakefilel` don't enable OpenMP), and both ifort setups give a big stack:
- macOS (`tools/genmakefilel`, `makefiles/genmakefilel_mac`): the OpenMP flags include
  `-Wl,-stack_size,0x40000000`, which links the executable with a **1 GB** main-thread stack.
- Linux/HPC (`makefiles/genmakefilel_linux`): no stack flag, but clusters usually run jobs with
  `ulimit -s unlimited`.
The crash in our tests came from a gfortran OpenMP build without such a flag. It is unrelated to the
OpenMP race fixed by D3 (which gave wrong answers, not crashes) and to F5 (an uninitialised read that only
showed up once local arrays were on the stack).

**What to do.**
- Linux: run with `ulimit -s unlimited` (in the job script, before starting the model).
- macOS: the stack can't be raised beyond 64 MB with `ulimit`, so link with `-Wl,-stack_size,0x40000000`
  (as the ifort macOS flags already do; gfortran accepts the same option).
- Permanent alternative (not done): make the large local arrays `ALLOCATABLE`, so they live on the heap
  and the stack size stops mattering.

The Julia port is not affected: Julia arrays always live on the heap.

## 5. How the changes were verified

The tests used the tools in `Julia_Port/validation/` (see `Julia_Port/VALIDATION.md` for details). The
Fortran was built with gfortran (`-O2 -fdefault-real-8 -fdefault-double-8`, i.e. the gfortran options
listed in `tools/_values_reference` plus `-O2`), serial unless noted.
- **Fortran vs Julia, whole model.** wiggle (default, Shchepetkin pressure gradient, hydrostatic) and
  waves (on a reduced 16×192×96 grid) were run with both codes from the same namelist, with the full
  model state written every step, and compared field by field. After all fixes they agree to round-off:
  kinetic energy to ≤1e-13 relative and all fields to ≤2e-10 relative over 3–100 steps (the exceptions are
  fields that are themselves pure round-off noise in these flows).
- **Fortran vs Julia, routine by routine.** The Fortran was built as a shared library and loaded into the
  Julia process. Each routine of the time step was run in both codes on identical inputs and every model
  variable compared afterwards (`difftest.jl`).
- **Uninitialised memory.** Builds with `-frecursive -finit-real=snan -ffpe-trap=invalid` (every
  uninitialised real is a signalling NaN that stops the program when used) run clean for wiggle (all three
  configurations) and waves. Builds with `-frecursive` (locals on the stack) are bit-identical to normal
  builds.
- **OpenMP.** After D3, OpenMP builds on 1, 3 and 4 threads are bit-identical to the serial build, and
  repeated 4-thread runs are identical (wiggle 10 steps, waves 3 steps).
- **This tree vs `Julia_Port/code`.** After the bug fixes F1–F7, the Fortran here gave bit-identical
  results to the Fortran copies in `Julia_Port/code` (wiggle, waves); D1–D3 were then applied identically
  to both.

Note: `tools/_values_reference` lists kinetic energies for wiggle from 2013 (V0.5, KE(0) = 2.92e-3). The
current wiggle setup gives KE(0) = 1.02e-3 with or without these changes, so those numbers come from a
different configuration and can't be used as a check.
