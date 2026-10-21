# TMT chemistry with MUSICA MIEM

The MPAS Makefile accepts `TMT=true` as an alternative to `MUSICA=true`.
TMT supplies the MUSICA-compatible interfaces to TurboMICM's Fable solver and
TurboTUV. MUSICA supplies MIEM, including the existing selected-cell, layered
emissions, and grid-metadata interfaces. Namelist settings and host call sites
are unchanged. Selecting TMT does not select a different Fable solver profile;
the existing Rosenbrock selection and configured tolerances still apply.

## Build the dependency interface

Use the same Fortran compiler as MPAS and the pinned MUSICA installation listed
in the main [README](../../README.md). Set `PKG_CONFIG_PATH` so that
`pkg-config --modversion musica-fortran` reports that installation, including its
NetCDF dependencies. TMT also needs a C compiler with OpenMP and a Python
interpreter with PyYAML. Its dependency snapshots come from the sibling
`TurboMICM`, `TurboTUV`, and `MUSICA` Git repositories at the revisions in TMT's
`dependencies.lock`.

From the CheMPAS-A root:

```sh
cmake -S cmake/tmt -B ../chempas-tmt-build \
  -DTMT_SOURCE_DIR="$PWD/../TMT" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/path/to/gcc \
  -DCMAKE_Fortran_COMPILER=/path/to/gfortran \
  -DPython3_EXECUTABLE=/path/to/python-with-pyyaml \
  -DCMAKE_INSTALL_PREFIX="$PWD/../chempas-tmt-install"
cmake --build ../chempas-tmt-build -j8
ctest --test-dir ../chempas-tmt-build --output-on-failure
cmake --install ../chempas-tmt-build

export PKG_CONFIG_PATH="$PWD/../chempas-tmt-install/lib/pkgconfig:$PKG_CONFIG_PATH"
make clean CORE=atmosphere
make -j8 gfortran CORE=atmosphere OPENMP=false PRECISION=double TMT=true
```

The usual MPI, NetCDF, PNetCDF, and PIO environment settings remain necessary.
Use `lib64/pkgconfig` instead if that is your installation's library directory.
The model's OpenMP option is independent of the runtime used by TMT/Fable.
Use the same compilers and OpenMP runtime throughout the dependency build.

The build checks the pinned MUSICA/MIEM metadata for both backends and links a
small probe against the selected chemistry and photolysis modules plus MIEM.
An ordinary `tmt-fortran` installation alone is insufficient: the additional
`chempas-tmt-fortran` package identifies the isolated interface.

## Why the isolated interface is necessary

TMT and MUSICA export identical Fortran module names and C entry points.
Their opaque handles and allocation routines cannot be mixed: for example,
TMT strings use `malloc/free`, while MUSICA strings use `new[]/delete[]`.
Choosing libraries by link order does not provide a safe ownership boundary.

This build adds `libchempas_tmt` and `libchempas_tmt_fortran` alongside the normal
TMT installation. Generated Fortran sources import `tmt_*` modules, and their
C bindings call `ChempasTmt_*` symbols. A generated preprocessor header applies
the same names to TMT's C definitions and internal calls to those entry points.
The dependency checkouts and numerical implementations are untouched. The
additional library uses TMT's existing Fable runtime and prepared-mechanism
format. The coexistence test exercises both libraries' strings/errors, a TMT
chemistry solve, a TurboTUV photolysis calculation, and MUSICA's MIEM lifecycle
in one process.

In CheMPAS-A, `MPAS_USE_TMT` selects conditional `use` statements in
`mpas_musica.F`, `mpas_tuvx.F`, and the framework's version reporting. The
existing `MPAS_USE_MUSICA` feature flag
remains defined for both builds because it enables the common host coupling and
MIEM. `mpas_miem.F` continues to import MUSICA's own utility and emissions types.

## Runtime configuration

TMT can prepare chemistry from the existing configuration path or load a
previously prepared mechanism directory. Source preparation needs the configured
Python interpreter and C compiler at runtime; prepared loading does not. Set
`TMT_CACHE` to a writable mechanism cache, especially for MPI jobs. Consult the
TMT checkout's `docs/BUILD.md` and `docs/COMPATIBILITY.md` for preparation and
supported mechanisms. The interface preserves call signatures; numerical
equivalence and performance across the research cases require separate
validation.

The CheMPAS YAML mechanisms require a TMT configuration reader that accepts
plain phase-species names and uses YAML 1.2 Boolean resolution, so unquoted `NO`
remains a species name. TMT's public photolysis interface must accept solar
zenith angles in radians and convert them to degrees for TurboTUV, matching
MUSICA's public API. Rebuild and reinstall the interface after updating TMT.

See [VALIDATION.md](VALIDATION.md) for the matched full-case numerical results,
compatibility fixes, and measurement status.

## Library timings

Both builds report `[CHEM TIMER]` entries at finalization:

- `MICM solve` measures each `micm%solve` call, including adaptive sub-calls.
  `MICM reference solve` separately measures the optional uncoupled diagnostic.
- `TUV run` measures each daytime `tuvx_solver%run` call. Nighttime columns that
  return zero photolysis without calling the library contribute no calls.

These boundaries are identical for MUSICA and TMT. They include the selected
library's Fortran/C interface overhead, but exclude MPAS state conversion,
profile preparation, output extraction, statistics queries, and logging.
Initialization and mechanism preparation are outside these solve timers.

The report gives the minimum, mean, and maximum accumulated wall seconds per
MPI rank, plus the total call count across ranks. Use the maximum when comparing
with elapsed model time; the mean helps identify workload imbalance. The
existing `chem MICM solve` timer encloses a broader host operation and should
not be confused with the new library timer. The existing `chemistry` timer
does not include photolysis.

The counters use 64-bit `system_clock` ticks locally and perform MPI reductions
only at finalization. This allows ranks to execute different numbers of
photolysis calls without introducing synchronization inside the column loop.
Run each timed model alone with the same MPI layout, threads, compiler options,
inputs, and mechanism-cache state. Concurrent validation runtimes are not
controlled performance measurements.
