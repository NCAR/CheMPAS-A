CheMPAS-A
=========

CheMPAS-A (Chemistry for MPAS - Atmosphere) is an ACOM integration pilot that couples
MUSICA/MICM atmospheric chemistry to MPAS-Atmosphere on its native
unstructured Voronoi mesh. The repository serves as a rapid-prototyping
ground for the chemistry coupling: runtime tracer allocation, MUSICA/MICM
state transfer, TUV-x photolysis, and idealized chemistry test cases are
developed here ahead of any upstream integration.

CheMPAS-A is decoupled from both
[MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model) and
[NCAR/MPAS-Model-ACOM-dev](https://github.com/NCAR/MPAS-Model-ACOM-dev) —
there is no sync or fork-tracking mechanism. Mature pieces are contributed
back to `MPAS-Model-ACOM-dev` as deliberate, hand-crafted pull requests,
which is the staging ground for eventual upstream integration into
`MPAS-Dev/MPAS-Model`.

## Upstream Integration Path

CheMPAS-A does not push its history to `MPAS-Model-ACOM-dev`. The two repos
share ancestry but have no merge or rebase relationship. Integration
contributions take the form of focused pull requests that reimplement a
mature CheMPAS-A feature against the current `MPAS-Model-ACOM-dev` tree.

**Baseline:** CheMPAS-A derives from MPAS-Model v8.3.1 (commit `b9090a143`).
At the file level, every v8.3.1 source file is present at the same path in
CheMPAS-A (1197 / 1197), 97.6% of them byte-identical; the divergence is
entirely additive (chemistry coupling, MUSICA/MICM/TUV-x integration, build
fixes for macOS LLVM/flang) plus per-file modifications. See
[docs/chempas/upstream/2026-04-19-vs-mpas-v8.3.1.md](docs/chempas/upstream/2026-04-19-vs-mpas-v8.3.1.md)
for the full systematic comparison.

In practice this means each upstream PR is:

- **Scoped to one capability** — e.g., runtime tracer allocation, the
  MUSICA/MICM coupler, or the TUV-x cloud radiator — rather than a bulk
  port of CheMPAS-A state.
- **Re-derived against the current upstream tree**, so the diff is clean
  against `MPAS-Model-ACOM-dev`'s `Registry.xml`, build system, and module
  layout at PR time.
- **Backed by the prototype evidence in this repo** — test cases,
  validation runs, and design notes — but not a literal port of every
  CheMPAS-A file.

This separation keeps CheMPAS-A free to iterate quickly while keeping
upstream PRs reviewable as standalone changes.

## What's Working

The pieces below are demonstrated end-to-end in this repo and are the
primary candidates for upstream integration PRs. Each entry points to the
code and the relevant deeper documentation.

- **Runtime chemistry tracer system** — chemistry tracers are removed from
  `Registry.xml` and discovered at startup from the active MICM YAML
  configuration. The `atm_extend_scalars_for_chemistry` hook, called from
  `atm_setup_block` before field-array allocation, queries MICM for the
  species list and extends the `scalars` pool's metadata in place
  (`num_scalars` dimension, `constituentNames`, `attLists`, index dims) for
  both time levels and the tend pool. The framework's normal allocation
  path then sizes field arrays from the updated `num_scalars`. Two known
  constraints: incompatible with `config_apply_lbcs=true` (because
  `lbc_scalars` is statically sized from Registry metadata), and currently
  hardcoded to the chemistry use case rather than a generic runtime
  pool-extension mechanism. **Stopgap intent:** this implementation is
  meant to be superseded by the generic runtime data-pool-variables
  infrastructure planned for a future MPAS-Dev release; CheMPAS-A will
  migrate when that lands upstream.
  Code: `src/core_atmosphere/mpas_atm_core_interface.F:657`
  (`atm_extend_scalars_for_chemistry`),
  `src/core_atmosphere/chemistry/mpas_atm_chemistry.F:999`
  (`chemistry_query_species`),
  `src/core_atmosphere/chemistry/musica/mpas_musica.F:232`
  (`musica_query_species`),
  `src/framework/mpas_block_creator.F` (allocates arrays from
  `num_scalars` after extension).

- **MUSICA/MICM coupler** — state transfer between MPAS and MICM in
  mol/m³, unit conversion, and external rate-parameter wiring.
  Code: `src/core_atmosphere/chemistry/musica/mpas_musica.F`,
  `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`.
  See [docs/chempas/musica/MUSICA_INTEGRATION.md](docs/chempas/musica/MUSICA_INTEGRATION.md).

- **TUV-x clear-sky photolysis** — j_NO2 computed via the delta-Eddington
  solver with a from-host wavelength grid (CAM 102-bin grid).
  See [docs/chempas/guides/TUVX_INTEGRATION.md](docs/chempas/guides/TUVX_INTEGRATION.md).

- **TUV-x cloud opacity from host** — cloud radiator built from MPAS LWC
  (`tau = 3·LWC·dz / (2·r_eff·ρ_water)`, SSA = 0.999999, g = 0.85),
  attached to the TUV-x solver before construction.
  See [docs/chempas/guides/TUVX_INTEGRATION.md](docs/chempas/guides/TUVX_INTEGRATION.md).

- **Idealized chemistry test cases** — supercell (with LNOx source),
  mountain wave, and Jablonowski–Williamson baroclinic wave configurations,
  all wired for MUSICA tracer transport.
  See [test_cases/](test_cases/) and
  [docs/chempas/results/TEST_RUNS.md](docs/chempas/results/TEST_RUNS.md).

- **Chemistry visualization tooling** — Python scripts for chemistry
  tracers, LNOx/O₃ evolution, and chemistry profile plotting (frame
  selection and time-series modes).
  See [scripts/](scripts/) and
  [docs/chempas/guides/VISUALIZE.md](docs/chempas/guides/VISUALIZE.md).

## Building

CheMPAS-A builds on macOS (LLVM/flang) and Ubuntu (GCC/gfortran via conda).
Both paths use the same preflight script to detect the toolchain and
export the required environment.

External dependencies: MPI, NetCDF, PnetCDF, PIO, and (for chemistry)
MUSICA-Fortran built with the matching Fortran compiler.

```bash
scripts/check_build_env.sh                       # report mode
eval "$(scripts/check_build_env.sh --export)"    # export mode
```

macOS (LLVM/flang):

```bash
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
  CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true
```

Ubuntu (GCC/gfortran, requires `conda activate mpas`):

```bash
eval "$(scripts/check_build_env.sh --export)" && make -j8 gfortran \
  CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true
```

Notes:
- `PKG_CONFIG_PATH` must be present in the same shell invocation as
  `make` — the Makefile invokes `pkg-config` at parse time.
- Do not mix flang and gfortran `.mod` files; rebuild MUSICA-Fortran with
  the same compiler as CheMPAS-A.

See [BUILD.md](BUILD.md) for the full preflight, troubleshooting, and
dependency-build notes, and [RUN.md](RUN.md) for executing test cases.

## Code Layout

```
CheMPAS-A/
├── src/
│   ├── driver                       -- Standalone driver
│   ├── framework                    -- MPAS framework (pools, fields, I/O, MPI)
│   ├── operators                    -- Mesh operators
│   ├── external                     -- Vendored external dependencies
│   ├── tools/
│   │   ├── registry                 -- Registry.xml parser / code gen
│   │   └── input_gen                -- Stream and namelist generators
│   ├── core_atmosphere/
│   │   ├── dynamics                 -- Dynamical core
│   │   ├── physics                  -- Physics parameterizations
│   │   ├── diagnostics              -- Diagnostic outputs
│   │   └── chemistry/
│   │       ├── mpas_atm_chemistry.F -- Generic chemistry interface
│   │       └── musica/              -- MUSICA/MICM coupler (mpas_musica.F)
│   ├── core_init_atmosphere         -- Initialization / preprocessing core
│   └── core_{ocean,seaice,landice,sw,test}  -- Inherited from upstream MPAS;
│                                              not actively maintained here
├── micm_configs/                    -- MICM YAML mechanism configs
├── test_cases/                      -- Idealized test case configurations
├── default_inputs/                  -- Default streams and namelists
├── scripts/                         -- Visualization and analysis tools
└── docs/                            -- Documentation tree
```

## Documentation

The most relevant docs for evaluating or extracting pieces of this pilot:

**Architecture & integration**

| Document | Description |
|----------|-------------|
| [docs/chempas/architecture/ARCHITECTURE.md](docs/chempas/architecture/ARCHITECTURE.md) | System architecture overview |
| [docs/chempas/architecture/COMPONENTS.md](docs/chempas/architecture/COMPONENTS.md) | Component-level details |
| [docs/chempas/musica/MUSICA_INTEGRATION.md](docs/chempas/musica/MUSICA_INTEGRATION.md) | MUSICA/MICM coupling design and implementation |
| [docs/chempas/musica/MUSICA_API.md](docs/chempas/musica/MUSICA_API.md) | MUSICA-Fortran API reference |
| [docs/chempas/guides/TUVX_INTEGRATION.md](docs/chempas/guides/TUVX_INTEGRATION.md) | TUV-x photolysis integration and test case |

**Build, run, and validation**

| Document | Description |
|----------|-------------|
| [BUILD.md](BUILD.md) | Build instructions and preflight notes |
| [RUN.md](RUN.md) | Running test cases |
| [docs/chempas/results/TEST_RUNS.md](docs/chempas/results/TEST_RUNS.md) | Recorded run outcomes and validation evidence |
| [docs/chempas/guides/VISUALIZE.md](docs/chempas/guides/VISUALIZE.md) | Chemistry visualization tools |

The full topic-organized index lives at
[docs/chempas/README.md](docs/chempas/README.md), and the MPAS user guide
chapters are imported under [docs/users-guide/](docs/users-guide/).

## License

See [LICENSE](LICENSE).
