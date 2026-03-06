# Claude Code Project Configuration

## Project Overview

CheMPAS (Chemistry for MPAS) is a standalone project derived from [NCAR/MPAS-Model-ACOM-dev](https://github.com/NCAR/MPAS-Model-ACOM-dev), itself a fork of [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model). It extends the Model for Prediction Across Scales (MPAS) with integrated MUSICA/MICM atmospheric chemistry, enabling coupled atmospheric-chemistry modeling on MPAS's unstructured mesh. CheMPAS is an independent project and is not intended to sync with the upstream repositories.

## Key Documentation

- `ARCHITECTURE.md` - High-level system architecture
- `BUILD.md` - Build system documentation (includes LLVM/macOS instructions)
- `RUN.md` - Test case execution instructions
- `VISUALIZE.md` - Chemistry visualization tools and workflows
- `COMPONENTS.md` - Detailed component descriptions
- `MUSICA_INTEGRATION.md` - MUSICA/MICM coupling details
- `MUSICA_API.md` - MUSICA Fortran API reference
- `AGENTS.md` - Agent roles, workflow, manifesto, and operational details
- `BENCHMARKS.md` - Agent model benchmark comparison
- `PURPOSE.md` - Project motivation and goals
- `TEST_RUNS.md` - Test run documentation
- `TODO.md` - Development task list

## Sister Project Reference

- `~/EarthSystem/DAVINCI-MPAS/` - Sister project with prior lightning-NOx/O3
  implementation work.
- Start with:
  - `PLAN.md` (Phase 6 LNOx-O3 details)
  - `SCIENCE.md` (physics assumptions and validation framing)
  - `DC3.md` (observational context and targets)

## Build Configuration (macOS with LLVM)

This project is configured to build with LLVM compilers (flang/clang) via Homebrew:

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"

# Build atmosphere core
make -j8 llvm \
  CORE=atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double

# Build init_atmosphere core
make -j8 llvm \
  CORE=init_atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double
```

**Important:** Always use `-j8` for parallel compilation.

## Key Source Locations

| Component | Path |
|-----------|------|
| Framework | `src/framework/` |
| Operators | `src/operators/` |
| Atmosphere Core | `src/core_atmosphere/` |
| Init Atmosphere | `src/core_init_atmosphere/` |
| Dynamics | `src/core_atmosphere/dynamics/` |
| Physics | `src/core_atmosphere/physics/` |
| Chemistry | `src/core_atmosphere/chemistry/` |
| MUSICA Integration | `src/core_atmosphere/chemistry/musica/` |
| Registry (metadata) | `src/core_atmosphere/Registry.xml` |

## MUSICA Integration Files

| File | Purpose |
|------|---------|
| `mpas_atm_chemistry.F` | Generic chemistry interface |
| `mpas_musica.F` | MUSICA/MICM coupler implementation |

## Build Artifacts

| File | Description |
|------|-------------|
| `atmosphere_model` | Main atmosphere executable |
| `init_atmosphere_model` | Initialization/preprocessing tool |
| `build_tables` | Physics lookup table generator |

## Available Skills

Skills are defined in `.claude/commands/` and can be invoked with `/skillname`:

| Skill | Command | Description |
|-------|---------|-------------|
| Build MPAS | `/build-mpas` | Build atmosphere or init_atmosphere with LLVM |
| Test MPAS | `/test-mpas` | Run a test case and verify output |
| Plot Chemistry | `/plot-chemistry` | Generate and open chemistry visualization plots |
| MUSICA Dev | `/musica-dev` | Work on MUSICA/MICM chemistry integration |

## Development Manifesto

CheMPAS is an agent-driven development project. See `AGENTS.md` for the full manifesto, agent roles, workflow, and human review gates.

## Multi-Agent Workflow

CheMPAS uses three agents from three vendors (see `AGENTS.md`). Codex 5.3 reviews Claude's work and may push doc edits or review findings directly to the repo.

**Convention for cross-agent communication:**
- Codex writes review findings to `CODEX_REVIEW.md` and pushes to `develop`
- When the user says "check the review", pull and read `CODEX_REVIEW.md`
- Implement fixes, then remove `CODEX_REVIEW.md` (it's a transient artifact, not a permanent doc)
- Codex may also edit `.md` documentation files directly — always `git pull` before committing if Codex has been active
- Before merging to main, ensure no concurrent Codex edits are in flight

## Notes for AI Assistants

1. **Fortran Standards**: This is a Fortran 2008 codebase using MPI for parallelism
2. **Build System**: Uses legacy Makefile (not CMake) with `make llvm` target for LLVM compilers
3. **Module Files**: LLVM flang `.mod` files are incompatible with gfortran - don't mix compilers
4. **MPI**: Uses `include 'mpif.h'` via `NOMPIMOD` flag (not `use mpi` module)
5. **Registry System**: Most variable definitions are in `Registry.xml` (build-time), but MUSICA chemistry tracers are injected at runtime from the MICM config
6. **Physics**: Many physics schemes from WRF, NoahMP, and other sources
7. **MUSICA**: MUSICA-Fortran must be built with flang; `musica-fortran.pc` may need yaml-cpp path fix. See `BUILD.md` for details.
8. **Testing**: Always run with 8 MPI ranks (`mpiexec -n 8`). A mismatched rank count with no matching partition file causes segfaults in the dynamics solver. See `RUN.md` for details.

## Common Tasks

### Full Rebuild
```bash
make clean CORE=atmosphere
find . -name "*.mod" -delete
find . -name "*.o" -delete
make -j8 llvm CORE=atmosphere PIO=$HOME/software NETCDF=/opt/homebrew PNETCDF=$HOME/software PRECISION=double
```

### Build with MUSICA
```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"

make -j8 llvm \
  CORE=atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double \
  MUSICA=true
```

**Note:** The `PKG_CONFIG_PATH` export is required for MUSICA builds.

### Check Build Status
```bash
ls -la atmosphere_model init_atmosphere_model build_tables
```

## Python Environment and Scripts

### Conda Environment
A conda environment `mpas` is available for visualization and analysis:
```bash
~/miniconda3/envs/mpas/bin/python  # Direct path
# Or activate: conda activate mpas
```

Packages: numpy, matplotlib, netcdf4

### Visualization Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `plot_chemistry.py` | Visualize chemistry tracer output (qA, qB, qAB) |
| `init_tracer_sine.py` | Initialize tracers with sine wave patterns |

**plot_chemistry.py usage:**
```bash
cd ~/Data/MPAS/supercell
~/miniconda3/envs/mpas/bin/python plot_chemistry.py -o output.png
~/miniconda3/envs/mpas/bin/python plot_chemistry.py --time-series  # Spatial evolution
~/miniconda3/envs/mpas/bin/python plot_chemistry.py --diff          # Difference from t=0
```

**init_tracer_sine.py usage:**
```bash
# Apply horizontal sine wave to qAB (for advection studies)
python init_tracer_sine.py -i supercell_init.nc -t qAB --waves-x 2 --amplitude 0.4 --offset 0.6
```

## Test Run Directory

The supercell test case is located at `~/Data/MPAS/supercell/` (shared with MPAS-Model-ACOM-dev). See `RUN.md` for execution instructions.

**Important:** The `streams.atmosphere` file must use `io_type="netcdf"` to avoid PnetCDF compatibility issues on macOS/LLVM builds.
