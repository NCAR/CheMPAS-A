# CheMPAS-A Ubuntu/GCC Build Plan

**Status: COMPLETE (2026-04-05)**

Ported CheMPAS-A from macOS/LLVM (flang/clang) to Ubuntu/GCC (gfortran/gcc).

## Current State

**Available:**
- gcc/gfortran 15.2.0 (conda-forge)
- MUSICA source at `~/EarthSystem/MUSICA/`

**Missing:**
- MPI (no mpifort/mpicc)
- CMake
- NetCDF C and Fortran libraries
- PnetCDF
- PIO (Parallel I/O)

**Codebase issues:**
- Makefile MUSICA integration hardcodes LLVM paths (`MUSICA-LLVM/flang-deps/...`)
- `scripts/check_build_env.sh` hardcodes `flang`/`clang` compilers
- `-DNOMPIMOD` flag was an LLVM/Homebrew MPI workaround (likely unnecessary with gfortran)

---

## Phase 1: Create Conda Environment and Install Dependencies

Create a dedicated `mpas` conda environment with the full GCC/MPI/NetCDF toolchain
from conda-forge. This avoids the gcc 15 (base conda) vs g++ 13 (system) mismatch
and keeps the build stack isolated.

```bash
conda create -n mpas -c conda-forge \
  gcc gxx gfortran cmake openmpi netcdf-c netcdf-fortran parallel-netcdf \
  pkg-config make
conda activate mpas
```

Verify after install:
- `gfortran --version` — consistent version
- `mpifort --version` shows gfortran
- `nc-config --all` shows paths
- `cmake --version` >= 3.21
- `pkg-config --libs netcdf-fortran` resolves

**Note:** All subsequent phases assume the `mpas` environment is active.

---

## Phase 2: Build and Install MUSICA

Build `~/EarthSystem/MUSICA` with gfortran and install to a local prefix.

```bash
cd ~/EarthSystem/MUSICA
mkdir -p build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/software \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx
make -j8
make install
```

Verify:
- `pkg-config --cflags musica-fortran` returns include paths
- `pkg-config --libs musica-fortran` returns link flags
- `.mod` files are gfortran-compatible

---

## Phase 3: Adapt CheMPAS-A Build System

### 3a. Makefile (Lines 884-902)

Remove the hardcoded LLVM netcdf-fortran path:

```makefile
# Current (line 896):
MUSICA_LIBS += $(HOME)/EarthSystem/MUSICA-LLVM/flang-deps/netcdf-fortran-install/lib/libnetcdff.a

# Replace with: use system/conda libnetcdff (gfortran-built, no mangling conflict)
# Delete the line or make it conditional on LLVM builds
```

### 3b. scripts/check_build_env.sh

Update to support GCC toolchain:
- Replace hardcoded `flang` checks with generic Fortran compiler detection
- Replace `OMPI_FC=flang` / `OMPI_CC=clang` / `OMPI_CXX=clang++` with
  `OMPI_FC=gfortran` / `OMPI_CC=gcc` / `OMPI_CXX=g++` (or auto-detect)
- Update MUSICA path search from `MUSICA-LLVM/build` to `MUSICA/build`
- Update suggested make target from `llvm` to `gfortran`

### 3c. Make Target

The Makefile already has a `gfortran` target (~line 375). Verify it has:
- `FC_SERIAL = gfortran`
- `CC_SERIAL = gcc`
- `CXX_SERIAL = g++`
- Compatible flags (`-fdefault-real-8`, `-fdefault-double-8`, `-fconvert=big-endian`)

Check whether `-DNOMPIMOD` is needed. With gfortran-built OpenMPI, the `mpi.mod`
should be compatible, so `NOMPIMOD` may not be required.

---

## Phase 4: Build CheMPAS-A

### 4a. Build without MUSICA first

```bash
export PKG_CONFIG_PATH=$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH
make clean CORE=atmosphere
make -j8 gfortran \
  CORE=atmosphere \
  PIO="$PIO" \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double
```

### 4b. Build with MUSICA

```bash
make clean CORE=atmosphere
make -j8 gfortran \
  CORE=atmosphere \
  PIO="$PIO" \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double \
  MUSICA=true
```

---

## Phase 5: Update Documentation

- `BUILD.md` — Add Ubuntu/GCC build section alongside existing LLVM/macOS section
- `CLAUDE.md` — Update build configuration for Ubuntu paths and GCC toolchain
- `scripts/check_build_env.sh` — Document new compiler-suite support
- `README.md` — Note multi-platform support

---

## Key Risks and Notes

1. **Compiler consistency**: MUSICA, NetCDF-Fortran, and CheMPAS-A must all be built
   with the same Fortran compiler. Mixing gfortran and flang `.mod` files causes
   link failures.

2. **PIO**: May not be available via conda. May need to build from source
   ([E3SM PIO](https://github.com/NCAR/ParallelIO)).

3. **NOMPIMOD**: Test without `-DNOMPIMOD` first. If gfortran's `mpi.mod` is
   compatible with the installed OpenMPI, it should work without this flag.

4. **Physics MMM submodule**: The build may try to `git fetch` MMM-physics.
   Ensure network access or pre-populate the submodule.

5. **Conda environment isolation**: All builds should happen within the same
   conda environment to ensure consistent compiler/library paths.

---

## Completion Notes (2026-04-05)

All five phases completed. Additional issues encountered and resolved:

1. **Conda package names differ from plan**: `netcdf-c` → `libnetcdf`,
   `parallel-netcdf` → `libpnetcdf`
2. **MUSICA needed BLAS**: Added `libblas`/`liblapack` to conda install
3. **`-lstdc++` required**: MUSICA's C++ code needs the C++ standard library
   at link time; added to `musica-fortran.pc`
4. **`checkout_data_files.sh` missing**: Created script to download WRF physics
   lookup tables from MPAS-Data GitHub repo
5. **NoahMP source files missing**: Copied from upstream MPAS-Dev/MPAS-Model
   (`physics_noahmp/{utility,src,drivers/mpas}/*.F90`)
6. **GSL orography files missing**: Copied `mpas_gsl_oro_data_{sm,lg}_scale.F`
   from upstream for `init_atmosphere` build
7. **NOMPIMOD not needed**: Confirmed — gfortran + conda OpenMPI uses `mpi_f08`
   natively
8. **Makefile LLVM path removed**: Deleted hardcoded
   `MUSICA-LLVM/flang-deps/netcdf-fortran-install/lib/libnetcdff.a` line
9. **`check_build_env.sh` updated**: Auto-detects flang vs gfortran, searches
   conda paths, supports `.so` (Linux) in addition to `.dylib` (macOS)
