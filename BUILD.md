# MPAS Build System

This document describes the build systems available for MPAS and their configuration options.

## Build Systems

MPAS supports two build systems:

1. **CMake** (Primary, recommended)
2. **Legacy Makefile** (Platform-specific targets)

## CMake Build

### Basic Usage

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `DO_PHYSICS` | ON | Enable built-in physics schemes |
| `MPAS_DOUBLE_PRECISION` | ON | Use 64-bit floating point |
| `MPAS_PROFILE` | OFF | Enable GPTL profiling |
| `MPAS_OPENMP` | OFF | Enable OpenMP parallelization |
| `BUILD_SHARED_LIBS` | ON | Build shared libraries |
| `MPAS_USE_PIO` | OFF | Use Parallel I/O library |
| `MPAS_USE_MUSICA` | OFF | Enable MUSICA/MICM chemistry |

### Selecting Cores

By default, both `atmosphere` and `init_atmosphere` cores are built. To build specific cores:

```bash
cmake -DCORES="atmosphere" ..
```

### Example Configurations

**Basic atmosphere build:**
```bash
cmake -DCORES="atmosphere;init_atmosphere" \
      -DMPAS_DOUBLE_PRECISION=ON \
      -DDO_PHYSICS=ON \
      ..
```

**With MUSICA chemistry:**
```bash
cmake -DMPAS_USE_MUSICA=ON \
      -DMUSICA_ROOT=/path/to/musica \
      ..
```

**Debug build:**
```bash
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DMPAS_PROFILE=ON \
      ..
```

## Legacy Makefile Build

The `Makefile` supports various platform-specific targets:

### Available Targets

| Target | Compiler Suite |
|--------|---------------|
| `gnu` | GNU Fortran/C/C++ |
| `llvm` | LLVM flang/clang/clang++ |
| `xlf` | IBM XL compilers |
| `xlf-summit-omp-offload` | IBM XL with OpenMP offload (Summit) |
| `ftn` | Cray compilers |
| `intel` | Intel compilers |
| `pgi` | PGI/NVIDIA compilers |
| `nag` | NAG Fortran |

### Usage

```bash
make -j8 gnu CORE=atmosphere
make -j8 gnu CORE=init_atmosphere
```

**Note:** Always use `-j8` (or appropriate number for your system) for parallel compilation.

### LLVM Build on macOS (Homebrew)

Building MPAS with LLVM compilers (flang/clang) on macOS requires special configuration due to compiler and library compatibility issues.

#### Prerequisites

Install LLVM compilers and dependencies via Homebrew:
```bash
brew install llvm flang open-mpi netcdf netcdf-fortran
```

**Important:** Homebrew's Open MPI is compiled with gfortran, so its Fortran module files (`.mod`) are incompatible with flang. The build system handles this via the `NOMPIMOD` flag.

#### PnetCDF Requirement

Homebrew's PnetCDF is built with gfortran, making its Fortran interface incompatible with flang. Build PnetCDF from source (C-only):

```bash
# Download pnetcdf-1.14.1 from https://parallel-netcdf.github.io/
cd /tmp/pnetcdf-1.14.1
export OMPI_FC=flang
export OMPI_CC=clang
export OMPI_CXX=clang++

./configure --prefix=$HOME/software \
  --disable-cxx \
  --disable-shared \
  --disable-fortran \
  MPICC=mpicc \
  CFLAGS="-O3"

make -j8 && make install
```

#### PIO Library

Build PIO with LLVM compilers (see PIO documentation). Install to `$HOME/software`.

#### Building MPAS

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"

make -j8 llvm \
  CORE=atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double
```

For `init_atmosphere`:
```bash
make -j8 llvm \
  CORE=init_atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double
```

#### LLVM Build Technical Details

The `llvm` target configures:
- **MPI Wrappers:** Sets `OMPI_FC=flang`, `OMPI_CC=clang`, `OMPI_CXX=clang++`
- **Fortran Flags:** Uses gfortran-compatible flags (`-fconvert=big-endian`, `-ffree-form`, `-fdefault-real-8`)
- **MPI Module:** Uses `-DNOMPIMOD` to use `include 'mpif.h'` instead of gfortran's binary `mpi.mod`

| NVIDIA/PGI Flag | LLVM flang Equivalent |
|-----------------|----------------------|
| `-Mbyteswapio` | `-fconvert=big-endian` |
| `-Mfreeform` | `-ffree-form` |
| `-r8` | `-fdefault-real-8` |
| `-Mstandard` | (not needed) |

#### Building with MUSICA Chemistry

To enable MUSICA/MICM atmospheric chemistry, MUSICA-Fortran must be built with flang (not gfortran) since `.mod` files are compiler-specific.

**Both `libmusica.a` and `libnetcdff.a` must be flang-built.** With TUV-x
integration (`-DMUSICA_ENABLE_TUVX=ON`), `libmusica.a` contains
`netcdf.F90.o` objects that reference flang-mangled `__QMnetcdfP*` symbols.
Homebrew's `netcdf-fortran` is built with gfortran (`___netcdf_MOD_*`
mangling) and cannot resolve those references, producing link errors like:

```
Undefined symbols for architecture arm64:
  "__QMnetcdfPnf90_put_var_eightbytereal", referenced from:
      __QMmusica_io_netcdfPappend_2d_double in libmusica.a[...]
```

`MUSICA-LLVM` bundles a flang-built `netcdf-fortran` at
`flang-deps/netcdf-fortran-install/`. Point CheMPAS-A at it with `NETCDFF=`.
The Makefile already knows how to split `NETCDF` (C library / headers) from
`NETCDFF` (Fortran library / `.mod`). The preflight script detects this
tree automatically and exports `NETCDFF=`.

Before running `make ... MUSICA=true`, verify these three things:

```bash
scripts/check_build_env.sh
eval "$(scripts/check_build_env.sh --export)"

pkg-config --modversion musica-fortran
pkg-config --libs musica-fortran
test -n "$PNETCDF"
```

If any of those fail, the top-level Makefile aborts before it reaches the atmosphere sources.

The script auto-detects the local `NETCDF`, `PNETCDF`, `PIO`, and sibling `../MUSICA-LLVM/build` layout used in this repo. If it finds a local MUSICA build tree whose `musica-fortran.pc` points at the wrong prefix, it writes a temporary override in `/tmp/chemps-musica-fortran.pc` and exports a corrected `PKG_CONFIG_PATH`.

`PKG_CONFIG_PATH` must be present in the environment that launches `make`. The Makefile uses `$(shell pkg-config ...)` during parse, so exporting it in a separate shell and then invoking `make` later is not sufficient.

**`musica-fortran.pc` must describe the real build/install tree.**

This is the recurring failure mode: `pkg-config` finds a `musica-fortran.pc`, but that file points at the wrong prefix (often `/usr/local`) instead of the actual MUSICA build tree. In that case, `make llvm ... MUSICA=true` fails in `musica_fortran_test` before compiling MPAS.

If you are using an uninstalled local MUSICA build tree, create a local override:

```bash
cat > /tmp/musica-fortran.pc <<'EOF'
prefix=/path/to/MUSICA-LLVM/build
exec_prefix=${prefix}
libdir=${prefix}/lib
includedir=${prefix}/mod_fortran

yaml_lib=yaml-cpp

Name: musica-fortran
Description: Fortran wrapper for the MUSICA library for modeling atmospheric chemistry
Version: 0.14.5
Cflags: -I${includedir}
Libs: -L${libdir} -lmusica-fortran -lmusica -lmechanism_configuration -l${yaml_lib}
EOF

export PKG_CONFIG_PATH="/tmp:$PKG_CONFIG_PATH"
pkg-config --cflags --libs musica-fortran
```

If yaml-cpp is installed in a non-default location, add its library directory to the `Libs:` line in the override (for example `-L/opt/homebrew/opt/yaml-cpp/lib`).

**Build command with MUSICA (recommended — uses preflight):**

```bash
eval "$(scripts/check_build_env.sh --export)" && \
  make -j8 llvm \
    CORE=atmosphere \
    PIO="$PIO" \
    NETCDF="$NETCDF" \
    NETCDFF="$NETCDFF" \
    PNETCDF="$PNETCDF" \
    PRECISION=double \
    MUSICA=true
```

The preflight exports `NETCDFF=$HOME/EarthSystem/MUSICA-LLVM/flang-deps/netcdf-fortran-install`
on flang hosts. The `eval ... && make` must run in a single shell
invocation (the Makefile's parse-time `$(shell pkg-config ...)` depends
on `PKG_CONFIG_PATH` being inherited directly).

The build will report `MPAS was linked with the MUSICA-Fortran library version X.Y.Z` on success.

**Current versions:** MUSICA-Fortran 0.13.0, MICM 3.10.0

#### MUSICA Build Pitfalls

1. **`PKG_CONFIG_PATH` must resolve a working `musica-fortran.pc`.** If `pkg-config` cannot find the package, or if the `.pc` file points at the wrong prefix, the build aborts in `musica_fortran_test` before MPAS source compilation starts.

2. **The MUSICA modules must be built with flang for the `llvm` target.** If `mpif90` resolves to gfortran while MUSICA was built with flang, or vice versa, module loading fails with errors like `... is not a GNU Fortran module file`.

3. **`netcdf-fortran` must also be flang-built for MUSICA+TUV-x.** Homebrew's `libnetcdff` is gfortran-mangled and cannot satisfy flang-mangled `__QMnetcdfP*` references from `libmusica.a` (pulled in by TUV-x's `netcdf.F90` object). Set `NETCDFF=` to the flang-built tree at `$HOME/EarthSystem/MUSICA-LLVM/flang-deps/netcdf-fortran-install`. The preflight detects and exports this automatically.

4. **`PNETCDF` is mandatory for the normal top-level build.** If `PNETCDF` is unset, the build aborts in `pnetcdf_test` before the atmosphere core is compiled.

5. **`-fdefault-real-8` affects `kind=` specifiers.** Under the LLVM flags, `real(kind=4)` becomes 8-byte. Use `RKIND` (from `mpas_kind_types`) for MPAS log calls (`realArgs`) and any real literals that interface with the framework. Do not use `kind=4`.

6. **Chemistry tracers are runtime-injected (not registry-defined).** The atmosphere registry no longer hardcodes chemistry tracers (`qAB/qA/qB`). With MUSICA enabled, `atm_extend_scalars_for_chemistry()` discovers MICM species at startup and extends `scalars`/`scalars_tend` metadata dynamically.

7. **MICM API returns scalars.** `micm%get_species_property_double(name, property, error)` returns a scalar `real(real64)`, not an array. Assigning to an allocatable array causes undefined behavior.

8. **Species ordering index vs iteration order.** `state%species_ordering%name(i)` gives the i-th species name in iteration order, but `state%species_ordering%index(name)` gives the stride-based index for concentration array access. These are not the same — always use `%index(name)` for array indexing.

#### Known Warnings

These warnings are expected and harmless:
- `-Wl,-flat_namespace: 'linker' input unused` - MPI wrapper passes unused linker flag
- `NUMERIC_STORAGE_SIZE from ISO_FORTRAN_ENV is not well-defined` - Due to `-fdefault-real-8`
- `Unrecognized compiler directive was ignored [-Wignored-directive]` - Intel `!DIR$` directives
- `Reference to procedure has an implicit interface` - Due to `include 'mpif.h'` usage

### GCC Build on Ubuntu (conda-forge)

Building MPAS with GCC/gfortran on Ubuntu using a conda-forge toolchain.

#### Prerequisites

Create a conda environment with the full toolchain:
```bash
conda create -n mpas -c conda-forge \
  gcc gxx gfortran cmake openmpi libnetcdf netcdf-fortran libpnetcdf \
  libblas liblapack pkg-config make
conda activate mpas
```

Verify:
- `gfortran --version` — GCC 15.x (conda-forge)
- `mpifort --version` — wraps gfortran
- `nc-config --prefix` — points to conda env
- `pkg-config --libs netcdf-fortran` — resolves

#### PIO Library

PIO is not available via conda and must be built from source:

```bash
git clone https://github.com/NCAR/ParallelIO.git
cd ParallelIO && mkdir build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/software \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DPIO_ENABLE_TIMING=OFF
make -j8 && make install
```

#### MUSICA Library

Build MUSICA from source with gfortran:

```bash
cd ~/EarthSystem/MUSICA
mkdir -p build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/software \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DMUSICA_BUILD_FORTRAN_INTERFACE=ON \
  -DMUSICA_ENABLE_MPI=ON \
  -DMUSICA_ENABLE_TESTS=OFF
make -j8 && make install
```

**Note:** The installed `musica-fortran.pc` must include `-lstdc++` in its
`Libs:` line for Fortran-driven linking to resolve C++ standard library
symbols.

#### Building MPAS

The preflight script auto-detects the GCC toolchain:

```bash
eval "$(scripts/check_build_env.sh --export)"

make -j8 gfortran \
  CORE=atmosphere \
  PIO=$HOME/software \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double \
  MUSICA=true
```

For `init_atmosphere`:
```bash
make -j8 gfortran \
  CORE=init_atmosphere \
  PIO=$HOME/software \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double
```

#### GCC Build Technical Details

The `gfortran` target configures:
- **MPI Wrappers:** Uses `mpif90`/`mpicc`/`mpicxx` (gfortran is the default wrapper compiler)
- **Fortran Flags:** `-fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-form -ffree-line-length-none`
- **MPI Module:** Uses `mpi_f08` module natively (no `NOMPIMOD` needed)

#### GCC-Specific Notes

1. **No `NOMPIMOD` needed.** With gfortran-built OpenMPI, the `mpi_f08` module
   is compatible. The `gfortran` make target does not pass `-DNOMPIMOD`.

2. **Conda environment isolation.** All builds (PIO, MUSICA, CheMPAS-A) must
   happen within the same conda environment to ensure consistent
   compiler/library paths.

3. **NoahMP and GSL orography source files.** These are fetched from upstream
   MPAS-Dev/MPAS-Model and are not tracked in the CheMPAS-A repo. If they are
   missing, copy from upstream or the sister repo.

4. **Physics data files.** The `checkout_data_files.sh` script in
   `src/core_atmosphere/physics/` downloads WRF lookup tables from
   [MPAS-Data](https://github.com/MPAS-Dev/MPAS-Data) on first build.

### Environment Variables

The Makefile reads several environment variables:

- `NETCDF` - Path to NetCDF installation
- `PNETCDF` - Path to Parallel NetCDF
- `PIO` - Path to PIO library
- `GPTL` - Path to GPTL profiling library

## Build Pipeline

The CMake build follows this sequence:

1. **Framework compilation** - `MPAS::framework`
2. **Operators compilation** - `MPAS::operators`
3. **Tool generation** - Registry parser, namelist generator
4. **External libraries** - ezxml, SMIOL, ESMF time
5. **Per-core compilation** - Atmosphere, ocean, etc.
6. **Physics libraries** - NoahMP, WRF physics, MMM physics, UGWP
7. **Executable linking** - Merge all targets
8. **Data file installation** - Physics tables, coefficient files

## Dependencies

### Required

| Dependency | Purpose |
|------------|---------|
| MPI | Distributed memory parallelism |
| PnetCDF | Parallel NetCDF I/O |

### Optional

| Dependency | Purpose |
|------------|---------|
| PIO | Parallel I/O abstraction |
| NetCDF | Additional I/O format |
| GPTL | Performance profiling |
| ESMF | Earth System Modeling Framework |
| MUSICA/MICM | Atmospheric chemistry |

## Physics Data Files

Physics parameterizations require lookup tables and coefficient files. The build system can automatically fetch these from the MPAS-Data repository:

```bash
# Data is fetched from:
# https://github.com/MPAS-Dev/MPAS-Data (v8.0 tag)
```

Data categories:
- Microphysics lookup tables (Thompson, WSM6)
- Radiation coefficient files (RRTMG)
- Land surface parameters (Noah, NoahMP)
- Ozone climatology
- Aerosol optical properties

## CMake Modules

Custom CMake modules are in `cmake/`:

- `FindNetCDF.cmake` - Locate NetCDF
- `FindPnetCDF.cmake` - Locate Parallel NetCDF
- `FindPIO.cmake` - Locate PIO
- `MPASFunctions.cmake` - MPAS-specific CMake functions

## Troubleshooting

### Common Issues

**MPI not found:**
```bash
export MPI_HOME=/path/to/mpi
cmake -DMPI_Fortran_COMPILER=mpif90 ..
```

**NetCDF not found:**
```bash
export NETCDF=/path/to/netcdf
cmake -DNETCDF_ROOT=/path/to/netcdf ..
```

**PnetCDF not found:**
```bash
export PNETCDF=/path/to/pnetcdf
cmake -DPNETCDF_ROOT=/path/to/pnetcdf ..
```

### Verbose Build

```bash
make VERBOSE=1
```

### Clean Build

```bash
rm -rf build
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Cross-Compilation

For HPC systems with compute nodes different from login nodes:

```bash
cmake -DCMAKE_SYSTEM_NAME=Linux \
      -DCMAKE_Fortran_COMPILER=ftn \
      -DCMAKE_C_COMPILER=cc \
      ..
```

## Related Documentation

- [ARCHITECTURE.md](docs/chempas/architecture/ARCHITECTURE.md) - System architecture
- [COMPONENTS.md](docs/chempas/architecture/COMPONENTS.md) - Component details
- `INSTALL` - Installation instructions
