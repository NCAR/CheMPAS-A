# Build MPAS

Build MPAS-Atmosphere or init_atmosphere. The preflight script auto-detects
the compiler toolchain (LLVM/flang on macOS, GCC/gfortran on Ubuntu) and
sets paths accordingly. See BUILD.md for full details.

## Instructions

1. Check if a clean build is needed by looking for stale `.mod` files or build option mismatches.

2. Run the repo build preflight. It detects the compiler toolchain, library
   paths, and prints the exports needed by `make`.

   ```bash
   scripts/check_build_env.sh
   ```

   The script will report the detected make target (`llvm` or `gfortran`)
   and all resolved paths.

3. Build with the preflight exports in a single shell invocation:

   ```bash
   eval "$(scripts/check_build_env.sh --export)" && make -j8 TARGET \
     CORE=atmosphere \
     PIO="$PIO" \
     NETCDF="$NETCDF" \
     PNETCDF="$PNETCDF" \
     PRECISION=double
   ```

   Replace `TARGET` with the make target reported by the preflight script:
   - macOS: `llvm`
   - Ubuntu: `gfortran`

4. To enable MUSICA chemistry, add `MUSICA=true`:

   ```bash
   eval "$(scripts/check_build_env.sh --export)" && make -j8 TARGET \
     CORE=atmosphere \
     PIO="$PIO" \
     NETCDF="$NETCDF" \
     PNETCDF="$PNETCDF" \
     PRECISION=double \
     MUSICA=true
   ```

5. Verify the executable was created:
   ```bash
   ls -la atmosphere_model
   ```

## Build Environment

The preflight script resolves paths per platform:

**macOS (LLVM):**
- `NETCDF=/opt/homebrew`
- `PNETCDF=$HOME/software`, `PIO=$HOME/software`
- `OMPI_FC=flang`, `OMPI_CC=clang`, `OMPI_CXX=clang++`

**Ubuntu (GCC/conda):**
- `NETCDF=$CONDA_PREFIX` (miniconda3/envs/mpas)
- `PNETCDF=$CONDA_PREFIX`, `PIO=$HOME/software`
- `OMPI_FC=gfortran`, `OMPI_CC=gcc`, `OMPI_CXX=g++`

The `eval ... && make ...` **must be in the same shell invocation**. GNU Make
evaluates `$(shell pkg-config ...)` while parsing the Makefile, so `make` must
inherit the correct `PKG_CONFIG_PATH` directly.

On Ubuntu, the `mpas` conda environment must be active (`conda activate mpas`).

The installed MUSICA package in `~/software/lib/pkgconfig/musica-fortran.pc`
is the preferred source of truth on both platforms.

## Options

- `CORE=atmosphere` (default) or `CORE=init_atmosphere`
- `MUSICA=true` to enable MUSICA/MICM chemistry
- `AUTOCLEAN=true` to automatically clean incompatible builds
- `DEBUG=true` for debug build

## Clean Build

If needed, perform a full clean:
```bash
make clean CORE=atmosphere
find . -name "*.mod" -delete
find . -name "*.o" -delete
eval "$(scripts/check_build_env.sh --export)" && make -j8 TARGET CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```
