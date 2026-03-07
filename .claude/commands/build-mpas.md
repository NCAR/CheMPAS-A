# Build MPAS

Build MPAS-Atmosphere with LLVM compilers on macOS. See BUILD.md for full details.

## Instructions

1. Check if a clean build is needed by looking for stale `.mod` files or build option mismatches.

2. Run the repo build preflight first. It detects the working local LLVM/MUSICA/PNetCDF
   environment and prints the exports needed by `make`.

   ```bash
   scripts/check_build_env.sh
   eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
     CORE=atmosphere \
     PIO="$PIO" \
     NETCDF="$NETCDF" \
     PNETCDF="$PNETCDF" \
     PRECISION=double
   ```

3. To enable MUSICA chemistry, add `MUSICA=true`:

   ```bash
   eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
     CORE=atmosphere \
     PIO="$PIO" \
     NETCDF="$NETCDF" \
     PNETCDF="$PNETCDF" \
     PRECISION=double \
     MUSICA=true
   ```

4. Verify the executable was created:
   ```bash
   ls -la atmosphere_model
   ```

## Important: Build Environment

The preflight script currently resolves the working local environment to:
- `NETCDF=/opt/homebrew`
- `PNETCDF=/Users/fillmore/software`
- `PIO=/Users/fillmore/software`
- `PKG_CONFIG_PATH=/Users/fillmore/software/lib/pkgconfig`

The `eval ... && make ...` **must be in the same shell invocation**. GNU Make evaluates
`$(shell pkg-config ...)` while parsing the Makefile, so `make` must inherit the correct
`PKG_CONFIG_PATH` directly.

```bash
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double MUSICA=true
```

Use the installed MUSICA package in `/Users/fillmore/software/lib/pkgconfig/musica-fortran.pc`.
Do not prefer the sibling `~/EarthSystem/MUSICA-LLVM/build` tree unless the installed package
is missing or broken.

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
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```
