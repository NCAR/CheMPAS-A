# Build MPAS

Build MPAS-Atmosphere with LLVM compilers on macOS. See BUILD.md for full details.

## Instructions

1. Check if a clean build is needed by looking for stale `.mod` files or build option mismatches.

2. Build the requested core (default: atmosphere). **Critical:** `PKG_CONFIG_PATH` must be
   exported in the shell environment before invoking `make`, because GNU Make `$(shell ...)`
   directives are evaluated at parse time and inherit the parent process environment.

   ```bash
   export PKG_CONFIG_PATH="/Users/fillmore/software/lib/pkgconfig:$PKG_CONFIG_PATH"
   make -j8 llvm \
     CORE=atmosphere \
     PIO=/Users/fillmore/software \
     NETCDF=/opt/homebrew \
     PNETCDF=/Users/fillmore/software \
     PRECISION=double
   ```

3. To enable MUSICA chemistry, add `MUSICA=true`:

   ```bash
   export PKG_CONFIG_PATH="/Users/fillmore/software/lib/pkgconfig:$PKG_CONFIG_PATH"
   make -j8 llvm \
     CORE=atmosphere \
     PIO=/Users/fillmore/software \
     NETCDF=/opt/homebrew \
     PNETCDF=/Users/fillmore/software \
     PRECISION=double \
     MUSICA=true
   ```

4. Verify the executable was created:
   ```bash
   ls -la atmosphere_model
   ```

## Important: PKG_CONFIG_PATH

The `export` and `make` **must be in the same shell invocation** (same Bash tool call,
joined with `&&`). If they are separate Bash calls, the export is lost because each
Bash tool call runs in a fresh shell.

```bash
export PKG_CONFIG_PATH="/Users/fillmore/software/lib/pkgconfig:$PKG_CONFIG_PATH" && make -j8 llvm CORE=atmosphere PIO=/Users/fillmore/software NETCDF=/opt/homebrew PNETCDF=/Users/fillmore/software PRECISION=double MUSICA=true
```

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
```
