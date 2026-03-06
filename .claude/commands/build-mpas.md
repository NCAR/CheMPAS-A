# Build MPAS

Build MPAS-Atmosphere with LLVM compilers on macOS.

## Instructions

1. Set up the environment (use absolute paths — `$HOME` doesn't expand reliably inside make):
   ```bash
   export PKG_CONFIG_PATH="/Users/fillmore/software/lib/pkgconfig:$PKG_CONFIG_PATH"
   ```

2. Check if a clean build is needed by looking for stale `.mod` files or build option mismatches.

3. Build the requested core (default: atmosphere):
   ```bash
   PKG_CONFIG_PATH="/Users/fillmore/software/lib/pkgconfig:$PKG_CONFIG_PATH" \
   make -j8 llvm \
     CORE=atmosphere \
     PIO=/Users/fillmore/software \
     NETCDF=/opt/homebrew \
     PNETCDF=/Users/fillmore/software \
     PRECISION=double
   ```

4. To enable MUSICA chemistry, add `MUSICA=true` to the build command.

5. Verify the executable was created:
   ```bash
   ls -la atmosphere_model
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
