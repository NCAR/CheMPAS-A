# MUSICA Development

Work on MUSICA/MICM atmospheric chemistry integration with MPAS.

## Key Files

| File | Purpose |
|------|---------|
| `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` | Generic chemistry interface |
| `src/core_atmosphere/chemistry/musica/mpas_musica.F` | MUSICA/MICM coupler |
| `src/core_atmosphere/Registry.xml` | Core metadata (non-chemistry tracers); MUSICA tracers are runtime-injected |

## Build with MUSICA

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

## MUSICA Requirements

- MUSICA-Fortran must be built with flang (not gfortran)
- The `musica-fortran.pc` file may need yaml-cpp path: `-L/opt/homebrew/opt/yaml-cpp/lib`
- Current version: MUSICA-Fortran 0.13.0, MICM 3.10.0

## State Coupling

The MPAS-MICM state coupler (`mpas_musica.F`) handles:
- Mapping MPAS atmospheric state to MICM species
- Converting units between MPAS and MICM conventions
- Managing chemistry timesteps within physics driver

## Testing Chemistry

1. Build with `MUSICA=true`
2. Verify MUSICA symbols in executable: `nm atmosphere_model | grep -i musica`
3. Configure chemistry options in `namelist.atmosphere`
4. Check `log.atmosphere.*.out` for MICM initialization messages

## Documentation

- `MUSICA_INTEGRATION.md` - Integration architecture
- `MUSICA_API.md` - Fortran API reference
- `~/Documents/MPAS-MUSICA/` - Downloaded MPAS documentation PDFs
