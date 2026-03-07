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
scripts/check_build_env.sh
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
  CORE=atmosphere \
  PIO="$PIO" \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double \
  MUSICA=true
```

## MUSICA Requirements

- MUSICA-Fortran must be built with flang (not gfortran)
- Use `scripts/check_build_env.sh` before building
- The preferred pkg-config file is `/Users/fillmore/software/lib/pkgconfig/musica-fortran.pc`
- `make` must inherit `PKG_CONFIG_PATH` in the same shell invocation
- Current version: MUSICA-Fortran 0.13.0, MICM 3.10.0
- Full atmosphere builds may still stop later if `physics_mmm` tries to fetch `MMM-physics`
  and network access is unavailable

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
