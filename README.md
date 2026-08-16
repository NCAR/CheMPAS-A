# CheMPAS-A

CheMPAS-A extends the MPAS-Atmosphere core with online, mechanism-configurable
atmospheric chemistry. It couples MPAS transport and meteorology to
[MUSICA](https://github.com/NCAR/MUSICA), using MICM for gas-phase chemistry,
TUV-x for photolysis, and MIEM for offline emissions.

The [CheMPAS-A wiki](https://github.com/NCAR/CheMPAS-A/wiki) contains build and
configuration references, worked cases, downloadable namelists and mechanism
files, and the data/provenance contract for reconstructing the example runs.

## Features

- Runtime discovery and registration of MICM species, transported tracers,
  photolysis diagnostics, and MIEM diagnostics.
- Dry-air mass mixing ratio as the authoritative transport and restart state,
  with optional output-only `fraction`, `percent`, `ppmv`, `ppbv`, and `pptv`
  diagnostics.
- Configurable chemistry cadence, MICM substeps and tolerance, reference solves,
  and bounded recovery from solver failures.
- TUV-x photolysis driven by MPAS columns and solar geometry, with static or
  exact-grid climatological upper-column extensions. A cosine-of-solar-zenith
  fallback is available for `jNO2` when TUV-x is disabled.
- MIEM surface and layered emissions on the exact MPAS mesh, including multiple
  inventories, signed exchange for selected species, and bounded sector,
  category, and layer diagnostics.
- Parameterized lightning NOx with altitude/vertical-velocity or isotherm
  gating.
- MPI-aware initialization, exact-grid validation, mass bookkeeping, and
  restart-aware chemistry state handling.

The implementation is an MVP for chemistry coupling and process integration.
The supplied mechanisms and initial conditions are research examples; they are
not a production air-quality forecast or a validated global composition
configuration.

## Build overview

CheMPAS-A uses the normal MPAS dependencies (MPI, NetCDF-C, NetCDF-Fortran,
PNetCDF, and PIO) plus MUSICA-Fortran 0.16.5 at source revision
`1403e3d22717bc87f3bf9d0aa591caf039c92bbc`, built with MIEM enabled. The
Makefile checks the MUSICA version/revision and MIEM feature metadata before
compilation.

After setting installation prefixes for the dependencies, a GNU build is:

```bash
export NETCDF=/path/to/netcdf
export NETCDFF=/path/to/netcdf-fortran
export PNETCDF=/path/to/pnetcdf
export PIO=/path/to/pio
export MUSICA_PREFIX=/path/to/musica-install
export PKG_CONFIG_PATH="${MUSICA_PREFIX}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"

make clean CORE=atmosphere
make -j8 gfortran CORE=atmosphere OPENMP=false \
  PIO="$PIO" NETCDF="$NETCDF" NETCDFF="$NETCDFF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true
```

This produces `atmosphere_model`. The wiki's
[building guide](https://github.com/NCAR/CheMPAS-A/wiki/Building) covers the
pinned MUSICA build, `init_atmosphere_model`, ABI checks, and the tested Ubuntu
toolchain.

## Run and configuration guides

- [Getting started](https://github.com/NCAR/CheMPAS-A/wiki/Getting-Started)
- [Features and architecture](https://github.com/NCAR/CheMPAS-A/wiki/Features-and-Architecture)
- [Namelist reference](https://github.com/NCAR/CheMPAS-A/wiki/Configuration-Reference)
- [Reconstructable examples](https://github.com/NCAR/CheMPAS-A/wiki/Examples)
- [Global chemistry and emissions MVP](https://github.com/NCAR/CheMPAS-A/wiki/Global-Chemistry-and-Emissions)
- [Data and provenance](https://github.com/NCAR/CheMPAS-A/wiki/Data-and-Provenance)

Chemistry is enabled only in a `MUSICA=true` atmosphere build. An empty
`config_micm_file` leaves the chemistry coupling inactive. Chemistry currently
requires one MPAS block per MPI task; multi-block chemistry is not supported.

## MPAS

CheMPAS-A retains the broader MPAS repository layout and upstream model cores.
General MPAS-Atmosphere concepts, mesh resources, and user guidance are
available from the [MPAS project](https://www2.mmm.ucar.edu/projects/mpas/).

## Contributing and attribution

Before contributing or building on CheMPAS-A, read the
[contributor and attribution guide](CONTRIBUTING.md) for licensing,
provenance, credit, authorship, and agentic-tool expectations.

See [LICENSE](LICENSE) for the repository license.
