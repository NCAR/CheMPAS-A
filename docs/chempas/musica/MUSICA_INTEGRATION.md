# MUSICA/MICM Integration in MPAS

This document describes the integration of the MUSICA (Multi-Scale
Infrastructure for Chemistry and Aerosols) chemistry stack into
MPAS-Atmosphere as implemented in CheMPAS-A.

## Overview

The MUSICA integration enables coupled atmospheric-chemistry modeling on
MPAS's unstructured mesh. CheMPAS-A uses MICM (Model Independent Chemistry
Module) as the chemical ODE solver and TUV-x as the optional photolysis
solver.

## Architecture

```
  MPAS Atmosphere Core
  ====================

  Dynamics --> Physics --> Chemistry
   (SRK3)    (Radiation,   (mpas_atm_chemistry.F)
             Convection)          |
                                  v
                  +------- mpas_musica.F
                  |        =============
                  |        MICM Solver (Rosenbrock)
                  |        State (Coupled + Reference)
                  |
                  +------- mpas_tuvx.F / mpas_solar_geometry.F
                  |        Photolysis rates
                  |
                  +------- mpas_lightning_nox.F
                           Operator-split NO source
```

## Source Files

| File | Location | Purpose |
|------|----------|---------|
| `mpas_atm_chemistry.F` | `src/core_atmosphere/chemistry/` | Top-level chemistry init/step/finalize driver, MPAS gather/scatter, photolysis update gating |
| `mpas_musica.F` | `src/core_atmosphere/chemistry/musica/` | MUSICA/MICM state management, species mapping, unit conversion, photolysis-rate writes |
| `mpas_tuvx.F` | `src/core_atmosphere/chemistry/` | TUV-x setup, host-column profile updates, optional upper-atmosphere extension, cloud optical depth |
| `mpas_solar_geometry.F` | `src/core_atmosphere/chemistry/` | Fallback solar zenith-angle calculation |
| `mpas_lightning_nox.F` | `src/core_atmosphere/chemistry/` | Operator-split lightning NO source |

## Conditional Compilation

The integration is conditionally compiled using the `MPAS_USE_MUSICA` preprocessor flag:

```fortran
#ifdef MPAS_USE_MUSICA
    use mpas_musica, only: musica_init
    ...
#endif
```

Enable at build time with the Makefile workflow used by this repository:
```bash
eval "$(scripts/check_build_env.sh --export)" && make -j8 TARGET \
    CORE=atmosphere MUSICA=true \
    PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```

Replace `TARGET` with the compiler target for the host, typically
`gfortran`, `llvm`, or `ifort`.

## API Overview

### Chemistry Interface (`mpas_atm_chemistry`)

| Routine | Purpose |
|---------|---------|
| `chemistry_init()` | Initialize chemistry packages |
| `chemistry_step()` | Advance chemistry one timestep |
| `chemistry_finalize()` | Clean up resources |
| `chemistry_from_MPAS()` | Extract MPAS state for chemistry |
| `chemistry_to_MPAS()` | Update MPAS state from chemistry |
| `chemistry_seed_chem()` | Seed MPAS chemistry tracers from MICM initial state |
| `chemistry_query_species()` | Query MICM config for runtime chemistry species |

### MUSICA Interface (`mpas_musica`)

| Routine | Purpose |
|---------|---------|
| `musica_init()` | Initialize MICM solver and state |
| `musica_query_species()` | Lightweight species discovery for runtime tracer allocation |
| `musica_step()` | Solve chemistry (coupled state) |
| `musica_step_ref()` | Solve chemistry (reference state) |
| `musica_finalize()` | Clean up MICM resources |
| `MICM_from_chemistry()` | Copy MPAS tracers to MICM |
| `MICM_to_chemistry()` | Copy MICM results to MPAS |
| `resolve_mpas_indices()` | Resolve `index_q*` dimensions for all chemistry species |
| `micm_to_mpas_chem()` | Seed MPAS with MICM initial state (generic species loop) |
| `log_column_comparison()` | Diagnostic logging |
| `copy_state_to_ref()` | Sync reference state |
| `musica_cache_photo_indices()` | Cache MICM `PHOTO.<name>` rate-parameter indices |
| `musica_set_photolysis_rates()` | Write the current photolysis-rate field into MICM state |

## Data Flow

### Each Chemistry Timestep

1. **Photolysis update** (`tuvx_compute_photolysis` or fallback `cos(SZA)`)
   - Compute rate parameters such as `PHOTO.jNO2`
   - Write them to MICM through `musica_set_photolysis_rates`
   - Mirror available rates to diagnostics such as `j_jNO2`

2. **Lightning NOx source** (`lightning_nox_inject`)
   - Adds operator-split NO to `qNO` when the active mechanism contains NO

3. **MPAS -> MICM** (`chemistry_from_MPAS`)
   - Extract scalars (tracers), temperature, pressure, density from MPAS pools
   - Convert mixing ratios [kg/kg] to concentrations [mol/m³]
   - Set MICM environmental conditions

4. **MICM Solve** (`musica_step`)
   - Rosenbrock ODE integration
   - Solve coupled state (advected by MPAS)
   - Optionally solve reference state (chemistry only, for diagnostics)
   - Optionally split the MPAS timestep into `config_chem_substeps` calls

5. **MICM -> MPAS** (`chemistry_to_MPAS`)
   - Convert concentrations [mol/m³] to mixing ratios [kg/kg]
   - Update MPAS scalar tracers

### Unit Conversion

```
MPAS → MICM:  C = q × ρ_air / M_species
MICM → MPAS:  q = C × M_species / ρ_air

where:
  C = concentration [mol/m³]
  q = mixing ratio [kg/kg]
  ρ_air = air density [kg/m³]
  M_species = molar mass [kg/mol]
```

## Mechanism Configuration

The coupling is mechanism-agnostic. Chemistry species are read from the MICM
configuration at runtime, and MPAS tracer names follow the convention:

`MICM species X -> MPAS tracer qX`

The shipped mechanisms include ABBA (`AB`, `A`, `B`), LNOx-O3 (`NO`, `NO2`,
`O3`), Chapman, and Chapman + NOx variants. MICM species map to MPAS tracers
by prefixing `q`, e.g. `NO2 -> qNO2` and `O3 -> qO3`.

Molar masses are read per-species from MICM properties (`__molar mass`) via
`micm%get_species_property_double(...)` during `musica_init`.

## Grid Cell Mapping

MICM processes a 1D array of grid cells. The mapping between MPAS's 2D (cell, level) indexing and MICM's 1D indexing:

```fortran
! MPAS: scalars(tracer, level, cell)
! MICM: state%concentrations(flat_index)

! MICM cell index = (iCell - 1) * nVertLevels + k
micm_cell = (iCell - 1) * nVertLevels + k

! MICM array index (strided)
idx = 1 + (micm_cell - 1) * cell_stride + (species - 1) * var_stride
```

## Dual State Tracking

The integration maintains two MICM states:

1. **Coupled State** (`state`)
   - Updated each timestep from MPAS
   - Experiences advection through MPAS tracer transport
   - Results written back to MPAS

2. **Reference State** (`state_ref`)
   - Chemistry only (no advection)
   - Same initial conditions as coupled state
   - Used to diagnose advection effects

Comparing these states reveals advection contributions to tracer evolution.

## Runtime Tracer Resolution

Chemistry tracers are not statically defined in `Registry.xml`. During
`atm_setup_block`, MPAS queries MICM species and extends `scalars` and
`scalars_tend` metadata dynamically.

During `chemistry_init`, `resolve_mpas_indices()` then resolves each
`index_q*` dimension from the state pool. Missing/invalid indices are treated
as initialization errors.

Runtime chemistry tracers are currently guarded against `config_apply_lbcs=true`
because `lbc_scalars` remains statically sized from registry metadata.

## Configuration

### Namelist Options

```fortran
&musica
    config_micm_file = 'lnox_o3.yaml'
    config_tuvx_config_file = 'tuvx_no2.json'
/
```

### MICM Configuration File

The MICM solver reads its mechanism from the YAML file named by
`config_micm_file`. The repository's current mechanism files live in
`micm_configs/` and specify:
- Chemical species
- Reactions
- Rate constants
- Initial concentrations

## Diagnostic Output

The integration logs diagnostic information:

```
[MUSICA] Initializing MICM chemistry package...
MICM version: X.Y.Z
MICM number of grid cells: 163840

[MUSICA] MICM species: AB
[MUSICA] MICM species: A
[MUSICA] MICM species: B

[MUSICA] Stepping MICM solver...
[MUSICA] MICM Solver statistics ...
  MICM function calls: 5
  MICM jacobian updates: 2
  MICM number of steps: 1
  ...

[COMPARE] Probe cell=512: Coupled vs Reference (chemistry-only)
[COMPARE] Level | AB_coupled AB_ref | A_coupled A_ref | B_coupled B_ref
```

## Dependencies

| Dependency | Module | Purpose |
|------------|--------|---------|
| `musica_micm` | `micm_t` | MICM solver type |
| `musica_micm` | `solver_stats_t` | Solver statistics |
| `musica_micm` | `RosenbrockStandardOrder` | Solver type constant |
| `musica_state` | `state_t` | Chemical state type |
| `musica_util` | `error_t` | Error handling |
| `musica_util` | `string_t` | String utilities |

## Error Handling

Errors from MUSICA are captured and propagated:

```fortran
type(error_t) :: error

call micm%solve(time_step, state, solver_state, solver_stats, error)
if (has_error_occurred(error, error_message, error_code)) return
```

The `has_error_occurred()` helper converts MUSICA errors to MPAS-compatible format.

## Current Status And Follow-On Work

Implemented in this branch:
- Runtime species discovery and MPAS tracer allocation from MICM YAML files
- MPAS <-> MICM state transfer and unit conversion
- Optional reference-state solve for chemistry-only comparison
- LNOx-O3, Chapman, and Chapman + NOx idealized mechanisms
- TUV-x photolysis rates, upper-atmosphere column extension, and cloud optical
  depth from MPAS `qc`/`qr`

Remaining follow-on work:
- Larger atmospheric chemistry mechanisms and aerosol chemistry
- More production-oriented validation beyond the idealized cases
- Performance work for expensive photolysis configurations

## Related Documentation

- [ARCHITECTURE.md](../architecture/ARCHITECTURE.md) - Overall system architecture
- `BUILD.md` - Build configuration for MUSICA
- [COMPONENTS.md](../architecture/COMPONENTS.md) - Atmosphere component details
