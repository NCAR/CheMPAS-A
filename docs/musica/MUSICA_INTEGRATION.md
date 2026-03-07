# MUSICA/MICM Integration in MPAS

This document describes the integration of the MUSICA (Multi-Scale Infrastructure for Chemistry and Aerosols) chemistry package into MPAS-Atmosphere. This is the primary focus of the `musica_micm_state` development branch.

## Overview

The MUSICA integration enables coupled atmospheric-chemistry modeling on MPAS's unstructured mesh. The integration uses MICM (Model Independent Chemistry Module) as the ODE solver for chemical kinetics.

## Architecture

```
  MPAS Atmosphere Core
  ====================

  Dynamics --> Physics --> Chemistry
   (SRK3)    (Radiation,   (mpas_atm_chemistry.F)
             Convection)          |
                                  v
                           mpas_musica.F
                           =============
                           MICM Solver (Rosenbrock)
                           State (Coupled + Reference)
```

## Source Files

| File | Location | Purpose |
|------|----------|---------|
| `mpas_atm_chemistry.F` | `src/core_atmosphere/chemistry/` | Generic chemistry interface |
| `mpas_musica.F` | `src/core_atmosphere/chemistry/musica/` | MUSICA/MICM-specific implementation |

## Conditional Compilation

The integration is conditionally compiled using the `MPAS_USE_MUSICA` preprocessor flag:

```fortran
#ifdef MPAS_USE_MUSICA
    use mpas_musica, only: musica_init
    ...
#endif
```

Enable at build time:
```bash
cmake -DMPAS_USE_MUSICA=ON -DMUSICA_ROOT=/path/to/musica ..
```

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

## Data Flow

### Each Chemistry Timestep

1. **MPAS → MICM** (`chemistry_from_MPAS`)
   - Extract scalars (tracers), temperature, pressure, density from MPAS pools
   - Convert mixing ratios [kg/kg] to concentrations [mol/m³]
   - Set MICM environmental conditions

2. **MICM Solve** (`musica_step`)
   - Rosenbrock ODE integration
   - Solve coupled state (advected by MPAS)
   - Solve reference state (chemistry only, for diagnostics)

3. **MICM → MPAS** (`chemistry_to_MPAS`)
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

Regression testing currently uses the ABBA mechanism (`AB`, `A`, `B`), which
maps to tracers `qAB`, `qA`, and `qB`.

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
config_micm_file = 'path/to/micm_config.json'
```

### MICM Configuration File

The MICM solver reads its mechanism from a JSON configuration file specifying:
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

## Future Development

Current status (as of this branch):
- Basic MPAS↔MICM coupling implemented
- ABBA test chemistry working
- Dual state tracking for advection diagnostics

Planned enhancements:
- Full atmospheric chemistry mechanisms (e.g., tropospheric ozone)
- TUV-x photolysis rate calculations
- Aerosol chemistry
- Parallel processing optimization

## Related Documentation

- [ARCHITECTURE.md](../architecture/ARCHITECTURE.md) - Overall system architecture
- [BUILD.md](../../BUILD.md) - Build configuration for MUSICA
- [COMPONENTS.md](../architecture/COMPONENTS.md) - Atmosphere component details
