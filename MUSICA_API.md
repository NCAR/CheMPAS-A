# MUSICA API Reference

This document provides a reference for the MUSICA (Multi-Scale Infrastructure for Chemistry and Aerosols) Fortran API used in MPAS-Atmosphere. The MUSICA library is located at `~/EarthSystem/MUSICA-LLVM`.

## Table of Contents

1. [Module Overview](#module-overview)
2. [MICM Solver API](#micm-solver-api)
3. [State Management API](#state-management-api)
4. [Utility Types](#utility-types)
5. [Error Handling](#error-handling)
6. [Configuration Files](#configuration-files)
7. [Solver Types](#solver-types)
8. [Usage Patterns](#usage-patterns)
9. [Unit Conversions](#unit-conversions)

---

## Module Overview

### Required Modules

```fortran
use musica_micm   ! MICM solver interface
use musica_state  ! State management
use musica_util   ! Utility types (error_t, string_t, mappings_t)
```

### Key Types

| Type | Module | Purpose |
|------|--------|---------|
| `micm_t` | `musica_micm` | MICM solver instance |
| `state_t` | `musica_state` | Chemical state container |
| `conditions_t` | `musica_state` | Environmental conditions (T, P, rho) |
| `solver_stats_t` | `musica_micm` | Solver performance statistics |
| `error_t` | `musica_util` | Error handling |
| `string_t` | `musica_util` | String wrapper |
| `mappings_t` | `musica_util` | Name-to-index mappings |
| `strides_t` | `musica_state` | Array stride information |

---

## MICM Solver API

### `micm_t` Type

```fortran
type :: micm_t
contains
  procedure :: solve
  procedure :: get_state
  procedure :: get_maximum_number_of_grid_cells
  procedure :: get_species_property_string
  procedure :: get_species_property_double
  procedure :: get_species_property_int
  procedure :: get_species_property_bool
end type micm_t
```

### Constructor

```fortran
type(micm_t), pointer :: micm
type(error_t) :: error

micm => micm_t(config_path, solver_type, error)
```

**Parameters:**
- `config_path` (character) - Path to configuration file or directory
- `solver_type` (integer) - Solver type constant (see [Solver Types](#solver-types))
- `error` (error_t, inout) - Error status

### solve()

Integrates the chemical system forward in time.

```fortran
call micm%solve(time_step, state, solver_state, solver_stats, error)
```

**Parameters:**
- `time_step` (real64, in) - Time step in seconds
- `state` (state_t, inout) - Chemical state
- `solver_state` (string_t, out) - Solver status ("Converged" or error)
- `solver_stats` (solver_stats_t, out) - Performance statistics
- `error` (error_t, inout) - Error status

### get_state()

Creates a new chemical state for the specified number of grid cells.

```fortran
type(state_t), pointer :: state
state => micm%get_state(number_of_grid_cells, error)
```

**Parameters:**
- `number_of_grid_cells` (integer, in) - Number of grid cells
- `error` (error_t, inout) - Error status

**Returns:** Pointer to new `state_t` instance

### get_species_property_*()

Query species properties from configuration.

```fortran
! String property
character(len=:), allocatable :: name
name = micm%get_species_property_string(species_name, property_name, error)

! Double property
real(real64) :: value
value = micm%get_species_property_double(species_name, property_name, error)

! Integer property
integer :: ivalue
ivalue = micm%get_species_property_int(species_name, property_name, error)

! Boolean property
logical :: flag
flag = micm%get_species_property_bool(species_name, property_name, error)
```

**Common Property Names:**
- `"__molar mass"` - Molar mass (used by CheMPAS runtime coupling)
- `"molecular weight [kg mol-1]"` - Legacy/alternate molar-mass key in some MICM examples
- `"__long name"` - Full species name
- `"__atoms"` - Number of atoms
- `"__do advect"` - Advection flag
- `"__initial concentration"` - Initial concentration
- `"__absolute tolerance"` - Solver tolerance

### get_maximum_number_of_grid_cells()

```fortran
integer :: max_cells
max_cells = micm%get_maximum_number_of_grid_cells()
```

**Returns:** Maximum grid cells supported by solver type
- Vector-ordered solvers: Limited by SIMD width
- Standard-ordered solvers: Essentially unlimited (>10^8)

### Version Information

```fortran
use musica_micm, only: get_micm_version
type(string_t) :: version
version = get_micm_version()
print *, version%value_
```

### CUDA Availability

```fortran
use musica_micm, only: is_cuda_available
logical :: cuda_ok
cuda_ok = is_cuda_available(error)
```

---

## State Management API

### `state_t` Type

```fortran
type :: state_t
  type(conditions_t), pointer :: conditions(:)     ! Environmental conditions
  real(real64), pointer :: concentrations(:)       ! Species concentrations
  real(real64), pointer :: rate_parameters(:)      ! User-defined rate constants
  type(mappings_t), pointer :: species_ordering    ! Species name→index
  type(mappings_t), pointer :: rate_parameters_ordering
  integer :: number_of_grid_cells
  integer :: number_of_species
  integer :: number_of_rate_parameters
  type(strides_t) :: species_strides
  type(strides_t) :: rate_parameters_strides
contains
  procedure :: update_references
end type state_t
```

### `conditions_t` Type

Environmental conditions for each grid cell.

```fortran
type, bind(c) :: conditions_t
  real(c_double) :: temperature   ! [K]
  real(c_double) :: pressure      ! [Pa]
  real(c_double) :: air_density   ! [mol/m³]
end type conditions_t
```

### `strides_t` Type

Array stride information for multi-cell access.

```fortran
type :: strides_t
  integer :: grid_cell  ! Stride between grid cells
  integer :: variable   ! Stride between variables (species/params)
end type strides_t
```

### Accessing Concentrations

Concentrations use strided array access:

```fortran
integer :: cell_stride, var_stride, idx, species_idx, cell_idx
real(real64) :: concentration

cell_stride = state%species_strides%grid_cell
var_stride = state%species_strides%variable

! Get species index by name
species_idx = state%species_ordering%index("O3", error)

! Compute array index for cell_idx, species_idx
idx = 1 + (cell_idx - 1) * cell_stride + (species_idx - 1) * var_stride

! Read or write concentration
concentration = state%concentrations(idx)
state%concentrations(idx) = new_value
```

### Accessing Rate Parameters

Same strided pattern as concentrations:

```fortran
integer :: param_idx, idx
real(real64) :: rate

param_idx = state%rate_parameters_ordering%index("PHOTO.jO3", error)

cell_stride = state%rate_parameters_strides%grid_cell
var_stride = state%rate_parameters_strides%variable

idx = 1 + (cell_idx - 1) * cell_stride + (param_idx - 1) * var_stride
state%rate_parameters(idx) = rate_value
```

### Setting Environmental Conditions

```fortran
do i = 1, state%number_of_grid_cells
  state%conditions(i)%temperature = T(i)      ! [K]
  state%conditions(i)%pressure = P(i)         ! [Pa]
  state%conditions(i)%air_density = rho(i)    ! [mol/m³]
end do
```

### update_references()

**Important:** C++ may reallocate arrays during solve. Call after `micm%solve()`:

```fortran
call state%update_references(error)
```

### Cleanup

```fortran
deallocate(state)  ! Calls finalizer automatically
deallocate(micm)
```

---

## Utility Types

### `mappings_t` - Name/Index Mappings

```fortran
type :: mappings_t
contains
  procedure :: name   ! Get name by position
  procedure :: index  ! Get index by name or position
  procedure :: size   ! Number of entries
end type mappings_t
```

**Usage:**

```fortran
integer :: n, i, idx
character(len=:), allocatable :: species_name

! Number of species
n = state%species_ordering%size()

! Iterate over all species
do i = 1, n
  species_name = state%species_ordering%name(i)
  idx = state%species_ordering%index(i)
  print *, "Species ", i, ": ", species_name, " at index ", idx
end do

! Lookup by name
idx = state%species_ordering%index("O3", error)
```

### `string_t` - String Wrapper

```fortran
type :: string_t
  character(len=:), allocatable :: value_
contains
  procedure :: get_char_array
end type string_t
```

**Usage:**

```fortran
type(string_t) :: solver_state
character(len=:), allocatable :: status

! After solve
if (solver_state%get_char_array() == "Converged") then
  ! Success
end if

! Or access directly
status = solver_state%value_
```

### `solver_stats_t` - Solver Statistics

```fortran
type :: solver_stats_t
contains
  procedure :: function_calls     ! Number of RHS evaluations
  procedure :: jacobian_updates   ! Jacobian recalculations
  procedure :: number_of_steps    ! Integration steps taken
  procedure :: accepted           ! Accepted steps
  procedure :: rejected           ! Rejected steps
  procedure :: decompositions     ! LU decompositions
  procedure :: solves             ! Linear system solves
  procedure :: final_time         ! Actual integration time [s]
end type solver_stats_t
```

**Usage:**

```fortran
call micm%solve(dt, state, solver_state, stats, error)

print *, "Function calls: ", stats%function_calls()
print *, "Steps: ", stats%number_of_steps()
print *, "Accepted/Rejected: ", stats%accepted(), "/", stats%rejected()
print *, "Final time: ", stats%final_time()
```

---

## Error Handling

### `error_t` Type

```fortran
type :: error_t
contains
  procedure :: code        ! Error code (integer)
  procedure :: category    ! Error category (string)
  procedure :: message     ! Error message (string)
  procedure :: is_success  ! True if no error
  procedure :: is_error    ! Check for specific error
end type error_t
```

### Error Codes

**MUSICA Error Category:**
| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Species not found |
| 2 | Solver type not found |
| 3 | Mapping not found |
| 4 | Mapping options undefined |
| 5 | Unknown error |
| 6 | Unsupported solver/state pair |

**MUSICA Parsing Category:**
| Code | Meaning |
|------|---------|
| 1 | Parsing failed |
| 2 | Invalid config file |
| 3 | Unsupported version |
| 4 | Failed to cast to version |

### Usage Pattern

```fortran
type(error_t) :: error

call micm%solve(dt, state, solver_state, stats, error)

if (.not. error%is_success()) then
  print *, "Error category: ", error%category()
  print *, "Error code: ", error%code()
  print *, "Message: ", error%message()
  return
end if
```

### Recommended Helper Function

```fortran
logical function has_error_occurred(error, error_message, error_code)
  use musica_util, only: error_t
  type(error_t), intent(in) :: error
  character(len=:), allocatable, intent(out) :: error_message
  integer, intent(out) :: error_code
  character(len=30) :: code_str

  if (error%is_success()) then
    error_code = 0
    error_message = ''
    has_error_occurred = .false.
  else
    error_code = error%code()
    write(code_str, '(I30)') error_code
    error_message = '[MUSICA Error]: ' // error%category() // &
                    '[' // trim(adjustl(code_str)) // ']: ' // error%message()
    has_error_occurred = .true.
  end if
end function
```

---

## Configuration Files

### v0 Format (Legacy - Directory)

Directory containing multiple JSON files:

```
config_directory/
├── config.json      # Points to other files
├── species.json     # Species definitions
└── reactions.json   # Reaction mechanisms
```

**config.json:**
```json
{
  "camp-files": ["species.json", "reactions.json"]
}
```

**species.json:**
```json
{
  "camp-data": [
    {
      "name": "O3",
      "type": "CHEM_SPEC",
      "molecular weight [kg mol-1]": 0.048,
      "__long name": "ozone",
      "__do advect": true,
      "__initial concentration": 8.1e-6
    }
  ]
}
```

**reactions.json:**
```json
{
  "camp-data": [
    {
      "type": "ARRHENIUS",
      "A": 2.9e19,
      "reactants": {"O": {"qty": 1}, "O2": {"qty": 1}},
      "products": {"O3": {"qty": 1}}
    },
    {
      "type": "PHOTOLYSIS",
      "name": "jO2",
      "reactants": [{"species name": "O2"}],
      "products": [{"species name": "O", "coefficient": 2.0}]
    }
  ]
}
```

**Path:** `micm => micm_t("configs/v0/chapman", solver_type, error)`

### v1 Format (Modern - Single File)

Single JSON file with all definitions:

```json
{
  "version": "1.0.0",
  "name": "Chapman",
  "species": [
    {
      "name": "O3",
      "molecular weight [kg mol-1]": 0.048,
      "__initial concentration": 8.1e-6
    }
  ],
  "phases": [
    {"name": "gas", "species": ["O2", "O", "O3"]}
  ],
  "reactions": [
    {
      "type": "ARRHENIUS",
      "A": 2.9e19,
      "reactants": [{"species name": "O"}, {"species name": "O2"}],
      "products": [{"species name": "O3"}]
    }
  ]
}
```

**Path:** `micm => micm_t("configs/v1/chapman/config.json", solver_type, error)`

### Reaction Types

| Type | Description | Rate Source |
|------|-------------|-------------|
| `ARRHENIUS` | Temperature-dependent | Computed from A, B, C, D, E params |
| `PHOTOLYSIS` | Light-dependent | Via `rate_parameters` array |
| `USER_DEFINED` | External rate | Via `rate_parameters` array |
| `TROE` | Pressure-dependent | Computed internally |
| `TERNARY_CHEMICAL_ACTIVATION` | Three-body | Computed internally |

---

## Solver Types

```fortran
use musica_micm, only: Rosenbrock, RosenbrockStandardOrder, &
                       BackwardEuler, BackwardEulerStandardOrder, &
                       CudaRosenbrock
```

| Solver | Description | Max Grid Cells |
|--------|-------------|----------------|
| `Rosenbrock` | Vectorized Rosenbrock | SIMD-limited (~100-1000) |
| `RosenbrockStandardOrder` | Non-vectorized Rosenbrock | Unlimited (>10^8) |
| `BackwardEuler` | Vectorized implicit Euler | SIMD-limited |
| `BackwardEulerStandardOrder` | Non-vectorized implicit Euler | Unlimited |
| `CudaRosenbrock` | GPU-accelerated | GPU memory limited |

**Recommendation for MPAS:** Use `RosenbrockStandardOrder` for unlimited grid cells.

---

## Usage Patterns

### Complete Integration Example

```fortran
subroutine run_chemistry(dt, nCells, nLevels, T, P, rho, tracers)
  use iso_fortran_env, only: real64
  use musica_micm, only: micm_t, solver_stats_t, RosenbrockStandardOrder
  use musica_state, only: state_t
  use musica_util, only: error_t, string_t

  real(real64), intent(in) :: dt
  integer, intent(in) :: nCells, nLevels
  real(real64), intent(in) :: T(:,:), P(:,:), rho(:,:)
  real(real64), intent(inout) :: tracers(:,:,:)

  type(micm_t), pointer :: micm
  type(state_t), pointer :: state
  type(error_t) :: error
  type(string_t) :: solver_state
  type(solver_stats_t) :: stats

  integer :: total_cells, spec_idx, cell_stride, var_stride
  integer :: iCell, k, micm_cell, idx

  total_cells = nCells * nLevels

  ! Initialize
  micm => micm_t("micm_config.json", RosenbrockStandardOrder, error)
  state => micm%get_state(total_cells, error)

  ! Get strides and species index
  cell_stride = state%species_strides%grid_cell
  var_stride = state%species_strides%variable
  spec_idx = state%species_ordering%index("O3", error)

  ! Copy data to MICM state
  do iCell = 1, nCells
    do k = 1, nLevels
      micm_cell = (iCell - 1) * nLevels + k

      ! Set conditions
      state%conditions(micm_cell)%temperature = T(k, iCell)
      state%conditions(micm_cell)%pressure = P(k, iCell)
      state%conditions(micm_cell)%air_density = rho(k, iCell) / 0.0289644_real64

      ! Set concentration (convert kg/kg to mol/m³)
      idx = 1 + (micm_cell - 1) * cell_stride + (spec_idx - 1) * var_stride
      state%concentrations(idx) = tracers(spec_idx, k, iCell) * rho(k, iCell) / 0.048_real64
    end do
  end do

  ! Solve
  call micm%solve(dt, state, solver_state, stats, error)

  ! Copy results back (convert mol/m³ to kg/kg)
  do iCell = 1, nCells
    do k = 1, nLevels
      micm_cell = (iCell - 1) * nLevels + k
      idx = 1 + (micm_cell - 1) * cell_stride + (spec_idx - 1) * var_stride
      tracers(spec_idx, k, iCell) = state%concentrations(idx) * 0.048_real64 / rho(k, iCell)
    end do
  end do

  ! Cleanup
  deallocate(state)
  deallocate(micm)

end subroutine
```

### Grid Cell Indexing (MPAS to MICM)

```fortran
! MPAS: 2D grid (iCell, k) where iCell=1..nCells, k=1..nVertLevels
! MICM: 1D array with nCells*nVertLevels elements

! Forward mapping
micm_cell = (iCell - 1) * nVertLevels + k

! Reverse mapping
iCell = (micm_cell - 1) / nVertLevels + 1
k = mod(micm_cell - 1, nVertLevels) + 1
```

### Setting Rate Parameters (Photolysis)

```fortran
integer :: jO3_idx, cell_stride, var_stride, idx
real(real64) :: j_rate

jO3_idx = state%rate_parameters_ordering%index("PHOTO.jO3", error)
cell_stride = state%rate_parameters_strides%grid_cell
var_stride = state%rate_parameters_strides%variable

do micm_cell = 1, total_cells
  idx = 1 + (micm_cell - 1) * cell_stride + (jO3_idx - 1) * var_stride
  state%rate_parameters(idx) = j_rate(micm_cell)
end do
```

---

## Unit Conversions

### MPAS to MICM

```fortran
! Mixing ratio [kg/kg] to concentration [mol/m³]
! C = q × rho_air / M_species

concentration = mixing_ratio * air_density_kg_m3 / molar_mass_kg_mol
```

### MICM to MPAS

```fortran
! Concentration [mol/m³] to mixing ratio [kg/kg]
! q = C × M_species / rho_air

mixing_ratio = concentration * molar_mass_kg_mol / air_density_kg_m3
```

### Air Density Conversion

```fortran
! MPAS provides: rho_air [kg/m³]
! MICM expects: air_density [mol/m³]

real(real64), parameter :: M_AIR = 0.0289644_real64  ! kg/mol

state%conditions(i)%air_density = rho_air_kg_m3 / M_AIR
```

### Common Molar Masses

| Species | Molar Mass [kg/mol] |
|---------|---------------------|
| Dry air | 0.0289644 |
| O₃ | 0.048 |
| O₂ | 0.032 |
| N₂ | 0.028 |
| H₂O | 0.018 |
| CO₂ | 0.044 |

---

## Related Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - MPAS system architecture
- [MUSICA_INTEGRATION.md](MUSICA_INTEGRATION.md) - MPAS-MUSICA integration details
- [MUSICA-LLVM Source](~/EarthSystem/MUSICA-LLVM) - Full library source code
