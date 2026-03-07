# CheMPAS Components

This document provides detailed information about the major components of
CheMPAS, with focus on the atmosphere core and the chemistry extensions in this
repository.

## Core Atmosphere (`src/core_atmosphere/`)

The atmosphere core is the primary component of this repository.

### Directory Structure

```
src/core_atmosphere/
├── Registry.xml                    # Main metadata and code-generation source
├── mpas_atm_core.F                 # Core main module
├── mpas_atm_core_interface.F       # Core interface (driver integration)
├── mpas_atm_dimensions.F           # Dimension definitions
├── mpas_atm_halos.F                # Halo exchange definitions
├── mpas_atm_threading.F            # Threading setup
├── CMakeLists.txt                  # Build configuration
├── dynamics/                       # Time stepping & dynamics
├── physics/                        # Physics parameterizations
├── chemistry/                      # Chemistry driver, MUSICA/MICM, TUV-x
├── diagnostics/                    # Post-processing diagnostics
└── utils/                          # Utility programs
```

### Dynamics (`dynamics/`)

The dynamics component handles atmospheric motion and time integration.

| File | Size | Purpose |
|------|------|---------|
| `mpas_atm_time_integration.F` | 360KB | Main dynamics solver (SRK3 scheme) |
| `mpas_atm_boundaries.F` | - | Boundary condition handling |
| `mpas_atm_iau.F` | - | Incremental Analysis Update |

**Time Integration Scheme:**
- Split-explicit Runge-Kutta 3rd order (SRK3)
- Non-hydrostatic equations
- Acoustic-mode substeps for stability
- Horizontal advection with monotonic limiters

### Physics (`physics/`)

Physics parameterizations from multiple sources (~50 modules).

#### Main Physics Drivers

| Driver | Purpose |
|--------|---------|
| `mpas_atmphys_driver.F` | Main physics driver (orchestration) |
| `mpas_atmphys_driver_convection.F` | Convection schemes |
| `mpas_atmphys_driver_microphysics.F` | Cloud microphysics |
| `mpas_atmphys_driver_gwdo.F` | Gravity wave drag over orography |
| `mpas_atmphys_driver_pbl.F` | Planetary boundary layer |
| `mpas_atmphys_driver_radiation_lw.F` | Long-wave radiation |
| `mpas_atmphys_driver_radiation_sw.F` | Short-wave radiation |
| `mpas_atmphys_driver_sfclayer.F` | Surface layer |
| `mpas_atmphys_driver_lsm.F` | Land surface model |
| `mpas_atmphys_driver_lsm_noahmp.F` | NoahMP land surface model |

#### Physics Scheme Options

**Convection:**
- Kain-Fritsch (cumulus parameterization)
- Tiedtke (mass flux scheme)
- NTIEDTKE (new Tiedtke)

**Microphysics:**
- Thompson (aerosol-aware)
- WSM6 (WRF Single-Moment 6-class)
- Kessler (simple warm rain)

**Radiation:**
- RRTMG (long-wave and short-wave)
- CAM radiation (optional)

**Land Surface:**
- Noah LSM
- NoahMP (multi-physics, ~280 files)

**PBL (Planetary Boundary Layer):**
- MYNN (Mellor-Yamada-Nakanishi-Niino)
- YSU (Yonsei University scheme)

**Gravity Waves:**
- UGWP (Unified Gravity Wave Physics)

#### Physics Submodules

| Directory | Source | Contents |
|-----------|--------|----------|
| `physics_noahmp/` | NoahMP | Full land surface model (~280 files) |
| `physics_wrf/` | WRF | Legacy WRF physics schemes |
| `physics_mmm/` | MMM | Git submodule for MMM physics |
| `physics_noaa/UGWP/` | NOAA | Git submodule for UGWP |

### Chemistry (`chemistry/`)

CheMPAS chemistry is integrated directly into the atmosphere core. The
chemistry layer manages:

- runtime tracer discovery from the MICM configuration
- MPAS <-> MICM state transfer
- operator-split lightning NOx source injection
- fallback `cos(SZA)` photolysis and TUV-x photolysis
- chemistry diagnostics such as `j_no2`

```
chemistry/
├── Makefile                        # Chemistry build integration
├── mpas_atm_chemistry.F            # Top-level chemistry manager and MPAS coupling
├── mpas_lightning_nox.F            # Operator-split lightning NOx source
├── mpas_solar_geometry.F           # Fallback solar zenith-angle calculation
├── mpas_tuvx.F                     # TUV-x photolysis wrapper
└── musica/
    └── mpas_musica.F               # MUSICA/MICM coupler and rate-parameter updates
```

See [MUSICA_INTEGRATION.md](../musica/MUSICA_INTEGRATION.md) for details.

#### Chemistry Responsibilities

| File | Purpose |
|------|---------|
| `mpas_atm_chemistry.F` | Initializes chemistry, resolves tracer indices, updates photolysis, runs source/solver coupling |
| `musica/mpas_musica.F` | Owns MICM state, species mapping, and photolysis-rate updates into MUSICA |
| `mpas_lightning_nox.F` | Applies the altitude- and updraft-dependent NO source before the chemistry solve |
| `mpas_solar_geometry.F` | Computes fallback solar geometry for Phase 1-style photolysis |
| `mpas_tuvx.F` | Computes profile-dependent `j_no2` with TUV-x, including cloud-opacity support |

#### Chemistry Call Path

At runtime, the chemistry-specific flow is:

```
chemistry_init
  |-> musica_init()
  |-> resolve_mpas_indices()
  |-> chemistry_lightning_nox_init()
  |-> tuvx_init()                  # when config_tuvx_config_file is set

chemistry_step
  |-> update photolysis (TUV-x or fallback cos(SZA))
  |-> lightning_nox_inject()
  |-> chemistry_from_MPAS()
  |-> musica_step()
  |-> chemistry_to_MPAS()
```

### Diagnostics (`diagnostics/`)

Post-processing and diagnostic calculations.

| File | Purpose |
|------|---------|
| `mpas_atm_diagnostics_manager.F` | Diagnostics orchestration |
| `mpas_isobaric_diagnostics.F` | Pressure-level interpolation |
| `mpas_pv_diagnostics.F` | Potential vorticity |
| `mpas_convective_diagnostics.F` | Convection diagnostics |
| `mpas_cloud_diagnostics.F` | Cloud properties |

Each diagnostic category has a `Registry_*.xml` file defining its variables.

---

## Core Init Atmosphere (`src/core_init_atmosphere/`)

Preprocessing tool for initializing atmosphere simulations.

### Key Files

| File | Size | Purpose |
|------|------|---------|
| `mpas_init_atm_cases.F` | 301KB | Idealized test case initialization |
| `mpas_init_atm_static.F` | 93KB | Static field initialization |
| `mpas_init_atm_llxy.F` | 79KB | Lat/lon to x/y coordinate conversion |
| `mpas_init_atm_hinterp.F` | 43KB | Horizontal interpolation |

### Capabilities

- Initialize from GFS/ERA5/CFSR analysis
- Set up idealized test cases (Jablonowski-Williamson, baroclinic wave, etc.)
- Generate terrain-following coordinates
- Interpolate to MPAS mesh
- Initialize surface fields

---

## Framework (`src/framework/`)

Shared infrastructure used by all cores.

### Core Services

| Module | Purpose |
|--------|---------|
| `mpas_dmpar` | MPI communication, halo exchanges |
| `mpas_io` | NetCDF/PnetCDF I/O |
| `mpas_stream_manager` | Stream handling and configuration |
| `mpas_pool_routines` | Variable pool management |
| `mpas_field_routines` | Field manipulation utilities |
| `mpas_timekeeping` | Time and calendar management |
| `mpas_log` | Logging infrastructure |
| `mpas_block_creator` | Domain decomposition |
| `mpas_hash_table` | Hash table implementation |
| `mpas_forcing` | External forcing interface |

### Data Structures

**Pools:** Named collections of variables, hierarchically organized.

**Fields:** Multi-dimensional arrays with metadata (dimensions, I/O info, halo status).

**Blocks:** Portions of the domain assigned to MPI ranks.

**Streams:** I/O channels for input, output, restart files.

---

## Operators (`src/operators/`)

Mathematical operations on mesh data.

| Operator | Purpose |
|----------|---------|
| `mpas_rbf_interpolation` | Radial basis function interpolation |
| `mpas_geometry_utils` | Geometric calculations |
| `mpas_tensor_operations` | Tensor computations on spherical mesh |
| `mpas_vector_operations` | Vector field operations |
| `mpas_tracer_advection` | Conservative tracer advection |
| `mpas_matrix_operations` | Sparse matrix utilities |
| `mpas_spline_interpolation` | Spline interpolation |

---

## Other Cores

### Ocean (`src/core_ocean/`)

Ocean circulation modeling with:
- Barotropic/baroclinic splitting
- Isopycnal/z-coordinate options
- Sea ice coupling interface
- BGC (biogeochemistry) tracers

### Sea Ice (`src/core_seaice/`)

Sea ice dynamics and thermodynamics:
- EVP (elastic-viscous-plastic) rheology
- Multi-category ice thickness distribution
- Column physics (vertical thermodynamics)

### Land Ice (`src/core_landice/`)

Ice sheet modeling:
- Shallow ice approximation
- First-order Stokes
- Albany/LI coupling (higher-order dynamics)

### Shallow Water (`src/core_sw/`)

Test cases for:
- Williamson test suite
- Numerical method verification
- Performance benchmarking

---

## Registry System

Each core has a `Registry.xml` file defining:

### Dimensions
```xml
<dims>
    <dim name="nCells"/>
    <dim name="nEdges"/>
    <dim name="nVertices"/>
    <dim name="nVertLevels" definition="80"/>
</dims>
```

### Variables
```xml
<var name="theta" type="real" dimensions="nVertLevels nCells Time"
     persistence="persistent" streams="input;output;restart"/>
```

### Namelists
```xml
<nml_record name="nhyd_model">
    <nml_option name="config_dt" type="real" default_value="450.0"/>
</nml_record>
```

### Streams
```xml
<stream name="output" type="output" filename_template="output.$Y-$M-$D_$h.$m.$s.nc">
    <var name="theta"/>
</stream>
```

---

## Related Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - High-level architecture
- [BUILD.md](../../BUILD.md) - Build system
- [MUSICA_INTEGRATION.md](../musica/MUSICA_INTEGRATION.md) - Chemistry integration
- [TEST_RUNS.md](../results/TEST_RUNS.md) - Runtime validation notes and case summaries
