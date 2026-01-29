# MPAS-Model Architecture

This document describes the high-level architecture of the MPAS (Model for Prediction Across Scales) framework and its components.

## Overview

MPAS is a modular, unstructured mesh framework for Earth system modeling. It supports multiple specialized "cores" for different physical domains:

- **core_atmosphere** - Atmospheric modeling (primary focus of this branch)
- **core_init_atmosphere** - Initialization preprocessing
- **core_ocean** - Ocean modeling
- **core_seaice** - Sea ice modeling
- **core_landice** - Land ice/glacier modeling
- **core_sw** - Shallow water test cases
- **core_test** - Framework testing suite

## Architecture Diagram

```
                    MPAS Main Driver
                     (src/driver/)
                          |
            +-------------+-------------+
            |                           |
      Atmosphere Core             Ocean/Other Core
            |
            +---------------------------+
            |                           |
        Dynamics                    Physics
         (SRK3)                  (Convection,
                                  Radiation,
                                  Microphysics,
                                  LSM)
                                      |
                                  Chemistry
                                  (MUSICA/MICM)

   ================================================
        MPAS Framework (Shared)
        I/O, DMpar, Logging, Halo, etc.
   ================================================
```

## Directory Structure

```
CheMPAS/
├── CMakeLists.txt          # Main CMake build configuration
├── Makefile                # Legacy Makefile build system
├── cmake/                  # CMake modules and functions
├── docs/                   # Sphinx-based documentation (RST format)
├── src/                    # Main source code
│   ├── driver/             # Main execution entry points
│   ├── framework/          # Shared infrastructure (~65 files)
│   ├── operators/          # Mathematical operations
│   ├── tools/              # Code generation tools
│   ├── external/           # External libraries (ezxml, SMIOL, ESMF)
│   ├── core_atmosphere/    # Atmosphere modeling (~254 files)
│   ├── core_init_atmosphere/
│   ├── core_ocean/
│   ├── core_seaice/
│   ├── core_landice/
│   ├── core_sw/
│   └── core_test/
└── testing_and_setup/      # Test configurations
```

## Core Components

### 1. Framework (`src/framework/`)

The framework provides shared infrastructure used by all cores:

| File | Size | Purpose |
|------|------|---------|
| `mpas_dmpar.F` | 415KB | Distributed memory parallel communication |
| `mpas_io_streams.F` | 180KB | I/O streaming infrastructure |
| `mpas_io.F` | 268KB | Core I/O routines |
| `mpas_pool_routines.F` | 235KB | Data pool management |
| `mpas_stream_manager.F` | 280KB | Stream handling |
| `mpas_field_routines.F` | 118KB | Field manipulation |
| `mpas_block_creator.F` | 94KB | Block decomposition |

Additional services: logging, timekeeping, halo exchanges, forcing, bootstrapping, hash tables.

### 2. Operators (`src/operators/`)

Specialized mathematical kernels:

- `mpas_rbf_interpolation.F` - Radial basis function interpolation
- `mpas_geometry_utils.F` - Geometric calculations
- `mpas_tensor_operations.F` - Tensor computations
- `mpas_vector_operations.F` - Vector routines
- Tracer advection schemes (monotonic, standard)
- Matrix and spline operations

### 3. Driver (`src/driver/`)

Main execution entry points:

- `mpas.F` - Main program
- `mpas_subdriver.F` - Common subdriver for all cores

### 4. External Libraries (`src/external/`)

- `ezxml/` - XML parsing library
- `SMIOL/` - Scalable Modeling I/O Library
- `esmf_time_f90/` - ESMF time management (internal copy)

## Registry System

The `Registry.xml` files define model metadata:

- **Dimensions**: Mesh sizes (nCells, nEdges, nVertices, nVertLevels)
- **Variables**: State variables, diagnostics, tracers
- **Namelists**: Configuration parameters
- **Streams**: I/O definitions (input, output, restart files)

The registry is processed by tools in `src/tools/registry/` to generate Fortran code.

## Data Flow

1. **Initialization**: Driver reads configuration, initializes framework
2. **Domain Decomposition**: Mesh partitioned across MPI ranks
3. **Time Integration**: Core-specific dynamics/physics stepping
4. **Halo Exchange**: Framework handles inter-processor communication
5. **I/O**: Stream manager handles checkpointing and output

## Dependencies

### Required
- MPI (Fortran bindings)
- PnetCDF (Fortran bindings)

### Optional
- PIO (Parallel I/O library)
- NetCDF
- GPTL (profiling)
- ESMF
- MUSICA/MICM (chemistry, this branch)

## Code Statistics

| Component | Files | Description |
|-----------|-------|-------------|
| Framework | 65 | Core infrastructure |
| Operators | 12 | Mathematical kernels |
| Atmosphere | 254 | Dynamics + physics |
| Init Atmosphere | 33 | Preprocessing |
| Chemistry | 2 | MUSICA/MICM interface |
| Ocean | 50+ | Ocean modeling |

Total: ~600+ Fortran/C source files

## Related Documentation

- [BUILD.md](BUILD.md) - Build system documentation
- [COMPONENTS.md](COMPONENTS.md) - Detailed component documentation
- [MUSICA_INTEGRATION.md](MUSICA_INTEGRATION.md) - Chemistry integration details

## External Resources

- MPAS-Atmosphere User Guide: http://mpas-dev.github.io/atmosphere/
- MPAS GitHub: https://github.com/MPAS-Dev
