CheMPAS
=======

CheMPAS (Chemistry for MPAS) is a coupled atmospheric-chemistry model built on
the [MPAS](https://mpas-dev.github.io/) framework. It integrates
[MUSICA/MICM](https://github.com/NCAR/musica) atmospheric chemistry into
MPAS-Atmosphere, enabling chemical transport modeling on MPAS's unstructured
Voronoi meshes.

CheMPAS is derived from
[NCAR/MPAS-Model-ACOM-dev](https://github.com/NCAR/MPAS-Model-ACOM-dev)
(a fork of [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model)).
It is an independent project and does not sync with the upstream repositories.

## Agent-Driven Development

CheMPAS uses an agent-driven development model. AI agents handle code
development, review, and CI, with human oversight reserved for scientific
correctness, architectural decisions, and physical validation. See
[PURPOSE.md](PURPOSE.md) for the full rationale and
[AGENTS.md](AGENTS.md) for operational details.

## Building

CheMPAS builds with LLVM compilers (flang/clang) on macOS:

```bash
scripts/check_build_env.sh
eval "$(scripts/check_build_env.sh --export)"

make -j8 llvm \
  CORE=atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double \
  MUSICA=true
```

See [BUILD.md](BUILD.md) for the MUSICA/pkg-config preflight notes and
[RUN.md](RUN.md) for test case execution.

## Code Layout

```
CheMPAS/
├── src/
│   ├── driver              -- Main driver (standalone mode)
│   ├── external            -- External dependencies
│   ├── framework           -- MPAS framework (data types, comms, I/O)
│   ├── operators           -- Mesh operators
│   ├── tools/
│   │   ├── registry        -- Registry.xml parser
│   │   └── input_gen       -- Stream and namelist generators
│   ├── core_atmosphere/
│   │   ├── dynamics        -- Dynamical core
│   │   ├── physics         -- Physics parameterizations
│   │   └── chemistry/
│   │       └── musica      -- MUSICA/MICM coupling
│   └── core_init_atmosphere
├── scripts/                -- Visualization and analysis tools
├── testing_and_setup/      -- Test case configuration
└── default_inputs/         -- Default stream and namelist files
```

## Documentation

| Document | Description |
|----------|-------------|
| [PURPOSE.md](PURPOSE.md) | Project motivation and development philosophy |
| [AGENTS.md](AGENTS.md) | Agent roles, workflow, and review gates |
| [docs/README.md](docs/README.md) | Documentation index by topic |
| [ARCHITECTURE.md](docs/architecture/ARCHITECTURE.md) | System architecture |
| [BUILD.md](BUILD.md) | Build instructions |
| [RUN.md](RUN.md) | Running test cases |
| [MUSICA_INTEGRATION.md](docs/musica/MUSICA_INTEGRATION.md) | MUSICA/MICM coupling details |
| [MUSICA_API.md](docs/musica/MUSICA_API.md) | MUSICA Fortran API reference |
| [TEST_RUNS.md](docs/results/TEST_RUNS.md) | Recorded run outcomes and validation notes |
| [BENCHMARKS.md](docs/results/BENCHMARKS.md) | Agent and model benchmark comparison |
| [TUVX_INTEGRATION.md](docs/guides/TUVX_INTEGRATION.md) | TUV-x integration summary and development test case |
| [VISUALIZE.md](docs/guides/VISUALIZE.md) | Chemistry visualization tools |
| [User Guide](docs/users-guide/00-foreword.md) | Imported user guide chapters |

## License

See [LICENSE](LICENSE).
