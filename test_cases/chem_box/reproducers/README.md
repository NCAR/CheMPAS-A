# chem_box Chapman + NOx Reproducers

This directory contains small standalone diagnostics for the tracked
`chem_box` case. They are meant to isolate MICM/MUSICA chemistry behaviour
without rerunning the full MPAS atmosphere executable.

## `chapman_nox_level2_reproducer.F90`

Replays the smallest currently known failing state: one local
`chem_box` grid cell at level 2, immediately before the first 3 s
Chapman + NOx MICM solve.

The default run uses the raw concentration and photolysis-rate slots
captured from CheMPAS-A. It is expected to reproduce the bad solve, not
pass a physical acceptance test:

- solver: `RosenbrockStandardOrder`
- timestep: `3.0 s`
- relative tolerance: `1.0e-15`
- expected result: O3 jumps from about `1.04e-6 mol m-3` to about
  `3.13e-1 mol m-3` in one solve.

Run from the repository root:

```bash
test_cases/chem_box/reproducers/run_chapman_nox_level2_reproducer.sh
```

Optional arguments are passed through to the executable:

```bash
test_cases/chem_box/reproducers/run_chapman_nox_level2_reproducer.sh 0.001 rosenbrock 1e-15
test_cases/chem_box/reproducers/run_chapman_nox_level2_reproducer.sh 3.0 backward_euler 1e-15
test_cases/chem_box/reproducers/run_chapman_nox_level2_reproducer.sh 3.0 rosenbrock 1e-15 name_species_raw_rates
test_cases/chem_box/reproducers/run_chapman_nox_level2_reproducer.sh 3.0 rosenbrock 1e-15 name_mapped
```

The runner uses `scripts/check_build_env.sh --export`, `mpifort`, and
`pkg-config --libs musica-fortran`. It writes the temporary executable
under `${TMPDIR:-/tmp}/chempas_chapman_nox_reproducer`.
