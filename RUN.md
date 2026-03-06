# Running MPAS Test Cases

This document describes how to run MPAS atmosphere test cases.

## Supercell Test Case

The supercell thunderstorm is an idealized convection test case located at `~/Data/MPAS/supercell`. This directory is shared with the upstream MPAS-Model-ACOM-dev project.

### Important: I/O Configuration

The `streams.atmosphere` file must use `io_type="netcdf"` for input and output streams to avoid PnetCDF compatibility issues:

```xml
<immutable_stream name="input"
                  type="input"
                  io_type="netcdf"
                  filename_template="supercell_init.nc"
                  input_interval="initial_only"/>

<stream name="output"
        type="output"
        io_type="netcdf"
        ...>
```

### Prerequisites

1. Build the atmosphere model (see [BUILD.md](BUILD.md)):
   ```bash
   make -j8 llvm CORE=atmosphere PIO=$HOME/software NETCDF=/opt/homebrew PNETCDF=$HOME/software PRECISION=double
   ```

2. Verify the executable exists:
   ```bash
   ls -la atmosphere_model
   ```

### Running the Test

**Important:** You must remove or move any existing `output.nc` before running. MPAS defaults to `clobber_mode = never_modify`, so if `output.nc` already exists the model will silently skip all output writes — the run completes but produces no new data.

```bash
cd ~/Data/MPAS/supercell

# Archive previous run output (REQUIRED — model won't overwrite existing output.nc)
timestamp=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${timestamp}.nc
[ -f log.atmosphere.0000.out ] && mv log.atmosphere.0000.out log.atmosphere.0000.${timestamp}.out

# Run with 8 MPI ranks (recommended for 10-core machine)
mpiexec -n 8 ~/EarthSystem/CheMPAS/atmosphere_model
```

### MPI Rank Selection

Choose based on available cores and partition files:

| Ranks | Partition File | Notes |
|-------|----------------|-------|
| 8 | `supercell.graph.info.part.8` | Recommended for 10-core Mac |
| 12 | `supercell.graph.info.part.12` | Slight oversubscription |
| 16 | `supercell.graph.info.part.16` | Oversubscribed on 10 cores |

The partition file prefix is set in `namelist.atmosphere`:
```
config_block_decomp_file_prefix = 'supercell.graph.info.part.'
```

### Configuration

Key settings in `namelist.atmosphere`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `config_dt` | 3.0 | Timestep (seconds) |
| `config_run_duration` | '00:2:00' | Run length (2 minutes) |
| `config_micm_file` | 'abba.yaml' | MICM chemistry config |
| `config_microp_scheme` | 'mp_kessler' | Microphysics scheme |

### Output Files

| File | Description |
|------|-------------|
| `output.nc` | Model output (written every 1 minute) |
| `log.atmosphere.0000.out` | Run log with diagnostics |

### Verifying the Run

Check the log for successful completion:
```bash
tail -20 log.atmosphere.0000.out
```

Look for final timestep and timing statistics.

### Quick Run Script

Create a helper script `run.sh`:
```bash
#!/bin/bash
cd ~/Data/MPAS/supercell

# Archive previous output
ts=$(date +%Y%m%d_%H%M%S)
for f in output.nc log.atmosphere.*.out; do
    [ -f "$f" ] && mv "$f" "${f%.nc}.${ts}.nc" 2>/dev/null || mv "$f" "${f%.out}.${ts}.out" 2>/dev/null
done

# Run
mpiexec -n 8 ~/EarthSystem/CheMPAS/atmosphere_model 2>&1 | tee run.out
```

## Advection Studies

To visualize advection effects on chemistry tracers, use custom initial conditions with spatial gradients.

### Initialize with Sine Wave Pattern

```bash
cd ~/Data/MPAS/supercell

# Backup uniform initial conditions
cp supercell_init.nc supercell_init_uniform.nc

# Apply 2-wave horizontal sine pattern to qAB
~/miniconda3/envs/mpas/bin/python init_tracer_sine.py \
  -i supercell_init.nc -t qAB --waves-x 2 --amplitude 0.4 --offset 0.6
```

This creates qAB varying from 0.2 to 1.0 across the domain.

### Longer Runs for Advection

For visible advection displacement, run 10-15 minutes:

1. Edit `namelist.atmosphere`:
   ```
   config_run_duration = '00:15:00'
   ```

2. Adjust output interval in `streams.atmosphere` (e.g., 30 seconds):
   ```xml
   output_interval="00:00:30"
   ```

3. Run as normal. At ~15 m/s wind, expect ~10-15 km displacement over 15 minutes.

### Visualization

```bash
cd ~/Data/MPAS/supercell

# Main summary plot
~/miniconda3/envs/mpas/bin/python plot_chemistry.py -o chemistry.png

# Spatial time evolution (pattern at multiple times)
~/miniconda3/envs/mpas/bin/python plot_chemistry.py -o chemistry.png --time-series

# Difference from initial (shows advection + chemistry effects)
~/miniconda3/envs/mpas/bin/python plot_chemistry.py -o chemistry.png --diff

# Open plots
open chemistry.png chemistry_timeseries.png chemistry_diff.png
```

## LNOx-O3 Chemistry Test Cases

The LNOx-O3 mechanism simulates lightning-NOx-driven ozone photochemistry:
- `NO + O3 → NO2` (Arrhenius, temperature-dependent)
- `NO2 + hv → NO + O3` (photolysis, prescribed rate)
- First-order NOx loss (configurable timescale)

### Setup

1. Copy the MICM config to the run directory:
   ```bash
   cp ~/EarthSystem/CheMPAS/micm_configs/lnox_o3.yaml ~/Data/MPAS/supercell/
   ```

2. Initialize tracers (O3=50 ppbv, NO=NO2=0):
   ```bash
   cd ~/Data/MPAS/supercell
   ~/miniconda3/envs/mpas/bin/python ~/EarthSystem/CheMPAS/scripts/init_lnox_o3.py -i supercell_init.nc
   ```

3. Configure `namelist.atmosphere` `&musica` section:
   ```
   &musica
       config_micm_file = 'lnox_o3.yaml'
       config_lnox_source_rate = 0.5      ! NO source [ppbv/s] when w - w_threshold = w_ref (0 = disabled)
       config_lnox_w_threshold = 5.0      ! min updraft for injection [m/s]
       config_lnox_w_ref = 10.0           ! excess updraft used in source scaling [m/s]
       config_lnox_z_min = 5000.0         ! min altitude for injection [m]
       config_lnox_z_max = 12000.0        ! max altitude for injection [m]
       config_lnox_j_no2 = 0.01           ! NO2 photolysis rate [s⁻¹] (0 = disabled)
       config_lnox_nox_tau = 0.0          ! NOx sink timescale [s] (0 = no sink)
   /
   ```

The source scaling follows `S = rate * max(0, w - w_threshold) / w_ref`, so
`config_lnox_source_rate` is reached when `w - w_threshold = config_lnox_w_ref`.

### Case B: Storm Chemistry (15-min supercell)

The primary test case with lightning NOx source active:

```bash
cd ~/Data/MPAS/supercell
ts=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${ts}.nc
[ -f log.atmosphere.0000.out ] && mv log.atmosphere.0000.out log.atmosphere.0000.${ts}.out
mpiexec -n 8 ~/EarthSystem/CheMPAS/atmosphere_model
```

**Expected behavior:** NO injected in updraft core (w > 5 m/s, 5–12 km altitude),
O3 titrated where NO is present, NO2 produced via Arrhenius reaction, O3 recovery
downwind via NO2 photolysis.

### Case A: Equilibrium Diagnostic

Set `config_lnox_source_rate = 0.0` and initialize with uniform NO=5 ppbv,
NO2=5 ppbv, O3=50 ppbv to test Leighton equilibrium approach. Use a short
run duration (2 minutes). Note: Ox conservation is approximate in coupled runs
due to advection; a rigorous test requires transport disabled.

### Verification

| Check | What to look for |
|-------|-----------------|
| Non-negativity | `qNO >= 0`, `qNO2 >= 0`, `qO3 >= 0` everywhere |
| O3 background | O3 ≈ 50 ppbv away from storm |
| O3 titration | O3 depressed in updraft core where NO is injected |
| NO2 production | NO2 present where NO+O3 reaction occurs |
| Photolysis | NO2 converted back to NO+O3 (prevents NO2 accumulation) |

## Other Test Cases

Additional test cases can be added to `~/Data/MPAS/` following the same pattern.
