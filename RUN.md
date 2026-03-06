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

## Other Test Cases

Additional test cases can be added to `~/Data/MPAS/` following the same pattern.
