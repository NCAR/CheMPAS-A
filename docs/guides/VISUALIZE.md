# MPAS Chemistry Visualization

This document describes tools for visualizing MPAS-MUSICA chemistry output.

## Python Environment

A conda environment `mpas` provides the required packages:

```bash
# Create environment (if not already done)
conda create -n mpas python=3.11 numpy matplotlib netcdf4 -y

# Use directly
~/miniconda3/envs/mpas/bin/python <script>

# Or activate
conda activate mpas
```

**Required packages:** numpy, matplotlib, netcdf4

## Scripts

All scripts are in the `scripts/` directory with symlinks in `~/Data/CheMPAS/supercell/`.

### plot_chemistry.py

Visualize chemistry tracer output (currently `qA`, `qB`, `qAB` for ABBA tests).

**Basic usage:**
```bash
cd ~/Data/CheMPAS/supercell
~/miniconda3/envs/mpas/bin/python plot_chemistry.py -o chemistry.png
```

**Options:**

| Option | Description |
|--------|-------------|
| `-i, --input` | Input file (default: output.nc) |
| `-o, --output` | Output figure filename |
| `-l, --level` | Vertical level for slices (default: 10) |
| `-t, --time` | Time index (default: -1 = last) |
| `--time-series` | Generate spatial maps at multiple times |
| `--diff` | Generate difference plots (t - t0) |
| `--diff-consecutive` | Generate consecutive diffs (t - t-1) |
| `--n-times` | Number of time steps for series (default: 6) |
| `--show` | Display plot interactively |

**Examples:**

```bash
# Quick summary plot
python plot_chemistry.py -o quick.png

# Spatial time evolution (6 panels showing pattern evolution)
python plot_chemistry.py -o advection.png --time-series

# Difference from initial (reveals advection + chemistry)
python plot_chemistry.py -o advection.png --diff

# Consecutive differences (instantaneous changes)
python plot_chemistry.py -o advection.png --diff-consecutive

# Specific level and time
python plot_chemistry.py -o level20.png --level 20 --time 5

# More time panels
python plot_chemistry.py -o detailed.png --time-series --n-times 9
```

**Output figures:**

| Suffix | Content |
|--------|---------|
| `.png` | Main 3x3 summary (horizontal slices, vertical cross-sections, time evolution) |
| `_timeseries.png` | Spatial maps at multiple times |
| `_diff.png` | Difference from initial conditions |
| `_diff_consecutive.png` | Consecutive time differences |

### init_tracer_sine.py

Initialize tracers with sine wave patterns for advection studies.

**Basic usage:**
```bash
python init_tracer_sine.py -i supercell_init.nc -t qAB --waves-x 2 --amplitude 0.4 --offset 0.6
```

**Options:**

| Option | Description |
|--------|-------------|
| `-i, --input` | Input init file (default: supercell_init.nc) |
| `-o, --output` | Output file (default: edit in place) |
| `-t, --tracer` | Tracer variable name (default: qAB) |
| `--amplitude` | Sine wave amplitude (default: 0.5) |
| `--offset` | Baseline value (default: 1.0) |
| `--waves-x` | Number of waves in x direction (default: 1) |
| `--waves-y` | Number of waves in y direction (default: 1) |

**Examples:**

```bash
# 2x2 wave pattern, values 0.2 to 1.0
python init_tracer_sine.py -t qAB --waves-x 2 --waves-y 2 --amplitude 0.4 --offset 0.6

# Single wave in x only
python init_tracer_sine.py -t qAB --waves-x 1 --waves-y 0 --amplitude 0.5 --offset 0.5

# Save to new file instead of editing in place
python init_tracer_sine.py -i supercell_init.nc -o supercell_init_sine.nc -t qAB
```

## Workflows

### Quick Look

After a run, generate a quick summary:

```bash
cd ~/Data/CheMPAS/supercell
~/miniconda3/envs/mpas/bin/python plot_chemistry.py -o quick.png
open quick.png
```

### Advection Study

1. **Set up initial conditions with gradients:**
   ```bash
   cp supercell_init.nc supercell_init_uniform.nc  # Backup
   python init_tracer_sine.py -t qAB --waves-x 2 --amplitude 0.4 --offset 0.6
   ```

2. **Configure longer run** (edit `namelist.atmosphere`):
   ```
   config_run_duration = '00:15:00'
   ```

3. **Adjust output interval** (edit `streams.atmosphere`):
   ```xml
   output_interval="00:00:30"
   ```

4. **Run the model:**
   ```bash
   mpiexec -n 8 ~/EarthSystem/CheMPAS-A/atmosphere_model
   ```

5. **Generate all visualizations:**
   ```bash
   python plot_chemistry.py -o advection.png --time-series --diff
   open advection.png advection_timeseries.png advection_diff.png
   ```

### Interpreting Results

**Time series plots:**
- Show spatial pattern evolution
- Chemistry decay: overall values decrease as AB → A + B
- Advection: pattern distortion/displacement

**Difference plots (t - t0):**
- Blue: tracer decreased (chemistry decay)
- Red: tracer increased (advection brought higher values)
- Symmetric pattern: chemistry-dominated
- Asymmetric pattern: advection effects visible

**Consecutive diffs (t - t-1):**
- Show instantaneous rate of change
- Useful for seeing where chemistry is most active
- Large values at concentration peaks (more reactant available)

## Technical Notes

### Unstructured Mesh Handling

The scripts use matplotlib's `Triangulation` to visualize the MPAS unstructured mesh:

```python
from matplotlib.tri import Triangulation
tri = Triangulation(xCell, yCell)  # Delaunay triangulation
ax.tricontourf(tri, values, ...)
```

**Limitations:**
- Uses Delaunay triangulation of cell centers (not actual MPAS Voronoi topology)
- May have minor artifacts at domain edges
- Future: consider uxarray for proper mesh handling (see ../project/TODO.md)

### Output Variables

Chemistry tracers in `output.nc` for the ABBA test mechanism:

| Variable | Description | Units |
|----------|-------------|-------|
| `qAB` | Molecular AB mixing ratio | kg/kg |
| `qA` | Atomic A mixing ratio | kg/kg |
| `qB` | Atomic B mixing ratio | kg/kg |

Wind fields for understanding advection:

| Variable | Description |
|----------|-------------|
| `uReconstructZonal` | Zonal (east-west) wind at cell centers |
| `uReconstructMeridional` | Meridional (north-south) wind at cell centers |
| `w` | Vertical velocity |
