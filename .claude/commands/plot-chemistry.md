# Plot Chemistry Output

Generate and open chemistry visualization plots from MPAS output.

## Arguments

- `$ARGUMENTS` - Optional: plot types and options
  - `all` - Generate all plot types (default)
  - `temporal` - Time series plots only
  - `spatial` - Spatial evolution plots only
  - `single-cell` - Single cell time series only
  - `vertical` - Vertical cross-section only
  - `full` - Use full time range instead of default 3 min

## Instructions

### 1. Set Up

```bash
TEST_DIR="$HOME/Data/MPAS/supercell"
SCRIPT="$HOME/EarthSystem/CheMPAS/scripts/plot_chemistry.py"
PYTHON="$HOME/miniconda3/envs/mpas/bin/python"
cd "$TEST_DIR"
```

Verify `output.nc` exists. If not, inform user to run the model first.

### 2. Parse Arguments

Default behavior generates all plots at both levels (5 km and 10 km).

Common options for all runs:
- `--wind --wind-skip 150 --wind-scale 200` for wind vectors on spatial plots
- Levels 10 (5 km) and 20 (10 km)

If `full` is in arguments, add `--time-end 60` to show full 30-minute run.

### 3. Generate Plots

**For `all` or default:**
```bash
$PYTHON $SCRIPT -o chemistry_L20.png --level 20 --temporal --multi-species --diff-consecutive --wind --wind-skip 150 --wind-scale 200
$PYTHON $SCRIPT -o chemistry_L10.png --level 10 --temporal --multi-species --diff-consecutive --wind --wind-skip 150 --wind-scale 200
$PYTHON $SCRIPT -o chemistry.png --single-cell --time-end 60
$PYTHON $SCRIPT -o chemistry.png --vertical
```

**For `temporal` only:**
```bash
$PYTHON $SCRIPT -o chemistry_L20.png --level 20 --temporal
$PYTHON $SCRIPT -o chemistry_L10.png --level 10 --temporal
```

**For `spatial` only:**
```bash
$PYTHON $SCRIPT -o chemistry_L20.png --level 20 --multi-species --diff-consecutive --wind --wind-skip 150 --wind-scale 200
$PYTHON $SCRIPT -o chemistry_L10.png --level 10 --multi-species --diff-consecutive --wind --wind-skip 150 --wind-scale 200
```

**For `single-cell` only:**
```bash
$PYTHON $SCRIPT -o chemistry.png --single-cell --time-end 60
```

**For `vertical` only:**
```bash
$PYTHON $SCRIPT -o chemistry.png --vertical
```

### 4. Clean Up

Remove the main summary plots (chemistry_L10.png, chemistry_L20.png, chemistry.png) as they contain vertical slices we don't need:

```bash
rm -f chemistry_L10.png chemistry_L10.pdf chemistry_L20.png chemistry_L20.pdf chemistry.png chemistry.pdf 2>/dev/null
```

### 5. Open Plots

```bash
open *.png
```

### 6. Report

List the generated plots with descriptions:

| Plot | Description |
|------|-------------|
| `*_temporal.png` | Domain-mean time series (first 3 min) |
| `*_multispecies.png` | qA + qAB spatial evolution with wind vectors |
| `*_diff_consecutive.png` | Frame-to-frame differences |
| `*_single_cell.png` | Single cell at 5 km and 10 km (full run) |
| `*_vertical.png` | Vertical cross-section through updraft |

## Plot Script Options Reference

| Option | Description |
|--------|-------------|
| `--level N` | Vertical level (10=5km, 20=10km) |
| `--temporal` | Generate time series plot |
| `--multi-species` | Generate qA + qAB spatial evolution |
| `--diff-consecutive` | Generate frame-to-frame diff plots |
| `--single-cell` | Generate single-cell time series |
| `--vertical` | Generate vertical cross-section |
| `--y-slice N` | Y coordinate for vertical slice (default: max updraft) |
| `--wind` | Add wind vectors to spatial plots |
| `--wind-skip N` | Wind vector spacing (default: 50) |
| `--wind-scale N` | Wind arrow scale (default: 300) |
| `--time-end N` | End time index (default: 6 = 3 min) |
| `--n-times N` | Number of time panels (default: 6) |
