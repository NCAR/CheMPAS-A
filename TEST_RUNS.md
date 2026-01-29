# MPAS-MUSICA Test Runs

Record of test runs for chemistry-atmosphere coupling validation.

## Test Case: Supercell

Idealized supercell thunderstorm with ABBA chemistry (AB ↔ A + B reactions).

**Location:** `~/MPAS/supercell`

### Run Configuration

| Parameter | Value |
|-----------|-------|
| Domain | 84 km × 84 km × 20 km |
| Horizontal resolution | ~1 km (28,080 cells) |
| Vertical levels | 40 |
| Timestep | 3 s |
| Output interval | 30 s |
| MPI ranks | 8 |

### Chemistry Configuration (abba.yaml)

| Reaction | Scaling Factor | Notes |
|----------|----------------|-------|
| AB → A + B | 2.0e-3 | Forward (dissociation) |
| A + B → AB | 1.0e-3 | Reverse (recombination) |

**Note:** Scaling factors reduced 10x from original (2.0e-2, 1.0e-2) to slow equilibration and better visualize advection effects.

### Initial Conditions

- **qAB:** Sine wave pattern (0.2 to 1.0 kg/kg), 2 waves in X direction
- **qA, qB:** Zero (produced by chemistry)
- **Wind:** ~15 m/s background zonal flow

```bash
# Apply sine wave to qAB
python init_tracer_sine.py -t qAB --waves-x 2 --amplitude 0.4 --offset 0.6
```

---

## Run History

### 2025-01-15: 30-minute run with slow chemistry

**Duration:** 30 minutes
**Chemistry:** Slow (scaling factors 2.0e-3, 1.0e-3)
**Output:** 61 time steps, 11 GB

**Key observations:**
- Supercell develops near domain center by t=15 min
- Max updraft ~30 m/s at t=30 min (x=38.7, y=39.1 km)
- Chemistry equilibration takes ~10 min (vs ~2 min with fast chemistry)
- Advection effects visible in single-cell time series
- Vertical transport of tracers by convective updraft

**Convection development:**

| Time | Max updraft (5 km) |
|------|-------------------|
| 0 min | 0 m/s |
| 8 min | 0.8 m/s |
| 15 min | 19 m/s |
| 22 min | 25 m/s |
| 30 min | 29 m/s |

**Generated plots:**
- `chemistry_L10_temporal.png` - Time series at 5 km
- `chemistry_L20_temporal.png` - Time series at 10 km
- `chemistry_L10_multispecies.png` - qA + qAB spatial at 5 km
- `chemistry_L20_multispecies.png` - qA + qAB spatial at 10 km
- `chemistry_L10_diff_consecutive.png` - Frame diffs at 5 km
- `chemistry_L20_diff_consecutive.png` - Frame diffs at 10 km
- `chemistry_single_cell.png` - Single cell at both heights
- `chemistry_vertical.png` - Vertical cross-section through updraft

---

### 2025-01-15: 15-minute run with slow chemistry

**Duration:** 15 minutes
**Chemistry:** Slow (scaling factors 2.0e-3, 1.0e-3)
**Output:** 31 time steps, 5.8 GB

**Key observations:**
- Initial chemistry equilibration visible
- Convection beginning to develop
- Sine wave pattern still partially visible

---

### 2025-01-15: 15-minute run with fast chemistry

**Duration:** 15 minutes
**Chemistry:** Fast (scaling factors 2.0e-2, 1.0e-2)
**Output:** 31 time steps, 5.8 GB

**Key observations:**
- Chemistry equilibrates in ~2 minutes
- Advection effects masked by rapid chemistry
- Single-cell time series shows flat profile after equilibration

---

## Visualization Commands

```bash
# Generate all plots
/plot-chemistry

# Specific plot types
/plot-chemistry temporal
/plot-chemistry spatial
/plot-chemistry single-cell
/plot-chemistry vertical

# Full time range
/plot-chemistry full
```

## Notes

1. **PnetCDF compatibility:** Use `io_type="netcdf"` in streams.atmosphere on macOS/LLVM builds
2. **Level mapping:** Level 10 ≈ 5 km, Level 20 ≈ 10 km
3. **Wind vectors:** Use `--wind-skip 150 --wind-scale 200` for readable arrows
4. **Vertical cross-section:** Auto-detects max updraft location; use `--y-slice N` to override
