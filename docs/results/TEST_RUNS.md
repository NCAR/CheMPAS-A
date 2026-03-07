# MPAS-MUSICA Test Runs

Record of test runs for chemistry-atmosphere coupling validation.

## Test Case: Supercell

Idealized supercell thunderstorm with ABBA chemistry (AB ↔ A + B reactions).

**Location:** `~/Data/MPAS/supercell`

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

---

### 2026-03-06: LNOx-O3 Case B — 15-minute storm with realistic source

**Duration:** 15 minutes
**Chemistry:** LNOx-O3 (NO+O3→NO2 Arrhenius, NO2+hv→NO+O3 photolysis)
**Source:** 0.5 ppbv/s when w - w_threshold = w_ref, w_threshold=5 m/s, w_ref=10 m/s, altitude 5–12 km
**Photolysis:** j_NO2 = 0.01 s⁻¹ (constant)
**Sink:** Disabled (tau = 0)
**Initial conditions:** NO=0, NO2=0, O3=50 ppbv

**Key observations:**
- Lightning source activates ~12 min (updraft reaches w > 5 m/s)
- NO peaks at 27.8 ppbv at t=15 min
- NO2 peaks at 5.9 ppbv (from Arrhenius titration)
- O3 depleted from 50 to 41.5 ppbv in updraft core (8.5 ppbv depletion)
- O3 background stable at 50 ppbv away from storm
- All tracers non-negative
- Photolysis active — NO2 recycled back to NO+O3

**Verification passed:**
- [x] Non-negativity for all species
- [x] O3 background preserved at 50 ppbv
- [x] O3 titration in updraft core
- [x] NO2 produced by Arrhenius reaction
- [x] Photolysis prevents NO2 accumulation
- [x] Physically realistic concentrations with 0.5 ppbv/s source

---

### 2026-03-06: LNOx-O3 Case A — 2-minute equilibrium diagnostic

**Duration:** 2 minutes
**Chemistry:** Same mechanism (`lnox_o3.yaml`, no LOSS reactions), source disabled
**Initial conditions:** NO=5, NO2=5, O3=50 ppbv (uniform)

**Key observations:**
- Chemistry drives system toward photostationary equilibrium
- Domain-integrated Ox (O3+NO2) conserved to machine precision (0.0000% drift)
- Domain-integrated NOx (NO+NO2) conserved to machine precision (0.0000% drift)
- Single-cell Ox varies due to advection (expected — advection redistributes, not destroys)
- All tracers non-negative

**MICM FIRST_ORDER_LOSS bug discovered:**
- With `lnox_o3.yaml` containing LOSS reactions (rate set to 0), Ox drifted 8.67%
  and NOx dropped 48% in 2 minutes — spurious mass loss
- MICM applies a nonzero loss even when rate parameters are explicitly set to 0
- Fix: default config (`lnox_o3.yaml`) omits LOSS reactions; separate config
  (`lnox_o3_sink.yaml`) available when finite sink timescale is intended

---

### 2026-03-06: Phase 1 — Solar geometry / SZA-dependent j_NO2

**Duration:** 30 minutes
**Chemistry:** LNOx-O3 with SZA-scaled photolysis: j_NO2 = j_max * max(0, cos_sza)
**Solar geometry:** Spencer (1971) algorithm in `mpas_solar_geometry.F`
**Location:** Kingfisher, OK (35.86°N, 97.93°W)
**Start time:** 0000-01-01 18:00 UTC (daytime at Kingfisher)
**j_NO2 max:** 0.01 s⁻¹

**Key observations:**
- cos_sza = 0.508 at 18:00 UTC (matches analytical prediction exactly)
- cos_sza evolves to 0.516 by 18:30 UTC (sun rising slightly)
- j_NO2 = 0.00508 → 0.00516 s⁻¹ (= j_max × cos_sza)
- Night test (midnight UTC): cos_sza = -0.9198, j_NO2 = 0 (correct)
- NO peak: 285.5 ppbv, NO2 peak: 32.6 ppbv, O3 background: 50.0 ppbv
- All tracers non-negative

**Verification passed:**
- [x] SZA matches analytical value for time/location
- [x] j_NO2 = 0 when sun below horizon
- [x] j_NO2 = j_max × cos_sza when sun above horizon
- [x] Non-negativity for all species
- [x] Physically plausible tracer concentrations

---

### 2026-03-06: Phase 2 — TUV-x clear-sky photolysis

**Duration:** 15 minutes
**Chemistry:** LNOx-O3 with TUV-x-computed j_NO2
**Photolysis config:** `config_tuvx_config_file = 'tuvx_no2.json'`
**Location / time:** Kingfisher, OK (35.86°N, 97.93°W), 18:00 UTC, SZA ≈ 59°
**Source:** 0.5 ppbv/s when w - w_threshold = w_ref, w_threshold=5 m/s, w_ref=10 m/s, altitude 5–12 km
**Sink:** Disabled (tau = 0)
**Initial conditions:** NO=0, NO2=0, O3=50 ppbv

**Key observations:**
- j_NO2 is no longer uniform with height
- j_NO2 increases from 7.2e-3 s^-1 at the surface to 1.2e-2 s^-1 at 20 km
- surface/top j_NO2 ratio = 0.61
- NO peak: 29.9 ppbv
- NO2 peak: 6.5 ppbv
- O3 minimum: 43.5 ppbv
- O3 background remains 50.000 ppbv away from the storm
- all tracers remain non-negative

**Phase comparison (15-minute Case B):**
- Phase 1 surface j_NO2: 5.1e-3 s^-1 (uniform with height)
- Phase 2 surface j_NO2: 7.2e-3 s^-1
- Phase 1 10 km j_NO2: 5.1e-3 s^-1
- Phase 2 10 km j_NO2: 1.0e-2 s^-1
- Phase 1 20 km j_NO2: 5.1e-3 s^-1
- Phase 2 20 km j_NO2: 1.2e-2 s^-1
- Phase 1 NO peak: 27.9 ppbv
- Phase 2 NO peak: 29.9 ppbv
- Phase 1 NO2 peak: 8.5 ppbv
- Phase 2 NO2 peak: 6.5 ppbv
- Phase 1 O3 minimum: 41.5 ppbv
- Phase 2 O3 minimum: 43.5 ppbv

**Verification passed:**
- [x] Non-negativity for all species
- [x] Nighttime shortcut sets j_NO2 = 0 when SZA >= 90 deg
- [x] j_NO2 shows plausible clear-sky vertical structure
- [x] Surface j_NO2 is in the expected clear-sky order of magnitude
- [x] O3 background preserved at 50 ppbv
- [x] 15-minute Case B remains stable with TUV-x photolysis
- [x] Empty `config_tuvx_config_file` falls back to the Phase 1 path (manual comparison)

**Deferred follow-up checks:**
- [ ] transition-smooth for extended runs spanning sunset
- [ ] decomposition compare across MPI decompositions
- [ ] Ox conservation with transport/source disabled

---

### 2026-03-06: Phase 3 — TUV-x cloud opacity

**Duration:** 15 minutes
**Chemistry:** LNOx-O3 with TUV-x + cloud radiator
**Photolysis config:** `config_tuvx_config_file = 'tuvx_no2.json'`
**Cloud radiator:** From-host, cloud water (qc) + rain (qr) from Kessler microphysics
**Cloud optical depth:** `tau = 3*LWC*dz / (2*r_eff*rho_water)`, r_eff=10 um (cloud), 500 um (rain)
**Cloud optical properties:** SSA=0.999999, g=0.85 (wavelength-independent)
**Location / time:** Kingfisher, OK (35.86°N, 97.93°W), 18:00 UTC, SZA ≈ 59°
**Source:** 0.5 ppbv/s when w - w_threshold = w_ref
**Initial conditions:** NO=0, NO2=0, O3=50 ppbv

**Key observations:**
- Cloud develops ~5 min into the simulation; max qc = 2.78e-3 kg/kg
- 346 cloudy columns (qc_max > 1e-6) vs 27,734 clear columns at t=15 min
- j_NO2 strongly attenuated below/inside optically thick cloud
- j_NO2 enhanced ~1.5x above cloud (delta-Eddington cloud albedo effect)
- NO2 accumulates inside cloud (photolysis recycling shut off)
- Clear-sky j_NO2 unchanged from Phase 2

**Phase 2 vs Phase 3 comparison (15-minute Case B):**

| Metric | Phase 2 (clear-sky) | Phase 3 (cloud) |
|--------|-------------------|----------------|
| Surface j_NO2 (clear) | 7.2e-3 s⁻¹ | 7.16e-3 s⁻¹ (unchanged) |
| Surface j_NO2 (cloudy min) | — | 6.27e-5 s⁻¹ (114x attenuation) |
| Surface j_NO2 (cloudy mean) | — | 6.68e-4 s⁻¹ |
| Above-cloud j_NO2 | ~1.0e-2 s⁻¹ | ~1.5e-2 s⁻¹ (1.5x enhanced) |
| NO2 max (cloudy) | 6.5 ppbv | 5.9 ppbv |
| NO2 max (clear) | — | 0.012 ppbv |
| O3 background | 50.0 ppbv | 50.0 ppbv |

**Verification passed:**
- [x] Cloud attenuation: j_NO2 = 6.3e-5 s⁻¹ inside thick cloud (< 1e-4 criterion)
- [x] Clear-sky unchanged: 7.16e-3 s⁻¹ (matches Phase 2)
- [x] Above-cloud enhancement: ~1.5x (physically correct cloud albedo effect)
- [x] Non-negativity for all species
- [x] NO2 higher in cloudy columns (less photolysis recycling)
- [x] O3 background preserved at 50 ppbv
- [x] 15-minute Case B stable with cloud-attenuated photolysis
- [x] Empty config_tuvx_config_file still falls back to Phase 1 path

**Deferred:**
- [ ] Clear-sky/cloudy column split optimization (per-column TUV-x for all cells currently)
- [ ] Ice-phase hydrometeors for non-Kessler microphysics
- [ ] Phase 2 vs Phase 3 chemistry response comparison plots

---

## Notes

1. **PnetCDF compatibility:** Use `io_type="netcdf"` in streams.atmosphere on macOS/LLVM builds
2. **Level mapping:** Level 10 ≈ 5 km, Level 20 ≈ 10 km
3. **Wind vectors:** Use `--wind-skip 150 --wind-scale 200` for readable arrows
4. **Vertical cross-section:** Auto-detects max updraft location; use `--y-slice N` to override
