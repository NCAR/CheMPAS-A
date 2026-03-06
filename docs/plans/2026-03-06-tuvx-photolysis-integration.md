# Photolysis and Tropospheric Chemistry Integration Plan

## Document Status

- `Historical Context:` Adapted from ancestor project plans — the TUV-x
  photolysis plan (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`) and the DAVINCI
  lightning-NOx/O3 mechanism (`DAVINCI-MPAS/PLAN.md` Phase 6, `SCIENCE.md`).
- `Current State:` Phase 0 complete (2026-03-06). LNOx-O3 mechanism runs
  end-to-end. Domain-integrated Ox/NOx conservation verified to machine
  precision. Ready for Phase 1 (solar geometry) and Phase 2 (TUV-x coupling).
- `Use This As:` Primary reference for post-ABBA chemistry development.

## Locked Decisions

1. **NOx sink path:** Carry over NOx loss via MICM-configured sink terms (not
   hard-coded Fortran sink tendencies). Use separate MICM config
   (`lnox_o3_sink.yaml`) when sink is enabled — MICM FIRST_ORDER_LOSS applies
   nonzero loss even with rate=0 (tracked for upstream report).
2. **Lightning source path:** Operator-split pre-MICM injection via standalone
   `mpas_lightning_nox.F` module. MICM doesn't support spatially varying source
   terms, so injection happens before solve.
3. **Verification:** Include a dedicated controlled-NOx pulse/equilibration
   case for Leighton-ratio validation.
4. **MICM units:** Concentrations in mol/m³. Arrhenius `A` must be in m³/mol/s
   (convert from cm³/molecule/s via `A × Nₐ × 10⁻⁶`). Photolysis rate
   parameters in s⁻¹, set externally.
5. **Initial conditions:** Set via MPAS init file (kg/kg), not MICM config.
   MICM `__initial concentration` is required by the parser but overwritten
   by state transfer; set to 0.
6. **Keep LNOx-O3 mechanism through TUV-x integration.** The ancestor plan
   targets Chapman oxygen photolysis (stratospheric), but our supercell domain
   is tropospheric (0–20 km). j_NO2 is the only photolysis rate we need now.
   Chapman can be added as a later phase.
7. **DC3 reference case:** Use the May 29, 2012 Kingfisher, Oklahoma supercell
   (35.86°N, 97.93°W) as the idealized test case reference. Set
   `config_start_time = '2012-05-29_21:00:00'` for realistic SZA.
8. **TUV-x runs every chemistry timestep.** No caching or interval — keep it
   simple, optimize later if needed. Our runs are short (≤30 min).

## Strategic Direction

The ABBA mechanism validated the MPAS-MICM coupling infrastructure and the
Phase 2 runtime tracer allocation. The LNOx-O3 mechanism (Phase 0) replaced
ABBA with scientifically meaningful tropospheric photochemistry. Now we need
to replace the prescribed constant j_NO2 with physically realistic,
spatially varying photolysis rates.

**Approach:** Rather than following the ancestor plan's Chapman-first path, we
keep the working LNOx-O3 mechanism and add TUV-x to compute j_NO2 as a
function of solar zenith angle, altitude, and atmospheric profiles. This gives
us:

1. Diurnal cycle — photolysis shuts off at night, NO2 accumulates
2. Altitude dependence — j_NO2 increases with height (less atmospheric
   absorption)
3. Physical realism for DC3 validation — correct illumination geometry for
   the May 29 Kingfisher storm (late afternoon → evening transition)

## Phase 0: LNOx-O3 with Fixed Rates — COMPLETE

**Status:** All items complete and merged to main. See implementation progress
below.

### Chemistry

The tropospheric photostationary system (from DAVINCI Phase 6):

```
NO + O3 → NO2 + O2     (k: Arrhenius, temperature-dependent)
NO2 + hv → NO + O3     (j_NO2: prescribed constant, ~0.01 s⁻¹ daytime)
```

The net O3 reaction is `NO2 + hv → NO + O + (O + O2 + M → O3)`, written as a
single step because atomic O reaches steady state in microseconds in the
troposphere.

### ODE System

```
d[NO]/dt  =  j*[NO2] - k*[NO]*[O3] + S_ltg - [NO]/tau
d[NO2]/dt = -j*[NO2] + k*[NO]*[O3] - [NO2]/tau
d[O3]/dt  =  j*[NO2] - k*[NO]*[O3]
```

where `S_ltg` is the lightning NO source term and `tau` is a configurable NOx
sink timescale represented in the MICM mechanism.

**Conservation:** Ox = [O3] + [NO2] is conserved when both lightning source and
NOx sink are disabled. Verified to machine precision (0.0000% domain-integrated
drift) using `scripts/verify_ox_conservation.py`.

### Source/Sink Coupling Strategy

1. **Operator-split lightning source** — `mpas_lightning_nox.F` injects NO into
   MPAS scalars before MICM runs each timestep. Source scales linearly:
   `S = rate * max(0, w - w_thr) / w_ref` in cells where altitude is in range.

2. **Chemistry in MICM** — Arrhenius titration + photolysis defined in MICM
   config. MICM handles the ODE solve.

3. **Sink in MICM config** — Separate config `lnox_o3_sink.yaml` with
   `FIRST_ORDER_LOSS` reactions. Default `lnox_o3.yaml` omits LOSS reactions
   (MICM bug: nonzero loss even with rate=0).

4. **Namelist control** — Seven parameters in `&musica`.

### Implementation Progress

- [x] `mpas_lightning_nox.F` — lightning source module (init + inject)
- [x] Registry.xml — namelist parameters for lightning source
- [x] `mpas_atm_chemistry.F` — hook lightning init and inject into chemistry pipeline
- [x] `chemistry/Makefile` — build integration with dependency ordering
- [x] Build passes with MUSICA=true
- [x] `micm_configs/lnox_o3.yaml` — MICM config (Arrhenius + photolysis, no LOSS)
- [x] `micm_configs/lnox_o3_sink.yaml` — config with FIRST_ORDER_LOSS for NO/NO2
- [x] Runtime tracer discovery: `qNO`, `qNO2`, `qO3` created automatically
- [x] `scripts/init_lnox_o3.py` — initialize tracers in supercell_init.nc
- [x] `scripts/plot_lnox_o3.py` + `scripts/style.py` — visualization suite
- [x] `scripts/verify_ox_conservation.py` — domain-integrated Ox/NOx verification
- [x] Arrhenius A parameter corrected from cm³/molecule/s to m³/mol/s
- [x] Discovered MICM FIRST_ORDER_LOSS bug (nonzero loss with rate=0)
- [x] Case B (storm): 30-min run, 0.5 ppbv/s source, O3 titration verified
- [x] Case A (equilibrium): Ox/NOx conserved to machine precision
- [x] Unit conversions verified (ppbv ↔ kg/kg ↔ mol/m³)
- [x] `test_cases/supercell/` — reference namelists tracked in repo

### Verification Results

| Check | Result |
|-------|--------|
| Non-negativity | PASS — qNO, qNO2, qO3 >= 0 everywhere |
| O3 background | PASS — 50 ppbv away from storm |
| O3 titration | PASS — depleted to near-zero in updraft core |
| Ox conservation | PASS — 0.0000% domain-integrated drift (Case A) |
| NOx conservation | PASS — 0.0000% domain-integrated drift (Case A) |
| Unit consistency | PASS — ppbv conversions verified |

---

## Phase 1: Solar Geometry and Day/Night Photolysis

**Goal:** Compute per-cell solar zenith angle (SZA) from MPAS model time and
namelist-specified coordinates. Replace the constant j_NO2 with a
SZA-dependent scaling: `j_NO2 = j_max * max(0, cos(SZA))`. Validate day/night
behavior.

**Rationale:** SZA computation is a prerequisite for TUV-x (Phase 2), which
requires SZA as input. Testing with a simple cosine scaling first validates
the solar geometry plumbing before adding TUV-x radiative transfer complexity.

### Test Case Configuration

The idealized supercell uses the DC3 May 29, 2012 Kingfisher, Oklahoma storm
as its reference case:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Latitude | 35.86°N | Kingfisher, OK |
| Longitude | 97.93°W | |
| Start time | 2012-05-29 21:00 UTC | ~4 PM CDT, convection initiation |
| Duration | 30 min (Case B) | SZA changes ~3° over this window |
| j_NO2 max | ~0.01 s⁻¹ | Daytime peak (surface, clear sky) |

At 21:00 UTC on May 29 at Kingfisher, the SZA ≈ 45° (late afternoon). Over a
30-minute run, SZA increases by ~3°. Over a 3-hour run, SZA would reach ~90°
(sunset around 01:30 UTC / 8:30 PM CDT).

### New Namelist Parameters

Add to `&musica` in `Registry.xml`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config_chemistry_latitude` | real | 0.0 | Reference latitude for SZA [degrees N] |
| `config_chemistry_longitude` | real | 0.0 | Reference longitude for SZA [degrees E] |

These are used for idealized cases where `latCell`/`lonCell` may not represent
real geographic coordinates. For real-data cases, per-cell lat/lon from the
mesh should be used instead.

### New Source Files

**`src/core_atmosphere/chemistry/mpas_solar_geometry.F`**

New module computing SZA from calendar time and geographic position:

```fortran
module mpas_solar_geometry
    ! Compute solar zenith angle from model time + geographic coordinates.
    ! Algorithm: Spencer (1971) solar declination + hour angle.
    ! Consistent with MPAS physics coszr computation.
contains
    subroutine solar_geometry_init(latitude, longitude)
    subroutine solar_zenith_angle(year, month, day, hour, minute, second, &
                                   latitude, longitude, sza, cos_sza)
end module
```

Inputs: model time (from MPAS clock), latitude, longitude.
Outputs: SZA [radians], cos(SZA).

The Spencer (1971) algorithm computes:
1. Day of year → fractional year angle
2. Solar declination from fractional year
3. Hour angle from UTC + longitude
4. `cos(SZA) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(ha)`

### Modified Files

| File | Change |
|------|--------|
| `Registry.xml` | Add `config_chemistry_latitude`, `config_chemistry_longitude` |
| `mpas_atm_chemistry.F` | Read lat/lon from namelist, call solar geometry each step, pass cos_sza to musica |
| `mpas_musica.F` | Modify `assign_rate_parameters` → per-cell j_NO2 = j_max * max(0, cos_sza); update rate params each step (not just init) |
| `chemistry/Makefile` | Add `mpas_solar_geometry.o` to build |
| `test_cases/supercell/namelist.atmosphere` | Add lat/lon, change start time to `2012-05-29_21:00:00` |

### Key Architecture Change

Currently `assign_rate_parameters` is called once at init. In Phase 1, rate
parameters must be updated **every chemistry timestep** because j_NO2 varies
with SZA. This means:

- Factor out rate parameter update from init into a per-step routine
- The per-step routine takes cos_sza and sets `PHOTO.no2_photolysis` =
  `config_lnox_j_no2 * max(0, cos_sza)` for each cell
- For uniform-SZA (idealized, single lat/lon), all cells get the same value
- Infrastructure supports per-cell SZA for future real-data cases

### Verification (Phase 1 Gate)

| Check | Criterion | Script |
|-------|-----------|--------|
| Non-negativity | qNO, qNO2, qO3 >= 0 | `check_tuvx_phase.py nonnegative` |
| Night j-zero | j_NO2 = 0 when SZA >= 90° | `check_tuvx_phase.py night-jzero` |
| SZA correctness | cos(SZA) matches analytical value for given time/location | Manual / unit test |
| Day/night response | NO2 accumulates at night (no photolysis), NO2→NO+O3 during day | Visual inspection |

### Exit Criteria

- Build passes with new solar geometry module.
- SZA computation matches expected value for Kingfisher, OK at 21:00 UTC May 29
  (SZA ≈ 45°, cos_sza ≈ 0.707).
- 30-min Case B run produces physically plausible results with SZA-scaled j_NO2.
- Extended run (3+ hours) shows j_NO2 → 0 at sunset (~01:30 UTC).
- `night-jzero` gate check passes.

---

## Phase 2: TUV-x Coupled Photolysis

**Goal:** Replace the simple `j_max * cos(SZA)` scaling with TUV-x radiative
transfer, giving physically accurate j_NO2 as a function of altitude, SZA,
and atmospheric profiles (temperature, pressure, O3 column).

**Rationale:** The cos(SZA) scaling from Phase 1 gives correct day/night
behavior but not the altitude dependence or atmospheric absorption effects
that make photolysis physically realistic. TUV-x computes the full actinic
flux → cross-section × quantum-yield integration.

### TUV-x Availability

TUV-x is already compiled into MUSICA-Fortran v0.13.0 (132 `.mod` files in
`~/software/include/musica/fortran/`). No additional build dependencies.

Key modules:
- `musica_tuvx.mod` — high-level MUSICA wrapper
- `tuvx_core.mod` — core TUV-x engine
- `tuvx_grid_from_host.mod` — runtime grid updates from host model
- `tuvx_cross_section_no2_tint.mod` — NO2 cross-section (temperature-interpolated)
- `tuvx_spherical_geometry.mod` — SZA + slant-path calculations

### TUV-x API Pattern

```fortran
! Initialization (once)
tuvx_core => core_t("tuvx_config.json")
height_updater => tuvx_core%get_updater(height_grid)
air_updater    => tuvx_core%get_updater(air_profile)
temp_updater   => tuvx_core%get_updater(temp_profile)

! Each timestep, per column
call height_updater%set_edges(zgrid_km)       ! MPAS zgrid [m] → [km]
call air_updater%set_midpoint_values(air_density)  ! [molecule/cm³]
call temp_updater%set_midpoint_values(temperature) ! [K]
call tuvx_core%run(sza, earth_sun_distance, j_values)
! j_values contains j_NO2(level) [s⁻¹]
```

### New Source Files

**`src/core_atmosphere/chemistry/mpas_tuvx.F`**

New module wrapping TUV-x for CheMPAS:

```fortran
module mpas_tuvx
contains
    subroutine tuvx_init(config_file, nVertLevels)
        ! Create TUV-x core from JSON config
        ! Get profile updaters for height, air density, temperature
        ! Cache j_NO2 reaction index in TUV-x output

    subroutine tuvx_compute_photolysis(sza, zgrid, rho_air, temperature, &
                                        nVertLevels, j_no2_column)
        ! Update TUV-x profiles from MPAS column data
        ! Run TUV-x radiative transfer
        ! Extract j_NO2(level) from output

    subroutine tuvx_finalize()
end module
```

**TUV-x config file** (`micm_configs/tuvx_no2.json`)

Minimal TUV-x configuration for NO2 photolysis only:
- Height grid: from-host (MPAS zgrid)
- Atmospheric profiles: from-host (air density, temperature, O3)
- Cross-section: NO2 (temperature-interpolated)
- Quantum yield: NO2 → NO + O (= 1.0 below 420 nm)
- Wavelength grid: standard or fast_tuv grid
- Radiative transfer: delta-Eddington (2-stream)
- Surface albedo: 0.1 (grassland, configurable)

### Modified Files

| File | Change |
|------|--------|
| `Registry.xml` | Add `config_tuvx_config_file` |
| `mpas_atm_chemistry.F` | Call `tuvx_init` during chemistry init, call `tuvx_compute_photolysis` per column each step |
| `mpas_musica.F` | Accept j_NO2 array from TUV-x, inject into rate parameters per cell/level |
| `chemistry/Makefile` | Add `mpas_tuvx.o`, link TUV-x modules |

### Key Design Decisions

1. **Per-column computation.** TUV-x solves radiative transfer per column. Loop
   over nCells, calling TUV-x for each column with that cell's atmospheric
   profile. SZA is uniform (idealized case) but profiles vary by cell.

2. **Profile sources.** Height edges from `zgrid`, air density and temperature
   from the same fields already extracted in `chemistry_from_MPAS`. O3 profile
   from the model's own `qO3` field (self-consistent).

3. **Rate parameter update.** The j_NO2 array (nVertLevels × nCells) replaces
   the scalar `config_lnox_j_no2`. Rate parameters updated every step via the
   same `PHOTO.no2_photolysis` index, but now with per-cell/level values.

4. **Fallback.** If `config_tuvx_config_file` is empty, fall back to Phase 1
   behavior (cos_sza scaling with `config_lnox_j_no2`). This preserves the
   simpler mode for quick testing.

### Verification (Phase 2 Gate)

| Check | Criterion | Script |
|-------|-----------|--------|
| Non-negativity | qNO, qNO2, qO3 >= 0 | `check_tuvx_phase.py nonnegative` |
| Night j-zero | j_NO2 = 0 when SZA >= 90° | `check_tuvx_phase.py night-jzero` |
| Transition smooth | Dawn/dusk j_NO2 varies smoothly | `check_tuvx_phase.py transition-smooth` |
| Altitude profile | j_NO2 increases with height | Log inspection / diagnostic output |
| Magnitude | j_NO2 ≈ 0.005–0.01 s⁻¹ at surface, higher aloft | Literature comparison |
| Ox conservation | Domain-integrated Ox conserved (source/sink off) | `verify_ox_conservation.py` |
| Decomp compare | Identical results across MPI decompositions | `check_tuvx_phase.py decomp-compare` |

### Exit Criteria

- Build passes with TUV-x module linked.
- j_NO2 profile shows expected altitude dependence.
- j_NO2 magnitude consistent with literature (~0.008 s⁻¹ surface, clear sky).
- 30-min Case B produces physically plausible storm chemistry with TUV-x j_NO2.
- All Phase 1 gate checks still pass.
- `transition-smooth` passes for extended (3+ hour) runs spanning sunset.

---

## Phase Gate Scripts

Adapted from ancestor project (`MPAS-Model-ACOM-dev/scripts/`):

| Script | Purpose |
|--------|---------|
| `scripts/check_tuvx_phase.py` | Suite of physics-based checks (nonnegative, night-jzero, transition-smooth, Ox-budget, decomp-compare) |
| `scripts/run_tuvx_phase_gate.sh` | Phase orchestrator — runs correct combination of checks per phase |

These scripts operate on MPAS NetCDF output and exit non-zero on failure.
Adapted for LNOx-O3 species (qNO, qNO2, qO3) rather than ancestor Chapman
species (qO, qO2, qO3).

### Phase Gate Matrix

| Phase | Checks |
|-------|--------|
| Phase 0 | `nonnegative`, `verify_ox_conservation.py` |
| Phase 1 | Phase 0 + `night-jzero` |
| Phase 2 | Phase 1 + `transition-smooth`, `decomp-compare` |

---

## Later Phases

- Phase 3: Solver robustness under real-world forcing (6–24h stress runs)
- Phase 4: Real-world robustness and reproducibility
- Phase 5: Full Chapman (O2/O3/O photolysis) — extends mechanism, adds j_O2/j_O3
- Phase 6: Performance optimization
- Phase 7 (Optional): Extended NOx chemistry (PAN, HNO3, organic nitrates)

---

## Key Constraints

1. **Vertical grid from MPAS** — TUV-x height edges from `zgrid(:, iCell)`.
2. **Profiles from MPAS state** — No static atmosphere data files.
3. **Domain-top limitation** — 20 km top only samples troposphere/UTLS.
4. **Source/sink representation split** — Lightning-NOx source is operator-split
   pre-MICM; NOx sink remains mechanism-defined within MICM.
5. **Keep `state_ref`** — Useful for diagnosing advection effects on chemistry.

## Reference Material

### DC3 Campaign

The May 29, 2012 Kingfisher, Oklahoma supercell is the benchmark DC3 case:

- DiGangi et al. (2016) — "An overview of the 29 May 2012 Kingfisher supercell
  during DC3" (JGR Atmospheres)
- Pickering et al. (2024) — "Lightning NOx in the 29–30 May 2012 DC3 Severe
  Storm and Its Downwind Chemical Consequences" (JGR Atmospheres)
- Cummings et al. (2024) — "Evaluation of Lightning Flash Rate
  Parameterizations in a Cloud-Resolved WRF-Chem Simulation of the 29–30 May
  2012 Oklahoma Severe Supercell System Observed During DC3" (JGR Atmospheres)

### DAVINCI Sister Project

The DAVINCI project (`~/EarthSystem/DAVINCI-MPAS/`) contains:

- `SCIENCE.md` — Lightning NOx physics, Leighton framework, DC3 findings
- `PLAN.md` Phase 6 — LNOx-O3 mechanism details, ODE system, Jacobian,
  verification criteria
- `DC3.md` — Deep Convective Clouds and Chemistry campaign reference
- `TUV.md` — TUV-x algorithm and data file reference

### Ancestor TUV-x Plan

`MPAS-Model-ACOM-dev/PLAN_TUVx.md` contains:
- Full 9-phase TUV-x integration plan with physical verification gates
- Phase gate runbook (namelist, runtime settings, pass/fail scripts)
- Fortran implementation sketches
- Photolysis-to-MICM rate parameter mapping
- `scripts/check_tuvx_phase.py` — phase gate check script (to be adapted)
- `scripts/run_tuvx_phase_gate.sh` — phase gate orchestrator (to be adapted)

## Dependencies

- MUSICA-Fortran v0.13.0 with MICM + TUV-x support (already linked)
- MICM LNOx-O3 mechanism config (`micm_configs/lnox_o3.yaml`)
- TUV-x config for NO2 photolysis (`micm_configs/tuvx_no2.json`, Phase 2)
- Python `netCDF4` for verification scripts
