# Photolysis and Tropospheric Chemistry Integration Plan

## Document Status

- `Historical Context:` Adapted from ancestor project plans — the TUV-x
  photolysis plan (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`) and the DAVINCI
  lightning-NOx/O3 mechanism (`DAVINCI-MPAS/PLAN.md` Phase 6, `SCIENCE.md`).
- `Current State:` Phase 0 complete. Phase 1 is partially implemented for the
  idealized supercell path: fallback solar geometry (`mpas_solar_geometry.F`)
  drives per-step scalar `j_NO2` updates in chemistry, and the chemistry-only
  build path is in place. The broader `coszr` / mesh-coordinate path is
  explicitly deferred because the needed mesh coordinates are not available in
  the current workflow. Current Phase 1 validation only requires daytime at the
  test coordinates from the namelist; DC3 calendar-time validation is deferred
  until there is a proper grid with grid coordinates. Remaining Phase 1 work:
  update stale metadata and add the planned analytical-SZA gate. Phase 2
  (TUV-x coupling) is not implemented yet.
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
7. **Phase 1 timing:** Use the Kingfisher test coordinates
   (35.86°N, 97.93°W) from the namelist for fallback SZA. For current Phase 1
   validation, any UTC start time that yields daytime at those coordinates is
   acceptable; the tracked namelist uses `0000-01-01_18:00:00`. Defer the
   DC3-specific timestamp (`2012-05-29_21:00:00`) to later validation once a
   proper grid with grid coordinates is available.
8. **TUV-x runs every chemistry timestep.** No caching or interval — keep it
   simple, optimize later if needed. Our runs are short (≤30 min).

## Strategic Direction

The ABBA mechanism validated the MPAS-MICM coupling infrastructure and the
Phase 2 runtime tracer allocation. The LNOx-O3 mechanism (Phase 0) replaced
ABBA with scientifically meaningful tropospheric photochemistry. Now we need
to replace the prescribed constant j_NO2 with time- and profile-dependent
photolysis rates.

**Approach:** Rather than following the ancestor plan's Chapman-first path, we
keep the working LNOx-O3 mechanism and add TUV-x to compute j_NO2 as a
function of solar zenith angle, altitude, and atmospheric profiles. This gives
us:

1. Diurnal cycle — photolysis shuts off at night, NO2 accumulates
2. Altitude dependence — j_NO2 increases with height (less atmospheric
   absorption)
3. Clear-sky illumination realism for DC3 validation — correct late-afternoon
   to evening geometry for the May 29 Kingfisher storm

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
geographic coordinates. Prefer the existing MPAS radiation `coszr` diagnostic
when it is available; otherwise compute the same solar geometry in chemistry.
Replace the constant j_NO2 with a SZA-dependent scaling:
`j_NO2 = j_max * max(0, cos(SZA))`. Validate day/night behavior.

**Rationale:** SZA computation is a prerequisite for TUV-x (Phase 2), which
requires SZA as input. Testing with a simple cosine scaling first validates
the solar-geometry plumbing before adding TUV-x radiative transfer complexity,
and keeps chemistry aligned with MPAS radiation geometry.

### Test Case Configuration

The idealized supercell uses the Kingfisher, Oklahoma test coordinates for
Phase 1 fallback-SZA validation:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Latitude | 35.86°N | Kingfisher, OK |
| Longitude | 97.93°W | |
| Start time | Any daytime UTC at the test coordinates | Tracked namelist currently uses `0000-01-01 18:00 UTC` |
| Duration | 30 min (Case B) | SZA changes modestly over this window |
| j_NO2 max | ~0.01 s⁻¹ | Daytime peak (surface, clear sky) |

For the currently tracked synthetic daytime case (`0000-01-01 18:00 UTC`) at
Kingfisher, `cos_sza ≈ 0.508`. That is sufficient for Phase 1 plumbing because
the fallback path only needs a reproducible daytime solar angle at the test
coordinates. The DC3-specific timing reference (`2012-05-29 21:00 UTC`,
`cos_sza ≈ 0.812`) remains useful later, but is deferred until grid-aware
validation is in scope.

Implementation status note (2026-03-06 review): the tracked idealized
namelist currently uses `0000-01-01_18:00:00` for a deterministic daytime
Phase 1 check. That is intentionally sufficient for fallback-SZA plumbing; the
DC3 calendar timestamp is deferred to later grid-aware validation.

### New Namelist Parameters

Add to `&musica` in `Registry.xml`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config_chemistry_latitude` | real | 0.0 | Reference latitude for SZA [degrees N] |
| `config_chemistry_longitude` | real | 0.0 | Reference longitude for SZA [degrees E] |

These are used for idealized cases where `latCell`/`lonCell` may not represent
real geographic coordinates. For real-data cases, per-cell lat/lon from the
mesh should be used instead.

### New Source Files / Fallback Helper

Prefer reusing MPAS radiation's `coszr` diagnostic when chemistry can access
it. If `coszr` is unavailable in the chemistry call path (for example,
idealized cases with radiation disabled), add:

**`src/core_atmosphere/chemistry/mpas_solar_geometry.F`**

New module computing SZA from calendar time and geographic position:

```fortran
module mpas_solar_geometry
    ! Fallback solar geometry helper for chemistry-only use.
    ! Algorithm: Spencer (1971) solar declination + hour angle.
    ! Must match the formula already used by MPAS radiation.
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
| `mpas_atm_chemistry.F` | Prefer `coszr` from MPAS diagnostics when available; otherwise call fallback solar geometry each step; pass cos_sza to musica |
| `mpas_musica.F` | Modify `assign_rate_parameters` → per-cell j_NO2 = j_max * max(0, cos_sza); update rate params each step (not just init) |
| `chemistry/Makefile` | Add `mpas_solar_geometry.o` only if the fallback helper is kept |
| `test_cases/supercell/namelist.atmosphere` | Add fallback lat/lon; current tracked Phase 1 namelist uses temporary `0000-01-01_18:00:00` start while DC3 timestamp revalidation remains pending |

### Key Architecture Change

Currently `assign_rate_parameters` is called once at init. In Phase 1, rate
parameters must be updated **every chemistry timestep** because j_NO2 varies
with SZA. This means:

- Prefer `coszr` from MPAS radiation diagnostics when it exists; only fall
  back to chemistry-side solar geometry when that diagnostic is unavailable
- Factor out rate parameter update from init into a per-step routine
- The per-step routine takes cos_sza and sets `PHOTO.no2_photolysis` =
  `config_lnox_j_no2 * max(0, cos_sza)` for each cell
- For uniform-SZA (idealized, single lat/lon), all cells get the same value
- Infrastructure supports per-cell SZA for future real-data cases

Implementation status note (2026-03-06 review): the committed code currently
implements only the fallback helper path. Reuse of MPAS `coszr` and real-data
SZA from mesh `latCell/lonCell` is deferred until the required mesh-coordinate
data is available in the chemistry workflow.

### Implementation Checklist

- [x] In `mpas_atm_core.F`, thread the current model time (`currTime`) into
  the chemistry call path for per-step SZA updates.
- [ ] Deferred: thread `diag_physics` and/or mesh `latCell/lonCell` into
  chemistry so Phase 1 can use existing MPAS solar geometry for real-data
  cases instead of relying only on fallback namelist coordinates. Blocked on
  mesh-coordinate availability.
- [x] Create `mpas_solar_geometry.F` with Spencer (1971) algorithm:
  `solar_cos_sza(DoY, hour_utc, lat_deg, lon_deg)` → cos(SZA).
- [x] In `chemistry_init`, cache `config_chemistry_latitude`,
  `config_chemistry_longitude`, and `config_lnox_j_no2` as module-level
  variables for per-step use.
- [x] In `mpas_musica.F`, add `musica_set_photolysis(j_max, cos_sza, ...)` that
  updates `PHOTO.no2_photolysis` rate parameter every timestep.
- [x] Cache the MICM rate-parameter index for `PHOTO.no2_photolysis` at init
  (`cache_photo_rp_index`) so per-step update is a direct array fill.
- [ ] Remove init-time `PHOTO.no2_photolysis` setup from `musica_init` /
  `assign_rate_parameters`; current implementation still carries the old
  constant-photolysis initialization path and then overwrites it per step.
- [x] Phase 1 uses scalar `j_NO2` update (uniform SZA across all cells).
- [x] Update supercell test namelist with fallback lat/lon and a temporary
  `0000-01-01_18:00:00` start for deterministic daytime verification.
- [ ] Deferred: restore the tracked Phase 1 reference case to the
  DC3-aligned `2012-05-29_21:00:00` timestamp only after there is a proper
  grid with usable grid coordinates for location-aware validation.
- [ ] Update `Registry.xml` metadata for `config_lnox_j_no2` to describe it as
  daytime `j_max`, not a constant prescribed photolysis rate.
- [ ] Add the analytical SZA check to the adapted `check_tuvx_phase.py`.
  The imported gate scripts are now wired to CheMPAS defaults
  (`qNO,qNO2,qO3`, `j_no2`) and include `night-jzero`, but this specific
  analytical-SZA gate is still missing.

### High-Level Fortran Design

Current committed path (implemented and accepted for now):

1. Thread `currTime` into `chemistry_step`.
2. Compute `cos_sza` from cached fallback lat/lon via `solar_cos_sza(...)`.
3. Call `musica_set_photolysis(...)` before `musica_step(...)`.

The broader Phase 1 target below remains deferred future work once
`diag_physics` / mesh-coordinate support is available.

The least disruptive Phase 1 path is:

1. Extend the chemistry driver call so `diag_physics` is available.
2. Resolve one `cos_sza` value for the idealized case.
3. Convert that into a scalar `j_NO2`.
4. Push the scalar photolysis rate into MICM before `musica_step`.

```fortran
! mpas_atm_core.F
call mpas_pool_get_subpool(block_chem % structs, 'diag_physics', diag_physics_chem)
call chemistry_step(dt, current_time, mesh_chem, state_chem, diag_chem,      &
                    diag_physics_chem, block_chem % dimensions, time_lev_chem)
```

```fortran
! mpas_atm_chemistry.F
call chemistry_get_cos_sza(diag_physics, current_time, cos_sza, &
                           error_code, error_message)
if (error_code /= 0) return

call musica_set_photolysis_scalar(config_lnox_j_no2, cos_sza, &
                                  error_code, error_message)
if (error_code /= 0) return
```

```fortran
! mpas_musica.F
subroutine musica_set_photolysis_scalar(j_no2_max, cos_sza, error_code, error_message)
    real(kind=RKIND), intent(in) :: j_no2_max, cos_sza
    real(kind=RKIND) :: j_no2_val

    j_no2_val = j_no2_max * max(0.0_RKIND, cos_sza)
    call fill_rate_parameter(state, photo_no2_rp_index, j_no2_val, &
                             error_code, error_message)
end subroutine
```

This keeps the Phase 1 code change narrow: `mpas_atm_core.F` gains one more
subpool argument, `mpas_atm_chemistry.F` gets a solar-geometry helper, and
`mpas_musica.F` gets a dedicated photolysis-update entry point.

### Verification (Phase 1 Gate)

| Check | Criterion | Script |
|-------|-----------|--------|
| Non-negativity | qNO, qNO2, qO3 >= 0 | `check_tuvx_phase.py nonnegative` |
| Night j-zero | j_NO2 = 0 when SZA >= 90° | `check_tuvx_phase.py night-jzero` |
| SZA correctness | cos(SZA) matches analytical value for given time/location | Manual / unit test |
| Day/night response | NO2 accumulates at night (no photolysis), NO2→NO+O3 during day | Visual inspection |

### Exit Criteria

- [x] Build passes with solar-geometry plumbing (fallback Spencer helper).
- [x] Fallback-SZA computation matches the synthetic idealized check used by
  the tracked namelist (`0000-01-01 18:00 UTC`, cos_sza = 0.508 actual vs
  0.508 predicted).
- [x] 30-min Case B run produces physically plausible results with SZA-scaled
  j_NO2 (j = 0.00508 at start, evolves to 0.00516 over 30 min).
- [x] Night test: j_NO2 = 0 when cos_sza < 0 (verified at midnight UTC,
  cos_sza = -0.9198).
- [x] All tracers non-negative throughout run.
- [ ] Deferred: re-run the Phase 1 checks on the DC3-aligned May 29, 2012
  timestamp after grid-aware validation becomes meaningful.
- [ ] Extended run (5+ hours) shows j_NO2 → 0 at sunset (deferred).
- [x] Add j_NO2 as diagnostic field in output.nc (3D field in diag pool,
  added to stream_list, verified in 30-min run: uniform 5.08–5.16e-3 s⁻¹).
- [x] Add j_NO2 photolysis plot to `plot_lnox_o3.py` (`--photolysis` flag).

---

## Phase 2: TUV-x Coupled Photolysis

**Goal:** Replace the simple `j_max * cos(SZA)` scaling with TUV-x radiative
transfer, giving clear-sky, in-domain j_NO2 as a function of altitude, SZA,
and atmospheric profiles (temperature, pressure, O3 column).

**Rationale:** The cos(SZA) scaling from Phase 1 gives correct day/night
behavior but not the altitude dependence or atmospheric absorption effects
captured by radiative transfer. TUV-x computes the full actinic flux →
cross-section × quantum-yield integration. This phase is explicitly a
clear-sky coupling milestone, not the final storm-radiation treatment.

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
call air_updater%set_midpoint_values(air_number_density)  ! [molecule/cm³]
call temp_updater%set_midpoint_values(temperature) ! [K]
call tuvx_core%run(sza, earth_sun_distance, j_values)
! j_values contains j_NO2(level) [s⁻¹]
```

where `air_number_density = (rho_air / M_AIR) * N_A * 1.0e-6`, converting
MPAS `rho_air` [kg/m³] to molecules/cm³.

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
- Clouds/aerosols: deferred in Phase 2 (clear-sky only)

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

2. **Profile sources.** Height edges from `zgrid`, air number density derived
   from `rho_air`, and temperature from the same fields already extracted in
   `chemistry_from_MPAS`. Convert `rho_air` [kg/m³] to molecules/cm³ before
   passing it to TUV-x. O3 profile comes from the model's own `qO3` field
   (self-consistent).

3. **Rate parameter update.** The j_NO2 array (nVertLevels × nCells) replaces
   the scalar `config_lnox_j_no2`. Rate parameters updated every step via the
   same `PHOTO.no2_photolysis` index, but now with per-cell/level values.

4. **Fallback.** If `config_tuvx_config_file` is empty, fall back to Phase 1
   behavior (cos_sza scaling with `config_lnox_j_no2`). This preserves the
   simpler mode for quick testing.

5. **Validation scope.** Phase 2 remains clear-sky and in-domain. It does not
   yet include cloud optics, aerosol optics, or an atmosphere above the MPAS
   model top, so absolute magnitude checks are sanity checks rather than hard
   pass/fail literature matches.

### Implementation Checklist

- [ ] Add `config_tuvx_config_file` to `Registry.xml` and read it in
  `chemistry_init`.
- [ ] Add `mpas_tuvx.F` plus the required makefile objects/modules.
- [ ] Initialize TUV-x once during `chemistry_init` when
  `config_tuvx_config_file` is non-empty; otherwise leave the TUV-x pointer
  unassociated and use the Phase 1 fallback branch.
- [ ] Refactor `chemistry_from_MPAS` so the environment-building loop
  (`rho_air`, `temperature`, `pressure`) is reusable by both MICM coupling and
  TUV-x, rather than recomputing those fields twice.
- [ ] Add a chemistry-side workspace for `j_no2(:,:)` with dimensions
  `(nVertLevels, nCells)` for the current block.
- [ ] Derive `air_number_density` from `rho_air` using the documented
  `kg/m³ -> molecules/cm³` conversion before each TUV-x column solve.
- [ ] Extract the model `qO3` profile for each column and convert it to the
  units required by the chosen TUV-x O3 updater.
- [ ] In `mpas_musica.F`, add a per-field photolysis update routine that writes
  the full `j_no2(level, cell)` field into the `PHOTO.no2_photolysis`
  rate-parameter slice.
- [ ] Keep the fallback behavior in the same driver path so Phase 1 and Phase 2
  both use the same `musica_set_photolysis_*` routines.
- [x] Extend the phase-gate scripts with `fallback-compare`,
  `transition-smooth`, and decomposition checks. The imported ancestor scripts
  are now adapted to the CheMPAS Phase 0/1/2 matrix.

### High-Level Fortran Design

Phase 2 is easiest to keep coherent if the chemistry driver owns the column
work arrays and chooses between two photolysis providers:

1. Phase 1 scalar fallback: `j = j_max * max(0, cos_sza)`
2. Phase 2 TUV-x provider: `j = TUV-x(z, rho, T, O3, SZA)`

That keeps MICM blind to where the photolysis field came from.

```fortran
! mpas_atm_chemistry.F
call chemistry_build_environment(mesh, state, diag, time_lev, &
                                 rho_air, temperature, pressure, zgrid, qo3, &
                                 error_code, error_message)
if (error_code /= 0) return

call chemistry_get_cos_sza(diag_physics, current_time, cos_sza, &
                           error_code, error_message)
if (error_code /= 0) return

if (tuvx_is_enabled()) then
    do iCell = 1, nCells
        call tuvx_compute_photolysis(cos_sza, zgrid(:, iCell), rho_air(:, iCell), &
                                     temperature(:, iCell), qo3(:, iCell),        &
                                     j_no2(:, iCell), error_code, error_message)
        if (error_code /= 0) return
    end do
    call musica_set_photolysis_field(j_no2, nCells, nVertLevels, &
                                     error_code, error_message)
else
    call musica_set_photolysis_scalar(config_lnox_j_no2, cos_sza, &
                                      error_code, error_message)
end if

call MICM_from_chemistry(scalars, rho_air, temperature, pressure, &
                         nCells, nVertLevels, error_code, error_message)
call musica_step(time_step, error_code, error_message)
```

```fortran
! mpas_tuvx.F
subroutine tuvx_compute_photolysis(cos_sza, z_edges_m, rho_air, temperature, qo3, &
                                   j_no2_column, error_code, error_message)
    real(kind=RKIND), intent(in)  :: cos_sza
    real(kind=RKIND), intent(in)  :: z_edges_m(:)
    real(kind=RKIND), intent(in)  :: rho_air(:), temperature(:), qo3(:)
    real(kind=RKIND), intent(out) :: j_no2_column(:)

    call update_height_profile(z_edges_m * 1.0e-3_RKIND)
    call update_air_profile(kgm3_to_moleccm3(rho_air))
    call update_temperature_profile(temperature)
    call update_o3_profile(qo3)
    call run_tuvx(acos(max(0.0_RKIND, cos_sza)), j_no2_column)
end subroutine
```

```fortran
! mpas_musica.F
subroutine musica_set_photolysis_field(j_no2, nCells, nVertLevels, error_code, error_message)
    real(kind=RKIND), intent(in) :: j_no2(nVertLevels, nCells)
    integer, intent(in)          :: nCells, nVertLevels

    call fill_rate_parameter_field(state, photo_no2_rp_index, j_no2, &
                                   nCells, nVertLevels, error_code, error_message)
end subroutine
```

The key refactor is to turn "set rate parameters" into a stable API with two
front doors:

- `musica_set_photolysis_scalar(...)` for Phase 1 and fallback mode
- `musica_set_photolysis_field(...)` for TUV-x mode

That prevents the chemistry driver from knowing MICM's stride details, and it
keeps the Phase 2 fallback path cheap to maintain.

### Verification (Phase 2 Gate)

| Check | Criterion | Script |
|-------|-----------|--------|
| Non-negativity | qNO, qNO2, qO3 >= 0 | `check_tuvx_phase.py nonnegative` |
| Night j-zero | j_NO2 = 0 when SZA >= 90° | `check_tuvx_phase.py night-jzero` |
| Transition smooth | Dawn/dusk j_NO2 varies smoothly | `check_tuvx_phase.py transition-smooth` |
| Altitude profile | j_NO2 increases with height | Log inspection / diagnostic output |
| Magnitude sanity | Surface j_NO2 is in the expected clear-sky order of magnitude; document bias from 20 km top and omitted cloud/aerosol optics | Literature sanity check / offline reference |
| Ox conservation | Domain-integrated Ox conserved (source/sink off) | `verify_ox_conservation.py` |
| Decomp compare | Identical results across MPI decompositions | `check_tuvx_phase.py decomp-compare` |
| Fallback compare | Empty `config_tuvx_config_file` reproduces Phase 1 behavior within roundoff | `check_tuvx_phase.py fallback-compare` |

### Exit Criteria

- Build passes with TUV-x module linked.
- j_NO2 profile shows expected altitude dependence.
- Surface j_NO2 remains in the expected clear-sky order of magnitude, with
  documented caveats from the 20 km top and omitted cloud/aerosol optics.
- 30-min Case B remains stable with TUV-x-provided j_NO2.
- All Phase 1 gate checks still pass.
- `fallback-compare` passes when `config_tuvx_config_file` is empty.
- `transition-smooth` passes for extended (5+ hour) runs spanning sunset.

---

## Phase Gate Scripts

Adapted from ancestor project (`MPAS-Model-ACOM-dev/scripts/`):

| Script | Purpose |
|--------|---------|
| `scripts/check_tuvx_phase.py` | Imported ancestor phase-gate checker, now adapted to CheMPAS defaults (`qNO,qNO2,qO3`, `j_no2`) and extended with `fallback-compare` |
| `scripts/run_tuvx_phase_gate.sh` | Phase wrapper adapted to the current CheMPAS Phase 0/1/2 matrix (`nonnegative`, `verify_ox_conservation.py`, `night-jzero`, `transition-smooth`, `decomp-compare`, `fallback-compare`) |

These scripts operate on MPAS NetCDF output and exit non-zero on failure.
Adapted for LNOx-O3 species (qNO, qNO2, qO3) rather than ancestor Chapman
species (qO, qO2, qO3).

Status note (2026-03-06 update): both scripts have now been copied from the
ancestor repo into `scripts/` and adapted to the current CheMPAS species and
phase matrix. Remaining gaps are narrower: the analytical-SZA Phase 1 check is
still not implemented, and runtime use still depends on the Python `netCDF4`
module plus whichever `j_*` diagnostics are written to output files.

### Phase Gate Matrix

| Phase | Checks |
|-------|--------|
| Phase 0 | `nonnegative`, `verify_ox_conservation.py` |
| Phase 1 | Phase 0 + `night-jzero` |
| Phase 2 | Phase 1 + `transition-smooth`, `decomp-compare`, `fallback-compare` |

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
3. **Domain-top limitation** — 20 km top omits part of the overhead ozone
   column, so Phase 2 magnitude checks are sanity checks rather than strict
   literature targets.
4. **Source/sink representation split** — Lightning-NOx source is operator-split
   pre-MICM; NOx sink remains mechanism-defined within MICM.
5. **Cloud/aerosol optics deferred** — Phase 2 is clear-sky only until cloud
   and aerosol optical inputs are coupled into TUV-x.

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
- `scripts/check_tuvx_phase.py` — copied into this repo from the ancestor and
  adapted for CheMPAS species defaults and `fallback-compare`
- `scripts/run_tuvx_phase_gate.sh` — copied into this repo from the ancestor
  and adapted to the current CheMPAS phase matrix

## Dependencies

- MUSICA-Fortran v0.13.0 with MICM + TUV-x support (already linked)
- MICM LNOx-O3 mechanism config (`micm_configs/lnox_o3.yaml`)
- TUV-x config for NO2 photolysis (`micm_configs/tuvx_no2.json`, Phase 2)
- Python `netCDF4` for verification scripts
