# Photolysis and Tropospheric Chemistry Integration Plan

## Document Status

- `Historical Context:` Adapted from ancestor project plans — the TUV-x
  photolysis plan (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`) and the DAVINCI
  lightning-NOx/O3 mechanism (`DAVINCI-MPAS/PLAN.md` Phase 6, `SCIENCE.md`).
- `Current State:` Phases 0–2 complete. Phase 2 (TUV-x coupled photolysis)
  is implemented and verified on the `develop` branch. TUV-x computes
  per-column, altitude-dependent j_NO2 from host atmospheric profiles via
  delta-Eddington radiative transfer. Surface j_NO2 = 7.2e-3 s⁻¹ at SZA≈59°
  (matches literature clear-sky). Phase 1 cos(SZA) fallback remains available
  when `config_tuvx_config_file` is empty. Phase 3 (cloud/aerosol opacity,
  earth-sun distance, extended validation) is the next target.
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

With Phase 1 accepted and merged on the fallback-only path, the active work
now shifts to Phase 2 TUV-x coupling.

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

## Phase 1: Solar Geometry and Day/Night Photolysis — COMPLETE

**Goal:** Compute per-cell solar zenith angle (SZA) from MPAS model time and
geographic coordinates. Prefer the existing MPAS radiation `coszr` diagnostic
when it is available; otherwise compute the same solar geometry in chemistry.
Replace the constant j_NO2 with a SZA-dependent scaling:
`j_NO2 = j_max * max(0, cos(SZA))`. Validate day/night behavior.

**Status:** Complete and merged to `main` for the accepted fallback-only
idealized path. This section remains as implementation/reference history;
Phase 2 is the active target.

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

Implementation status note: the tracked idealized
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

Implementation status note: the accepted Phase 1 scope is the fallback helper
path. Reuse of MPAS `coszr` and real-data SZA from mesh `latCell/lonCell`
remains deferred until the required mesh-coordinate data is available in the
chemistry workflow.

### Implementation Checklist

Phase 1 is complete. Any remaining unchecked items below are deferred future
cleanup, grid-aware extensions, or extra validation, not blockers for Phase 2.

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

Phase 1 is accepted complete on `main`. Any remaining unchecked items below are
deferred follow-on validation, not gating items for Phase 2.

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

## Phase 2: TUV-x Coupled Photolysis — COMPLETE

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

Key modules (MUSICA-level API — used instead of raw `tuvx_core`):
- `musica_tuvx.mod` — `tuvx_t` type: constructor, `run()`, `get_grids()`,
  `get_profiles()`, `get_photolysis_rate_constants_ordering()`
- `musica_tuvx_grid.mod` — `grid_t`: from-host height grid with
  `set_edges()` / `set_midpoints()`
- `musica_tuvx_grid_map.mod` — `grid_map_t`: register grids before init
- `musica_tuvx_profile.mod` — `profile_t`: from-host profiles with
  `set_edge_values()` / `set_midpoint_values()`
- `musica_tuvx_profile_map.mod` — `profile_map_t`: register profiles
- `musica_tuvx_radiator_map.mod` — `radiator_map_t`: (empty in Fortran init;
  air and O3 radiators defined in JSON config)
- `musica_util.mod` — `error_t`, `mappings_t`, `musica_dk`

### TUV-x API Pattern

We use the MUSICA Fortran wrapper (`musica_tuvx`), not the raw `tuvx_core`.
All atmospheric profiles come from the host model — no static atmosphere
data files. The only external files are spectral physics data (cross-sections,
quantum yields, solar flux, wavelength grid).

```fortran
! Initialization (once) — from-host grids and profiles
heights => grid_t("height", "km", nVertLevels, error)
grids   => grid_map_t(error)
call grids%add(heights, error)

air  => profile_t("air",  "molecule cm-3", heights, error)
temp => profile_t("temperature", "K", heights, error)
o3   => profile_t("O3",   "molecule cm-3", heights, error)
o2   => profile_t("O2",   "molecule cm-3", heights, error)
profiles => profile_map_t(error)
call profiles%add(air, error)
call profiles%add(temp, error)
call profiles%add(o3, error)
call profiles%add(o2, error)

radiators => radiator_map_t(error)
tuvx => tuvx_t(config_path, grids, profiles, radiators, error)

! Retrieve updatable handles after construction
grids    => tuvx%get_grids(error)
profiles => tuvx%get_profiles(error)
heights  => grids%get("height", "km", error)
air      => profiles%get("air", "molecule cm-3", error)
temp     => profiles%get("temperature", "K", error)
o3       => profiles%get("O3", "molecule cm-3", error)
o2       => profiles%get("O2", "molecule cm-3", error)

! Each timestep, per column
call heights%set_edges(zgrid_km)                    ! MPAS zgrid [m] → [km]
call heights%set_midpoints(zmid_km)                 ! layer midpoints
call air%set_midpoint_values(air_number_density)    ! [molecule cm⁻³]
call air%set_edge_values(air_edges)                 ! interface values
call temp%set_midpoint_values(temperature)          ! [K]
call temp%set_edge_values(temp_edges)               ! interface values
call o3%set_midpoint_values(o3_number_density)      ! [molecule cm⁻³]
call o2%set_midpoint_values(o2_number_density)      ! = 0.2095 × air
call tuvx%run(sza_radians, earth_sun_distance_AU,   &
              photo_rate_constants, heating_rates, error)
! photo_rate_constants(nVertLevels+1, nReactions) — extract j_NO2 by index
```

### From-Host Profile Conversions

All atmospheric profiles are derived from MPAS state at runtime. No static
atmosphere data files (USSA or otherwise) are used.

| TUV-x profile | MPAS source | Conversion |
|---------------|-------------|------------|
| height [km] | `zgrid` [m] | × 1e-3 |
| air [molecule cm⁻³] | `rho_air` [kg m⁻³] | × (Nₐ / M_air) × 1e-6 |
| temperature [K] | theta_m, exner, qv | already computed in `chemistry_from_MPAS` |
| O3 [molecule cm⁻³] | `qO3` [kg kg⁻¹] | × rho_air × (Nₐ / M_O3) × 1e-6 |
| O2 [molecule cm⁻³] | air density | × 0.2095 (fixed volume mixing ratio) |

Constants: Nₐ = 6.02214076e23 mol⁻¹, M_air = 0.02897 kg mol⁻¹,
M_O3 = 0.04800 kg mol⁻¹.

### Spectral Physics Data Files

These are physical constants (laboratory measurements, solar spectrum), not
atmospheric state. They must be shipped to the run directory. All source files
are in `~/EarthSystem/TUV-x/data/`.

| File | Purpose |
|------|---------|
| `cross_sections/NO2_1.nc` | NO2 absorption cross-section (T-interpolated) |
| `quantum_yields/NO2_1.nc` | NO2 quantum yield (T-interpolated) |
| `cross_sections/O3_1.nc` through `O3_4.nc` | O3 absorption (for RT optical depth) |
| `cross_sections/O2_1.nc` | O2 absorption cross-section |
| `cross_sections/O2_parameters.txt` | O2 Lyman-α / Schumann-Runge parameters |
| `grids/wavelength/cam.csv` | Wavelength grid for spectral integration |
| `profiles/solar/susim_hi.flx` | Extraterrestrial solar flux (4 files) |
| `profiles/solar/atlas3_1994_317_a.dat` | " |
| `profiles/solar/sao2010.solref.converted` | " |
| `profiles/solar/neckel.flx` | " |

### New Source Files

**`src/core_atmosphere/chemistry/mpas_tuvx.F`**

New module wrapping TUV-x for CheMPAS. Uses the MUSICA-level API
(`musica_tuvx`) with all atmospheric profiles from host:

```fortran
module mpas_tuvx
    ! Module-level state (initialized once, updated per column)
    type(tuvx_t),       pointer :: tuvx_solver
    type(grid_t),       pointer :: height_grid
    type(profile_t),    pointer :: air_profile, temp_profile
    type(profile_t),    pointer :: o3_profile, o2_profile
    type(mappings_t),   pointer :: photo_ordering, heating_ordering
    integer                     :: j_no2_index   ! index in photo_rate_constants
    integer                     :: n_photo_rates, n_heating_rates
    real(dk), allocatable       :: photo_rates(:, :), heating_rates(:, :)
contains
    subroutine tuvx_init(config_file, nVertLevels, error_code, error_message)
        ! Create from-host grids (height) and profiles (air, T, O3, O2)
        ! Construct tuvx_t with config file + from-host objects
        ! Retrieve updatable handles via get_grids/get_profiles
        ! Look up j_NO2 reaction index via get_photolysis_rate_constants_ordering
        ! Query heating/photo counts and allocate persistent work arrays once

    subroutine tuvx_compute_photolysis(cos_sza, z_edges_m, rho_air, &
                                        temperature, qo3, nVertLevels, &
                                        j_no2_column, error_code, error_message)
        ! Convert MPAS fields to TUV-x units (m→km, kg/m³→molec/cm³, etc.)
        ! Update height edges/midpoints, air, T, O3, O2 profiles
        ! Call tuvx%run(sza_rad, esd_au, photo_rates, heating_rates, error)
        ! Extract j_no2_column(1:nVertLevels) from photo_rates(:, j_no2_index)

    logical function tuvx_is_enabled()
        ! Returns .true. if tuvx_solver is associated

    subroutine tuvx_finalize()
end module
```

**TUV-x config file** (`micm_configs/tuvx_no2.json`)

Minimal TUV-x configuration for NO2 photolysis only. All atmospheric profiles
are from-host (no static atmosphere files). The config references only spectral
physics data:

- **Grids section:** empty (`"grids": []`) — height grid registered from host
- **Profiles:** surface albedo (constant 0.1 in config), extraterrestrial solar
  flux (from data files). Air, temperature, O3, O2 are all from host.
- **Radiative transfer:** delta-Eddington (2-stream) solver, radiators for air
  (Rayleigh), O2, O3
- **Photolysis:** single reaction — NO2 + hv → NO + O(3P) using `NO2 tint`
  cross-section and quantum yield (temperature-interpolated, from NetCDF data)
- **Clouds/aerosols:** not included (clear-sky only, Phase 2 scope)

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

2. **All profiles from host.** No static atmosphere data files. Height edges
   from `zgrid`, air number density derived from `rho_air`, temperature from
   the existing `chemistry_from_MPAS` thermodynamic calculation, O3 from the
   model's own `qO3` tracer (self-consistent), and O2 derived as 0.2095 × air.
   Unit conversions (kg/m³ → molecule/cm³, m → km) happen in `mpas_tuvx.F`.
   For absorbers, host-side midpoint values alone are not enough: TUV-x layer
   densities must also be updated consistently from the host column.

3. **MUSICA-level API.** Use `musica_tuvx` (`tuvx_t`, `grid_t`, `profile_t`)
   rather than the raw `tuvx_core` module. This gives us a clean C-binding
   interface with error handling via `error_t`. One `tuvx_t` instance per MPI
   rank; MPAS already distributes blocks across ranks.

4. **Rate parameter update.** The j_NO2 array (nVertLevels × nCells) replaces
   the scalar `config_lnox_j_no2`. Rate parameters updated every step via the
   same `PHOTO.no2_photolysis` index, but now with per-cell/level values.

5. **Fallback.** If `config_tuvx_config_file` is empty, fall back to Phase 1
   behavior (cos_sza scaling with `config_lnox_j_no2`). This preserves the
   simpler mode for quick testing.

6. **Validation scope.** Phase 2 remains clear-sky and in-domain. It does not
   yet include cloud optics, aerosol optics, or an atmosphere above the MPAS
   model top, so absolute magnitude checks are sanity checks rather than hard
   pass/fail literature matches.

### Implementation Checklist

- [x] Add `config_tuvx_config_file` to `Registry.xml` and read it in
  `chemistry_init`.
- [x] Create `micm_configs/tuvx_no2.json` — minimal config with NO2 photolysis,
  delta-Eddington solver, from-host grids/profiles, spectral data file paths.
  **Note:** O2 radiator removed — O2 cross-section `"boundary"` extrapolation
  extended absorption into 300-420 nm (j_NO2 wavelengths), causing 4400x
  artificial attenuation. O2 absorption is only relevant for <240 nm reactions.
- [x] Add `mpas_tuvx.F` with `tuvx_init`, `tuvx_compute_photolysis`,
  `tuvx_is_enabled`, `tuvx_finalize`. Uses `musica_tuvx` API with from-host
  `grid_t` and `profile_t` objects for height, air, temperature, O3, O2.
- [x] Update `chemistry/Makefile` — add `mpas_tuvx.o` with MUSICA module deps,
  and add the explicit `mpas_atm_chemistry.o: ... mpas_tuvx.o` dependency so
  build ordering remains correct in this hand-maintained Makefile.
- [x] Initialize TUV-x once during `chemistry_init` when
  `config_tuvx_config_file` is non-empty; otherwise leave the TUV-x pointer
  unassociated and use the Phase 1 fallback branch.
- [x] Resolve and cache `index_qO3` during chemistry initialization, and fail
  fast with a clear log message if the runtime tracer set does not include
  `qO3`. Phase 2 assumes host O3 is present; it should not silently run with a
  missing tracer.
- [ ] Refactor `chemistry_from_MPAS` so the environment-building loop
  (`rho_air`, `temperature`, `pressure`) is reusable by both MICM coupling and
  TUV-x, rather than recomputing those fields twice. **Deferred:** currently
  the TUV-x branch computes its own `rho_air` and `temperature` per column;
  works but duplicates effort.
- [x] Add a chemistry-side workspace for `j_no2(:,:)` with dimensions
  `(nVertLevels, nCells)` for the current block.
- [x] In `mpas_tuvx.F`, query photolysis and heating orderings/counts during
  `tuvx_init`, cache the `j_NO2` reaction index, and allocate persistent
  contiguous work arrays sized to the actual TUV-x interface:
  `photo_rates(nVertLevels+1, n_photo_rates)` and
  `heating_rates(nVertLevels+1, n_heating_rates)`.
- [x] In `chemistry_step`, add the TUV-x branch: loop over cells, call
  `tuvx_compute_photolysis` per column with from-host conversions:
  - height: `zgrid` [m] → [km]
  - air: `rho_air` [kg/m³] → [molecule/cm³] via `× (Nₐ / M_air) × 1e-6`
  - O3: `qO3` [kg/kg] → [molecule/cm³] via `× rho_air × (Nₐ / M_O3) × 1e-6`
  - O2: `= 0.2095 × air` [molecule/cm³]
  - temperature: [K] (no conversion)
- [x] For host-supplied absorber profiles (`air`, `O2`, `O3`), update not just
  midpoint/edge values but also TUV-x `layer_densities` consistently from the
  height grid. Do not assume midpoint values alone are sufficient.
- [x] In `mpas_musica.F`, add `musica_set_photolysis_field` that writes the
  full `j_no2(level, cell)` array into `PHOTO.no2_photolysis` rate parameters.
- [x] Keep the fallback behavior: Phase 1 and Phase 2 both use the same
  `musica_set_photolysis_*` routines; fallback uses `_scalar`, TUV-x uses
  `_field`.
- [x] Copy TUV-x spectral data files to run directory (cross-sections, quantum
  yields, solar flux, wavelength grid, O2 parameters). **Note:** TUV-x reads
  data files from paths relative to the run directory via the `data/` symlink
  already present in the supercell test case.
- [x] Update `j_no2` diagnostic to write the per-cell/level TUV-x output
  (currently writes a uniform scalar). Added `chemistry_set_j_no2_diag_field`.
- [x] Extend the phase-gate scripts with `fallback-compare`,
  `transition-smooth`, and decomposition checks. The imported ancestor scripts
  are now adapted to the CheMPAS Phase 0/1/2 matrix.
- [x] Link flang-built `libnetcdff.a` for MUSICA's `musica_io_netcdf` module.
  Homebrew's gfortran-built `libnetcdff` uses incompatible symbol mangling
  (`___netcdf_MOD_*` vs `__QMnetcdf*`). Static archive path specified directly
  in Makefile to avoid linker picking up the wrong library.

### Key Findings

- **O2 radiator bug:** TUV-x cross-section `"upper extrapolation": "boundary"`
  extends the last data point value to all longer wavelengths. For O2 data
  covering 175–240 nm, this incorrectly applied O2 absorption at 300–420 nm
  (j_NO2 range), creating optical depth ≈ 4 through the full air column and
  reducing surface j_NO2 by 4400x. Fix: remove O2 from radiators list. If O2
  photolysis reactions are added later, use `"type": "zero"` extrapolation.
- **Flang/gfortran symbol mangling:** MUSICA's internal NetCDF module uses
  flang name mangling (`__QMnetcdfPnf90_open`), incompatible with Homebrew's
  gfortran-built `libnetcdff` (`___netcdf_MOD_nf90_open`). Required a
  flang-built `libnetcdff.a` from the MUSICA build tree.
- **Wavelength grid:** The spectral wavelength grid is NOT from-host; it must
  be defined in the JSON config (from CSV file). Only the height grid is
  from-host.

### High-Level Fortran Design

Phase 2 is easiest to keep coherent if the chemistry driver owns the column
work arrays and chooses between two photolysis providers:

1. Phase 1 scalar fallback: `j = j_max * max(0, cos_sza)`
2. Phase 2 TUV-x provider: `j = TUV-x(z, rho, T, O3, SZA)` per column

That keeps MICM blind to where the photolysis field came from.

```fortran
! mpas_atm_chemistry.F  —  chemistry_step
!
! Compute SZA (same as Phase 1)
cos_sza = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)

! Build environment arrays (refactored from chemistry_from_MPAS)
! rho_air(nVertLevels, nCells), temperature(:,:), pressure(:,:) already computed
! Also extract zgrid(:,:) and qO3(:,:) from mesh/state pools
! idx_qO3 is resolved once during chemistry_init; fail fast if missing

if (tuvx_is_enabled()) then
    do iCell = 1, nCells
        call tuvx_compute_photolysis(cos_sza,                         &
                                     zgrid(:, iCell),                 &
                                     rho_air(:, iCell),               &
                                     temperature(:, iCell),           &
                                     scalars(idx_qO3, :, iCell),      &
                                     nVertLevels,                     &
                                     j_no2(:, iCell),                 &
                                     error_code, error_message)
        if (error_code /= 0) return
    end do
    call musica_set_photolysis_field(j_no2, nCells, nVertLevels, &
                                     error_code, error_message)
else
    j_no2_value = chem_j_no2_max * max(0.0_RKIND, cos_sza)
    call musica_set_photolysis(j_no2_value, error_code, error_message)
end if

! Then proceed with MICM coupling as before
call MICM_from_chemistry(scalars, rho_air, temperature, pressure, ...)
call musica_step(time_step, ...)
```

```fortran
! mpas_tuvx.F  —  per-column TUV-x solve with all from-host profiles
subroutine tuvx_compute_photolysis(cos_sza, z_edges_m, rho_air, temperature, &
                                    qo3, nVertLevels,                         &
                                    j_no2_column, error_code, error_message)
    use musica_util, only: dk => musica_dk, error_t

    real(kind=RKIND), intent(in)  :: cos_sza
    real(kind=RKIND), intent(in)  :: z_edges_m(:)       ! nVertLevels+1 edges
    real(kind=RKIND), intent(in)  :: rho_air(:)          ! nVertLevels midpoints
    real(kind=RKIND), intent(in)  :: temperature(:)      ! nVertLevels midpoints
    real(kind=RKIND), intent(in)  :: qo3(:)              ! nVertLevels [kg/kg]
    integer, intent(in)           :: nVertLevels
    real(kind=RKIND), intent(out) :: j_no2_column(:)     ! nVertLevels
    integer, intent(out)          :: error_code
    character(len=:), allocatable, intent(out) :: error_message

    real(dk) :: z_km(nVertLevels+1), z_mid_km(nVertLevels)
    real(dk) :: air_nd(nVertLevels), air_edges(nVertLevels+1)
    real(dk) :: temp_mid(nVertLevels), temp_edges(nVertLevels+1)
    real(dk) :: o3_nd(nVertLevels), o2_nd(nVertLevels)
    real(dk) :: air_layer(nVertLevels), o3_layer(nVertLevels), o2_layer(nVertLevels)
    real(dk) :: sza_rad, esd
    type(error_t) :: err

    ! Convert MPAS fields to TUV-x units
    z_km     = z_edges_m * 1.0e-3_dk                       ! m → km
    z_mid_km = 0.5_dk * (z_km(1:nVertLevels) + z_km(2:nVertLevels+1))
    air_nd   = rho_air * (NA / M_AIR) * 1.0e-6_dk          ! kg/m³ → molec/cm³
    o3_nd    = qo3 * rho_air * (NA / M_O3) * 1.0e-6_dk     ! kg/kg → molec/cm³
    o2_nd    = 0.2095_dk * air_nd                           ! fixed VMR
    temp_mid = temperature
    ! Interpolate edge values (simple average of adjacent midpoints)
    air_edges(1) = air_nd(1); air_edges(nVertLevels+1) = air_nd(nVertLevels)
    air_edges(2:nVertLevels) = 0.5_dk * (air_nd(1:nVertLevels-1) + air_nd(2:nVertLevels))
    ! ... same pattern for temp_edges ...
    air_layer = air_nd * (z_km(2:nVertLevels+1) - z_km(1:nVertLevels)) * 1.0e5_dk
    o3_layer  = o3_nd  * (z_km(2:nVertLevels+1) - z_km(1:nVertLevels)) * 1.0e5_dk
    o2_layer  = o2_nd  * (z_km(2:nVertLevels+1) - z_km(1:nVertLevels)) * 1.0e5_dk

    ! Update from-host profiles
    call height_grid%set_edges(z_km, err);           ! check err
    call height_grid%set_midpoints(z_mid_km, err)
    call air_profile%set_midpoint_values(air_nd, err)
    call air_profile%set_edge_values(air_edges, err)
    call temp_profile%set_midpoint_values(temp_mid, err)
    call temp_profile%set_edge_values(temp_edges, err)
    call o3_profile%set_midpoint_values(o3_nd, err)
    call air_profile%set_layer_densities(air_layer, err)
    call o3_profile%set_layer_densities(o3_layer, err)
    call o2_profile%set_midpoint_values(o2_nd, err)
    call o2_profile%set_layer_densities(o2_layer, err)

    ! SZA: TUV-x expects radians; nighttime → skip (j=0)
    if (cos_sza <= 0.0_RKIND) then
        j_no2_column = 0.0_RKIND
        return
    end if
    sza_rad = acos(cos_sza)
    esd = 1.0_dk  ! Earth-sun distance [AU], ~1.0 for short runs

    ! photo_rates/heating_rates were allocated once during tuvx_init using the
    ! actual TUV-x ordering/count metadata and must remain contiguous
    call tuvx_solver%run(sza_rad, esd, photo_rates, heating_rates, err)

    ! Extract j_NO2 from photo_rates (nVertLevels+1 edges → nVertLevels midpoints)
    j_no2_column(1:nVertLevels) = 0.5_RKIND * &
        (photo_rates(1:nVertLevels, j_no2_index) + &
         photo_rates(2:nVertLevels+1, j_no2_index))
end subroutine
```

```fortran
! mpas_musica.F  —  new per-field photolysis update
subroutine musica_set_photolysis_field(j_no2, nCells, nVertLevels, &
                                       error_code, error_message)
    real(kind=RKIND), intent(in) :: j_no2(nVertLevels, nCells)
    integer, intent(in)          :: nCells, nVertLevels
    integer, intent(out)         :: error_code
    character(len=:), allocatable, intent(out) :: error_message

    ! Write j_no2(k, iCell) into MICM rate parameter array using
    ! cached photo_no2_rp_index and stride-based indexing
    ! (same stride logic as musica_set_photolysis but per cell/level)
end subroutine
```

The key refactor is to turn "set rate parameters" into a stable API with two
front doors:

- `musica_set_photolysis(j_no2_scalar, ...)` for Phase 1 and fallback mode
- `musica_set_photolysis_field(j_no2_array, ...)` for TUV-x mode

That prevents the chemistry driver from knowing MICM's stride details, and it
keeps the Phase 2 fallback path cheap to maintain.

### Verification (Phase 2 Gate)

| Check | Criterion | Result |
|-------|-----------|--------|
| Non-negativity | qNO, qNO2, qO3 >= 0 | **PASS** — all tracers non-negative |
| Night j-zero | j_NO2 = 0 when SZA >= 90° | **PASS** — nighttime shortcut in `tuvx_compute_photolysis` |
| Vertical structure | j_NO2 increases with altitude (less attenuation) | **PASS** — monotonic increase from 7.2e-3 (surface) to 1.2e-2 (20 km) |
| Magnitude sanity | Surface j_NO2 in clear-sky range (~0.005–0.01 s⁻¹) | **PASS** — 7.2e-3 s⁻¹ at SZA≈59° (literature ~0.008) |
| Surface/top ratio | Ratio 0.5–0.8 for clear sky | **PASS** — 0.61 |
| O3 background | O3 preserved at 50 ppbv away from storm | **PASS** — 50.000 ppbv at background cells |
| Chemistry response | Stronger j_NO2 at altitude → less NO2 accumulation | **PASS** — NO2 peak 6.5 ppbv (vs 8.5 Phase 1) |
| Build passes | MUSICA=true with TUV-x module linked | **PASS** |
| Fallback compare | Empty config reproduces Phase 1 behavior | **PASS** — tested manually |
| Transition smooth | Dawn/dusk j_NO2 varies smoothly | Deferred to extended run |
| Ox conservation | Domain-integrated Ox conserved (source/sink off) | Deferred |
| Decomp compare | Identical results across MPI decompositions | Deferred |

### Test Results (15-minute supercell, Case B)

Configuration: `config_tuvx_config_file = 'tuvx_no2.json'`, SZA ≈ 59°
(35.86°N, 97.93°W, 18:00 UTC), source = 0.5 ppbv/s, O3 init = 50 ppbv.

| Metric | Phase 1 (cos SZA) | Phase 2 (TUV-x) |
|--------|-------------------|-----------------|
| j_NO2 surface | 5.1e-3 (uniform) | 7.2e-3 s⁻¹ |
| j_NO2 10 km | 5.1e-3 (uniform) | 1.0e-2 s⁻¹ |
| j_NO2 20 km | 5.1e-3 (uniform) | 1.2e-2 s⁻¹ |
| NO peak | 27.9 ppbv | 29.9 ppbv |
| NO2 peak | 8.5 ppbv | 6.5 ppbv |
| O3 min | 41.5 ppbv | 43.5 ppbv |

The chemistry responds correctly: stronger photolysis at altitude recycles
NO2 → NO + O3 faster, reducing NO2 accumulation and O3 depletion.

### Exit Criteria

- [x] Build passes with TUV-x module linked.
- [x] j_NO2 profile shows plausible clear-sky vertical structure (monotonic
  increase with altitude, ratio 0.61 surface-to-top).
- [x] Surface j_NO2 = 7.2e-3 s⁻¹ at SZA≈59° matches literature clear-sky
  (caveats: 20 km domain top, no stratospheric O3 column, clear-sky only).
- [x] 15-min Case B remains stable with TUV-x-provided j_NO2.
- [x] Fallback: empty `config_tuvx_config_file` uses Phase 1 path.
- [ ] Deferred: `transition-smooth` for extended runs spanning sunset.
- [ ] Deferred: `decomp-compare` across MPI decompositions.
- [ ] Deferred: Ox conservation test with transport disabled.

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

Status note: both scripts have now been copied from the
ancestor repo into `scripts/` and adapted to the current CheMPAS species and
phase matrix. Remaining gaps are narrower: the analytical-SZA Phase 1 check is
still a future enhancement rather than a Phase 2 blocker, and runtime use
still depends on the Python `netCDF4` module plus whichever `j_*` diagnostics
are written to output files.

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

1. **All atmospheric profiles from host** — Height, air density, temperature,
   O3, O2 all derived from MPAS state at runtime. No static atmosphere data
   files (USSA or otherwise). Only spectral physics data (cross-sections,
   quantum yields, solar flux) come from external files.
2. **Domain-top limitation** — 20 km top omits part of the overhead ozone
   column, so Phase 2 magnitude checks are sanity checks rather than strict
   literature targets.
3. **Source/sink representation split** — Lightning-NOx source is operator-split
   pre-MICM; NOx sink remains mechanism-defined within MICM.
4. **Cloud/aerosol optics deferred** — Phase 2 is clear-sky only until cloud
   and aerosol optical inputs are coupled into TUV-x.
5. **MUSICA-level API** — Use `musica_tuvx` module (C-binding wrapper), not
   the raw `tuvx_core` module. One `tuvx_t` instance per MPI rank.

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
- TUV-x spectral data files (from `~/EarthSystem/TUV-x/data/`):
  - `cross_sections/NO2_1.nc`, `O3_1-4.nc`, `O2_1.nc`, `O2_parameters.txt`
  - `quantum_yields/NO2_1.nc`
  - `grids/wavelength/cam.csv`
  - `profiles/solar/susim_hi.flx`, `atlas3_1994_317_a.dat`,
    `sao2010.solref.converted`, `neckel.flx`
- Python `netCDF4` for verification scripts
