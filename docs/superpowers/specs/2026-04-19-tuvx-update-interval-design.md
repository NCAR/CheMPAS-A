# TUV-x Update Interval — Design

Date: 2026-04-19
Status: Design (awaiting user review)
Target files: `src/core_atmosphere/Registry.xml`,
              `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`

## Goal

Throttle TUV-x updates to a configurable simulated-time interval. A new
namelist option `config_tuvx_update_interval` (real, seconds, default
`0.0`) gates the entire photolysis block in `chemistry_step` —
`cos_sza_cell` allocation + per-cell column extraction +
`tuvx_compute_photolysis` calls + `musica_set_photolysis_rates` +
`chemistry_set_photolysis_diag`.

When the throttle is "off" (default `0.0`), behavior is bit-identical to
the current per-cell-lat/lon baseline at commit `d7601e4`. When the
throttle is "on" (positive interval), TUV-x runs only when at least
`interval` seconds of simulated time have accumulated since the last
update; on skipped chemistry steps, MICM keeps using the rate parameters
last set (MICM stores them internally), and the `j_<rxn>` diagnostic
fields retain their last-written values.

## Scope

**In scope:**
- New namelist option in `&musica`: `config_tuvx_update_interval`
- Two module-scope save variables for the accumulator:
  `tuvx_update_interval` and `tuvx_time_since_last`
- A single `if/else` gate around the existing photolysis block in
  `chemistry_step`
- Doc update describing the new option

**Out of scope:**
- Adaptive cadence (e.g., faster updates near sunrise/sunset)
- Per-rate or per-mechanism throttling (one global interval covers all
  photolysis reactions)
- Restart-file persistence of the accumulator (module `save` re-init to
  `huge` correctly forces an update on the first chemistry step after
  restart, matching MICM's zeroed rate parameters)
- Diagnostic field "freshness" handling beyond "leave at last-written
  values" (no separate update timestamp on the diag)
- Wall-clock timer instrumentation beyond what already exists in
  `chemistry_step`

## Pre-existing State

After commit `d7601e4` (per-cell lat/lon switch), `chemistry_step` in
`src/core_atmosphere/chemistry/mpas_atm_chemistry.F`:

- Calls `mpas_get_time` and computes `hour_utc` once per step
- Allocates `cos_sza_cell(nCells)`, fills it (per-cell or broadcast),
  allocates `photo_rates(n_photo_rp, nVertLevels, nCells)`, runs the
  TUV-x branch or fallback branch, calls
  `musica_set_photolysis_rates`, calls `chemistry_set_photolysis_diag`,
  deallocates both arrays
- Does the above on **every** chemistry step

`musica_set_photolysis_rates` overwrites MICM's internal photolysis
rate parameters. If not called, MICM keeps the last-set values across
subsequent `musica_step` calls.

The `j_<rxn>` diagnostic fields (e.g., `j_jNO2`) are MPAS history
variables. Their values persist between writes if not overwritten.

## Design

### Namelist & Registry

Add one new namelist option to the existing `&musica` record in
`src/core_atmosphere/Registry.xml`, immediately after the
`config_chemistry_use_grid_coords` block (which currently lives just
after `config_chemistry_longitude`):

```xml
                <nml_option name="config_tuvx_update_interval" type="real" default_value="0.0"
                     units="s"
                     description="TUV-x update interval in simulated seconds. 0.0 means update every chemistry step (default; preserves bit-reproducibility). Positive values gate the photolysis block: TUV-x runs only when at least this many simulated seconds have accumulated since the last update; MICM reuses the last-set photolysis rates on skipped steps."
                     possible_values="Any non-negative real"/>
```

Default `0.0` → "update every step" → exact preservation of current
behavior. Negative values are not explicitly checked; the implementation
treats `interval <= 0.0` as "always update" because the accumulator
comparison `huge >= interval` is always true on the first call and
`0.0 >= interval` is true thereafter.

### Module-scope state

Add two module-scope save variables to
`src/core_atmosphere/chemistry/mpas_atm_chemistry.F`, alongside the
existing `chem_lat` / `chem_lon` / `chem_j_no2_max` /
`chem_use_grid_coords` group (currently around lines 37–40):

```fortran
    ! TUV-x update interval state (set in chemistry_init, accumulated in chemistry_step)
    real(kind=RKIND), save :: tuvx_update_interval = 0.0_RKIND
    real(kind=RKIND), save :: tuvx_time_since_last = huge(1.0_RKIND)
```

`tuvx_time_since_last` initialized to `huge` so the very first
`chemistry_step` call always fires an update. After process restart, the
`save` re-initializes to `huge` (Fortran `save` semantics within a
process), which is the desired behavior: MICM's rate parameters are
zeroed at init, so the first post-restart chemistry step must run TUV-x.

### `chemistry_init` — read the new config

In `chemistry_init`, after the existing read of
`config_chemistry_use_grid_coords` (around lines 167–169), add:

```fortran
        block
            real(kind=RKIND), pointer :: tuvx_update_interval_ptr
            nullify(tuvx_update_interval_ptr)
            call mpas_pool_get_config(configs, 'config_tuvx_update_interval', tuvx_update_interval_ptr)
            if (associated(tuvx_update_interval_ptr)) tuvx_update_interval = tuvx_update_interval_ptr
        end block
```

The `block` construct (rather than adding a top-of-routine pointer
declaration) follows the existing pattern in this file for
`config_chemistry_ref_solve`, `config_chem_substeps`, and
`config_micm_relative_tolerance` (lines 133–155).

### `chemistry_step` — gate the photolysis block

The throttle wraps the existing photolysis block in a single
`if/else` block. The hard-error `n_photo_rp < 1` check stays outside
the gate (it's a process-fatal condition). The `mpas_get_time` and
`hour_utc` computation moves inside the gate (only used inside).

After the change, the top of `chemistry_step` reads:

```fortran
        if (n_photo_rp < 1) then
            call mpas_log_write('[Chemistry] photolysis rate indices not cached.', &
                messageType=MPAS_LOG_CRIT)
            return
        end if

        ! TUV-x throttle: run photolysis only if at least tuvx_update_interval
        ! simulated seconds have accumulated since the last update. On skip
        ! steps, MICM keeps using the rate parameters last set, and the
        ! j_<rxn> diag fields keep their last-written values.
        if (tuvx_time_since_last >= tuvx_update_interval) then
            call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
            hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                     + real(second, RKIND) / 3600.0_RKIND

            allocate(cos_sza_cell(nCells))
            if (chem_use_grid_coords) then
                ! ... per-cell fill (unchanged) ...
            else
                cos_sza_cell(:) = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)
            end if

            allocate(photo_rates(n_photo_rp, nVertLevels, nCells))
            photo_rates = 0.0_RKIND

            if (use_tuvx) then
                ! ... per-cell TUV-x loop (unchanged) ...
            else
                ! ... fallback loop (unchanged) ...
            end if

            call musica_set_photolysis_rates(photo_rates, ...)
            if (error_code /= 0) then
                deallocate(photo_rates, cos_sza_cell)
                return
            end if

            call chemistry_set_photolysis_diag(diag, photo_rates, ...)
            deallocate(photo_rates, cos_sza_cell)

            tuvx_time_since_last = 0.0_RKIND
        else
            tuvx_time_since_last = tuvx_time_since_last + dt
        end if

        ! ... lightning_nox, MPAS->MICM, MICM step, MICM->MPAS (unchanged) ...
```

### Correctness analysis

1. **Default behavior preserved.** With `tuvx_update_interval = 0.0` and
   `tuvx_time_since_last = huge` initially, the comparison
   `huge >= 0.0` is true on the first call → update fires →
   `tuvx_time_since_last = 0.0`. Next call: `0.0 >= 0.0` is true →
   update fires again. The default path runs the photolysis block every
   chemistry step, exactly as today.

2. **First call always updates.** The `huge` initial value of
   `tuvx_time_since_last` guarantees this regardless of `interval`.

3. **Restart correctness.** Module `save` re-initializes to `huge` in a
   fresh process → first chemistry_step call after restart updates →
   matches the fact that MICM's rate parameters are zeroed at init.

4. **Early-return paths inside the photolysis block do not corrupt the
   accumulator.** The TUV-x per-column error path and the
   `musica_set_photolysis_rates` error path both `return` from
   `chemistry_step` without resetting `tuvx_time_since_last`. So a
   transient TUV-x error on what would have been an update step leaves
   the accumulator unchanged (still `>= interval`), and the next call
   will re-attempt the update.

5. **The accumulator only ticks up on skip steps.** It does *not* tick
   on update steps (we set it to `0.0` at the end of the update block).
   This matches the natural interpretation of "time since the last
   update."

### Validation

1. **Default behavior preserved (`interval = 0.0`):** Run the supercell
   Phase 2/3 case with the new code and the default
   `config_tuvx_update_interval = 0.0`. Spot-check `j_jNO2`, `qNO`,
   `qNO2`, `qO3` against the values recorded after the lat/lon-work run
   (clean reference at `output.nc` from `d7601e4`). Expected:
   bit-for-bit identical.

2. **Throttle actually skips work (`interval = 60.0`):** Run the same
   supercell case with `config_tuvx_update_interval = 60.0` (one update
   every 20 dynamics steps at `dt = 3 s`). Verify two things:

   - **Wall-clock improvement.** TUV-x dominates per-step chemistry
     cost; the throttled run should show a noticeably lower per-step
     wall time in the chemistry timer.
   - **Update cadence.** Add temporary instrumentation: a single
     `mpas_log_write('[TUV-x] photolysis update fired at xtime ' // ...)`
     line inside the update branch. Run for 3 simulated min and confirm
     exactly **4** update messages: at t=0 (initial), t=60s, t=120s,
     t=180s. Remove the temporary log line before commit.

3. **Physical sanity (`interval = 60.0`):** `j_jNO2` min/max at the end
   of the throttled run should be within ~1% of the unthrottled run's
   values at the same simulated time. Slow SZA evolution over 60 s
   means the cached value is a very small staleness error — sub-percent
   at any time of day except dawn/dusk transitions.

## Files Touched

**Modified (2):**
- `src/core_atmosphere/Registry.xml` (one new namelist option in
  `&musica`)
- `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (two
  module-scope save vars, one `block` config-read in `chemistry_init`,
  one `if/else` wrapper around the photolysis block in
  `chemistry_step`, plus moving `mpas_get_time` + `hour_utc` lines
  inside the gate)

**Doc update:**
- `docs/guides/TUVX_INTEGRATION.md` — add a new `### TUV-x update
  interval` subsection in "Multi-Photolysis Plumbing" (parallel in
  shape to the recently-added `### Per-cell solar geometry`
  subsection), describing the new switch, the default behavior, and a
  typical operational value (60 s for short-duration / fine-dt cases;
  longer for global runs).

**Not touched:**
- `mpas_solar_geometry.F` (signature unchanged)
- `mpas_tuvx.F` (signature unchanged)
- `mpas_musica.F`, `mpas_lightning_nox.F` (no interaction with
  throttling)
- `mpas_atm_core.F` (calls chemistry_step every dynamics step
  unchanged)
- MICM YAML configs, TUV-x JSON configs, namelist defaults in
  `test_cases/`

## Commit Shape

Single commit on `develop`:

```
feat(chemistry): throttle TUV-x to a configurable update interval
```

Body explains the new `config_tuvx_update_interval` switch, default
`0.0` = every step, and the bit-reproducibility guarantee at default.
The doc update lands in the same commit.

## Sequencing

This is the second of two sequential TUV-x PRs (per the user's choice
in the prior brainstorming, option B). The first (per-cell lat/lon)
landed in commit `d7601e4` and is currently on `develop` awaiting
Ubuntu cross-platform validation. Implementation of this spec should
wait until the Ubuntu reproducibility check confirms the prior change
is clean on both platforms — landing throttling on a confirmed-working
baseline keeps debugging simpler if anything regresses.
