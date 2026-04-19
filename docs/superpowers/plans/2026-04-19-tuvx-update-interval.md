# TUV-x Update Interval Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `config_tuvx_update_interval` namelist switch (real, seconds, default `0.0`) so chemistry_step gates the entire photolysis block — only running TUV-x and updating MICM rate parameters when at least `interval` simulated seconds have accumulated since the last update.

**Architecture:** Two new module-scope save variables track the configured interval and the accumulator. A single `if/else` block in `chemistry_step` wraps the existing photolysis section: on update steps, run as today and reset accumulator to zero; on skip steps, add `dt` to accumulator. MICM keeps the last-set photolysis rates internally between updates; the `j_<rxn>` diag fields keep their last-written values. Default `interval = 0.0` plus the `huge` initial value of the accumulator preserve exact bit-reproducibility of the per-cell-lat/lon baseline at commit `d7601e4`.

**Tech Stack:** Fortran 2008 (LLVM/flang on macOS or GCC/gfortran on Ubuntu), MPAS Registry XML, MPI. Spec: `docs/superpowers/specs/2026-04-19-tuvx-update-interval-design.md`.

---

## File Structure

- Modify: `src/core_atmosphere/Registry.xml` (add namelist option)
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (two save vars, one block in `chemistry_init`, one if/else wrapper in `chemistry_step`)
- Modify: `docs/guides/TUVX_INTEGRATION.md` (add usage subsection)

No new files. No deletions. No changes to other source files, MICM YAMLs, TUV-x JSONs, or `test_cases/` namelists.

---

### Task 1: Add `config_tuvx_update_interval` to Registry

**Files:**
- Modify: `src/core_atmosphere/Registry.xml` (insert after the `config_chemistry_use_grid_coords` block, which currently lives at lines 451–453)

- [ ] **Step 1: Locate the existing `config_chemistry_use_grid_coords` block**

Run: `grep -n 'config_chemistry_use_grid_coords' src/core_atmosphere/Registry.xml`
Expected: one line number (currently 451) showing the existing namelist option.

- [ ] **Step 2: Edit `Registry.xml` — insert new option after `config_chemistry_use_grid_coords`**

Find this exact block (currently around lines 451–453):

```xml
                <nml_option name="config_chemistry_use_grid_coords" type="logical" default_value=".false."
                     description="If .true., compute per-cell solar geometry from latCell/lonCell. If .false., broadcast namelist config_chemistry_latitude/longitude to every cell (preserves idealized-case behavior)."
                     possible_values=".true. or .false."/>
```

Add the following block immediately after its closing `/>` (with the same indentation):

```xml
                <nml_option name="config_tuvx_update_interval" type="real" default_value="0.0"
                     units="s"
                     description="TUV-x update interval in simulated seconds. 0.0 means update every chemistry step (default; preserves bit-reproducibility). Positive values gate the photolysis block: TUV-x runs only when at least this many simulated seconds have accumulated since the last update; MICM reuses the last-set photolysis rates on skipped steps."
                     possible_values="Any non-negative real"/>
```

- [ ] **Step 3: Verify the new option is present and well-formed**

Run: `grep -n 'config_tuvx_update_interval' src/core_atmosphere/Registry.xml`
Expected: at least one line showing the new option name.

Run: `xmllint --noout src/core_atmosphere/Registry.xml; echo "exit=$?"`
Expected: `exit=0` (Registry XML is well-formed). If `xmllint` is not installed, skip; the build will catch malformed XML.

---

### Task 2: Add module-scope state and read in `chemistry_init`

**Files:**
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (module-scope save vars around line 40; chemistry_init body around line 169)

- [ ] **Step 1: Add the module-scope save variables**

Locate the existing solar-geometry / lat-lon-switch save variables (currently around lines 37–40):

```fortran
    ! Solar geometry config values (cached at init, used per-step)
    real(kind=RKIND), save :: chem_lat = 0.0_RKIND
    real(kind=RKIND), save :: chem_lon = 0.0_RKIND
    real(kind=RKIND), save :: chem_j_no2_max = 0.0_RKIND
    logical, save :: chem_use_grid_coords = .false.
```

Add the following two-line block IMMEDIATELY AFTER the `chem_use_grid_coords` line (same indentation):

```fortran

    ! TUV-x update interval state (set in chemistry_init, accumulated in chemistry_step)
    real(kind=RKIND), save :: tuvx_update_interval = 0.0_RKIND
    real(kind=RKIND), save :: tuvx_time_since_last = huge(1.0_RKIND)
```

(Leading blank line for readability; mirrors the existing structure between the solar-geometry group and the `use_ref_solve` group around lines 41–48.)

- [ ] **Step 2: Read the new config in `chemistry_init`**

Locate the existing `chem_use_grid_coords` config-read block in `chemistry_init` (currently around lines 167–169):

```fortran
        nullify(use_grid_coords_ptr)
        call mpas_pool_get_config(configs, 'config_chemistry_use_grid_coords', use_grid_coords_ptr)
        if (associated(use_grid_coords_ptr)) chem_use_grid_coords = use_grid_coords_ptr
```

Add the following `block` construct IMMEDIATELY AFTER the last `if (associated(use_grid_coords_ptr)) ...` line, with the same indentation:

```fortran

        block
            real(kind=RKIND), pointer :: tuvx_update_interval_ptr
            nullify(tuvx_update_interval_ptr)
            call mpas_pool_get_config(configs, 'config_tuvx_update_interval', tuvx_update_interval_ptr)
            if (associated(tuvx_update_interval_ptr)) tuvx_update_interval = tuvx_update_interval_ptr
        end block
```

The `block` pattern (rather than adding a top-of-routine pointer declaration) follows the existing precedent in this file at lines 133–155 (`config_chemistry_ref_solve`, `config_chem_substeps`, `config_micm_relative_tolerance`).

- [ ] **Step 3: Verify the additions**

Run:
```bash
grep -n 'tuvx_update_interval\|tuvx_time_since_last' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: 4 matches —
- Line ~42: `tuvx_update_interval` save declaration
- Line ~43: `tuvx_time_since_last` save declaration
- Line ~172 (or thereabouts): `mpas_pool_get_config` call for `config_tuvx_update_interval`
- Line ~173 (or thereabouts): `if (associated(tuvx_update_interval_ptr))` line

(Plus one more once Task 3 lands — but that's a future state. After Task 2, expect 4.)

---

### Task 3: Wrap the photolysis block in `chemistry_step` with the throttle gate

**Files:**
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (chemistry_step body, currently lines 383–493)

This task is one large structural edit. The `mpas_get_time` + `hour_utc` lines move inside a new `if/else` block, which also wraps the entire existing photolysis section through the trailing `deallocate`.

- [ ] **Step 1: Replace the photolysis section with the gated version**

Locate this exact block in `chemistry_step` (currently lines 383–493). Note: the block starts with the `! Compute the time-of-day inputs` comment and ends with the `deallocate(photo_rates, cos_sza_cell)` line that immediately precedes the `call mpas_log_write('Stepping chemistry packages...')` line:

```fortran
        ! Compute the time-of-day inputs once; per-cell SZA is filled below
        ! after the n_photo_rp check (so we can clean up cos_sza_cell on error
        ! paths cleanly alongside photo_rates).
        call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
        hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                 + real(second, RKIND) / 3600.0_RKIND

        ! Update photolysis rates. Unified multi-rate path: TUV-x fills
        ! photo_rates(ir,k,iCell) for every reaction it reports; when TUV-x
        ! is disabled, the fallback fills a single 'jNO2' slot with
        ! j_max * max(0, cos_sza).
        if (n_photo_rp < 1) then
            call mpas_log_write('[Chemistry] photolysis rate indices not cached.', &
                messageType=MPAS_LOG_CRIT)
            return
        end if

        ! Per-cell solar zenith angle. Filled per cell when the grid-coords
        ! switch is enabled, broadcast from the namelist scalars otherwise.
        allocate(cos_sza_cell(nCells))
        if (chem_use_grid_coords) then
            call mpas_pool_get_array(mesh, 'latCell', latCell)
            call mpas_pool_get_array(mesh, 'lonCell', lonCell)
            do iCell = 1, nCells
                cos_sza_cell(iCell) = solar_cos_sza(DoY, hour_utc, &
                                                    latCell(iCell) * RAD2DEG, &
                                                    lonCell(iCell) * RAD2DEG)
            end do
        else
            cos_sza_cell(:) = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)
        end if

        allocate(photo_rates(n_photo_rp, nVertLevels, nCells))
        photo_rates = 0.0_RKIND

        if (use_tuvx) then
            if (idx_qO3 < 1) then
                call mpas_log_write('[Chemistry] TUV-x enabled but index_qO3 is unresolved.', &
                    messageType=MPAS_LOG_CRIT)
                deallocate(photo_rates, cos_sza_cell)
                return
            end if

            allocate(rho_air_col(nVertLevels))
            allocate(temperature_col(nVertLevels))
            allocate(qc_col(nVertLevels))
            allocate(qr_col(nVertLevels))

            call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
            call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
            call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
            call mpas_pool_get_array(mesh, 'zz', zz)
            call mpas_pool_get_array(mesh, 'zgrid', zgrid)
            call mpas_pool_get_array(diag, 'exner', exner)
            call mpas_pool_get_dimension(state, 'index_qv', index_qv)

            do iCell = 1, nCells
                do k = 1, nVertLevels
                    rho_air_col(k) = zz(k, iCell) * rho_zz(k, iCell) * &
                                     (1.0_RKIND + scalars(index_qv, k, iCell))
                    temperature_col(k) = (theta_m(k, iCell) / &
                        (1.0_RKIND + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell)
                end do

                if (idx_qc > 0) then
                    qc_col(:) = scalars(idx_qc, :, iCell)
                else
                    qc_col(:) = 0.0_RKIND
                end if
                if (idx_qr > 0) then
                    qr_col(:) = scalars(idx_qr, :, iCell)
                else
                    qr_col(:) = 0.0_RKIND
                end if

                call tuvx_compute_photolysis(cos_sza_cell(iCell), zgrid(:, iCell), &
                                             rho_air_col, temperature_col, &
                                             scalars(idx_qO3, :, iCell), &
                                             qc_col, qr_col, &
                                             nVertLevels, &
                                             photo_rates(:, :, iCell), &
                                             error_code, error_message)
                if (error_code /= 0) then
                    call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                    deallocate(photo_rates, rho_air_col, temperature_col, qc_col, qr_col, cos_sza_cell)
                    return
                end if
            end do

            deallocate(rho_air_col, temperature_col, qc_col, qr_col)
        else
            ! Phase-1 fallback: single rate 'jNO2' = j_max * max(0, cos_sza).
            ! Assumes n_photo_rp == 1 and the cached rate name is 'jNO2' —
            ! musica_cache_photo_indices was called with that name at init.
            do iCell = 1, nCells
                photo_rates(1, :, iCell) = chem_j_no2_max * &
                                            max(0.0_RKIND, cos_sza_cell(iCell))
            end do
        end if

        call musica_set_photolysis_rates(photo_rates, nCells, nVertLevels, &
                                         n_photo_rp, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            deallocate(photo_rates, cos_sza_cell)
            return
        end if

        call chemistry_set_photolysis_diag(diag, photo_rates, nCells, nVertLevels)

        deallocate(photo_rates, cos_sza_cell)
```

REPLACE that entire block with this version (the `n_photo_rp` check stays before the gate; everything else moves inside; the closing `else` adds `dt` to the accumulator):

```fortran
        ! Update photolysis rates. Unified multi-rate path: TUV-x fills
        ! photo_rates(ir,k,iCell) for every reaction it reports; when TUV-x
        ! is disabled, the fallback fills a single 'jNO2' slot with
        ! j_max * max(0, cos_sza). The whole block is gated by
        ! tuvx_update_interval; on skip steps the accumulator just ticks
        ! forward and MICM keeps using the last-set rate parameters.
        if (n_photo_rp < 1) then
            call mpas_log_write('[Chemistry] photolysis rate indices not cached.', &
                messageType=MPAS_LOG_CRIT)
            return
        end if

        if (tuvx_time_since_last >= tuvx_update_interval) then
            ! Compute the time-of-day inputs (only needed when an update fires).
            call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
            hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                     + real(second, RKIND) / 3600.0_RKIND

            ! Per-cell solar zenith angle. Filled per cell when the grid-coords
            ! switch is enabled, broadcast from the namelist scalars otherwise.
            allocate(cos_sza_cell(nCells))
            if (chem_use_grid_coords) then
                call mpas_pool_get_array(mesh, 'latCell', latCell)
                call mpas_pool_get_array(mesh, 'lonCell', lonCell)
                do iCell = 1, nCells
                    cos_sza_cell(iCell) = solar_cos_sza(DoY, hour_utc, &
                                                        latCell(iCell) * RAD2DEG, &
                                                        lonCell(iCell) * RAD2DEG)
                end do
            else
                cos_sza_cell(:) = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)
            end if

            allocate(photo_rates(n_photo_rp, nVertLevels, nCells))
            photo_rates = 0.0_RKIND

            if (use_tuvx) then
                if (idx_qO3 < 1) then
                    call mpas_log_write('[Chemistry] TUV-x enabled but index_qO3 is unresolved.', &
                        messageType=MPAS_LOG_CRIT)
                    deallocate(photo_rates, cos_sza_cell)
                    return
                end if

                allocate(rho_air_col(nVertLevels))
                allocate(temperature_col(nVertLevels))
                allocate(qc_col(nVertLevels))
                allocate(qr_col(nVertLevels))

                call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
                call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
                call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
                call mpas_pool_get_array(mesh, 'zz', zz)
                call mpas_pool_get_array(mesh, 'zgrid', zgrid)
                call mpas_pool_get_array(diag, 'exner', exner)
                call mpas_pool_get_dimension(state, 'index_qv', index_qv)

                do iCell = 1, nCells
                    do k = 1, nVertLevels
                        rho_air_col(k) = zz(k, iCell) * rho_zz(k, iCell) * &
                                         (1.0_RKIND + scalars(index_qv, k, iCell))
                        temperature_col(k) = (theta_m(k, iCell) / &
                            (1.0_RKIND + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell)
                    end do

                    if (idx_qc > 0) then
                        qc_col(:) = scalars(idx_qc, :, iCell)
                    else
                        qc_col(:) = 0.0_RKIND
                    end if
                    if (idx_qr > 0) then
                        qr_col(:) = scalars(idx_qr, :, iCell)
                    else
                        qr_col(:) = 0.0_RKIND
                    end if

                    call tuvx_compute_photolysis(cos_sza_cell(iCell), zgrid(:, iCell), &
                                                 rho_air_col, temperature_col, &
                                                 scalars(idx_qO3, :, iCell), &
                                                 qc_col, qr_col, &
                                                 nVertLevels, &
                                                 photo_rates(:, :, iCell), &
                                                 error_code, error_message)
                    if (error_code /= 0) then
                        call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                        deallocate(photo_rates, rho_air_col, temperature_col, qc_col, qr_col, cos_sza_cell)
                        return
                    end if
                end do

                deallocate(rho_air_col, temperature_col, qc_col, qr_col)
            else
                ! Phase-1 fallback: single rate 'jNO2' = j_max * max(0, cos_sza).
                ! Assumes n_photo_rp == 1 and the cached rate name is 'jNO2' —
                ! musica_cache_photo_indices was called with that name at init.
                do iCell = 1, nCells
                    photo_rates(1, :, iCell) = chem_j_no2_max * &
                                                max(0.0_RKIND, cos_sza_cell(iCell))
                end do
            end if

            call musica_set_photolysis_rates(photo_rates, nCells, nVertLevels, &
                                             n_photo_rp, error_code, error_message)
            if (error_code /= 0) then
                call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                deallocate(photo_rates, cos_sza_cell)
                return
            end if

            call chemistry_set_photolysis_diag(diag, photo_rates, nCells, nVertLevels)

            deallocate(photo_rates, cos_sza_cell)

            tuvx_time_since_last = 0.0_RKIND
        else
            tuvx_time_since_last = tuvx_time_since_last + dt
        end if
```

The structural changes versus the original:

1. The first comment block ("Compute the time-of-day inputs once...") and the `mpas_get_time` + `hour_utc` lines are **removed** from the top of the section and **re-added** inside the new `if (tuvx_time_since_last >= tuvx_update_interval) then` branch (with a slightly updated leading comment).
2. The leading comment for the photolysis block is updated to mention the throttle.
3. A new `if (tuvx_time_since_last >= tuvx_update_interval) then` wraps everything from the `cos_sza_cell` allocation through the existing `deallocate(photo_rates, cos_sza_cell)` line.
4. Inside the gate, after the existing `deallocate`, a new line `tuvx_time_since_last = 0.0_RKIND` resets the accumulator.
5. A new `else` branch contains a single line `tuvx_time_since_last = tuvx_time_since_last + dt`.

Indentation of the gated body increased by one level (4 spaces). All early-return paths (`return` after MPAS_LOG_CRIT) inside the gated body still work — they exit `chemistry_step` entirely without touching the accumulator, leaving it `>= tuvx_update_interval` so the next call will re-attempt the update.

- [ ] **Step 2: Verify the gate structure**

Run:
```bash
grep -n 'tuvx_time_since_last\|tuvx_update_interval' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: 7 matches —
- Two save declarations (Task 2 Step 1)
- One `mpas_pool_get_config` call (Task 2 Step 2)
- One `if (associated(...))` line (Task 2 Step 2)
- One `if (tuvx_time_since_last >= tuvx_update_interval)` gate condition (Task 3)
- One `tuvx_time_since_last = 0.0_RKIND` reset (Task 3)
- One `tuvx_time_since_last = tuvx_time_since_last + dt` accumulate (Task 3)

Run:
```bash
grep -c 'mpas_get_time(currTime' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: `1` (the call moved into the gate; not duplicated).

Run:
```bash
grep -c 'allocate(cos_sza_cell' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
grep -c 'allocate(photo_rates' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: `1` and `1` (allocations not duplicated by the structural change).

---

### Task 4: Build verification

**Files:** none modified.

- [ ] **Step 1: Clean the chemistry artifacts**

Run from the repo root:
```bash
rm -f src/core_atmosphere/chemistry/mpas_atm_chemistry.o
rm -f src/core_atmosphere/chemistry/mpas_atm_chemistry.mod
```

- [ ] **Step 2: Build atmosphere with MUSICA**

Use the platform-appropriate command (per `BUILD.md`). On macOS:

```bash
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
  CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true 2>&1 | tail -40
```

On Ubuntu (with `conda activate mpas`):

```bash
eval "$(scripts/check_build_env.sh --export)" && make -j8 gfortran \
  CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true 2>&1 | tail -40
```

Expected: build completes, `atmosphere_model` produced. No error messages mentioning `tuvx_update_interval`, `tuvx_time_since_last`, `cos_sza_cell`, or `photo_rates`.

If the build stops in `physics_mmm` with a network-fetch error, that's an unrelated upstream issue documented in `BUILD.md`. Confirm chemistry compiled successfully:

```bash
test -f src/core_atmosphere/chemistry/mpas_atm_chemistry.o && echo "chemistry: ok"
```
Expected: `chemistry: ok`.

- [ ] **Step 3: Confirm the executable links**

Run: `ls -l atmosphere_model`
Expected: file exists, recently modified, ~15 MB.

---

### Task 5: Update TUV-x integration doc

**Files:**
- Modify: `docs/guides/TUVX_INTEGRATION.md` (add subsection in "Multi-Photolysis Plumbing")

- [ ] **Step 1: Locate the existing `### Per-cell solar geometry` subsection**

Run: `grep -n '### Per-cell solar geometry' docs/guides/TUVX_INTEGRATION.md`
Expected: one line number (currently around line 263).

- [ ] **Step 2: Add a new `### TUV-x update interval` subsection IMMEDIATELY AFTER the existing per-cell subsection**

Find the `### Per-cell solar geometry` subsection and its closing paragraph. The full subsection currently looks like:

```
### Per-cell solar geometry

By default chemistry uses the namelist scalars
`config_chemistry_latitude` / `config_chemistry_longitude` to compute a
single `cos_sza` shared by every column — appropriate for idealized
Cartesian-plane test cases (supercell, mountain wave, baroclinic wave).
For real spherical-mesh runs, set `config_chemistry_use_grid_coords =
.true.` in the `&musica` namelist; chemistry will then read `latCell`
and `lonCell` from the mesh and compute `cos_sza` per cell. The
per-cell path is used by both the TUV-x photolysis call and the
fallback `cos(SZA)` `j_NO2` rate. Default is `.false.` to preserve
exact bit-reproducibility of all idealized cases.
```

Insert this new subsection IMMEDIATELY AFTER it (with one blank line separator, before the next existing `## Column Extension Above MPAS Top` header):

```
### TUV-x update interval

By default chemistry calls TUV-x every dynamics step. For long
integrations or fine `dt` values, this is more often than necessary —
photolysis rates evolve on minute-to-hour timescales, while the
dynamics `dt` may be a few seconds. Set
`config_tuvx_update_interval` (real, simulated seconds) in the
`&musica` namelist to throttle TUV-x to a coarser cadence. When the
configured interval has not yet elapsed since the last update,
`chemistry_step` skips the entire photolysis block; MICM keeps using
the rate parameters last set, and the `j_<rxn>` diagnostic fields
keep their last-written values. Default is `0.0` (update every step)
to preserve exact bit-reproducibility of all existing runs.

Typical operational values: 60 s for short-duration / fine-`dt`
idealized cases; several minutes for global production runs. A 60 s
interval at `dt = 3 s` runs TUV-x once every 20 chemistry steps.
```

- [ ] **Step 3: Verify the doc edit**

Run: `grep -c '^### ' docs/guides/TUVX_INTEGRATION.md`
Expected: count up by exactly 1 from the pre-edit value (one new `### TUV-x update interval` subsection).

Run: `grep -n 'config_tuvx_update_interval' docs/guides/TUVX_INTEGRATION.md`
Expected: one match (the new namelist option in the new paragraph).

---

### Task 6: Commit code + doc together

**Files:** stages all of the changes from Tasks 1, 2, 3, and 5.

- [ ] **Step 1: Inspect what will be committed**

Run:
```bash
git status
git diff --stat src/core_atmosphere/Registry.xml src/core_atmosphere/chemistry/mpas_atm_chemistry.F docs/guides/TUVX_INTEGRATION.md
```

Expected: three files modified — `src/core_atmosphere/Registry.xml`, `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`, and `docs/guides/TUVX_INTEGRATION.md`. No untracked files in `src/`.

- [ ] **Step 2: Stage and commit**

Run:
```bash
git add \
  src/core_atmosphere/Registry.xml \
  src/core_atmosphere/chemistry/mpas_atm_chemistry.F \
  docs/guides/TUVX_INTEGRATION.md
git commit -m "$(cat <<'EOF'
feat(chemistry): throttle TUV-x to a configurable update interval

Adds config_tuvx_update_interval to the &musica namelist record
(real, seconds, default 0.0). When 0.0, behavior is bit-identical to
the per-cell-lat/lon baseline at d7601e4 — TUV-x runs every chemistry
step. When positive, the entire photolysis block in chemistry_step
(time-of-day computation, cos_sza_cell allocation/fill, photo_rates
allocation, per-cell TUV-x or fallback loop, musica_set_photolysis_rates,
chemistry_set_photolysis_diag, deallocations) is gated by an
accumulator. On skip steps, MICM keeps using the rate parameters
last set, and the j_<rxn> diag fields keep their last-written values.

The accumulator is a module-scope save initialized to huge() so the
first chemistry_step call always fires an update — also covers the
restart case, since save re-init in a fresh process matches MICM's
zeroed rate parameters.

Also updates docs/guides/TUVX_INTEGRATION.md with a new
"TUV-x update interval" subsection in the "Multi-Photolysis Plumbing"
section.

Spec: docs/superpowers/specs/2026-04-19-tuvx-update-interval-design.md
Plan: docs/superpowers/plans/2026-04-19-tuvx-update-interval.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 3: Verify the commit**

Run: `git log -1 --stat`
Expected: commit subject `feat(chemistry): throttle TUV-x to a configurable update interval`; three files in the stat line.

---

### Task 7: Reproducibility validation (`interval = 0.0` on supercell)

**Files:** none modified — runtime check.

This task verifies that the default (`config_tuvx_update_interval = 0.0`) preserves bit-reproducibility of the existing supercell run. Compares against the reference `output.nc` produced after the lat/lon work landed at commit `d7601e4`.

- [ ] **Step 1: Stage the new build into the supercell run directory**

Run:
```bash
cp atmosphere_model ~/Data/CheMPAS/supercell/
```

- [ ] **Step 2: Confirm the namelist does NOT set `config_tuvx_update_interval`**

Run: `grep 'config_tuvx_update_interval' ~/Data/CheMPAS/supercell/namelist.atmosphere`
Expected: no match (option absent → defaults to `0.0`).

If the option is present in the namelist, remove it (or set it to `0.0`) for this validation; the reproducibility check requires the default behavior path.

- [ ] **Step 3: Run the supercell case (3-min sim is sufficient)**

If a `config_run_duration` of `00:15:00` is in the namelist, temporarily reduce it to `00:03:00` for this check (and restore afterward). 3 minutes of simulated time covers the first output write (output_interval `00:03:00`) and is plenty for spot-check.

Run:
```bash
cd ~/Data/CheMPAS/supercell/
rm -f log.atmosphere.* output.nc
mpiexec -n 8 ./atmosphere_model 2>&1 | tail -10
```

Expected: clean run completion (no `CRITICAL ERROR` lines in `log.atmosphere.0000.out`).

- [ ] **Step 4: Spot-check key values**

Run:
```bash
~/miniconda3/envs/mpas/bin/python -c "
import netCDF4 as nc
import numpy as np
f = nc.Dataset('/Users/fillmore/Data/CheMPAS/supercell/output.nc')
for var in ['j_jNO2', 'qNO', 'qNO2', 'qO3']:
    if var in f.variables:
        a = f.variables[var][-1]
        print(f'{var:8s} min={np.min(a):.4e}  max={np.max(a):.4e}')
"
```

Expected (from the reference recorded after commit `d7601e4` at `t=3 min`):
- `j_jNO2`: min ≈ `5.671e-05`, max ≈ `1.704e-02`
- `qNO`: min ≈ `0.000e+00`, max ≈ `1.115e-07`
- `qNO2`: min ≈ `0.000e+00`, max ≈ `2.985e-08`
- `qO3`: min ≈ `5.162e-08`, max ≈ `8.276e-08`

All values must match within at most a few units in the last decimal place (floating-point identical, since the broadcast path through both the lat/lon switch and the new throttle should be mathematically identical to the original scalar path).

- [ ] **Step 5: Restore the namelist run duration if changed**

If you reduced `config_run_duration` for Step 3, restore it now:
```bash
sed -i.bak "s/config_run_duration = '00:03:00'/config_run_duration = '00:15:00'/" ~/Data/CheMPAS/supercell/namelist.atmosphere
rm -f ~/Data/CheMPAS/supercell/namelist.atmosphere.bak
```

- [ ] **Step 6: Report**

If values match the reference within the stated tolerance: report DONE.
If values diverge significantly: report DONE_WITH_CONCERNS, paste the actual values, and stop. Do not proceed to Task 8 — investigate the discrepancy first (likely cause: bug in the gate condition or accumulator logic).

---

### Task 8: Throttle validation (`interval = 60.0` on supercell)

**Files:** Stages a one-off namelist override in the run directory and a temporary log line in source for verification.

This task confirms the throttle actually skips work and the cached photolysis rates are reused on skip steps.

- [ ] **Step 1: Add temporary instrumentation to count update fires**

Edit `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`. Find the new line (added in Task 3):

```fortran
            tuvx_time_since_last = 0.0_RKIND
```

Add a `mpas_log_write` immediately BEFORE it:

```fortran
            call mpas_log_write('[TUV-x] photolysis update fired')
            tuvx_time_since_last = 0.0_RKIND
```

This log line will be printed once per update fire. It will be removed in Step 6.

- [ ] **Step 2: Rebuild**

Run from the repo root:
```bash
rm -f src/core_atmosphere/chemistry/mpas_atm_chemistry.o src/core_atmosphere/chemistry/mpas_atm_chemistry.mod
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
  CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true 2>&1 | tail -10
```
(Use `gfortran` instead of `llvm` on Ubuntu.)

Expected: clean build.

```bash
cp atmosphere_model ~/Data/CheMPAS/supercell/
```

- [ ] **Step 3: Set the throttle interval and run for 3 simulated minutes**

Edit `~/Data/CheMPAS/supercell/namelist.atmosphere` and add (or set, if already present) inside the `&musica` block:

```
config_tuvx_update_interval = 60.0
```

Confirm `config_run_duration = '00:03:00'` is set in the `&nhyd_model` block (temporarily reduce from `00:15:00` if needed; the throttle test only needs 3 simulated minutes).

Run:
```bash
cd ~/Data/CheMPAS/supercell/
rm -f log.atmosphere.* output.nc
time mpiexec -n 8 ./atmosphere_model 2>&1 | tail -10
```

Expected: clean run completion, noticeably lower wall time than the unthrottled Task 7 run (per-step chemistry cost dominated by TUV-x; throttling by 20× reduces total chemistry cost roughly proportionally).

- [ ] **Step 4: Verify the update count**

Run:
```bash
grep -c '\[TUV-x\] photolysis update fired' ~/Data/CheMPAS/supercell/log.atmosphere.0000.out
```

Expected: `4` exactly. The four firings should be at simulated times t=0 (initial, accumulator is `huge`), t=60 s, t=120 s, and t=180 s (run end at t=180 s).

If the count is significantly different (e.g., 60 like the unthrottled case, or 0 = never): the throttle is not working as designed. Stop and investigate.

- [ ] **Step 5: Verify cached `j_jNO2` is physically reasonable**

Run:
```bash
~/miniconda3/envs/mpas/bin/python -c "
import netCDF4 as nc
import numpy as np
f = nc.Dataset('/Users/fillmore/Data/CheMPAS/supercell/output.nc')
a = f.variables['j_jNO2'][-1]
print(f'j_jNO2 throttled: min={np.min(a):.4e}  max={np.max(a):.4e}')
"
```

Expected: min and max within ~1% of the unthrottled values from Task 7 Step 4. (Slow SZA evolution over 60 s means the cached value is sub-percent stale at most times of day.)

- [ ] **Step 6: Remove temporary instrumentation, restore namelist, and rebuild**

Edit `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` and REMOVE the temporary line added in Step 1:

```fortran
            call mpas_log_write('[TUV-x] photolysis update fired')
```

(Keep the `tuvx_time_since_last = 0.0_RKIND` line that follows it.)

Verify removal:
```bash
grep -c 'photolysis update fired' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: `0`.

Edit `~/Data/CheMPAS/supercell/namelist.atmosphere` to remove the `config_tuvx_update_interval = 60.0` line (return to default). Restore `config_run_duration = '00:15:00'` if you changed it.

Rebuild to confirm the source is back to the committed state and still compiles:
```bash
rm -f src/core_atmosphere/chemistry/mpas_atm_chemistry.o src/core_atmosphere/chemistry/mpas_atm_chemistry.mod
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
  CORE=atmosphere PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" \
  PRECISION=double MUSICA=true 2>&1 | tail -5
```

Expected: clean build.

Run: `git diff src/core_atmosphere/chemistry/mpas_atm_chemistry.F`
Expected: no diff (the source is back to the committed state from Task 6).

- [ ] **Step 7: Report**

If Steps 4 and 5 pass: report DONE. The throttle is validated.

---

## Self-Review Notes

**Spec coverage:**

- Spec § *Namelist & Registry* → Task 1.
- Spec § *Module-scope state* → Task 2 Step 1.
- Spec § *`chemistry_init` — read the new config* → Task 2 Step 2.
- Spec § *`chemistry_step` — gate the photolysis block* → Task 3.
- Spec § *Correctness analysis* → preserved by Task 3's structure (default `0.0` + `huge` initial covered; first call always updates; restart correctness; early-return paths don't corrupt accumulator; accumulator only ticks on skip steps).
- Spec § *Validation* item 1 (default behavior preserved) → Task 7.
- Spec § *Validation* item 2 (throttle skips work, update cadence) → Task 8.
- Spec § *Validation* item 3 (physical sanity) → Task 8 Step 5.
- Spec § *Files Touched* → Tasks 1, 3, 5 cover exactly the listed files.
- Spec § *Doc update* → Task 5.
- Spec § *Commit Shape* → Task 6.

**Placeholder scan:** No TBD/TODO/vague-instruction patterns. Every step shows the exact code or command and an expected outcome. The temporary instrumentation in Task 8 is explicitly added and explicitly removed in the same task.

**Type / name consistency:** New names — `config_tuvx_update_interval` (namelist), `tuvx_update_interval` (module save), `tuvx_time_since_last` (module save), `tuvx_update_interval_ptr` (init local pointer). Each appears in the same form in every task that references it.

**Note on line numbers:** Line numbers throughout the plan are anchored to the current state of `develop` at planning time (after commit `a6b029e`). If new commits land between planning and execution, prefer the surrounding-context anchors (e.g., "the existing `chem_use_grid_coords` block") over the literal line numbers.
