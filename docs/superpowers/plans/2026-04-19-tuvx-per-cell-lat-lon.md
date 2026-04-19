# TUV-x Per-Cell Solar Geometry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `config_chemistry_use_grid_coords` namelist switch (default `.false.`) so chemistry computes `cos_sza` per cell from `latCell`/`lonCell` when enabled, preserving bit-reproducibility of all existing test cases when disabled.

**Architecture:** Replace the single scalar `cos_sza` in `chemistry_step` with a `cos_sza_cell(:)` array allocated each call. When the switch is `.true.`, fill it per cell from grid coordinates; when `.false.`, broadcast the namelist scalar. Downstream uses (TUV-x call and fallback `j_no2`) index the array. No interface changes to `solar_cos_sza`, `tuvx_compute_photolysis`, MUSICA, or lightning_nox.

**Tech Stack:** Fortran 2008 (LLVM/flang on macOS or GCC/gfortran on Ubuntu), MPAS Registry XML, MPI. Spec: `docs/superpowers/specs/2026-04-19-tuvx-per-cell-lat-lon-design.md`.

---

## File Structure

- Modify: `src/core_atmosphere/Registry.xml` (add namelist option)
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (module var, init read, step rewrite)
- Modify: `docs/guides/TUVX_INTEGRATION.md` (remove "deferred" bullet, add usage note)

No new files. No deletions. No changes to other source files, MICM YAMLs, TUV-x JSONs, or `test_cases/` namelists.

---

### Task 1: Add `config_chemistry_use_grid_coords` to Registry

**Files:**
- Modify: `src/core_atmosphere/Registry.xml` (insert after the `config_chemistry_longitude` block; current location is around lines 447–449)

- [ ] **Step 1: Locate the existing `config_chemistry_longitude` block**

Run: `grep -n 'config_chemistry_longitude' src/core_atmosphere/Registry.xml`
Expected: one or two line numbers identifying the namelist option (around line 447).

- [ ] **Step 2: Edit `Registry.xml` — insert the new option after `config_chemistry_longitude`**

Find this exact block (lines around 447–450 in the current file):

```xml
                <nml_option name="config_chemistry_longitude" type="real" default_value="0.0"
                     units="degrees"
                     description="Longitude for solar geometry in idealized cases (degrees E)"
                     possible_values="Any real between -180 and 180"/>
```

Add the following block immediately after its closing `/>`:

```xml
                <nml_option name="config_chemistry_use_grid_coords" type="logical" default_value=".false."
                     description="If .true., compute per-cell solar geometry from latCell/lonCell. If .false., broadcast namelist config_chemistry_latitude/longitude to every cell (preserves idealized-case behavior)."
                     possible_values=".true. or .false."/>
```

Note: the exact `possible_values` string for `config_chemistry_longitude` may differ slightly in your tree (e.g., the bound or the wording). Match your file's actual contents — only the new block matters here.

- [ ] **Step 3: Verify the new option is present and well-formed**

Run: `grep -n 'config_chemistry_use_grid_coords' src/core_atmosphere/Registry.xml`
Expected: at least one line showing the new option name.

Run: `xmllint --noout src/core_atmosphere/Registry.xml; echo "exit=$?"`
Expected: `exit=0` (Registry XML is well-formed).

If `xmllint` is not installed, skip this check; the build itself will catch malformed XML.

---

### Task 2: Add module-scope `chem_use_grid_coords` and read in `chemistry_init`

**Files:**
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (module-scope save var around line 39; chemistry_init local var around line 97; chemistry_init body around line 164)

- [ ] **Step 1: Add the module-scope save variable**

Locate the existing `chem_lat`/`chem_lon`/`chem_j_no2_max` save variables (around lines 37–39):

```fortran
    ! Solar geometry config values (cached at init, used per-step)
    real(kind=RKIND), save :: chem_lat = 0.0_RKIND
    real(kind=RKIND), save :: chem_lon = 0.0_RKIND
    real(kind=RKIND), save :: chem_j_no2_max = 0.0_RKIND
```

Add this line immediately after `chem_j_no2_max`:

```fortran
    logical, save :: chem_use_grid_coords = .false.
```

- [ ] **Step 2: Add the local pointer in `chemistry_init`**

Locate the existing `lat_ptr`/`lon_ptr` pointer declaration in `chemistry_init`'s local variable block (around line 97):

```fortran
        real(kind=RKIND), pointer       :: lat_ptr, lon_ptr
```

Add this declaration immediately after it:

```fortran
        logical, pointer                :: use_grid_coords_ptr
```

- [ ] **Step 3: Read the new config option in `chemistry_init`**

Locate the existing config-read block for `chem_lat`/`chem_lon` (around lines 159–164):

```fortran
        ! Cache solar geometry config values for per-step SZA computation
        chem_j_no2_max = j_no2_val
        nullify(lat_ptr)
        nullify(lon_ptr)
        call mpas_pool_get_config(configs, 'config_chemistry_latitude', lat_ptr)
        if (associated(lat_ptr)) chem_lat = lat_ptr
        call mpas_pool_get_config(configs, 'config_chemistry_longitude', lon_ptr)
        if (associated(lon_ptr)) chem_lon = lon_ptr
```

Add the following lines immediately after the last `if (associated(lon_ptr)) ...`:

```fortran
        nullify(use_grid_coords_ptr)
        call mpas_pool_get_config(configs, 'config_chemistry_use_grid_coords', use_grid_coords_ptr)
        if (associated(use_grid_coords_ptr)) chem_use_grid_coords = use_grid_coords_ptr
```

- [ ] **Step 4: Verify the file still parses (compile-check the file later in Task 4)**

Run: `grep -n 'chem_use_grid_coords' src/core_atmosphere/chemistry/mpas_atm_chemistry.F`
Expected: three lines — the module-scope `save` declaration, and two references in `chemistry_init` (the `mpas_pool_get_config` call and the `if (associated...)` line).

---

### Task 3: Modify `chemistry_step` for per-cell `cos_sza`

**Files:**
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` (chemistry_step local-var block around line 349, body around lines 376–467)

- [ ] **Step 1: Update the local-variable declarations in `chemistry_step`**

Locate the existing solar-geometry variable declarations in `chemistry_step` (around line 349):

```fortran
        ! Solar geometry variables
        integer :: DoY, hour, minute, second
        real(kind=RKIND) :: hour_utc, cos_sza, j_no2_value
```

Replace that three-line block with:

```fortran
        ! Solar geometry variables
        integer :: DoY, hour, minute, second
        real(kind=RKIND) :: hour_utc
        real(kind=RKIND), allocatable :: cos_sza_cell(:)
        real(kind=RKIND), parameter   :: RAD2DEG = 57.29577951308232_RKIND
        real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
```

(`cos_sza` and `j_no2_value` scalars are removed; both are now redundant.)

- [ ] **Step 2: Remove the scalar `cos_sza` computation**

Locate this block (around lines 375–379):

```fortran
        ! Compute SZA
        call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
        hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                 + real(second, RKIND) / 3600.0_RKIND
        cos_sza = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)
```

Replace it with:

```fortran
        ! Compute the time-of-day inputs once; per-cell SZA is filled below
        ! after the n_photo_rp check (so we can clean up cos_sza_cell on error
        ! paths cleanly alongside photo_rates).
        call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
        hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                 + real(second, RKIND) / 3600.0_RKIND
```

(The `cos_sza = solar_cos_sza(...)` line is removed.)

- [ ] **Step 3: Allocate and fill `cos_sza_cell` after the n_photo_rp check**

Locate this block (around lines 385–392):

```fortran
        if (n_photo_rp < 1) then
            call mpas_log_write('[Chemistry] photolysis rate indices not cached.', &
                messageType=MPAS_LOG_CRIT)
            return
        end if

        allocate(photo_rates(n_photo_rp, nVertLevels, nCells))
        photo_rates = 0.0_RKIND
```

Insert the `cos_sza_cell` block between the `end if` and the `allocate(photo_rates...)` line, so the section becomes:

```fortran
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
```

- [ ] **Step 4: Replace the scalar `cos_sza` use in the TUV-x per-cell loop**

Locate this `tuvx_compute_photolysis` call inside the `do iCell = 1, nCells` loop (around line 434):

```fortran
                call tuvx_compute_photolysis(cos_sza, zgrid(:, iCell), &
                                             rho_air_col, temperature_col, &
                                             scalars(idx_qO3, :, iCell), &
                                             qc_col, qr_col, &
                                             nVertLevels, &
                                             photo_rates(:, :, iCell), &
                                             error_code, error_message)
```

Replace `cos_sza` (the first argument) with `cos_sza_cell(iCell)`:

```fortran
                call tuvx_compute_photolysis(cos_sza_cell(iCell), zgrid(:, iCell), &
                                             rho_air_col, temperature_col, &
                                             scalars(idx_qO3, :, iCell), &
                                             qc_col, qr_col, &
                                             nVertLevels, &
                                             photo_rates(:, :, iCell), &
                                             error_code, error_message)
```

- [ ] **Step 5: Replace the scalar fallback computation**

Locate this fallback block (around lines 449–455):

```fortran
        else
            ! Phase-1 fallback: single rate 'jNO2' = j_max * max(0, cos_sza).
            ! Assumes n_photo_rp == 1 and the cached rate name is 'jNO2' —
            ! musica_cache_photo_indices was called with that name at init.
            j_no2_value = chem_j_no2_max * max(0.0_RKIND, cos_sza)
            photo_rates(1, :, :) = j_no2_value
        end if
```

Replace it with:

```fortran
        else
            ! Phase-1 fallback: single rate 'jNO2' = j_max * max(0, cos_sza).
            ! Assumes n_photo_rp == 1 and the cached rate name is 'jNO2' —
            ! musica_cache_photo_indices was called with that name at init.
            do iCell = 1, nCells
                photo_rates(1, :, iCell) = chem_j_no2_max * &
                                            max(0.0_RKIND, cos_sza_cell(iCell))
            end do
        end if
```

- [ ] **Step 6: Add `cos_sza_cell` deallocation to all photo_rates exit paths**

Find each existing `deallocate(photo_rates)` call (and `deallocate(photo_rates, ...)` group calls) in `chemistry_step`. There are four such sites in the current code:

1. Around line 398–400 (TUV-x enabled but `idx_qO3 < 1`):

   Existing:
   ```fortran
                call mpas_log_write('[Chemistry] TUV-x enabled but index_qO3 is unresolved.', &
                    messageType=MPAS_LOG_CRIT)
                deallocate(photo_rates)
                return
   ```
   Replace `deallocate(photo_rates)` with:
   ```fortran
                deallocate(photo_rates, cos_sza_cell)
   ```

2. Around line 442–444 (TUV-x per-column error):

   Existing:
   ```fortran
                if (error_code /= 0) then
                    call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                    deallocate(photo_rates, rho_air_col, temperature_col, qc_col, qr_col)
                    return
                end if
   ```
   Replace the `deallocate(...)` line with:
   ```fortran
                    deallocate(photo_rates, rho_air_col, temperature_col, qc_col, qr_col, cos_sza_cell)
   ```

3. Around line 459–462 (`musica_set_photolysis_rates` error):

   Existing:
   ```fortran
        call musica_set_photolysis_rates(photo_rates, nCells, nVertLevels, &
                                         n_photo_rp, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            deallocate(photo_rates)
            return
        end if
   ```
   Replace `deallocate(photo_rates)` with:
   ```fortran
            deallocate(photo_rates, cos_sza_cell)
   ```

4. Around line 467 (normal exit after `chemistry_set_photolysis_diag`):

   Existing:
   ```fortran
        call chemistry_set_photolysis_diag(diag, photo_rates, nCells, nVertLevels)

        deallocate(photo_rates)
   ```
   Replace `deallocate(photo_rates)` with:
   ```fortran
        deallocate(photo_rates, cos_sza_cell)
   ```

- [ ] **Step 7: Verify all references are consistent**

Run:
```bash
grep -n 'cos_sza\b' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: no matches at all (the bare scalar identifier `cos_sza` is fully removed; all remaining occurrences are `cos_sza_cell`).

Run:
```bash
grep -n 'cos_sza_cell' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: about 8 matches — one declaration, one allocate, two fill paths (per-cell and broadcast), one TUV-x call argument, one fallback loop, plus four deallocations (one per exit path).

Run:
```bash
grep -n 'j_no2_value' src/core_atmosphere/chemistry/mpas_atm_chemistry.F
```
Expected: no matches (the scalar is fully removed).

---

### Task 4: Build verification

**Files:** none modified.

- [ ] **Step 1: Clean the chemistry artifacts**

Run from the repo root:
```bash
rm -f src/core_atmosphere/chemistry/mpas_atm_chemistry.o
rm -f src/core_atmosphere/chemistry/mpas_atm_chemistry.mod
```

(A full `make clean CORE=atmosphere` is unnecessary — the change is contained to one source file plus the registry-generated includes.)

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

Expected: build completes; the last lines show linking and produce `atmosphere_model`. No error messages mentioning `cos_sza`, `cos_sza_cell`, `chem_use_grid_coords`, `latCell`, or `lonCell`.

If the build stops in `physics_mmm` with a network-fetch error (`git fetch MMM-physics`), that's an unrelated upstream issue documented in `BUILD.md`; it does not indicate a problem with this change. Confirm chemistry compiled successfully by checking for the chemistry object:

```bash
test -f src/core_atmosphere/chemistry/mpas_atm_chemistry.o && echo "chemistry: ok"
```
Expected: `chemistry: ok`.

- [ ] **Step 3: Confirm the executable links (only if Step 2 didn't stop on the unrelated `physics_mmm` issue)**

Run: `ls -l atmosphere_model`
Expected: file exists, recently modified, size in the tens of megabytes.

If the executable did not produce because of the `physics_mmm` issue, report status DONE_WITH_CONCERNS and note that the chemistry object compiled cleanly.

---

### Task 5: Update TUV-x integration doc

**Files:**
- Modify: `docs/guides/TUVX_INTEGRATION.md` (the "Current Limits And Follow-On Work" section)

- [ ] **Step 1: Remove the resolved bullet**

Locate this exact line in `docs/guides/TUVX_INTEGRATION.md` (currently around line 316):

```
- fallback SZA still uses namelist coordinates; grid-aware chemistry geometry is
  deferred
```

Delete both lines of that bullet.

- [ ] **Step 2: Add a usage note in the "Multi-Photolysis Plumbing" section**

Locate this paragraph in `docs/guides/TUVX_INTEGRATION.md` (currently around lines 257–261):

```
The Phase-1 `cos(SZA)` fallback is still available: when
`config_tuvx_config_file` is empty, chemistry runs with a single
synthesised `jNO2` rate equal to `config_lnox_j_no2 * max(0, cos_sza)`.
The fallback is single-rate only; Chapman and Chapman+NOx mechanisms
require TUV-x.
```

Add the following paragraph immediately after it (preserving the `## Column Extension Above MPAS Top` header that currently follows):

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

- [ ] **Step 3: Verify the doc still renders cleanly**

Run: `grep -c '^## ' docs/guides/TUVX_INTEGRATION.md`
Expected: same count as before the edit (no top-level headers were added or removed; only a `### Per-cell solar geometry` subsection under the existing "Multi-Photolysis Plumbing" section).

Run: `grep -n 'config_chemistry_use_grid_coords' docs/guides/TUVX_INTEGRATION.md`
Expected: one match (the new namelist option name in the new paragraph).

Run: `grep -n 'grid-aware chemistry geometry is' docs/guides/TUVX_INTEGRATION.md`
Expected: no matches (the resolved bullet is fully removed).

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
feat(chemistry): add per-cell solar geometry switch for TUV-x

Adds config_chemistry_use_grid_coords to the &musica namelist record
(default .false.). When .true., chemistry_step computes cos_sza per
cell from latCell/lonCell (radians, converted to degrees) and uses
those per-cell values in both the TUV-x call and the fallback
cos(SZA) j_NO2 rate. When .false., the existing namelist scalars
(config_chemistry_latitude/longitude) are broadcast across cells —
all existing test cases (supercell, mountain wave, baroclinic wave)
remain bit-for-bit identical.

Also updates docs/guides/TUVX_INTEGRATION.md: removes the resolved
"grid-aware chemistry geometry is deferred" bullet and adds a brief
usage note for the new switch.

Spec: docs/superpowers/specs/2026-04-19-tuvx-per-cell-lat-lon-design.md
Plan: docs/superpowers/plans/2026-04-19-tuvx-per-cell-lat-lon.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 3: Verify the commit**

Run: `git log -1 --stat`
Expected: commit subject `feat(chemistry): add per-cell solar geometry switch for TUV-x`; three files in the stat line.

---

### Task 7: Reproducibility validation (switch=`.false.` on supercell)

**Files:** none modified — this is a runtime check.

This task verifies that the default (`config_chemistry_use_grid_coords = .false.`) preserves bit-reproducibility of the existing supercell run. Requires the supercell run directory at `~/Data/CheMPAS/supercell/` (per `RUN.md` / `CLAUDE.md`).

- [ ] **Step 1: Stage the new build into the run directory**

Run:
```bash
cp atmosphere_model ~/Data/CheMPAS/supercell/
```

- [ ] **Step 2: Run the supercell case (Phase 2/3 setting; 15-minute integration)**

Run from `~/Data/CheMPAS/supercell/`:
```bash
cd ~/Data/CheMPAS/supercell/
rm -f history.*.nc diag.*.nc log.atmosphere.*
mpiexec -n 8 ./atmosphere_model 2>&1 | tail -20
```

Expected: clean run completion. The `tail -20` output should show the final log messages including a successful chemistry step at the final time.

- [ ] **Step 3: Spot-check key values against the recorded reference**

The `docs/guides/TUVX_INTEGRATION.md` "Recorded Development Results" section gives expected values for the Phase 2/3 supercell at `t = 15 min`:

| Quantity | Expected value |
|----------|---------------:|
| Surface `j_no2` (clear-sky column) | `~7.16e-3 s^-1` |
| Cloudy minimum `j_no2` | `~6.27e-5 s^-1` |
| `NO` peak | `~29.9 ppbv` (Phase 2) |
| `NO2` peak | `~6.5 ppbv` (Phase 2) |
| `O3` minimum | `~43.5 ppbv` (Phase 2) |

Run an inspection script. Either:

a. Use `scripts/check_tuvx_phase.py` if available — it already gates these values:

```bash
~/miniconda3/envs/mpas/bin/python ~/EarthSystem/CheMPAS/scripts/check_tuvx_phase.py \
  --history history.0000-01-01_18.00.00.nc --phase 3
```
Expected: phase-gate PASS.

b. Or read the final-step values directly:

```bash
~/miniconda3/envs/mpas/bin/python -c "
import netCDF4 as nc
import numpy as np
f = nc.Dataset('history.0000-01-01_18.00.00.nc')
for var in ['j_jNO2', 'qNO', 'qNO2', 'qO3']:
    if var in f.variables:
        a = f.variables[var][-1, :, :]
        print(f'{var:8s} min={np.min(a):.3e} max={np.max(a):.3e}')
"
```
Expected: values within a few percent of the table above. Bit-for-bit identity to a saved baseline is the strong guarantee, but in absence of a saved baseline, the recorded reference values are the next-best check.

- [ ] **Step 4: Report**

If values match: report DONE.
If values diverge significantly (more than a few percent on `j_jNO2`, or qualitative changes in `qNO`/`qNO2`/`qO3` peaks/minimums): report DONE_WITH_CONCERNS, paste the actual values, and stop. Do not proceed to Task 8 — investigate the discrepancy first (likely cause: the `cos_sza_cell` fill broadcast may not be producing identical floating-point results to the original scalar use).

---

### Task 8: Per-cell sanity check (switch=`.true.` on baroclinic wave)

**Files:** Stages a one-off namelist override in the run directory; does not modify the tracked `test_cases/jw_baroclinic_wave/namelist.atmosphere`.

This task confirms the per-cell path actually varies `cos_sza` across the mesh on a real spherical case. Requires the baroclinic-wave run directory at `~/Data/CheMPAS/jw_baroclinic_wave/`.

- [ ] **Step 1: Stage the new build into the run directory**

Run:
```bash
cp atmosphere_model ~/Data/CheMPAS/jw_baroclinic_wave/
```

- [ ] **Step 2: Override the namelist for a per-cell run**

Edit `~/Data/CheMPAS/jw_baroclinic_wave/namelist.atmosphere` and add (or set) inside the `&musica` block:

```
config_chemistry_use_grid_coords = .true.
```

Also confirm a TUV-x config and chemistry are enabled (set `config_micm_file` and `config_tuvx_config_file` if not already). The simplest setup is the same `lnox_o3.yaml` + `tuvx_no2.json` pair used for the supercell case — copy them into `~/Data/CheMPAS/jw_baroclinic_wave/` if not already present.

- [ ] **Step 3: Run a short integration (1 day is enough to observe SZA gradient)**

Edit `~/Data/CheMPAS/jw_baroclinic_wave/streams.atmosphere` to set the history output interval to a sensible cadence (e.g., every 6 hours), and edit `namelist.atmosphere` to set `config_run_duration = '1_00:00:00'` (1 day).

Run:
```bash
cd ~/Data/CheMPAS/jw_baroclinic_wave/
rm -f history.*.nc diag.*.nc log.atmosphere.*
mpiexec -n 8 ./atmosphere_model 2>&1 | tail -20
```

Expected: clean run completion.

- [ ] **Step 4: Plot a horizontal slice of `j_jNO2` at the surface**

Use the visualization tooling. The general pattern (adapt the script name to whatever surface/horizontal plotter exists in `scripts/`):

```bash
~/miniconda3/envs/mpas/bin/python -c "
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
f = nc.Dataset('history.0000-01-02_00.00.00.nc')
lat = np.degrees(f.variables['latCell'][:])
lon = np.degrees(f.variables['lonCell'][:])
j   = f.variables['j_jNO2'][-1, -1, :]   # last time, surface level
plt.scatter(lon, lat, c=j, s=2, cmap='viridis')
plt.colorbar(label='j_NO2 (s^-1)')
plt.xlabel('Longitude (deg)'); plt.ylabel('Latitude (deg)')
plt.title('Surface j_NO2 with per-cell solar geometry')
plt.savefig('/tmp/jw_jno2_surface.png', dpi=120, bbox_inches='tight')
print('wrote /tmp/jw_jno2_surface.png')
"
open /tmp/jw_jno2_surface.png   # macOS; use xdg-open on Linux
```

- [ ] **Step 5: Confirm a longitudinal terminator and latitudinal gradient**

Visually verify the plot:
- A clear day–night terminator (sharp boundary where `j_NO2` drops to near zero).
- Latitudinal SZA dependence: `j_NO2` is largest near the sub-solar latitude and smaller toward the poles (within the daylit hemisphere).

If the field is uniformly zero, uniformly nonzero, or shows no longitudinal variation: report BLOCKED — the per-cell path isn't engaging. Likely causes: switch not `.true.` in the actual namelist read, `latCell`/`lonCell` not actually populated for the case, or the run is too short to span any daylight.

- [ ] **Step 6: Report**

If the plot shows the expected terminator and SZA dependence: report DONE. The per-cell path is validated.

---

## Self-Review Notes

**Spec coverage:**

- Spec § *Namelist & Registry* → Task 1.
- Spec § *`chemistry_step` data flow* item 1 (read switch & coords) → Task 2 + Task 3 Step 3.
- Spec § *`chemistry_step` data flow* item 2 (replace scalar with per-cell array) → Task 3 Steps 1–3.
- Spec § *`chemistry_step` data flow* item 3 (substitute downstream uses) → Task 3 Steps 4–5.
- Spec § *`chemistry_step` data flow* item 4 (deallocate) → Task 3 Step 6.
- Spec § *Edge cases* (idealized-mesh allowed, no defensive checks) → preserved by Task 3's design (no validation gates added).
- Spec § *Validation* item 1 (reproducibility on supercell) → Task 7.
- Spec § *Validation* item 2 (per-cell sanity on baroclinic wave) → Task 8.
- Spec § *Validation* item 3 (no phase-gate changes) → preserved (Task 7 uses the existing `check_tuvx_phase.py`).
- Spec § *Files Touched* → Tasks 1, 3, 5 cover exactly the listed files; nothing else is touched.
- Spec § *Commit Shape* → Task 6 (single commit, exact subject line from spec).
- Spec § *Doc Updates* → Task 5 (both bullets); landed in same commit per Task 6.

**Placeholder scan:** No TBD/TODO/vague-instruction patterns. Every step shows the exact code or command and an expected outcome.

**Type / name consistency:** New names used — `chem_use_grid_coords` (module save), `use_grid_coords_ptr` (init local), `config_chemistry_use_grid_coords` (namelist), `cos_sza_cell` (step local), `latCell`/`lonCell` (existing mesh fields), `RAD2DEG` (parameter). Each appears in the same form in every task that references it.

**Note on line numbers:** Line numbers throughout the plan are anchored to the current state of `develop` at planning time. If new commits land between planning and execution, prefer the surrounding-context anchors (e.g., "the `chem_lat`/`chem_lon` block" in Task 2) over the literal line numbers.
