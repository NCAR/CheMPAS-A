# Phase 1: Solar Geometry for Photolysis — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the constant `config_lnox_j_no2` photolysis rate with `j_max * max(0, cos_sza)`, where SZA is computed from model time and namelist lat/lon coordinates.

**Architecture:** Add `mpas_solar_geometry.F` as a standalone SZA calculator. Thread model time through the chemistry call path. In `mpas_musica.F`, cache the photolysis rate-parameter index at init, then update it every chemistry timestep via a new `musica_set_photolysis` entry point. `config_lnox_j_no2` becomes the peak daytime rate (j_max).

**Tech Stack:** Fortran 2008, MPAS timekeeping (`mpas_get_time`), MUSICA-Fortran rate parameters API.

---

### Task 1: Add namelist parameters for lat/lon

**Files:**
- Modify: `src/core_atmosphere/Registry.xml:397-431`
- Modify: `test_cases/supercell/namelist.atmosphere:56-65`

**Step 1: Add `config_chemistry_latitude` and `config_chemistry_longitude` to Registry.xml**

In `Registry.xml`, inside the `<nml_record name="musica">` block, add two new options after the existing `config_lnox_nox_tau` entry:

```xml
                <nml_option name="config_chemistry_latitude" type="real" default_value="0.0"
                     units="degrees"
                     description="Latitude for solar geometry in idealized cases (degrees N)"
                     possible_values="Any real between -90 and 90"/>
                <nml_option name="config_chemistry_longitude" type="real" default_value="0.0"
                     units="degrees"
                     description="Longitude for solar geometry in idealized cases (degrees E)"
                     possible_values="Any real between -180 and 360"/>
```

**Step 2: Update supercell namelist with Kingfisher, OK coordinates and DC3 start time**

Edit `test_cases/supercell/namelist.atmosphere`:

```
&nhyd_model
    config_start_time = '2012-05-29_21:00:00'
```

```
&musica
    ...
    config_chemistry_latitude = 35.86
    config_chemistry_longitude = -97.93
/
```

**Step 3: Commit**

```bash
git add src/core_atmosphere/Registry.xml test_cases/supercell/namelist.atmosphere
git commit -m "feat(phase1): add lat/lon namelist params and DC3 start time"
```

---

### Task 2: Create `mpas_solar_geometry.F`

**Files:**
- Create: `src/core_atmosphere/chemistry/mpas_solar_geometry.F`
- Modify: `src/core_atmosphere/chemistry/Makefile`

**Step 1: Write the solar geometry module**

This module computes `cos(SZA)` from year, day-of-year, hour (UTC), latitude, and longitude using the Spencer (1971) algorithm for solar declination.

```fortran
! src/core_atmosphere/chemistry/mpas_solar_geometry.F
module mpas_solar_geometry

    use mpas_kind_types, only : RKIND

    implicit none
    private
    public :: solar_cos_sza

    real(kind=RKIND), parameter :: PI = 3.14159265358979323846_RKIND
    real(kind=RKIND), parameter :: DEG2RAD = PI / 180.0_RKIND

contains

    !> Compute cos(SZA) from date/time and location.
    !>
    !> Uses Spencer (1971) formula for solar declination.
    !> Returns cos(SZA); negative values mean sun is below horizon.
    !>
    !> @param[in] DoY        Day of year (1-366)
    !> @param[in] hour_utc   Hour of day in UTC (fractional, 0-24)
    !> @param[in] lat_deg    Latitude in degrees North
    !> @param[in] lon_deg    Longitude in degrees East
    !> @return cos_sza       Cosine of solar zenith angle
    function solar_cos_sza(DoY, hour_utc, lat_deg, lon_deg) result(cos_sza)

        integer, intent(in)           :: DoY
        real(kind=RKIND), intent(in)  :: hour_utc
        real(kind=RKIND), intent(in)  :: lat_deg
        real(kind=RKIND), intent(in)  :: lon_deg
        real(kind=RKIND)              :: cos_sza

        real(kind=RKIND) :: gamma, decl, eqtime, hour_angle, lat_rad

        ! Fractional year [radians] — Spencer (1971)
        gamma = 2.0_RKIND * PI * (real(DoY - 1, RKIND) + (hour_utc - 12.0_RKIND) / 24.0_RKIND) / 365.0_RKIND

        ! Solar declination [radians]
        decl = 0.006918_RKIND &
             - 0.399912_RKIND * cos(gamma) &
             + 0.070257_RKIND * sin(gamma) &
             - 0.006758_RKIND * cos(2.0_RKIND * gamma) &
             + 0.000907_RKIND * sin(2.0_RKIND * gamma) &
             - 0.002697_RKIND * cos(3.0_RKIND * gamma) &
             + 0.001480_RKIND * sin(3.0_RKIND * gamma)

        ! Equation of time [minutes]
        eqtime = 229.18_RKIND * ( &
                   0.000075_RKIND &
                 + 0.001868_RKIND * cos(gamma) &
                 - 0.032077_RKIND * sin(gamma) &
                 - 0.014615_RKIND * cos(2.0_RKIND * gamma) &
                 - 0.040849_RKIND * sin(2.0_RKIND * gamma) )

        ! Hour angle [radians]
        ! Solar noon occurs when hour_angle = 0
        hour_angle = DEG2RAD * ((hour_utc * 60.0_RKIND + eqtime + 4.0_RKIND * lon_deg) / 4.0_RKIND - 180.0_RKIND)

        ! Latitude in radians
        lat_rad = lat_deg * DEG2RAD

        ! Cosine of solar zenith angle
        cos_sza = sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(hour_angle)

    end function solar_cos_sza

end module mpas_solar_geometry
```

**Step 2: Add to Makefile**

In `src/core_atmosphere/chemistry/Makefile`, add `mpas_solar_geometry.o` to OBJS and make `mpas_atm_chemistry.o` depend on it:

```makefile
OBJS = \
	mpas_solar_geometry.o \
	mpas_lightning_nox.o \
	mpas_atm_chemistry.o
```

And update the dependency:

```makefile
mpas_atm_chemistry.o: mpas_lightning_nox.o mpas_solar_geometry.o
```

**Step 3: Commit**

```bash
git add src/core_atmosphere/chemistry/mpas_solar_geometry.F src/core_atmosphere/chemistry/Makefile
git commit -m "feat(phase1): add mpas_solar_geometry module with Spencer SZA"
```

---

### Task 3: Thread model time into chemistry call path

**Files:**
- Modify: `src/core_atmosphere/mpas_atm_core.F:1020-1057`
- Modify: `src/core_atmosphere/chemistry/mpas_atm_chemistry.F:138`

**Step 1: Pass `currTime` to `chemistry_step` from `atm_do_timestep`**

In `mpas_atm_core.F`, the `currTime` variable already exists (line 1030). Add it to the `chemistry_step` call:

```fortran
! In the DO_CHEMISTRY block declarations (after line 1023):
      type (MPAS_pool_type), pointer :: configs_chem
```

```fortran
! In the DO_CHEMISTRY block (replace lines 1046-1057):
      block_chem => domain % blocklist
      call mpas_pool_get_subpool(domain % blocklist % structs, 'configs', configs_chem)
      do while (associated(block_chem))
         call mpas_pool_get_subpool(block_chem % structs, 'mesh', mesh_chem)
         call mpas_pool_get_subpool(block_chem % structs, 'state', state_chem)
         call mpas_pool_get_subpool(block_chem % structs, 'diag', diag_chem)

         time_lev_chem = 1
         call chemistry_step(dt, currTime, configs_chem, &
                             mesh_chem, state_chem, diag_chem, &
                             block_chem % dimensions, time_lev_chem)

         block_chem => block_chem % next
      end do
```

**Step 2: Update `chemistry_step` signature**

In `mpas_atm_chemistry.F`, update the subroutine signature to accept `currTime` and `configs`:

```fortran
    subroutine chemistry_step(dt, currTime, configs, mesh, state, diag, dimensions, time_lev)

#ifdef MPAS_USE_MUSICA
        use iso_fortran_env, only: real64
        use mpas_musica, only: musica_step, musica_step_ref, log_column_comparison, &
                               copy_state_to_ref, musica_set_photolysis
        use mpas_lightning_nox, only: lightning_nox_inject
        use mpas_solar_geometry, only: solar_cos_sza
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT, MPAS_Time_Type
        use mpas_pool_routines, only: mpas_pool_get_dimension, mpas_pool_get_array, &
                                      mpas_pool_get_config
        use mpas_timekeeping, only: mpas_get_time

        real (kind=RKIND), intent(in) :: dt
        type (MPAS_Time_Type), intent(in) :: currTime
        type (mpas_pool_type), intent(in) :: configs
        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(inout) :: state
        type (mpas_pool_type), intent(in) :: diag
        type (mpas_pool_type), intent(in) :: dimensions
        integer, intent(in) :: time_lev
```

Add SZA computation and photolysis update before lightning injection (after the `first_call` block and dimension reads, before Step 0):

```fortran
#ifdef MPAS_USE_MUSICA
        ...existing declarations...

        ! Solar geometry variables
        integer :: year, DoY, hour, minute, second
        real(kind=RKIND) :: hour_utc, cos_sza, j_no2_max
        real(kind=RKIND), pointer :: lat_ptr, lon_ptr, j_no2_ptr
        real(kind=RKIND), save :: chem_lat = 0.0_RKIND
        real(kind=RKIND), save :: chem_lon = 0.0_RKIND
        real(kind=RKIND), save :: chem_j_no2_max = 0.0_RKIND
        logical, save :: solar_params_cached = .false.
```

After dimension reads, before Step 0:

```fortran
        ! Cache solar geometry parameters on first call
        if (.not. solar_params_cached) then
            nullify(lat_ptr, lon_ptr, j_no2_ptr)
            call mpas_pool_get_config(configs, 'config_chemistry_latitude', lat_ptr)
            if (associated(lat_ptr)) chem_lat = lat_ptr
            call mpas_pool_get_config(configs, 'config_chemistry_longitude', lon_ptr)
            if (associated(lon_ptr)) chem_lon = lon_ptr
            call mpas_pool_get_config(configs, 'config_lnox_j_no2', j_no2_ptr)
            if (associated(j_no2_ptr)) chem_j_no2_max = j_no2_ptr
            solar_params_cached = .true.
        end if

        ! Compute SZA and update photolysis rate
        call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
        hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                 + real(second, RKIND) / 3600.0_RKIND
        cos_sza = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)

        call musica_set_photolysis(chem_j_no2_max, cos_sza, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            return
        end if
```

**Step 3: Commit**

```bash
git add src/core_atmosphere/mpas_atm_core.F src/core_atmosphere/chemistry/mpas_atm_chemistry.F
git commit -m "feat(phase1): thread model time into chemistry, compute SZA each step"
```

---

### Task 4: Add `musica_set_photolysis` to `mpas_musica.F`

**Files:**
- Modify: `src/core_atmosphere/chemistry/musica/mpas_musica.F`

**Step 1: Add module-level cached index and public routine**

Add to the module's private data (after `state_ref` declaration):

```fortran
    integer, save :: photo_no2_rp_index = -1   ! cached rate-parameter index for PHOTO.no2_photolysis
```

Add to the `public` list:

```fortran
    public :: musica_set_photolysis
```

**Step 2: Cache the rate parameter index at end of `musica_init`**

At the end of `musica_init`, before `deallocate(init_conc)`:

```fortran
        ! Cache photolysis rate parameter index for per-step updates
        call cache_photo_rp_index()
```

**Step 3: Write `cache_photo_rp_index` subroutine**

```fortran
    subroutine cache_photo_rp_index()
        use mpas_log, only : mpas_log_write
        integer :: i, n_rp
        character(len=:), allocatable :: rp_name

        if (.not. associated(state)) return
        if (state%number_of_rate_parameters == 0) return

        n_rp = state%rate_parameters_ordering%size()
        do i = 1, n_rp
            rp_name = trim(state%rate_parameters_ordering%name(i))
            if (rp_name == 'PHOTO.no2_photolysis') then
                photo_no2_rp_index = i
                call mpas_log_write('[MUSICA] Cached PHOTO.no2_photolysis rate param index = $i', &
                    intArgs=[i])
                return
            end if
        end do

        call mpas_log_write('[MUSICA] PHOTO.no2_photolysis not found in rate parameters (SZA scaling disabled)')

    end subroutine cache_photo_rp_index
```

**Step 4: Write `musica_set_photolysis` subroutine**

```fortran
    subroutine musica_set_photolysis(j_no2_max, cos_sza, error_code, error_message)
        use iso_fortran_env, only : real64
        use mpas_kind_types, only : RKIND
        use mpas_log, only : mpas_log_write

        real(kind=RKIND), intent(in) :: j_no2_max   ! peak daytime j_NO2 [s-1]
        real(kind=RKIND), intent(in) :: cos_sza      ! cosine of solar zenith angle
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        real(real64) :: j_val
        integer :: cell_stride, var_stride, cell, idx, n_cells
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        if (photo_no2_rp_index < 1) return  ! no photolysis in mechanism
        if (.not. associated(state)) return

        j_val = real(j_no2_max, real64) * max(0.0_real64, real(cos_sza, real64))

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride  = state%rate_parameters_strides%variable
        n_cells = state%number_of_grid_cells

        do cell = 1, n_cells
            idx = 1 + (cell - 1) * cell_stride + (photo_no2_rp_index - 1) * var_stride
            state%rate_parameters(idx) = j_val
        end do

        ! Also update reference state
        if (associated(state_ref)) then
            cell_stride = state_ref%rate_parameters_strides%grid_cell
            var_stride  = state_ref%rate_parameters_strides%variable
            do cell = 1, state_ref%number_of_grid_cells
                idx = 1 + (cell - 1) * cell_stride + (photo_no2_rp_index - 1) * var_stride
                state_ref%rate_parameters(idx) = j_val
            end do
        end if

        write(msg, '(A,F8.4,A,ES10.3)') '[MUSICA] cos_sza=', cos_sza, ' j_NO2=', j_val
        call mpas_log_write(trim(msg))

    end subroutine musica_set_photolysis
```

**Step 5: Remove `j_no2_prescribed` from `musica_init` signature**

Since photolysis is now set per-step, `musica_init` no longer needs it. Change:
- `musica_init` signature: remove `j_no2_prescribed`, keep `nox_tau`
- `assign_rate_parameters`: remove `j_no2_val` parameter, skip the PHOTO block (let `musica_set_photolysis` handle it)
- `mpas_atm_chemistry.F`: update `chemistry_init` to match (stop reading/passing `config_lnox_j_no2` to `musica_init`)

Actually — **keep it simple**. Leave `assign_rate_parameters` setting the initial value at init. The per-step `musica_set_photolysis` will overwrite it immediately on the first timestep. No signature changes needed — this avoids unnecessary churn.

**Step 6: Commit**

```bash
git add src/core_atmosphere/chemistry/musica/mpas_musica.F
git commit -m "feat(phase1): add musica_set_photolysis for per-step rate updates"
```

---

### Task 5: Build and verify

**Step 1: Build**

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"
make clean CORE=atmosphere && find . -name "*.mod" -delete && find . -name "*.o" -delete
make -j8 llvm CORE=atmosphere PIO=$HOME/software NETCDF=/opt/homebrew PNETCDF=$HOME/software PRECISION=double MUSICA=true
```

Expected: Clean build, `atmosphere_model` produced.

**Step 2: Commit if any fixups needed**

---

### Task 6: Run supercell test and verify SZA

**Step 1: Update run directory namelist**

Copy the updated reference namelist to the run directory:

```bash
cp test_cases/supercell/namelist.atmosphere ~/Data/MPAS/supercell/namelist.atmosphere
```

**Step 2: Re-initialize tracers for Case B**

```bash
cd ~/Data/MPAS/supercell
~/miniconda3/envs/mpas/bin/python ~/EarthSystem/CheMPAS/scripts/init_lnox_o3.py -i supercell_init.nc
```

**Step 3: Run 30-minute Case B**

```bash
cd ~/Data/MPAS/supercell
ts=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${ts}.nc
[ -f log.atmosphere.0000.out ] && mv log.atmosphere.0000.out log.atmosphere.0000.${ts}.out
mpiexec -n 8 ~/EarthSystem/CheMPAS/atmosphere_model
```

**Step 4: Verify SZA in log**

```bash
grep 'cos_sza' log.atmosphere.0000.out | head -5
```

Expected: `cos_sza ≈ 0.81` (Kingfisher, OK at 21:00 UTC May 29 = SZA ≈ 35.7°), slowly decreasing over 30 minutes as sun descends.

**Step 5: Verify j_NO2 is non-zero during daytime**

```bash
grep 'j_NO2' log.atmosphere.0000.out | head -5
```

Expected: `j_NO2 ≈ 0.0081` (= 0.01 * 0.81).

**Step 6: Check for non-negative tracers**

```bash
~/miniconda3/envs/mpas/bin/python -c "
import netCDF4 as nc
ds = nc.Dataset('output.nc')
for v in ['qNO','qNO2','qO3']:
    d = ds.variables[v][:]
    print(f'{v}: min={d.min():.6e} max={d.max():.6e}')
"
```

Expected: All min values >= 0.

---

### Task 7: Update plan and commit final Phase 1

**Step 1: Update the TUV-x plan document**

Mark Phase 1 checklist items as complete in `docs/plans/2026-03-06-tuvx-photolysis-integration.md`.

**Step 2: Update docs/project/TODO.md**

Mark "Solar geometry / day-night physics for j_NO2 (Phase 1)" as complete.

**Step 3: Commit**

```bash
git add -A
git commit -m "feat(phase1): solar geometry for SZA-scaled photolysis - complete"
```

---

## Summary of Changes

| File | Change |
|------|--------|
| `Registry.xml` | Add `config_chemistry_latitude`, `config_chemistry_longitude` |
| `mpas_solar_geometry.F` | **NEW** — Spencer SZA calculator |
| `chemistry/Makefile` | Add `mpas_solar_geometry.o` |
| `mpas_atm_core.F` | Pass `currTime` and `configs` to `chemistry_step` |
| `mpas_atm_chemistry.F` | Accept `currTime`/`configs`, compute SZA, call `musica_set_photolysis` |
| `mpas_musica.F` | Add `musica_set_photolysis`, cache `photo_no2_rp_index` |
| `test_cases/supercell/namelist.atmosphere` | DC3 start time, Kingfisher lat/lon |

## Key Design Decisions

1. **`config_lnox_j_no2` becomes j_max** — the peak daytime rate. Actual rate = j_max × max(0, cos_sza).
2. **Per-step update** — `musica_set_photolysis` overwrites `PHOTO.no2_photolysis` every timestep for both coupled and reference states.
3. **Fallback SZA** — Chemistry-side `mpas_solar_geometry.F` computes SZA from namelist lat/lon. No dependency on physics radiation.
4. **Uniform SZA** — All cells get the same cos_sza (appropriate for idealized cases with single lat/lon).
