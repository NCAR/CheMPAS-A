# Runtime Chemistry Tracer Allocation — Implementation Plan

## Document Status

- `Historical Context:` Task-by-task implementation plan used to deliver runtime
  chemistry tracer allocation.
- `Current State:` Complete and merged to `main` (as of March 6, 2026).
- `Use This As:` Implementation history and design rationale, not active task
  instructions.
- `Original Authoring Note:` "For Claude: REQUIRED SUB-SKILL: Use
  superpowers:executing-plans to implement this plan task-by-task."

**Goal:** Remove chemistry tracers from Registry.xml and inject them dynamically at runtime based on the MICM config, so switching chemistry mechanisms requires zero source code changes.

**Architecture:** Two-phase init. A lightweight early scan in `atm_setup_block` creates a temporary MICM instance to discover species names/count, then extends the existing `scalars` and `scalars_tend` var_arrays. The full MICM solver initialization happens later in `chemistry_init` as before.

**Tech Stack:** Fortran 2008, MPAS framework pool system, MUSICA-Fortran/MICM API

---

## Task 1: Add `musica_query_species` to `mpas_musica.F`

**Files:**
- Modify: `src/core_atmosphere/chemistry/musica/mpas_musica.F:18-31`

**Step 1: Add the new public subroutine**

Add `musica_query_species` to the public list (line 31) and implement it after `resolve_mpas_indices`. This routine creates a temporary `micm_t`, queries species names, and destroys it.

```fortran
public :: musica_query_species
```

Implement below `resolve_mpas_indices` (after line 215):

```fortran
!------------------------------------------------------------------------
!  subroutine musica_query_species
!
!> \brief Query MICM config for chemistry species names without full init
!> \details
!>  Creates a temporary MICM instance to read species_ordering from
!>  the config file. Returns species count and MPAS tracer names
!>  (with 'q' prefix). The temporary MICM instance is destroyed.
!>  Used by atm_setup_block to extend scalars before chemistry_init.
!------------------------------------------------------------------------
subroutine musica_query_species(config_path, num_species, species_names, error_code, error_message)

    use musica_micm,  only: micm_t, RosenbrockStandardOrder
    use musica_state, only: state_t
    use musica_util,  only: error_t
    use mpas_log,     only: mpas_log_write

    character(len=*), intent(in)               :: config_path
    integer, intent(out)                       :: num_species
    character(len=64), allocatable, intent(out) :: species_names(:)
    integer, intent(out)                       :: error_code
    character(len=:), allocatable, intent(out)  :: error_message

    type(micm_t), pointer  :: tmp_micm
    type(state_t), pointer :: tmp_state
    type(error_t)          :: error
    integer                :: i

    error_code = 0
    error_message = ''
    num_species = 0

    if (len_trim(config_path) == 0) then
        call mpas_log_write('[MUSICA] No MICM config path; skipping chemistry tracer discovery')
        return
    end if

    call mpas_log_write('[MUSICA] Querying MICM config for species list...')

    tmp_micm => micm_t(trim(config_path), RosenbrockStandardOrder, error)
    if (has_error_occurred(error, error_message, error_code)) return

    ! We need a state just to access species_ordering
    tmp_state => tmp_micm%get_state(1, error)
    if (has_error_occurred(error, error_message, error_code)) return

    num_species = tmp_state%species_ordering%size()
    allocate(species_names(num_species))

    do i = 1, num_species
        species_names(i) = 'q' // trim(tmp_state%species_ordering%name(i))
        call mpas_log_write('[MUSICA] Discovered chemistry tracer: ' // trim(species_names(i)))
    end do

    ! Clean up temporary instances
    deallocate(tmp_state)
    deallocate(tmp_micm)

end subroutine musica_query_species
```

**Step 2: Build to verify compilation**

```bash
make -j8 llvm CORE=atmosphere PIO=$HOME/software NETCDF=/opt/homebrew PNETCDF=$HOME/software PRECISION=double MUSICA=true 2>&1 | tail -20
```

Expected: Compiles successfully (new subroutine not yet called).

**Step 3: Commit**

```bash
git add src/core_atmosphere/chemistry/musica/mpas_musica.F
git commit -m "Add musica_query_species for lightweight MICM config scan"
```

---

## Task 2: Add `atm_extend_scalars_for_chemistry` to `mpas_atm_core_interface.F`

**Files:**
- Modify: `src/core_atmosphere/mpas_atm_core_interface.F:461-498` (atm_setup_block)
- Modify: `src/core_atmosphere/mpas_atm_core_interface.F:517-632` (new routine after atm_allocate_scalars)

**Step 1: Add the extension routine**

Add this new function after `atm_allocate_scalars` (after line 632):

```fortran
#ifdef MPAS_USE_MUSICA
!-----------------------------------------------------------------------
!
!  function atm_extend_scalars_for_chemistry
!
!> \brief   Extend scalars/scalars_tend with chemistry tracers from MICM
!> \author  CheMPAS-A Developers
!> \date    March 2026
!> \details
!>  Queries MICM config for species names, then extends the existing
!>  scalars and scalars_tend var_arrays by appending chemistry tracers.
!>  Called from atm_setup_block after atm_generate_structs.
!
!-----------------------------------------------------------------------
function atm_extend_scalars_for_chemistry(block) result(ierr)

    use mpas_derived_types, only : block_type, mpas_pool_type, field3dReal, &
                                   MPAS_LOG_ERR, MPAS_POOL_SILENT
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_config, &
                                   mpas_pool_get_field, mpas_pool_get_dimension, &
                                   mpas_pool_add_dimension, &
                                   mpas_pool_set_error_level, mpas_pool_get_error_level
    use mpas_kind_types,    only : StrKIND, RKIND
    use mpas_log,           only : mpas_log_write
    use mpas_musica,        only : musica_query_species

    implicit none

    type (block_type), pointer :: block
    integer :: ierr

    ! Local variables
    character(len=StrKIND), pointer :: config_micm_file
    integer :: num_chem_species
    character(len=64), allocatable :: chem_names(:)
    integer :: error_code
    character(len=:), allocatable :: error_message

    type (mpas_pool_type), pointer :: statePool, tendPool
    type (field3dReal), dimension(:), pointer :: scalarsField, tendField
    integer, pointer :: num_scalars_ptr, num_tend_ptr

    real (kind=RKIND), dimension(:,:,:), pointer :: old_array, new_array
    character(len=StrKIND), dimension(:), pointer :: old_names, new_names
    type (att_lists_type), dimension(:), pointer :: old_atts, new_atts

    integer :: old_count, new_count, timeLevs, t, i, j
    integer :: nVertLevels, nCells
    integer, pointer :: nVertLevels_ptr, nCells_ptr
    integer :: err_level

    ierr = 0

    ! Get MICM config path from namelist
    nullify(config_micm_file)
    err_level = mpas_pool_get_error_level()
    call mpas_pool_set_error_level(MPAS_POOL_SILENT)
    call mpas_pool_get_config(block % domain % configs, 'config_micm_file', config_micm_file)
    call mpas_pool_set_error_level(err_level)

    if (.not. associated(config_micm_file)) then
        call mpas_log_write('[Chemistry] No config_micm_file found; skipping dynamic tracer allocation')
        return
    end if

    if (len_trim(config_micm_file) == 0) then
        call mpas_log_write('[Chemistry] config_micm_file is empty; skipping dynamic tracer allocation')
        return
    end if

    ! Query MICM for species names
    call musica_query_species(config_micm_file, num_chem_species, chem_names, error_code, error_message)
    if (error_code /= 0) then
        call mpas_log_write('[Chemistry] Error querying MICM species: ' // trim(error_message), &
                            messageType=MPAS_LOG_ERR)
        ierr = error_code
        return
    end if

    if (num_chem_species == 0) then
        call mpas_log_write('[Chemistry] MICM config has no species; no tracers to add')
        return
    end if

    call mpas_log_write('[Chemistry] Extending scalars with $i chemistry tracers', &
                        intArgs=[num_chem_species])

    ! Get dimensions
    call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels_ptr)
    call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells_ptr)
    nVertLevels = nVertLevels_ptr
    nCells = nCells_ptr

    ! --- Extend scalars in state pool ---
    nullify(statePool)
    call mpas_pool_get_subpool(block % structs, 'state', statePool)
    if (.not. associated(statePool)) then
        call mpas_log_write('[Chemistry] No state pool found', messageType=MPAS_LOG_ERR)
        ierr = 1
        return
    end if

    timeLevs = 2
    nullify(scalarsField)
    call mpas_pool_get_field(statePool, 'scalars', scalarsField, 1)
    if (.not. associated(scalarsField)) then
        call mpas_log_write('[Chemistry] No scalars field found in state pool', messageType=MPAS_LOG_ERR)
        ierr = 1
        return
    end if

    ! Current size of first dimension
    old_count = size(scalarsField(1) % array, 1)
    new_count = old_count + num_chem_species

    call mpas_log_write('[Chemistry] scalars: old_count=$i new_count=$i', &
                        intArgs=[old_count, new_count])

    ! Update num_scalars dimension
    call mpas_pool_add_dimension(statePool, 'num_scalars', new_count)

    do t = 1, timeLevs
        ! Save old array pointer
        old_array => scalarsField(t) % array

        ! Allocate new larger array
        allocate(new_array(new_count, nVertLevels, nCells))
        new_array = 0.0_RKIND

        ! Copy existing tracer data
        new_array(1:old_count, :, :) = old_array(1:old_count, :, :)

        ! Swap pointers
        scalarsField(t) % array => new_array
        deallocate(old_array)

        ! Update dimSizes
        scalarsField(t) % dimSizes(1) = new_count

        ! Extend constituentNames
        old_names => scalarsField(t) % constituentNames
        allocate(new_names(new_count))
        new_names(1:old_count) = old_names(1:old_count)
        do i = 1, num_chem_species
            new_names(old_count + i) = trim(chem_names(i))
        end do
        scalarsField(t) % constituentNames => new_names
        deallocate(old_names)

        ! Extend attLists
        old_atts => scalarsField(t) % attLists
        allocate(new_atts(new_count))
        do j = 1, old_count
            new_atts(j) = old_atts(j)
        end do
        do j = old_count + 1, new_count
            allocate(new_atts(j) % attList)
        end do
        scalarsField(t) % attLists => new_atts
        ! Don't deallocate old_atts elements — they were shallow-copied
        deallocate(old_atts)
    end do

    ! Add index dimensions for each chemistry tracer
    do i = 1, num_chem_species
        call mpas_pool_add_dimension(statePool, 'index_' // trim(chem_names(i)), old_count + i)
        call mpas_log_write('[Chemistry] Added index_' // trim(chem_names(i)) // ' = $i', &
                            intArgs=[old_count + i])
    end do

    ! --- Extend scalars_tend in tend pool ---
    nullify(tendPool)
    call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
    if (.not. associated(tendPool)) then
        call mpas_log_write('[Chemistry] No tend pool found', messageType=MPAS_LOG_ERR)
        ierr = 1
        return
    end if

    nullify(tendField)
    call mpas_pool_get_field(tendPool, 'scalars_tend', tendField, 1)
    if (.not. associated(tendField)) then
        call mpas_log_write('[Chemistry] No scalars_tend field found in tend pool', messageType=MPAS_LOG_ERR)
        ierr = 1
        return
    end if

    old_count = size(tendField(1) % array, 1)
    new_count = old_count + num_chem_species

    call mpas_pool_add_dimension(tendPool, 'num_scalars_tend', new_count)

    ! scalars_tend has 1 time level
    old_array => tendField(1) % array
    allocate(new_array(new_count, nVertLevels, nCells))
    new_array = 0.0_RKIND
    new_array(1:old_count, :, :) = old_array(1:old_count, :, :)
    tendField(1) % array => new_array
    deallocate(old_array)
    tendField(1) % dimSizes(1) = new_count

    ! Extend constituentNames for tend
    old_names => tendField(1) % constituentNames
    allocate(new_names(new_count))
    new_names(1:old_count) = old_names(1:old_count)
    do i = 1, num_chem_species
        new_names(old_count + i) = trim(chem_names(i))
    end do
    tendField(1) % constituentNames => new_names
    deallocate(old_names)

    ! Extend attLists for tend
    old_atts => tendField(1) % attLists
    allocate(new_atts(new_count))
    do j = 1, old_count
        new_atts(j) = old_atts(j)
    end do
    do j = old_count + 1, new_count
        allocate(new_atts(j) % attList)
    end do
    tendField(1) % attLists => new_atts
    deallocate(old_atts)

    ! Add index dimensions for tend pool too
    do i = 1, num_chem_species
        call mpas_pool_add_dimension(tendPool, 'index_' // trim(chem_names(i)), old_count + i)
    end do

    if (allocated(chem_names)) deallocate(chem_names)

    call mpas_log_write('[Chemistry] Dynamic tracer allocation complete')

end function atm_extend_scalars_for_chemistry
#endif
```

**Step 2: Call it from `atm_setup_block`**

In `atm_setup_block` (line 496), add the call after the CAM block:

```fortran
      ! After the existing cam_pcnst block, add:
#ifdef MPAS_USE_MUSICA
      ierr = atm_extend_scalars_for_chemistry(block)
      if (ierr /= 0) then
         call mpas_log_write('** Error extending scalars for chemistry', messageType=MPAS_LOG_ERR)
         return
      end if
#endif
```

Note: The `use` for `att_lists_type` may be needed. Check `mpas_derived_types` for the correct import.

**Step 3: Build to verify compilation**

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"
make clean CORE=atmosphere && find . -name "*.mod" -delete && find . -name "*.o" -delete
make -j8 llvm CORE=atmosphere PIO=$HOME/software NETCDF=/opt/homebrew PNETCDF=$HOME/software PRECISION=double MUSICA=true 2>&1 | tail -30
```

Expected: Compiles successfully.

**Step 4: Commit**

```bash
git add src/core_atmosphere/mpas_atm_core_interface.F
git commit -m "Add atm_extend_scalars_for_chemistry for runtime tracer injection"
```

---

## Task 3: Remove chemistry tracers from `Registry.xml`

**Files:**
- Modify: `src/core_atmosphere/Registry.xml:1666-1675` (scalars block)
- Modify: `src/core_atmosphere/Registry.xml:2028-2037` (scalars_tend block)

**Step 1: Remove the `#ifdef MPAS_USE_MUSICA` blocks**

Delete lines 1666-1675 (the qAB/qA/qB scalars block including the `#ifdef`/`#endif`):

```xml
<!-- DELETE THIS ENTIRE BLOCK -->
#ifdef MPAS_USE_MUSICA
                        <var name="qAB" array_group="passive" units="kg kg^{-1}"
                             description="Molecular AB mixing ratio"/>

                        <var name="qA" array_group="passive" units="kg kg^{-1}"
                             description="Atomic A mixing ratio"/>

                        <var name="qB" array_group="passive" units="kg kg^{-1}"
                             description="Atomic B mixing ratio"/>
#endif
```

Delete lines 2028-2037 (the tend_qAB/tend_qA/tend_qB block including the `#ifdef`/`#endif`):

```xml
<!-- DELETE THIS ENTIRE BLOCK -->
#ifdef MPAS_USE_MUSICA
                        <var name="tend_qAB" name_in_code="qAB" array_group="passive" units="kg m^{-3} s^{-1}"
                             description="Tendency of molecular AB mass per unit volume divided by d(zeta)/dz"/>

                        <var name="tend_qA" name_in_code="qA" array_group="passive" units="kg m^{-3} s^{-1}"
                             description="Tendency of atomic A mass per unit volume divided by d(zeta)/dz"/>

                        <var name="tend_qB" name_in_code="qB" array_group="passive" units="kg m^{-3} s^{-1}"
                             description="Tendency of atomic B mass per unit volume divided by d(zeta)/dz"/>
#endif
```

**Step 2: Remove the consistency check from `musica_init`**

In `mpas_musica.F`, the `check_registry_tracer_consistency` call (lines 118-124) checks that MICM species match registry tracers. With dynamic allocation, registry tracers won't exist yet at `musica_init` time. Remove or convert to an informational log. Also remove the `registry_tracer_names` parameter from `musica_init` if no longer needed, and update the caller in `mpas_atm_chemistry.F`.

**Step 3: Full rebuild and verify**

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"
make clean CORE=atmosphere && find . -name "*.mod" -delete && find . -name "*.o" -delete
make -j8 llvm CORE=atmosphere PIO=$HOME/software NETCDF=/opt/homebrew PNETCDF=$HOME/software PRECISION=double MUSICA=true 2>&1 | tail -30
```

Expected: Clean build. No references to qAB/qA/qB in generated code.

**Step 4: Commit**

```bash
git add src/core_atmosphere/Registry.xml src/core_atmosphere/chemistry/musica/mpas_musica.F src/core_atmosphere/chemistry/mpas_atm_chemistry.F
git commit -m "Remove static chemistry tracers from Registry.xml, rely on runtime injection"
```

---

## Task 4: Run regression test

**Files:** None (runtime test only)

**Step 1: Run the supercell test**

```bash
cd ~/Data/MPAS/supercell
cp ~/EarthSystem/CheMPAS/atmosphere_model .
mpiexec -n 8 ./atmosphere_model 2>&1 | tee log_phase2.txt
```

**Step 2: Check logs for dynamic tracer allocation**

```bash
grep '\[Chemistry\]' log_phase2.txt
grep '\[MUSICA\]' log_phase2.txt
```

Expected output should include:
- `[MUSICA] Discovered chemistry tracer: qAB`
- `[MUSICA] Discovered chemistry tracer: qA`
- `[MUSICA] Discovered chemistry tracer: qB`
- `[Chemistry] Extending scalars with 3 chemistry tracers`
- `[Chemistry] scalars: old_count=N new_count=N+3`
- `[Chemistry] Added index_qAB = N+1`
- `[Chemistry] Added index_qA = N+2`
- `[Chemistry] Added index_qB = N+3`
- `[MUSICA] Resolved AB -> MPAS index N+1` (etc.)

**Step 3: Compare output against reference**

```bash
cd ~/Data/MPAS/supercell
~/miniconda3/envs/mpas/bin/python -c "
import netCDF4 as nc
import numpy as np

ref = nc.Dataset('reference_phase1_output.nc')
new = nc.Dataset('output.nc')

# Check chemistry tracers exist and match
for var in ['qv', 'theta']:
    r = ref.variables[var][:]
    n = new.variables[var][:]
    diff = np.max(np.abs(r - n))
    print(f'{var}: max_diff = {diff}')

ref.close()
new.close()
"
```

Note: Chemistry tracers (qAB, qA, qB) may not appear in output.nc since I/O registration is deferred. The key check is that met fields are unchanged and the run completes without errors.

**Step 4: Verify chemistry log output matches reference**

```bash
# Compare MUSICA species resolution and chemistry diagnostics
grep 'MUSICA.*Resolved' log_phase2.txt
grep 'column_comparison' log_phase2.txt | head -5
```

**Step 5: Commit test results to PLAN.md**

Update PLAN.md with verification results and commit.

---

## Task 5: Cleanup and finalize

**Step 1: Remove `check_registry_tracer_consistency` subroutine entirely**

If it was only converted to a log in Task 3, now fully remove the dead code from `mpas_musica.F` (lines 547-619).

**Step 2: Update PLAN.md with Phase 2 completion status**

**Step 3: Final commit**

```bash
git add -A
git commit -m "Complete Phase 2: runtime chemistry tracer allocation"
git push origin develop
```
