! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_atm_chemistry
!
!> \brief Manages interactions with chemistry packages
!> \author CheMPAS-A Developers
!> \date   20 August 2025
!> \details
!>  This module manages the interactions with chemistry packages,
!>  including initialization, time-stepping, and finalization.
!>  It provides a framework for integrating various chemistry models
!>  into the CheMPAS-A system.
!
!-------------------------------------------------------------------------
module mpas_atm_chemistry

    use mpas_kind_types, only : StrKIND, RKIND

    implicit none

    private

    public :: chemistry_init, chemistry_step, chemistry_finalize, chemistry_seed_chem
    public :: chemistry_from_MPAS, chemistry_to_MPAS

    contains

    !------------------------------------------------------------------------
    !  routine chemistry_init
    !
    !> \brief Initializes the chemistry packages
    !> \author CheMPAS-A Developers
    !> \date   20 August 2025
    !> \details
    !>  This routine initializes the chemistry packages, setting up
    !>  necessary parameters and data structures for the simulation.
    !------------------------------------------------------------------------
    subroutine chemistry_init(configs, dimensions, structs)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_init, resolve_mpas_indices
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT, field3dReal
        use mpas_pool_routines, only: mpas_pool_get_config, mpas_pool_get_dimension, &
                                      mpas_pool_get_subpool, mpas_pool_get_field, mpas_pool_get_array

        type (mpas_pool_type), intent(in) :: configs
        type (mpas_pool_type), intent(in) :: dimensions
        type (mpas_pool_type), intent(in) :: structs

#ifdef MPAS_USE_MUSICA
        character(len=StrKIND), pointer :: filepath_ptr
        integer                         :: error_code
        character(len=:), allocatable   :: error_message
        integer                         :: nCells
        integer                         :: nVertLevels
        integer, pointer                :: nCells_ptr
        integer, pointer                :: nVertLevels_ptr
        type (mpas_pool_type), pointer  :: statePool
        type (field3dReal), pointer     :: scalarsField
        character(len=StrKIND), dimension(:), pointer :: tracer_names
#endif

        call mpas_log_write('Initializing chemistry packages...')

#ifdef MPAS_USE_MUSICA
        call mpas_pool_get_dimension(dimensions, 'nVertLevels', nVertLevels_ptr)
        nVertLevels = nVertLevels_ptr

        call mpas_pool_get_dimension(dimensions, 'nCells', nCells_ptr)
        if (associated(nCells_ptr)) then
            nCells = nCells_ptr
        else
            nCells = 0
        end if

        call mpas_pool_get_config(configs, 'config_micm_file', filepath_ptr)

        nullify(statePool)
        nullify(scalarsField)
        nullify(tracer_names)
        nullify(tracer_names)

        call mpas_pool_get_subpool(structs, 'state', statePool)
        if (associated(statePool)) then
            call mpas_pool_get_field(statePool, 'scalars', scalarsField, 1)
            if (associated(scalarsField)) then
                tracer_names => scalarsField%constituentNames
            end if
        end if

        ! Initialize MICM with total grid cells = nCells × nVertLevels
        ! This allows processing all cells in parallel rather than column-by-column
        if (associated(tracer_names)) then
            call musica_init(filepath_ptr, nCells * nVertLevels, error_code, error_message, tracer_names)
        else
            call musica_init(filepath_ptr, nCells * nVertLevels, error_code, error_message)
        end if

        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if

        ! Resolve MPAS tracer indices for all chemistry species
        if (associated(statePool)) then
            call resolve_mpas_indices(statePool, error_code, error_message)
            if (error_code /= 0) then
                call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            end if
        end if

        ! Seed MPAS tracers from MICM initial state (both time levels)
        call chemistry_seed_chem(dimensions, structs)
#endif

    end subroutine chemistry_init


    !------------------------------------------------------------------------
    !  routine chemistry_step
    !
    !> \brief Steps the chemistry packages
    !> \author CheMPAS-A Developers
    !> \date   20 August 2025
    !> \details
    !>  This routine steps the chemistry packages, updating their state
    !>  based on the current simulation time and conditions.
    !>  Now includes bidirectional coupling: MPAS -> MICM -> MPAS
    !------------------------------------------------------------------------
    subroutine chemistry_step(dt, mesh, state, diag, dimensions, time_lev)

#ifdef MPAS_USE_MUSICA
        use iso_fortran_env, only: real64
        use mpas_musica, only: musica_step, musica_step_ref, log_column_comparison, &
                               copy_state_to_ref
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_dimension

        real (kind=RKIND), intent(in) :: dt
        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(inout) :: state
        type (mpas_pool_type), intent(in) :: diag
        type (mpas_pool_type), intent(in) :: dimensions
        integer, intent(in) :: time_lev

#ifdef MPAS_USE_MUSICA
        real(real64)                  :: time_step
        integer                       :: error_code
        character(len=:), allocatable :: error_message
        integer                       :: nCells, nVertLevels
        integer, pointer              :: nCells_ptr, nVertLevels_ptr
        integer                       :: probe_cell
        logical, save                 :: first_call = .true.

        time_step = dt

        ! Get dimensions
        call mpas_pool_get_dimension(dimensions, 'nCells', nCells_ptr)
        call mpas_pool_get_dimension(dimensions, 'nVertLevels', nVertLevels_ptr)
        nCells = nCells_ptr
        nVertLevels = nVertLevels_ptr

        ! Probe cell near center of domain for advection diagnostics
        probe_cell = nCells / 2

        call mpas_log_write('Stepping chemistry packages...')

        ! Step 1: Copy MPAS state to MICM (MPAS -> MICM)
        call chemistry_from_MPAS(mesh, state, diag, dimensions, time_lev, &
                                  nCells, nVertLevels, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            return
        end if

        ! On first call, copy coupled state to reference so both start with
        ! identical initial conditions from MPAS. Any subsequent divergence
        ! will indicate advection effects.
        if (first_call) then
            call copy_state_to_ref()
            call mpas_log_write('[Chemistry] Copied coupled state to reference state (first timestep).')
            first_call = .false.
        end if

        ! Step 2: Solve chemistry (coupled state)
        call musica_step(time_step, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            return
        end if

        ! Step 2b: Solve chemistry on reference state (uncoupled, for advection diagnostics)
        call musica_step_ref(time_step, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if

        ! Log comparison of coupled vs reference at probe cell
        call log_column_comparison(probe_cell, nVertLevels)

        ! Step 3: Copy MICM state back to MPAS (MICM -> MPAS)
        call chemistry_to_MPAS(mesh, state, diag, dimensions, time_lev, &
                                nCells, nVertLevels, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if
#endif

    end subroutine chemistry_step


    !------------------------------------------------------------------------
    !  routine chemistry_finalize
    !
    !> \brief Finalizes the chemistry packages
    !> \author CheMPAS-A Developers
    !> \date   20 August 2025
    !> \details
    !>  This routine finalizes the chemistry packages, cleaning up
    !>  any resources and data structures used during the simulation.
    !------------------------------------------------------------------------
    subroutine chemistry_finalize()

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_finalize
#endif
        use mpas_log, only : mpas_log_write

#ifdef MPAS_USE_MUSICA
        call mpas_log_write('Finalizing chemistry packages...')
        call musica_finalize()
#endif

    end subroutine chemistry_finalize


    !------------------------------------------------------------------------
    !  routine chemistry_from_MPAS
    !
    !> \brief Extracts MPAS state for chemistry processing
    !> \author CheMPAS-A Developers
    !> \date   January 2026
    !> \details
    !>  This routine extracts MPAS state variables (scalars, thermodynamic
    !>  fields) and passes them to the chemistry package (MICM) for processing.
    !>  Following the physics interface pattern from mpas_atmphys_interface.F
    !------------------------------------------------------------------------
    subroutine chemistry_from_MPAS(mesh, state, diag, dimensions, time_lev, &
                                    nCells, nVertLevels, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: MICM_from_chemistry
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_array, mpas_pool_get_dimension

        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(in) :: state
        type (mpas_pool_type), intent(in) :: diag
        type (mpas_pool_type), intent(in) :: dimensions
        integer, intent(in) :: time_lev
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        ! MPAS state arrays
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz, theta_m
        real (kind=RKIND), dimension(:,:), pointer :: zz, exner
        real (kind=RKIND), dimension(:,:), pointer :: pressure_p, pressure_b
        integer, pointer :: index_qv

        ! Local work arrays
        real (kind=RKIND), allocatable :: rho_air(:,:)
        real (kind=RKIND), allocatable :: temperature(:,:)
        real (kind=RKIND), allocatable :: pressure(:,:)

        real (kind=RKIND), parameter :: rvord = 461.5_RKIND / 287.0_RKIND  ! Rv/Rd

        integer :: k, iCell

        error_code = 0
        error_message = ''

        ! Get MPAS arrays from pools
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
        call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
        call mpas_pool_get_array(mesh, 'zz', zz)
        call mpas_pool_get_array(diag, 'exner', exner)
        call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)

        ! Get water vapor index (needed for thermodynamic calculations)
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)

        ! Validate pointers
        if (.not. associated(scalars) .or. .not. associated(rho_zz) .or. &
            .not. associated(theta_m) .or. .not. associated(zz) .or. &
            .not. associated(exner)) then
            error_code = 1
            error_message = '[Chemistry] Required MPAS arrays not available in chemistry_from_MPAS.'
            return
        end if

        ! Allocate work arrays
        allocate(rho_air(nVertLevels, nCells))
        allocate(temperature(nVertLevels, nCells))
        allocate(pressure(nVertLevels, nCells))

        ! Calculate derived quantities following physics interface pattern
        do iCell = 1, nCells
            do k = 1, nVertLevels
                ! Moist air density [kg/m³]: rho = zz * rho_zz * (1 + qv)
                rho_air(k, iCell) = zz(k, iCell) * rho_zz(k, iCell) * &
                                    (1.0_RKIND + scalars(index_qv, k, iCell))

                ! Temperature [K]: T = theta * exner, where theta = theta_m / (1 + Rv/Rd * qv)
                temperature(k, iCell) = (theta_m(k, iCell) / &
                    (1.0_RKIND + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell)

                ! Total pressure [Pa]
                if (associated(pressure_p) .and. associated(pressure_b)) then
                    pressure(k, iCell) = pressure_p(k, iCell) + pressure_b(k, iCell)
                else
                    ! Fallback: estimate from ideal gas law
                    pressure(k, iCell) = rho_air(k, iCell) * 287.0_RKIND * temperature(k, iCell)
                end if
            end do
        end do

        ! Call MICM-specific coupling routine (uses chem_species table internally)
        call MICM_from_chemistry(scalars, rho_air, temperature, pressure, &
                                  nCells, nVertLevels, &
                                  error_code, error_message)

        deallocate(rho_air)
        deallocate(temperature)
        deallocate(pressure)
#else
        error_code = 0
        error_message = ''
#endif

    end subroutine chemistry_from_MPAS


    !------------------------------------------------------------------------
    !  routine chemistry_to_MPAS
    !
    !> \brief Updates MPAS state from chemistry results
    !> \author CheMPAS-A Developers
    !> \date   January 2026
    !> \details
    !>  This routine copies the updated chemistry state back to MPAS scalars
    !>  after the chemistry solver has completed.
    !------------------------------------------------------------------------
    subroutine chemistry_to_MPAS(mesh, state, diag, dimensions, time_lev, &
                                  nCells, nVertLevels, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: MICM_to_chemistry
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_array, mpas_pool_get_dimension

        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(inout) :: state
        type (mpas_pool_type), intent(in) :: diag
        type (mpas_pool_type), intent(in) :: dimensions
        integer, intent(in) :: time_lev
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        ! MPAS state arrays
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
        integer, pointer :: index_qv

        ! Local work arrays
        real (kind=RKIND), allocatable :: rho_air(:,:)

        integer :: k, iCell

        error_code = 0
        error_message = ''

        ! Get MPAS arrays from pools
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
        call mpas_pool_get_array(mesh, 'zz', zz)

        ! Get water vapor index
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)

        ! Validate pointers
        if (.not. associated(scalars) .or. .not. associated(rho_zz) .or. &
            .not. associated(zz)) then
            error_code = 1
            error_message = '[Chemistry] Required MPAS arrays not available in chemistry_to_MPAS.'
            return
        end if

        ! Allocate work arrays
        allocate(rho_air(nVertLevels, nCells))

        ! Calculate air density for unit conversion
        do iCell = 1, nCells
            do k = 1, nVertLevels
                rho_air(k, iCell) = zz(k, iCell) * rho_zz(k, iCell) * &
                                    (1.0_RKIND + scalars(index_qv, k, iCell))
            end do
        end do

        ! Call MICM-specific coupling routine (uses chem_species table internally)
        call MICM_to_chemistry(scalars, rho_air, nCells, nVertLevels, &
                                error_code, error_message)

        deallocate(rho_air)
#else
        error_code = 0
        error_message = ''
#endif

    end subroutine chemistry_to_MPAS


    !------------------------------------------------------------------------
    !  routine chemistry_seed_chem
    !
    !> \brief One-time (or repeated) MICM -> MPAS chemistry tracer seed
    !------------------------------------------------------------------------
    subroutine chemistry_seed_chem(dimensions, structs)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: micm_to_mpas_chem, n_chem_species, chem_species
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT, field3dReal
        use mpas_pool_routines, only: mpas_pool_get_dimension, mpas_pool_get_subpool, &
                                      mpas_pool_get_field, mpas_pool_get_array

        type (mpas_pool_type), intent(in) :: dimensions
        type (mpas_pool_type), intent(in) :: structs

#ifdef MPAS_USE_MUSICA
        integer                         :: nCells
        integer                         :: nVertLevels
        integer, pointer                :: nCells_ptr
        integer, pointer                :: nVertLevels_ptr
        type (mpas_pool_type), pointer  :: statePool
        type (field3dReal), pointer     :: scalarsField
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars_time1
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars_time2
        integer                         :: error_code, i_sp, mpas_idx
        character(len=:), allocatable   :: error_message

        call mpas_log_write('[MUSICA] Seeding chemistry tracers from MICM into MPAS scalars (both time levels).')

        if (n_chem_species == 0) then
            call mpas_log_write('[MUSICA] No chemistry species to seed.')
            return
        end if

        call mpas_pool_get_dimension(dimensions, 'nVertLevels', nVertLevels_ptr)
        nVertLevels = nVertLevels_ptr

        call mpas_pool_get_dimension(dimensions, 'nCells', nCells_ptr)
        if (associated(nCells_ptr)) then
            nCells = nCells_ptr
        else
            nCells = 0
        end if

        nullify(statePool)
        nullify(scalarsField)
        nullify(scalars_time1)
        nullify(scalars_time2)

        call mpas_pool_get_subpool(structs, 'state', statePool)
        if (associated(statePool)) then
            call mpas_pool_get_field(statePool, 'scalars', scalarsField, 1)
            call mpas_pool_get_array(statePool, 'scalars', scalars_time1, 1)
            call mpas_pool_get_array(statePool, 'scalars', scalars_time2, 2)
        end if

        if (.not. associated(scalars_time1)) then
            call mpas_log_write('[MUSICA] Skipping MICM->MPAS seed: scalars time level 1 not available.')
            return
        else if (.not. associated(scalars_time2)) then
            call mpas_log_write('[MUSICA] Skipping MICM->MPAS seed: scalars time level 2 not available.')
            return
        end if

        ! Check if any chemistry tracer has MPAS index not yet resolved
        do i_sp = 1, n_chem_species
            if (chem_species(i_sp)%mpas_index < 1) then
                call mpas_log_write('[MUSICA] Skipping MICM->MPAS seed: MPAS indices not resolved.')
                return
            end if
        end do

        ! Skip seeding if first chemistry tracer already has spatial gradients (initialized from file)
        mpas_idx = chem_species(1)%mpas_index
        if (maxval(scalars_time1(mpas_idx,:,:)) - minval(scalars_time1(mpas_idx,:,:)) > 1.0e-10_RKIND) then
            call mpas_log_write('[MUSICA] Tracers already initialized with gradients, skipping MICM seed.')
            return
        end if

        error_code = 0
        call micm_to_mpas_chem(scalars_time1, nCells, nVertLevels, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if

        error_code = 0
        call micm_to_mpas_chem(scalars_time2, nCells, nVertLevels, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if
#endif

    end subroutine chemistry_seed_chem

end module mpas_atm_chemistry
