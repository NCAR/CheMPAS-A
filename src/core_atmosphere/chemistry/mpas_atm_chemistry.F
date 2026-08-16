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
    use mpas_chem_constants, only : pii, rgas, rvord, CHEM_M_AIR
    use mpas_timer, only : mpas_timer_start, mpas_timer_stop

    implicit none

    private

    public :: chemistry_init, chemistry_step, chemistry_finalize, chemistry_seed_chem
    public :: chemistry_from_MPAS, chemistry_to_MPAS
    public :: chemistry_query_species, chemistry_query_photo_rates
    public :: chemistry_query_emission_species
    public :: chemistry_prepare_miem_diagnostic_selection
    public :: chemistry_prepare_tracer_outputs, chemistry_clear_tracer_outputs
    public :: chemistry_add_tracer_output_vars
    public :: chemistry_configure_tracer_output_streams
    public :: chemistry_compute_tracer_outputs
    public :: chem_env_t, chem_env_fill

    ! Lightning NOx source (no-op if qNO not in mechanism)
    public :: chemistry_lightning_nox_init

    type :: chem_env_t
        real(kind=RKIND), allocatable :: temperature(:,:)  ! (nVertLevels, nCellsSolve) [K]
        real(kind=RKIND), allocatable :: pressure(:,:)     ! [Pa]
        real(kind=RKIND), allocatable :: rho_dry(:,:)      ! [kg m-3]
        real(kind=RKIND), dimension(:,:), pointer :: zgrid => null()
    end type chem_env_t

    integer, parameter :: CHEM_OUTPUT_MMR      = 0
    integer, parameter :: CHEM_OUTPUT_FRACTION = 1
    integer, parameter :: CHEM_OUTPUT_PERCENT  = 2
    integer, parameter :: CHEM_OUTPUT_PPMV     = 3
    integer, parameter :: CHEM_OUTPUT_PPBV     = 4
    integer, parameter :: CHEM_OUTPUT_PPTV     = 5

    type :: chem_output_descriptor_t
        character(len=64)      :: micm_name = ''
        character(len=64)      :: mpas_name = ''
        character(len=64)      :: output_name = ''
        character(len=64)      :: units = '-'
        character(len=StrKIND) :: description = ''
        real(kind=RKIND)       :: molar_mass = 0.0_RKIND
        real(kind=RKIND)       :: unit_scale = 1.0_RKIND
        real(kind=RKIND)       :: fixed_vmr = -1.0_RKIND
        integer                :: mpas_index = -1
    end type chem_output_descriptor_t

    integer, save :: n_chem_tracer_outputs = 0
    type(chem_output_descriptor_t), allocatable, save :: chem_tracer_outputs(:)

    ! Solar geometry config values (cached at init, used per-step)
    real(kind=RKIND), save :: chem_lat = 0.0_RKIND
    real(kind=RKIND), save :: chem_lon = 0.0_RKIND
    real(kind=RKIND), save :: chem_j_no2_max = 0.0_RKIND
    logical, save :: chem_use_grid_coords = .false.
    integer, save :: chem_consecutive_failures = 0
    integer, save :: chem_max_consecutive_failures = 3

    ! TUV-x update interval state (set in chemistry_init, accumulated in chemistry_step)
    real(kind=RKIND), save :: tuvx_update_interval = 0.0_RKIND
    real(kind=RKIND), save :: tuvx_time_since_last = huge(1.0_RKIND)

    ! Chemistry solve interval state (0 = every MPAS step). The clock carries
    ! cadence remainder; accumulated_dt is elapsed time since the last solve.
    real(kind=RKIND), save :: chemistry_interval = 0.0_RKIND
    real(kind=RKIND), save :: chemistry_interval_clock = 0.0_RKIND
    real(kind=RKIND), save :: chemistry_accumulated_dt = 0.0_RKIND

        ! Reference solve (doubles chemistry cost; off by default)
        logical, save :: use_ref_solve = .false.

        ! Number of MICM sub-steps per MPAS dt (1 = no sub-stepping)
        integer, save :: chem_substeps_val = 1

        ! MICM state relative tolerance (MICM default is 1e-6)
        real(kind=RKIND), save :: micm_relative_tolerance_val = 1.0e-6_RKIND

        ! TUV-x state
        logical, save :: use_tuvx = .false.
        integer, save :: idx_qO3 = -1
        integer, save :: idx_qc  = -1
        integer, save :: idx_qr  = -1

    contains

    !------------------------------------------------------------------------
    !  routine chem_fatal
    !
    !> \brief Log an init-time chemistry fatal error.
    !------------------------------------------------------------------------
    subroutine chem_fatal(msg)

        use mpas_log, only: mpas_log_write
        use mpas_derived_types, only: MPAS_LOG_CRIT

        character(len=*), intent(in) :: msg

        call mpas_log_write('[Chemistry] ' // trim(msg), messageType=MPAS_LOG_CRIT)

    end subroutine chem_fatal


    pure function chem_lowercase(value) result(lower)

        character(len=*), intent(in) :: value
        character(len=len(value)) :: lower
        integer :: i, code

        lower = value
        do i = 1, len(value)
            code = iachar(value(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do

    end function chem_lowercase


    subroutine parse_chem_output_unit(value, unit_kind, valid)

        character(len=*), intent(in) :: value
        integer, intent(out) :: unit_kind
        logical, intent(out) :: valid

        valid = .true.
        select case (trim(chem_lowercase(adjustl(value))))
        case ('mmr')
            unit_kind = CHEM_OUTPUT_MMR
        case ('fraction')
            unit_kind = CHEM_OUTPUT_FRACTION
        case ('percent')
            unit_kind = CHEM_OUTPUT_PERCENT
        case ('ppmv')
            unit_kind = CHEM_OUTPUT_PPMV
        case ('ppbv')
            unit_kind = CHEM_OUTPUT_PPBV
        case ('pptv')
            unit_kind = CHEM_OUTPUT_PPTV
        case default
            unit_kind = CHEM_OUTPUT_MMR
            valid = .false.
        end select

    end subroutine parse_chem_output_unit


    subroutine set_chem_output_properties(unit_kind, scale, units)

        integer, intent(in) :: unit_kind
        real(kind=RKIND), intent(out) :: scale
        character(len=*), intent(out) :: units

        select case (unit_kind)
        case (CHEM_OUTPUT_FRACTION)
            scale = 1.0_RKIND
            units = 'mol mol^{-1}'
        case (CHEM_OUTPUT_PERCENT)
            scale = 1.0e2_RKIND
            units = 'percent'
        case (CHEM_OUTPUT_PPMV)
            scale = 1.0e6_RKIND
            units = 'ppmv'
        case (CHEM_OUTPUT_PPBV)
            scale = 1.0e9_RKIND
            units = 'ppbv'
        case (CHEM_OUTPUT_PPTV)
            scale = 1.0e12_RKIND
            units = 'pptv'
        case default
            scale = 1.0_RKIND
            units = 'kg kg^{-1}'
        end select

    end subroutine set_chem_output_properties


    subroutine chemistry_clear_tracer_outputs()

        n_chem_tracer_outputs = 0
        if (allocated(chem_tracer_outputs)) deallocate(chem_tracer_outputs)

    end subroutine chemistry_clear_tracer_outputs


    !------------------------------------------------------------------------
    !  routine chemistry_prepare_tracer_outputs
    !
    !> \brief Parse tracer output units and cache output-only VMR descriptors.
    !------------------------------------------------------------------------
    subroutine chemistry_prepare_tracer_outputs(global_unit, overrides, num_species, &
            micm_names, mpas_names, molar_masses, fixed_vmrs, error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_log, only: mpas_log_write

        character(len=*), intent(in) :: global_unit
        character(len=*), intent(in) :: overrides
        integer, intent(in) :: num_species
        character(len=64), intent(in) :: micm_names(:)
        character(len=64), intent(in) :: mpas_names(:)
        real(real64), intent(in) :: molar_masses(:)
        real(real64), intent(in) :: fixed_vmrs(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer, allocatable :: selected_units(:)
        logical, allocatable :: overridden(:)
        character(len=StrKIND) :: work, token, species_name, unit_name
        character(len=64) :: output_name
        character(len=64) :: output_units
        integer :: global_kind, override_kind
        integer :: i, j, start_pos, stop_pos, comma_pos, equal_pos, match_idx
        integer :: output_idx
        logical :: valid
        real(kind=RKIND) :: output_scale

        error_code = 0
        error_message = ''
        call chemistry_clear_tracer_outputs()

        if (size(micm_names) /= num_species .or. size(mpas_names) /= num_species .or. &
            size(molar_masses) /= num_species .or. size(fixed_vmrs) /= num_species) then
            error_code = 1
            error_message = 'Internal species metadata size mismatch while preparing tracer output units.'
            return
        end if

        call parse_chem_output_unit(global_unit, global_kind, valid)
        if (.not. valid) then
            error_code = 1
            error_message = 'Invalid config_chem_tracer_output_unit=''' // &
                trim(global_unit) // '''. Expected mmr, fraction, percent, ppmv, ppbv, or pptv.'
            return
        end if

        allocate(selected_units(num_species))
        allocate(overridden(num_species))
        selected_units = global_kind
        overridden = .false.

        work = trim(overrides)
        if (len_trim(work) > 0) then
            start_pos = 1
            do
                comma_pos = index(work(start_pos:len_trim(work)), ',')
                if (comma_pos == 0) then
                    stop_pos = len_trim(work)
                else
                    stop_pos = start_pos + comma_pos - 2
                end if

                if (stop_pos < start_pos) then
                    error_code = 1
                    error_message = 'Empty entry in config_chem_tracer_output_overrides.'
                    return
                end if

                token = adjustl(work(start_pos:stop_pos))
                equal_pos = index(token, '=')
                if (equal_pos <= 1 .or. equal_pos >= len_trim(token)) then
                    error_code = 1
                    error_message = 'Malformed chemistry tracer output override ''' // &
                        trim(token) // '''. Expected SPECIES=UNIT.'
                    return
                end if
                if (index(token(equal_pos + 1:len_trim(token)), '=') > 0) then
                    error_code = 1
                    error_message = 'Malformed chemistry tracer output override ''' // &
                        trim(token) // ''': too many equals signs.'
                    return
                end if

                species_name = trim(adjustl(token(1:equal_pos - 1)))
                unit_name = trim(adjustl(token(equal_pos + 1:len_trim(token))))
                match_idx = 0
                do i = 1, num_species
                    if (trim(micm_names(i)) == trim(species_name)) then
                        match_idx = i
                        exit
                    end if
                end do
                if (match_idx == 0) then
                    error_code = 1
                    error_message = 'Unknown MICM species ''' // trim(species_name) // &
                        ''' in config_chem_tracer_output_overrides.'
                    return
                end if
                if (overridden(match_idx)) then
                    error_code = 1
                    error_message = 'Duplicate output-unit override for MICM species ''' // &
                        trim(species_name) // '''.'
                    return
                end if

                call parse_chem_output_unit(unit_name, override_kind, valid)
                if (.not. valid) then
                    error_code = 1
                    error_message = 'Invalid output unit ''' // trim(unit_name) // &
                        ''' for MICM species ''' // trim(species_name) // '''.'
                    return
                end if
                selected_units(match_idx) = override_kind
                overridden(match_idx) = .true.

                if (comma_pos == 0) exit
                start_pos = stop_pos + 2
                if (start_pos > len_trim(work)) then
                    error_code = 1
                    error_message = 'Trailing comma in config_chem_tracer_output_overrides.'
                    return
                end if
            end do
        end if

        n_chem_tracer_outputs = count(selected_units /= CHEM_OUTPUT_MMR)
        do i = 1, num_species
            if (selected_units(i) == CHEM_OUTPUT_MMR) cycle
            if (molar_masses(i) <= 0.0_real64) then
                error_code = 1
                error_message = 'MICM species ''' // trim(micm_names(i)) // &
                    ''' has a non-positive molar mass; cannot produce VMR output.'
                call chemistry_clear_tracer_outputs()
                return
            end if
            if (len_trim(micm_names(i)) + len('vmr_') > len(output_name)) then
                error_code = 1
                error_message = 'VMR output name for MICM species ''' // trim(micm_names(i)) // &
                    ''' exceeds the MPAS 64-character field-name limit.'
                call chemistry_clear_tracer_outputs()
                return
            end if
        end do

        if (n_chem_tracer_outputs == 0) then
            call mpas_log_write('[Chemistry] MICM tracer history output remains in mass mixing ratio units')
            return
        end if

        allocate(chem_tracer_outputs(n_chem_tracer_outputs))
        output_idx = 0
        do i = 1, num_species
            if (selected_units(i) == CHEM_OUTPUT_MMR) cycle
            output_idx = output_idx + 1
            output_name = 'vmr_' // trim(micm_names(i))
            call set_chem_output_properties(selected_units(i), output_scale, output_units)

            do j = 1, output_idx - 1
                if (trim(chem_tracer_outputs(j) % output_name) == trim(output_name)) then
                    error_code = 1
                    error_message = 'Duplicate VMR output field name ''' // trim(output_name) // '''.'
                    call chemistry_clear_tracer_outputs()
                    return
                end if
            end do

            chem_tracer_outputs(output_idx) % micm_name = trim(micm_names(i))
            chem_tracer_outputs(output_idx) % mpas_name = trim(mpas_names(i))
            chem_tracer_outputs(output_idx) % output_name = trim(output_name)
            chem_tracer_outputs(output_idx) % molar_mass = real(molar_masses(i), RKIND)
            chem_tracer_outputs(output_idx) % unit_scale = output_scale
            chem_tracer_outputs(output_idx) % fixed_vmr = real(fixed_vmrs(i), RKIND)
            chem_tracer_outputs(output_idx) % units = trim(output_units)
            chem_tracer_outputs(output_idx) % description = 'Dry-air volume mixing ratio of MICM species ' // &
                trim(micm_names(i)) // ' in ' // trim(output_units)

            if (fixed_vmrs(i) >= 0.0_real64) then
                call mpas_log_write('[Chemistry] Configured ' // trim(output_name) // &
                    ' from its fixed host dry-air mole fraction in ' // trim(output_units))
            else
                call mpas_log_write('[Chemistry] Configured ' // trim(output_name) // &
                    ' from ' // trim(mpas_names(i)) // ' in ' // trim(output_units))
            end if
        end do

    end subroutine chemistry_prepare_tracer_outputs


    subroutine chemistry_resolve_tracer_output_indices(state_pool, error_code, error_message)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_dimension

        type(mpas_pool_type), intent(in) :: state_pool
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: i
        integer, pointer :: index_ptr

        error_code = 0
        error_message = ''
        do i = 1, n_chem_tracer_outputs
            if (chem_tracer_outputs(i) % fixed_vmr >= 0.0_RKIND) cycle
            nullify(index_ptr)
            call mpas_pool_get_dimension(state_pool, &
                'index_' // trim(chem_tracer_outputs(i) % mpas_name), index_ptr)
            if (.not. associated(index_ptr) .or. index_ptr < 1) then
                error_code = 1
                error_message = 'No valid MPAS scalar index for VMR source field ''' // &
                    trim(chem_tracer_outputs(i) % mpas_name) // '''.'
                return
            end if
            chem_tracer_outputs(i) % mpas_index = index_ptr
        end do

    end subroutine chemistry_resolve_tracer_output_indices


    !------------------------------------------------------------------------
    !  routine chemistry_add_tracer_output_vars
    !
    !> \brief Add configured output-only VMR diagnostics to MPAS pools.
    !------------------------------------------------------------------------
    subroutine chemistry_add_tracer_output_vars(block, diag_pool, package_pool)

        use mpas_derived_types, only: block_type, mpas_pool_type, MPAS_LOG_CRIT
        use mpas_runtime_vars, only: MPAS_var, MPAS_VAR_REAL, MPAS_VAR_OK
        use mpas_log, only: mpas_log_write

        type(block_type), target, intent(inout) :: block
        type(mpas_pool_type), intent(inout) :: diag_pool
        type(mpas_pool_type), intent(in) :: package_pool

        type(MPAS_var) :: new_var
        integer :: i, add_err

        do i = 1, n_chem_tracer_outputs
            new_var % name = trim(chem_tracer_outputs(i) % output_name)
            new_var % hasTimeDimension = .true.
            new_var % vartype = MPAS_VAR_REAL
            new_var % dimensions = 'nVertLevels nCells'
            new_var % units = trim(chem_tracer_outputs(i) % units)
            new_var % description = trim(chem_tracer_outputs(i) % description)
            new_var % isPersistent = .true.
            new_var % nTimeLevels = 1

            call new_var % add_to_pool(block, diag_pool, package_pool, &
                useNameInCode=.true., ierr=add_err)
            if (add_err /= MPAS_VAR_OK) then
                call mpas_log_write('Failed to add runtime chemistry output diagnostic ' // &
                    trim(chem_tracer_outputs(i) % output_name) // ' to diag', &
                    messageType=MPAS_LOG_CRIT)
                return
            end if
            call new_var % add_to_pool(block, block % allFields, package_pool, ierr=add_err)
            if (add_err /= MPAS_VAR_OK) then
                call mpas_log_write('Failed to add runtime chemistry output diagnostic ' // &
                    trim(chem_tracer_outputs(i) % output_name) // ' to allFields', &
                    messageType=MPAS_LOG_CRIT)
                return
            end if

            call mpas_log_write('[Chemistry] Added runtime output diagnostic ' // &
                trim(chem_tracer_outputs(i) % output_name) // ' to diag')
            call new_var % reset()
        end do

    end subroutine chemistry_add_tracer_output_vars


    !------------------------------------------------------------------------
    !  routine chemistry_configure_tracer_output_streams
    !
    !> \brief Mirror configured VMR diagnostics into mutable scalar outputs.
    !------------------------------------------------------------------------
    subroutine chemistry_configure_tracer_output_streams(stream_manager, error_code)

        use mpas_derived_types, only: MPAS_streamManager_type, MPAS_STREAM_OUTPUT, &
            MPAS_STREAM_INPUT_OUTPUT, MPAS_STREAM_MGR_NOERR
        use mpas_stream_manager, only: MPAS_stream_mgr_begin_iteration, &
            MPAS_stream_mgr_get_next_stream, MPAS_stream_mgr_get_next_field, &
            MPAS_stream_mgr_add_field
        use mpas_log, only: mpas_log_write

        type(MPAS_streamManager_type), intent(inout) :: stream_manager
        integer, intent(out) :: error_code

        character(len=StrKIND) :: stream_name, field_name
        integer :: stream_direction, add_err, i
        logical :: stream_immutable, field_active, has_scalars
        logical, allocatable :: already_present(:)

        error_code = 0
        if (n_chem_tracer_outputs == 0) return

        allocate(already_present(n_chem_tracer_outputs))
        call MPAS_stream_mgr_begin_iteration(stream_manager)
        do while (MPAS_stream_mgr_get_next_stream(stream_manager, streamID=stream_name, &
                directionProperty=stream_direction, immutableProperty=stream_immutable))
            if (stream_immutable) cycle
            if (stream_direction /= MPAS_STREAM_OUTPUT .and. &
                stream_direction /= MPAS_STREAM_INPUT_OUTPUT) cycle

            has_scalars = .false.
            already_present = .false.
            call MPAS_stream_mgr_begin_iteration(stream_manager, streamID=stream_name)
            do while (MPAS_stream_mgr_get_next_field(stream_manager, stream_name, &
                    field_name, isActive=field_active))
                if (field_active .and. trim(field_name) == 'scalars') has_scalars = .true.
                do i = 1, n_chem_tracer_outputs
                    if (trim(field_name) == trim(chem_tracer_outputs(i) % output_name)) then
                        already_present(i) = .true.
                    end if
                end do
            end do
            if (.not. has_scalars) cycle

            do i = 1, n_chem_tracer_outputs
                if (already_present(i)) cycle
                call MPAS_stream_mgr_add_field(stream_manager, trim(stream_name), &
                    trim(chem_tracer_outputs(i) % output_name), ierr=add_err)
                if (add_err /= MPAS_STREAM_MGR_NOERR) then
                    error_code = add_err
                    return
                end if
                call mpas_log_write('[Chemistry] Added ' // &
                    trim(chem_tracer_outputs(i) % output_name) // &
                    ' to mutable scalar output stream ' // trim(stream_name))
            end do
        end do

    end subroutine chemistry_configure_tracer_output_streams


    !------------------------------------------------------------------------
    !  routine chemistry_compute_tracer_outputs
    !
    !> \brief Convert authoritative dry-air MMR scalars into VMR diagnostics.
    !------------------------------------------------------------------------
    subroutine chemistry_compute_tracer_outputs(state_pool, time_lev, diag_pool)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_array

        type(mpas_pool_type), intent(in) :: state_pool
        integer, intent(in) :: time_lev
        type(mpas_pool_type), intent(inout) :: diag_pool

        real(kind=RKIND), dimension(:,:,:), pointer :: scalars
        real(kind=RKIND), dimension(:,:), pointer :: output_field
        real(kind=RKIND) :: conversion_factor
        integer :: i

        if (n_chem_tracer_outputs == 0) return

        nullify(scalars)
        call mpas_pool_get_array(state_pool, 'scalars', scalars, time_lev)
        if (.not. associated(scalars)) then
            call chem_fatal('State/scalars is unavailable while computing chemistry tracer output diagnostics.')
            return
        end if

        do i = 1, n_chem_tracer_outputs
            nullify(output_field)
            call mpas_pool_get_array(diag_pool, trim(chem_tracer_outputs(i) % output_name), &
                output_field)
            if (.not. associated(output_field)) then
                call chem_fatal('Diagnostic field ' // &
                    trim(chem_tracer_outputs(i) % output_name) // &
                    ' is unavailable during chemistry tracer output conversion.')
                return
            end if

            if (chem_tracer_outputs(i) % fixed_vmr >= 0.0_RKIND) then
                output_field(:,:) = chem_tracer_outputs(i) % fixed_vmr * &
                    chem_tracer_outputs(i) % unit_scale
                cycle
            end if
            if (chem_tracer_outputs(i) % mpas_index < 1 .or. &
                chem_tracer_outputs(i) % mpas_index > size(scalars, 1)) then
                call chem_fatal('Invalid source scalar index for chemistry output diagnostic ' // &
                    trim(chem_tracer_outputs(i) % output_name) // '.')
                return
            end if

            conversion_factor = real(CHEM_M_AIR, RKIND) / &
                chem_tracer_outputs(i) % molar_mass * chem_tracer_outputs(i) % unit_scale
            output_field(:,:) = scalars(chem_tracer_outputs(i) % mpas_index,:,:) * &
                conversion_factor
        end do

    end subroutine chemistry_compute_tracer_outputs


    !------------------------------------------------------------------------
    !  routine chem_step_error
    !
    !> \brief Log a step-time chemistry error and count consecutive failures.
    !------------------------------------------------------------------------
    subroutine chem_step_error(msg, nfail)

        use mpas_log, only: mpas_log_write
        use mpas_derived_types, only: MPAS_LOG_ERR, MPAS_LOG_CRIT

        character(len=*), intent(in) :: msg
        integer, intent(inout) :: nfail

        nfail = nfail + 1
        if (nfail >= chem_max_consecutive_failures) then
            call mpas_log_write('[Chemistry] ' // trim(msg) // ' ($i consecutive failures)', &
                intArgs=[nfail], messageType=MPAS_LOG_CRIT)
        else
            call mpas_log_write('[Chemistry] ' // trim(msg) // &
                ' - keeping previous chemistry state', messageType=MPAS_LOG_ERR)
        end if

    end subroutine chem_step_error


#ifdef MPAS_USE_MUSICA
    !------------------------------------------------------------------------
    !  routine chemistry_miem_epoch
    !
    !> \brief Convert the chemistry interval start to Unix epoch seconds
    !> \details
    !>  Uses MPAS timekeeping for Gregorian calendar arithmetic. For an
    !>  accumulated chemistry interval, curr_time is rewound by
    !>  ``time_step - model_dt`` so MIEM samples the interval start.
    !------------------------------------------------------------------------
    subroutine chemistry_miem_epoch(curr_time, time_step, model_dt, epoch_seconds, &
                                    error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use iso_fortran_env, only: real64
        use mpas_timekeeping, only: MPAS_Time_type, MPAS_TimeInterval_type, &
                                    mpas_set_time, mpas_set_timeInterval, &
                                    mpas_get_timeInterval, operator(-)

        type(MPAS_Time_type), intent(in) :: curr_time
        real(real64), intent(in) :: time_step
        real(kind=RKIND), intent(in) :: model_dt
        real(real64), intent(out) :: epoch_seconds
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(MPAS_Time_type) :: unix_epoch
        type(MPAS_Time_type) :: interval_start
        type(MPAS_TimeInterval_type) :: rewind_interval
        real(kind=RKIND) :: elapsed_seconds
        real(real64) :: rewind_seconds, tolerance
        integer :: ierr

        epoch_seconds = 0.0_real64
        error_code = 0
        error_message = ''

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_real64 .or. &
            .not. ieee_is_finite(real(model_dt, real64)) .or. model_dt <= 0.0_RKIND) then
            error_code = 1
            error_message = '[MIEM] chemistry and model timesteps must be finite and positive'
            return
        end if

        rewind_seconds = time_step - real(model_dt, real64)
        tolerance = 32.0_real64 * epsilon(time_step) * &
                    max(1.0_real64, abs(time_step), abs(real(model_dt, real64)))
        if (rewind_seconds < -tolerance) then
            error_code = 1
            error_message = '[MIEM] accumulated chemistry timestep is shorter than the MPAS timestep'
            return
        end if
        rewind_seconds = max(0.0_real64, rewind_seconds)

        ierr = 0
        call mpas_set_timeInterval(rewind_interval, dt=real(rewind_seconds, RKIND), ierr=ierr)
        if (ierr /= 0) then
            error_code = 1
            error_message = '[MIEM] could not construct the chemistry-interval rewind'
            return
        end if
        interval_start = curr_time - rewind_interval

        ierr = 0
        call mpas_set_time(unix_epoch, dateTimeString='1970-01-01_00:00:00', ierr=ierr)
        if (ierr /= 0) then
            error_code = 1
            error_message = '[MIEM] could not construct the Unix epoch with MPAS timekeeping'
            return
        end if

        ierr = 0
        call mpas_get_timeInterval(interval_start - unix_epoch, dt=elapsed_seconds, ierr=ierr)
        if (ierr /= 0) then
            error_code = 1
            error_message = '[MIEM] could not convert the chemistry interval start to epoch seconds'
            return
        end if
        epoch_seconds = real(elapsed_seconds, real64)
        if (.not. ieee_is_finite(epoch_seconds)) then
            error_code = 1
            error_message = '[MIEM] MPAS timekeeping returned nonfinite epoch seconds'
            return
        end if

    end subroutine chemistry_miem_epoch
#endif


    !------------------------------------------------------------------------
    !  routine chem_env_fill
    !
    !> \brief Derive one thermodynamic chemistry environment for a step.
    !------------------------------------------------------------------------
    subroutine chem_env_fill(env, mesh, state, diag, time_lev, nVertLevels, nCellsSolve, &
                             error_code, error_message)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_array, mpas_pool_get_dimension

        type(chem_env_t), intent(inout) :: env
        type(mpas_pool_type), intent(in) :: mesh
        type(mpas_pool_type), intent(in) :: state
        type(mpas_pool_type), intent(in) :: diag
        integer, intent(in) :: time_lev
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: nCellsSolve
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        real(kind=RKIND), dimension(:,:,:), pointer :: scalars
        real(kind=RKIND), dimension(:,:), pointer :: rho_zz, theta_m, zz, exner
        real(kind=RKIND), dimension(:,:), pointer :: pressure_p, pressure_b
        integer, pointer :: index_qv
        integer :: k, iCell
        logical :: have_pressure_fields

        error_code = 0
        error_message = ''

        if (allocated(env%temperature)) deallocate(env%temperature)
        if (allocated(env%pressure)) deallocate(env%pressure)
        if (allocated(env%rho_dry)) deallocate(env%rho_dry)
        nullify(env%zgrid)

        nullify(scalars)
        nullify(rho_zz)
        nullify(theta_m)
        nullify(zz)
        nullify(exner)
        nullify(pressure_p)
        nullify(pressure_b)
        nullify(index_qv)

        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
        call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
        call mpas_pool_get_array(mesh, 'zz', zz)
        call mpas_pool_get_array(mesh, 'zgrid', env%zgrid)
        call mpas_pool_get_array(diag, 'exner', exner)
        call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)

        if (.not. associated(scalars) .or. .not. associated(rho_zz) .or. &
            .not. associated(theta_m) .or. .not. associated(zz) .or. &
            .not. associated(exner) .or. .not. associated(index_qv)) then
            error_code = 1
            error_message = 'Required MPAS arrays not available in chem_env_fill.'
            return
        end if

        allocate(env%temperature(nVertLevels, nCellsSolve))
        allocate(env%pressure(nVertLevels, nCellsSolve))
        allocate(env%rho_dry(nVertLevels, nCellsSolve))

        ! Evaluate the pressure-source test once, not per (k,iCell).
        have_pressure_fields = associated(pressure_p) .and. associated(pressure_b)

        do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
                ! MPAS scalars are dry-air mass mixing ratios; unit conversion
                ! must use dry density.
                env%rho_dry(k, iCell) = zz(k, iCell) * rho_zz(k, iCell)
                env%temperature(k, iCell) = (theta_m(k, iCell) / &
                    (1.0_RKIND + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell)
            end do
        end do

        if (have_pressure_fields) then
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    env%pressure(k, iCell) = pressure_p(k, iCell) + pressure_b(k, iCell)
                end do
            end do
        else
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    ! Ideal-gas fallback. rho_dry is dry-air density, so the
                    ! moist-air factor is (1 + rvord*qv), matching how
                    ! temperature is derived from theta_m above.
                    env%pressure(k, iCell) = env%rho_dry(k, iCell) * rgas * &
                        env%temperature(k, iCell) * &
                        (1.0_RKIND + rvord * scalars(index_qv, k, iCell))
                end do
            end do
        end if

    end subroutine chem_env_fill


    !------------------------------------------------------------------------
    !  routine chemistry_prepare_miem_diagnostic_selection
    !
    !> \brief Expose the shared MIEM diagnostic parser to core setup
    !------------------------------------------------------------------------
    subroutine chemistry_prepare_miem_diagnostic_selection(sector_text, category_text, &
                                                           number_of_vertical_levels, &
                                                           number_of_species, &
                                                           layered_diagnostics, &
                                                           max_diagnostic_fields, sectors, &
                                                           sector_labels, categories, &
                                                           error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_miem, only: miem_prepare_diagnostic_selection
#endif

        character(len=*), intent(in) :: sector_text, category_text
        integer, intent(in) :: number_of_vertical_levels, number_of_species
        logical, intent(in) :: layered_diagnostics
        integer, intent(in) :: max_diagnostic_fields
        character(len=64), allocatable, intent(out) :: sectors(:)
        character(len=64), allocatable, intent(out) :: sector_labels(:)
        integer, allocatable, intent(out) :: categories(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        call miem_prepare_diagnostic_selection(sector_text, category_text, &
                                               number_of_vertical_levels, number_of_species, &
                                               layered_diagnostics, max_diagnostic_fields, &
                                               sectors, sector_labels, categories, &
                                               error_code, error_message)
#else
        allocate(sectors(0), sector_labels(0), categories(0))
        error_code = 1
        error_message = '[MIEM] diagnostic selection requires a MUSICA-enabled build'
#endif

    end subroutine chemistry_prepare_miem_diagnostic_selection


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
    subroutine chemistry_init(configs, dimensions, structs, dminfo, blocklist)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_init, resolve_mpas_indices, &
                               musica_cache_photo_indices, musica_cache_emission_indices
        use mpas_miem, only: miem_init, miem_is_enabled, miem_get_species, &
                             miem_validate_owned_ids, miem_prepare_diagnostic_selection, &
                             miem_configure_net_flux_species
        use mpas_tuvx, only: tuvx_init, tuvx_is_enabled, &
                             tuvx_n_rates, tuvx_rate_names
        use mpas_dmpar, only: mpas_dmpar_max_int
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, dm_info, block_type
        use mpas_pool_routines, only: mpas_pool_get_array, mpas_pool_get_config, &
                                      mpas_pool_get_dimension, mpas_pool_get_subpool

        type (mpas_pool_type), intent(in) :: configs
        type (mpas_pool_type), intent(in) :: dimensions
        type (mpas_pool_type), intent(inout) :: structs
        type (dm_info), pointer :: dminfo
        type (block_type), pointer :: blocklist

#ifdef MPAS_USE_MUSICA
        character(len=StrKIND), pointer :: filepath_ptr
        integer                         :: error_code
        character(len=:), allocatable   :: error_message
        integer                         :: nCells
        integer                         :: nCellsSolve
        integer                         :: nVertLevels
        integer, pointer                :: nCells_ptr
        integer, pointer                :: nCellsSolve_ptr
        integer, pointer                :: nVertLevels_ptr
        type (mpas_pool_type), pointer  :: statePool
        type (mpas_pool_type), pointer  :: diagPool
        type (mpas_pool_type), pointer  :: meshPool
        real(kind=RKIND), pointer       :: j_no2_ptr, nox_tau_ptr
        real(kind=RKIND)                :: j_no2_val, nox_tau_val
        real(kind=RKIND), pointer       :: lat_ptr, lon_ptr
        logical, pointer                :: use_grid_coords_ptr
        character(len=StrKIND), pointer :: tuvx_filepath_ptr
        logical, pointer                :: tuvx_top_extension_ptr
        character(len=StrKIND), pointer :: tuvx_extension_file_ptr
        character(len=StrKIND), pointer :: tuvx_upper_column_mode_ptr
        character(len=StrKIND), pointer :: tuvx_prescribed_field_file_ptr
        logical                         :: tuvx_top_extension
        logical                         :: tuvx_configured
        character(len=StrKIND)          :: tuvx_extension_file
        character(len=StrKIND)          :: tuvx_upper_column_mode
        character(len=StrKIND)          :: tuvx_prescribed_field_file
        integer, pointer                :: idx_qO3_ptr
        integer, pointer                :: idx_qc_ptr, idx_qr_ptr
        character(len=StrKIND), pointer :: miem_filepath_ptr
        character(len=StrKIND), pointer :: miem_net_flux_species_ptr
        character(len=StrKIND), pointer :: miem_diag_sectors_ptr
        character(len=StrKIND), pointer :: miem_diag_categories_ptr
        logical, pointer                :: miem_layered_diagnostics_ptr
        integer, pointer                :: miem_max_diagnostic_fields_ptr
        character(len=StrKIND), pointer :: calendar_type_ptr
        character(len=StrKIND)          :: miem_filepath
        character(len=StrKIND)          :: miem_net_flux_species_text
        character(len=StrKIND)          :: miem_diag_sector_text
        character(len=StrKIND)          :: miem_diag_category_text
        character(len=StrKIND)          :: calendar_type
        character(len=64), allocatable  :: miem_diag_sectors(:)
        character(len=64), allocatable  :: miem_diag_sector_labels(:)
        integer, allocatable            :: miem_diag_categories(:)
        logical                         :: miem_layered_diagnostics
        integer                         :: miem_max_diagnostic_fields
        integer, dimension(:), pointer  :: indexToCellID
        real(kind=RKIND), dimension(:), pointer :: miem_areaCell
        real(kind=RKIND), dimension(:), pointer :: miem_latCell, miem_lonCell
        real(kind=RKIND), dimension(:), pointer :: miem_xCell, miem_yCell, miem_zCell
        logical, pointer                :: miem_on_a_sphere_ptr
        logical, pointer                :: miem_is_periodic_ptr
        real(kind=RKIND), pointer       :: miem_sphere_radius_ptr
        integer                         :: local_max_cell_id
        integer                         :: global_nCells
        integer                         :: n_miem_species
        character(len=64), allocatable  :: miem_species(:)
#endif

        call mpas_log_write('Initializing chemistry packages...')

        if (associated(blocklist)) then
            if (associated(blocklist % next)) then
                call chem_fatal('Chemistry requires one block per MPI task ' // &
                                '(config_number_of_blocks); multi-block chemistry is not supported.')
            end if
        end if

#ifdef MPAS_USE_MUSICA
        call mpas_pool_get_dimension(dimensions, 'nVertLevels', nVertLevels_ptr)
        nVertLevels = nVertLevels_ptr

        call mpas_pool_get_dimension(dimensions, 'nCells', nCells_ptr)
        if (associated(nCells_ptr)) then
            nCells = nCells_ptr
        else
            nCells = 0
        end if

        call mpas_pool_get_dimension(dimensions, 'nCellsSolve', nCellsSolve_ptr)
        if (associated(nCellsSolve_ptr)) then
            nCellsSolve = nCellsSolve_ptr
        else
            nCellsSolve = nCells
        end if

        call mpas_pool_get_config(configs, 'config_micm_file', filepath_ptr)

        ! Read prescribed photolysis and sink parameters
        j_no2_val = 0.0_RKIND
        nox_tau_val = 0.0_RKIND
        nullify(j_no2_ptr)
        nullify(nox_tau_ptr)
        call mpas_pool_get_config(configs, 'config_j_no2_max', j_no2_ptr)
        if (associated(j_no2_ptr)) j_no2_val = j_no2_ptr
        call mpas_pool_get_config(configs, 'config_lnox_nox_tau', nox_tau_ptr)
        if (associated(nox_tau_ptr)) nox_tau_val = nox_tau_ptr

        ! Read reference solve flag
        block
            logical, pointer :: ref_solve_ptr
            nullify(ref_solve_ptr)
            call mpas_pool_get_config(configs, 'config_chemistry_ref_solve', ref_solve_ptr)
            if (associated(ref_solve_ptr)) use_ref_solve = ref_solve_ptr
        end block

        ! Read chemistry sub-step count (split MICM dt for stiff mechanisms)
        block
            integer, pointer :: chem_substeps_ptr
            nullify(chem_substeps_ptr)
            call mpas_pool_get_config(configs, 'config_chem_substeps', chem_substeps_ptr)
            if (associated(chem_substeps_ptr)) chem_substeps_val = chem_substeps_ptr
            if (chem_substeps_val < 1) chem_substeps_val = 1
        end block

        ! Read chemistry solve interval. The default 0.0 preserves every-step chemistry.
        block
            real(kind=RKIND), pointer :: chemistry_interval_ptr
            nullify(chemistry_interval_ptr)
            call mpas_pool_get_config(configs, 'config_chemistry_interval', chemistry_interval_ptr)
            if (associated(chemistry_interval_ptr)) chemistry_interval = chemistry_interval_ptr
            if (chemistry_interval < 0.0_RKIND) chemistry_interval = 0.0_RKIND
            chemistry_interval_clock = 0.0_RKIND
            chemistry_accumulated_dt = 0.0_RKIND
        end block

        ! Read maximum consecutive step failures before chemistry aborts
        block
            integer, pointer :: chem_max_failures_ptr
            nullify(chem_max_failures_ptr)
            call mpas_pool_get_config(configs, 'config_chem_max_failures', chem_max_failures_ptr)
            if (associated(chem_max_failures_ptr)) chem_max_consecutive_failures = chem_max_failures_ptr
            if (chem_max_consecutive_failures < 1) chem_max_consecutive_failures = 1
        end block

        ! Read MICM relative tolerance (tighten from default 1e-6 for stiff mechanisms)
        block
            real(kind=RKIND), pointer :: reltol_ptr
            nullify(reltol_ptr)
            call mpas_pool_get_config(configs, 'config_micm_relative_tolerance', reltol_ptr)
            if (associated(reltol_ptr)) micm_relative_tolerance_val = reltol_ptr
        end block

        ! Cache solar geometry config values for per-step SZA computation
        chem_j_no2_max = j_no2_val
        nullify(lat_ptr)
        nullify(lon_ptr)
        call mpas_pool_get_config(configs, 'config_chemistry_latitude', lat_ptr)
        if (associated(lat_ptr)) chem_lat = lat_ptr
        call mpas_pool_get_config(configs, 'config_chemistry_longitude', lon_ptr)
        if (associated(lon_ptr)) chem_lon = lon_ptr
        nullify(use_grid_coords_ptr)
        call mpas_pool_get_config(configs, 'config_chemistry_use_grid_coords', use_grid_coords_ptr)
        if (associated(use_grid_coords_ptr)) chem_use_grid_coords = use_grid_coords_ptr

        block
            real(kind=RKIND), pointer :: tuvx_update_interval_ptr
            nullify(tuvx_update_interval_ptr)
            call mpas_pool_get_config(configs, 'config_tuvx_update_interval', tuvx_update_interval_ptr)
            if (associated(tuvx_update_interval_ptr)) tuvx_update_interval = tuvx_update_interval_ptr
        end block

        ! Initialize MICM (species already injected into scalars by atm_setup_block)
        call musica_init(filepath_ptr, nCellsSolve * nVertLevels, &
                         nox_tau_val, &
                         real(micm_relative_tolerance_val, kind(1.0d0)), &
                         use_ref_solve, &
                         error_code, error_message)

        if (error_code /= 0) then
            call chem_fatal(error_message)
        end if

        ! Initialize MIEM only after MICM is available. An empty path is a
        ! deliberate no-op: no mesh collective, mechanism object, or inventory
        ! access occurs on the disabled path.
        miem_filepath = ''
        miem_net_flux_species_text = ''
        miem_diag_sector_text = ''
        miem_diag_category_text = ''
        miem_layered_diagnostics = .false.
        miem_max_diagnostic_fields = 256
        calendar_type = 'gregorian'
        nullify(miem_filepath_ptr)
        nullify(miem_net_flux_species_ptr)
        nullify(miem_diag_sectors_ptr, miem_diag_categories_ptr)
        nullify(miem_layered_diagnostics_ptr, miem_max_diagnostic_fields_ptr)
        nullify(calendar_type_ptr)
        call mpas_pool_get_config(configs, 'config_miem_file', miem_filepath_ptr)
        call mpas_pool_get_config(configs, 'config_miem_net_flux_species', &
                                  miem_net_flux_species_ptr)
        call mpas_pool_get_config(configs, 'config_miem_diagnostic_sectors', miem_diag_sectors_ptr)
        call mpas_pool_get_config(configs, 'config_miem_diagnostic_categories', miem_diag_categories_ptr)
        call mpas_pool_get_config(configs, 'config_miem_layered_diagnostics', miem_layered_diagnostics_ptr)
        call mpas_pool_get_config(configs, 'config_miem_max_diagnostic_fields', &
                                  miem_max_diagnostic_fields_ptr)
        call mpas_pool_get_config(configs, 'config_calendar_type', calendar_type_ptr)
        if (associated(miem_filepath_ptr)) miem_filepath = miem_filepath_ptr
        if (associated(miem_net_flux_species_ptr)) then
            miem_net_flux_species_text = miem_net_flux_species_ptr
        end if
        if (associated(miem_diag_sectors_ptr)) miem_diag_sector_text = miem_diag_sectors_ptr
        if (associated(miem_diag_categories_ptr)) miem_diag_category_text = miem_diag_categories_ptr
        if (associated(miem_layered_diagnostics_ptr)) then
            miem_layered_diagnostics = miem_layered_diagnostics_ptr
        end if
        if (associated(miem_max_diagnostic_fields_ptr)) then
            miem_max_diagnostic_fields = miem_max_diagnostic_fields_ptr
        end if
        if (associated(calendar_type_ptr)) calendar_type = calendar_type_ptr

        global_nCells = 0
        if (len_trim(miem_filepath) > 0) then
            if (trim(calendar_type) /= 'gregorian') then
                call chem_fatal('[MIEM] config_calendar_type must be gregorian when emissions are enabled; got ' // &
                                trim(calendar_type))
            end if

            nullify(meshPool)
            nullify(indexToCellID)
            nullify(miem_areaCell, miem_latCell, miem_lonCell)
            nullify(miem_xCell, miem_yCell, miem_zCell)
            nullify(miem_on_a_sphere_ptr, miem_is_periodic_ptr, miem_sphere_radius_ptr)
            call mpas_pool_get_subpool(structs, 'mesh', meshPool)
            if (.not. associated(meshPool)) then
                call chem_fatal('[MIEM] mesh pool is unavailable during initialization.')
            end if
            call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
            if (.not. associated(indexToCellID)) then
                call chem_fatal('[MIEM] indexToCellID is unavailable during initialization.')
            end if
            call mpas_pool_get_array(meshPool, 'areaCell', miem_areaCell)
            call mpas_pool_get_array(meshPool, 'latCell', miem_latCell)
            call mpas_pool_get_array(meshPool, 'lonCell', miem_lonCell)
            call mpas_pool_get_array(meshPool, 'xCell', miem_xCell)
            call mpas_pool_get_array(meshPool, 'yCell', miem_yCell)
            call mpas_pool_get_array(meshPool, 'zCell', miem_zCell)
            if (.not. associated(miem_areaCell) .or. .not. associated(miem_latCell) .or. &
                .not. associated(miem_lonCell) .or. .not. associated(miem_xCell) .or. &
                .not. associated(miem_yCell) .or. .not. associated(miem_zCell)) then
                call chem_fatal('[MIEM] complete cell-center mesh geometry is unavailable during initialization.')
            end if
            call mpas_pool_get_config(meshPool, 'on_a_sphere', miem_on_a_sphere_ptr)
            call mpas_pool_get_config(meshPool, 'is_periodic', miem_is_periodic_ptr)
            call mpas_pool_get_config(meshPool, 'sphere_radius', miem_sphere_radius_ptr)
            if (.not. associated(miem_on_a_sphere_ptr) .or. &
                .not. associated(miem_is_periodic_ptr) .or. &
                .not. associated(miem_sphere_radius_ptr)) then
                call chem_fatal('[MIEM] mesh geometry attributes are unavailable during initialization.')
            end if

            call miem_validate_owned_ids(indexToCellID, nCellsSolve, local_max_cell_id, &
                                         error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
            if (.not. associated(dminfo)) then
                call chem_fatal('[MIEM] MPI decomposition metadata is unavailable.')
            end if
            call mpas_dmpar_max_int(dminfo, local_max_cell_id, global_nCells)
        end if

        if (len_trim(miem_filepath) > 0) then
            call miem_prepare_diagnostic_selection(miem_diag_sector_text, &
                                                   miem_diag_category_text, nVertLevels, 0, &
                                                   miem_layered_diagnostics, &
                                                   miem_max_diagnostic_fields, miem_diag_sectors, &
                                                   miem_diag_sector_labels, miem_diag_categories, &
                                                   error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
                return
            end if
            call miem_init(miem_filepath, global_nCells, nVertLevels, &
                           indexToCellID(1:nCellsSolve), miem_on_a_sphere_ptr, &
                           miem_is_periodic_ptr, miem_sphere_radius_ptr, &
                           miem_areaCell(1:nCellsSolve), miem_latCell(1:nCellsSolve), &
                           miem_lonCell(1:nCellsSolve), miem_xCell(1:nCellsSolve), &
                           miem_yCell(1:nCellsSolve), miem_zCell(1:nCellsSolve), &
                           error_code, error_message, &
                           diagnostic_sectors=miem_diag_sectors, &
                           diagnostic_categories=miem_diag_categories, &
                           layered_diagnostics=miem_layered_diagnostics, &
                           max_diagnostic_fields=miem_max_diagnostic_fields)
            if (allocated(miem_diag_sectors)) deallocate(miem_diag_sectors)
            if (allocated(miem_diag_sector_labels)) deallocate(miem_diag_sector_labels)
            if (allocated(miem_diag_categories)) deallocate(miem_diag_categories)
        else
            call miem_init(miem_filepath, error_code=error_code, error_message=error_message)
        end if
        if (error_code /= 0) then
            call chem_fatal(error_message)
        end if
        if (len_trim(miem_filepath) > 0 .and. .not. miem_is_enabled()) then
            call chem_fatal('[MIEM] configured emissions did not initialize.')
        end if

        if (miem_is_enabled()) then
            call miem_configure_net_flux_species(miem_net_flux_species_text, &
                                                 error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
            call miem_get_species(n_miem_species, miem_species)
            call musica_cache_emission_indices(miem_species, error_code, error_message, &
                                               net_flux_species=miem_net_flux_species_text)
            if (allocated(miem_species)) deallocate(miem_species)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
        else if (len_trim(miem_net_flux_species_text) > 0) then
            call chem_fatal('[MIEM] config_miem_net_flux_species requires config_miem_file.')
        end if

        ! Resolve MPAS tracer indices for all chemistry species
        nullify(statePool)
        call mpas_pool_get_subpool(structs, 'state', statePool)
        if (associated(statePool)) then
            call resolve_mpas_indices(statePool, error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
            call chemistry_resolve_tracer_output_indices(statePool, error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
        end if

        ! Initialize lightning NOx source (no-op if qNO not in mechanism)
        call chemistry_lightning_nox_init(configs, statePool)

        ! Initialize TUV-x if config file is specified
        nullify(tuvx_filepath_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_config_file', tuvx_filepath_ptr)
        tuvx_configured = .false.
        if (associated(tuvx_filepath_ptr)) then
            tuvx_configured = len_trim(tuvx_filepath_ptr) > 0
        end if
        nullify(tuvx_top_extension_ptr)
        nullify(tuvx_extension_file_ptr)
        nullify(tuvx_upper_column_mode_ptr)
        nullify(tuvx_prescribed_field_file_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_top_extension', tuvx_top_extension_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_extension_file', tuvx_extension_file_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_upper_column_mode', tuvx_upper_column_mode_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_prescribed_field_file', &
                                  tuvx_prescribed_field_file_ptr)
        tuvx_top_extension = .false.
        tuvx_extension_file = ''
        tuvx_upper_column_mode = 'none'
        tuvx_prescribed_field_file = ''
        if (associated(tuvx_top_extension_ptr)) tuvx_top_extension = tuvx_top_extension_ptr
        if (associated(tuvx_extension_file_ptr)) tuvx_extension_file = tuvx_extension_file_ptr
        if (associated(tuvx_upper_column_mode_ptr)) then
            tuvx_upper_column_mode = trim(adjustl(tuvx_upper_column_mode_ptr))
        end if
        if (associated(tuvx_prescribed_field_file_ptr)) then
            tuvx_prescribed_field_file = tuvx_prescribed_field_file_ptr
        end if

        select case (trim(tuvx_upper_column_mode))
        case ('none')
            if (tuvx_top_extension .or. len_trim(tuvx_extension_file) > 0 .or. &
                len_trim(tuvx_prescribed_field_file) > 0) then
                call chem_fatal('[TUV-x] upper-column mode none is inconsistent with configured extension inputs.')
            end if
        case ('legacy_static')
            if (.not. tuvx_top_extension .or. len_trim(tuvx_extension_file) == 0 .or. &
                len_trim(tuvx_prescribed_field_file) > 0) then
                call chem_fatal('[TUV-x] legacy_static requires the compatibility flag and only a CSV file.')
            end if
        case ('spatial_climatology')
            if (.not. tuvx_top_extension .or. len_trim(tuvx_extension_file) > 0 .or. &
                len_trim(tuvx_prescribed_field_file) == 0) then
                call chem_fatal('[TUV-x] spatial_climatology requires the compatibility flag and only a NetCDF package.')
            end if
            if (trim(calendar_type) /= 'gregorian') then
                call chem_fatal('[Prescribed fields] spatial climatology requires config_calendar_type=gregorian.')
            end if
        case default
            call chem_fatal('[TUV-x] config_tuvx_upper_column_mode must be none, legacy_static, or spatial_climatology.')
        end select

        if (.not. tuvx_configured .and. trim(tuvx_upper_column_mode) /= 'none') then
            call chem_fatal('[TUV-x] an upper-column mode requires config_tuvx_config_file.')
        end if

        if (tuvx_configured) then
            if (trim(tuvx_upper_column_mode) == 'spatial_climatology') then
                block
                    use mpas_prescribed_fields, only: prescribed_fields_init

                    nullify(meshPool)
                    nullify(indexToCellID)
                    nullify(miem_areaCell, miem_latCell, miem_lonCell)
                    call mpas_pool_get_subpool(structs, 'mesh', meshPool)
                    if (.not. associated(meshPool)) then
                        call chem_fatal('[Prescribed fields] mesh pool is unavailable during initialization.')
                    end if
                    call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
                    call mpas_pool_get_array(meshPool, 'areaCell', miem_areaCell)
                    call mpas_pool_get_array(meshPool, 'latCell', miem_latCell)
                    call mpas_pool_get_array(meshPool, 'lonCell', miem_lonCell)
                    if (.not. associated(indexToCellID) .or. .not. associated(miem_areaCell) .or. &
                        .not. associated(miem_latCell) .or. .not. associated(miem_lonCell)) then
                        call chem_fatal('[Prescribed fields] exact-grid mesh metadata is unavailable.')
                    end if
                    if (.not. associated(dminfo)) then
                        call chem_fatal('[Prescribed fields] MPI decomposition metadata is unavailable.')
                    end if
                    if (nCellsSolve < 1) then
                        call chem_fatal('[Prescribed fields] each chemistry rank must own at least one cell.')
                    end if
                    local_max_cell_id = maxval(indexToCellID(1:nCellsSolve))
                    call mpas_dmpar_max_int(dminfo, local_max_cell_id, global_nCells)
                    call prescribed_fields_init(tuvx_prescribed_field_file, global_nCells, &
                                                indexToCellID(1:nCellsSolve), &
                                                miem_areaCell(1:nCellsSolve), &
                                                miem_latCell(1:nCellsSolve), &
                                                miem_lonCell(1:nCellsSolve), &
                                                error_code, error_message)
                    if (error_code /= 0) call chem_fatal(error_message)
                end block
            end if

            call tuvx_init(tuvx_filepath_ptr, nVertLevels, tuvx_upper_column_mode, &
                           tuvx_extension_file, error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
            use_tuvx = tuvx_is_enabled()

            ! Resolve qO3 index (required for TUV-x)
            if (use_tuvx .and. associated(statePool)) then
                nullify(idx_qO3_ptr)
                call mpas_pool_get_dimension(statePool, 'index_qO3', idx_qO3_ptr)
                if (associated(idx_qO3_ptr) .and. idx_qO3_ptr > 0) then
                    idx_qO3 = idx_qO3_ptr
                    call mpas_log_write('[Chemistry] Resolved index_qO3 = $i for TUV-x', &
                        intArgs=[idx_qO3])
                else
                    use_tuvx = .false.
                    idx_qO3 = -1
                    call chem_fatal('qO3 tracer not found; TUV-x requires O3 and cannot continue.')
                end if
                ! Resolve cloud water index (required for Phase 3 cloud opacity)
                nullify(idx_qc_ptr)
                call mpas_pool_get_dimension(statePool, 'index_qc', idx_qc_ptr)
                if (associated(idx_qc_ptr) .and. idx_qc_ptr > 0) then
                    idx_qc = idx_qc_ptr
                    call mpas_log_write('[Chemistry] Resolved index_qc = $i for TUV-x cloud opacity', &
                        intArgs=[idx_qc])
                else
                    idx_qc = -1
                    call mpas_log_write('[Chemistry] qc tracer not found; TUV-x will use clear-sky only.')
                end if

                ! Resolve rain water index (optional — zero fallback if absent)
                nullify(idx_qr_ptr)
                call mpas_pool_get_dimension(statePool, 'index_qr', idx_qr_ptr)
                if (associated(idx_qr_ptr) .and. idx_qr_ptr > 0) then
                    idx_qr = idx_qr_ptr
                    call mpas_log_write('[Chemistry] Resolved index_qr = $i for TUV-x cloud opacity', &
                        intArgs=[idx_qr])
                else
                    idx_qr = -1
                    call mpas_log_write('[Chemistry] qr tracer not found; rain opacity will be zero.')
                end if
            else if (use_tuvx) then
                use_tuvx = .false.
                idx_qO3 = -1
                call chem_fatal('State pool unavailable; cannot resolve qO3 for TUV-x.')
            end if
        else
            call mpas_log_write('[Chemistry] TUV-x disabled (no config_tuvx_config_file); using Phase 1 fallback.')
        end if

        ! Cache MICM rate-parameter indices for every photolysis rate we plan
        ! to write. When TUV-x is active, names come from its photolysis
        ! ordering. When TUV-x is disabled, a positive config_j_no2_max
        ! enables the single-rate Phase-1 fallback ('jNO2').
        block
            character(len=64), allocatable :: photo_names(:)
            integer                        :: n_rates
            if (use_tuvx) then
                n_rates = tuvx_n_rates()
                if (n_rates < 1) then
                    call chem_fatal('TUV-x reported zero photolysis rates.')
                end if
                allocate(photo_names(n_rates))
                call tuvx_rate_names(photo_names)
            else if (chem_j_no2_max > 0.0_RKIND) then
                n_rates = 1
                allocate(photo_names(1))
                photo_names(1) = 'jNO2'
            else
                n_rates = 0
                call mpas_log_write('[Chemistry] TUV-x disabled and config_j_no2_max = 0; ' // &
                    'no photolysis rate parameters will be driven.')
            end if
            if (n_rates > 0) then
                call musica_cache_photo_indices(photo_names, error_code, error_message)
                if (error_code /= 0) then
                    call chem_fatal(error_message)
                end if
                deallocate(photo_names)
            end if
        end block

        ! Seed MPAS tracers only after every chemistry, emissions, photolysis,
        ! and prescribed-field configuration has passed. A fail-closed input
        ! error must not partially mutate prognostic chemistry state.
        call chemistry_seed_chem(dimensions, structs, dminfo)

        ! Initialize photolysis diagnostics before any initial history write.
        ! Runtime-injected j_<reaction> fields are required for each cached
        ! photolysis rate parameter.
        nullify(diagPool)
        call mpas_pool_get_subpool(structs, 'diag', diagPool)
        if (associated(diagPool)) then
            call zero_photolysis_diag(diagPool)
            call zero_emissions_diag(diagPool, error_code, error_message)
            if (error_code /= 0) then
                call chem_fatal(error_message)
            end if
        end if
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
    subroutine chemistry_step(dt, currTime, mesh, state, diag, dimensions, time_lev)

#ifdef MPAS_USE_MUSICA
        use iso_fortran_env, only: real64
        use mpas_musica, only: musica_step, musica_step_ref, log_column_comparison, &
                               copy_state_to_ref, musica_set_photolysis_rates, &
                               musica_set_emission_fluxes, &
                               n_photo_rp, n_chem_species, chem_species
        use mpas_miem, only: miem_is_enabled, miem_run, miem_commit_mass, &
                             miem_get_diagnostic_fluxes
        use mpas_lightning_nox, only: lightning_nox_inject
        use mpas_solar_geometry, only: solar_cos_sza
        use mpas_tuvx, only: tuvx_compute_photolysis, tuvx_rate_names, &
                             tuvx_uses_spatial_climatology, tuvx_n_extension_layers
        use mpas_prescribed_fields, only: prescribed_fields_update_o3
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_dimension, mpas_pool_get_array
        use mpas_timekeeping, only: mpas_get_time, MPAS_Time_type

        real (kind=RKIND), intent(in) :: dt
        type (MPAS_Time_type), intent(in) :: currTime
        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(inout) :: state
        type (mpas_pool_type), intent(inout) :: diag
        type (mpas_pool_type), intent(in) :: dimensions
        integer, intent(in) :: time_lev

#ifdef MPAS_USE_MUSICA
        real(real64)                  :: time_step
        real(kind=RKIND)              :: tuvx_time_before_step
        integer                       :: error_code
        character(len=:), allocatable :: error_message
        integer                       :: nCells, nCellsSolve, nVertLevels
        integer, pointer              :: nCells_ptr, nCellsSolve_ptr, nVertLevels_ptr
        integer                       :: probe_cell
        logical, save                 :: first_call = .true.
        logical                       :: photolysis_updated
        logical                       :: photolysis_cadence_active
        type(chem_env_t)              :: env
        real(real64)                  :: miem_epoch_seconds
        real(real64), allocatable     :: local_emission_flux(:,:)
        real(real64), allocatable     :: local_emission_layer_flux(:,:,:)
        real(real64), allocatable     :: local_emission_sector_flux(:,:,:)
        real(real64), allocatable     :: local_emission_category_flux(:,:,:)
        real(real64), allocatable     :: local_emission_sector_layer_flux(:,:,:,:)
        real(real64), allocatable     :: local_emission_category_layer_flux(:,:,:,:)
        real(real64), allocatable     :: owned_cell_area(:)
        real(kind=RKIND), dimension(:), pointer :: areaCell

        ! Solar geometry variables
        integer :: year, month, day, DoY, hour, minute, second
        real(kind=RKIND) :: hour_utc
        real(kind=RKIND), allocatable :: cos_sza_cell(:)
        real(kind=RKIND), dimension(:), pointer :: latCell, lonCell

        ! Lightning NOx injection arrays
        real(kind=RKIND), dimension(:,:,:), pointer :: scalars
        real(kind=RKIND), allocatable :: scalars_before_chem(:,:,:)
        real(kind=RKIND), dimension(:,:), pointer :: w

        ! Photolysis variables
        real(kind=RKIND), allocatable :: photo_rates(:,:,:)  ! (n_rates, nVertLevels, nCells)
        real(kind=RKIND), allocatable :: qc_col(:), qr_col(:)
        real(kind=RKIND), allocatable :: prescribed_o3_density(:,:)
        integer :: n_upper_levels
        integer :: iCell

        time_step = dt
        if (chemistry_interval > 0.0_RKIND) then
            chemistry_interval_clock = chemistry_interval_clock + dt
            chemistry_accumulated_dt = chemistry_accumulated_dt + dt
            if (chemistry_interval_clock < chemistry_interval) return
            time_step = real(chemistry_accumulated_dt, real64)
        end if

        ! Get dimensions
        call mpas_pool_get_dimension(dimensions, 'nCells', nCells_ptr)
        call mpas_pool_get_dimension(dimensions, 'nCellsSolve', nCellsSolve_ptr)
        call mpas_pool_get_dimension(dimensions, 'nVertLevels', nVertLevels_ptr)
        nCells = nCells_ptr
        if (associated(nCellsSolve_ptr)) then
            nCellsSolve = nCellsSolve_ptr
        else
            nCellsSolve = nCells
        end if
        nVertLevels = nVertLevels_ptr

        ! Probe cell near center of domain for advection diagnostics
        if (nCellsSolve > 0) then
            probe_cell = max(1, nCellsSolve / 2)
        else
            probe_cell = 0
        end if
        photolysis_updated = .false.
        photolysis_cadence_active = .false.

        call chem_env_fill(env, mesh, state, diag, time_lev, nVertLevels, nCellsSolve, &
                           error_code, error_message)
        if (error_code /= 0) then
            call chem_step_error(error_message, chem_consecutive_failures)
            if (use_ref_solve) first_call = .true.
            return
        end if

        ! Update photolysis rates. Unified multi-rate path: TUV-x fills
        ! photo_rates(ir,k,iCell) for every reaction it reports; when TUV-x
        ! is disabled, the fallback fills a single 'jNO2' slot with
        ! j_max * max(0, cos_sza). The whole block is gated by
        ! tuvx_update_interval; on skip steps the accumulator just ticks
        ! forward and MICM keeps using the last-set rate parameters. If no
        ! photolysis rate parameters were cached, MICM keeps its defaults.
        if (n_photo_rp > 0) then
            photolysis_cadence_active = .true.
            tuvx_time_before_step = tuvx_time_since_last
            tuvx_time_since_last = tuvx_time_since_last + real(time_step, RKIND)
        end if
        if (photolysis_cadence_active .and. tuvx_time_since_last >= tuvx_update_interval) then
            ! Compute the time-of-day inputs (only needed when an update fires).
            call mpas_get_time(currTime, YYYY=year, MM=month, DD=day, DoY=DoY, &
                               H=hour, M=minute, S=second)
            hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
                     + real(second, RKIND) / 3600.0_RKIND

            ! Per-cell solar zenith angle. Filled per cell when the grid-coords
            ! switch is enabled, broadcast from the namelist scalars otherwise.
            allocate(cos_sza_cell(nCellsSolve))
            if (chem_use_grid_coords) then
                call mpas_pool_get_array(mesh, 'latCell', latCell)
                call mpas_pool_get_array(mesh, 'lonCell', lonCell)
                do iCell = 1, nCellsSolve
                    cos_sza_cell(iCell) = solar_cos_sza(DoY, hour_utc, &
                                                        latCell(iCell) * 180.0_RKIND / pii, &
                                                        lonCell(iCell) * 180.0_RKIND / pii)
                end do
            else
                cos_sza_cell(:) = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)
            end if

            allocate(photo_rates(n_photo_rp, nVertLevels, nCellsSolve))
            photo_rates = 0.0_RKIND

            if (use_tuvx) then
                if (idx_qO3 < 1) then
                    call mpas_log_write('[Chemistry] TUV-x enabled but index_qO3 is unresolved.', &
                        messageType=MPAS_LOG_CRIT)
                    deallocate(photo_rates, cos_sza_cell)
                    return
                end if
                if (.not. associated(env%zgrid)) then
                    call mpas_log_write('[Chemistry] TUV-x enabled but zgrid is unavailable.', &
                        messageType=MPAS_LOG_CRIT)
                    deallocate(photo_rates, cos_sza_cell)
                    return
                end if

                allocate(qc_col(nVertLevels))
                allocate(qr_col(nVertLevels))

                if (tuvx_uses_spatial_climatology()) then
                    n_upper_levels = tuvx_n_extension_layers()
                    if (n_upper_levels < 1) then
                        call mpas_log_write('[Chemistry] Spatial climatology has no upper layers.', &
                                            messageType=MPAS_LOG_CRIT)
                        deallocate(photo_rates, qc_col, qr_col, cos_sza_cell)
                        return
                    end if
                    allocate(prescribed_o3_density(n_upper_levels, nCellsSolve))
                    call prescribed_fields_update_o3(year, month, day, hour, minute, second, &
                                                     prescribed_o3_density, &
                                                     error_code, error_message)
                    if (error_code /= 0) then
                        call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                        deallocate(photo_rates, qc_col, qr_col, cos_sza_cell, &
                                   prescribed_o3_density)
                        return
                    end if
                end if

                call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

                do iCell = 1, nCellsSolve
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

                    if (tuvx_uses_spatial_climatology()) then
                        call tuvx_compute_photolysis(cos_sza_cell(iCell), env%zgrid(:, iCell), &
                                                     env%rho_dry(:, iCell), &
                                                     env%temperature(:, iCell), &
                                                     env%pressure(:, iCell), &
                                                     scalars(idx_qO3, :, iCell), &
                                                     qc_col, qr_col, &
                                                     nVertLevels, &
                                                     photo_rates(:, :, iCell), &
                                                     error_code, error_message, &
                                                     prescribed_o3_density(:, iCell))
                    else
                        call tuvx_compute_photolysis(cos_sza_cell(iCell), env%zgrid(:, iCell), &
                                                     env%rho_dry(:, iCell), &
                                                     env%temperature(:, iCell), &
                                                     env%pressure(:, iCell), &
                                                     scalars(idx_qO3, :, iCell), &
                                                     qc_col, qr_col, &
                                                     nVertLevels, &
                                                     photo_rates(:, :, iCell), &
                                                     error_code, error_message)
                    end if
                    if (error_code /= 0) then
                        call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                        deallocate(photo_rates, qc_col, qr_col, cos_sza_cell)
                        if (allocated(prescribed_o3_density)) deallocate(prescribed_o3_density)
                        return
                    end if
                end do

                deallocate(qc_col, qr_col)
                if (allocated(prescribed_o3_density)) deallocate(prescribed_o3_density)
            else
                ! Phase-1 fallback: single rate 'jNO2' = j_max * max(0, cos_sza).
                ! Assumes n_photo_rp == 1 and the cached rate name is 'jNO2' —
                ! musica_cache_photo_indices was called with that name at init.
                do iCell = 1, nCellsSolve
                    photo_rates(1, :, iCell) = chem_j_no2_max * &
                                                max(0.0_RKIND, cos_sza_cell(iCell))
                end do
            end if

            call musica_set_photolysis_rates(photo_rates, nCellsSolve, &
                                             nVertLevels, n_photo_rp, error_code, error_message)
            if (error_code /= 0) then
                call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                deallocate(photo_rates, cos_sza_cell)
                return
            end if

            photolysis_updated = .true.
            deallocate(cos_sza_cell)
        end if

        ! MIEM is sampled exactly once for this accumulated chemistry interval.
        ! Rank-local column and layer mass fluxes already use owned MPAS cell
        ! order. They are converted to EMIS.* concentration tendencies before either MICM
        ! state is synchronized or solved. The rates remain fixed across every
        ! configured MICM substep. Existing lightning NOx is applied afterward
        ! as a separate operator-split volumetric source; configuring both
        ! sources intentionally combines them and may double-count NOx if the
        ! offline inventory already includes lightning.
        if (miem_is_enabled()) then
            if (.not. associated(env%zgrid)) then
                if (allocated(photo_rates)) deallocate(photo_rates)
                call mpas_log_write('[Chemistry] [MIEM] zgrid is unavailable for surface-flux conversion.', &
                                    messageType=MPAS_LOG_CRIT)
                return
            end if

            call chemistry_miem_epoch(currTime, time_step, dt, miem_epoch_seconds, &
                                      error_code, error_message)
            if (error_code /= 0) then
                if (allocated(photo_rates)) deallocate(photo_rates)
                call mpas_log_write('[Chemistry] ' // trim(error_message), messageType=MPAS_LOG_CRIT)
                return
            end if

            call mpas_timer_start('chem MIEM')
            call miem_run(miem_epoch_seconds, time_step, local_emission_flux, &
                          local_emission_layer_flux, error_code, error_message)
            call mpas_timer_stop('chem MIEM')
            if (error_code /= 0) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                call mpas_log_write('[Chemistry] ' // trim(error_message), messageType=MPAS_LOG_CRIT)
                return
            end if

            call miem_get_diagnostic_fluxes(local_emission_sector_flux, &
                                            local_emission_category_flux, &
                                            local_emission_sector_layer_flux, &
                                            local_emission_category_layer_flux, &
                                            error_code, error_message)
            if (error_code /= 0) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                call mpas_log_write('[Chemistry] ' // trim(error_message), messageType=MPAS_LOG_CRIT)
                return
            end if

            call musica_set_emission_fluxes(local_emission_layer_flux, env%zgrid, &
                                             nCellsSolve, nVertLevels, &
                                             error_code, error_message)
            if (error_code /= 0) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                call mpas_log_write('[Chemistry] ' // trim(error_message), messageType=MPAS_LOG_CRIT)
                return
            end if
        end if

        call mpas_log_write('Stepping chemistry packages...')
        call mpas_timer_start('chemistry')

        ! Step 0: Inject lightning NOx source (pre-MICM, operator split)
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
        call mpas_pool_get_array(state, 'w', w, time_lev)

        ! Snapshot only the chemistry-species rows (not all scalars) so a failed
        ! step can be rolled back. This includes qNO when LNOx is active, since
        ! qNO is itself a chemistry species.
        if (associated(scalars)) then
            block
                integer :: i_snap
                allocate(scalars_before_chem(n_chem_species, nVertLevels, nCellsSolve))
                do i_snap = 1, n_chem_species
                    scalars_before_chem(i_snap, :, :) = &
                        scalars(chem_species(i_snap)%mpas_index, 1:nVertLevels, 1:nCellsSolve)
                end do
            end block
        end if

        if (associated(scalars) .and. associated(w) .and. associated(env%zgrid)) then
            call lightning_nox_inject(real(time_step, RKIND), scalars, w, env%zgrid, &
                                      nVertLevels, nCellsSolve, time_lev, &
                                      temperature=env%temperature)
        end if

        ! Step 1: Copy MPAS state to MICM (MPAS -> MICM)
        call mpas_timer_start('chem MPAS->MICM')
        call chemistry_from_MPAS(env, state, time_lev, nCells, nCellsSolve, nVertLevels, &
                                 error_code, error_message)
        call mpas_timer_stop('chem MPAS->MICM')
        if (error_code /= 0) then
            call chem_step_error(error_message, chem_consecutive_failures)
            call chem_step_rollback()
            return
        end if

        ! On first call, copy coupled state to reference so both start with
        ! identical initial conditions from MPAS. Any subsequent divergence
        ! will indicate advection effects.
        if (first_call .and. use_ref_solve) then
            call copy_state_to_ref()
            call mpas_log_write('[Chemistry] Copied coupled state to reference state (first timestep).')
        end if

        ! Step 2: Solve chemistry (coupled state). Split the chemistry timestep into
        ! config_chem_substeps MICM calls. TUV-x photolysis rates were set
        ! once above and stay frozen across sub-steps. Sub-stepping is needed
        ! for stiff mechanisms where Backward Euler at the full dt lands on
        ! a non-physical implicit fixed point (see
        ! docs/plans/2026-04-18-chapman-nox-chem-box-issue.md).
        block
            real(kind=RKIND) :: sub_dt
            integer          :: i_sub
            sub_dt = time_step / real(chem_substeps_val, RKIND)
            call mpas_timer_start('chem MICM solve')
            do i_sub = 1, chem_substeps_val
                call musica_step(real(sub_dt, kind(time_step)), error_code, error_message)
                if (error_code /= 0) then
                    call mpas_timer_stop('chem MICM solve')
                    call chem_step_error(error_message, chem_consecutive_failures)
                    call chem_step_rollback()
                    return
                end if
            end do
            call mpas_timer_stop('chem MICM solve')

            ! Step 2b: Reference state (uncoupled, for advection diagnostics)
            if (use_ref_solve) then
                call mpas_timer_start('chem MICM ref solve')
                do i_sub = 1, chem_substeps_val
                    call musica_step_ref(real(sub_dt, kind(time_step)), error_code, error_message)
                    if (error_code /= 0) then
                        call mpas_timer_stop('chem MICM ref solve')
                        call chem_step_error(error_message, chem_consecutive_failures)
                        call chem_step_rollback()
                        return
                    end if
                end do
                call mpas_timer_stop('chem MICM ref solve')

                if (nCellsSolve > 0) call log_column_comparison(probe_cell, nVertLevels)
            end if
        end block

        ! Step 3: Copy MICM state back to MPAS (MICM -> MPAS)
        call mpas_timer_start('chem MICM->MPAS')
        call chemistry_to_MPAS(env, state, time_lev, nCells, nCellsSolve, nVertLevels, &
                               error_code, error_message)
        call mpas_timer_stop('chem MICM->MPAS')
        if (error_code /= 0) then
            call chem_step_error(error_message, chem_consecutive_failures)
            call chem_step_rollback()
            return
        end if

        ! A sampled inventory is not considered applied until the complete
        ! chemistry transaction succeeds. Publish the exact owned-cell fluxes
        ! used by this solve and add their area/time integral to rank-local
        ! mass totals only after MICM state has been copied back to MPAS.
        ! Rollback paths above therefore leave both diagnostics and budgets at
        ! their previous successful values.
        if (miem_is_enabled()) then
            nullify(areaCell)
            call mpas_pool_get_array(mesh, 'areaCell', areaCell)
            if (.not. associated(areaCell)) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                if (allocated(scalars_before_chem)) deallocate(scalars_before_chem)
                call mpas_timer_stop('chemistry')
                call mpas_log_write('[Chemistry] [MIEM] areaCell is unavailable for mass accounting.', &
                                    messageType=MPAS_LOG_CRIT)
                return
            end if
            if (size(areaCell) < nCellsSolve) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                if (allocated(scalars_before_chem)) deallocate(scalars_before_chem)
                call mpas_timer_stop('chemistry')
                call mpas_log_write('[Chemistry] [MIEM] areaCell is undersized for mass accounting.', &
                                    messageType=MPAS_LOG_CRIT)
                return
            end if

            allocate(owned_cell_area(nCellsSolve))
            owned_cell_area(:) = real(areaCell(1:nCellsSolve), real64)
            call miem_commit_mass(local_emission_flux, owned_cell_area, time_step, &
                                  error_code, error_message)
            deallocate(owned_cell_area)
            if (error_code /= 0) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                if (allocated(scalars_before_chem)) deallocate(scalars_before_chem)
                call mpas_timer_stop('chemistry')
                call mpas_log_write('[Chemistry] ' // trim(error_message), messageType=MPAS_LOG_CRIT)
                return
            end if

            call chemistry_set_emissions_diag(diag, local_emission_flux, &
                                               local_emission_layer_flux, &
                                               local_emission_sector_flux, &
                                               local_emission_category_flux, &
                                               local_emission_sector_layer_flux, &
                                               local_emission_category_layer_flux, &
                                               nCellsSolve, nVertLevels, &
                                               error_code, error_message)
            if (error_code /= 0) then
                call release_miem_step_arrays()
                if (allocated(photo_rates)) deallocate(photo_rates)
                if (allocated(scalars_before_chem)) deallocate(scalars_before_chem)
                call mpas_timer_stop('chemistry')
                call mpas_log_write('[Chemistry] ' // trim(error_message), messageType=MPAS_LOG_CRIT)
                return
            end if
        end if

        if (photolysis_updated) then
            call chemistry_set_photolysis_diag(diag, photo_rates, nCells, nCellsSolve, nVertLevels)
            deallocate(photo_rates)

            if (tuvx_time_since_last > 1.0e30_RKIND) then
                tuvx_time_since_last = 0.0_RKIND
            else
                tuvx_time_since_last = tuvx_time_since_last - tuvx_update_interval
                if (tuvx_time_since_last >= tuvx_update_interval) tuvx_time_since_last = 0.0_RKIND
            end if
        end if
        if (chemistry_interval > 0.0_RKIND) then
            chemistry_interval_clock = chemistry_interval_clock - chemistry_interval
            if (chemistry_interval_clock >= chemistry_interval) chemistry_interval_clock = 0.0_RKIND
            chemistry_accumulated_dt = 0.0_RKIND
        end if
        if (allocated(scalars_before_chem)) deallocate(scalars_before_chem)
        call release_miem_step_arrays()
        if (use_ref_solve) first_call = .false.
        chem_consecutive_failures = 0
        call mpas_timer_stop('chemistry')

    contains

        !--------------------------------------------------------------------
        !  routine release_miem_step_arrays
        !--------------------------------------------------------------------
        subroutine release_miem_step_arrays()

            if (allocated(local_emission_flux)) deallocate(local_emission_flux)
            if (allocated(local_emission_layer_flux)) deallocate(local_emission_layer_flux)
            if (allocated(local_emission_sector_flux)) deallocate(local_emission_sector_flux)
            if (allocated(local_emission_category_flux)) deallocate(local_emission_category_flux)
            if (allocated(local_emission_sector_layer_flux)) deallocate(local_emission_sector_layer_flux)
            if (allocated(local_emission_category_layer_flux)) deallocate(local_emission_category_layer_flux)

        end subroutine release_miem_step_arrays

        !--------------------------------------------------------------------
        !  routine chem_step_rollback
        !
        !> \brief Undo a failed chemistry step and stop its timer.
        !> \details
        !>  Restores the chemistry-species rows from the pre-step snapshot,
        !>  frees the photolysis-rate buffer, rewinds the TUV-x update clock,
        !>  and (when the reference solve is active) forces a reference
        !>  re-sync on the next successful step. Stops the 'chemistry' timer.
        !>  Callers issue their own `return` after this returns.
        !--------------------------------------------------------------------
        subroutine chem_step_rollback()

            integer :: i_rb

            if (allocated(scalars_before_chem) .and. associated(scalars)) then
                do i_rb = 1, n_chem_species
                    scalars(chem_species(i_rb)%mpas_index, 1:nVertLevels, 1:nCellsSolve) = &
                        scalars_before_chem(i_rb, :, :)
                end do
                deallocate(scalars_before_chem)
            end if
            if (allocated(photo_rates)) deallocate(photo_rates)
            call release_miem_step_arrays()
            if (use_ref_solve) first_call = .true.
            if (photolysis_cadence_active) tuvx_time_since_last = tuvx_time_before_step
            call mpas_timer_stop('chemistry')

        end subroutine chem_step_rollback
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
    subroutine chemistry_finalize(dminfo)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_finalize
        use mpas_miem, only: miem_finalize
        use mpas_tuvx, only: tuvx_finalize_sub => tuvx_finalize
        use mpas_prescribed_fields, only: prescribed_fields_finalize
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: dm_info

        ! dminfo is the run-level global sync point: the TUV-x and MUSICA
        ! finalizers reduce and log accumulated negative-clip totals here (the
        ! only place a collective is safe, since per-step reductions could deadlock).
        type(dm_info), intent(in) :: dminfo

#ifdef MPAS_USE_MUSICA
        call mpas_log_write('Finalizing chemistry packages...')
        call tuvx_finalize_sub(dminfo)
        call prescribed_fields_finalize()
        call miem_finalize(dminfo)
        call musica_finalize(dminfo)
#endif

    end subroutine chemistry_finalize


    !------------------------------------------------------------------------
    !  routine chemistry_set_photolysis_diag
    !
    !> \brief Write per-cell/level photolysis rates into diagnostic pool
    !> \details
    !>  For every cached photolysis rate name (j_<name>), look up the
    !>  matching diag array and copy the corresponding slice from
    !>  photo_rates. Missing diag variables are logged as critical errors
    !>  because runtime injection should have registered every rate.
    !------------------------------------------------------------------------
    subroutine chemistry_set_photolysis_diag(diag, photo_rates, nCells, nCellsSolve, nVertLevels)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: photo_rp_names, n_photo_rp
#endif
        use mpas_pool_routines, only: mpas_pool_get_array
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_log, only: mpas_log_write

        type(mpas_pool_type), intent(inout) :: diag
        real(kind=RKIND), intent(in) :: photo_rates(:,:,:)  ! (n_rates, nVertLevels, nCellsSolve)
        integer, intent(in) :: nCells, nCellsSolve, nVertLevels

#ifdef MPAS_USE_MUSICA
        real(kind=RKIND), dimension(:,:), pointer :: diag_field
        integer :: ir, k, iCell

        do ir = 1, n_photo_rp
            nullify(diag_field)
            call mpas_pool_get_array(diag, 'j_' // trim(photo_rp_names(ir)), diag_field)
            if (.not. associated(diag_field)) then
                call mpas_log_write('[Chemistry] Missing photolysis diagnostic j_' // &
                                    trim(photo_rp_names(ir)), messageType=MPAS_LOG_CRIT)
                return
            end if

            do iCell = 1, min(nCellsSolve, size(diag_field, 2))
                do k = 1, min(nVertLevels, size(diag_field, 1))
                    diag_field(k, iCell) = photo_rates(ir, k, iCell)
                end do
            end do
        end do
#endif

    end subroutine chemistry_set_photolysis_diag


    !------------------------------------------------------------------------
    !  routine zero_photolysis_diag
    !
    !> \brief Zero cached photolysis diagnostic fields at init time
    !> \details
    !>  Clears every runtime-injected j_<reaction> diagnostic corresponding
    !>  to a cached MICM photolysis rate parameter.
    !------------------------------------------------------------------------
    subroutine zero_photolysis_diag(diag)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: photo_rp_names, n_photo_rp
#endif
        use mpas_pool_routines, only: mpas_pool_get_array
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_log, only: mpas_log_write

        type(mpas_pool_type), intent(inout) :: diag

#ifdef MPAS_USE_MUSICA
        real(kind=RKIND), dimension(:,:), pointer :: diag_field
        integer :: ir

        do ir = 1, n_photo_rp
            nullify(diag_field)
            call mpas_pool_get_array(diag, 'j_' // trim(photo_rp_names(ir)), diag_field)
            if (.not. associated(diag_field)) then
                call mpas_log_write('[Chemistry] Missing photolysis diagnostic j_' // &
                                    trim(photo_rp_names(ir)), messageType=MPAS_LOG_CRIT)
                return
            end if
            diag_field = 0.0_RKIND
        end do
#endif

    end subroutine zero_photolysis_diag


    !------------------------------------------------------------------------
    !  routine chemistry_set_emissions_diag
    !
    !> \brief Publish total, layered, and selected-group MIEM fluxes
    !> \details
    !>  Only owned cells are changed. Halo cells retain their prior values,
    !>  and this routine is called only after the chemistry transaction has
    !>  completed successfully.
    !------------------------------------------------------------------------
    subroutine chemistry_set_emissions_diag(diag, local_flux, local_layer_flux, &
                                             local_sector_flux, local_category_flux, &
                                             local_sector_layer_flux, local_category_layer_flux, &
                                             nCellsSolve, nVertLevels, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use iso_fortran_env, only: real64
        use mpas_miem, only: miem_is_enabled, miem_get_species, &
                             miem_get_diagnostic_selection
#endif
        use mpas_pool_routines, only: mpas_pool_get_array
        use mpas_derived_types, only: mpas_pool_type

        type(mpas_pool_type), intent(inout) :: diag
#ifdef MPAS_USE_MUSICA
        real(real64), intent(in) :: local_flux(:,:)
        real(real64), intent(in) :: local_layer_flux(:,:,:)
        real(real64), intent(in) :: local_sector_flux(:,:,:)
        real(real64), intent(in) :: local_category_flux(:,:,:)
        real(real64), intent(in) :: local_sector_layer_flux(:,:,:,:)
        real(real64), intent(in) :: local_category_layer_flux(:,:,:,:)
#else
        real(kind=RKIND), intent(in) :: local_flux(:,:)
        real(kind=RKIND), intent(in) :: local_layer_flux(:,:,:)
        real(kind=RKIND), intent(in) :: local_sector_flux(:,:,:)
        real(kind=RKIND), intent(in) :: local_category_flux(:,:,:)
        real(kind=RKIND), intent(in) :: local_sector_layer_flux(:,:,:,:)
        real(kind=RKIND), intent(in) :: local_category_layer_flux(:,:,:,:)
#endif
        integer, intent(in) :: nCellsSolve, nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        character(len=64), allocatable :: emission_names(:)
        character(len=64), allocatable :: sector_names(:), sector_labels(:)
        integer, allocatable :: category_ids(:)
        character(len=128) :: field_name
        character(len=32) :: category_text
        integer :: n_emission_species, n_sectors, n_categories
        integer :: i, group
        logical :: layered_diagnostics

        error_code = 0
        error_message = ''
        if (.not. miem_is_enabled()) return

        call miem_get_species(n_emission_species, emission_names)
        call miem_get_diagnostic_selection(n_sectors, sector_names, sector_labels, &
                                           n_categories, category_ids, layered_diagnostics)
        if (size(local_flux, 1) /= n_emission_species .or. &
            size(local_flux, 2) /= nCellsSolve .or. &
            size(local_layer_flux, 1) /= n_emission_species .or. &
            size(local_layer_flux, 2) /= nVertLevels .or. &
            size(local_layer_flux, 3) /= nCellsSolve .or. &
            size(local_sector_flux, 1) /= n_emission_species .or. &
            size(local_sector_flux, 2) /= n_sectors .or. &
            size(local_sector_flux, 3) /= nCellsSolve .or. &
            size(local_category_flux, 1) /= n_emission_species .or. &
            size(local_category_flux, 2) /= n_categories .or. &
            size(local_category_flux, 3) /= nCellsSolve) then
            error_code = 1
            error_message = '[MIEM] emissions diagnostic dimensions do not match the sampled inventory'
            return
        end if
        if (layered_diagnostics) then
            if (size(local_sector_layer_flux, 1) /= n_emission_species .or. &
                size(local_sector_layer_flux, 2) /= nVertLevels .or. &
                size(local_sector_layer_flux, 3) /= n_sectors .or. &
                size(local_sector_layer_flux, 4) /= nCellsSolve .or. &
                size(local_category_layer_flux, 1) /= n_emission_species .or. &
                size(local_category_layer_flux, 2) /= nVertLevels .or. &
                size(local_category_layer_flux, 3) /= n_categories .or. &
                size(local_category_layer_flux, 4) /= nCellsSolve) then
                error_code = 1
                error_message = '[MIEM] layered diagnostic dimensions do not match the sampled inventory'
                return
            end if
        end if

        ! Validate the complete destination set before changing any field.
        do i = 1, n_emission_species
            field_name = 'emis_' // trim(emission_names(i))
            call require_emissions_diag_1d(trim(field_name))
            if (error_code /= 0) return
            if (layered_diagnostics) then
                field_name = 'emis_' // trim(emission_names(i)) // '_layer'
                call require_emissions_diag_2d(trim(field_name))
                if (error_code /= 0) return
            end if
            do group = 1, n_sectors
                field_name = 'emis_' // trim(emission_names(i)) // '__sector_' // &
                             trim(sector_labels(group))
                call require_emissions_diag_1d(trim(field_name))
                if (error_code /= 0) return
                if (layered_diagnostics) then
                    field_name = 'emis_' // trim(emission_names(i)) // '__sector_' // &
                                 trim(sector_labels(group)) // '_layer'
                    call require_emissions_diag_2d(trim(field_name))
                    if (error_code /= 0) return
                end if
            end do
            do group = 1, n_categories
                write(category_text, '(i0)') category_ids(group)
                field_name = 'emis_' // trim(emission_names(i)) // '__category_' // &
                             trim(category_text)
                call require_emissions_diag_1d(trim(field_name))
                if (error_code /= 0) return
                if (layered_diagnostics) then
                    field_name = 'emis_' // trim(emission_names(i)) // '__category_' // &
                                 trim(category_text) // '_layer'
                    call require_emissions_diag_2d(trim(field_name))
                    if (error_code /= 0) return
                end if
            end do
        end do

        do i = 1, n_emission_species
            field_name = 'emis_' // trim(emission_names(i))
            call write_emissions_diag_1d(trim(field_name), local_flux(i, :))
            if (layered_diagnostics) then
                field_name = 'emis_' // trim(emission_names(i)) // '_layer'
                call write_emissions_diag_2d(trim(field_name), local_layer_flux(i, :, :))
            end if
            do group = 1, n_sectors
                field_name = 'emis_' // trim(emission_names(i)) // '__sector_' // &
                             trim(sector_labels(group))
                call write_emissions_diag_1d(trim(field_name), local_sector_flux(i, group, :))
                if (layered_diagnostics) then
                    field_name = 'emis_' // trim(emission_names(i)) // '__sector_' // &
                                 trim(sector_labels(group)) // '_layer'
                    call write_emissions_diag_2d(trim(field_name), &
                                                 local_sector_layer_flux(i, :, group, :))
                end if
            end do
            do group = 1, n_categories
                write(category_text, '(i0)') category_ids(group)
                field_name = 'emis_' // trim(emission_names(i)) // '__category_' // &
                             trim(category_text)
                call write_emissions_diag_1d(trim(field_name), local_category_flux(i, group, :))
                if (layered_diagnostics) then
                    field_name = 'emis_' // trim(emission_names(i)) // '__category_' // &
                                 trim(category_text) // '_layer'
                    call write_emissions_diag_2d(trim(field_name), &
                                                 local_category_layer_flux(i, :, group, :))
                end if
            end do
        end do
#else
        error_code = 0
        error_message = ''
#endif

#ifdef MPAS_USE_MUSICA
    contains

        subroutine require_emissions_diag_1d(name)

            character(len=*), intent(in) :: name
            real(kind=RKIND), dimension(:), pointer :: field

            nullify(field)
            call mpas_pool_get_array(diag, trim(name), field)
            if (.not. associated(field)) then
                error_code = 1
                error_message = '[MIEM] missing emissions diagnostic ' // trim(name)
            else if (size(field) < nCellsSolve) then
                error_code = 1
                error_message = '[MIEM] undersized emissions diagnostic ' // trim(name)
            end if

        end subroutine require_emissions_diag_1d

        subroutine require_emissions_diag_2d(name)

            character(len=*), intent(in) :: name
            real(kind=RKIND), dimension(:,:), pointer :: field

            nullify(field)
            call mpas_pool_get_array(diag, trim(name), field)
            if (.not. associated(field)) then
                error_code = 1
                error_message = '[MIEM] missing emissions diagnostic ' // trim(name)
            else if (size(field, 1) < nVertLevels .or. size(field, 2) < nCellsSolve) then
                error_code = 1
                error_message = '[MIEM] undersized emissions diagnostic ' // trim(name)
            end if

        end subroutine require_emissions_diag_2d

        subroutine write_emissions_diag_1d(name, values)

            character(len=*), intent(in) :: name
            real(real64), intent(in) :: values(:)
            real(kind=RKIND), dimension(:), pointer :: field

            nullify(field)
            call mpas_pool_get_array(diag, trim(name), field)
            field(1:nCellsSolve) = real(values(1:nCellsSolve), RKIND)

        end subroutine write_emissions_diag_1d

        subroutine write_emissions_diag_2d(name, values)

            character(len=*), intent(in) :: name
            real(real64), intent(in) :: values(:,:)
            real(kind=RKIND), dimension(:,:), pointer :: field

            nullify(field)
            call mpas_pool_get_array(diag, trim(name), field)
            field(1:nVertLevels, 1:nCellsSolve) = &
                real(values(1:nVertLevels, 1:nCellsSolve), RKIND)

        end subroutine write_emissions_diag_2d
#endif

    end subroutine chemistry_set_emissions_diag


    !------------------------------------------------------------------------
    !  routine zero_emissions_diag
    !
    !> \brief Zero all runtime MIEM diagnostics before initial output
    !------------------------------------------------------------------------
    subroutine zero_emissions_diag(diag, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_miem, only: miem_is_enabled, miem_get_species, &
                             miem_get_diagnostic_selection
#endif
        use mpas_pool_routines, only: mpas_pool_get_array
        use mpas_derived_types, only: mpas_pool_type

        type(mpas_pool_type), intent(inout) :: diag
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        character(len=64), allocatable :: emission_names(:)
        character(len=64), allocatable :: sector_names(:), sector_labels(:)
        integer, allocatable :: category_ids(:)
        character(len=128) :: field_name
        character(len=32) :: category_text
        integer :: n_emission_species, n_sectors, n_categories
        integer :: i, group
        logical :: layered_diagnostics

        error_code = 0
        error_message = ''
        if (.not. miem_is_enabled()) return

        call miem_get_species(n_emission_species, emission_names)
        call miem_get_diagnostic_selection(n_sectors, sector_names, sector_labels, &
                                           n_categories, category_ids, layered_diagnostics)
        do i = 1, n_emission_species
            field_name = 'emis_' // trim(emission_names(i))
            call zero_emissions_diag_1d(trim(field_name))
            if (error_code /= 0) return
            if (layered_diagnostics) then
                field_name = 'emis_' // trim(emission_names(i)) // '_layer'
                call zero_emissions_diag_2d(trim(field_name))
                if (error_code /= 0) return
            end if
            do group = 1, n_sectors
                field_name = 'emis_' // trim(emission_names(i)) // '__sector_' // &
                             trim(sector_labels(group))
                call zero_emissions_diag_1d(trim(field_name))
                if (error_code /= 0) return
                if (layered_diagnostics) then
                    field_name = 'emis_' // trim(emission_names(i)) // '__sector_' // &
                                 trim(sector_labels(group)) // '_layer'
                    call zero_emissions_diag_2d(trim(field_name))
                    if (error_code /= 0) return
                end if
            end do
            do group = 1, n_categories
                write(category_text, '(i0)') category_ids(group)
                field_name = 'emis_' // trim(emission_names(i)) // '__category_' // &
                             trim(category_text)
                call zero_emissions_diag_1d(trim(field_name))
                if (error_code /= 0) return
                if (layered_diagnostics) then
                    field_name = 'emis_' // trim(emission_names(i)) // '__category_' // &
                                 trim(category_text) // '_layer'
                    call zero_emissions_diag_2d(trim(field_name))
                    if (error_code /= 0) return
                end if
            end do
        end do
#else
        error_code = 0
        error_message = ''
#endif

#ifdef MPAS_USE_MUSICA
    contains

        subroutine zero_emissions_diag_1d(name)

            character(len=*), intent(in) :: name
            real(kind=RKIND), dimension(:), pointer :: field

            nullify(field)
            call mpas_pool_get_array(diag, trim(name), field)
            if (.not. associated(field)) then
                error_code = 1
                error_message = '[MIEM] missing emissions diagnostic ' // trim(name)
                return
            end if
            field = 0.0_RKIND

        end subroutine zero_emissions_diag_1d

        subroutine zero_emissions_diag_2d(name)

            character(len=*), intent(in) :: name
            real(kind=RKIND), dimension(:,:), pointer :: field

            nullify(field)
            call mpas_pool_get_array(diag, trim(name), field)
            if (.not. associated(field)) then
                error_code = 1
                error_message = '[MIEM] missing emissions diagnostic ' // trim(name)
                return
            end if
            field = 0.0_RKIND

        end subroutine zero_emissions_diag_2d
#endif

    end subroutine zero_emissions_diag


    !------------------------------------------------------------------------
    !  routine chemistry_from_MPAS
    !
    !> \brief Extracts MPAS state for chemistry processing
    !> \author CheMPAS-A Developers
    !> \date   January 2026
    !> \details
    !>  This routine extracts MPAS state variables (scalars, thermodynamic
    !>  fields) and passes them to the chemistry package (MICM) for processing.
    !>  Tracer unit conversion uses dry density because MPAS scalars are
    !>  dry-air mass mixing ratios. MICM third-body density uses total moist
    !>  molar density from p/(R*T), computed in MICM_from_chemistry.
    !------------------------------------------------------------------------
    subroutine chemistry_from_MPAS(env, state, time_lev, nCells, nCellsSolve, nVertLevels, &
                                    error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: MICM_from_chemistry
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_array

        type (chem_env_t), intent(in) :: env
        type (mpas_pool_type), intent(in) :: state
        integer, intent(in) :: time_lev
        integer, intent(in) :: nCells
        integer, intent(in) :: nCellsSolve
        integer, intent(in) :: nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        ! MPAS state arrays
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars

        error_code = 0
        error_message = ''

        ! Get MPAS arrays from pools
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

        ! Validate pointers
        if (.not. associated(scalars) .or. .not. allocated(env%rho_dry) .or. &
            .not. allocated(env%temperature) .or. .not. allocated(env%pressure)) then
            error_code = 1
            error_message = 'Required MPAS arrays not available in chemistry_from_MPAS.'
            return
        end if

        ! Call MICM-specific coupling routine (uses chem_species table internally).
        ! Negative-clip totals are accumulated rank-locally and reduced once at
        ! finalize, so no collective is issued here (this call sits downstream of
        ! rank-local recoverable early returns in chemistry_step).
        call MICM_from_chemistry(scalars, env%rho_dry, env%temperature, env%pressure, &
                                  nCellsSolve, nVertLevels, &
                                  error_code, error_message)
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
    !>  after the chemistry solver has completed. Tracer unit conversion uses
    !>  dry density because MPAS scalars are dry-air mass mixing ratios.
    !------------------------------------------------------------------------
    subroutine chemistry_to_MPAS(env, state, time_lev, nCells, nCellsSolve, nVertLevels, &
                                  error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: MICM_to_chemistry
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_array

        type (chem_env_t), intent(in) :: env
        type (mpas_pool_type), intent(inout) :: state
        integer, intent(in) :: time_lev
        integer, intent(in) :: nCells
        integer, intent(in) :: nCellsSolve
        integer, intent(in) :: nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        ! MPAS state arrays
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars

        error_code = 0
        error_message = ''

        ! Get MPAS arrays from pools
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

        ! Validate pointers
        if (.not. associated(scalars) .or. .not. allocated(env%rho_dry)) then
            error_code = 1
            error_message = 'Required MPAS arrays not available in chemistry_to_MPAS.'
            return
        end if

        ! Call MICM-specific coupling routine (uses chem_species table internally).
        ! Negative-clip totals are accumulated rank-locally and reduced once at
        ! finalize, so no collective is issued here.
        call MICM_to_chemistry(scalars, env%rho_dry, nCellsSolve, nVertLevels, &
                                error_code, error_message)
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
    subroutine chemistry_seed_chem(dimensions, structs, dminfo)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: micm_to_mpas_chem, n_chem_species, chem_species
#endif
        use mpas_log, only : mpas_log_write
        use mpas_dmpar, only: mpas_dmpar_max_real, mpas_dmpar_min_real
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT, field3dReal, dm_info
        use mpas_pool_routines, only: mpas_pool_get_dimension, mpas_pool_get_subpool, &
                                      mpas_pool_get_field, mpas_pool_get_array

        type (mpas_pool_type), intent(in) :: dimensions
        type (mpas_pool_type), intent(in) :: structs
        type (dm_info), pointer :: dminfo

#ifdef MPAS_USE_MUSICA
        integer                         :: nCells
        integer                         :: nCellsSolve
        integer                         :: nVertLevels
        integer, pointer                :: nCells_ptr
        integer, pointer                :: nCellsSolve_ptr
        integer, pointer                :: nVertLevels_ptr
        type (mpas_pool_type), pointer  :: statePool
        type (field3dReal), pointer     :: scalarsField
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars_time1
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars_time2
        real (kind=RKIND)               :: local_max, local_min
        real (kind=RKIND)               :: global_max, global_min
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

        call mpas_pool_get_dimension(dimensions, 'nCellsSolve', nCellsSolve_ptr)
        if (associated(nCellsSolve_ptr)) then
            nCellsSolve = nCellsSolve_ptr
        else
            nCellsSolve = nCells
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

        ! Skip seeding if any chemistry tracer already has spatial gradients.
        ! Host-bound species (e.g. H2O -> qv) are excluded: their host field
        ! (water vapor) always has gradients, which would otherwise wrongly
        ! suppress seeding of the genuine chemistry tracers.
        do i_sp = 1, n_chem_species
            if (chem_species(i_sp)%host_bound) cycle
            mpas_idx = chem_species(i_sp)%mpas_index
            if (nCellsSolve > 0) then
                local_max = maxval(scalars_time1(mpas_idx,:,1:nCellsSolve))
                local_min = minval(scalars_time1(mpas_idx,:,1:nCellsSolve))
            else
                local_max = -huge(local_max)
                local_min = huge(local_min)
            end if
            call mpas_dmpar_max_real(dminfo, local_max, global_max)
            call mpas_dmpar_min_real(dminfo, local_min, global_min)
            ! Skip MICM seeding when a tracer already carries initial values —
            ! either spatial structure (a gradient) OR a uniform nonzero
            ! background (e.g. the LNOx qO3 8.276e-8 field). Seeding must only
            ! fill tracers that are genuinely empty (all ~zero over owned
            ! cells); a gradient-only test wrongly treats a uniform nonzero
            ! file-initialized field as uninitialized and overwrites it with the
            ! mechanism's (often zero) initial concentration. (FABLE R11)
            if (global_max - global_min > 1.0e-10_RKIND .or. &
                max(abs(global_max), abs(global_min)) > 1.0e-30_RKIND) then
                call mpas_log_write('[MUSICA] Tracer ' // trim(chem_species(i_sp)%mpas_name) // &
                    ' already initialized (nonzero or nonuniform); skipping MICM seed.')
                return
            end if
        end do

        error_code = 0
        call micm_to_mpas_chem(scalars_time1, nCellsSolve, nVertLevels, &
                               error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if

        error_code = 0
        call micm_to_mpas_chem(scalars_time2, nCellsSolve, nVertLevels, &
                               error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
        end if
#endif

    end subroutine chemistry_seed_chem


    !------------------------------------------------------------------------
    !  subroutine chemistry_query_species
    !
    !> \brief Wrapper to query MICM config for chemistry species names
    !> \details
    !>  Delegates to musica_query_species. Provides a clean interface
    !>  for mpas_atm_core_interface without requiring direct musica access.
    !------------------------------------------------------------------------
    subroutine chemistry_query_species(config_path, num_species, micm_names, mpas_names, &
            molar_masses, host_bound, fixed_vmrs, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_query_species
#endif
        use iso_fortran_env, only: real64

        character(len=*), intent(in)                :: config_path
        integer, intent(out)                        :: num_species
        character(len=64), allocatable, intent(out)  :: micm_names(:)
        character(len=64), allocatable, intent(out)  :: mpas_names(:)
        real(real64), allocatable, intent(out)       :: molar_masses(:)
        logical, allocatable, intent(out)            :: host_bound(:)
        real(real64), allocatable, intent(out)       :: fixed_vmrs(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        error_code = 0
        error_message = ''
        num_species = 0

#ifdef MPAS_USE_MUSICA
        call musica_query_species(config_path, num_species, micm_names, mpas_names, &
            molar_masses, host_bound, fixed_vmrs, error_code, error_message)
#else
        allocate(micm_names(0), mpas_names(0), molar_masses(0), host_bound(0), fixed_vmrs(0))
#endif

    end subroutine chemistry_query_species


    !------------------------------------------------------------------------
    !  subroutine chemistry_query_emission_species
    !
    !> \brief Wrapper to query EMIS.* names from the MICM configuration
    !------------------------------------------------------------------------
    subroutine chemistry_query_emission_species(config_path, num_species, species_names, &
                                                error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_query_emission_species
#endif

        character(len=*), intent(in) :: config_path
        integer, intent(out) :: num_species
        character(len=64), allocatable, intent(out) :: species_names(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        error_code = 0
        error_message = ''
        num_species = 0

#ifdef MPAS_USE_MUSICA
        call musica_query_emission_species(config_path, num_species, species_names, &
                                           error_code, error_message)
#else
        allocate(species_names(0))
#endif

    end subroutine chemistry_query_emission_species


    !------------------------------------------------------------------------
    !  subroutine chemistry_query_photo_rates
    !
    !> \brief Wrapper to query TUV-x config for photolysis rate names
    !> \details
    !>  Delegates to mpas_tuvx without exposing the chemistry setup path to
    !>  the TUV-x module directly.
    !------------------------------------------------------------------------
    subroutine chemistry_query_photo_rates(config_path, num_rates, rate_names, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_tuvx, only: tuvx_query_photo_rates
#endif

        character(len=*), intent(in)                :: config_path
        integer, intent(out)                        :: num_rates
        character(len=64), allocatable, intent(out)  :: rate_names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        error_code = 0
        error_message = ''
        num_rates = 0

#ifdef MPAS_USE_MUSICA
        call tuvx_query_photo_rates(config_path, num_rates, rate_names, error_code, error_message)
#else
        allocate(rate_names(0))
#endif

    end subroutine chemistry_query_photo_rates


    !------------------------------------------------------------------------
    !  subroutine chemistry_lightning_nox_init
    !
    !> \brief Initialize lightning NOx source from namelist + tracer pool
    !------------------------------------------------------------------------
    subroutine chemistry_lightning_nox_init(configs, statePool)

        use mpas_lightning_nox, only: lightning_nox_init
        use mpas_derived_types, only: mpas_pool_type

        type(mpas_pool_type), intent(in) :: configs
        type(mpas_pool_type), intent(in) :: statePool

        call lightning_nox_init(configs, statePool)

    end subroutine chemistry_lightning_nox_init

end module mpas_atm_chemistry
