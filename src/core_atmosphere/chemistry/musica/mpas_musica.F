! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_musica
!
!> \brief Manages interactions with the MUSICA chemistry package
!> \author CheMPAS-A Developers
!> \date   20 August 2025
!> \details
!>  This module manages the interactions with the MUSICA chemistry package,
!>  including initialization, time-stepping, and finalization for
!>  MICM (ODE solver) and TUV-x (photolysis rate constant calculator).
!--------------------------------------------------------------------------
module mpas_musica

    use musica_micm,  only: micm_t, get_micm_version
    ! get_micm_version is here to avoid an ICE when within musica_init
    use musica_state, only: state_t

    implicit none

    private

    public :: musica_init, musica_step, musica_step_ref, musica_finalize, micm_to_mpas_chem
    public :: MICM_from_chemistry, MICM_to_chemistry
    public :: log_column_comparison, copy_state_to_ref
    public :: resolve_mpas_indices, n_chem_species, chem_species
    public :: musica_query_species
    public :: musica_cache_photo_indices, musica_set_photolysis_rates
    public :: n_photo_rp, photo_rp_names

    type(micm_t),  pointer :: micm  => null ( )  ! Pointer to the MICM ODE solver instance
    type(state_t), pointer :: state => null ( )  ! Pointer to the state of the MICM solver (coupled to MPAS)
    type(state_t), pointer :: state_ref => null ( )  ! Reference state (uncoupled, chemistry only)

    ! Photolysis rate-parameter bookkeeping. For each TUV-x reaction name we
    ! resolve the MICM rate_parameters_ordering key PHOTO.<name> once at
    ! init, then reuse the cached stride index on every step.
    integer, save                    :: n_photo_rp = 0
    integer,             allocatable :: photo_rp_indices(:)
    character(len=64),   allocatable :: photo_rp_names(:)

    type :: species_entry_t
        character(len=64) :: micm_name    ! e.g. 'AB'
        character(len=64) :: mpas_name    ! e.g. 'qAB'
        integer           :: micm_index   ! position in MICM species_ordering
        integer           :: mpas_index   ! position in MPAS scalars (1st dim)
        real(kind=8)      :: molar_mass   ! kg/mol from MICM config
    end type

    integer :: n_chem_species = 0
    type(species_entry_t), allocatable :: chem_species(:)
    real(kind=8), parameter :: M_AIR = 0.0289644_8  ! kg/mol (dry air molar mass)

    contains

    !------------------------------------------------------------------------
    !  subroutine musica_init
    !
    !> \brief Initializes the MUSICA chemistry package
    !> \author CheMPAS-A Developers
    !> \date   20 August 2025
    !> \details
    !>  This subroutine initializes the MUSICA chemistry package,
    !>  setting up necessary parameters and data structures for the simulation.
    !>  For now, we will load fixed configurations for MICM and TUV-x.
    !>  Later, we will gradually replace the fixed configuration elements
    !>  with runtime updates using the MUSICA API.
    !------------------------------------------------------------------------
    subroutine musica_init(filename_of_micm_configuration, &
            number_of_grid_cells, &
            nox_tau, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND

        use musica_micm, only : BackwardEulerStandardOrder
        use musica_util, only : error_t, string_t

        use mpas_log, only : mpas_log_write

        character(len=*),              intent(in)  :: filename_of_micm_configuration
        integer,                       intent(in)  :: number_of_grid_cells
        real(kind=RKIND),              intent(in)  :: nox_tau           ! [s] NOx sink timescale (0 = disabled)
        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=*), parameter :: INITIAL_CONCENTRATION_PROPERTY = "__initial concentration"
        character(len=*), parameter :: MOLAR_MASS_PROPERTY = "__molar mass"

        real(real64), allocatable :: init_conc(:)
        real(real64) :: molar_mass_val

        type(error_t)  :: error
        type(string_t) :: micm_version
        ! Try Backward Euler for chapman_nox to test robustness; see
        ! docs/plans/2026-04-18-chapman-nox-chem-box-issue.md.
        integer :: solver_type = BackwardEulerStandardOrder
        integer :: i_species, spec_idx
        character(len=:), allocatable :: sp_name

        micm_version = get_micm_version()

        call mpas_log_write('[MUSICA] Initializing MICM chemistry package...')
        call mpas_log_write('MICM version: ' // micm_version%value_)
        call mpas_log_write('MICM number of grid cells: $i', intArgs=[number_of_grid_cells])

        micm => micm_t(trim(filename_of_micm_configuration), solver_type, error)
        if (has_error_occurred(error, error_message, error_code)) return

        state => micm%get_state(number_of_grid_cells, error)
        if (has_error_occurred(error, error_message, error_code)) return

        ! Create reference state (uncoupled, for advection diagnostics)
        state_ref => micm%get_state(number_of_grid_cells, error)
        if (has_error_occurred(error, error_message, error_code)) return
        call mpas_log_write('[MUSICA] Created reference state for advection diagnostics')

        do i_species = 1, state%species_ordering%size()
            call mpas_log_write('[MUSICA] MICM species: ' // state%species_ordering%name(i_species))
        end do

        call assign_rate_parameters(state, number_of_grid_cells, nox_tau)
        call assign_rate_parameters(state_ref, number_of_grid_cells, nox_tau)

        ! Build species table from MICM config
        n_chem_species = state%species_ordering%size()
        if (allocated(chem_species)) deallocate(chem_species)
        allocate(chem_species(n_chem_species))
        allocate(init_conc(number_of_grid_cells))

        do i_species = 1, n_chem_species
            sp_name = trim(state%species_ordering%name(i_species))
            chem_species(i_species)%micm_name = sp_name
            chem_species(i_species)%mpas_name = 'q' // sp_name
            chem_species(i_species)%mpas_index = 0  ! resolved later by resolve_mpas_indices

            ! Get the MICM stride-based index (may differ from iteration order)
            spec_idx = state%species_ordering%index(sp_name, error)
            if (has_error_occurred(error, error_message, error_code)) return
            chem_species(i_species)%micm_index = spec_idx

            ! Read molar mass from MICM species properties
            molar_mass_val = micm%get_species_property_double(sp_name, MOLAR_MASS_PROPERTY, error)
            if (has_error_occurred(error, error_message, error_code)) return
            chem_species(i_species)%molar_mass = molar_mass_val

            call mpas_log_write('[MUSICA] Species ' // trim(sp_name) // &
                ' molar_mass=$r kg/mol', realArgs=[real(chem_species(i_species)%molar_mass, kind=RKIND)])

            ! Read initial concentration and seed both states
            init_conc = micm%get_species_property_double(sp_name, INITIAL_CONCENTRATION_PROPERTY, error)
            if (has_error_occurred(error, error_message, error_code)) return

            call set_species_profile(state, sp_name, init_conc, error)
            if (has_error_occurred(error, error_message, error_code)) return
            call set_species_profile(state_ref, sp_name, init_conc, error)
            if (has_error_occurred(error, error_message, error_code)) return
        end do

        ! Photolysis rate-parameter indices are cached later, after TUV-x
        ! init, by chemistry_init calling musica_cache_photo_indices with
        ! the reaction names TUV-x enumerated.

        deallocate(init_conc)

    end subroutine musica_init

    !------------------------------------------------------------------------
    !  subroutine resolve_mpas_indices
    !
    !> \brief Resolves MPAS scalar indices for all chemistry species
    !> \details
    !>  Loops over chem_species and looks up each species' MPAS tracer index
    !>  from the state pool dimension 'index_<mpas_name>'.
    !------------------------------------------------------------------------
    subroutine resolve_mpas_indices(statePool, error_code, error_message)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_dimension
        use mpas_log, only: mpas_log_write

        type(mpas_pool_type), intent(in) :: statePool
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: i_sp
        integer, pointer :: idx_ptr

        error_code = 0
        error_message = ''

        do i_sp = 1, n_chem_species
            nullify(idx_ptr)
            call mpas_pool_get_dimension(statePool, 'index_' // trim(chem_species(i_sp)%mpas_name), idx_ptr)

            if (.not. associated(idx_ptr)) then
                error_code = 1
                error_message = '[MUSICA] No MPAS tracer index for ' // trim(chem_species(i_sp)%mpas_name)
                call mpas_log_write(error_message)
                return
            end if

            if (idx_ptr < 1) then
                error_code = 1
                error_message = '[MUSICA] Invalid MPAS tracer index for ' // trim(chem_species(i_sp)%mpas_name)
                call mpas_log_write(error_message)
                return
            end if

            chem_species(i_sp)%mpas_index = idx_ptr
            call mpas_log_write('[MUSICA] Resolved ' // trim(chem_species(i_sp)%micm_name) // &
                ' -> MPAS index $i', intArgs=[chem_species(i_sp)%mpas_index])
        end do

    end subroutine resolve_mpas_indices

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

        use musica_micm,  only: RosenbrockStandardOrder
        use musica_state, only: state_t
        use musica_util,  only: error_t
        use mpas_log,     only: mpas_log_write

        character(len=*), intent(in)                :: config_path
        integer, intent(out)                        :: num_species
        character(len=64), allocatable, intent(out)  :: species_names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        type(micm_t), pointer  :: tmp_micm => null()
        type(state_t), pointer :: tmp_state => null()
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
        if (has_error_occurred(error, error_message, error_code)) goto 100

        tmp_state => tmp_micm%get_state(1, error)
        if (has_error_occurred(error, error_message, error_code)) goto 100

        num_species = tmp_state%species_ordering%size()
        allocate(species_names(num_species))

        do i = 1, num_species
            species_names(i) = 'q' // trim(tmp_state%species_ordering%name(i))
            call mpas_log_write('[MUSICA] Discovered chemistry tracer: ' // trim(species_names(i)))
        end do

100     continue
        if (associated(tmp_state)) deallocate(tmp_state)
        if (associated(tmp_micm)) deallocate(tmp_micm)

    end subroutine musica_query_species

    !------------------------------------------------------------------------
    !  subroutine musica_step
    !
    !> \brief Steps the MUSICA chemistry package
    !> \author David Fillmore
    !> \date   17 November 2025
    !> \details
    !>  This subroutine steps the MICM ODE solver.
    !------------------------------------------------------------------------
    subroutine musica_step(time_step, error_code, error_message)

        use iso_fortran_env, only: real64

        use musica_micm, only : solver_stats_t
        use musica_util, only : error_t, string_t

        use mpas_log, only : mpas_log_write

        real(real64), intent(in) :: time_step

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(string_t)       :: solver_state 
        type(solver_stats_t) :: solver_stats
        type(error_t)        :: error

        call mpas_log_write('[MUSICA] Stepping MICM solver...')

        call micm%solve(time_step, state, solver_state, solver_stats, error)
        if (has_error_occurred(error, error_message, error_code)) return

        if (solver_state%get_char_array() /= 'Converged') then
            call mpas_log_write('[MUSICA Warning] MICM solver failure')
        end if

        call log_solver_stats(solver_stats)

    end subroutine musica_step

    !------------------------------------------------------------------------
    !  subroutine musica_step_ref
    !
    !> \brief Steps the reference (uncoupled) MICM state
    !> \details
    !>  This subroutine steps the reference MICM state that is NOT coupled
    !>  to MPAS. Used for advection diagnostics - comparing this to the
    !>  coupled state shows the effect of advection.
    !------------------------------------------------------------------------
    subroutine musica_step_ref(time_step, error_code, error_message)

        use iso_fortran_env, only: real64

        use musica_micm, only : solver_stats_t
        use musica_util, only : error_t, string_t

        real(real64), intent(in) :: time_step

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(string_t)       :: solver_state
        type(solver_stats_t) :: solver_stats
        type(error_t)        :: error

        ! Solve reference state (not coupled to MPAS, chemistry only)
        call micm%solve(time_step, state_ref, solver_state, solver_stats, error)
        if (has_error_occurred(error, error_message, error_code)) return

        if (solver_state%get_char_array() /= 'Converged') then
            error_code = 1
            error_message = '[MUSICA Warning] Reference MICM solver failure'
        end if

    end subroutine musica_step_ref

    !------------------------------------------------------------------------
    !  subroutine musica_finalize
    !
    !> \brief Finalizes the MUSICA chemistry package
    !> \author CheMPAS-A Developers
    !> \date   20 August 2025
    !> \details
    !>  This subroutine finalizes the MUSICA chemistry package,
    !>  cleaning up any resources and data structures used during the simulation.
    !-------------------------------------------------------------------------
    subroutine musica_finalize()

        use mpas_log, only : mpas_log_write

        call mpas_log_write('[MUSICA] Finalizing MICM chemistry package...')

    end subroutine musica_finalize

    !------------------------------------------------------------------------
    !  subroutine micm_to_mpas_chem
    !
    !> \brief One-time copy of MICM initial state into MPAS tracers
    !> \details
    !>  Loops over chem_species to seed all chemistry tracers generically.
    !------------------------------------------------------------------------
    subroutine micm_to_mpas_chem(scalars, nCells, nVertLevels, &
            error_code, error_message)

        use mpas_kind_types, only : RKIND
        use musica_util, only : error_t
        use mpas_log, only : mpas_log_write

        real (kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: cell_stride, var_stride
        integer :: i_sp, k, iCell
        integer :: micm_levels, target_cells, target_levels
        integer :: offset
        real (kind=RKIND) :: conc_val
        logical :: logged_first
        real (kind=RKIND) :: q_min, q_max

        error_code = 0
        error_message = ''
        logged_first = .false.

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized; cannot seed MPAS tracers.'
            call mpas_log_write(error_message)
            return
        end if

        call mpas_log_write('[MUSICA] Seeding MPAS tracers from MICM initial state.')

        ! Validate all MPAS indices
        do i_sp = 1, n_chem_species
            if (chem_species(i_sp)%mpas_index < 1 .or. &
                chem_species(i_sp)%mpas_index > size(scalars, 1)) then
                error_code = 1
                error_message = '[MUSICA] Tracer index not valid for ' // trim(chem_species(i_sp)%mpas_name)
                call mpas_log_write(error_message)
                return
            end if
        end do

        if (nVertLevels > 0) then
            target_levels = min(nVertLevels, size(scalars, 2))
        else
            target_levels = size(scalars, 2)
        end if
        if (nCells > 0) then
            target_cells = min(nCells, size(scalars, 3))
        else
            target_cells = size(scalars, 3)
        end if
        micm_levels = min(state%number_of_grid_cells, target_levels)

        if (target_cells <= 0 .or. micm_levels <= 0) then
            error_code = 1
            error_message = '[MUSICA] No target cells/levels available for MICM -> MPAS seeding.'
            call mpas_log_write(error_message)
            return
        end if

        call mpas_log_write('[MUSICA] MICM->MPAS seed: micm_cells=$i target_levels=$i target_cells=$i', &
                intArgs=[state%number_of_grid_cells, micm_levels, target_cells])

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do i_sp = 1, n_chem_species
            do k = 1, micm_levels
                offset = 1 + (k - 1) * cell_stride + (chem_species(i_sp)%micm_index - 1) * var_stride
                conc_val = state%concentrations(offset)

                do iCell = 1, target_cells
                    scalars(chem_species(i_sp)%mpas_index, k, iCell) = conc_val
                end do
            end do

            q_min = minval(scalars(chem_species(i_sp)%mpas_index, 1:micm_levels, 1:target_cells))
            q_max = maxval(scalars(chem_species(i_sp)%mpas_index, 1:micm_levels, 1:target_cells))
            call mpas_log_write('[MUSICA] Seeded ' // trim(chem_species(i_sp)%mpas_name) // &
                ' min=$r max=$r', realArgs=[q_min, q_max])
        end do

    end subroutine micm_to_mpas_chem

    !------------------------------------------------------------------------
    !  subroutine MICM_from_chemistry
    !
    !> \brief Copies MPAS chemistry state to MICM state with unit conversion
    !> \details
    !>  Loops over chem_species to convert MPAS mixing ratios (kg/kg)
    !>  to MICM concentrations (mol/m³) and sets environmental conditions.
    !------------------------------------------------------------------------
    subroutine MICM_from_chemistry(scalars, rho_air, temperature, pressure, &
            nCells, nVertLevels, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND
        use mpas_log, only: mpas_log_write

        real(kind=RKIND), dimension(:,:,:), intent(in) :: scalars      ! (tracer, level, cell)
        real(kind=RKIND), dimension(:,:), intent(in) :: rho_air        ! (level, cell) [kg/m³]
        real(kind=RKIND), dimension(:,:), intent(in) :: temperature    ! (level, cell) [K]
        real(kind=RKIND), dimension(:,:), intent(in) :: pressure       ! (level, cell) [Pa]
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: cell_stride, var_stride
        integer :: i_sp, k, iCell, idx, micm_cell
        real(real64) :: q_val, rho, conc

        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized in MICM_from_chemistry.'
            return
        end if

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do iCell = 1, nCells
            do k = 1, nVertLevels
                micm_cell = (iCell - 1) * nVertLevels + k
                rho = rho_air(k, iCell)

                ! Set MICM environmental conditions
                state%conditions(micm_cell)%temperature = temperature(k, iCell)
                state%conditions(micm_cell)%pressure = pressure(k, iCell)
                state%conditions(micm_cell)%air_density = rho / M_AIR

                ! Convert each species: mixing ratio [kg/kg] -> concentration [mol/m³]
                do i_sp = 1, n_chem_species
                    q_val = scalars(chem_species(i_sp)%mpas_index, k, iCell)
                    conc = q_val * rho / chem_species(i_sp)%molar_mass
                    idx = 1 + (micm_cell - 1) * cell_stride + (chem_species(i_sp)%micm_index - 1) * var_stride
                    state%concentrations(idx) = conc
                end do
            end do
        end do

    end subroutine MICM_from_chemistry

    !------------------------------------------------------------------------
    !  subroutine MICM_to_chemistry
    !
    !> \brief Copies MICM state to MPAS chemistry state with unit conversion
    !> \details
    !>  Loops over chem_species to convert MICM concentrations (mol/m³)
    !>  back to MPAS mixing ratios (kg/kg) after the solver step.
    !------------------------------------------------------------------------
    subroutine MICM_to_chemistry(scalars, rho_air, nCells, nVertLevels, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND
        use mpas_log, only: mpas_log_write

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars   ! (tracer, level, cell)
        real(kind=RKIND), dimension(:,:), intent(in) :: rho_air        ! (level, cell) [kg/m³]
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: cell_stride, var_stride
        integer :: i_sp, k, iCell, idx, micm_cell
        real(real64) :: rho, conc, q_val

        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized in MICM_to_chemistry.'
            return
        end if

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do iCell = 1, nCells
            do k = 1, nVertLevels
                micm_cell = (iCell - 1) * nVertLevels + k
                rho = rho_air(k, iCell)

                ! Convert each species: concentration [mol/m³] -> mixing ratio [kg/kg]
                do i_sp = 1, n_chem_species
                    idx = 1 + (micm_cell - 1) * cell_stride + (chem_species(i_sp)%micm_index - 1) * var_stride
                    conc = state%concentrations(idx)
                    q_val = conc * chem_species(i_sp)%molar_mass / rho
                    scalars(chem_species(i_sp)%mpas_index, k, iCell) = q_val
                end do
            end do
        end do

    end subroutine MICM_to_chemistry

    subroutine assign_rate_parameters(state, number_of_grid_cells, nox_tau)
        use iso_fortran_env, only : real64
        use mpas_kind_types, only : RKIND
        use mpas_log, only : mpas_log_write

        type(state_t), intent(inout) :: state
        integer, intent(in)          :: number_of_grid_cells
        real(kind=RKIND), intent(in) :: nox_tau     ! [s] NOx sink timescale (0 = disabled)

        integer :: cell_stride, var_stride, rp_index, cell, idx, n_rp
        character(len=256) :: msg
        character(len=:), allocatable :: rp_name
        real(real64) :: rate_val

        if (state%number_of_rate_parameters == 0) return

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride  = state%rate_parameters_strides%variable
        n_rp = state%rate_parameters_ordering%size()

        do rp_index = 1, n_rp
            rp_name = trim(state%rate_parameters_ordering%name(rp_index))

            ! Set rate based on parameter name.
            ! USER_DEFINED reactions: host multiplier — default 1.0 so the
            ! YAML scaling_factor defines the effective rate constant.
            ! PHOTO and LOSS parameters stay at 0 unless explicitly configured.
            if (rp_name(1:min(5,len(rp_name))) == 'USER.') then
                rate_val = 1.0_real64
            else
                rate_val = 0.0_real64
            end if
            if (rp_name == 'LOSS.no_loss' .or. rp_name == 'LOSS.no2_loss') then
                if (nox_tau > 0.0_RKIND) then
                    rate_val = 1.0_real64 / real(nox_tau, real64)
                end if
            end if

            do cell = 1, number_of_grid_cells
                idx = 1 + (cell - 1) * cell_stride + (rp_index - 1) * var_stride
                state%rate_parameters(idx) = rate_val
            end do

            write(msg, '(A,A,A,ES10.3)') '[MUSICA] Rate param "', rp_name, '" = ', rate_val
            call mpas_log_write(trim(msg))
        end do

    end subroutine assign_rate_parameters

    !------------------------------------------------------------------------
    !  subroutine musica_cache_photo_indices
    !
    !> \brief Resolve MICM rate-parameter indices for a list of rate names
    !> \details
    !>  For each name, looks up ``PHOTO.<name>`` in state%rate_parameters_ordering
    !>  and caches the stride index in photo_rp_indices(:). Errors if a rate
    !>  name is missing from the MICM ordering; the chemistry mechanism
    !>  (yaml) and the TUV-x config must be consistent about which
    !>  photolysis reactions exist.
    !------------------------------------------------------------------------
    subroutine musica_cache_photo_indices(names, error_code, error_message)
        use mpas_log, only : mpas_log_write

        character(len=*), intent(in)                :: names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        integer :: i, j, n_rp
        character(len=:), allocatable :: rp_name, target_name
        logical :: found

        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] musica_cache_photo_indices called before state init'
            return
        end if

        n_photo_rp = size(names)

        if (allocated(photo_rp_names))   deallocate(photo_rp_names)
        if (allocated(photo_rp_indices)) deallocate(photo_rp_indices)
        allocate(photo_rp_names(n_photo_rp))
        allocate(photo_rp_indices(n_photo_rp))

        if (state%number_of_rate_parameters == 0) then
            error_code = 1
            error_message = '[MUSICA] MICM state has zero rate parameters'
            return
        end if

        n_rp = state%rate_parameters_ordering%size()

        do i = 1, n_photo_rp
            photo_rp_names(i) = trim(names(i))
            target_name = 'PHOTO.' // trim(names(i))
            found = .false.
            do j = 1, n_rp
                rp_name = trim(state%rate_parameters_ordering%name(j))
                if (rp_name == target_name) then
                    photo_rp_indices(i) = j
                    found = .true.
                    call mpas_log_write('[MUSICA] Cached ' // target_name // &
                        ' rate param index = $i', intArgs=[j])
                    exit
                end if
            end do
            if (.not. found) then
                error_code = 1
                error_message = '[MUSICA] Rate parameter ' // target_name // &
                    ' not found; TUV-x and MICM configs disagree'
                return
            end if
        end do

    end subroutine musica_cache_photo_indices

    !------------------------------------------------------------------------
    !  subroutine musica_set_photolysis_rates
    !
    !> \brief Write photolysis rates for every cached reaction to MICM
    !> \details
    !>  Takes a (n_rates, nVertLevels, nCells) array of photolysis rates and
    !>  writes each rate's (k,iCell) value into both state and state_ref
    !>  using the cached stride indices. The reference state exists for
    !>  uncoupled advection diagnostics.
    !------------------------------------------------------------------------
    subroutine musica_set_photolysis_rates(rates, nCells, nVertLevels, n_rates, &
                                           error_code, error_message)
        use iso_fortran_env, only : real64
        use mpas_kind_types, only : RKIND
        use mpas_log, only : mpas_log_write

        real(kind=RKIND), intent(in) :: rates(:,:,:) ! (n_rates, nVertLevels, nCells)
        integer, intent(in)          :: nCells, nVertLevels, n_rates
        integer, intent(out)         :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: cell_stride, var_stride, iCell, k, micm_cell, idx, ir, rp_idx
        character(len=256) :: msg
        real(real64) :: rmin, rmax

        error_code = 0
        error_message = ''

        if (n_photo_rp < 1) return
        if (.not. associated(state)) return

        if (n_rates /= n_photo_rp) then
            error_code = 1
            error_message = '[MUSICA] musica_set_photolysis_rates: rate-count mismatch'
            return
        end if

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride  = state%rate_parameters_strides%variable

        do ir = 1, n_photo_rp
            rp_idx = photo_rp_indices(ir)
            do iCell = 1, nCells
                do k = 1, nVertLevels
                    micm_cell = (iCell - 1) * nVertLevels + k
                    idx = 1 + (micm_cell - 1) * cell_stride + (rp_idx - 1) * var_stride
                    state%rate_parameters(idx) = real(rates(ir, k, iCell), real64)
                end do
            end do

            if (associated(state_ref)) then
                do iCell = 1, nCells
                    do k = 1, nVertLevels
                        micm_cell = (iCell - 1) * nVertLevels + k
                        idx = 1 + (micm_cell - 1) * cell_stride + (rp_idx - 1) * var_stride
                        state_ref%rate_parameters(idx) = real(rates(ir, k, iCell), real64)
                    end do
                end do
            end if

            rmin = real(minval(rates(ir, 1:nVertLevels, 1:nCells)), real64)
            rmax = real(maxval(rates(ir, 1:nVertLevels, 1:nCells)), real64)
            write(msg, '(A,A,A,ES10.3,A,ES10.3)') '[MUSICA] ', &
                trim(photo_rp_names(ir)), ' min=', rmin, ' max=', rmax
            call mpas_log_write(trim(msg))
        end do

    end subroutine musica_set_photolysis_rates

    subroutine set_species_profile(state, species_name, values, err)
        ! Write a 1-D profile into the contiguous concentrations array.
        use iso_fortran_env, only : real64
        use musica_util, only : error_t, string_t

        type(state_t), intent(inout) :: state
        character(len=*), intent(in) :: species_name
        real(real64), intent(in)     :: values(:)
        type(error_t), intent(inout) :: err
        integer :: species_index, cell_stride, var_stride, cell, idx

        species_index = state%species_ordering%index(species_name, err)
        if (.not. err%is_success()) return

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do cell = 1, size(values)
            idx = 1 + (cell - 1) * cell_stride + (species_index - 1) * var_stride
            state%concentrations(idx) = values(cell)
        end do
    end subroutine set_species_profile

    subroutine log_concentrations(state, number_of_grid_cells, species_name, err)
        use iso_fortran_env, only : real64
        use mpas_log, only        : mpas_log_write
        use mpas_kind_types, only : RKIND
        use musica_util, only : error_t, string_t

        type(state_t), intent(in)    :: state
        integer, intent(in)          :: number_of_grid_cells
        character(len=*), intent(in) :: species_name
        type(error_t), intent(inout) :: err

        integer :: species_index, cell_stride, var_stride, cell, idx
        real (kind=RKIND) :: concentration

        species_index = state%species_ordering%index(species_name, err)
        if (.not. err%is_success()) return

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do cell = 1, number_of_grid_cells
            idx = 1 + (cell - 1) * cell_stride + (species_index - 1) * var_stride
            concentration = state%concentrations(idx)
            call mpas_log_write('  MICM concentration: $r', realArgs=[concentration])
        end do
    end subroutine log_concentrations

    subroutine log_solver_stats(solver_stats)
        use mpas_log, only        : mpas_log_write
        use mpas_kind_types, only : RKIND
        use musica_micm, only     : solver_stats_t

        type(solver_stats_t), intent(in) :: solver_stats
 
        integer :: function_calls, jacobian_updates, number_of_steps
        integer :: accepted, rejected, decompositions, solves

        real (kind=RKIND) :: final_time

        function_calls = solver_stats%function_calls()
        jacobian_updates = solver_stats%jacobian_updates()
        number_of_steps = solver_stats%number_of_steps()
        accepted = solver_stats%accepted()
        rejected = solver_stats%rejected()
        decompositions = solver_stats%decompositions()
        solves = solver_stats%solves()
        final_time = solver_stats%final_time()

        call mpas_log_write('[MUSICA] MICM Solver statistics ...')
        call mpas_log_write('  MICM function calls: $i', intArgs=[function_calls])
        call mpas_log_write('  MICM jacobian updates: $i', intArgs=[jacobian_updates])
        call mpas_log_write('  MICM number of steps: $i', intArgs=[number_of_steps])
        call mpas_log_write('  MICM accepted: $i', intArgs=[accepted])
        call mpas_log_write('  MICM rejected: $i', intArgs=[rejected])
        call mpas_log_write('  MICM LU decompositions: $i', intArgs=[decompositions])
        call mpas_log_write('  MICM linear solves: $i', intArgs=[solves])
        call mpas_log_write('  MICM final time (s): $r', realArgs=[final_time])
    end subroutine log_solver_stats

    !-------------------------------------------------------------------------
    !  function has_error_occurred
    !
    !> \author CheMPAS-A Developers
    !> \date   20 August 2025
    !  \details
    !>  Evaluate a MUSICA error for failure and convert to error data
    !>  @param[in] error The error code to evaluate and convert.
    !>  @param[out] error_message The error message.
    !>  @param[out] error_code The error code.
    !>  @return True for an error, false for success.
    !-------------------------------------------------------------------------
    logical function has_error_occurred(error, error_message, error_code)
        use musica_util, only: error_t

        type(error_t),                 intent(in)  :: error
        character(len=:), allocatable, intent(out) :: error_message
        integer,                       intent(out) :: error_code

        character(len=30) :: error_code_str

        if ( error%is_success( ) ) then
          error_code = 0
          error_message = ''
          has_error_occurred = .false.
          return
        end if
        error_code = error%code( )
        write(error_code_str, '(I30)') error%code( )
        error_message = '[MUSICA Error]: ' // error%category( ) // '[' // &
                        trim( adjustl( error_code_str ) ) // ']: ' // error%message( )
        has_error_occurred = .true.
    end function has_error_occurred

    !------------------------------------------------------------------------
    !  subroutine copy_state_to_ref
    !
    !> \brief Copies coupled state to reference state
    !> \details
    !>  Copies concentrations and conditions from the coupled state to the
    !>  reference state. Called after first chemistry_from_MPAS to ensure
    !>  both states start with identical initial conditions.
    !------------------------------------------------------------------------
    subroutine copy_state_to_ref()

        use mpas_log, only: mpas_log_write

        integer :: i

        if (.not. associated(state) .or. .not. associated(state_ref)) return

        ! Copy concentrations array
        state_ref%concentrations(:) = state%concentrations(:)

        ! Copy conditions (T, P, rho) for each grid cell
        do i = 1, size(state%conditions)
            state_ref%conditions(i)%temperature = state%conditions(i)%temperature
            state_ref%conditions(i)%pressure = state%conditions(i)%pressure
            state_ref%conditions(i)%air_density = state%conditions(i)%air_density
        end do

        call mpas_log_write('[MUSICA] Copied coupled state to reference state for advection diagnostics')

    end subroutine copy_state_to_ref

    !------------------------------------------------------------------------
    !  subroutine log_column_comparison
    !
    !> \brief Logs comparison of coupled vs uncoupled state for a probe column
    !> \details
    !>  Logs concentrations from both the coupled state (which has been advected
    !>  by MPAS) and the reference state (chemistry only) for a single column.
    !>  Differences indicate advection is occurring.
    !------------------------------------------------------------------------
    subroutine log_column_comparison(probe_cell, nVertLevels)

        use iso_fortran_env, only: real64
        use mpas_log, only: mpas_log_write
        use mpas_kind_types, only: RKIND

        integer, intent(in) :: probe_cell
        integer, intent(in) :: nVertLevels

        integer :: cell_stride, var_stride
        integer :: i_sp, k, micm_cell, idx
        real(real64) :: conc_coupled, conc_ref

        if (.not. associated(state) .or. .not. associated(state_ref)) return
        if (n_chem_species == 0) return

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        call mpas_log_write('[COMPARE] Probe cell=$i: Coupled vs Reference (chemistry-only)', &
            intArgs=[probe_cell])

        do k = 1, nVertLevels, max(1, nVertLevels/4)
            micm_cell = (probe_cell - 1) * nVertLevels + k

            do i_sp = 1, n_chem_species
                idx = 1 + (micm_cell - 1) * cell_stride + (chem_species(i_sp)%micm_index - 1) * var_stride
                conc_coupled = state%concentrations(idx)
                conc_ref = state_ref%concentrations(idx)

                call mpas_log_write('[COMPARE] lev=$i ' // trim(chem_species(i_sp)%micm_name) // &
                    ' coupled=$r ref=$r', &
                    intArgs=[k], &
                    realArgs=[real(conc_coupled, kind=RKIND), real(conc_ref, kind=RKIND)])
            end do
        end do

    end subroutine log_column_comparison

end module mpas_musica
