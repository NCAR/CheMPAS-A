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

    public :: musica_init, musica_step, musica_step_ref, musica_finalize, micm_to_mpas_abba
    public :: MICM_from_chemistry, MICM_to_chemistry
    public :: log_column_comparison, copy_state_to_ref

    type(micm_t),  pointer :: micm  => null ( )  ! Pointer to the MICM ODE solver instance
    type(state_t), pointer :: state => null ( )  ! Pointer to the state of the MICM solver (coupled to MPAS)
    type(state_t), pointer :: state_ref => null ( )  ! Reference state (uncoupled, chemistry only)

    ! Molar masses for unit conversion [kg/mol]
    ! TEMPORARY: ABBA-specific values; future: read from MICM species properties
    real(kind=8), parameter :: M_AB  = 0.058_8   ! kg/mol (placeholder for AB molecule)
    real(kind=8), parameter :: M_A   = 0.029_8   ! kg/mol (placeholder for A atom)
    real(kind=8), parameter :: M_B   = 0.029_8   ! kg/mol (placeholder for B atom)
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
            error_code, error_message, &
            registry_tracer_names)

        use iso_fortran_env, only: real64

        use musica_micm, only : RosenbrockStandardOrder
        use musica_util, only : error_t, string_t

        use mpas_log, only : mpas_log_write

        character(len=*),              intent(in)  :: filename_of_micm_configuration
        integer,                       intent(in)  :: number_of_grid_cells
        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        character(len=*), dimension(:), intent(in), optional :: registry_tracer_names

        ! TEMPORARY: ABBA specific initial concentrations
        real(real64), allocatable :: init_ab(:)
        real(real64), allocatable :: init_a(:)
        real(real64), allocatable :: init_b(:)

        character(len=*), parameter :: INITIAL_CONCENTRATION_PROPERTY = "__initial concentration"

        type(error_t)  :: error
        type(string_t) :: micm_version
        ! TEMPORARY: Hard-coded options for the MICM solver
        integer :: solver_type = RosenbrockStandardOrder
        integer :: i_species
        integer :: check_error_code
        character(len=:), allocatable :: check_error_message

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

        call check_registry_tracer_consistency(state, registry_tracer_names, &
                check_error_code, check_error_message)
        if (check_error_code /= 0) then
            error_code = check_error_code
            error_message = check_error_message
            return
        end if

        ! TEMPORARY: ABBA specific initialization
        call assign_rate_parameters(state, number_of_grid_cells)
        call assign_rate_parameters(state_ref, number_of_grid_cells)

        allocate(init_ab(number_of_grid_cells))
        allocate(init_a(number_of_grid_cells))
        allocate(init_b(number_of_grid_cells))

        init_ab = micm%get_species_property_double('AB', INITIAL_CONCENTRATION_PROPERTY, error)
        if (has_error_occurred(error, error_message, error_code)) return
        init_a = micm%get_species_property_double('A', INITIAL_CONCENTRATION_PROPERTY, error)
        if (has_error_occurred(error, error_message, error_code)) return
        init_b = micm%get_species_property_double('B', INITIAL_CONCENTRATION_PROPERTY, error)
        if (has_error_occurred(error, error_message, error_code)) return

        ! Initialize coupled state
        call set_species_profile(state, 'AB', init_ab, error)
        if (has_error_occurred(error, error_message, error_code)) return
        call set_species_profile(state, 'A', init_a, error)
        if (has_error_occurred(error, error_message, error_code)) return
        call set_species_profile(state, 'B', init_b, error)
        if (has_error_occurred(error, error_message, error_code)) return

        ! Initialize reference state (same initial conditions)
        call set_species_profile(state_ref, 'AB', init_ab, error)
        if (has_error_occurred(error, error_message, error_code)) return
        call set_species_profile(state_ref, 'A', init_a, error)
        if (has_error_occurred(error, error_message, error_code)) return
        call set_species_profile(state_ref, 'B', init_b, error)
        if (has_error_occurred(error, error_message, error_code)) return

    end subroutine musica_init

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
    !  subroutine micm_to_mpas_abba
    !
    !> \brief One-time copy of ABBA MICM initial state into MPAS tracers
    !------------------------------------------------------------------------
    subroutine micm_to_mpas_abba(scalars, nCells, nVertLevels, &
            index_qAB, index_qA, index_qB, error_code, error_message)

        use mpas_kind_types, only : RKIND
        use musica_util, only : error_t
        use mpas_log, only : mpas_log_write

        real (kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: index_qAB
        integer, intent(in) :: index_qA
        integer, intent(in) :: index_qB

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err

        integer :: spec_ab, spec_a, spec_b
        integer :: cell_stride, var_stride
        integer :: k, iCell
        integer :: micm_levels, target_cells, target_levels
        integer :: offset_ab, offset_a, offset_b
        real (kind=RKIND) :: conc_ab, conc_a, conc_b
        logical :: logged_first
        real (kind=RKIND) :: qAB_min, qAB_max

        error_code = 0
        error_message = ''
        logged_first = .false.

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized; cannot seed MPAS tracers.'
            call mpas_log_write(error_message)
            return
        end if

        call mpas_log_write('[MUSICA] Seeding MPAS tracers from MICM ABBA initial state.')

        spec_ab = state%species_ordering%index('AB', err)
        if (has_error_occurred(err, error_message, error_code)) return
        spec_a = state%species_ordering%index('A', err)
        if (has_error_occurred(err, error_message, error_code)) return
        spec_b = state%species_ordering%index('B', err)
        if (has_error_occurred(err, error_message, error_code)) return

        if (index_qAB < 1 .or. index_qA < 1 .or. index_qB < 1) then
            error_code = 1
            error_message = '[MUSICA] ABBA tracer indices are not valid; skipping MICM -> MPAS seed.'
            call mpas_log_write(error_message)
            return
        end if

        if (index_qAB > size(scalars, 1) .or. index_qA > size(scalars, 1) .or. &
                index_qB > size(scalars, 1)) then
            error_code = 1
            error_message = '[MUSICA] ABBA tracer indices exceed MPAS scalars dimension.'
            call mpas_log_write(error_message)
            return
        end if

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
        micm_levels   = min(state%number_of_grid_cells, target_levels)

        if (target_cells <= 0 .or. micm_levels <= 0) then
            error_code = 1
            error_message = '[MUSICA] No target cells/levels available for MICM -> MPAS seeding.'
            call mpas_log_write(error_message)
            return
        end if

        call mpas_log_write('[MUSICA] MICM->MPAS ABBA seed: micm_cells=$i target_levels=$i target_cells=$i', &
                intArgs=[state%number_of_grid_cells, micm_levels, target_cells])

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do k = 1, micm_levels
            offset_ab = 1 + (k - 1) * cell_stride + (spec_ab - 1) * var_stride
            offset_a  = 1 + (k - 1) * cell_stride + (spec_a  - 1) * var_stride
            offset_b  = 1 + (k - 1) * cell_stride + (spec_b  - 1) * var_stride

            conc_ab = state%concentrations(offset_ab)
            conc_a  = state%concentrations(offset_a)
            conc_b  = state%concentrations(offset_b)

            do iCell = 1, target_cells
                scalars(index_qAB, k, iCell) = conc_ab
                scalars(index_qA,  k, iCell) = conc_a
                scalars(index_qB,  k, iCell) = conc_b
                if (.not. logged_first) then
                    call mpas_log_write('[MUSICA] Seeded ABBA cell=1 level=1 -> qAB=$r qA=$r qB=$r', &
                            realArgs=[conc_ab, conc_a, conc_b])
                    logged_first = .true.
                end if
            end do
        end do

        qAB_min = minval(scalars(index_qAB, 1:target_levels, 1:target_cells))
        qAB_max = maxval(scalars(index_qAB, 1:target_levels, 1:target_cells))
        call mpas_log_write('[MUSICA] ABBA seed summary: qAB min=$r max=$r', realArgs=[qAB_min, qAB_max])

    end subroutine micm_to_mpas_abba

    !------------------------------------------------------------------------
    !  subroutine MICM_from_chemistry
    !
    !> \brief Copies MPAS chemistry state to MICM state with unit conversion
    !> \author CheMPAS-A Developers
    !> \date   January 2026
    !> \details
    !>  This subroutine copies MPAS scalar tracers (mixing ratios in kg/kg)
    !>  to MICM concentrations (mol/m³) and sets environmental conditions.
    !>  Unit conversion: C = q × rho_air / M_species
    !------------------------------------------------------------------------
    subroutine MICM_from_chemistry(scalars, rho_air, temperature, pressure, &
            nCells, nVertLevels, &
            index_qAB, index_qA, index_qB, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND
        use musica_util, only: error_t
        use mpas_log, only: mpas_log_write

        real(kind=RKIND), dimension(:,:,:), intent(in) :: scalars      ! (tracer, level, cell)
        real(kind=RKIND), dimension(:,:), intent(in) :: rho_air        ! (level, cell) [kg/m³]
        real(kind=RKIND), dimension(:,:), intent(in) :: temperature    ! (level, cell) [K]
        real(kind=RKIND), dimension(:,:), intent(in) :: pressure       ! (level, cell) [Pa]
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: index_qAB
        integer, intent(in) :: index_qA
        integer, intent(in) :: index_qB
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err
        integer :: spec_ab, spec_a, spec_b
        integer :: cell_stride, var_stride
        integer :: k, iCell, idx, micm_cell
        real(real64) :: q_ab, q_a, q_b
        real(real64) :: rho, T, P
        real(real64) :: conc_ab, conc_a, conc_b

        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized in MICM_from_chemistry.'
            return
        end if

        ! Get MICM species indices
        spec_ab = state%species_ordering%index('AB', err)
        if (has_error_occurred(err, error_message, error_code)) return
        spec_a = state%species_ordering%index('A', err)
        if (has_error_occurred(err, error_message, error_code)) return
        spec_b = state%species_ordering%index('B', err)
        if (has_error_occurred(err, error_message, error_code)) return

        ! Get strides for MICM array access
        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        ! Loop over cells and levels
        ! MICM state has nCells × nVertLevels grid cells
        ! MICM cell index = (iCell - 1) * nVertLevels + k
        do iCell = 1, nCells
            do k = 1, nVertLevels
                ! Compute MICM grid cell index (1-based)
                micm_cell = (iCell - 1) * nVertLevels + k

                ! Get MPAS values
                q_ab = scalars(index_qAB, k, iCell)
                q_a  = scalars(index_qA,  k, iCell)
                q_b  = scalars(index_qB,  k, iCell)
                rho  = rho_air(k, iCell)
                T    = temperature(k, iCell)
                P    = pressure(k, iCell)

                ! Set MICM environmental conditions
                state%conditions(micm_cell)%temperature = T
                state%conditions(micm_cell)%pressure = P
                state%conditions(micm_cell)%air_density = rho / M_AIR  ! Convert kg/m³ to mol/m³

                ! Convert mixing ratio [kg/kg] to concentration [mol/m³]
                ! C = q × rho_air / M_species
                conc_ab = q_ab * rho / M_AB
                conc_a  = q_a  * rho / M_A
                conc_b  = q_b  * rho / M_B

                ! Store in MICM state using strided access
                idx = 1 + (micm_cell - 1) * cell_stride + (spec_ab - 1) * var_stride
                state%concentrations(idx) = conc_ab
                idx = 1 + (micm_cell - 1) * cell_stride + (spec_a - 1) * var_stride
                state%concentrations(idx) = conc_a
                idx = 1 + (micm_cell - 1) * cell_stride + (spec_b - 1) * var_stride
                state%concentrations(idx) = conc_b
            end do
        end do

        ! Probe cell logging: log one cell for debugging (cell 100, level 1)
        if (nCells >= 100 .and. nVertLevels >= 1) then
            call mpas_log_write('[MICM_from] Probe cell=100 lev=1: qAB=$r qA=$r qB=$r T=$r P=$r rho=$r', &
                realArgs=[real(scalars(index_qAB, 1, 100), kind=RKIND), &
                          real(scalars(index_qA, 1, 100), kind=RKIND), &
                          real(scalars(index_qB, 1, 100), kind=RKIND), &
                          temperature(1, 100), pressure(1, 100), rho_air(1, 100)])
        end if

    end subroutine MICM_from_chemistry

    !------------------------------------------------------------------------
    !  subroutine MICM_to_chemistry
    !
    !> \brief Copies MICM state to MPAS chemistry state with unit conversion
    !> \author CheMPAS-A Developers
    !> \date   January 2026
    !> \details
    !>  This subroutine copies MICM concentrations (mol/m³) back to
    !>  MPAS scalar tracers (mixing ratios in kg/kg) after the solver step.
    !>  Unit conversion: q = C × M_species / rho_air
    !------------------------------------------------------------------------
    subroutine MICM_to_chemistry(scalars, rho_air, nCells, nVertLevels, &
            index_qAB, index_qA, index_qB, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND
        use musica_util, only: error_t
        use mpas_log, only: mpas_log_write

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars   ! (tracer, level, cell)
        real(kind=RKIND), dimension(:,:), intent(in) :: rho_air        ! (level, cell) [kg/m³]
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: index_qAB
        integer, intent(in) :: index_qA
        integer, intent(in) :: index_qB
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err
        integer :: spec_ab, spec_a, spec_b
        integer :: cell_stride, var_stride
        integer :: k, iCell, idx, micm_cell
        real(real64) :: rho
        real(real64) :: conc_ab, conc_a, conc_b
        real(real64) :: q_ab, q_a, q_b

        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized in MICM_to_chemistry.'
            return
        end if

        ! Get MICM species indices
        spec_ab = state%species_ordering%index('AB', err)
        if (has_error_occurred(err, error_message, error_code)) return
        spec_a = state%species_ordering%index('A', err)
        if (has_error_occurred(err, error_message, error_code)) return
        spec_b = state%species_ordering%index('B', err)
        if (has_error_occurred(err, error_message, error_code)) return

        ! Get strides for MICM array access
        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        ! Loop over cells and levels
        ! MICM state has nCells × nVertLevels grid cells
        ! MICM cell index = (iCell - 1) * nVertLevels + k
        do iCell = 1, nCells
            do k = 1, nVertLevels
                ! Compute MICM grid cell index (1-based)
                micm_cell = (iCell - 1) * nVertLevels + k

                rho = rho_air(k, iCell)

                ! Read from MICM state using strided access
                idx = 1 + (micm_cell - 1) * cell_stride + (spec_ab - 1) * var_stride
                conc_ab = state%concentrations(idx)
                idx = 1 + (micm_cell - 1) * cell_stride + (spec_a - 1) * var_stride
                conc_a = state%concentrations(idx)
                idx = 1 + (micm_cell - 1) * cell_stride + (spec_b - 1) * var_stride
                conc_b = state%concentrations(idx)

                ! Convert concentration [mol/m³] to mixing ratio [kg/kg]
                ! q = C × M_species / rho_air
                q_ab = conc_ab * M_AB / rho
                q_a  = conc_a  * M_A  / rho
                q_b  = conc_b  * M_B  / rho

                ! Store in MPAS scalars
                scalars(index_qAB, k, iCell) = q_ab
                scalars(index_qA,  k, iCell) = q_a
                scalars(index_qB,  k, iCell) = q_b
            end do
        end do

        ! Probe cell logging: log one cell for debugging (cell 100, level 1)
        if (nCells >= 100 .and. nVertLevels >= 1) then
            call mpas_log_write('[MICM_to] Probe cell=100 lev=1: qAB=$r qA=$r qB=$r', &
                realArgs=[real(scalars(index_qAB, 1, 100), kind=RKIND), &
                          real(scalars(index_qA, 1, 100), kind=RKIND), &
                          real(scalars(index_qB, 1, 100), kind=RKIND)])
        end if

    end subroutine MICM_to_chemistry

    subroutine assign_rate_parameters(state, number_of_grid_cells)
        ! Provide a default value for any rate parameters (none for ABBA by default).
        use iso_fortran_env, only : real64

        type(state_t), intent(inout) :: state
        integer, intent(in)          :: number_of_grid_cells

        integer :: cell_stride, var_stride, rp_index, cell, idx

        if (state%number_of_rate_parameters == 0) return

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride  = state%rate_parameters_strides%variable

        do rp_index = 1, state%rate_parameters_ordering%size()
            do cell = 1, number_of_grid_cells
                idx = 1 + (cell - 1) * cell_stride + (rp_index - 1) * var_stride
                state%rate_parameters(idx) = 1.0_real64
            end do
        end do
    end subroutine assign_rate_parameters

    subroutine check_registry_tracer_consistency(state, registry_tracer_names, &
            error_code, error_message)
        use mpas_log, only : mpas_log_write

        type(state_t), intent(in) :: state
        character(len=*), dimension(:), intent(in), optional :: registry_tracer_names
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: i_species, i_tracer, registry_count
        logical :: found
        logical, allocatable :: registry_used(:)
        character(len=:), allocatable :: missing_list
        character(len=:), allocatable :: extra_list
        character(len=:), allocatable :: species_name
        character(len=:), allocatable :: expected_tracer

        error_code = 0
        error_message = ''

        if (.not. present(registry_tracer_names)) then
            call mpas_log_write('[MUSICA] Registry tracers not provided; skipping MICM/MPAS consistency check.')
            return
        end if

        registry_count = size(registry_tracer_names)
        if (registry_count == 0) then
            call mpas_log_write('[MUSICA] Registry tracer list is empty; skipping MICM/MPAS consistency check.')
            return
        end if

        allocate(registry_used(registry_count))
        registry_used = .false.
        missing_list = ''
        extra_list = ''

        do i_species = 1, state%species_ordering%size()
            species_name = trim(state%species_ordering%name(i_species))
            expected_tracer = 'q' // species_name
            found = .false.
            do i_tracer = 1, registry_count
                if (expected_tracer == trim(registry_tracer_names(i_tracer))) then
                    found = .true.
                    registry_used(i_tracer) = .true.
                    call mpas_log_write('[MUSICA] Matched MICM species ' // trim(species_name) // ' to MPAS tracer ' // trim(expected_tracer))
                    exit
                end if
            end do
            if (.not. found) then
                missing_list = missing_list // ' ' // expected_tracer
            end if
        end do

        do i_tracer = 1, registry_count
            if (.not. registry_used(i_tracer)) then
                extra_list = extra_list // ' ' // trim(registry_tracer_names(i_tracer))
            end if
        end do

        if (len_trim(missing_list) > 0) then
            error_code = 1
            error_message = '[MUSICA] MICM species missing in MPAS registry tracers:' // trim(adjustl(missing_list))
            call mpas_log_write(error_message)
        else
            call mpas_log_write('[MUSICA] MICM/MPAS tracer consistency check passed.')
        end if

        if (len_trim(extra_list) > 0) then
            call mpas_log_write('[MUSICA] MPAS registry tracers not present in MICM config:' // trim(adjustl(extra_list)))
        end if

        deallocate(registry_used)
    end subroutine check_registry_tracer_consistency

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
        use musica_util, only: error_t
        use mpas_log, only: mpas_log_write
        use mpas_kind_types, only: RKIND

        integer, intent(in) :: probe_cell   ! The cell index to log
        integer, intent(in) :: nVertLevels  ! Number of vertical levels

        type(error_t) :: err
        integer :: spec_ab, spec_a, spec_b
        integer :: cell_stride, var_stride
        integer :: k, micm_cell, idx
        real(real64) :: conc_ab, conc_a, conc_b
        real(real64) :: conc_ab_ref, conc_a_ref, conc_b_ref

        if (.not. associated(state) .or. .not. associated(state_ref)) return

        ! Get MICM species indices
        spec_ab = state%species_ordering%index('AB', err)
        spec_a = state%species_ordering%index('A', err)
        spec_b = state%species_ordering%index('B', err)

        ! Get strides
        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        call mpas_log_write('[COMPARE] Probe cell=$i: Coupled vs Reference (chemistry-only)', &
            intArgs=[probe_cell])
        call mpas_log_write('[COMPARE] Level | AB_coupled AB_ref | A_coupled A_ref | B_coupled B_ref')

        ! Log a few levels (1, middle, top)
        do k = 1, nVertLevels, max(1, nVertLevels/4)
            micm_cell = (probe_cell - 1) * nVertLevels + k

            ! Coupled state
            idx = 1 + (micm_cell - 1) * cell_stride + (spec_ab - 1) * var_stride
            conc_ab = state%concentrations(idx)
            idx = 1 + (micm_cell - 1) * cell_stride + (spec_a - 1) * var_stride
            conc_a = state%concentrations(idx)
            idx = 1 + (micm_cell - 1) * cell_stride + (spec_b - 1) * var_stride
            conc_b = state%concentrations(idx)

            ! Reference state
            idx = 1 + (micm_cell - 1) * cell_stride + (spec_ab - 1) * var_stride
            conc_ab_ref = state_ref%concentrations(idx)
            idx = 1 + (micm_cell - 1) * cell_stride + (spec_a - 1) * var_stride
            conc_a_ref = state_ref%concentrations(idx)
            idx = 1 + (micm_cell - 1) * cell_stride + (spec_b - 1) * var_stride
            conc_b_ref = state_ref%concentrations(idx)

            call mpas_log_write('[COMPARE] $i | $r $r | $r $r | $r $r', &
                intArgs=[k], &
                realArgs=[real(conc_ab, kind=RKIND), real(conc_ab_ref, kind=RKIND), &
                          real(conc_a, kind=RKIND), real(conc_a_ref, kind=RKIND), &
                          real(conc_b, kind=RKIND), real(conc_b_ref, kind=RKIND)])
        end do

    end subroutine log_column_comparison

end module mpas_musica
