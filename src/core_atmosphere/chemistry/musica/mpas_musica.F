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

    use iso_fortran_env, only: real64
    use musica_micm,  only: micm_t, get_micm_version
    ! get_micm_version is here to avoid an ICE when within musica_init
    use musica_state, only: state_t
    use mpas_chem_constants, only: CHEM_GAS_R, CHEM_M_AIR, CHEM_NA, &
                                   CHEM_O2_VMR, CHEM_N2_VMR

    implicit none

    private

    public :: musica_init, musica_step, musica_step_ref, musica_finalize, micm_to_mpas_chem
    public :: MICM_from_chemistry, MICM_to_chemistry
    public :: log_column_comparison, copy_state_to_ref
    public :: resolve_mpas_indices, n_chem_species, chem_species
    public :: musica_query_species, musica_query_emission_species
    public :: musica_get_solver_tolerances
    public :: musica_tracer_name
    public :: musica_cache_photo_indices, musica_set_photolysis_rates
    public :: musica_cache_emission_indices, musica_set_emission_fluxes
    public :: musica_mozart35_rate
    public :: n_photo_rp, photo_rp_names

    type(micm_t),  pointer :: micm  => null ( )  ! Pointer to the MICM ODE solver instance
    type(state_t), pointer :: state => null ( )  ! Pointer to the state of the MICM solver (coupled to MPAS)
    type(state_t), pointer :: state_ref => null ( )  ! Reference state (uncoupled, chemistry only)

    integer, parameter :: MAX_SUB_CALLS = 100
    real(real64), parameter :: TIME_EPS_REL = 1.0e-9_real64
    ! Photolysis rate-parameter bookkeeping. For each TUV-x reaction name we
    ! resolve the MICM rate_parameters_ordering key PHOTO.<name> once at
    ! init, then reuse the cached stride index on every step.
    integer, save                    :: n_photo_rp = 0
    integer,             allocatable :: photo_rp_indices(:)
    character(len=64),   allocatable :: photo_rp_names(:)

    type :: emission_entry_t
        character(len=64) :: name
        integer           :: species_index
        integer           :: rate_parameter_index
        real(real64)      :: molar_mass
        logical           :: allows_net_flux = .false.
    end type emission_entry_t

    integer, save :: n_emission_rp = 0
    type(emission_entry_t), allocatable, save :: emission_entries(:)

    type :: species_entry_t
        character(len=64) :: micm_name    ! e.g. 'AB'
        character(len=64) :: mpas_name    ! e.g. 'qAB' (or a host field name like 'qv')
        integer           :: micm_index   ! position in MICM species_ordering
        integer           :: mpas_index   ! position in MPAS scalars (1st dim)
        real(kind=8)      :: molar_mass   ! kg/mol from MICM config
        ! Host binding: when .true., this MICM species reads from an existing
        ! MPAS prognostic field (mpas_name above, e.g. qv for H2O) rather than a
        ! freshly injected qXXX chemistry tracer. write_back controls whether the
        ! MICM result is copied back into that field after the solve.
        logical           :: host_bound = .false.
        logical           :: write_back  = .true.
    end type

    integer :: n_chem_species = 0
    type(species_entry_t), allocatable :: chem_species(:)

    ! Negative-value clip accounting. Accumulated rank-locally every chemistry
    ! step (no collective), then reduced and logged once at finalize. Doing the
    ! reduction per step would deadlock: chemistry_step has rank-local
    ! recoverable early returns on solver failure, and a collective downstream
    ! of one rank's early return hangs the run.
    integer,      save :: clip_count_in  = 0
    integer,      save :: clip_count_out = 0
    real(real64), save :: clip_min_in    = huge(1.0_real64)
    real(real64), save :: clip_min_out   = huge(1.0_real64)

    interface musica_set_emission_fluxes
        module procedure musica_set_surface_emission_fluxes
        module procedure musica_set_layer_emission_fluxes
    end interface musica_set_emission_fluxes

    contains

    pure function musica_tracer_name(species_name) result(name)

        character(len=*), intent(in) :: species_name
        character(len=64) :: name

        name = 'q' // trim(species_name)

    end function musica_tracer_name

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
            relative_tolerance, &
            enable_ref, &
            error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND

        use musica_micm, only : BackwardEulerStandardOrder, RosenbrockStandardOrder, &
                                rosenbrock_solver_parameters_t
        use musica_util, only : error_t, string_t

        use mpas_log, only : mpas_log_write

        character(len=*),              intent(in)  :: filename_of_micm_configuration
        integer,                       intent(in)  :: number_of_grid_cells
        real(kind=RKIND),              intent(in)  :: nox_tau           ! [s] NOx sink timescale (0 = disabled)
        real(real64),                  intent(in)  :: relative_tolerance ! MICM adaptive solver tolerance
        logical,                       intent(in)  :: enable_ref
        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=*), parameter :: INITIAL_CONCENTRATION_PROPERTY = "__initial concentration"
        character(len=*), parameter :: MOLAR_MASS_PROPERTY = "__molar mass"
        character(len=*), parameter :: ABSOLUTE_TOLERANCE_PROPERTY = "__absolute tolerance"
        real(real64), parameter :: MICM_DEFAULT_ABSOLUTE_TOLERANCE = 1.0e-3_real64

        real(real64), allocatable :: init_conc(:)
        real(real64) :: molar_mass_val, absolute_tolerance

        type(error_t)  :: error
        type(string_t) :: micm_version
        type(rosenbrock_solver_parameters_t) :: ros_params
        type(state_t), pointer :: tolerance_state
        integer :: solver_type = RosenbrockStandardOrder
        integer :: i_species, spec_idx
        character(len=:), allocatable :: sp_name

        call reset_emission_cache()
        micm_version = get_micm_version()

        call mpas_log_write('[MUSICA] Initializing MICM chemistry package...')
        call mpas_log_write('MICM version: ' // micm_version%value_)
        call mpas_log_write('MICM number of grid cells: $i', intArgs=[number_of_grid_cells])

        micm => micm_t(trim(filename_of_micm_configuration), solver_type, error)
        if (has_error_occurred(error, error_message, error_code)) return

        ! mechanism_configuration v1 intentionally retains application-specific
        ! species keys with their leading "__". MICM's builder, however, only
        ! consumes its native "absolute tolerance" property and otherwise uses
        ! 1e-3. Translate the CheMPAS configuration property explicitly before
        ! allocating any production state so the configured science tolerances
        ! are the effective solver tolerances rather than documentation only.
        tolerance_state => micm%get_state(1, error)
        if (has_error_occurred(error, error_message, error_code)) return

        ros_params = micm%get_rosenbrock_solver_parameters(error)
        if (has_error_occurred(error, error_message, error_code)) then
            deallocate(tolerance_state)
            return
        end if
        ros_params%relative_tolerance = real(relative_tolerance, real64)
        if (allocated(ros_params%absolute_tolerances)) then
            deallocate(ros_params%absolute_tolerances)
        end if
        allocate(ros_params%absolute_tolerances(tolerance_state%species_ordering%size()))
        ros_params%absolute_tolerances = MICM_DEFAULT_ABSOLUTE_TOLERANCE

        do i_species = 1, tolerance_state%species_ordering%size()
            sp_name = trim(tolerance_state%species_ordering%name(i_species))
            spec_idx = tolerance_state%species_ordering%index(sp_name, error)
            if (has_error_occurred(error, error_message, error_code)) then
                deallocate(tolerance_state)
                return
            end if
            absolute_tolerance = micm%get_species_property_double( &
                sp_name, ABSOLUTE_TOLERANCE_PROPERTY, error)
            if (error%is_success()) then
                if (.not. ieee_is_finite(absolute_tolerance) .or. absolute_tolerance <= 0.0_real64) then
                    error_code = 1
                    error_message = '[MUSICA] Invalid absolute tolerance for species ' // &
                        trim(sp_name) // '; expected a finite value greater than zero.'
                    deallocate(tolerance_state)
                    return
                end if
                ros_params%absolute_tolerances(spec_idx) = absolute_tolerance
            else
                call mpas_log_write('[MUSICA] Species ' // trim(sp_name) // &
                    ' has no __absolute tolerance; using MICM default 1e-3')
            end if
            call mpas_log_write('[MUSICA] Effective absolute tolerance for ' // &
                trim(sp_name) // ': $r', &
                realArgs=[real(ros_params%absolute_tolerances(spec_idx), RKIND)])
        end do

        deallocate(tolerance_state)
        call micm%set_rosenbrock_solver_parameters(ros_params, error)
        if (has_error_occurred(error, error_message, error_code)) return
        call mpas_log_write('[MUSICA] MICM relative tolerance set to $r', &
            realArgs=[real(relative_tolerance, RKIND)])

        state => micm%get_state(number_of_grid_cells, error)
        if (has_error_occurred(error, error_message, error_code)) return

        nullify(state_ref)
        if (enable_ref) then
            state_ref => micm%get_state(number_of_grid_cells, error)
            if (has_error_occurred(error, error_message, error_code)) return
            call mpas_log_write('[MUSICA] Created reference state for advection diagnostics')
        end if

        do i_species = 1, state%species_ordering%size()
            call mpas_log_write('[MUSICA] MICM species: ' // state%species_ordering%name(i_species))
        end do

        call assign_rate_parameters(state, number_of_grid_cells, nox_tau)
        if (associated(state_ref)) then
            call assign_rate_parameters(state_ref, number_of_grid_cells, nox_tau)
        end if

        ! Build species table from MICM config
        n_chem_species = state%species_ordering%size()
        if (allocated(chem_species)) deallocate(chem_species)
        allocate(chem_species(n_chem_species))
        allocate(init_conc(number_of_grid_cells))

        do i_species = 1, n_chem_species
            sp_name = trim(state%species_ordering%name(i_species))
            chem_species(i_species)%micm_name = sp_name
            ! Host-bound species (e.g. H2O -> qv) read from an existing MPAS
            ! field instead of an injected qXXX tracer. host_binding sets the
            ! mpas_name to the host field and clears write_back; all other
            ! species keep the 'q'-prefixed name and are written back normally.
            block
                logical :: is_bound, do_wb
                character(len=64) :: host_name
                call host_binding(micm, sp_name, is_bound, host_name, do_wb)
                chem_species(i_species)%host_bound = is_bound
                chem_species(i_species)%write_back = do_wb
                if (is_bound) then
                    chem_species(i_species)%mpas_name = trim(host_name)
                    call mpas_log_write('[MUSICA] Species ' // trim(sp_name) // &
                        ' bound to host field ' // trim(host_name) // &
                        ' (read-only)')
                else
                    chem_species(i_species)%mpas_name = musica_tracer_name(sp_name)
                end if
            end block
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
            if (associated(state_ref)) then
                call set_species_profile(state_ref, sp_name, init_conc, error)
                if (has_error_occurred(error, error_message, error_code)) return
            end if
        end do

        ! Photolysis rate-parameter indices are cached later, after TUV-x
        ! init, by chemistry_init calling musica_cache_photo_indices with
        ! the reaction names TUV-x enumerated.

        deallocate(init_conc)

    end subroutine musica_init

    subroutine musica_get_solver_tolerances(relative_tolerance, absolute_tolerances, &
            species_names, error_code, error_message)

        use musica_micm, only: rosenbrock_solver_parameters_t
        use musica_util, only: error_t

        real(real64), intent(out) :: relative_tolerance
        real(real64), allocatable, intent(out) :: absolute_tolerances(:)
        character(len=64), allocatable, intent(out) :: species_names(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: error
        type(rosenbrock_solver_parameters_t) :: parameters
        integer :: i_species, n_species, spec_idx
        character(len=:), allocatable :: sp_name

        error_code = 0
        error_message = ''
        relative_tolerance = 0.0_real64
        if (.not. associated(micm) .or. .not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] Solver tolerances requested before musica_init.'
            return
        end if

        parameters = micm%get_rosenbrock_solver_parameters(error)
        if (has_error_occurred(error, error_message, error_code)) return
        n_species = state%species_ordering%size()
        if (.not. allocated(parameters%absolute_tolerances) .or. &
            size(parameters%absolute_tolerances) /= n_species) then
            error_code = 1
            error_message = '[MUSICA] Effective absolute-tolerance vector is missing or has the wrong size.'
            return
        end if

        relative_tolerance = parameters%relative_tolerance
        allocate(absolute_tolerances(n_species))
        allocate(species_names(n_species))
        absolute_tolerances = parameters%absolute_tolerances
        do i_species = 1, n_species
            sp_name = trim(state%species_ordering%name(i_species))
            spec_idx = state%species_ordering%index(sp_name, error)
            if (has_error_occurred(error, error_message, error_code)) return
            species_names(spec_idx) = sp_name
        end do

    end subroutine musica_get_solver_tolerances

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
    !>  Creates a temporary MICM instance to read species_ordering and molar
    !>  masses from the config file. Returns metadata for every MICM species,
    !>  including species that bind to existing MPAS host fields. The temporary
    !>  MICM instance is destroyed. Used by atm_setup_block to extend scalars
    !>  and define output-only chemistry diagnostics before chemistry_init.
    !------------------------------------------------------------------------
    subroutine musica_query_species(config_path, num_species, micm_names, mpas_names, &
            molar_masses, host_bound, fixed_vmrs, error_code, error_message)

        use iso_fortran_env, only: real64
        use musica_micm,  only: RosenbrockStandardOrder
        use musica_state, only: state_t
        use musica_util,  only: error_t
        use mpas_log,     only: mpas_log_write

        character(len=*), intent(in)                :: config_path
        integer, intent(out)                        :: num_species
        character(len=64), allocatable, intent(out)  :: micm_names(:)
        character(len=64), allocatable, intent(out)  :: mpas_names(:)
        real(real64), allocatable, intent(out)       :: molar_masses(:)
        logical, allocatable, intent(out)            :: host_bound(:)
        real(real64), allocatable, intent(out)       :: fixed_vmrs(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        character(len=*), parameter :: MOLAR_MASS_PROPERTY = "__molar mass"
        type(micm_t), pointer  :: tmp_micm => null()
        type(state_t), pointer :: tmp_state => null()
        type(error_t)          :: error
        integer                :: i
        logical                :: is_bound, do_wb
        real(real64)           :: fixed_vmr
        character(len=64)      :: host_name, micm_nm

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
        allocate(micm_names(num_species))
        allocate(mpas_names(num_species))
        allocate(molar_masses(num_species))
        allocate(host_bound(num_species))
        allocate(fixed_vmrs(num_species))

        do i = 1, num_species
            micm_nm = trim(tmp_state%species_ordering%name(i))
            micm_names(i) = micm_nm
            call host_binding(tmp_micm, micm_nm, is_bound, host_name, do_wb, fixed_vmr)
            host_bound(i) = is_bound
            fixed_vmrs(i) = fixed_vmr
            if (is_bound) then
                mpas_names(i) = trim(host_name)
                call mpas_log_write('[MUSICA] Species ' // trim(micm_nm) // &
                    ' bound to host field ' // trim(host_name) // &
                    '; no chemistry tracer injected')
            else
                mpas_names(i) = musica_tracer_name(micm_nm)
                call mpas_log_write('[MUSICA] Discovered chemistry tracer: ' // &
                    trim(mpas_names(i)))
            end if

            molar_masses(i) = tmp_micm%get_species_property_double( &
                micm_nm, MOLAR_MASS_PROPERTY, error)
            if (has_error_occurred(error, error_message, error_code)) goto 100
        end do

100     continue
        if (associated(tmp_state)) deallocate(tmp_state)
        if (associated(tmp_micm)) deallocate(tmp_micm)

    end subroutine musica_query_species

    !------------------------------------------------------------------------
    !  subroutine musica_query_emission_species
    !
    !> \brief Query EMIS.* names from a MICM chemistry configuration
    !> \details
    !>  Creates a temporary one-cell MICM state, selects rate parameters with
    !>  the ``EMIS.`` prefix, and returns their stripped species names. This
    !>  routine must receive config_micm_file; a MIEM emissions-only file is
    !>  not a source of MICM diagnostic names.
    !------------------------------------------------------------------------
    subroutine musica_query_emission_species(config_path, num_species, species_names, &
                                             error_code, error_message)

        use musica_micm, only: RosenbrockStandardOrder
        use musica_state, only: state_t
        use musica_util, only: error_t
        use mpas_log, only: mpas_log_write

        character(len=*), intent(in) :: config_path
        integer, intent(out) :: num_species
        character(len=64), allocatable, intent(out) :: species_names(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(micm_t), pointer :: tmp_micm => null()
        type(state_t), pointer :: tmp_state => null()
        type(error_t) :: error
        character(len=:), allocatable :: rate_name
        integer :: i, j, n_found, storage_index

        error_code = 0
        error_message = ''
        num_species = 0

        if (len_trim(config_path) == 0) then
            error_code = 1
            error_message = '[MUSICA] MICM config path is empty during emission diagnostic discovery'
            allocate(species_names(0))
            return
        end if

        call mpas_log_write('[MUSICA] Querying MICM config for EMIS.* diagnostics...')

        tmp_micm => micm_t(trim(config_path), RosenbrockStandardOrder, error)
        if (has_error_occurred(error, error_message, error_code)) goto 100
        if (.not. associated(tmp_micm)) then
            error_code = 1
            error_message = '[MUSICA] MICM diagnostic query returned no mechanism'
            goto 100
        end if

        tmp_state => tmp_micm%get_state(1, error)
        if (has_error_occurred(error, error_message, error_code)) goto 100
        if (.not. associated(tmp_state)) then
            error_code = 1
            error_message = '[MUSICA] MICM diagnostic query returned no state'
            goto 100
        end if

        n_found = 0
        do i = 1, tmp_state%rate_parameters_ordering%size()
            rate_name = trim(tmp_state%rate_parameters_ordering%name(i))
            if (index(rate_name, 'EMIS.') == 1) n_found = n_found + 1
        end do

        allocate(species_names(n_found))
        species_names = ''
        n_found = 0
        do i = 1, tmp_state%rate_parameters_ordering%size()
            rate_name = trim(tmp_state%rate_parameters_ordering%name(i))
            if (index(rate_name, 'EMIS.') /= 1) cycle

            if (len_trim(rate_name) <= 5) then
                error_code = 1
                error_message = '[MUSICA] MICM contains an empty EMIS. diagnostic suffix'
                goto 100
            end if
            if (len_trim(rate_name) - 5 > len(species_names)) then
                error_code = 1
                error_message = '[MUSICA] MICM emission diagnostic species name exceeds 64 characters'
                goto 100
            end if

            storage_index = tmp_state%rate_parameters_ordering%index(rate_name, error)
            if (.not. error%is_success() .or. storage_index < 1) then
                error_code = 1
                error_message = '[MUSICA] Could not resolve emission diagnostic parameter ' // rate_name
                goto 100
            end if

            n_found = n_found + 1
            species_names(n_found) = trim(rate_name(6:))
            do j = 1, n_found - 1
                if (trim(species_names(j)) == trim(species_names(n_found))) then
                    error_code = 1
                    error_message = '[MUSICA] Duplicate MICM emission diagnostic species: ' // &
                                    trim(species_names(n_found))
                    goto 100
                end if
            end do
        end do
        num_species = n_found

        do i = 1, num_species
            call mpas_log_write('[MUSICA] Discovered emissions diagnostic: emis_' // &
                                trim(species_names(i)))
        end do

100     continue
        if (error_code /= 0) then
            num_species = 0
            if (allocated(species_names)) deallocate(species_names)
            allocate(species_names(0))
        end if
        if (associated(tmp_state)) deallocate(tmp_state)
        if (associated(tmp_micm)) deallocate(tmp_micm)

    end subroutine musica_query_emission_species

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

        use mpas_log, only : mpas_log_write

        real(real64), intent(in) :: time_step

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        ! Re-entrant solve loop: MICM's adaptive controller returns early when
        ! it exhausts its internal step budget (max_number_of_steps, default 11)
        ! before reaching the requested interval. solver_stats%final_time()
        ! reports how far it actually advanced. Resubmit the remainder until
        ! the full time_step is covered or the safety cap is hit. This is the
        ! pattern used by MUSICA/fortran/examples/column_model.F90 and is
        ! required for stiff transients (e.g. Chapman chemistry at sunrise).

        call mpas_log_write('[MUSICA] Stepping MICM solver...')

        call micm_substep_loop(state, time_step, error_code, error_message)

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

        real(real64), intent(in) :: time_step

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        if (.not. associated(state_ref)) then
            error_code = 1
            error_message = '[MUSICA] musica_step_ref: reference state is not enabled'
            return
        end if

        call micm_substep_loop(state_ref, time_step, error_code, error_message)

    end subroutine musica_step_ref

    subroutine micm_substep_loop(a_state, time_step, error_code, error_message)

        use musica_micm, only : solver_stats_t
        use musica_util, only : error_t, string_t
        use mpas_log, only : mpas_log_write
        use mpas_kind_types, only : RKIND

        type(state_t), pointer, intent(inout) :: a_state
        real(real64), intent(in)              :: time_step

        integer,                       intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(string_t)       :: solver_state
        type(solver_stats_t) :: solver_stats
        type(error_t)        :: error

        real(real64) :: elapsed, remaining, advanced
        integer      :: sub_call
        logical      :: is_ref_state

        error_code = 0
        error_message = ''

        if (.not. associated(a_state)) then
            error_code = 1
            error_message = '[MUSICA] MICM state not initialized in sub-step loop.'
            return
        end if
        if (a_state%number_of_grid_cells < 1) return

        elapsed = 0.0_real64
        sub_call = 0
        is_ref_state = associated(a_state, state_ref)

        do while (time_step - elapsed > TIME_EPS_REL * time_step .and. sub_call < MAX_SUB_CALLS)
            sub_call = sub_call + 1
            remaining = time_step - elapsed

            call micm%solve(remaining, a_state, solver_state, solver_stats, error)
            if (has_error_occurred(error, error_message, error_code)) return

            advanced = solver_stats%final_time()

            if (advanced <= 0.0_real64) then
                if (remaining <= TIME_EPS_REL * time_step) exit

                if (is_ref_state) then
                    error_code = 1
                    error_message = '[MUSICA] Reference MICM stalled at sub-step'
                else
                    call mpas_log_write('[MUSICA Error] MICM advanced 0 s with state "' // &
                        solver_state%get_char_array() // '"; aborting sub-step loop')
                    error_code = 1
                    error_message = '[MUSICA] MICM stalled at sub-step'
                end if
                return
            end if

            elapsed = elapsed + advanced
        end do

        if (time_step - elapsed > TIME_EPS_REL * time_step) then
            if (is_ref_state) then
                error_code = 1
                error_message = '[MUSICA] Reference MICM sub-step budget exhausted'
            else
                call mpas_log_write('[MUSICA Error] MICM exhausted sub-call budget ($i calls); ' // &
                    'elapsed $r of $r s', &
                    intArgs=[MAX_SUB_CALLS], &
                    realArgs=[real(elapsed, RKIND), real(time_step, RKIND)])
                error_code = 1
                error_message = '[MUSICA] sub-step budget exhausted'
            end if
            return
        end if

        if (.not. is_ref_state) then
            if (sub_call > 1) then
                call mpas_log_write('[MUSICA] MICM completed time_step $r s in $i sub-calls', &
                    intArgs=[sub_call], realArgs=[real(time_step, RKIND)])
            end if

            if (sub_call > 0) call log_solver_stats(solver_stats)
        end if

    end subroutine micm_substep_loop

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
    subroutine musica_finalize(dminfo)

        use mpas_log, only : mpas_log_write
        use mpas_dmpar, only : mpas_dmpar_sum_int, mpas_dmpar_min_real
        use mpas_derived_types, only : dm_info
        use mpas_kind_types, only : RKIND

        type(dm_info), intent(in) :: dminfo

        integer          :: g_count_in, g_count_out
        real(kind=RKIND) :: l_min_in, l_min_out, g_min_in, g_min_out

        call mpas_log_write('[MUSICA] Finalizing MICM chemistry package...')
        call reset_emission_cache()

        ! Reduce the rank-local clip accumulators once, at this global sync
        ! point, and report the totals. This is the only safe place for a
        ! collective: every rank reaches finalize together.
        call mpas_dmpar_sum_int(dminfo, clip_count_in,  g_count_in)
        call mpas_dmpar_sum_int(dminfo, clip_count_out, g_count_out)

        if (g_count_in > 0 .or. g_count_out > 0) then
            l_min_in  = real(clip_min_in,  kind=RKIND)
            l_min_out = real(clip_min_out, kind=RKIND)
            call mpas_dmpar_min_real(dminfo, l_min_in,  g_min_in)
            call mpas_dmpar_min_real(dminfo, l_min_out, g_min_out)
            call mpas_log_write('[MUSICA] Clipped negative inputs over run: ' // &
                '$i (min $r); outputs: $i (min $r)', &
                intArgs=[g_count_in, g_count_out], &
                realArgs=[g_min_in, g_min_out])
        end if

    end subroutine musica_finalize

    !------------------------------------------------------------------------
    !  subroutine micm_to_mpas_chem
    !
    !> \brief One-time copy of MICM initial state into MPAS tracers
    !> \details
    !>  Loops over chem_species to seed all chemistry tracers generically.
    !------------------------------------------------------------------------
    subroutine micm_to_mpas_chem(scalars, nCellsSolve, nVertLevels, &
            error_code, error_message)

        use mpas_kind_types, only : RKIND
        use musica_util, only : error_t
        use mpas_log, only : mpas_log_write

        real (kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        integer, intent(in) :: nCellsSolve
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
        if (nCellsSolve > 0) then
            target_cells = min(nCellsSolve, size(scalars, 3))
        else
            target_cells = 0
        end if
        micm_levels = min(state%number_of_grid_cells, target_levels)

        if (target_cells <= 0 .or. micm_levels <= 0) then
            call mpas_log_write('[MUSICA] No owned cells/levels for MICM -> MPAS seeding; skipping.')
            return
        end if

        call mpas_log_write('[MUSICA] MICM->MPAS seed: micm_cells=$i target_levels=$i target_cells=$i', &
                intArgs=[state%number_of_grid_cells, micm_levels, target_cells])

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do i_sp = 1, n_chem_species
            ! Host-bound species (e.g. H2O -> qv) must never be seeded: their
            ! MPAS field is owned by the dynamics and seeding the MICM initial
            ! concentration (typically 0) would wipe out, e.g., water vapor.
            if (chem_species(i_sp)%host_bound) cycle

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
    !>  Loops over chem_species to convert MPAS dry-air mass mixing ratios
    !>  (kg/kg dry air) to MICM concentrations (mol/m3) with dry density.
    !>  MICM third-body density is total moist molar density from p/(R*T).
    !------------------------------------------------------------------------
    subroutine MICM_from_chemistry(scalars, rho_dry, temperature, pressure, &
            nCellsSolve, nVertLevels, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND

        real(kind=RKIND), dimension(:,:,:), intent(in) :: scalars      ! (tracer, level, cell)
        real(kind=RKIND), dimension(:,:), intent(in) :: rho_dry        ! (level, cell) [kg/m3]
        real(kind=RKIND), dimension(:,:), intent(in) :: temperature    ! (level, cell) [K]
        real(kind=RKIND), dimension(:,:), intent(in) :: pressure       ! (level, cell) [Pa]
        integer, intent(in) :: nCellsSolve
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
        if (nCellsSolve < 1 .or. nVertLevels < 1 .or. state%number_of_grid_cells < 1) return

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
                micm_cell = (iCell - 1) * nVertLevels + k
                rho = rho_dry(k, iCell)

                ! Set MICM environmental conditions
                state%conditions(micm_cell)%temperature = temperature(k, iCell)
                state%conditions(micm_cell)%pressure = pressure(k, iCell)
                state%conditions(micm_cell)%air_density = &
                    pressure(k, iCell) / (real(CHEM_GAS_R, real64) * temperature(k, iCell))

                ! Convert each species: mixing ratio [kg/kg] -> concentration [mol/m³]
                do i_sp = 1, n_chem_species
                    select case (trim(chem_species(i_sp)%micm_name))
                    case ('O2')
                        if (chem_species(i_sp)%host_bound) then
                            ! MOZART-35 declares O2 as a fixed fraction of dry
                            ! air. Other mechanisms may still transport qO2.
                            conc = real(CHEM_O2_VMR, real64) * rho / real(CHEM_M_AIR, real64)
                        else
                            q_val = scalars(chem_species(i_sp)%mpas_index, k, iCell)
                            if (q_val < 0.0_real64) then
                                clip_count_in = clip_count_in + 1
                                clip_min_in = min(clip_min_in, q_val)
                                q_val = 0.0_real64
                            end if
                            conc = q_val * rho / chem_species(i_sp)%molar_mass
                        end if
                    case ('N2')
                        if (chem_species(i_sp)%host_bound) then
                            conc = real(CHEM_N2_VMR, real64) * rho / real(CHEM_M_AIR, real64)
                        else
                            q_val = scalars(chem_species(i_sp)%mpas_index, k, iCell)
                            if (q_val < 0.0_real64) then
                                clip_count_in = clip_count_in + 1
                                clip_min_in = min(clip_min_in, q_val)
                                q_val = 0.0_real64
                            end if
                            conc = q_val * rho / chem_species(i_sp)%molar_mass
                        end if
                    case default
                        q_val = scalars(chem_species(i_sp)%mpas_index, k, iCell)
                        if (q_val < 0.0_real64) then
                            ! Accumulate rank-locally; reduced/logged once at finalize.
                            clip_count_in = clip_count_in + 1
                            clip_min_in = min(clip_min_in, q_val)
                            q_val = 0.0_real64
                        end if
                        conc = q_val * rho / chem_species(i_sp)%molar_mass
                    end select
                    idx = 1 + (micm_cell - 1) * cell_stride + (chem_species(i_sp)%micm_index - 1) * var_stride
                    state%concentrations(idx) = conc
                end do
            end do
        end do

        call update_mozart35_user_rates(state, error_code, error_message)

    end subroutine MICM_from_chemistry

    !------------------------------------------------------------------------
    !> \brief Evaluate a cell-local MOZART-35 rate that has no native MICM form
    !> \details
    !>  Inputs M and H2O are mol m-3. Returned coefficients use MICM's SI
    !>  molar concentration convention. The expressions are the exact pinned
    !>  source forms; unsupported names fail instead of falling back.
    !------------------------------------------------------------------------
    function musica_mozart35_rate(rate_name, temperature, pressure, air_molar_density, &
                                  h2o_molar_density, error_code) result(rate)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        character(len=*), intent(in) :: rate_name
        real(real64), intent(in) :: temperature, pressure
        real(real64), intent(in) :: air_molar_density, h2o_molar_density
        integer, intent(out) :: error_code
        real(real64) :: rate

        real(real64), parameter :: CGS_PER_SI = real(CHEM_NA, real64) * 1.0e-6_real64
        real(real64) :: mc_cgs, h2oc_cgs, dry_mc_cgs
        real(real64) :: association

        error_code = 0
        rate = 0.0_real64
        if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_real64 .or. &
            .not. ieee_is_finite(pressure) .or. pressure < 0.0_real64 .or. &
            .not. ieee_is_finite(air_molar_density) .or. air_molar_density < 0.0_real64 .or. &
            .not. ieee_is_finite(h2o_molar_density) .or. h2o_molar_density < 0.0_real64 .or. &
            h2o_molar_density > air_molar_density) then
            error_code = 1
            return
        end if

        mc_cgs = air_molar_density * CGS_PER_SI
        h2oc_cgs = h2o_molar_density * CGS_PER_SI
        dry_mc_cgs = max(air_molar_density - h2o_molar_density, 0.0_real64) * CGS_PER_SI

        select case (trim(rate_name))
        case ('R_HO2_HO2')
            rate = (2.3e-13_real64 * exp(600.0_real64 / temperature) + &
                    1.7e-33_real64 * mc_cgs * exp(1000.0_real64 / temperature)) * &
                   (1.0_real64 + 1.4e-21_real64 * h2oc_cgs * &
                    exp(2200.0_real64 / temperature)) * CGS_PER_SI
        case ('R_N2O5_M')
            association = mozart35_troe_rate(temperature, mc_cgs, &
                2.0e-30_real64, -4.4_real64, 1.4e-12_real64, -0.7_real64, 0.60_real64)
            rate = association * (3.333e26_real64 / CGS_PER_SI) * &
                   exp(-10990.0_real64 / temperature)
        case ('R_HNO3_OH')
            rate = (2.4e-14_real64 * exp(460.0_real64 / temperature) + &
                    (6.5e-34_real64 * exp(1335.0_real64 / temperature) * mc_cgs) / &
                    (1.0_real64 + (6.5e-34_real64 * exp(1335.0_real64 / temperature) * &
                     mc_cgs) / (2.7e-17_real64 * exp(2199.0_real64 / temperature)))) * &
                   CGS_PER_SI
        case ('R_HO2NO2_M')
            association = mozart35_troe_rate(temperature, mc_cgs, &
                1.8e-31_real64, -3.2_real64, 4.7e-12_real64, -1.4_real64, 0.60_real64)
            rate = association * (4.761904761904762e26_real64 / CGS_PER_SI) * &
                   exp(-10900.0_real64 / temperature)
        case ('R_PAN_M')
            association = mozart35_troe_rate(temperature, mc_cgs, &
                8.5e-29_real64, -6.5_real64, 1.1e-11_real64, -1.0_real64, 0.60_real64)
            rate = association * (1.111e28_real64 / CGS_PER_SI) * &
                   exp(-14000.0_real64 / temperature)
        case ('R_DMS_OH_B')
            rate = 1.7e-42_real64 * exp(7810.0_real64 / temperature) * &
                   dry_mc_cgs * 0.21_real64 / &
                   (1.0_real64 + 5.5e-31_real64 * exp(7460.0_real64 / temperature) * &
                    dry_mc_cgs * 0.21_real64) * CGS_PER_SI
        case default
            error_code = 1
            return
        end select

        if (.not. ieee_is_finite(rate) .or. rate < 0.0_real64) then
            error_code = 1
            rate = 0.0_real64
        end if

    end function musica_mozart35_rate

    pure function mozart35_troe_rate(temperature, mc_cgs, k0_a, k0_b, &
                                     kinf_a, kinf_b, broadening) result(rate)
        real(real64), intent(in) :: temperature, mc_cgs
        real(real64), intent(in) :: k0_a, k0_b, kinf_a, kinf_b, broadening
        real(real64), parameter :: CGS_PER_SI = real(CHEM_NA, real64) * 1.0e-6_real64
        real(real64) :: rate, k0, kinf, reduced_pressure

        k0 = k0_a * (temperature / 300.0_real64)**k0_b
        kinf = kinf_a * (temperature / 300.0_real64)**kinf_b
        reduced_pressure = max(k0 * mc_cgs / max(kinf, tiny(1.0_real64)), &
                               tiny(1.0_real64))
        rate = (k0 * mc_cgs / (1.0_real64 + reduced_pressure)) * &
               broadening**(1.0_real64 / &
               (1.0_real64 + log10(reduced_pressure)**2)) * CGS_PER_SI
    end function mozart35_troe_rate

    subroutine update_mozart35_user_rates(target_state, error_code, error_message)
        use musica_util, only: error_t

        type(state_t), intent(inout) :: target_state
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: rp, rp_index, rp_cell_stride, rp_var_stride
        integer :: sp_cell_stride, sp_var_stride, h2o_index
        integer :: cell, idx, h2o_storage, local_error
        character(len=:), allocatable :: parameter_name
        real(real64) :: rate, h2o
        logical :: has_mozart_rates
        type(error_t) :: error

        error_code = 0
        error_message = ''
        has_mozart_rates = .false.
        do rp = 1, target_state%rate_parameters_ordering%size()
            parameter_name = trim(target_state%rate_parameters_ordering%name(rp))
            select case (parameter_name)
            case ('USER.R_HO2_HO2', 'USER.R_N2O5_M', 'USER.R_HNO3_OH', &
                  'USER.R_HO2NO2_M', 'USER.R_PAN_M', 'USER.R_DMS_OH_B')
                has_mozart_rates = .true.
            end select
        end do
        if (.not. has_mozart_rates) return

        h2o_index = target_state%species_ordering%index('H2O', error)
        if (.not. error%is_success()) then
            error_code = 1
            error_message = '[MUSICA] MOZART-35 user rates require host-bound H2O'
            return
        end if
        rp_cell_stride = target_state%rate_parameters_strides%grid_cell
        rp_var_stride = target_state%rate_parameters_strides%variable
        sp_cell_stride = target_state%species_strides%grid_cell
        sp_var_stride = target_state%species_strides%variable

        do rp = 1, target_state%rate_parameters_ordering%size()
            parameter_name = trim(target_state%rate_parameters_ordering%name(rp))
            select case (parameter_name)
            case ('USER.R_HO2_HO2', 'USER.R_N2O5_M', 'USER.R_HNO3_OH', &
                  'USER.R_HO2NO2_M', 'USER.R_PAN_M', 'USER.R_DMS_OH_B')
                rp_index = target_state%rate_parameters_ordering%index(parameter_name, error)
                if (.not. error%is_success()) then
                    error_code = 1
                    error_message = '[MUSICA] Could not resolve MOZART-35 rate parameter ' // &
                                    trim(parameter_name)
                    return
                end if
                do cell = 1, target_state%number_of_grid_cells
                    h2o_storage = 1 + (cell - 1) * sp_cell_stride + &
                                  (h2o_index - 1) * sp_var_stride
                    h2o = target_state%concentrations(h2o_storage)
                    rate = musica_mozart35_rate(parameter_name(6:), &
                        target_state%conditions(cell)%temperature, &
                        target_state%conditions(cell)%pressure, &
                        target_state%conditions(cell)%air_density, h2o, local_error)
                    if (local_error /= 0) then
                        error_code = 1
                        error_message = '[MUSICA] Invalid cell-local MOZART-35 rate ' // &
                                        trim(parameter_name)
                        return
                    end if
                    idx = 1 + (cell - 1) * rp_cell_stride + (rp_index - 1) * rp_var_stride
                    target_state%rate_parameters(idx) = rate
                end do
            end select
        end do
    end subroutine update_mozart35_user_rates

    !------------------------------------------------------------------------
    !  subroutine MICM_to_chemistry
    !
    !> \brief Copies MICM state to MPAS chemistry state with unit conversion
    !> \details
    !>  Loops over chem_species to convert MICM concentrations (mol/m3)
    !>  back to MPAS dry-air mass mixing ratios (kg/kg dry air) after the
    !>  solver step, using dry density.
    !------------------------------------------------------------------------
    subroutine MICM_to_chemistry(scalars, rho_dry, nCellsSolve, nVertLevels, &
            error_code, error_message)

        use iso_fortran_env, only: real64
        use mpas_kind_types, only: RKIND

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars   ! (tracer, level, cell)
        real(kind=RKIND), dimension(:,:), intent(in) :: rho_dry        ! (level, cell) [kg/m3]
        integer, intent(in) :: nCellsSolve
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
        if (nCellsSolve < 1 .or. nVertLevels < 1 .or. state%number_of_grid_cells < 1) return

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
                micm_cell = (iCell - 1) * nVertLevels + k
                rho = rho_dry(k, iCell)

                ! Convert each species: concentration [mol/m³] -> mixing ratio [kg/kg]
                do i_sp = 1, n_chem_species
                    ! Read-only host-bound species (e.g. H2O -> qv) are not
                    ! written back: chemistry must not overwrite the dynamical
                    ! water-vapor field with its negligible H2O tendency.
                    if (.not. chem_species(i_sp)%write_back) cycle
                    idx = 1 + (micm_cell - 1) * cell_stride + (chem_species(i_sp)%micm_index - 1) * var_stride
                    conc = state%concentrations(idx)
                    if (conc < 0.0_real64) then
                        ! Accumulate rank-locally; reduced/logged once at finalize.
                        clip_count_out = clip_count_out + 1
                        clip_min_out = min(clip_min_out, conc)
                        conc = 0.0_real64
                    end if
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
        use musica_util, only : error_t

        type(state_t), intent(inout) :: state
        integer, intent(in)          :: number_of_grid_cells
        real(kind=RKIND), intent(in) :: nox_tau     ! [s] NOx sink timescale (0 = disabled)

        integer :: cell_stride, var_stride, rp_index, storage_index, cell, idx, n_rp
        character(len=256) :: msg
        character(len=:), allocatable :: rp_name
        real(real64) :: rate_val
        type(error_t) :: error

        if (state%number_of_rate_parameters == 0) return

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride  = state%rate_parameters_strides%variable
        n_rp = state%rate_parameters_ordering%size()

        do rp_index = 1, n_rp
            rp_name = trim(state%rate_parameters_ordering%name(rp_index))
            storage_index = state%rate_parameters_ordering%index(rp_name, error)
            if (.not. error%is_success()) then
                call mpas_log_write('[MUSICA] Could not resolve storage index for rate param "' // &
                    rp_name // '"')
                cycle
            end if

            ! Set rate based on parameter name.
            ! USER_DEFINED reactions: host multiplier — default 1.0 so the
            ! YAML scaling_factor defines the effective rate constant.
            ! PHOTO parameters stay at 0 unless explicitly configured.
            ! All LOSS.* parameters use config_lnox_nox_tau when enabled.
            if (rp_name(1:min(5,len(rp_name))) == 'USER.') then
                rate_val = 1.0_real64
            else if (rp_name(1:min(5,len(rp_name))) == 'LOSS.') then
                if (nox_tau > 0.0_RKIND) then
                    rate_val = 1.0_real64 / real(nox_tau, real64)
                else
                    rate_val = 0.0_real64
                end if
            else
                rate_val = 0.0_real64
            end if

            do cell = 1, number_of_grid_cells
                idx = 1 + (cell - 1) * cell_stride + (storage_index - 1) * var_stride
                state%rate_parameters(idx) = rate_val
            end do

            write(msg, '(A,A,A,ES10.3)') '[MUSICA] Rate param "', rp_name, '" = ', rate_val
            call mpas_log_write(trim(msg))
        end do

    end subroutine assign_rate_parameters

    !------------------------------------------------------------------------
    !  subroutine musica_cache_emission_indices
    !
    !> \brief Cross-match MIEM species with writable MICM emission reactions
    !> \details
    !>  Resolves the MICM species and ``EMIS.<species>`` rate-parameter
    !>  storage indices once during initialization. The configured MIEM
    !>  species and MICM EMIS parameter sets must match exactly.
    !------------------------------------------------------------------------
    subroutine musica_cache_emission_indices(names, error_code, error_message, net_flux_species)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use musica_util, only: error_t
        use mpas_log, only: mpas_log_write

        character(len=*), intent(in) :: names(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        character(len=*), optional, intent(in) :: net_flux_species

        type(error_t) :: error
        character(len=:), allocatable :: rate_name
        character(len=:), allocatable :: configured_name
        integer :: i, j, chem_index, configured_count, matches

        call reset_emission_cache()
        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            call set_emission_cache_error('emission index cache requested before MICM state initialization', &
                                          error_code, error_message)
            return
        end if
        if (size(names) < 1) then
            call set_emission_cache_error('MIEM reported no emissions species', &
                                          error_code, error_message)
            return
        end if
        if (state%number_of_rate_parameters < 1) then
            call set_emission_cache_error('MICM state has no rate parameters; expected EMIS.* entries', &
                                          error_code, error_message)
            return
        end if

        do i = 1, size(names)
            if (len_trim(names(i)) < 1) then
                call set_emission_cache_error('MIEM reported an empty species name', &
                                              error_code, error_message)
                return
            end if
            if (len_trim(names(i)) > 64) then
                call set_emission_cache_error('MIEM species name exceeds 64 characters: ' // trim(names(i)), &
                                              error_code, error_message)
                return
            end if
            do j = 1, i - 1
                if (trim(names(j)) == trim(names(i))) then
                    call set_emission_cache_error('duplicate MIEM species: ' // trim(names(i)), &
                                                  error_code, error_message)
                    return
                end if
            end do
        end do

        allocate(emission_entries(size(names)))
        n_emission_rp = size(names)

        do i = 1, n_emission_rp
            emission_entries(i)%name = trim(names(i))
            emission_entries(i)%species_index = &
                state%species_ordering%index(trim(names(i)), error)
            if (.not. error%is_success()) then
                call set_emission_cache_error('MIEM species is absent from MICM: ' // trim(names(i)), &
                                              error_code, error_message)
                return
            end if

            rate_name = 'EMIS.' // trim(names(i))
            emission_entries(i)%rate_parameter_index = &
                state%rate_parameters_ordering%index(rate_name, error)
            if (.not. error%is_success()) then
                call set_emission_cache_error('MICM rate parameter is missing: ' // rate_name, &
                                              error_code, error_message)
                return
            end if

            chem_index = 0
            matches = 0
            do j = 1, n_chem_species
                if (trim(chem_species(j)%micm_name) == trim(names(i))) then
                    chem_index = j
                    matches = matches + 1
                end if
            end do
            if (matches /= 1) then
                call set_emission_cache_error('MIEM species has no unique CheMPAS chemistry target: ' // &
                                              trim(names(i)), error_code, error_message)
                return
            end if
            if (chem_species(chem_index)%host_bound .or. &
                .not. chem_species(chem_index)%write_back) then
                call set_emission_cache_error('MIEM species targets a host-bound/read-only field: ' // &
                                              trim(names(i)), error_code, error_message)
                return
            end if

            emission_entries(i)%molar_mass = chem_species(chem_index)%molar_mass
            if (.not. ieee_is_finite(emission_entries(i)%molar_mass) .or. &
                emission_entries(i)%molar_mass <= 0.0_real64) then
                call set_emission_cache_error('invalid molar mass for MIEM species ' // trim(names(i)), &
                                              error_code, error_message)
                return
            end if

        end do

        ! Every MIEM name is now known to have a species and EMIS parameter.
        ! Scan the complete MICM parameter set afterward to reject extras.
        configured_count = 0
        do i = 1, state%rate_parameters_ordering%size()
            configured_name = trim(state%rate_parameters_ordering%name(i))
            if (index(configured_name, 'EMIS.') == 1) then
                if (len_trim(configured_name) <= 5) then
                    call set_emission_cache_error('MICM contains an empty EMIS. rate-parameter suffix', &
                                                  error_code, error_message)
                    return
                end if
                configured_count = configured_count + 1
                matches = 0
                do j = 1, size(names)
                    if (trim(configured_name(6:)) == trim(names(j))) matches = matches + 1
                end do
                if (matches /= 1) then
                    call set_emission_cache_error('MICM emission parameter ' // trim(configured_name) // &
                                                  ' has no one-to-one MIEM species match', &
                                                  error_code, error_message)
                    return
                end if
            end if
        end do
        if (configured_count /= size(names)) then
            call set_emission_cache_error('MIEM species count does not match the MICM EMIS.* parameter count', &
                                          error_code, error_message)
            return
        end if

        if (present(net_flux_species)) then
            call configure_musica_net_flux_species(net_flux_species, error_code, error_message)
            if (error_code /= 0) return
        end if

        do i = 1, n_emission_rp
            rate_name = 'EMIS.' // trim(emission_entries(i)%name)
            call mpas_log_write('[MUSICA] Cached ' // rate_name // &
                ' for MIEM species ' // trim(emission_entries(i)%name) // ' at index $i', &
                intArgs=[emission_entries(i)%rate_parameter_index])
        end do

    end subroutine musica_cache_emission_indices


    subroutine configure_musica_net_flux_species(text, error_code, error_message)

        character(len=*), intent(in) :: text
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=64) :: token
        integer :: comma, first, i, last, matches, stop

        error_code = 0
        error_message = ''
        if (allocated(emission_entries)) emission_entries%allows_net_flux = .false.
        last = len_trim(text)
        if (last == 0) return
        first = 1
        do
            if (first > last) then
                call set_emission_cache_error('config_miem_net_flux_species contains an empty name', &
                                              error_code, error_message)
                return
            end if
            comma = index(text(first:last), ',')
            if (comma == 0) then
                stop = last
            else
                stop = first + comma - 2
            end if
            if (stop < first .or. stop - first + 1 > len(token)) then
                call set_emission_cache_error('invalid config_miem_net_flux_species entry', &
                                              error_code, error_message)
                return
            end if
            token = trim(adjustl(text(first:stop)))
            if (len_trim(token) == 0) then
                call set_emission_cache_error('config_miem_net_flux_species contains an empty name', &
                                              error_code, error_message)
                return
            end if
            matches = 0
            do i = 1, n_emission_rp
                if (trim(token) == trim(emission_entries(i)%name)) then
                    if (emission_entries(i)%allows_net_flux) then
                        call set_emission_cache_error('duplicate config_miem_net_flux_species entry: ' // &
                                                      trim(token), error_code, error_message)
                        return
                    end if
                    emission_entries(i)%allows_net_flux = .true.
                    matches = matches + 1
                end if
            end do
            if (matches /= 1) then
                call set_emission_cache_error('config_miem_net_flux_species is not an exact emission species: ' // &
                                              trim(token), error_code, error_message)
                return
            end if
            if (comma == 0) exit
            first = stop + 2
        end do

    end subroutine configure_musica_net_flux_species

    !------------------------------------------------------------------------
    !  subroutine musica_set_surface_emission_fluxes
    !
    !> \brief Convert surface mass fluxes to lowest-layer MICM source rates
    !> \details
    !>  Validates the complete input and target layout before changing either
    !>  MICM state. Upper-level EMIS parameters are explicitly zeroed and the
    !>  same lowest-layer rates are installed in coupled and reference states.
    !------------------------------------------------------------------------
    subroutine musica_set_surface_emission_fluxes(local_flux, zgrid, nCellsSolve, nVertLevels, &
                                                  error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use mpas_kind_types, only: RKIND

        real(real64), intent(in) :: local_flux(:,:) ! (n_emission_rp, nCellsSolve)
        real(kind=RKIND), intent(in) :: zgrid(:,:)  ! (nVertLevels+1, nCells)
        integer, intent(in) :: nCellsSolve, nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        real(real64), allocatable :: level1_rate(:,:)
        real(real64) :: dz, denominator
        integer :: state_cell_stride, state_var_stride
        integer :: ref_cell_stride, ref_var_stride
        integer :: i, iCell, k, micm_cell, idx, ref_idx
        character(len=256) :: detail

        error_code = 0
        error_message = ''

        if (n_emission_rp < 1) return
        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] Emission flux setter called before MICM state initialization'
            return
        end if
        if (nCellsSolve < 0 .or. nVertLevels < 1) then
            error_code = 1
            error_message = '[MUSICA] Invalid local grid dimensions for emissions'
            return
        end if
        if (size(local_flux, 1) /= n_emission_rp .or. size(local_flux, 2) /= nCellsSolve) then
            error_code = 1
            error_message = '[MUSICA] MIEM local-flux dimensions do not match the cached species/local-cell layout'
            return
        end if
        if (size(zgrid, 1) < nVertLevels + 1 .or. size(zgrid, 2) < nCellsSolve) then
            error_code = 1
            error_message = '[MUSICA] zgrid is too small for MIEM lowest-layer conversion'
            return
        end if
        if (state%number_of_grid_cells /= nCellsSolve * nVertLevels) then
            error_code = 1
            error_message = '[MUSICA] MICM grid-cell count does not match the MPAS emissions layout'
            return
        end if
        if (associated(state_ref)) then
            if (state_ref%number_of_grid_cells /= state%number_of_grid_cells) then
                error_code = 1
                error_message = '[MUSICA] Reference MICM grid does not match the coupled emissions layout'
                return
            end if
        end if

        state_cell_stride = state%rate_parameters_strides%grid_cell
        state_var_stride = state%rate_parameters_strides%variable
        if (state_cell_stride < 1 .or. state_var_stride < 1) then
            error_code = 1
            error_message = '[MUSICA] MICM returned invalid emission rate-parameter strides'
            return
        end if
        if (associated(state_ref)) then
            ref_cell_stride = state_ref%rate_parameters_strides%grid_cell
            ref_var_stride = state_ref%rate_parameters_strides%variable
            if (ref_cell_stride < 1 .or. ref_var_stride < 1) then
                error_code = 1
                error_message = '[MUSICA] Reference MICM returned invalid emission rate-parameter strides'
                return
            end if
        else
            ref_cell_stride = 0
            ref_var_stride = 0
        end if

        allocate(level1_rate(n_emission_rp, nCellsSolve))
        level1_rate = 0.0_real64

        do iCell = 1, nCellsSolve
            dz = real(zgrid(2, iCell) - zgrid(1, iCell), real64)
            if (.not. ieee_is_finite(dz) .or. dz <= 0.0_real64) then
                write(detail, '(A,I0)') '[MUSICA] Invalid lowest-layer thickness for MIEM cell ', iCell
                error_code = 1
                error_message = trim(detail)
                deallocate(level1_rate)
                return
            end if
            do i = 1, n_emission_rp
                if (.not. ieee_is_finite(local_flux(i, iCell)) .or. &
                    (local_flux(i, iCell) < 0.0_real64 .and. &
                     .not. emission_entries(i)%allows_net_flux)) then
                    write(detail, '(A,A,A,I0)') '[MUSICA] Invalid MIEM surface flux for ', &
                        trim(emission_entries(i)%name), ' at local cell ', iCell
                    error_code = 1
                    error_message = trim(detail)
                    deallocate(level1_rate)
                    return
                end if
                denominator = dz * emission_entries(i)%molar_mass
                if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_real64) then
                    error_code = 1
                    error_message = '[MUSICA] Invalid MIEM conversion denominator for ' // &
                                    trim(emission_entries(i)%name)
                    deallocate(level1_rate)
                    return
                end if
                level1_rate(i, iCell) = local_flux(i, iCell) / denominator
                if (.not. ieee_is_finite(level1_rate(i, iCell))) then
                    error_code = 1
                    error_message = '[MUSICA] Nonfinite MICM emission rate for ' // &
                                    trim(emission_entries(i)%name)
                    deallocate(level1_rate)
                    return
                end if
            end do
        end do

        ! Bounds-check every target before changing either state.
        do i = 1, n_emission_rp
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    micm_cell = (iCell - 1) * nVertLevels + k
                    idx = 1 + (micm_cell - 1) * state_cell_stride + &
                          (emission_entries(i)%rate_parameter_index - 1) * state_var_stride
                    if (idx < 1 .or. idx > size(state%rate_parameters)) then
                        error_code = 1
                        error_message = '[MUSICA] Coupled MICM emission-rate buffer index is out of bounds'
                        deallocate(level1_rate)
                        return
                    end if
                    if (associated(state_ref)) then
                        ref_idx = 1 + (micm_cell - 1) * ref_cell_stride + &
                                  (emission_entries(i)%rate_parameter_index - 1) * ref_var_stride
                        if (ref_idx < 1 .or. ref_idx > size(state_ref%rate_parameters)) then
                            error_code = 1
                            error_message = '[MUSICA] Reference MICM emission-rate buffer index is out of bounds'
                            deallocate(level1_rate)
                            return
                        end if
                    end if
                end do
            end do
        end do

        do i = 1, n_emission_rp
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    micm_cell = (iCell - 1) * nVertLevels + k
                    idx = 1 + (micm_cell - 1) * state_cell_stride + &
                          (emission_entries(i)%rate_parameter_index - 1) * state_var_stride
                    state%rate_parameters(idx) = 0.0_real64
                    if (associated(state_ref)) then
                        ref_idx = 1 + (micm_cell - 1) * ref_cell_stride + &
                                  (emission_entries(i)%rate_parameter_index - 1) * ref_var_stride
                        state_ref%rate_parameters(ref_idx) = 0.0_real64
                    end if
                end do

                micm_cell = (iCell - 1) * nVertLevels + 1
                idx = 1 + (micm_cell - 1) * state_cell_stride + &
                      (emission_entries(i)%rate_parameter_index - 1) * state_var_stride
                state%rate_parameters(idx) = level1_rate(i, iCell)
                if (associated(state_ref)) then
                    ref_idx = 1 + (micm_cell - 1) * ref_cell_stride + &
                              (emission_entries(i)%rate_parameter_index - 1) * ref_var_stride
                    state_ref%rate_parameters(ref_idx) = level1_rate(i, iCell)
                end if
            end do
        end do

        deallocate(level1_rate)

    end subroutine musica_set_surface_emission_fluxes

    !------------------------------------------------------------------------
    !  subroutine musica_set_layer_emission_fluxes
    !
    !> \brief Convert per-layer mass fluxes to MICM source rates
    !> \details
    !>  Every source layer is converted independently using its MPAS layer
    !>  depth and species molar mass. The complete input and both MICM target
    !>  layouts are validated before either coupled or reference state changes.
    !------------------------------------------------------------------------
    subroutine musica_set_layer_emission_fluxes(local_layer_flux, zgrid, nCellsSolve, &
                                                nVertLevels, error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use mpas_kind_types, only: RKIND

        real(real64), intent(in) :: local_layer_flux(:,:,:) ! species, level, cell
        real(kind=RKIND), intent(in) :: zgrid(:,:)           ! level interface, cell
        integer, intent(in) :: nCellsSolve, nVertLevels
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        real(real64), allocatable :: layer_rate(:,:,:)
        real(real64) :: dz, denominator
        integer :: state_cell_stride, state_var_stride
        integer :: ref_cell_stride, ref_var_stride
        integer :: i, iCell, k, micm_cell, idx, ref_idx
        character(len=256) :: detail

        error_code = 0
        error_message = ''

        if (n_emission_rp < 1) return
        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] Emission flux setter called before MICM state initialization'
            return
        end if
        if (nCellsSolve < 0 .or. nVertLevels < 1) then
            error_code = 1
            error_message = '[MUSICA] Invalid local grid dimensions for emissions'
            return
        end if
        if (size(local_layer_flux, 1) /= n_emission_rp .or. &
            size(local_layer_flux, 2) /= nVertLevels .or. &
            size(local_layer_flux, 3) /= nCellsSolve) then
            error_code = 1
            error_message = '[MUSICA] MIEM layer-flux dimensions do not match the cached species/local-grid layout'
            return
        end if
        if (size(zgrid, 1) < nVertLevels + 1 .or. size(zgrid, 2) < nCellsSolve) then
            error_code = 1
            error_message = '[MUSICA] zgrid is too small for MIEM layer-flux conversion'
            return
        end if
        if (state%number_of_grid_cells /= nCellsSolve * nVertLevels) then
            error_code = 1
            error_message = '[MUSICA] MICM grid-cell count does not match the MPAS emissions layout'
            return
        end if
        if (associated(state_ref)) then
            if (state_ref%number_of_grid_cells /= state%number_of_grid_cells) then
                error_code = 1
                error_message = '[MUSICA] Reference MICM grid does not match the coupled emissions layout'
                return
            end if
        end if

        state_cell_stride = state%rate_parameters_strides%grid_cell
        state_var_stride = state%rate_parameters_strides%variable
        if (state_cell_stride < 1 .or. state_var_stride < 1) then
            error_code = 1
            error_message = '[MUSICA] MICM returned invalid emission rate-parameter strides'
            return
        end if
        if (associated(state_ref)) then
            ref_cell_stride = state_ref%rate_parameters_strides%grid_cell
            ref_var_stride = state_ref%rate_parameters_strides%variable
            if (ref_cell_stride < 1 .or. ref_var_stride < 1) then
                error_code = 1
                error_message = '[MUSICA] Reference MICM returned invalid emission rate-parameter strides'
                return
            end if
        else
            ref_cell_stride = 0
            ref_var_stride = 0
        end if

        allocate(layer_rate(n_emission_rp, nVertLevels, nCellsSolve))
        layer_rate = 0.0_real64

        do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
                dz = real(zgrid(k + 1, iCell) - zgrid(k, iCell), real64)
                if (.not. ieee_is_finite(dz) .or. dz <= 0.0_real64) then
                    write(detail, '(A,I0,A,I0)') '[MUSICA] Invalid layer thickness for MIEM cell ', &
                                                  iCell, ', level ', k
                    error_code = 1
                    error_message = trim(detail)
                    deallocate(layer_rate)
                    return
                end if
                do i = 1, n_emission_rp
                    if (.not. ieee_is_finite(local_layer_flux(i, k, iCell)) .or. &
                        (local_layer_flux(i, k, iCell) < 0.0_real64 .and. &
                         .not. emission_entries(i)%allows_net_flux)) then
                        write(detail, '(A,A,A,I0,A,I0)') '[MUSICA] Invalid MIEM layer flux for ', &
                            trim(emission_entries(i)%name), ' at local cell ', iCell, ', level ', k
                        error_code = 1
                        error_message = trim(detail)
                        deallocate(layer_rate)
                        return
                    end if
                    denominator = dz * emission_entries(i)%molar_mass
                    if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_real64) then
                        error_code = 1
                        error_message = '[MUSICA] Invalid MIEM conversion denominator for ' // &
                                        trim(emission_entries(i)%name)
                        deallocate(layer_rate)
                        return
                    end if
                    layer_rate(i, k, iCell) = local_layer_flux(i, k, iCell) / denominator
                    if (.not. ieee_is_finite(layer_rate(i, k, iCell))) then
                        error_code = 1
                        error_message = '[MUSICA] Nonfinite MICM emission rate for ' // &
                                        trim(emission_entries(i)%name)
                        deallocate(layer_rate)
                        return
                    end if
                end do
            end do
        end do

        ! Bounds-check every target before mutating either MICM state.
        do i = 1, n_emission_rp
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    micm_cell = (iCell - 1) * nVertLevels + k
                    idx = 1 + (micm_cell - 1) * state_cell_stride + &
                          (emission_entries(i)%rate_parameter_index - 1) * state_var_stride
                    if (idx < 1 .or. idx > size(state%rate_parameters)) then
                        error_code = 1
                        error_message = '[MUSICA] Coupled MICM emission-rate buffer index is out of bounds'
                        deallocate(layer_rate)
                        return
                    end if
                    if (associated(state_ref)) then
                        ref_idx = 1 + (micm_cell - 1) * ref_cell_stride + &
                                  (emission_entries(i)%rate_parameter_index - 1) * ref_var_stride
                        if (ref_idx < 1 .or. ref_idx > size(state_ref%rate_parameters)) then
                            error_code = 1
                            error_message = '[MUSICA] Reference MICM emission-rate buffer index is out of bounds'
                            deallocate(layer_rate)
                            return
                        end if
                    end if
                end do
            end do
        end do

        do i = 1, n_emission_rp
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    micm_cell = (iCell - 1) * nVertLevels + k
                    idx = 1 + (micm_cell - 1) * state_cell_stride + &
                          (emission_entries(i)%rate_parameter_index - 1) * state_var_stride
                    state%rate_parameters(idx) = layer_rate(i, k, iCell)
                    if (associated(state_ref)) then
                        ref_idx = 1 + (micm_cell - 1) * ref_cell_stride + &
                                  (emission_entries(i)%rate_parameter_index - 1) * ref_var_stride
                        state_ref%rate_parameters(ref_idx) = layer_rate(i, k, iCell)
                    end if
                end do
            end do
        end do

        deallocate(layer_rate)

    end subroutine musica_set_layer_emission_fluxes


    subroutine reset_emission_cache()

        if (allocated(emission_entries)) deallocate(emission_entries)
        n_emission_rp = 0

    end subroutine reset_emission_cache

    subroutine set_emission_cache_error(message, error_code, error_message)

        character(len=*), intent(in) :: message
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        call reset_emission_cache()
        error_code = 1
        error_message = '[MUSICA] ' // trim(message)

    end subroutine set_emission_cache_error

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
        use musica_util, only : error_t

        character(len=*), intent(in)                :: names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        integer :: i
        character(len=:), allocatable :: target_name
        type(error_t) :: error

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

        do i = 1, n_photo_rp
            photo_rp_names(i) = trim(names(i))
            target_name = 'PHOTO.' // trim(names(i))
            photo_rp_indices(i) = state%rate_parameters_ordering%index(target_name, error)
            if (.not. error%is_success()) then
                error_code = 1
                error_message = '[MUSICA] Rate parameter ' // target_name // &
                    ' not found; TUV-x and MICM configs disagree'
                return
            end if
            call mpas_log_write('[MUSICA] Cached ' // target_name // &
                ' rate param index = $i', intArgs=[photo_rp_indices(i)])
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
    subroutine musica_set_photolysis_rates(rates, nCellsSolve, nVertLevels, n_rates, &
                                           error_code, error_message)
        use iso_fortran_env, only : real64
        use mpas_kind_types, only : RKIND
        use mpas_log, only : mpas_log_write

        real(kind=RKIND), intent(in) :: rates(:,:,:) ! (n_rates, nVertLevels, nCellsSolve)
        integer, intent(in)          :: nCellsSolve, nVertLevels, n_rates
        integer, intent(out)         :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: cell_stride, var_stride, ref_cell_stride, ref_var_stride
        integer :: iCell, k, micm_cell, idx, ref_idx, ir, rp_idx
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
        if (nCellsSolve < 1) return

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride  = state%rate_parameters_strides%variable
        if (associated(state_ref)) then
            ref_cell_stride = state_ref%rate_parameters_strides%grid_cell
            ref_var_stride  = state_ref%rate_parameters_strides%variable
        end if

        do ir = 1, n_photo_rp
            rp_idx = photo_rp_indices(ir)
            do iCell = 1, nCellsSolve
                do k = 1, nVertLevels
                    micm_cell = (iCell - 1) * nVertLevels + k
                    idx = 1 + (micm_cell - 1) * cell_stride + (rp_idx - 1) * var_stride
                    state%rate_parameters(idx) = real(rates(ir, k, iCell), real64)
                end do
            end do

            if (associated(state_ref)) then
                do iCell = 1, nCellsSolve
                    do k = 1, nVertLevels
                        micm_cell = (iCell - 1) * nVertLevels + k
                        ref_idx = 1 + (micm_cell - 1) * ref_cell_stride + &
                                      (rp_idx - 1) * ref_var_stride
                        state_ref%rate_parameters(ref_idx) = real(rates(ir, k, iCell), real64)
                    end do
                end do
            end if

            rmin = real(minval(rates(ir, 1:nVertLevels, 1:nCellsSolve)), real64)
            rmax = real(maxval(rates(ir, 1:nVertLevels, 1:nCellsSolve)), real64)
            write(msg, '(A,A,A,ES10.3,A,ES10.3)') '[MUSICA] ', &
                trim(photo_rp_names(ir)), ' min=', rmin, ' max=', rmax
            call mpas_log_write(trim(msg))
        end do

    end subroutine musica_set_photolysis_rates

    !------------------------------------------------------------------------
    !  subroutine host_binding
    !
    !> \brief Map a MICM species to an existing MPAS prognostic field
    !> \details
    !>  Some MICM species are not independent chemistry tracers but are
    !>  already carried by MPAS as prognostic fields. Rather than inject a
    !>  redundant advected qXXX tracer, we bind such a species to the host
    !>  field: the coupler reads its concentration from that field each step
    !>  and (optionally) writes the result back.
    !>
    !>  H2O is always bound to qv. O2 and N2 are bound only when their MICM
    !>  species metadata explicitly declares a dry-air host relationship;
    !>  this preserves transported qO2 in the accepted Chapman mechanisms.
    !------------------------------------------------------------------------
    subroutine host_binding(mechanism, micm_name, is_bound, host_name, do_write_back, fixed_vmr)

        use musica_util, only: error_t, string_t

        type(micm_t), pointer, intent(inout) :: mechanism
        character(len=*), intent(in)  :: micm_name
        logical,          intent(out) :: is_bound
        character(len=*), intent(out) :: host_name
        logical,          intent(out) :: do_write_back
        real(real64), optional, intent(out) :: fixed_vmr

        type(error_t) :: property_error
        type(string_t) :: binding

        is_bound      = .false.
        host_name     = ''
        do_write_back = .true.
        if (present(fixed_vmr)) fixed_vmr = -1.0_real64

        select case (trim(micm_name))
        case ('H2O')
            is_bound      = .true.
            host_name     = 'qv'     ! resolve_mpas_indices looks up index_qv
            do_write_back = .false.  ! read-only environmental parameter
        case ('O2', 'N2')
            binding = mechanism%get_species_property_string( &
                trim(micm_name), '__host binding', property_error)
            if (property_error%is_success()) then
                if (index(trim(binding%value_), 'dry-air mole fraction') == 1) then
                    is_bound      = .true.
                    host_name     = 'qv'  ! index anchor; concentration derives from dry air
                    do_write_back = .false.
                    if (present(fixed_vmr)) then
                        if (trim(micm_name) == 'O2') then
                            fixed_vmr = real(CHEM_O2_VMR, real64)
                        else
                            fixed_vmr = real(CHEM_N2_VMR, real64)
                        end if
                    end if
                end if
            end if
        end select
    end subroutine host_binding

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
        state_ref%rate_parameters(:) = state%rate_parameters(:)

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
