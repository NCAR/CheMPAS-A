! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_tuvx
!
!> \brief TUV-x photolysis rate calculations for CheMPAS-A
!> \author CheMPAS-A Developers
!> \date   March 2026
!> \details
!>  Wraps the MUSICA TUV-x Fortran API to compute per-column photolysis
!>  rate constants from atmospheric profiles provided by the host model.
!>  All atmospheric profiles (height, air, temperature, O3, O2) come
!>  from MPAS state — no static data files for atmospheric state.
!>  Spectral physics data (cross-sections, quantum yields, solar flux)
!>  are read from external files by TUV-x.
!>
!>  When config_tuvx_config_file is empty, this module is a no-op and
!>  the chemistry driver falls back to Phase 1 cos(SZA) scaling.
!
!-----------------------------------------------------------------------
module mpas_tuvx

#ifdef MPAS_USE_MUSICA
    use musica_util, only: dk => musica_dk, error_t, mappings_t
    use musica_tuvx, only: tuvx_t
    use musica_tuvx_grid, only: grid_t
    use musica_tuvx_grid_map, only: grid_map_t
    use musica_tuvx_profile, only: profile_t
    use musica_tuvx_profile_map, only: profile_map_t
    use musica_tuvx_radiator, only: radiator_t
    use musica_tuvx_radiator_map, only: radiator_map_t
#endif

    use mpas_kind_types, only: RKIND
    use mpas_chem_constants, only: CHEM_GAS_R, CHEM_NA, CHEM_M_O3, CHEM_O2_VMR, CHEM_RHO_H2O
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private

    public :: tuvx_init, tuvx_compute_photolysis, tuvx_is_enabled, tuvx_finalize
    public :: tuvx_query_photo_rates
    public :: tuvx_n_rates, tuvx_rate_names
    public :: tuvx_uses_spatial_climatology, tuvx_n_extension_layers

#ifdef MPAS_USE_MUSICA
    ! Module-level state (initialized once, updated per column)
    type(tuvx_t),       pointer :: tuvx_solver   => null()
    type(grid_t),       pointer :: height_grid   => null()
    type(grid_t),       pointer :: cloud_wavelength_grid => null()
    type(profile_t),    pointer :: air_profile   => null()
    type(profile_t),    pointer :: temp_profile  => null()
    type(profile_t),    pointer :: o3_profile    => null()
    type(profile_t),    pointer :: o2_profile    => null()

    ! Cloud radiator (Phase 3)
    type(radiator_t),   pointer :: cloud_radiator => null()

    ! Wavelength-grid section count the from-host cloud radiator is sized for.
    ! The CAM wavelength grid (cam.csv) has 103 edges = 102 sections. The
    ! from-host cloud radiator cannot be resized after construction, so the
    ! configured wavelength grid must have exactly this many sections.
    integer, parameter :: N_CLOUD_WAVELENGTH_BINS = 102  ! CAM wavelength grid: 103 edges - 1
    integer, save :: n_wavelength_bins = 0

    real(dk), allocatable :: cloud_od(:,:)   ! (n_total_layers, n_wavelength_bins)
    real(dk), allocatable :: cloud_ssa(:,:)  ! (n_total_layers, n_wavelength_bins)
    real(dk), allocatable :: cloud_g(:,:,:)  ! (n_total_layers, n_wavelength_bins, 1)

    ! Cloud optical property constants (warm cloud, wavelength-independent)
    real(dk), parameter :: CLOUD_SSA_VAL = 0.999999_dk  ! pure scattering
    real(dk), parameter :: CLOUD_G_VAL   = 0.85_dk      ! strong forward scattering

    ! Cloud microphysics constants
    real(dk), parameter :: R_EFF_CLOUD = 10.0e-6_dk   ! cloud droplet r_eff [m]
    real(dk), parameter :: R_EFF_RAIN  = 500.0e-6_dk  ! rain drop r_eff [m]

    ! Photolysis rate bookkeeping: cache all rates TUV-x reports.
    ! tuvx_indices(i) is the TUV-x ordering index for photo_names(i).
    integer, save              :: n_photo_rates   = 0
    integer, save              :: n_heating_rates = 0
    character(len=64), allocatable :: photo_names(:)  ! one per TUV-x reaction
    integer,           allocatable :: tuvx_indices(:)

    ! Column decomposition:
    !   n_mpas_layers   : layers supplied by MPAS host state
    !   n_ext_layers    : climatology layers above MPAS top (0 when disabled)
    !   n_total_layers  : full column TUV-x sees (= n_mpas + n_ext)
    integer, save :: n_mpas_layers   = 0
    integer, save :: n_ext_layers    = 0
    integer, save :: n_total_layers  = 0
    character(len=32), save :: upper_column_mode = 'none'

    ! Upper-atmosphere climatology edge values (sized n_ext_layers+1 when enabled).
    ! Values are defined at edge heights; midpoints are computed on the fly.
    real(dk), allocatable :: ext_z_km_edges(:)
    real(dk), allocatable :: ext_T_K_edges(:)
    real(dk), allocatable :: ext_n_air_edges(:)
    real(dk), allocatable :: ext_n_o3_edges(:)

    real(dk), allocatable :: photo_rates(:,:)    ! (n_total_layers+1, n_photo_rates)
    real(dk), allocatable :: heating_rates(:,:)  ! (n_total_layers+1, n_heating_rates)

    ! Negative host mixing ratios are clipped only at the TUV-x boundary.
    ! Accumulate rank-local statistics here and reduce them at finalize, where
    ! every rank is guaranteed to participate in collectives.
    integer,          save :: clip_count_input = 0
    real(kind=RKIND), save :: clip_min_input = huge(1.0_RKIND)

#endif

contains

    !------------------------------------------------------------------------
    !  subroutine tuvx_init
    !
    !> \brief Initialize TUV-x with from-host grids and profiles
    !> \details
    !>  Creates height grid and atmospheric profiles (air, temperature, O3, O2)
    !>  as from-host objects, registers them, constructs the TUV-x solver,
    !>  then retrieves updatable handles. Looks up the j_NO2 reaction index
    !>  and allocates persistent work arrays.
    !------------------------------------------------------------------------
    subroutine tuvx_init(config_file, nVertLevels_in, requested_upper_column_mode, &
                         extension_file, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif

        character(len=*), intent(in)               :: config_file
        integer, intent(in)                        :: nVertLevels_in
        character(len=*), intent(in)               :: requested_upper_column_mode
        character(len=*), intent(in)               :: extension_file
        integer, intent(out)                       :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

#ifdef MPAS_USE_MUSICA
        type(grid_map_t),     pointer :: grids    => null()
        type(profile_map_t),  pointer :: profiles => null()
        type(radiator_map_t), pointer :: radiators => null()
        type(mappings_t),     pointer :: photo_ordering  => null()
        type(mappings_t),     pointer :: heat_ordering   => null()
        type(error_t) :: err
        integer :: i
        character(len=256) :: msg
        real(kind=RKIND), allocatable :: spatial_z_edges_m(:)
        real(kind=RKIND), allocatable :: spatial_temperature_edges(:)
        real(kind=RKIND), allocatable :: spatial_air_edges(:)

        error_code = 0
        error_message = ''

        clip_count_input = 0
        clip_min_input = huge(1.0_RKIND)

        if (len_trim(config_file) == 0) then
            call mpas_log_write('[TUV-x] No config file specified; TUV-x disabled (Phase 1 fallback).')
            return
        end if

        n_mpas_layers = nVertLevels_in
        n_ext_layers  = 0
        upper_column_mode = trim(adjustl(requested_upper_column_mode))

        call mpas_log_write('[TUV-x] Initializing with config: ' // trim(config_file))
        write(msg, '(A,I0)') '[TUV-x] MPAS nVertLevels = ', n_mpas_layers
        call mpas_log_write(trim(msg))

        ! Select exactly one upper-column mode. The prescribed-field provider
        ! owns spatial/time I/O; this TUV-x adapter owns only the column seam.
        select case (trim(upper_column_mode))
        case ('none')
            if (len_trim(extension_file) > 0) then
                error_code = 1
                error_message = '[TUV-x] upper-column mode none requires an empty legacy extension file.'
                return
            end if
        case ('legacy_static')
            if (len_trim(extension_file) == 0) then
                error_code = 1
                error_message = '[TUV-x] legacy_static requires config_tuvx_extension_file.'
                return
            end if
            call load_extension_csv(trim(extension_file), error_code, error_message)
            if (error_code /= 0) return
            write(msg, '(A,I0,A,A)') '[TUV-x] Extension enabled: n_ext_layers = ', n_ext_layers, &
                                     ' from ', trim(extension_file)
            call mpas_log_write(trim(msg))
        case ('spatial_climatology')
            block
                use mpas_prescribed_fields, only: prescribed_fields_is_enabled, &
                                                  prescribed_fields_get_structure
                if (.not. prescribed_fields_is_enabled()) then
                    error_code = 1
                    error_message = '[TUV-x] spatial_climatology provider is not initialized.'
                    return
                end if
                call prescribed_fields_get_structure(spatial_z_edges_m, &
                                                     spatial_temperature_edges, &
                                                     spatial_air_edges, &
                                                     error_code, error_message)
                if (error_code /= 0) return
                if (size(spatial_z_edges_m) < 2 .or. &
                    size(spatial_temperature_edges) /= size(spatial_z_edges_m) .or. &
                    size(spatial_air_edges) /= size(spatial_z_edges_m)) then
                    error_code = 1
                    error_message = '[TUV-x] spatial climatology structure has inconsistent dimensions.'
                    return
                end if
                n_ext_layers = size(spatial_z_edges_m) - 1
                allocate(ext_z_km_edges(n_ext_layers + 1), &
                         ext_T_K_edges(n_ext_layers + 1), &
                         ext_n_air_edges(n_ext_layers + 1))
                ext_z_km_edges = real(spatial_z_edges_m, dk) * 1.0e-3_dk
                ext_T_K_edges = real(spatial_temperature_edges, dk)
                ext_n_air_edges = real(spatial_air_edges, dk)
                deallocate(spatial_z_edges_m, spatial_temperature_edges, spatial_air_edges)
                write(msg, '(A,I0)') '[TUV-x] Spatial climatology extension layers = ', n_ext_layers
                call mpas_log_write(trim(msg))
            end block
        case default
            error_code = 1
            error_message = '[TUV-x] config_tuvx_upper_column_mode must be none, legacy_static, or spatial_climatology.'
            return
        end select

        n_total_layers = n_mpas_layers + n_ext_layers

        write(msg, '(A,I0)') '[TUV-x] Total column layers (grid sections) = ', n_total_layers
        call mpas_log_write(trim(msg))

        ! Create from-host height grid sized for the full composite column
        height_grid => grid_t("height", "km", n_total_layers, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Create grid map and register height grid
        grids => grid_map_t(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call grids%add(height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Create from-host profiles on the height grid
        air_profile  => profile_t("air",         "molecule cm-3", height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        temp_profile => profile_t("temperature", "K",             height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        o3_profile   => profile_t("O3",          "molecule cm-3", height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        o2_profile   => profile_t("O2",          "molecule cm-3", height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Create profile map and register all profiles
        profiles => profile_map_t(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call profiles%add(air_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call profiles%add(temp_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call profiles%add(o3_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call profiles%add(o2_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Create from-host cloud radiator (Phase 3)
        ! The cloud radiator needs a wavelength grid for array sizing.
        ! The CAM wavelength grid (cam.csv) has 103 edges = 102 sections.
        ! Keep this construction-time grid alive for the lifetime of the
        ! cloud radiator; freeing it during init leaves the radiator built
        ! on a deleted grid object.
        n_wavelength_bins = N_CLOUD_WAVELENGTH_BINS  ! CAM wavelength grid: 103 edges - 1

        cloud_wavelength_grid => grid_t("wavelength", "nm", n_wavelength_bins, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        cloud_radiator => radiator_t("clouds", height_grid, cloud_wavelength_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        write(msg, '(A,I0)') '[TUV-x] Cloud radiator created with n_wavelength_bins = ', n_wavelength_bins
        call mpas_log_write(trim(msg))

        ! Create radiator map with the cloud radiator
        radiators => radiator_map_t(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call radiators%add(cloud_radiator, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Construct TUV-x solver
        tuvx_solver => tuvx_t(trim(config_file), grids, profiles, radiators, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Retrieve updatable handles after construction
        grids    => tuvx_solver%get_grids(err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! NOTE (P2-F10 / FABLE R4): a runtime check that the configured
        ! wavelength grid actually has N_CLOUD_WAVELENGTH_BINS sections is not
        ! feasible with the current MUSICA API. The wavelength grid lives inside
        ! the cloud radiator, not the host grid map, so grid_map%get(
        ! "wavelength", "nm") raises MUSICA Error 105 (grid type mismatch). All
        ! in-tree TUV-x configs use the 102-section cam.csv wavelength grid; a
        ! mismatch would still surface inside the TUV-x core. Deriving the bin
        ! count remains open pending a MUSICA API to query the radiator grid.

        profiles => tuvx_solver%get_profiles(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        radiators => tuvx_solver%get_radiators(err)
        if (tuvx_has_error(err, error_code, error_message)) return

        height_grid  => grids%get("height", "km", err)
        if (tuvx_has_error(err, error_code, error_message)) return
        air_profile  => profiles%get("air", "molecule cm-3", err)
        if (tuvx_has_error(err, error_code, error_message)) return
        temp_profile => profiles%get("temperature", "K", err)
        if (tuvx_has_error(err, error_code, error_message)) return
        o3_profile   => profiles%get("O3", "molecule cm-3", err)
        if (tuvx_has_error(err, error_code, error_message)) return
        o2_profile   => profiles%get("O2", "molecule cm-3", err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Retrieve solver-owned cloud radiator handle
        cloud_radiator => radiators%get("clouds", err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call mpas_log_write('[TUV-x] Cloud radiator registered and retrieved from solver.')

        ! Allocate persistent cloud optical property arrays
        allocate(cloud_od(n_total_layers, n_wavelength_bins))
        allocate(cloud_ssa(n_total_layers, n_wavelength_bins))
        allocate(cloud_g(n_total_layers, n_wavelength_bins, 1))

        cloud_od  = 0.0_dk                ! extension rows stay zero for all time
        cloud_ssa = CLOUD_SSA_VAL
        cloud_g(:,:,1) = CLOUD_G_VAL

        ! Enumerate ALL photolysis reactions TUV-x reports; cache names and indices.
        photo_ordering => tuvx_solver%get_photolysis_rate_constants_ordering(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        n_photo_rates = photo_ordering%size()

        if (n_photo_rates < 1) then
            error_code = 1
            error_message = '[TUV-x] TUV-x config registered zero photolysis reactions'
            deallocate(photo_ordering)
            return
        end if

        if (allocated(photo_names))  deallocate(photo_names)
        if (allocated(tuvx_indices)) deallocate(tuvx_indices)
        allocate(photo_names(n_photo_rates))
        allocate(tuvx_indices(n_photo_rates))

        do i = 1, n_photo_rates
            photo_names(i)  = trim(photo_ordering%name(i))
            tuvx_indices(i) = photo_ordering%index(i)
            write(msg, '(A,A,A,I0)') '[TUV-x] Photolysis reaction ', &
                trim(photo_names(i)), ' index = ', tuvx_indices(i)
            call mpas_log_write(trim(msg))
        end do
        deallocate(photo_ordering)

        ! Query heating rates ordering
        heat_ordering => tuvx_solver%get_heating_rates_ordering(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        n_heating_rates = heat_ordering%size()
        deallocate(heat_ordering)

        ! Allocate persistent work arrays (n_total_layers+1 edges for TUV-x output)
        allocate(photo_rates(n_total_layers + 1, max(n_photo_rates, 1)))
        allocate(heating_rates(n_total_layers + 1, max(n_heating_rates, 1)))
        photo_rates   = 0.0_dk
        heating_rates = 0.0_dk

        write(msg, '(A,I0,A,I0)') '[TUV-x] Initialized: n_photo=', n_photo_rates, &
            ' n_heating=', n_heating_rates
        call mpas_log_write(trim(msg))
#else
        error_code = 0
        error_message = ''
#endif

    end subroutine tuvx_init


    !------------------------------------------------------------------------
    !  subroutine tuvx_compute_photolysis
    !
    !> \brief Compute photolysis rates for a single atmospheric column
    !> \details
    !>  Converts MPAS fields to TUV-x units, updates from-host profiles,
    !>  calls TUV-x run, and extracts every cached photolysis rate
    !>  (edge values averaged to midpoints) for the MPAS slice. The air
    !>  number-density profile uses total moist molar density from p/(R*T);
    !>  O3 number density uses dry density because MPAS tracers are dry-air
    !>  mass mixing ratios.
    !------------------------------------------------------------------------
    subroutine tuvx_compute_photolysis(cos_sza, z_edges_m, rho_dry, temperature, &
                                        pressure, qo3, qc, qr, nVertLevels_in, &
                                        photo_rates_out, error_code, error_message, &
                                        prescribed_o3_density)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif

        real(kind=RKIND), intent(in)  :: cos_sza
        real(kind=RKIND), intent(in)  :: z_edges_m(:)       ! nVertLevels+1 edges [m]
        real(kind=RKIND), intent(in)  :: rho_dry(:)          ! nVertLevels dry density [kg/m3]
        real(kind=RKIND), intent(in)  :: temperature(:)      ! nVertLevels midpoints [K]
        real(kind=RKIND), intent(in)  :: pressure(:)         ! nVertLevels total pressure [Pa]
        real(kind=RKIND), intent(in)  :: qo3(:)              ! nVertLevels [kg/kg]
        real(kind=RKIND), intent(in)  :: qc(:)               ! nVertLevels cloud water [kg/kg]
        real(kind=RKIND), intent(in)  :: qr(:)               ! nVertLevels rain water [kg/kg]
        integer, intent(in)           :: nVertLevels_in
        real(kind=RKIND), intent(out) :: photo_rates_out(:,:) ! (n_rates, nVertLevels)
        integer, intent(out)          :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        real(kind=RKIND), intent(in), optional :: prescribed_o3_density(:) ! extension layers [molec cm-3]

#ifdef MPAS_USE_MUSICA
        real(dk) :: z_km(n_total_layers + 1), z_mid_km(n_total_layers)
        real(dk) :: air_nd(n_total_layers), air_edges(n_total_layers + 1)
        real(dk) :: temp_mid(n_total_layers), temp_edges(n_total_layers + 1)
        real(dk) :: o3_nd(n_total_layers), o2_nd(n_total_layers)
        real(dk) :: o3_edges(n_total_layers + 1), o2_edges(n_total_layers + 1)
        real(dk) :: air_ld(n_total_layers), o3_ld(n_total_layers), o2_ld(n_total_layers)
        real(dk) :: dz_cm(n_total_layers), dz_m_mpas(n_mpas_layers)
        real(dk) :: cloud_tau_1d(n_mpas_layers)
        real(dk) :: sza_rad, esd
        type(error_t) :: err
        integer :: k, j, iWave, new_clip_count
        real(kind=RKIND) :: new_clip_min

        error_code = 0
        error_message = ''

        if (nVertLevels_in /= n_mpas_layers .or. size(z_edges_m) /= n_mpas_layers + 1 .or. &
            size(rho_dry) /= n_mpas_layers .or. size(temperature) /= n_mpas_layers .or. &
            size(pressure) /= n_mpas_layers .or. size(qo3) /= n_mpas_layers .or. &
            size(qc) /= n_mpas_layers .or. size(qr) /= n_mpas_layers .or. &
            size(photo_rates_out, 1) < n_photo_rates .or. &
            size(photo_rates_out, 2) /= n_mpas_layers) then
            error_code = 1
            error_message = '[TUV-x] host column arrays do not match the initialized solver dimensions.'
            return
        end if
        if (.not. ieee_is_finite(cos_sza) .or. cos_sza < -1.0_RKIND .or. cos_sza > 1.0_RKIND .or. &
            .not. all(ieee_is_finite(z_edges_m)) .or. &
            .not. all(ieee_is_finite(rho_dry)) .or. &
            .not. all(ieee_is_finite(temperature)) .or. &
            .not. all(ieee_is_finite(pressure)) .or. .not. all(ieee_is_finite(qo3)) .or. &
            .not. all(ieee_is_finite(qc)) .or. .not. all(ieee_is_finite(qr)) .or. &
            any(z_edges_m(2:) <= z_edges_m(:n_mpas_layers)) .or. any(rho_dry <= 0.0_RKIND) .or. &
            any(temperature <= 0.0_RKIND) .or. any(pressure <= 0.0_RKIND)) then
            error_code = 1
            error_message = '[TUV-x] host column contains invalid coordinates or state values.'
            return
        end if
        if (trim(upper_column_mode) == 'spatial_climatology') then
            if (.not. present(prescribed_o3_density)) then
                error_code = 1
                error_message = '[TUV-x] spatial_climatology requires a per-cell O3 profile.'
                return
            end if
            if (size(prescribed_o3_density) /= n_ext_layers .or. &
                .not. all(ieee_is_finite(prescribed_o3_density)) .or. &
                any(prescribed_o3_density < 0.0_RKIND)) then
                error_code = 1
                error_message = '[TUV-x] prescribed per-cell O3 profile is invalid.'
                return
            end if
        else if (present(prescribed_o3_density)) then
            error_code = 1
            error_message = '[TUV-x] prescribed O3 was supplied outside spatial_climatology mode.'
            return
        end if

        ! Transport and microphysics can leave roundoff-scale negative mixing
        ! ratios. Match the MICM coupling policy: retain the prognostic inputs,
        ! clip only the values passed across this package boundary, and report
        ! aggregate statistics once at finalize.
        new_clip_count = count(qo3 < 0.0_RKIND) + count(qc < 0.0_RKIND) + &
                         count(qr < 0.0_RKIND)
        if (new_clip_count > 0) then
            new_clip_min = min(minval(qo3, mask=qo3 < 0.0_RKIND), &
                               minval(qc,  mask=qc  < 0.0_RKIND), &
                               minval(qr,  mask=qr  < 0.0_RKIND))
            clip_count_input = clip_count_input + new_clip_count
            clip_min_input = min(clip_min_input, new_clip_min)
        end if

        ! Nighttime: skip TUV-x call, return zero for every cached rate
        if (cos_sza <= 0.0_RKIND) then
            photo_rates_out(1:n_photo_rates, 1:n_mpas_layers) = 0.0_RKIND
            return
        end if

        ! ---- Build the composite column (MPAS below + climatology above) ----

        ! Height edges [km]: MPAS edges 1..n_mpas+1, then climatology edges
        ! 2..n_ext+1 (skipping ext edge 1 which coincides with MPAS top).
        z_km(1:n_mpas_layers + 1) = real(z_edges_m(1:n_mpas_layers + 1), dk) * 1.0e-3_dk
        if (n_ext_layers > 0) then
            block
                real(dk) :: model_top_km
                character(len=256) :: msg

                model_top_km = z_km(n_mpas_layers + 1)
                if (abs(ext_z_km_edges(1) - model_top_km) > 5.0e-4_dk) then
                    error_code = 1
                    write(msg, '(A,F0.2,A,F0.2,A)') &
                        '[TUV-x] Extension file first edge (', ext_z_km_edges(1), &
                        ' km) does not match model top (', model_top_km, &
                        ' km); regenerate the upper-column input for this mesh.'
                    error_message = trim(msg)
                    return
                end if
            end block
        end if
        do j = 1, n_ext_layers
            z_km(n_mpas_layers + 1 + j) = ext_z_km_edges(j + 1)
        end do
        z_mid_km = 0.5_dk * (z_km(1:n_total_layers) + z_km(2:n_total_layers + 1))

        ! Layer thickness in cm (for layer densities)
        dz_cm = (z_km(2:n_total_layers + 1) - z_km(1:n_total_layers)) * 1.0e5_dk

        ! MPAS-layer midpoint values from host state
        temp_mid(1:n_mpas_layers) = real(temperature(1:n_mpas_layers), dk)
        air_nd (1:n_mpas_layers) = real(pressure(1:n_mpas_layers), dk) / &
                                   (real(CHEM_GAS_R, dk) * temp_mid(1:n_mpas_layers)) * &
                                   real(CHEM_NA, dk) * 1.0e-6_dk
        o3_nd  (1:n_mpas_layers) = real(max(qo3(1:n_mpas_layers), 0.0_RKIND), dk) &
                                   * real(rho_dry(1:n_mpas_layers), dk) &
                                   * (real(CHEM_NA, dk) / real(CHEM_M_O3, dk)) * 1.0e-6_dk

        ! Extension-layer midpoint values from CSV (average adjacent edges)
        ! Anchor the extension T profile to the MPAS top-midpoint T: keep
        ! the climatology's vertical gradient (stratopause warming,
        ! mesopause cooling) but shift so the join is continuous. Avoids
        ! an ~30 K step between the Weisman-Klemp sounding (no stratopause)
        ! and the US-Std-Atm-based CSV at 50 km.
        if (n_ext_layers > 0) then
        block
            real(dk) :: T_anchor_offset, ext_first_mid_T
            ext_first_mid_T = 0.5_dk * (ext_T_K_edges(1) + ext_T_K_edges(2))
            T_anchor_offset = temp_mid(n_mpas_layers) - ext_first_mid_T

            do j = 1, n_ext_layers
                air_nd (n_mpas_layers + j) = 0.5_dk * (ext_n_air_edges(j) + ext_n_air_edges(j + 1))
                temp_mid(n_mpas_layers + j) = 0.5_dk * (ext_T_K_edges(j)   + ext_T_K_edges(j + 1)) &
                                              + T_anchor_offset
            end do
            if (trim(upper_column_mode) == 'spatial_climatology') then
                o3_nd(n_mpas_layers + 1:n_total_layers) = real(prescribed_o3_density, dk)
            else
                do j = 1, n_ext_layers
                    o3_nd(n_mpas_layers + j) = 0.5_dk * (ext_n_o3_edges(j) + ext_n_o3_edges(j + 1))
                end do
            end if
        end block
        end if

        ! O2 = 0.2095 × air [molecule/cm³] across the full column
        o2_nd = real(CHEM_O2_VMR, dk) * air_nd

        ! Edge values (simple extrapolation from adjacent midpoints).
        ! At the MPAS/extension join (edge n_mpas+1), blend the top MPAS midpoint
        ! with the first extension midpoint instead of the plain extrapolation.
        ! Compute edges for all midpoint-defined profiles: air, T, O3, O2.
        air_edges (1) = air_nd (1)
        temp_edges(1) = temp_mid(1)
        o3_edges  (1) = o3_nd  (1)
        o2_edges  (1) = o2_nd  (1)
        do k = 2, n_total_layers
            air_edges (k) = 0.5_dk * (air_nd (k - 1) + air_nd (k))
            temp_edges(k) = 0.5_dk * (temp_mid(k - 1) + temp_mid(k))
            o3_edges  (k) = 0.5_dk * (o3_nd  (k - 1) + o3_nd  (k))
            o2_edges  (k) = 0.5_dk * (o2_nd  (k - 1) + o2_nd  (k))
        end do
        air_edges (n_total_layers + 1) = air_nd (n_total_layers)
        temp_edges(n_total_layers + 1) = temp_mid(n_total_layers)
        o3_edges  (n_total_layers + 1) = o3_nd  (n_total_layers)
        o2_edges  (n_total_layers + 1) = o2_nd  (n_total_layers)

        ! Layer thickness in metres (MPAS slice, for cloud OD only)
        dz_m_mpas = (z_km(2:n_mpas_layers + 1) - z_km(1:n_mpas_layers)) * 1.0e3_dk

        ! Layer densities [molecule/cm²] = number_density × thickness [cm]
        air_ld = air_nd * dz_cm
        o3_ld  = o3_nd  * dz_cm
        o2_ld  = o2_nd  * dz_cm

        ! Update from-host height grid
        call height_grid%set_edges(z_km, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call height_grid%set_midpoints(z_mid_km, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Update from-host profiles
        call air_profile%set_midpoint_values(air_nd, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call air_profile%set_edge_values(air_edges, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call air_profile%set_layer_densities(air_ld, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Populate the "exo" layer density above the top of the grid so
        ! spherical-geometry slant-path air mass accounts for air above
        ! the extension top (100 km) rather than treating it as vacuum.
        ! Use a scale-height approximation: layer_dens_exo = n_air_top *
        ! scale_height. Pick H = 7 km, representative of the mesosphere/
        ! lower-thermosphere around the extension top (T ≈ 195 K →
        ! H = R T / (M_air g) ≈ 5.7 km; rounded up).
        call air_profile%calculate_exo_layer_density(7.0_dk, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        call temp_profile%set_midpoint_values(temp_mid, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call temp_profile%set_edge_values(temp_edges, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        call o3_profile%set_midpoint_values(o3_nd, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call o3_profile%set_edge_values(o3_edges, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call o3_profile%set_layer_densities(o3_ld, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        call o2_profile%set_midpoint_values(o2_nd, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call o2_profile%set_edge_values(o2_edges, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call o2_profile%set_layer_densities(o2_ld, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Compute cloud optical depth per layer and update cloud radiator.
        ! Clouds only live in the MPAS slice; extension rows of cloud_od were
        ! zero-filled at init and are never touched here.
        call compute_cloud_optical_depth( &
            real(max(qc(1:n_mpas_layers), 0.0_RKIND), dk), &
            real(max(qr(1:n_mpas_layers), 0.0_RKIND), dk), &
            real(rho_dry(1:n_mpas_layers), dk), dz_m_mpas, n_mpas_layers, cloud_tau_1d)

        do iWave = 1, n_wavelength_bins
            cloud_od(1:n_mpas_layers, iWave) = cloud_tau_1d(:)
        end do

        call cloud_radiator%set_optical_depths(cloud_od, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call cloud_radiator%set_single_scattering_albedos(cloud_ssa, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call cloud_radiator%set_asymmetry_factors(cloud_g, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! SZA in radians; earth-sun distance ~1 AU
        sza_rad = acos(real(cos_sza, dk))
        esd = 1.0_dk

        ! Run TUV-x
        call tuvx_solver%run(sza_rad, esd, photo_rates, heating_rates, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Extract every cached photolysis rate on the MPAS slice; extension-
        ! layer values are computed by TUV-x but discarded (the host has no
        ! storage for them).
        block
            integer :: ir
            do ir = 1, n_photo_rates
                do k = 1, n_mpas_layers
                    photo_rates_out(ir, k) = real( &
                        0.5_dk * (photo_rates(k,     tuvx_indices(ir)) + &
                                  photo_rates(k + 1, tuvx_indices(ir))), RKIND)
                end do
            end do
        end block
#else
        error_code = 0
        error_message = ''
        photo_rates_out = 0.0_RKIND
#endif

    end subroutine tuvx_compute_photolysis


    !------------------------------------------------------------------------
    !  function tuvx_n_rates
    !
    !> \brief Returns the number of photolysis reactions TUV-x reported.
    !> \details Zero when TUV-x is disabled; otherwise equals n_photo_rates.
    !------------------------------------------------------------------------
    integer function tuvx_n_rates()

#ifdef MPAS_USE_MUSICA
        tuvx_n_rates = n_photo_rates
#else
        tuvx_n_rates = 0
#endif

    end function tuvx_n_rates


    logical function tuvx_uses_spatial_climatology()
#ifdef MPAS_USE_MUSICA
        tuvx_uses_spatial_climatology = associated(tuvx_solver) .and. &
                                        trim(upper_column_mode) == 'spatial_climatology'
#else
        tuvx_uses_spatial_climatology = .false.
#endif
    end function tuvx_uses_spatial_climatology


    integer function tuvx_n_extension_layers()
#ifdef MPAS_USE_MUSICA
        tuvx_n_extension_layers = n_ext_layers
#else
        tuvx_n_extension_layers = 0
#endif
    end function tuvx_n_extension_layers


    !------------------------------------------------------------------------
    !  subroutine tuvx_rate_names
    !
    !> \brief Return the cached list of TUV-x photolysis reaction names
    !> \details Caller-allocated character array must be at least
    !>  tuvx_n_rates() elements long. Fills entries 1..n_photo_rates.
    !------------------------------------------------------------------------
    subroutine tuvx_rate_names(names_out)

        character(len=*), intent(out) :: names_out(:)

#ifdef MPAS_USE_MUSICA
        integer :: i
        do i = 1, n_photo_rates
            names_out(i) = photo_names(i)
        end do
#endif

    end subroutine tuvx_rate_names


    !------------------------------------------------------------------------
    !  subroutine tuvx_query_photo_rates
    !
    !> \brief Query a TUV-x config for photolysis reaction names.
    !> \details
    !>  Used during MPAS pool setup, before the real TUV-x solver is
    !>  initialized, so runtime diagnostic fields can be injected from the
    !>  configured rate list.
    !------------------------------------------------------------------------
    subroutine tuvx_query_photo_rates(config_file, n_rates, rate_names, error_code, error_message)

        character(len=*), intent(in)                :: config_file
        integer, intent(out)                        :: n_rates
        character(len=64), allocatable, intent(out)  :: rate_names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

#ifdef MPAS_USE_MUSICA
        integer, parameter :: QUERY_N_LAYERS = 1
        type(grid_t),       pointer :: query_height_grid => null()
        type(grid_t),       pointer :: query_wavelength_grid => null()
        type(profile_t),    pointer :: query_air_profile => null()
        type(profile_t),    pointer :: query_temp_profile => null()
        type(profile_t),    pointer :: query_o3_profile => null()
        type(profile_t),    pointer :: query_o2_profile => null()
        type(radiator_t),   pointer :: query_cloud_radiator => null()
        type(grid_map_t),     pointer :: query_grids => null()
        type(profile_map_t),  pointer :: query_profiles => null()
        type(radiator_map_t), pointer :: query_radiators => null()
        type(tuvx_t),         pointer :: query_solver => null()
        type(mappings_t),     pointer :: photo_ordering => null()
        type(error_t) :: err
        integer :: i
#endif

        error_code = 0
        error_message = ''
        n_rates = 0

        if (len_trim(config_file) == 0) then
            allocate(rate_names(0))
            return
        end if

#ifdef MPAS_USE_MUSICA
        query_height_grid => grid_t("height", "km", QUERY_N_LAYERS, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        query_grids => grid_map_t(err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        call query_grids%add(query_height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        query_air_profile  => profile_t("air",         "molecule cm-3", query_height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        query_temp_profile => profile_t("temperature", "K",             query_height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        query_o3_profile   => profile_t("O3",          "molecule cm-3", query_height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        query_o2_profile   => profile_t("O2",          "molecule cm-3", query_height_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        query_profiles => profile_map_t(err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        call query_profiles%add(query_air_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        call query_profiles%add(query_temp_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        call query_profiles%add(query_o3_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        call query_profiles%add(query_o2_profile, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        query_wavelength_grid => grid_t("wavelength", "nm", N_CLOUD_WAVELENGTH_BINS, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        query_cloud_radiator => radiator_t("clouds", query_height_grid, query_wavelength_grid, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        query_radiators => radiator_map_t(err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100
        call query_radiators%add(query_cloud_radiator, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        query_solver => tuvx_t(trim(config_file), query_grids, query_profiles, query_radiators, err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        ! (P2-F10 / FABLE R4) No wavelength-grid bin-count validation here: the
        ! wavelength grid is not retrievable from the host grid map (see the
        ! note in tuvx_init); grid_map%get("wavelength","nm") raises MUSICA
        ! Error 105. This query path only needs the photolysis rate names.

        photo_ordering => query_solver%get_photolysis_rate_constants_ordering(err)
        if (tuvx_has_error(err, error_code, error_message)) goto 100

        n_rates = photo_ordering%size()
        allocate(rate_names(n_rates))
        do i = 1, n_rates
            rate_names(i) = trim(photo_ordering%name(i))
        end do

100     continue
        ! Deallocate the query-side scaffolding. Free the solver first so it no
        ! longer references the grids/profiles/radiators, then the maps and the
        ! from-host objects (all caller-owned deep copies; the tuvx_t
        ! constructor copies them, so freeing here does not double-free).
        if (associated(photo_ordering)) deallocate(photo_ordering)
        if (associated(query_solver)) deallocate(query_solver)
        if (associated(query_radiators)) deallocate(query_radiators)
        if (associated(query_cloud_radiator)) deallocate(query_cloud_radiator)
        if (associated(query_wavelength_grid)) deallocate(query_wavelength_grid)
        if (associated(query_profiles)) deallocate(query_profiles)
        if (associated(query_air_profile)) deallocate(query_air_profile)
        if (associated(query_temp_profile)) deallocate(query_temp_profile)
        if (associated(query_o3_profile)) deallocate(query_o3_profile)
        if (associated(query_o2_profile)) deallocate(query_o2_profile)
        if (associated(query_grids)) deallocate(query_grids)
        if (associated(query_height_grid)) deallocate(query_height_grid)

        if (error_code /= 0) then
            n_rates = 0
            if (allocated(rate_names)) deallocate(rate_names)
        end if
#else
        allocate(rate_names(0))
#endif

    end subroutine tuvx_query_photo_rates


    !------------------------------------------------------------------------
    !  function tuvx_is_enabled
    !
    !> \brief Returns .true. if TUV-x solver is initialized
    !------------------------------------------------------------------------
    logical function tuvx_is_enabled()

#ifdef MPAS_USE_MUSICA
        tuvx_is_enabled = associated(tuvx_solver)
#else
        tuvx_is_enabled = .false.
#endif

    end function tuvx_is_enabled


    !------------------------------------------------------------------------
    !  subroutine tuvx_finalize
    !
    !> \brief Clean up TUV-x resources
    !------------------------------------------------------------------------
    subroutine tuvx_finalize(dminfo)

        use mpas_derived_types, only: dm_info
#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
        use mpas_dmpar, only: mpas_dmpar_sum_int, mpas_dmpar_min_real
#endif

        type(dm_info), intent(in), optional :: dminfo

#ifdef MPAS_USE_MUSICA
        integer :: global_clip_count
        real(kind=RKIND) :: global_clip_min

        if (present(dminfo)) then
            call mpas_dmpar_sum_int(dminfo, clip_count_input, global_clip_count)
            if (global_clip_count > 0) then
                call mpas_dmpar_min_real(dminfo, clip_min_input, global_clip_min)
            end if
        else
            global_clip_count = clip_count_input
            global_clip_min = clip_min_input
        end if

        if (global_clip_count > 0) then
            call mpas_log_write('[TUV-x] Clipped negative host qO3/qc/qr inputs over run: ' // &
                                '$i (min $r)', intArgs=[global_clip_count], &
                                realArgs=[global_clip_min])
        end if

        if (allocated(photo_rates))     deallocate(photo_rates)
        if (allocated(heating_rates))   deallocate(heating_rates)
        if (allocated(cloud_od))        deallocate(cloud_od)
        if (allocated(cloud_ssa))       deallocate(cloud_ssa)
        if (allocated(cloud_g))         deallocate(cloud_g)
        if (allocated(ext_z_km_edges))  deallocate(ext_z_km_edges)
        if (allocated(ext_T_K_edges))   deallocate(ext_T_K_edges)
        if (allocated(ext_n_air_edges)) deallocate(ext_n_air_edges)
        if (allocated(ext_n_o3_edges))  deallocate(ext_n_o3_edges)
        if (allocated(photo_names))     deallocate(photo_names)
        if (allocated(tuvx_indices))    deallocate(tuvx_indices)

        if (associated(tuvx_solver)) then
            deallocate(tuvx_solver)
            nullify(tuvx_solver)
        end if

        if (associated(cloud_wavelength_grid)) then
            deallocate(cloud_wavelength_grid)
            nullify(cloud_wavelength_grid)
        end if

        nullify(height_grid)
        nullify(air_profile)
        nullify(temp_profile)
        nullify(o3_profile)
        nullify(o2_profile)
        nullify(cloud_radiator)
        upper_column_mode = 'none'
        clip_count_input = 0
        clip_min_input = huge(1.0_RKIND)

        call mpas_log_write('[TUV-x] Finalized.')
#endif

    end subroutine tuvx_finalize

#ifdef MPAS_USE_MUSICA
    !------------------------------------------------------------------------
    !  subroutine load_extension_csv
    !
    !> \brief Read upper-atmosphere climatology CSV for the TUV-x extension
    !> \details
    !>  Expects a header line followed by edge values, one row per edge,
    !>  ordered bottom-up. Columns:
    !>    z_km, T_K, n_air_molec_cm3, n_O3_molec_cm3
    !>  With N rows after the header, n_ext_layers = N - 1.
    !------------------------------------------------------------------------
    subroutine load_extension_csv(filename, error_code, error_message)

        use mpas_io_units, only: mpas_new_unit, mpas_release_unit
        use mpas_log,      only: mpas_log_write

        character(len=*), intent(in)               :: filename
        integer, intent(out)                       :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: iunit, ios, nrows, k
        character(len=512) :: line
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        call mpas_new_unit(iunit)
        open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            error_code = 1
            error_message = '[TUV-x] Failed to open extension file: ' // trim(filename)
            call mpas_release_unit(iunit)
            return
        end if

        ! Count data rows (skip the header)
        read(iunit, '(A)', iostat=ios) line    ! header
        nrows = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            nrows = nrows + 1
        end do

        if (nrows < 2) then
            error_code = 1
            error_message = '[TUV-x] Extension file must contain at least 2 data rows (1 layer)'
            close(iunit)
            call mpas_release_unit(iunit)
            return
        end if

        n_ext_layers = nrows - 1

        allocate(ext_z_km_edges (nrows))
        allocate(ext_T_K_edges  (nrows))
        allocate(ext_n_air_edges(nrows))
        allocate(ext_n_o3_edges (nrows))

        ! Re-read from the top to pull the numbers
        rewind(iunit)
        read(iunit, '(A)', iostat=ios) line    ! skip header
        do k = 1, nrows
            read(iunit, *, iostat=ios) ext_z_km_edges(k), ext_T_K_edges(k), &
                                       ext_n_air_edges(k), ext_n_o3_edges(k)
            if (ios /= 0) then
                error_code = 1
                write(msg, '(A,I0,A,A)') '[TUV-x] Failed parsing row ', k, &
                                          ' of ', trim(filename)
                error_message = trim(msg)
                close(iunit)
                call mpas_release_unit(iunit)
                return
            end if
        end do

        close(iunit)
        call mpas_release_unit(iunit)

        ! Sanity-check monotonic increasing heights
        do k = 2, nrows
            if (ext_z_km_edges(k) <= ext_z_km_edges(k - 1)) then
                error_code = 1
                error_message = '[TUV-x] Extension z_km values must be strictly increasing'
                return
            end if
        end do

        write(msg, '(A,F6.1,A,F6.1,A,I0,A)') '[TUV-x] Extension: z = ', &
             ext_z_km_edges(1), ' -> ', ext_z_km_edges(nrows), ' km (', &
             n_ext_layers, ' layers)'
        call mpas_log_write(trim(msg))

    end subroutine load_extension_csv


    !------------------------------------------------------------------------
    !  subroutine compute_cloud_optical_depth
    !
    !> \brief Compute cloud optical depth per layer from cloud/rain water
    !> \details
    !>  tau = 3 * LWC * dz / (2 * r_eff * rho_water)
    !>  Cloud droplets (r_eff ~10 um) dominate; rain drops (r_eff ~500 um)
    !>  contribute ~50x less per unit mass.
    !------------------------------------------------------------------------
    subroutine compute_cloud_optical_depth(qc, qr, rho_air, dz_m, nLev, cloud_od_1d)

        real(dk), intent(in)  :: qc(:), qr(:), rho_air(:), dz_m(:)
        integer,  intent(in)  :: nLev
        real(dk), intent(out) :: cloud_od_1d(:)

        real(dk) :: lwc_cloud, lwc_rain
        integer :: k

        do k = 1, nLev
            lwc_cloud = max(qc(k), 0.0_dk) * rho_air(k)  ! [kg/m^3]
            lwc_rain  = max(qr(k), 0.0_dk) * rho_air(k)  ! [kg/m^3]

            cloud_od_1d(k) = (3.0_dk * lwc_cloud * dz_m(k)) / &
                             (2.0_dk * R_EFF_CLOUD * real(CHEM_RHO_H2O, dk)) + &
                             (3.0_dk * lwc_rain  * dz_m(k)) / &
                             (2.0_dk * R_EFF_RAIN  * real(CHEM_RHO_H2O, dk))
        end do

    end subroutine compute_cloud_optical_depth

    !------------------------------------------------------------------------
    !  function tuvx_has_error (private helper)
    !------------------------------------------------------------------------
    logical function tuvx_has_error(err, error_code, error_message)

        type(error_t), intent(in) :: err
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=30) :: code_str

        if (err%is_success()) then
            error_code = 0
            error_message = ''
            tuvx_has_error = .false.
            return
        end if

        error_code = err%code()
        write(code_str, '(I30)') err%code()
        error_message = '[TUV-x Error]: ' // err%category() // '[' // &
                        trim(adjustl(code_str)) // ']: ' // err%message()
        tuvx_has_error = .true.

    end function tuvx_has_error
#endif

end module mpas_tuvx
