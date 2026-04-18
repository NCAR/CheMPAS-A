! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_tuvx
!
!> \brief TUV-x photolysis rate calculations for CheMPAS
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

    implicit none

    private

    public :: tuvx_init, tuvx_compute_photolysis, tuvx_is_enabled, tuvx_finalize

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
    real(dk), parameter :: RHO_WATER   = 1000.0_dk    ! water density [kg/m^3]

    integer, save :: j_no2_index     = -1  ! index in photolysis rate constants
    integer, save :: n_photo_rates   = 0
    integer, save :: n_heating_rates = 0

    ! Column decomposition:
    !   n_mpas_layers   : layers supplied by MPAS host state
    !   n_ext_layers    : climatology layers above MPAS top (0 when disabled)
    !   n_total_layers  : full column TUV-x sees (= n_mpas + n_ext)
    integer, save :: n_mpas_layers   = 0
    integer, save :: n_ext_layers    = 0
    integer, save :: n_total_layers  = 0

    ! Upper-atmosphere climatology edge values (sized n_ext_layers+1 when enabled).
    ! Values are defined at edge heights; midpoints are computed on the fly.
    real(dk), allocatable :: ext_z_km_edges(:)
    real(dk), allocatable :: ext_T_K_edges(:)
    real(dk), allocatable :: ext_n_air_edges(:)
    real(dk), allocatable :: ext_n_o3_edges(:)

    real(dk), allocatable :: photo_rates(:,:)    ! (n_total_layers+1, n_photo_rates)
    real(dk), allocatable :: heating_rates(:,:)  ! (n_total_layers+1, n_heating_rates)

    ! Physical constants
    real(dk), parameter :: NA    = 6.02214076e23_dk  ! Avogadro [mol^-1]
    real(dk), parameter :: M_AIR = 0.02897_dk        ! dry air molar mass [kg/mol]
    real(dk), parameter :: M_O3  = 0.04800_dk        ! ozone molar mass [kg/mol]
    real(dk), parameter :: VMR_O2 = 0.2095_dk        ! O2 volume mixing ratio
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
    subroutine tuvx_init(config_file, nVertLevels_in, enable_extension, &
                         extension_file, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif

        character(len=*), intent(in)               :: config_file
        integer, intent(in)                        :: nVertLevels_in
        logical, intent(in)                        :: enable_extension
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

        error_code = 0
        error_message = ''

        if (len_trim(config_file) == 0) then
            call mpas_log_write('[TUV-x] No config file specified; TUV-x disabled (Phase 1 fallback).')
            return
        end if

        n_mpas_layers = nVertLevels_in
        n_ext_layers  = 0

        call mpas_log_write('[TUV-x] Initializing with config: ' // trim(config_file))
        write(msg, '(A,I0)') '[TUV-x] MPAS nVertLevels = ', n_mpas_layers
        call mpas_log_write(trim(msg))

        ! Load upper-atmosphere climatology if the extension is enabled
        if (enable_extension .and. len_trim(extension_file) > 0) then
            call load_extension_csv(trim(extension_file), error_code, error_message)
            if (error_code /= 0) return
            write(msg, '(A,I0,A,A)') '[TUV-x] Extension enabled: n_ext_layers = ', n_ext_layers, &
                                     ' from ', trim(extension_file)
            call mpas_log_write(trim(msg))
        else if (enable_extension) then
            error_code = 1
            error_message = '[TUV-x] config_tuvx_top_extension is .true. but config_tuvx_extension_file is empty'
            return
        end if

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
        n_wavelength_bins = 102  ! CAM wavelength grid: 103 edges - 1

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

        ! Query photolysis ordering to find j_NO2 index
        photo_ordering => tuvx_solver%get_photolysis_rate_constants_ordering(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        n_photo_rates = photo_ordering%size()

        j_no2_index = -1
        do i = 1, n_photo_rates
            call mpas_log_write('[TUV-x] Photolysis reaction ' // trim(photo_ordering%name(i)))
            if (trim(photo_ordering%name(i)) == 'jNO2') then
                j_no2_index = photo_ordering%index(i)
            end if
        end do
        deallocate(photo_ordering)

        if (j_no2_index < 1) then
            error_code = 1
            error_message = '[TUV-x] jNO2 reaction not found in photolysis ordering'
            return
        end if

        write(msg, '(A,I0)') '[TUV-x] j_NO2 index in photo_rates = ', j_no2_index
        call mpas_log_write(trim(msg))

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
    !> \brief Compute j_NO2 for a single atmospheric column
    !> \details
    !>  Converts MPAS fields to TUV-x units, updates from-host profiles,
    !>  calls TUV-x run, and extracts j_NO2 (edge values averaged to
    !>  midpoints).
    !------------------------------------------------------------------------
    subroutine tuvx_compute_photolysis(cos_sza, z_edges_m, rho_air, temperature, &
                                        qo3, qc, qr, nVertLevels_in, &
                                        j_no2_column, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif

        real(kind=RKIND), intent(in)  :: cos_sza
        real(kind=RKIND), intent(in)  :: z_edges_m(:)       ! nVertLevels+1 edges [m]
        real(kind=RKIND), intent(in)  :: rho_air(:)          ! nVertLevels midpoints [kg/m³]
        real(kind=RKIND), intent(in)  :: temperature(:)      ! nVertLevels midpoints [K]
        real(kind=RKIND), intent(in)  :: qo3(:)              ! nVertLevels [kg/kg]
        real(kind=RKIND), intent(in)  :: qc(:)               ! nVertLevels cloud water [kg/kg]
        real(kind=RKIND), intent(in)  :: qr(:)               ! nVertLevels rain water [kg/kg]
        integer, intent(in)           :: nVertLevels_in
        real(kind=RKIND), intent(out) :: j_no2_column(:)     ! nVertLevels
        integer, intent(out)          :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        real(dk) :: z_km(n_total_layers + 1), z_mid_km(n_total_layers)
        real(dk) :: air_nd(n_total_layers), air_edges(n_total_layers + 1)
        real(dk) :: temp_mid(n_total_layers), temp_edges(n_total_layers + 1)
        real(dk) :: o3_nd(n_total_layers), o2_nd(n_total_layers)
        real(dk) :: air_ld(n_total_layers), o3_ld(n_total_layers), o2_ld(n_total_layers)
        real(dk) :: dz_cm(n_total_layers), dz_m_mpas(n_mpas_layers)
        real(dk) :: cloud_tau_1d(n_mpas_layers)
        real(dk) :: sza_rad, esd
        type(error_t) :: err
        integer :: k, j, iWave

        error_code = 0
        error_message = ''

        ! Nighttime: skip TUV-x call, return zero
        if (cos_sza <= 0.0_RKIND) then
            j_no2_column(1:n_mpas_layers) = 0.0_RKIND
            return
        end if

        ! ---- Build the composite column (MPAS below + climatology above) ----

        ! Height edges [km]: MPAS edges 1..n_mpas+1, then climatology edges
        ! 2..n_ext+1 (skipping ext edge 1 which coincides with MPAS top).
        z_km(1:n_mpas_layers + 1) = real(z_edges_m(1:n_mpas_layers + 1), dk) * 1.0e-3_dk
        do j = 1, n_ext_layers
            z_km(n_mpas_layers + 1 + j) = ext_z_km_edges(j + 1)
        end do
        z_mid_km = 0.5_dk * (z_km(1:n_total_layers) + z_km(2:n_total_layers + 1))

        ! Layer thickness in cm (for layer densities)
        dz_cm = (z_km(2:n_total_layers + 1) - z_km(1:n_total_layers)) * 1.0e5_dk

        ! MPAS-layer midpoint values from host state
        air_nd (1:n_mpas_layers) = real(rho_air(1:n_mpas_layers), dk) * (NA / M_AIR) * 1.0e-6_dk
        temp_mid(1:n_mpas_layers) = real(temperature(1:n_mpas_layers), dk)
        o3_nd  (1:n_mpas_layers) = real(qo3(1:n_mpas_layers), dk) &
                                   * real(rho_air(1:n_mpas_layers), dk) &
                                   * (NA / M_O3) * 1.0e-6_dk

        ! Extension-layer midpoint values from CSV (average adjacent edges)
        do j = 1, n_ext_layers
            air_nd (n_mpas_layers + j) = 0.5_dk * (ext_n_air_edges(j) + ext_n_air_edges(j + 1))
            temp_mid(n_mpas_layers + j) = 0.5_dk * (ext_T_K_edges(j)   + ext_T_K_edges(j + 1))
            o3_nd  (n_mpas_layers + j) = 0.5_dk * (ext_n_o3_edges(j)  + ext_n_o3_edges(j + 1))
        end do

        ! O2 = 0.2095 × air [molecule/cm³] across the full column
        o2_nd = VMR_O2 * air_nd

        ! Edge values (simple extrapolation from adjacent midpoints).
        ! At the MPAS/extension join (edge n_mpas+1), blend the top MPAS midpoint
        ! with the first extension midpoint instead of the plain extrapolation.
        air_edges (1) = air_nd (1)
        temp_edges(1) = temp_mid(1)
        do k = 2, n_total_layers
            air_edges (k) = 0.5_dk * (air_nd (k - 1) + air_nd (k))
            temp_edges(k) = 0.5_dk * (temp_mid(k - 1) + temp_mid(k))
        end do
        air_edges (n_total_layers + 1) = air_nd (n_total_layers)
        temp_edges(n_total_layers + 1) = temp_mid(n_total_layers)

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

        call temp_profile%set_midpoint_values(temp_mid, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call temp_profile%set_edge_values(temp_edges, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        call o3_profile%set_midpoint_values(o3_nd, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call o3_profile%set_layer_densities(o3_ld, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        call o2_profile%set_midpoint_values(o2_nd, err)
        if (tuvx_has_error(err, error_code, error_message)) return
        call o2_profile%set_layer_densities(o2_ld, err)
        if (tuvx_has_error(err, error_code, error_message)) return

        ! Compute cloud optical depth per layer and update cloud radiator.
        ! Clouds only live in the MPAS slice; extension rows of cloud_od were
        ! zero-filled at init and are never touched here.
        call compute_cloud_optical_depth( &
            real(qc(1:n_mpas_layers), dk), real(qr(1:n_mpas_layers), dk), &
            real(rho_air(1:n_mpas_layers), dk), dz_m_mpas, n_mpas_layers, cloud_tau_1d)

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

        ! Extract j_NO2 for the MPAS slice only; extension-layer values are
        ! computed by TUV-x but discarded (the host has no storage for them).
        do k = 1, n_mpas_layers
            j_no2_column(k) = real(0.5_dk * (photo_rates(k, j_no2_index) + &
                                              photo_rates(k + 1, j_no2_index)), RKIND)
        end do
#else
        error_code = 0
        error_message = ''
        j_no2_column = 0.0_RKIND
#endif

    end subroutine tuvx_compute_photolysis


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
    subroutine tuvx_finalize()

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write

        if (allocated(photo_rates))     deallocate(photo_rates)
        if (allocated(heating_rates))   deallocate(heating_rates)
        if (allocated(cloud_od))        deallocate(cloud_od)
        if (allocated(cloud_ssa))       deallocate(cloud_ssa)
        if (allocated(cloud_g))         deallocate(cloud_g)
        if (allocated(ext_z_km_edges))  deallocate(ext_z_km_edges)
        if (allocated(ext_T_K_edges))   deallocate(ext_T_K_edges)
        if (allocated(ext_n_air_edges)) deallocate(ext_n_air_edges)
        if (allocated(ext_n_o3_edges))  deallocate(ext_n_o3_edges)

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
                             (2.0_dk * R_EFF_CLOUD * RHO_WATER) + &
                             (3.0_dk * lwc_rain  * dz_m(k)) / &
                             (2.0_dk * R_EFF_RAIN  * RHO_WATER)
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
