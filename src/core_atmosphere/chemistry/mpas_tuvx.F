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
    type(profile_t),    pointer :: air_profile   => null()
    type(profile_t),    pointer :: temp_profile  => null()
    type(profile_t),    pointer :: o3_profile    => null()
    type(profile_t),    pointer :: o2_profile    => null()

    ! Cloud radiator (Phase 3)
    type(radiator_t),   pointer :: cloud_radiator => null()
    integer, save :: n_wavelength_bins = 0

    real(dk), allocatable :: cloud_od(:,:)   ! (n_layers, n_wavelength_bins)
    real(dk), allocatable :: cloud_ssa(:,:)  ! (n_layers, n_wavelength_bins)
    real(dk), allocatable :: cloud_g(:,:,:)  ! (n_layers, n_wavelength_bins, 1)

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
    integer, save :: n_layers        = 0   ! nVertLevels (= number of grid sections)

    real(dk), allocatable :: photo_rates(:,:)    ! (n_layers+1, n_photo_rates)
    real(dk), allocatable :: heating_rates(:,:)  ! (n_layers+1, n_heating_rates)

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
    subroutine tuvx_init(config_file, nVertLevels_in, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif

        character(len=*), intent(in)               :: config_file
        integer, intent(in)                        :: nVertLevels_in
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

        n_layers = nVertLevels_in

        call mpas_log_write('[TUV-x] Initializing with config: ' // trim(config_file))
        write(msg, '(A,I0)') '[TUV-x] nVertLevels (grid sections) = ', n_layers
        call mpas_log_write(trim(msg))

        ! Create from-host height grid
        height_grid => grid_t("height", "km", n_layers, err)
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
        ! We create a standalone from-host wavelength grid here, used only
        ! to dimension the radiator arrays. The solver's actual wavelength
        ! grid is loaded from the JSON config's CSV file.
        block
            type(grid_t), pointer :: wl_grid => null()

            n_wavelength_bins = 102  ! CAM wavelength grid: 103 edges - 1

            wl_grid => grid_t("wavelength", "nm", n_wavelength_bins, err)
            if (tuvx_has_error(err, error_code, error_message)) return

            cloud_radiator => radiator_t("clouds", height_grid, wl_grid, err)
            if (tuvx_has_error(err, error_code, error_message)) return

            write(msg, '(A,I0)') '[TUV-x] Cloud radiator created with n_wavelength_bins = ', n_wavelength_bins
            call mpas_log_write(trim(msg))

            ! wl_grid is not added to the grid map — it's only used for sizing
            deallocate(wl_grid)
        end block

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
        allocate(cloud_od(n_layers, n_wavelength_bins))
        allocate(cloud_ssa(n_layers, n_wavelength_bins))
        allocate(cloud_g(n_layers, n_wavelength_bins, 1))

        cloud_od  = 0.0_dk
        cloud_ssa = CLOUD_SSA_VAL
        cloud_g(:,:,1) = CLOUD_G_VAL

        ! Query photolysis ordering to find j_NO2 index
        photo_ordering => tuvx_solver%get_photolysis_rate_constants_ordering(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        n_photo_rates = photo_ordering%size()

        j_no2_index = -1
        do i = 1, n_photo_rates
            call mpas_log_write('[TUV-x] Photolysis reaction ' // trim(photo_ordering%name(i)))
            if (trim(photo_ordering%name(i)) == 'jno2') then
                j_no2_index = photo_ordering%index(i)
            end if
        end do
        deallocate(photo_ordering)

        if (j_no2_index < 1) then
            error_code = 1
            error_message = '[TUV-x] jno2 reaction not found in photolysis ordering'
            return
        end if

        write(msg, '(A,I0)') '[TUV-x] j_NO2 index in photo_rates = ', j_no2_index
        call mpas_log_write(trim(msg))

        ! Query heating rates ordering
        heat_ordering => tuvx_solver%get_heating_rates_ordering(err)
        if (tuvx_has_error(err, error_code, error_message)) return
        n_heating_rates = heat_ordering%size()
        deallocate(heat_ordering)

        ! Allocate persistent work arrays (n_layers+1 edges for TUV-x output)
        allocate(photo_rates(n_layers + 1, max(n_photo_rates, 1)))
        allocate(heating_rates(n_layers + 1, max(n_heating_rates, 1)))
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
        real(dk) :: z_km(nVertLevels_in + 1), z_mid_km(nVertLevels_in)
        real(dk) :: air_nd(nVertLevels_in), air_edges(nVertLevels_in + 1)
        real(dk) :: temp_mid(nVertLevels_in), temp_edges(nVertLevels_in + 1)
        real(dk) :: o3_nd(nVertLevels_in), o2_nd(nVertLevels_in)
        real(dk) :: air_ld(nVertLevels_in), o3_ld(nVertLevels_in), o2_ld(nVertLevels_in)
        real(dk) :: dz_cm(nVertLevels_in), dz_m(nVertLevels_in)
        real(dk) :: cloud_tau_1d(nVertLevels_in)
        real(dk) :: sza_rad, esd
        type(error_t) :: err
        integer :: k, iWave

        error_code = 0
        error_message = ''

        ! Nighttime: skip TUV-x call, return zero
        if (cos_sza <= 0.0_RKIND) then
            j_no2_column(1:nVertLevels_in) = 0.0_RKIND
            return
        end if

        ! Convert MPAS height edges [m] to [km]
        z_km = real(z_edges_m(1:nVertLevels_in + 1), dk) * 1.0e-3_dk
        z_mid_km = 0.5_dk * (z_km(1:nVertLevels_in) + z_km(2:nVertLevels_in + 1))

        ! Layer thickness in cm (for layer densities)
        dz_cm = (z_km(2:nVertLevels_in + 1) - z_km(1:nVertLevels_in)) * 1.0e5_dk

        ! Air number density: rho [kg/m³] -> [molecule/cm³]
        air_nd = real(rho_air(1:nVertLevels_in), dk) * (NA / M_AIR) * 1.0e-6_dk

        ! Temperature midpoints [K]
        temp_mid = real(temperature(1:nVertLevels_in), dk)

        ! O3 number density: qO3 [kg/kg] * rho [kg/m³] -> [molecule/cm³]
        o3_nd = real(qo3(1:nVertLevels_in), dk) * real(rho_air(1:nVertLevels_in), dk) &
                * (NA / M_O3) * 1.0e-6_dk

        ! O2 = 0.2095 × air [molecule/cm³]
        o2_nd = VMR_O2 * air_nd

        ! Edge values (simple extrapolation from adjacent midpoints)
        air_edges(1) = air_nd(1)
        air_edges(nVertLevels_in + 1) = air_nd(nVertLevels_in)
        do k = 2, nVertLevels_in
            air_edges(k) = 0.5_dk * (air_nd(k-1) + air_nd(k))
        end do

        temp_edges(1) = temp_mid(1)
        temp_edges(nVertLevels_in + 1) = temp_mid(nVertLevels_in)
        do k = 2, nVertLevels_in
            temp_edges(k) = 0.5_dk * (temp_mid(k-1) + temp_mid(k))
        end do

        ! Layer thickness in metres (for cloud OD)
        dz_m = z_km(2:nVertLevels_in + 1) - z_km(1:nVertLevels_in)
        dz_m = dz_m * 1.0e3_dk  ! km -> m

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

        ! Compute cloud optical depth per layer and update cloud radiator
        call compute_cloud_optical_depth( &
            real(qc(1:nVertLevels_in), dk), real(qr(1:nVertLevels_in), dk), &
            real(rho_air(1:nVertLevels_in), dk), dz_m, nVertLevels_in, cloud_tau_1d)

        ! Broadcast 1-D cloud OD to all wavelength bins
        do iWave = 1, n_wavelength_bins
            cloud_od(:, iWave) = cloud_tau_1d(:)
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

        ! Extract j_NO2: average edge values to midpoints
        do k = 1, nVertLevels_in
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

        if (allocated(photo_rates))   deallocate(photo_rates)
        if (allocated(heating_rates)) deallocate(heating_rates)
        if (allocated(cloud_od))      deallocate(cloud_od)
        if (allocated(cloud_ssa))     deallocate(cloud_ssa)
        if (allocated(cloud_g))       deallocate(cloud_g)

        if (associated(tuvx_solver)) then
            deallocate(tuvx_solver)
            nullify(tuvx_solver)
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
