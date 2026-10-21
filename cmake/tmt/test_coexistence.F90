! Exercise TMT and MUSICA ownership in the same process, including MIEM.
program test_coexistence
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use musica_util, only: miem_error_t => error_t, miem_string_t => string_t
  use musica, only: get_musica_version
  use musica_emissions, only: mechanism_t, emissions_t
  use tmt_micm, only: micm_t, RosenbrockStandardOrder, solver_stats_t
  use tmt_state, only: state_t
  use tmt_tuvx, only: tuvx_t, get_tuvx_version
  use tmt_tuvx_grid, only: grid_t
  use tmt_tuvx_grid_map, only: grid_map_t
  use tmt_tuvx_profile_map, only: profile_map_t
  use tmt_tuvx_radiator_map, only: radiator_map_t
  use tmt_util, only: error_t, string_t, mappings_t
  implicit none
  type(micm_t), pointer :: solver
  type(state_t), pointer :: state
  type(grid_t), pointer :: grid
  type(tuvx_t), pointer :: photolysis
  type(grid_map_t), pointer :: grids, model_grids
  type(profile_map_t), pointer :: profiles
  type(radiator_map_t), pointer :: radiators
  type(mappings_t), pointer :: photo_names, heating_names
  type(mechanism_t), pointer :: mechanism
  type(emissions_t), pointer :: emissions
  type(error_t) :: error
  type(miem_error_t) :: miem_error
  type(string_t) :: version, result
  type(miem_string_t) :: miem_version
  type(solver_stats_t) :: stats
  character(4096) :: path, photo_path
  real(real64), allocatable :: rates(:,:), heating(:,:)
  integer :: a, iteration, levels

  call get_command_argument(1, path)
  call get_command_argument(2, photo_path)
  do iteration = 1, 3
    version = get_tuvx_version()
    miem_version = get_musica_version()
    if (index(version%value_, 'TurboTUV') == 0) stop 1
    if (index(miem_version%value_, 'TMT') /= 0) stop 2

    ! Repeated failures exercise each library's error allocation/finalization.
    mechanism => mechanism_t('nonexistent-coexistence-config.yaml', miem_error)
    if (miem_error%is_success()) stop 3
    if (associated(mechanism)) deallocate(mechanism)
    solver => micm_t('nonexistent-coexistence-config.json', RosenbrockStandardOrder, error)
    if (error%is_success()) stop 4
    if (associated(solver)) deallocate(solver)

    mechanism => mechanism_t(trim(path), miem_error)
    if (.not. miem_error%is_success()) stop 5
    emissions => emissions_t(mechanism, 1, 1, miem_error)
    if (.not. miem_error%is_success()) stop 6
    call emissions%run(0.0_real64, 1.0_real64, miem_error)
    if (.not. miem_error%is_success()) stop 7

    grid => grid_t('height', 'km', 1, error)
    if (.not. error%is_success()) stop 8
    call grid%set_edges([0.0_real64, 1.0_real64], error)
    if (.not. error%is_success()) stop 9

    solver => micm_t(trim(path), RosenbrockStandardOrder, error)
    if (.not. error%is_success()) then
      result = error%message()
      print *, result%value_
      stop 10
    end if
    state => solver%get_state(1, error)
    if (.not. error%is_success()) stop 11
    state%conditions(1)%temperature = 300.0_real64
    state%conditions(1)%pressure = 100000.0_real64
    state%conditions(1)%air_density = 40.0_real64
    a = state%species_ordering%index('A', error)
    if (.not. error%is_success()) stop 12
    state%concentrations = 0.0_real64
    state%concentrations(a) = 1.0_real64
    call solver%solve(1.0_real64, state, result, stats, error)
    if (.not. error%is_success()) stop 13
    if (.not. all(ieee_is_finite(state%concentrations))) stop 14
    if (abs(state%concentrations(a) - exp(-0.2_real64)) > 1.0e-6_real64) stop 14
    deallocate(state, solver, grid, emissions, mechanism)
  end do

  grids => grid_map_t(error)
  if (.not. error%is_success()) stop 15
  profiles => profile_map_t(error)
  if (.not. error%is_success()) stop 16
  radiators => radiator_map_t(error)
  if (.not. error%is_success()) stop 17
  photolysis => tuvx_t(trim(photo_path), grids, profiles, radiators, error)
  if (.not. error%is_success()) then
    result = error%message()
    print *, result%value_
    stop 18
  end if
  model_grids => photolysis%get_grids(error)
  if (.not. error%is_success()) stop 19
  grid => model_grids%get('height', 'km', error)
  if (.not. error%is_success()) stop 20
  levels = grid%number_of_sections(error) + 1
  if (.not. error%is_success()) stop 21
  photo_names => photolysis%get_photolysis_rate_constants_ordering(error)
  if (.not. error%is_success()) stop 22
  heating_names => photolysis%get_heating_rates_ordering(error)
  if (.not. error%is_success()) stop 23
  allocate(rates(levels, photo_names%size()), heating(levels, max(1, heating_names%size())))
  call photolysis%run(acos(-1.0_real64) / 6.0_real64, 1.0_real64, rates, heating, error)
  if (.not. error%is_success()) stop 24
  if (.not. all(ieee_is_finite(rates))) stop 25
  if (any(rates < 0.0_real64) .or. .not. any(rates > 0.0_real64)) stop 26
  deallocate(photolysis, grids, model_grids, profiles, radiators, grid, photo_names, heating_names)
  print *, 'TMT chemistry/photolysis and MUSICA MIEM coexist successfully'
end program test_coexistence
