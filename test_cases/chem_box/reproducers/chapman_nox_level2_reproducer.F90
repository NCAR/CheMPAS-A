program chapman_nox_level2_reproducer
  use iso_fortran_env, only: real64
  use musica_util, only: dk => musica_dk, error_t, string_t
  use musica_micm, only: micm_t, solver_stats_t, RosenbrockStandardOrder, &
      BackwardEulerStandardOrder
  use musica_state, only: state_t

  implicit none

  integer, parameter :: NUM_CELLS = 1
  real(dk), parameter :: DEFAULT_DT = 3.0_dk
  real(real64), parameter :: DEFAULT_REL_TOL = 1.0e-15_real64

  type(micm_t), pointer :: micm
  type(state_t), pointer :: state
  type(error_t) :: error
  type(string_t) :: solver_state
  type(solver_stats_t) :: solver_stats

  character(len=512) :: yaml_path
  character(len=64) :: arg
  character(len=64) :: solver_name
  character(len=64) :: state_mode
  integer :: solver_type
  integer :: ios
  integer :: var_stride, rp_var_stride
  integer :: idx_o2, idx_o, idx_o1d, idx_o3, idx_no, idx_no2
  integer :: rp_jo2, rp_jo3_o, rp_jo3_o1d, rp_jno2
  integer :: idx
  real(dk) :: dt
  real(real64) :: rel_tol
  real(real64) :: elapsed, remaining, advanced
  real(real64) :: target_time, time_epsilon

  call parse_arguments()

  micm => micm_t(trim(yaml_path), solver_type, error)
  call check(error, "create MICM")

  state => micm%get_state(NUM_CELLS, error)
  call check(error, "get state")
  call state%set_relative_tolerance(rel_tol, error)
  call check(error, "set relative tolerance")

  call lookup_indices()

  var_stride = state%species_strides%variable
  rp_var_stride = state%rate_parameters_strides%variable

  state%concentrations = 0.0_dk
  state%rate_parameters = 0.0_dk

  state%conditions(1)%temperature = real(296.542333817352528_real64, dk)
  state%conditions(1)%pressure = real(94417.3244009938644_real64, dk)
  state%conditions(1)%air_density = real(37.9828906716704537_real64, dk)

  select case (trim(state_mode))
  case ("chempas_raw")
    call set_chempas_raw_state()
  case ("name_species_raw_rates")
    call set_name_species_raw_rates_state()
  case ("name_mapped")
    call set_name_mapped_state()
  case default
    write(*,'(A,A)') "invalid state_mode: ", trim(state_mode)
    stop 2
  end select

  write(*,'(A)') "chapman_nox level-2 one-cell reproducer"
  write(*,'(A,A)') "  config = ", trim(yaml_path)
  write(*,'(A,A)') "  solver = ", trim(solver_name)
  write(*,'(A,A)') "  state_mode = ", trim(state_mode)
  write(*,'(A,ES24.16E3)') "  dt_s = ", real(dt, real64)
  write(*,'(A,ES24.16E3)') "  relative_tolerance = ", rel_tol
  write(*,'(A,ES24.16E3)') "  temperature_K = ", real(state%conditions(1)%temperature, real64)
  write(*,'(A,ES24.16E3)') "  pressure_Pa = ", real(state%conditions(1)%pressure, real64)
  write(*,'(A,ES24.16E3)') "  air_density_mol_m3 = ", real(state%conditions(1)%air_density, real64)
  call print_ordering()
  call print_rate_parameters()

  call print_state("before")

  target_time = real(dt, real64)
  time_epsilon = max(1.0e-7_real64 * target_time, 1.0e-10_real64)

  elapsed = 0.0_real64
  do while (elapsed < target_time - time_epsilon)
    remaining = target_time - elapsed
    call micm%solve(remaining, state, solver_state, solver_stats, error)
    call check(error, "solve")
    advanced = solver_stats%final_time()
    if (advanced <= 0.0_real64) then
      write(*,*) "ERROR: solver advanced <= 0"
      stop 2
    end if
    elapsed = elapsed + advanced
    write(*,'(A,A,A,ES24.16E3,A,I0,A,I0,A,I0,A,I0)') &
        "solver_state=", solver_state%get_char_array(), &
        " advanced_s=", advanced, &
        " accepted=", solver_stats%accepted(), &
        " rejected=", solver_stats%rejected(), &
        " steps=", solver_stats%number_of_steps(), &
        " function_calls=", solver_stats%function_calls()
  end do

  call print_state("after")

contains

  subroutine parse_arguments()
    call get_command_argument(1, yaml_path)
    if (len_trim(yaml_path) == 0) then
      write(*,'(A)') "usage: chapman_nox_level2_reproducer <chapman_nox.yaml> [dt_s] [solver] [relative_tolerance] [state_mode]"
      write(*,'(A)') "       solver: rosenbrock or backward_euler"
      write(*,'(A)') "       state_mode: chempas_raw, name_species_raw_rates, or name_mapped"
      stop 2
    end if

    dt = DEFAULT_DT
    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) then
      read(arg, *, iostat=ios) dt
      if (ios /= 0 .or. dt <= 0.0_dk) then
        write(*,'(A,A)') "invalid dt_s: ", trim(arg)
        stop 2
      end if
    end if

    solver_name = "rosenbrock"
    solver_type = RosenbrockStandardOrder
    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) then
      select case (trim(arg))
      case ("rosenbrock", "Rosenbrock", "ROSENBROCK")
        solver_name = "rosenbrock"
        solver_type = RosenbrockStandardOrder
      case ("backward_euler", "BackwardEuler", "BACKWARD_EULER", "be", "BE")
        solver_name = "backward_euler"
        solver_type = BackwardEulerStandardOrder
      case default
        write(*,'(A,A)') "invalid solver: ", trim(arg)
        stop 2
      end select
    end if

    rel_tol = DEFAULT_REL_TOL
    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) then
      read(arg, *, iostat=ios) rel_tol
      if (ios /= 0 .or. rel_tol <= 0.0_real64) then
        write(*,'(A,A)') "invalid relative_tolerance: ", trim(arg)
        stop 2
      end if
    end if

    state_mode = "chempas_raw"
    call get_command_argument(5, arg)
    if (len_trim(arg) > 0) state_mode = trim(arg)
  end subroutine parse_arguments

  subroutine set_chempas_raw_state()
    call set_conc_raw(1, 0.0_real64)
    call set_conc_raw(2, 7.95587740008720612_real64)
    call set_conc_raw(3, 1.03579990236078073e-6_real64)
    call set_conc_raw(4, 4.06074210565908480e-22_real64)
    call set_conc_raw(5, 5.69633226736554358e-10_real64)
    call set_conc_raw(6, 1.32914419571862701e-9_real64)

    call set_rate_raw(1, 6.91433685646463889e-3_real64)
    call set_rate_raw(2, 3.54345945700791284e-37_real64)
    call set_rate_raw(3, 3.66372021584529420e-4_real64)
    call set_rate_raw(4, 9.25367724139163325e-6_real64)
  end subroutine set_chempas_raw_state

  subroutine set_name_mapped_state()
    call set_conc(idx_o,   0.0_real64)
    call set_conc(idx_o2,  7.95587740008720612_real64)
    call set_conc(idx_o3,  1.03579990236078073e-6_real64)
    call set_conc(idx_o1d, 4.06074210565908480e-22_real64)
    call set_conc(idx_no,  5.69633226736554358e-10_real64)
    call set_conc(idx_no2, 1.32914419571862701e-9_real64)

    call set_rate(rp_jno2,    6.91433685646463889e-3_real64)
    call set_rate(rp_jo2,     3.54345945700791284e-37_real64)
    call set_rate(rp_jo3_o,   3.66372021584529420e-4_real64)
    call set_rate(rp_jo3_o1d, 9.25367724139163325e-6_real64)
  end subroutine set_name_mapped_state

  subroutine set_name_species_raw_rates_state()
    call set_conc(idx_o,   0.0_real64)
    call set_conc(idx_o2,  7.95587740008720612_real64)
    call set_conc(idx_o3,  1.03579990236078073e-6_real64)
    call set_conc(idx_o1d, 4.06074210565908480e-22_real64)
    call set_conc(idx_no,  5.69633226736554358e-10_real64)
    call set_conc(idx_no2, 1.32914419571862701e-9_real64)

    call set_rate_raw(1, 6.91433685646463889e-3_real64)
    call set_rate_raw(2, 3.54345945700791284e-37_real64)
    call set_rate_raw(3, 3.66372021584529420e-4_real64)
    call set_rate_raw(4, 9.25367724139163325e-6_real64)
  end subroutine set_name_species_raw_rates_state

  subroutine lookup_indices()
    idx_o2  = state%species_ordering%index("O2", error)
    call check(error, "index O2")
    idx_o   = state%species_ordering%index("O", error)
    call check(error, "index O")
    idx_o1d = state%species_ordering%index("O1D", error)
    call check(error, "index O1D")
    idx_o3  = state%species_ordering%index("O3", error)
    call check(error, "index O3")
    idx_no  = state%species_ordering%index("NO", error)
    call check(error, "index NO")
    idx_no2 = state%species_ordering%index("NO2", error)
    call check(error, "index NO2")

    rp_jo2     = state%rate_parameters_ordering%index("PHOTO.jO2", error)
    call check(error, "index PHOTO.jO2")
    rp_jo3_o   = state%rate_parameters_ordering%index("PHOTO.jO3_O", error)
    call check(error, "index PHOTO.jO3_O")
    rp_jo3_o1d = state%rate_parameters_ordering%index("PHOTO.jO3_O1D", error)
    call check(error, "index PHOTO.jO3_O1D")
    rp_jno2    = state%rate_parameters_ordering%index("PHOTO.jNO2", error)
    call check(error, "index PHOTO.jNO2")
  end subroutine lookup_indices

  subroutine check(err, context)
    type(error_t), intent(in) :: err
    character(len=*), intent(in) :: context
    if (.not. err%is_success()) then
      write(*,*) "ERROR in ", trim(context), ": ", err%message()
      stop 1
    end if
  end subroutine check

  subroutine set_conc(spec_idx, value)
    integer, intent(in) :: spec_idx
    real(real64), intent(in) :: value
    idx = 1 + (spec_idx - 1) * var_stride
    state%concentrations(idx) = real(value, dk)
  end subroutine set_conc

  subroutine set_conc_raw(slot_idx, value)
    integer, intent(in) :: slot_idx
    real(real64), intent(in) :: value
    state%concentrations(slot_idx) = real(value, dk)
  end subroutine set_conc_raw

  function get_conc(spec_idx) result(value)
    integer, intent(in) :: spec_idx
    real(real64) :: value
    idx = 1 + (spec_idx - 1) * var_stride
    value = real(state%concentrations(idx), real64)
  end function get_conc

  subroutine set_rate(rp_idx, value)
    integer, intent(in) :: rp_idx
    real(real64), intent(in) :: value
    idx = 1 + (rp_idx - 1) * rp_var_stride
    state%rate_parameters(idx) = real(value, dk)
  end subroutine set_rate

  subroutine set_rate_raw(slot_idx, value)
    integer, intent(in) :: slot_idx
    real(real64), intent(in) :: value
    state%rate_parameters(slot_idx) = real(value, dk)
  end subroutine set_rate_raw

  function get_rate(rp_idx) result(value)
    integer, intent(in) :: rp_idx
    real(real64) :: value
    idx = 1 + (rp_idx - 1) * rp_var_stride
    value = real(state%rate_parameters(idx), real64)
  end function get_rate

  subroutine print_ordering()
    write(*,'(A)') "  MICM species indices:"
    write(*,'(A,I0)') "    O1D = ", idx_o1d
    write(*,'(A,I0)') "    O2  = ", idx_o2
    write(*,'(A,I0)') "    O3  = ", idx_o3
    write(*,'(A,I0)') "    O   = ", idx_o
    write(*,'(A,I0)') "    NO  = ", idx_no
    write(*,'(A,I0)') "    NO2 = ", idx_no2
    write(*,'(A)') "  MICM rate-parameter indices:"
    write(*,'(A,I0)') "    PHOTO.jO3_O1D = ", rp_jo3_o1d
    write(*,'(A,I0)') "    PHOTO.jO3_O   = ", rp_jo3_o
    write(*,'(A,I0)') "    PHOTO.jO2     = ", rp_jo2
    write(*,'(A,I0)') "    PHOTO.jNO2    = ", rp_jno2
  end subroutine print_ordering

  subroutine print_rate_parameters()
    write(*,'(A)') "rate parameters interpreted by MICM"
    write(*,'(A,ES24.16E3)') "  PHOTO.jO2     = ", get_rate(rp_jo2)
    write(*,'(A,ES24.16E3)') "  PHOTO.jO3_O   = ", get_rate(rp_jo3_o)
    write(*,'(A,ES24.16E3)') "  PHOTO.jO3_O1D = ", get_rate(rp_jo3_o1d)
    write(*,'(A,ES24.16E3)') "  PHOTO.jNO2    = ", get_rate(rp_jno2)
  end subroutine print_rate_parameters

  subroutine print_state(label)
    character(len=*), intent(in) :: label
    real(real64) :: o_atoms, n_atoms

    o_atoms = 2.0_real64 * get_conc(idx_o2) + get_conc(idx_o) + &
        get_conc(idx_o1d) + 3.0_real64 * get_conc(idx_o3)
    n_atoms = get_conc(idx_no) + get_conc(idx_no2)

    write(*,'(A)') trim(label)
    write(*,'(A,ES24.16E3)') "  O   = ", get_conc(idx_o)
    write(*,'(A,ES24.16E3)') "  O2  = ", get_conc(idx_o2)
    write(*,'(A,ES24.16E3)') "  O3  = ", get_conc(idx_o3)
    write(*,'(A,ES24.16E3)') "  O1D = ", get_conc(idx_o1d)
    write(*,'(A,ES24.16E3)') "  NO  = ", get_conc(idx_no)
    write(*,'(A,ES24.16E3)') "  NO2 = ", get_conc(idx_no2)
    write(*,'(A,ES24.16E3)') "  O atoms = ", o_atoms
    write(*,'(A,ES24.16E3)') "  N atoms = ", n_atoms
  end subroutine print_state

end program chapman_nox_level2_reproducer
