! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
! SPDX-License-Identifier: BSD-3-Clause

module mpas_chem_timing
    use iso_fortran_env, only: int64, real64
    use mpas_kind_types, only: RKIND
    use mpas_derived_types, only: dm_info
    use mpas_dmpar, only: mpas_dmpar_sum_real, mpas_dmpar_min_real, mpas_dmpar_max_real
    use mpas_log, only: mpas_log_write

    implicit none
    private
    public :: chem_library_timer

    ! Rank-local wall time. No MPI operations or logging in start/stop: ranks
    ! can make different numbers of calls (including none on nighttime ranks).
    type :: chem_library_timer
        private
        integer(int64) :: started = 0_int64, ticks = 0_int64, calls = 0_int64
    contains
        procedure :: reset => timer_reset
        procedure :: start => timer_start
        procedure :: stop => timer_stop
        procedure :: report => timer_report
    end type chem_library_timer

contains

    subroutine timer_reset(timer)
        class(chem_library_timer), intent(inout) :: timer
        timer%started = 0_int64
        timer%ticks = 0_int64
        timer%calls = 0_int64
    end subroutine timer_reset

    subroutine timer_start(timer)
        class(chem_library_timer), intent(inout) :: timer
        call system_clock(count=timer%started)
    end subroutine timer_start

    subroutine timer_stop(timer)
        class(chem_library_timer), intent(inout) :: timer
        integer(int64) :: finished, maximum, elapsed
        call system_clock(count=finished, count_max=maximum)
        elapsed = finished - timer%started
        if (finished < timer%started) elapsed = (maximum - timer%started) + finished + 1_int64
        timer%ticks = timer%ticks + elapsed
        timer%calls = timer%calls + 1_int64
    end subroutine timer_stop

    ! Call once at a common finalization point, even on ranks with no calls.
    ! Seconds are accumulated per rank, then reduced; they are not CPU time
    ! and do not include waiting for other ranks in these reporting reductions.
    subroutine timer_report(timer, name, dminfo)
        class(chem_library_timer), intent(in) :: timer
        character(len=*), intent(in) :: name
        type(dm_info), intent(in), optional :: dminfo
        integer(int64) :: rate
        real(kind=RKIND) :: seconds, minimum, mean, maximum, local_calls, total_calls
        character(len=512) :: message

        call system_clock(count_rate=rate)
        if (rate <= 0_int64) then
            call mpas_log_write('[CHEM TIMER] ' // trim(name) // ': system_clock unavailable')
            return
        end if
        seconds = real(real(timer%ticks, real64) / real(rate, real64), RKIND)
        local_calls = real(timer%calls, RKIND)
        if (present(dminfo)) then
            call mpas_dmpar_min_real(dminfo, seconds, minimum)
            call mpas_dmpar_sum_real(dminfo, seconds, mean)
            mean = mean / dminfo%nprocs
            call mpas_dmpar_max_real(dminfo, seconds, maximum)
            call mpas_dmpar_sum_real(dminfo, local_calls, total_calls)
        else
            minimum = seconds
            mean = seconds
            maximum = seconds
            total_calls = local_calls
        end if
        write(message, '(a,a,a,3(1x,es16.8),a,i0)') &
            '[CHEM TIMER] ', trim(name), ': rank seconds min/mean/max =', &
            minimum, mean, maximum, '; total calls = ', nint(total_calls, int64)
        call mpas_log_write(trim(message))
    end subroutine timer_report

end module mpas_chem_timing
