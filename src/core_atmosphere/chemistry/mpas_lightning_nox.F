! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_lightning_nox
!
!> \brief Lightning NOx source parameterization
!> \author CheMPAS Developers
!> \date   March 2026
!> \details
!>  Injects NO into MPAS scalars where vertical velocity exceeds a
!>  threshold and altitude is within a prescribed range. This is a
!>  source-only module; the chemistry (NO + O3 -> NO2, NO2 + hv -> NO + O3)
!>  is handled by MICM.
!>
!>  The source rate scales linearly with updraft strength:
!>    S = source_rate * max(0, w - w_threshold) / w_ref  [ppbv/s]
!>
!>  At init time the module looks up index_qNO in the scalars pool.
!>  If qNO is not present (mechanism doesn't use it), the module is a no-op.
!>
!>  Adapted from DAVINCI-MPAS mpas_davinci_lightning_nox.F
!
!-----------------------------------------------------------------------
module mpas_lightning_nox

    use mpas_kind_types, only: RKIND, StrKIND

    implicit none
    private

    public :: lightning_nox_init
    public :: lightning_nox_inject

    ! Module state
    logical, save :: lnox_active = .false.
    integer, save :: idx_qNO = -1

    ! Parameters (set from namelist at init)
    real(kind=RKIND), save :: source_rate_ppbv = 0.0_RKIND   ! peak rate [ppbv/s]
    real(kind=RKIND), save :: w_threshold      = 5.0_RKIND   ! [m/s]
    real(kind=RKIND), save :: w_ref            = 10.0_RKIND  ! reference w for scaling [m/s]
    real(kind=RKIND), save :: z_min            = 5000.0_RKIND ! [m]
    real(kind=RKIND), save :: z_max            = 12000.0_RKIND ! [m]

contains

    !-----------------------------------------------------------------------
    !  routine lightning_nox_init
    !
    !> \brief Initialize lightning NOx source; no-op if qNO not present
    !-----------------------------------------------------------------------
    subroutine lightning_nox_init(configs, state)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_config, mpas_pool_get_dimension
        use mpas_log, only: mpas_log_write

        type(mpas_pool_type), intent(in) :: configs
        type(mpas_pool_type), intent(in) :: state

        integer, pointer :: index_qNO_ptr
        real(kind=RKIND), pointer :: cfg_ptr

        ! Look for qNO in the scalars pool
        nullify(index_qNO_ptr)
        call mpas_pool_get_dimension(state, 'index_qNO', index_qNO_ptr)

        if (.not. associated(index_qNO_ptr)) then
            call mpas_log_write('[LNOx] qNO not found in scalars; lightning source disabled.')
            lnox_active = .false.
            return
        end if

        if (index_qNO_ptr < 1) then
            call mpas_log_write('[LNOx] qNO index invalid; lightning source disabled.')
            lnox_active = .false.
            return
        end if

        idx_qNO = index_qNO_ptr

        ! Read namelist parameters
        nullify(cfg_ptr)
        call mpas_pool_get_config(configs, 'config_lnox_source_rate', cfg_ptr)
        if (associated(cfg_ptr)) source_rate_ppbv = cfg_ptr

        nullify(cfg_ptr)
        call mpas_pool_get_config(configs, 'config_lnox_w_threshold', cfg_ptr)
        if (associated(cfg_ptr)) w_threshold = cfg_ptr

        nullify(cfg_ptr)
        call mpas_pool_get_config(configs, 'config_lnox_w_ref', cfg_ptr)
        if (associated(cfg_ptr)) w_ref = cfg_ptr

        nullify(cfg_ptr)
        call mpas_pool_get_config(configs, 'config_lnox_z_min', cfg_ptr)
        if (associated(cfg_ptr)) z_min = cfg_ptr

        nullify(cfg_ptr)
        call mpas_pool_get_config(configs, 'config_lnox_z_max', cfg_ptr)
        if (associated(cfg_ptr)) z_max = cfg_ptr

        if (source_rate_ppbv <= 0.0_RKIND) then
            call mpas_log_write('[LNOx] source_rate <= 0; lightning source disabled.')
            lnox_active = .false.
            return
        end if

        lnox_active = .true.

        call mpas_log_write('[LNOx] Lightning NOx source enabled:')
        call mpas_log_write('[LNOx]   qNO index    = $i', intArgs=(/idx_qNO/))
        call mpas_log_write('[LNOx]   source_rate  = $r ppbv/s', realArgs=(/source_rate_ppbv/))
        call mpas_log_write('[LNOx]   w_threshold  = $r m/s', realArgs=(/w_threshold/))
        call mpas_log_write('[LNOx]   w_ref        = $r m/s', realArgs=(/w_ref/))
        call mpas_log_write('[LNOx]   z_min        = $r m', realArgs=(/z_min/))
        call mpas_log_write('[LNOx]   z_max        = $r m', realArgs=(/z_max/))

    end subroutine lightning_nox_init


    !-----------------------------------------------------------------------
    !  routine lightning_nox_inject
    !
    !> \brief Inject NO source into scalars where w > threshold
    !> \details
    !>  Modifies scalars(index_qNO, k, iCell) in place. Called before
    !>  MICM solves the chemistry. Source is in mixing ratio (ppbv * 1e-9).
    !>
    !>  Uses layer-midpoint height from zgrid for altitude gating and
    !>  w interpolated to layer midpoints from interface values.
    !-----------------------------------------------------------------------
    subroutine lightning_nox_inject(dt, scalars, w, zgrid, &
                                    nVertLevels, nCells, time_lev)

        use mpas_log, only: mpas_log_write

        real(kind=RKIND), intent(in) :: dt
        real(kind=RKIND), intent(inout) :: scalars(:,:,:)
        real(kind=RKIND), intent(in) :: w(:,:)             ! (nVertLevelsP1, nCells)
        real(kind=RKIND), intent(in) :: zgrid(:,:)         ! (nVertLevelsP1, nCells)
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: nCells
        integer, intent(in) :: time_lev

        integer :: iCell, k, n_activated
        real(kind=RKIND) :: w_mid, z_mid, delta_q
        logical, save :: first_call = .true.

        if (.not. lnox_active) return

        n_activated = 0

        do iCell = 1, nCells
            do k = 1, nVertLevels

                ! Layer-midpoint vertical velocity (average of interfaces)
                w_mid = 0.5_RKIND * (w(k, iCell) + w(k+1, iCell))

                ! Layer-midpoint height
                z_mid = 0.5_RKIND * (zgrid(k, iCell) + zgrid(k+1, iCell))

                if (w_mid > w_threshold .and. &
                    z_mid >= z_min .and. z_mid <= z_max) then

                    ! Source in mixing ratio: rate [ppbv/s] * scale * dt * 1e-9
                    delta_q = source_rate_ppbv * (w_mid - w_threshold) / w_ref &
                              * dt * 1.0e-9_RKIND

                    scalars(idx_qNO, k, iCell) = scalars(idx_qNO, k, iCell) + delta_q
                    n_activated = n_activated + 1
                end if

            end do
        end do

        if (first_call) then
            call mpas_log_write('[LNOx] First injection: $i cell-levels activated', &
                                intArgs=(/n_activated/))
            first_call = .false.
        end if

    end subroutine lightning_nox_inject

end module mpas_lightning_nox
