! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_chem_constants
!
!> \brief Shared chemistry constants
!-----------------------------------------------------------------------
module mpas_chem_constants

    use iso_fortran_env, only: real64
    use mpas_constants, only: pii, rgas, rv, rvord

    implicit none

    private

    public :: CHEM_GAS_R
    public :: CHEM_NA
    public :: CHEM_M_AIR
    public :: CHEM_M_O3
    public :: CHEM_O2_VMR
    public :: CHEM_N2_VMR
    public :: CHEM_RHO_H2O
    public :: pii
    public :: rgas
    public :: rv
    public :: rvord

    real(real64), parameter :: CHEM_GAS_R   = 8.314462618_real64   ! J mol-1 K-1
    real(real64), parameter :: CHEM_NA      = 6.02214076e23_real64 ! mol-1
    real(real64), parameter :: CHEM_M_AIR   = 0.0289644_real64     ! kg mol-1 dry air
    real(real64), parameter :: CHEM_M_O3    = 0.04800_real64       ! kg mol-1 ozone
    real(real64), parameter :: CHEM_O2_VMR  = 0.2095_real64        ! O2 volume mixing ratio
    real(real64), parameter :: CHEM_N2_VMR  = 0.7808_real64        ! N2 volume mixing ratio
    real(real64), parameter :: CHEM_RHO_H2O = 1000.0_real64        ! kg m-3

end module mpas_chem_constants
