! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_runtime_vars
!
!> \brief Provides facilities for defining new fields / variables at runtime
!> \author Michael Duda
!> \date   17 July 2025
!> \details
!>  This module provides facilities for defining new fields (e.g., field2DReal)
!>  at runtime.
!>
!
!-----------------------------------------------------------------------
module mpas_runtime_vars

    use mpas_kind_types, only : ShortStrKIND
    use mpas_derived_types, only : mpas_pool_type

    private

    public :: MPAS_var

    integer, parameter, public :: MPAS_VAR_INT = 1000, &
                                  MPAS_VAR_REAL = 1001

    type MPAS_var
        integer :: ndims = 2
        integer :: vartype = MPAS_VAR_REAL
        character(len=ShortStrKIND) :: name = ''
        character(len=ShortStrKIND) :: units = '-'
        character(len=ShortStrKIND) :: description = 'unknown'
        character(len=ShortStrKIND), dimension(:), allocatable :: dimensions
        character(len=ShortStrKIND) :: group = 'default'
        character(len=ShortStrKIND):: packages = ''
    contains
        procedure :: instantiate_in_pools
        final :: destroy_mpas_var
    end type MPAS_var


contains


    subroutine instantiate_in_pools(self, block, pool1, pool2)

        use mpas_pool_routines, only : mpas_pool_add_field
        use mpas_derived_types, only : field2DReal, block_type
        use mpas_attlist, only : mpas_add_att

        class (MPAS_var), intent(inout) :: self
        type (block_type), target, intent(in) :: block
        ! type (mpas_pool_type), dimension(:), target, intent(in) :: pools
        type (mpas_pool_type), target, intent(inout) :: pool1, pool2

        type (mpas_pool_type), pointer :: p
        type (field2DReal), dimension(:), pointer :: f
        integer :: i


        allocate(f(1))

        f(1) % block => block
        f(1) % fieldName = self % name
        f(1) % dimNames(1) = self % dimensions(1)
        f(1) % dimNames(2) = self % dimensions(2)
        f(1) % defaultValue = 0.0
        f(1) % missingValue = 0.0
        f(1) % isDecomposed = .true.
        f(1) % hasTimeDimension = .false.
        f(1) % isActive = .true.
        f(1) % isVarArray = .false.
        f(1) % isPersistent = .true.

        allocate(f(1) % attLists(1))
        allocate(f(1) % attLists(1) % attList)

        call mpas_add_att(f(1) % attLists(1) % attList, 'units', self % units)
        call mpas_add_att(f(1) % attLists(1) % attList, 'long_name', self % description)

        ! do i = 1, self % ndims
        !     p => pools(i)
        !     call mpas_pool_add_field(p, f(1) % fieldName, f)
        ! end do
        p => pool1
        call mpas_pool_add_field(p, f(1) % fieldName, f)
        p => pool2
        call mpas_pool_add_field(p, f(1) % fieldName, f)

    end subroutine instantiate_in_pools


    subroutine destroy_mpas_var(self)

        type (MPAS_var), intent(inout) :: self

        deallocate(self % dimensions)

    end subroutine destroy_mpas_var

end module mpas_runtime_vars
