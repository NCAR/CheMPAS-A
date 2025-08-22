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
    use mpas_derived_types, only : mpas_pool_type, &
                                   field1DReal, &
                                   field2DReal, &
                                   field3DReal

    private

    public :: MPAS_var

    integer, parameter, public :: MPAS_VAR_OK            =     0, &
                                  MPAS_VAR_INVALID_TYPE  = -1000, &
                                  MPAS_VAR_INVALID_NDIMS = -1001

    integer, parameter, public :: MPAS_VAR_INTEGER = 1000, &
                                  MPAS_VAR_REAL    = 1001, &
                                  MPAS_VAR_LOGICAL = 1002, &
                                  MPAS_VAR_CHAR    = 1003

! TO DO:
!        - handle packages
!        - maybe dimensions can be a string?
!        - defaultValue
!        - missingValue
!        - time dimension

    type MPAS_var
        ! Public members
        integer :: ndims = 2
        integer :: vartype = MPAS_VAR_REAL
        character(len=ShortStrKIND) :: name = ''
        character(len=ShortStrKIND) :: units = '-'
        character(len=ShortStrKIND) :: description = 'unknown'
        character(len=ShortStrKIND), dimension(:), allocatable :: dimensions
        character(len=ShortStrKIND) :: group = 'default'  ! MGD maybe group only matters
                                                          !   when adding constitutent to
                                                          !   a var_array?
        character(len=ShortStrKIND):: packages = ''

        ! Private members
        type (field1DReal), dimension(:), pointer, private :: real1d => null()
        type (field2DReal), dimension(:), pointer, private :: real2d => null()
        type (field3DReal), dimension(:), pointer, private :: real3d => null()
    contains
        procedure :: add_to_pool
        final :: destroy_mpas_var
    end type MPAS_var

!    type MPAS_var_array
!        character(len=ShortStrKIND) :: name = ''
!
!        !MGD everything else private?
!
!        integer :: ndims = 2
!        integer :: vartype = MPAS_VAR_REAL
!        character(len=ShortStrKIND) :: units = '-'
!        character(len=ShortStrKIND) :: description = 'unknown'
!        character(len=ShortStrKIND), dimension(:), allocatable :: dimensions
!        character(len=ShortStrKIND) :: group = 'default'  ! MGD maybe group only matters
!                                                          !   when adding constitutent to
!                                                          !   a var_array?
!        character(len=ShortStrKIND):: packages = ''
!    contains
!        procedure :: add_constituent
!        procedure :: instantiate_in_pools
!        final :: destroy_mpas_var
!    end type MPAS_var_array


contains


    subroutine add_to_pool(self, block, pool, ierr)

        use mpas_pool_routines, only : mpas_pool_add_field
        use mpas_derived_types, only : block_type
        use mpas_attlist, only : mpas_add_att

        class (MPAS_var), intent(inout) :: self
        type (block_type), target, intent(in) :: block
        type (mpas_pool_type), target, intent(inout) :: pool
        integer, intent(out), optional :: ierr


        if (present(ierr)) ierr = MPAS_VAR_OK

        select case (self % vartype)
        case (MPAS_VAR_REAL)

            select case (self % ndims)
            case (1)
                if (.not. associated(self % real1d)) then
                    allocate(self % real1d(1))

                    self % real1d(1) % block => block
                    self % real1d(1) % fieldName = self % name
                    self % real1d(1) % dimNames(1) = self % dimensions(1)
                    self % real1d(1) % defaultValue = 0.0
                    self % real1d(1) % missingValue = 0.0
                    self % real1d(1) % isDecomposed = .true.
                    self % real1d(1) % hasTimeDimension = .false.
                    self % real1d(1) % isActive = .false. ! MGD should depend on packages
                    self % real1d(1) % isVarArray = .false.
                    self % real1d(1) % isPersistent = .true.

                    allocate(self % real1d(1) % attLists(1))
                    allocate(self % real1d(1) % attLists(1) % attList)

                    call mpas_add_att(self % real1d(1) % attLists(1) % attList, 'units', self % units)
                    call mpas_add_att(self % real1d(1) % attLists(1) % attList, 'long_name', self % description)
                end if

                call mpas_pool_add_field(pool, trim(self % real1d(1) % fieldName), self % real1d)

            case (2)
                if (.not. associated(self % real2d)) then
                    allocate(self % real2d(1))

                    self % real2d(1) % block => block
                    self % real2d(1) % fieldName = self % name
                    self % real2d(1) % dimNames(1) = self % dimensions(1)
                    self % real2d(1) % dimNames(2) = self % dimensions(2)
                    self % real2d(1) % defaultValue = 0.0
                    self % real2d(1) % missingValue = 0.0
                    self % real2d(1) % isDecomposed = .true.
                    self % real2d(1) % hasTimeDimension = .false.
                    self % real2d(1) % isActive = .true.
                    self % real2d(1) % isVarArray = .false.
                    self % real2d(1) % isPersistent = .true.

                    allocate(self % real2d(1) % attLists(1))
                    allocate(self % real2d(1) % attLists(1) % attList)

                    call mpas_add_att(self % real2d(1) % attLists(1) % attList, 'units', self % units)
                    call mpas_add_att(self % real2d(1) % attLists(1) % attList, 'long_name', self % description)
                end if

                call mpas_pool_add_field(pool, trim(self % real2d(1) % fieldName), self % real2d)

            case (3)
                if (.not. associated(self % real3d)) then
                    allocate(self % real3d(1))

                    self % real3d(1) % block => block
                    self % real3d(1) % fieldName = self % name
                    self % real3d(1) % dimNames(1) = self % dimensions(1)
                    self % real3d(1) % dimNames(2) = self % dimensions(2)
                    self % real3d(1) % dimNames(3) = self % dimensions(3)
                    self % real3d(1) % defaultValue = 0.0
                    self % real3d(1) % missingValue = 0.0
                    self % real3d(1) % isDecomposed = .true.
                    self % real3d(1) % hasTimeDimension = .false.
                    self % real3d(1) % isActive = .true.
                    self % real3d(1) % isVarArray = .false.
                    self % real3d(1) % isPersistent = .true.

                    allocate(self % real3d(1) % attLists(1))
                    allocate(self % real3d(1) % attLists(1) % attList)

                    call mpas_add_att(self % real3d(1) % attLists(1) % attList, 'units', self % units)
                    call mpas_add_att(self % real3d(1) % attLists(1) % attList, 'long_name', self % description)
                end if

                call mpas_pool_add_field(pool, trim(self % real3d(1) % fieldName), self % real3d)

            case default
                if (present(ierr)) ierr = MPAS_VAR_INVALID_NDIMS
            end select

        case default
            if (present(ierr)) ierr = MPAS_VAR_INVALID_TYPE
        end select

    end subroutine add_to_pool


    subroutine destroy_mpas_var(self)

        type (MPAS_var), intent(inout) :: self

        if (allocated(self % dimensions)) deallocate(self % dimensions)

    end subroutine destroy_mpas_var


    subroutine mpas_split_string_new(string, delimiter, subStrings)

        implicit none

        ! Arguments
        character(len=*), intent(in) :: string
        character, intent(in) :: delimiter
        character(len=:), dimension(:), allocatable, intent(inout) :: subStrings

        ! Local variables
        integer :: i, j, n_strs, strlen, max_strlen


        i = 1
        n_strs = 0
        strlen = 0
        max_strlen = 1
        PARSE_LOOP: do while (i <= len(string))
            do while (string(i:i) == delimiter)
                i = i + 1
                if (i > len(string)) exit PARSE_LOOP
            end do

            n_strs = n_strs + 1
            strlen = 0
            do while (string(i:i) /= delimiter)
                i = i + 1
                strlen = strlen + 1
                if (i > len(string)) exit
            end do
            max_strlen = max(strlen, max_strlen)
        end do PARSE_LOOP

        if (allocated(subStrings)) deallocate(subStrings)
        allocate(character(len=max_strlen) :: subStrings(n_strs))

        i = 1
        n_strs = 0
        COPY_LOOP: do while (i <= len(string))
            do while (string(i:i) == delimiter)
                i = i + 1
                if (i > len(string)) exit COPY_LOOP
            end do

            n_strs = n_strs + 1
            j = 1
            do while (string(i:i) /= delimiter)
                subStrings(n_strs)(j:j) = string(i:i)
                i = i + 1
                j = j + 1
                if (i > len(string)) exit
            end do
        end do COPY_LOOP

    end subroutine mpas_split_string_new

end module mpas_runtime_vars
