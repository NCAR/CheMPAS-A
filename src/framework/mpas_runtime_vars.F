! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_runtime_vars
!
!> \brief Provides facilities for defining fields / variables at runtime
!> \author Michael Duda
!> \date   17 July 2025
!> \details
!>  This module provides facilities for defining fields (e.g., field2DReal)
!>  at runtime.
!>
!
!-----------------------------------------------------------------------

#define PTR_CAST(PTR, RES, TYPE) if (associated(PTR)) then; select type(PTR); type is (TYPE); RES => PTR; end select; end if

module mpas_runtime_vars

    use mpas_kind_types, only : StrKIND, ShortStrKIND
    use mpas_derived_types, only : mpas_pool_type, &
                                   field1DInteger, &
                                   field2DInteger, &
                                   field1DReal, &
                                   field2DReal, &
                                   field3DReal

    implicit none

    private

    public :: MPAS_var, &
              MPAS_var_array

    integer, parameter, public :: MPAS_VAR_OK            =     0, &
                                  MPAS_VAR_INVALID_TYPE  = -1000, &
                                  MPAS_VAR_INVALID_NDIMS = -1001

    integer, parameter, public :: MPAS_VAR_INTEGER   = 1000, &
                                  MPAS_VAR_REAL      = 1001, &
                                  MPAS_VAR_LOGICAL   = 1002, &
                                  MPAS_VAR_CHARACTER = 1003

! TO DO:
!        - defaultValue
!        - missingValue

    type linked_list_node
        character(len=ShortStrKIND) :: key
        class(*), pointer :: item => null()
        type (linked_list_node), pointer :: next => null()
    end type linked_list_node


    type linked_list
        ! Private members
        integer, private :: nitems = 0
        type (linked_list_node), pointer, private :: head => null()
        type (linked_list_node), pointer, private :: iter_next => null()
    contains
        procedure, pass :: length => ll_length
        procedure, pass :: insert => ll_insert
        procedure, pass :: remove => ll_remove
        procedure, pass :: find => ll_find
        procedure, pass :: begin_iteration => ll_begin_iteration
        procedure, pass :: get_next_item => ll_get_next_item
        final :: ll_destroy
    end type linked_list


    type MPAS_var
        ! Public members
        integer :: vartype = MPAS_VAR_REAL
        character(len=ShortStrKIND) :: name = ''
        character(len=ShortStrKIND) :: name_in_code = ''
        character(len=ShortStrKIND) :: units = '-'
        character(len=StrKIND) :: description = 'unknown'
        character(len=StrKIND) :: dimensions = 'undefined'
        logical :: hasTimeDimension = .true.
        integer :: nTimeLevels = 1
        logical :: isPersistent = .true.
        character(len=ShortStrKIND) :: group = ''  ! Only needed for constituents of an MPAS_var_array
        character(len=StrKIND) :: packages = ''

        ! Private members
        type (field1DInteger), pointer, private :: int1d => null()
        type (field1DInteger), dimension(:), pointer, private :: int1d_arr => null()
        type (field2DInteger), pointer, private :: int2d => null()
        type (field2DInteger), dimension(:), pointer, private :: int2d_arr => null()
        type (field1DReal), pointer, private :: real1d => null()
        type (field1DReal), dimension(:), pointer, private :: real1d_arr => null()
        type (field2DReal), pointer, private :: real2d => null()
        type (field2DReal), dimension(:), pointer, private :: real2d_arr => null()
        type (field3DReal), pointer, private :: real3d => null()
        type (field3DReal), dimension(:), pointer, private :: real3d_arr => null()
    contains
        procedure :: add_to_pool => add_to_pool_var
        procedure :: reset => reset_var
        final :: destroy_mpas_var
    end type MPAS_var


    type MPAS_var_array
        ! Public members
        character(len=ShortStrKIND) :: name = ''
        character(len=ShortStrKIND) :: name_in_code = ''
        character(len=StrKIND) :: packages = ''
        integer :: nTimeLevels = 1
        logical :: isPersistent = .true.

        ! Private members
        integer, private :: ndims = -1
        integer, private :: vartype = MPAS_VAR_INVALID_TYPE
        character(len=ShortStrKIND), dimension(:), allocatable, private :: dimensions
        type (linked_list), allocatable :: const_list
        logical :: hasTimeDimension = .true.
        type (field3DReal), pointer, private :: real3d => null()
        type (field3DReal), dimension(:), pointer, private :: real3d_arr => null()
    contains
        procedure :: add_constituent
        procedure :: add_to_pool => add_to_pool_var_array
        procedure :: reset => reset_var_array
        final :: destroy_mpas_var_array
    end type MPAS_var_array


contains


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine add_to_pool_var(self, block, fieldPool, packagePool, useNameInCode, ierr)

        use mpas_pool_routines, only : mpas_pool_add_field
        use mpas_derived_types, only : block_type
        use mpas_kind_types, only : RKIND
        use mpas_attlist, only : mpas_add_att

        implicit none

        ! Arguments
        class (MPAS_var), intent(inout) :: self
        type (block_type), target, intent(in) :: block
        type (mpas_pool_type), target, intent(inout) :: fieldPool
        type (mpas_pool_type), intent(in) :: packagePool
        logical, intent(in), optional :: useNameInCode
        integer, intent(out), optional :: ierr

        ! Local variables
        integer :: i, ndims
        character(len=:), dimension(:), allocatable :: dimensions
        type (field1DInteger), pointer :: int1d_ptr
        type (field2DInteger), pointer :: int2d_ptr
        type (field1DReal), pointer :: real1d_ptr
        type (field2DReal), pointer :: real2d_ptr
        type (field3DReal), pointer :: real3d_ptr
        character(len=ShortStrKIND) :: name_in_pool

        if (present(ierr)) ierr = MPAS_VAR_OK

        if (len_trim(self % name_in_code) == 0) then
            self % name_in_code = self % name
        end if

        if (present(useNameInCode)) then
            name_in_pool = self % name_in_code
        else
            name_in_pool = self % name
        end if

        call mpas_split_string_new(self % dimensions, ' ', dimensions)
        ndims = size(dimensions)

        select case (self % vartype)
        case (MPAS_VAR_INTEGER)

            select case (ndims)
            case (1)
                if ((.not. associated(self % int1d) .and. self % nTimeLevels == 1) &
                    .or. (.not. associated(self % int1d_arr) .and. self % nTimeLevels > 1)) then

                    if (self % nTimeLevels > 1) then
                        allocate(self % int1d_arr(self % nTimeLevels))
                    else
                        allocate(self % int1d)
                    end if

                    do i = 1, self % nTimeLevels
                        if (self % nTimeLevels > 1) then
                            int1d_ptr => self % int1d_arr(i)
                        else
                            int1d_ptr => self % int1d
                        end if

                        int1d_ptr % block => block
                        int1d_ptr % fieldName = self % name
                        int1d_ptr % dimNames(1) = dimensions(1)
                        int1d_ptr % defaultValue = 0
                        int1d_ptr % missingValue = 0
                        int1d_ptr % isDecomposed = is_decomposed_dim(dimensions(1))
                        int1d_ptr % hasTimeDimension = self % hasTimeDimension
                        int1d_ptr % isActive = contains_active_package(self % packages, packagePool)
                        int1d_ptr % isVarArray = .false.
                        int1d_ptr % isPersistent = self % isPersistent

                        allocate(int1d_ptr % attLists(1))
                        allocate(int1d_ptr % attLists(1) % attList)

                        call mpas_add_att(int1d_ptr % attLists(1) % attList, 'units', self % units)
                        call mpas_add_att(int1d_ptr % attLists(1) % attList, 'long_name', self % description)
                    end do

                end if

                if (self % nTimeLevels > 1) then
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % int1d_arr)
                else
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % int1d)
                end if

            case (2)
                if ((.not. associated(self % int2d) .and. self % nTimeLevels == 1) &
                    .or. (.not. associated(self % int2d_arr) .and. self % nTimeLevels > 1)) then

                    if (self % nTimeLevels > 1) then
                        allocate(self % int2d_arr(self % nTimeLevels))
                    else
                        allocate(self % int2d)
                    end if

                    do i = 1, self % nTimeLevels
                        if (self % nTimeLevels > 1) then
                            int2d_ptr => self % int2d_arr(i)
                        else
                            int2d_ptr => self % int2d
                        end if

                        int2d_ptr % block => block
                        int2d_ptr % fieldName = self % name
                        int2d_ptr % dimNames(1) = dimensions(1)
                        int2d_ptr % dimNames(2) = dimensions(2)
                        int2d_ptr % defaultValue = 0
                        int2d_ptr % missingValue = 0
                        int2d_ptr % isDecomposed = is_decomposed_dim(dimensions(2))
                        int2d_ptr % hasTimeDimension = self % hasTimeDimension
                        int2d_ptr % isActive = contains_active_package(self % packages, packagePool)
                        int2d_ptr % isVarArray = .false.
                        int2d_ptr % isPersistent = self % isPersistent

                        allocate(int2d_ptr % attLists(1))
                        allocate(int2d_ptr % attLists(1) % attList)

                        call mpas_add_att(int2d_ptr % attLists(1) % attList, 'units', self % units)
                        call mpas_add_att(int2d_ptr % attLists(1) % attList, 'long_name', self % description)
                    end do

                end if

                if (self % nTimeLevels > 1) then
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % int2d_arr)
                else
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % int2d)
                end if

            case default
                if (present(ierr)) ierr = MPAS_VAR_INVALID_NDIMS
            end select

        case (MPAS_VAR_REAL)

            select case (ndims)
            case (1)
                if ((.not. associated(self % real1d) .and. self % nTimeLevels == 1) &
                    .or. (.not. associated(self % real1d_arr) .and. self % nTimeLevels > 1)) then

                    if (self % nTimeLevels > 1) then
                        allocate(self % real1d_arr(self % nTimeLevels))
                    else
                        allocate(self % real1d)
                    end if

                    do i = 1, self % nTimeLevels
                        if (self % nTimeLevels > 1) then
                            real1d_ptr => self % real1d_arr(i)
                        else
                            real1d_ptr => self % real1d
                        end if

                        real1d_ptr % block => block
                        real1d_ptr % fieldName = self % name
                        real1d_ptr % dimNames(1) = dimensions(1)
                        real1d_ptr % defaultValue = 0.0_RKIND
                        real1d_ptr % missingValue = 0.0_RKIND
                        real1d_ptr % isDecomposed = is_decomposed_dim(dimensions(1))
                        real1d_ptr % hasTimeDimension = self % hasTimeDimension
                        real1d_ptr % isActive = contains_active_package(self % packages, packagePool)
                        real1d_ptr % isVarArray = .false.
                        real1d_ptr % isPersistent = self % isPersistent

                        allocate(real1d_ptr % attLists(1))
                        allocate(real1d_ptr % attLists(1) % attList)

                        call mpas_add_att(real1d_ptr % attLists(1) % attList, 'units', self % units)
                        call mpas_add_att(real1d_ptr % attLists(1) % attList, 'long_name', self % description)
                    end do

                end if

                if (self % nTimeLevels > 1) then
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real1d_arr)
                else
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real1d)
                end if

            case (2)
                if ((.not. associated(self % real2d) .and. self % nTimeLevels == 1) &
                    .or. (.not. associated(self % real2d_arr) .and. self % nTimeLevels > 1)) then

                    if (self % nTimeLevels > 1) then
                        allocate(self % real2d_arr(self % nTimeLevels))
                    else
                        allocate(self % real2d)
                    end if

                    do i = 1, self % nTimeLevels
                        if (self % nTimeLevels > 1) then
                            real2d_ptr => self % real2d_arr(i)
                        else
                            real2d_ptr => self % real2d
                        end if

                        real2d_ptr % block => block
                        real2d_ptr % fieldName = self % name
                        real2d_ptr % dimNames(1) = dimensions(1)
                        real2d_ptr % dimNames(2) = dimensions(2)
                        real2d_ptr % defaultValue = 0.0_RKIND
                        real2d_ptr % missingValue = 0.0_RKIND
                        real2d_ptr % isDecomposed = is_decomposed_dim(dimensions(2))
                        real2d_ptr % hasTimeDimension = self % hasTimeDimension
                        real2d_ptr % isActive = contains_active_package(self % packages, packagePool)
                        real2d_ptr % isVarArray = .false.
                        real2d_ptr % isPersistent = self % isPersistent

                        allocate(real2d_ptr % attLists(1))
                        allocate(real2d_ptr % attLists(1) % attList)

                        call mpas_add_att(real2d_ptr % attLists(1) % attList, 'units', self % units)
                        call mpas_add_att(real2d_ptr % attLists(1) % attList, 'long_name', self % description)
                    end do

                end if

                if (self % nTimeLevels > 1) then
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real2d_arr)
                else
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real2d)
                end if

            case (3)
                if ((.not. associated(self % real3d) .and. self % nTimeLevels == 1) &
                    .or. (.not. associated(self % real3d_arr) .and. self % nTimeLevels > 1)) then

                    if (self % nTimeLevels > 1) then
                        allocate(self % real3d_arr(self % nTimeLevels))
                    else
                        allocate(self % real3d)
                    end if

                    do i = 1, self % nTimeLevels
                        if (self % nTimeLevels > 1) then
                            real3d_ptr => self % real3d_arr(i)
                        else
                            real3d_ptr => self % real3d
                        end if

                        real3d_ptr % block => block
                        real3d_ptr % fieldName = self % name
                        real3d_ptr % dimNames(1) = dimensions(1)
                        real3d_ptr % dimNames(2) = dimensions(2)
                        real3d_ptr % dimNames(3) = dimensions(3)
                        real3d_ptr % defaultValue = 0.0_RKIND
                        real3d_ptr % missingValue = 0.0_RKIND
                        real3d_ptr % isDecomposed = is_decomposed_dim(dimensions(3))
                        real3d_ptr % hasTimeDimension = self % hasTimeDimension
                        real3d_ptr % isActive = contains_active_package(self % packages, packagePool)
                        real3d_ptr % isVarArray = .false.
                        real3d_ptr % isPersistent = self % isPersistent

                        allocate(real3d_ptr % attLists(1))
                        allocate(real3d_ptr % attLists(1) % attList)

                        call mpas_add_att(real3d_ptr % attLists(1) % attList, 'units', self % units)
                        call mpas_add_att(real3d_ptr % attLists(1) % attList, 'long_name', self % description)
                    end do

                end if

                if (self % nTimeLevels > 1) then
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real3d_arr)
                else
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real3d)
                end if

            case default
                if (present(ierr)) ierr = MPAS_VAR_INVALID_NDIMS
            end select

        case default
            if (present(ierr)) ierr = MPAS_VAR_INVALID_TYPE
        end select

        deallocate(dimensions)

    end subroutine add_to_pool_var


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine reset_var(self)

        implicit none

        ! Arguments
        class (MPAS_var), intent(inout) :: self

        ! Local variables
        type (MPAS_var) :: clean

        call destroy_mpas_var(self)
        select type (var => self)
        type is (MPAS_var)
           var = clean
        end select

    end subroutine reset_var


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine destroy_mpas_var(self)

        ! Arguments
        type (MPAS_var), intent(inout) :: self

        ! if (allocated(self % dimensions)) deallocate(self % dimensions)

    end subroutine destroy_mpas_var


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine add_constituent(self, constituent, group)

        ! Arguments
        class (MPAS_var_array), intent(inout) :: self
        type (MPAS_var), intent(in) :: constituent
        character(len=*), intent(in) :: group

        ! Local variables
        integer :: i, ndims
        character(len=:), dimension(:), allocatable :: dimensions
        type (MPAS_var), pointer :: const_copy

        call mpas_split_string_new(constituent % dimensions, ' ', dimensions)
        ndims = size(dimensions)

        if (self % ndims == -1) then

            ! write(0,*) 'Adding first constituent '//trim(constituent % name)
            self % ndims = ndims + 1
            self % vartype = constituent % vartype
            allocate(self % dimensions(self % ndims))
            do i = 1, ndims
                self % dimensions(i+1) = dimensions(i)
            end do
            ! self % dimensions(1) will be set when instantiating in a pool
            allocate(self % const_list)
            ! Assume that the first added constituent will determine whether var_array as
            ! a whole has a time dimension or not
            self % hasTimeDimension = constituent % hasTimeDimension
        else
            ! write(0,*) 'Adding constituent '//trim(constituent % name)
            if ((self % ndims - 1) /= ndims) then
                write(0,*) 'Constituent dimensions incompatible...'
            end if
            if (self % vartype /= constituent % vartype) then
                write(0,*) 'Constituent type incompatible...'
            end if
            do i = 1, ndims
                if (trim(self % dimensions(i+1)) /= trim(dimensions(i))) then
                    write(0,*) 'Constituent dimension incompatible...'
                end if
            end do
        end if

        deallocate(dimensions)

        allocate(const_copy)
        const_copy = constituent

        if (len_trim(const_copy % name_in_code) == 0) then
            const_copy % name_in_code = const_copy % name
        end if

        const_copy % group = trim(group)
        call self % const_list % insert(trim(const_copy % name), const_copy)

    end subroutine add_constituent


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine add_to_pool_var_array(self, block, fieldPool, packagePool, useNameInCode, addIndices, ierr)

        use mpas_pool_routines, only : mpas_pool_add_field, mpas_pool_add_dimension
        use mpas_derived_types, only : block_type
        use mpas_attlist, only : mpas_add_att
        use mpas_kind_types, only : RKIND

        implicit none

        ! Arguments
        class (MPAS_var_array), intent(inout) :: self
        type (block_type), target, intent(in) :: block
        type (mpas_pool_type), target, intent(inout) :: fieldPool
        type (mpas_pool_type), intent(in) :: packagePool
        logical, intent(in), optional :: useNameInCode
        logical, intent(in), optional :: addIndices
        integer, intent(out), optional :: ierr

        ! Local variables
        class (*), pointer :: current, cursor
        type (MPAS_var), pointer :: current_cnst, cursor_cnst
        type (linked_list) :: sorted_list
        integer :: idx, n, ncnst, i
        type (field3DReal), pointer :: real3d_ptr
        character(len=ShortStrKIND) :: name_in_pool
        logical :: addIndices_local


        if (present(ierr)) ierr = MPAS_VAR_OK

        if (len_trim(self % name_in_code) == 0) then
            self % name_in_code = self % name
        end if

        if (present(useNameInCode)) then
            name_in_pool = self % name_in_code
        else
            name_in_pool = self % name
        end if

        if (present(addIndices)) then
            addIndices_local = addIndices
        else
            addIndices_local = .true.
        end if

        call sort_constituents(self % const_list, sorted_list, packagePool)

        write(0,*) 'num_'//trim(self % name)//' = ', sorted_list % length()
        call mpas_pool_add_dimension(fieldPool, 'num_'//trim(self % name), sorted_list % length())

        call sorted_list % begin_iteration()
        cursor => sorted_list % get_next_item()
        current => cursor
        nullify(current_cnst)
        PTR_CAST(current, current_cnst, MPAS_var)
        n = 1
        ncnst = sorted_list % length()

        if (addIndices_local) then

            if (ncnst > 0) then
                write(0,*) trim(current_cnst % group)//'_start = ', 1
                call mpas_pool_add_dimension(fieldPool, trim(current_cnst % group)//'_start', 1)
            end if

            do idx = 1, ncnst
                nullify(cursor_cnst)
                PTR_CAST(cursor, cursor_cnst, MPAS_var)
                write(0,*) 'MGD idx = ', idx, &
                           ' cnst = '//trim(cursor_cnst % name_in_code)//' group = '//trim(cursor_cnst % group)
                call mpas_pool_add_dimension(fieldPool, 'index_'//trim(cursor_cnst % name_in_code), idx)

                if (trim(cursor_cnst % group) /= trim(current_cnst % group)) then
                    write(0,*) trim(current_cnst % group)//'_end = ', idx - 1
                    call mpas_pool_add_dimension(fieldPool, trim(current_cnst % group)//'_end', idx - 1)
                    write(0,*) 'num_'//trim(current_cnst % group)//' = ', n - 1
                    call mpas_pool_add_dimension(fieldPool, 'num_'//trim(current_cnst % group), n - 1)

                    n = 1
                    current => cursor
                    nullify(current_cnst)
                    PTR_CAST(current, current_cnst, MPAS_var)
                    write(0,*) trim(current_cnst % group)//'_start = ', idx
                    call mpas_pool_add_dimension(fieldPool, trim(current_cnst % group)//'_start', idx)
                end if

                n = n + 1
                cursor => sorted_list % get_next_item()
            end do

            if (ncnst > 0) then
                write(0,*) trim(current_cnst % group)//'_end = ', idx - 1
                call mpas_pool_add_dimension(fieldPool, trim(current_cnst % group)//'_end', idx - 1)
                write(0,*) 'num_'//trim(current_cnst % group)//' = ', n - 1
                call mpas_pool_add_dimension(fieldPool, 'num_'//trim(current_cnst % group), n - 1)
            end if

            !
            ! Add index_* dimensions for inactive constituents
            !
            call self % const_list % begin_iteration()
            cursor => self % const_list % get_next_item()
            do while (associated(cursor))
                PTR_CAST(cursor, cursor_cnst, MPAS_var)
                if (.not. contains_active_package(cursor_cnst % packages, packagePool)) then
                    call mpas_pool_add_dimension(fieldPool, 'index_'//trim(cursor_cnst % name_in_code), -1)
                end if
                cursor => self % const_list % get_next_item()
            end do

        end if

        select case (self % vartype)
        case (MPAS_VAR_REAL)

            select case (self % ndims)
            case (3)
                if ((.not. associated(self % real3d) .and. self % nTimeLevels == 1) &
                    .or. (.not. associated(self % real3d_arr) .and. self % nTimeLevels > 1)) then
                    if (self % nTimeLevels > 1) then
                        allocate(self % real3d_arr(self % nTimeLevels))
                    else
                        allocate(self % real3d)
                    end if

                    do i = 1, self % nTimeLevels
                        if (self % nTimeLevels > 1) then
                            real3d_ptr => self % real3d_arr(i)
                        else
                            real3d_ptr => self % real3d
                        end if

                        real3d_ptr % block => block
                        real3d_ptr % fieldName = self % name
                        real3d_ptr % dimNames(1) = 'num_'//trim(self % name)
                        real3d_ptr % dimNames(2) = self % dimensions(2)
                        real3d_ptr % dimNames(3) = self % dimensions(3)
                        real3d_ptr % defaultValue = 0.0_RKIND
                        real3d_ptr % missingValue = 0.0_RKIND
                        real3d_ptr % isDecomposed = is_decomposed_dim(self % dimensions(3))
                        real3d_ptr % hasTimeDimension = self % hasTimeDimension
                        real3d_ptr % isActive = contains_active_package(self % packages, packagePool)
                        real3d_ptr % isVarArray = .true.
                        real3d_ptr % isPersistent = self % isPersistent

                        allocate(real3d_ptr % attLists(ncnst))
                        if (ncnst > 0) allocate(real3d_ptr % attLists(ncnst) % attList)

                        allocate(real3d_ptr % constituentNames(ncnst))

                        call sorted_list % begin_iteration()
                        do idx = 1, ncnst
                            cursor => sorted_list % get_next_item()
                            nullify(cursor_cnst)
                            PTR_CAST(cursor, cursor_cnst, MPAS_var)
                            real3d_ptr % constituentNames(idx) = trim(cursor_cnst % name)
                            call mpas_add_att(real3d_ptr % attLists(idx) % attList, 'units', &
                                              cursor_cnst % units)
                            call mpas_add_att(real3d_ptr % attLists(idx) % attList, 'long_name', &
                                              cursor_cnst % description)
                        end do
                    end do
                end if

                if (self % nTimeLevels > 1) then
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real3d_arr)
                else
                    call mpas_pool_add_field(fieldPool, trim(name_in_pool), self % real3d)
                end if

            case default
                if (present(ierr)) ierr = MPAS_VAR_INVALID_NDIMS
            end select

        case default
            if (present(ierr)) ierr = MPAS_VAR_INVALID_TYPE
        end select

    end subroutine add_to_pool_var_array


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine reset_var_array(self)

        implicit none

        ! Arguments
        class (MPAS_var_array), intent(inout) :: self

        ! Local variables
        type (MPAS_var_array) :: clean

        call destroy_mpas_var_array(self)
        select type (var_arr => self)
        type is (MPAS_var_array)
           var_arr = clean
        end select

    end subroutine reset_var_array


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine destroy_mpas_var_array(self)

        ! Arguments
        type (MPAS_var_array), intent(inout) :: self

        ! Local variables
        class (*), pointer :: cursor
        type (MPAS_var), pointer :: cursor_cnst

        if (allocated(self % dimensions)) deallocate(self % dimensions)

        if (allocated(self % const_list)) then
            call self % const_list % begin_iteration()
            cursor => self % const_list % get_next_item()
            do while (associated(cursor))
                nullify(cursor_cnst)
                PTR_CAST(cursor, cursor_cnst, MPAS_var)

                write(0,*) 'Deallocating '//trim(cursor_cnst % name)
                deallocate(cursor_cnst)

                cursor => self % const_list % get_next_item()
            end do

            deallocate(self % const_list)
        end if

    end subroutine destroy_mpas_var_array


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine sort_constituents(orig_list, sorted_list, packagePool)

        implicit none

        ! Arguments
        type (linked_list), intent(inout) :: orig_list
        type (linked_list), intent(inout) :: sorted_list
        type (mpas_pool_type), intent(in) :: packagePool

        ! Local variables
        class (*), pointer :: x, cursor, current
        type (MPAS_var), pointer :: cursor_cnst, current_cnst
        type (linked_list) :: tmp_list


        !
        ! Copy orig_list to tmp_list, removing all inactive variables
        !
        call orig_list % begin_iteration()
        cursor => orig_list % get_next_item()
        do while (associated(cursor))
            nullify(cursor_cnst)
            PTR_CAST(cursor, cursor_cnst, MPAS_var)

            if (contains_active_package(cursor_cnst % packages, packagePool)) then
                ! write(0,*) trim(cnst % name) // ' contains active packages and will be used'
                call tmp_list % insert(trim(cursor_cnst % name), cursor_cnst)
            end if

            cursor => orig_list % get_next_item()
        end do

        !
        ! Sort items from tmp_list into sorted_list, placing variables with matching groups
        ! next to each other
        !
        do while (tmp_list % length() > 0)

            ! Get item at beginning of tmp_list
            call tmp_list % begin_iteration()
            current => tmp_list % get_next_item()
            nullify(current_cnst)
            PTR_CAST(current, current_cnst, MPAS_var)

            ! Add first item to sorted_list and remove it from tmp_list
            write(0,*) 'Adding '//trim(current_cnst % name)//' to sorted_list ('//trim(current_cnst % group)//')'
            call sorted_list % insert(trim(current_cnst % name), current_cnst)
            cursor => tmp_list % get_next_item()
            x => tmp_list % remove(trim(current_cnst % name))

            do while (associated(cursor))
                nullify(cursor_cnst)
                PTR_CAST(cursor, cursor_cnst, MPAS_var)
                if (trim(cursor_cnst % group) == trim(current_cnst % group)) then
                    write(0,*) 'Adding '//trim(cursor_cnst % name)//' to sorted_list ('//trim(cursor_cnst % group)//')'
                    call sorted_list % insert(trim(cursor_cnst % name), cursor_cnst)
                    cursor => tmp_list % get_next_item()
                    x => tmp_list % remove(trim(cursor_cnst % name))
                else
                    cursor => tmp_list % get_next_item()
                end if
            end do

        end do

    end subroutine sort_constituents


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
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
            subStrings(n_strs)(j:max_strlen) = ''
        end do COPY_LOOP

    end subroutine mpas_split_string_new


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    logical function contains_active_package(packageList, packagePool)

        use mpas_pool_routines, only : mpas_pool_get_package

        implicit none

        ! Arguments
        character(len=*), intent(in) :: packageList
        type (mpas_pool_type), intent(in) :: packagePool

        ! Local variables
        integer :: i
        logical, pointer :: active
        character(len=:), dimension(:), allocatable :: packages


        contains_active_package = .false.

        !
        ! If there are no packages at all, consider that equivalent to having
        ! an active package
        !
        if (len_trim(packageList) == 0) then
            contains_active_package = .true.
            return
        end if

        !
        ! If there are packages, check each one
        !
        call mpas_split_string_new(packageList, ';', packages)

        do i = 1, size(packages)
            ! write(0,*) 'Considering package '//trim(packages(i))
            nullify(active)
            call mpas_pool_get_package(packagePool, trim(packages(i))//'Active', active)
            if (associated(active)) then
                if (active) then
                    contains_active_package = .true.
                    exit
                end if
            else
                write(0,*) 'Package does not exist: '//trim(packages(i))
            end if
        end do

        deallocate(packages)

    end function contains_active_package


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    integer function ll_length(self) result(length)

        ! Arguments
        class (linked_list), intent(in) :: self

        length = self % nitems

    end function ll_length


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine ll_insert(self, key, item)

        ! Arguments
        class (linked_list), intent(inout) :: self
        character(len=*), intent(in) :: key
        class(*), target, intent(in) :: item

        ! Local variables
        type (linked_list_node), pointer :: cursor


        if (.not. associated(self % head)) then
            allocate(self % head)
            self % head % key = key
            self % head % item => item
        else
            cursor => self % head
            do while (associated(cursor % next))
                cursor => cursor % next
            end do
            allocate(cursor % next)
            cursor => cursor % next
            cursor % key = key
            cursor % item => item
        end if

        self % nitems = self % nitems + 1

    end subroutine ll_insert


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    function ll_remove(self, key) result(item)

        ! Arguments
        class (linked_list), intent(inout) :: self
        character(len=*), intent(in) :: key

        ! Return value
        class(*), pointer :: item

        ! Local variables
        type (linked_list_node), pointer :: cursor, prev


        nullify(item)

        nullify(prev)
        cursor => self % head
        do while (associated(cursor))
            if (trim(cursor % key) == trim(key)) then
                item => cursor % item
                if (associated(prev)) then
                    prev % next => cursor % next
                else
                    self % head => cursor % next
                end if
                self % nitems = self % nitems - 1
                deallocate(cursor)
                exit
            end if
            prev => cursor
            cursor => cursor % next
        end do

    end function ll_remove


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    function ll_find(self, key) result(item)

        ! Arguments
        class (linked_list), intent(inout) :: self
        character(len=*), intent(in) :: key

        ! Return value
        class(*), pointer :: item

        ! Local variables
        type (linked_list_node), pointer :: cursor


        nullify(item)

        cursor => self % head
        do while (associated(cursor))
            if (trim(cursor % key) == trim(key)) then
                item => cursor % item
                exit
            end if
            cursor => cursor % next
        end do

    end function ll_find


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine ll_begin_iteration(self)

        ! Arguments
        class (linked_list), intent(inout) :: self

        self % iter_next => self % head

    end subroutine ll_begin_iteration


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    function ll_get_next_item(self) result(item)

        ! Arguments
        class (linked_list), intent(inout) :: self

        ! Return value
        class(*), pointer :: item


        if (associated(self % iter_next)) then
            item => self % iter_next % item
            self % iter_next => self % iter_next % next
        else
            nullify(item)
        end if

    end function ll_get_next_item


    !-----------------------------------------------------------------------
    !  routine foo_bar
    !
    !> \brief Foos the bar
    !> \author Michael Duda
    !> \date   22 August 2028
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  MPAS_streamInfo_type. The new instance has valid methods init(), query(),
    !>  and finalize() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    subroutine ll_destroy(self)

        ! Arguments
        type (linked_list), intent(inout) :: self

        ! Local variables
        type (linked_list_node), pointer :: cursor


        do while (associated(self % head))
            cursor => self % head % next
            deallocate(self % head)
            self % head => cursor
        end do

        nullify(self % iter_next)
        self % nitems = 0

    end subroutine ll_destroy


    !-----------------------------------------------------------------------
    !  routine is_decomposed_dim
    !
    !> \brief Returns logical indicating whether a dimension is decomposed
    !> \author Michael Duda
    !> \date   13 October 2025
    !> \details
    !>  This function returns a logical value indicating whether the given
    !>  dimension name corresponds to a decomposed dimension or not.
    !
    !-----------------------------------------------------------------------
    function is_decomposed_dim(dimname) result(is_decomposed)

        ! Arguments
        character(len=*), intent(in) :: dimname

        ! Return value
        logical :: is_decomposed


        is_decomposed = .false.

        if (trim(dimname) == 'nCells'         &
            .or. trim(dimname) == 'nVertices' &
            .or. trim(dimname) == 'nEdges') then

            is_decomposed = .true.

        end if

    end function is_decomposed_dim

end module mpas_runtime_vars
