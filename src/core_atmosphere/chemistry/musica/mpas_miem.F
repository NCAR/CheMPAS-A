! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_miem
!
!> \brief Owns the MUSICA MIEM emissions lifecycle for MPAS
!> \details
!>  Every MPI rank owns a selected-cell MIEM object ordered exactly like its
!>  owned MPAS cells. The adapter validates immutable inventory grid metadata
!>  against rank-local MPAS geometry after the first file read and before any
!>  flux is returned for MICM mutation. Raw MUSICA pointers are never cached
!>  outside the emissions object because MUSICA refreshes them after every run.
!-----------------------------------------------------------------------
module mpas_miem

    use iso_fortran_env, only: int64, real64
    use musica_emissions, only: mechanism_t, emissions_t, &
        MUSICA_EMISSIONS_GRID_GEOMETRY_PLANAR, &
        MUSICA_EMISSIONS_GRID_GEOMETRY_SPHERICAL, &
        MUSICA_EMISSIONS_GRID_FIELD_AREA_CELL, &
        MUSICA_EMISSIONS_GRID_FIELD_LAT_CELL, &
        MUSICA_EMISSIONS_GRID_FIELD_LON_CELL, &
        MUSICA_EMISSIONS_GRID_FIELD_X_CELL, &
        MUSICA_EMISSIONS_GRID_FIELD_Y_CELL, &
        MUSICA_EMISSIONS_GRID_FIELD_Z_CELL
    use musica_util, only: error_t

    implicit none

    private

    integer, parameter :: MIEM_NAME_LEN = 64

    type(mechanism_t), pointer, save :: mechanism => null()
    type(emissions_t), pointer, save :: emissions => null()

    logical, save :: enabled = .false.
    logical, save :: grid_metadata_validated = .false.
    logical, save :: diagnostic_layered_output = .false.
    logical, save :: expected_on_a_sphere = .false.
    logical, save :: expected_is_periodic = .false.
    integer, save :: global_ncells = 0
    integer, save :: local_ncells = 0
    integer, save :: n_vert_levels = 0
    integer, save :: n_species = 0
    integer(int64), save :: cell_stride = 0_int64
    integer(int64), save :: species_stride = 0_int64
    integer(int64), save :: layer_cell_stride = 0_int64
    integer(int64), save :: layer_level_stride = 0_int64
    integer(int64), save :: layer_species_stride = 0_int64
    real(real64), save :: expected_sphere_radius = 0.0_real64
    real(real64), save :: first_run_seconds = -1.0_real64
    real(real64), save :: warm_run_seconds = 0.0_real64
    integer, save :: warm_run_calls = 0
    integer, allocatable, save :: species_storage_indices(:)
    logical, allocatable, save :: net_flux_species(:)
    integer, allocatable, save :: expected_global_ids(:)
    integer, allocatable, save :: diagnostic_category_ids(:)
    character(len=MIEM_NAME_LEN), allocatable, save :: species_names(:)
    character(len=MIEM_NAME_LEN), allocatable, save :: diagnostic_sector_names(:)
    character(len=MIEM_NAME_LEN), allocatable, save :: diagnostic_sector_labels(:)
    real(real64), allocatable, save :: local_cumulative_mass_kg(:)
    real(real64), allocatable, save :: expected_area_cell(:)
    real(real64), allocatable, save :: expected_lat_cell(:)
    real(real64), allocatable, save :: expected_lon_cell(:)
    real(real64), allocatable, save :: expected_x_cell(:)
    real(real64), allocatable, save :: expected_y_cell(:)
    real(real64), allocatable, save :: expected_z_cell(:)

    public :: miem_init, miem_get_species, miem_run
    public :: miem_prepare_diagnostic_selection, miem_get_diagnostic_selection
    public :: miem_get_diagnostic_fluxes
    public :: miem_commit_mass, miem_finalize, miem_is_enabled
    public :: miem_validate_owned_ids
    public :: miem_configure_net_flux_species

contains

    !-------------------------------------------------------------------
    !  routine miem_prepare_diagnostic_selection
    !
    !> \brief Parse and validate bounded runtime diagnostic selections
    !-------------------------------------------------------------------
    subroutine miem_prepare_diagnostic_selection(sector_text, category_text, &
                                                 number_of_vertical_levels, number_of_species, &
                                                 layered_diagnostics, max_diagnostic_fields, &
                                                 sectors, sector_labels, categories, &
                                                 error_code, error_message)

        character(len=*), intent(in) :: sector_text, category_text
        integer, intent(in) :: number_of_vertical_levels, number_of_species
        logical, intent(in) :: layered_diagnostics
        integer, intent(in) :: max_diagnostic_fields
        character(len=MIEM_NAME_LEN), allocatable, intent(out) :: sectors(:)
        character(len=MIEM_NAME_LEN), allocatable, intent(out) :: sector_labels(:)
        integer, allocatable, intent(out) :: categories(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: i

        error_code = 0
        error_message = ''
        call parse_sector_list(sector_text, sectors, error_code, error_message)
        if (error_code /= 0) return
        call parse_category_list(category_text, categories, error_code, error_message)
        if (error_code /= 0) return
        call validate_diagnostic_inputs(sectors, categories, number_of_vertical_levels, &
                                        number_of_species, layered_diagnostics, &
                                        max_diagnostic_fields, error_code, error_message)
        if (error_code /= 0) return

        allocate(sector_labels(size(sectors)))
        do i = 1, size(sectors)
            call sanitize_sector_label(sectors(i), sector_labels(i), error_code, error_message)
            if (error_code /= 0) return
        end do

    end subroutine miem_prepare_diagnostic_selection


    !-------------------------------------------------------------------
    !  routine miem_configure_net_flux_species
    !
    !> \brief Allow explicit signed surface exchange for named species
    !> \details Positive values remain upward sources. Negative values are
    !>  downward uptake and are accepted only for exact names in the
    !>  comma-separated configuration. The default is source-only for every
    !>  species, preserving the existing emissions safety contract.
    !-------------------------------------------------------------------
    subroutine miem_configure_net_flux_species(text, error_code, error_message)

        character(len=*), intent(in) :: text
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=MIEM_NAME_LEN), allocatable :: requested(:)
        integer :: i, j, matches

        error_code = 0
        error_message = ''
        if (allocated(net_flux_species)) deallocate(net_flux_species)
        allocate(net_flux_species(n_species))
        net_flux_species = .false.

        if (len_trim(text) == 0) return
        if (.not. enabled) then
            call set_error('config_miem_net_flux_species requires a configured MIEM file', &
                           error_code, error_message)
            return
        end if
        call parse_sector_list(text, requested, error_code, error_message)
        if (error_code /= 0) then
            error_message = 'invalid config_miem_net_flux_species: ' // trim(error_message)
            return
        end if
        do i = 1, size(requested)
            matches = 0
            do j = 1, n_species
                if (trim(requested(i)) == trim(species_names(j))) then
                    net_flux_species(j) = .true.
                    matches = matches + 1
                end if
            end do
            if (matches /= 1) then
                call set_error('config_miem_net_flux_species is not an exact MIEM species: ' // &
                               trim(requested(i)), error_code, error_message)
                return
            end if
            do j = 1, i - 1
                if (trim(requested(i)) == trim(requested(j))) then
                    call set_error('duplicate config_miem_net_flux_species entry: ' // &
                                   trim(requested(i)), error_code, error_message)
                    return
                end if
            end do
        end do

    end subroutine miem_configure_net_flux_species


    logical function allows_net_flux(species)

        integer, intent(in) :: species

        allows_net_flux = .false.
        if (.not. allocated(net_flux_species)) return
        if (species < 1 .or. species > size(net_flux_species)) return
        allows_net_flux = net_flux_species(species)

    end function allows_net_flux


    !-------------------------------------------------------------------
    !  routine parse_sector_list
    !-------------------------------------------------------------------
    subroutine parse_sector_list(text, values, error_code, error_message)

        character(len=*), intent(in) :: text
        character(len=MIEM_NAME_LEN), allocatable, intent(out) :: values(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=:), allocatable :: token
        integer :: comma, count, first, i, last, stop

        error_code = 0
        error_message = ''
        last = len_trim(text)
        if (last == 0) then
            allocate(values(0))
            return
        end if

        count = 1
        do i = 1, last
            if (text(i:i) == ',') count = count + 1
        end do
        allocate(values(count))
        values = ''

        first = 1
        do i = 1, count
            if (first > last) then
                call set_error('config_miem_diagnostic_sectors contains an empty label', &
                               error_code, error_message)
                return
            end if
            comma = index(text(first:last), ',')
            if (comma == 0) then
                stop = last
            else
                stop = first + comma - 2
            end if
            if (stop < first) then
                call set_error('config_miem_diagnostic_sectors contains an empty label', &
                               error_code, error_message)
                return
            end if
            token = trim(adjustl(text(first:stop)))
            if (len(token) == 0) then
                call set_error('config_miem_diagnostic_sectors contains an empty label', &
                               error_code, error_message)
                return
            end if
            if (len(token) > MIEM_NAME_LEN) then
                call set_error('diagnostic sector label exceeds 64 characters: ' // token, &
                               error_code, error_message)
                return
            end if
            values(i) = token
            first = stop + 2
        end do

    end subroutine parse_sector_list


    !-------------------------------------------------------------------
    !  routine parse_category_list
    !-------------------------------------------------------------------
    subroutine parse_category_list(text, values, error_code, error_message)

        character(len=*), intent(in) :: text
        integer, allocatable, intent(out) :: values(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=:), allocatable :: token
        integer :: comma, count, first, i, ios, last, start_digit, stop

        error_code = 0
        error_message = ''
        last = len_trim(text)
        if (last == 0) then
            allocate(values(0))
            return
        end if

        count = 1
        do i = 1, last
            if (text(i:i) == ',') count = count + 1
        end do
        allocate(values(count))
        values = 0

        first = 1
        do i = 1, count
            if (first > last) then
                call set_error('config_miem_diagnostic_categories contains an empty ID', &
                               error_code, error_message)
                return
            end if
            comma = index(text(first:last), ',')
            if (comma == 0) then
                stop = last
            else
                stop = first + comma - 2
            end if
            if (stop < first) then
                call set_error('config_miem_diagnostic_categories contains an empty ID', &
                               error_code, error_message)
                return
            end if
            token = trim(adjustl(text(first:stop)))
            if (len(token) == 0) then
                call set_error('config_miem_diagnostic_categories contains an empty ID', &
                               error_code, error_message)
                return
            end if
            start_digit = 1
            if (token(1:1) == '+') start_digit = 2
            if (start_digit > len(token)) then
                call set_error('invalid diagnostic category ID: ' // token, &
                               error_code, error_message)
                return
            end if
            do ios = start_digit, len(token)
                if (token(ios:ios) < '0' .or. token(ios:ios) > '9') then
                    call set_error('invalid diagnostic category ID: ' // token, &
                                   error_code, error_message)
                    return
                end if
            end do
            read(token, *, iostat=ios) values(i)
            if (ios /= 0 .or. values(i) < 0) then
                call set_error('diagnostic category IDs must be nonnegative integers: ' // token, &
                               error_code, error_message)
                return
            end if
            first = stop + 2
        end do

    end subroutine parse_category_list


    !-------------------------------------------------------------------
    !  routine validate_diagnostic_inputs
    !-------------------------------------------------------------------
    subroutine validate_diagnostic_inputs(sectors, categories, number_of_vertical_levels, &
                                          number_of_species, layered_diagnostics, &
                                          max_diagnostic_fields, error_code, error_message)

        character(len=*), intent(in) :: sectors(:)
        integer, intent(in) :: categories(:)
        integer, intent(in) :: number_of_vertical_levels, number_of_species
        logical, intent(in) :: layered_diagnostics
        integer, intent(in) :: max_diagnostic_fields
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=MIEM_NAME_LEN), allocatable :: labels(:)
        integer(int64) :: field_count, level_factor
        integer :: i, j

        error_code = 0
        error_message = ''
        if (number_of_vertical_levels < 1) then
            call set_error('diagnostic validation requires a positive vertical-level count', &
                           error_code, error_message)
            return
        end if
        if (number_of_species < 0) then
            call set_error('diagnostic validation received a negative species count', &
                           error_code, error_message)
            return
        end if
        if (max_diagnostic_fields < 1) then
            call set_error('config_miem_max_diagnostic_fields must be positive', &
                           error_code, error_message)
            return
        end if

        allocate(labels(size(sectors)))
        labels = ''
        do i = 1, size(sectors)
            if (len_trim(sectors(i)) < 1 .or. len_trim(sectors(i)) > MIEM_NAME_LEN) then
                call set_error('diagnostic sector names must contain 1..64 characters', &
                               error_code, error_message)
                return
            end if
            call sanitize_sector_label(sectors(i), labels(i), error_code, error_message)
            if (error_code /= 0) return
            do j = 1, i - 1
                if (trim(sectors(j)) == trim(sectors(i))) then
                    call set_error('duplicate diagnostic sector: ' // trim(sectors(i)), &
                                   error_code, error_message)
                    return
                end if
                if (trim(labels(j)) == trim(labels(i))) then
                    call set_error('diagnostic sector labels collide after sanitization: ' // &
                                   trim(sectors(j)) // ' and ' // trim(sectors(i)), &
                                   error_code, error_message)
                    return
                end if
            end do
        end do

        do i = 1, size(categories)
            if (categories(i) < 0) then
                call set_error('diagnostic category IDs must be nonnegative', &
                               error_code, error_message)
                return
            end if
            do j = 1, i - 1
                if (categories(j) == categories(i)) then
                    call set_error('duplicate diagnostic category: ' // &
                                   integer_string(categories(i)), error_code, error_message)
                    return
                end if
            end do
        end do

        if (number_of_species > 0) then
            level_factor = 1_int64
            if (layered_diagnostics) level_factor = int(number_of_vertical_levels, int64)
            field_count = int(number_of_species, int64) * &
                          int(size(sectors) + size(categories), int64) * level_factor
            if (field_count > int(max_diagnostic_fields, int64)) then
                call set_error('requested disaggregated emissions diagnostics require ' // &
                               int64_string(field_count) // ' fields, exceeding cap ' // &
                               integer_string(max_diagnostic_fields), error_code, error_message)
                return
            end if
        end if

    end subroutine validate_diagnostic_inputs


    !-------------------------------------------------------------------
    !  routine sanitize_sector_label
    !
    !> \brief Lowercase alphanumerics and collapse other runs to underscores
    !-------------------------------------------------------------------
    subroutine sanitize_sector_label(name, label, error_code, error_message)

        character(len=*), intent(in) :: name
        character(len=MIEM_NAME_LEN), intent(out) :: label
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character :: ch
        logical :: pending_separator
        integer :: code, i, output_length

        error_code = 0
        error_message = ''
        label = ''
        output_length = 0
        pending_separator = .false.
        do i = 1, len_trim(name)
            ch = name(i:i)
            code = iachar(ch)
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                ch = achar(code + iachar('a') - iachar('A'))
                code = iachar(ch)
            end if
            if ((code >= iachar('a') .and. code <= iachar('z')) .or. &
                (code >= iachar('0') .and. code <= iachar('9'))) then
                if (pending_separator .and. output_length > 0) then
                    output_length = output_length + 1
                    label(output_length:output_length) = '_'
                end if
                output_length = output_length + 1
                label(output_length:output_length) = ch
                pending_separator = .false.
            else
                pending_separator = .true.
            end if
        end do
        if (output_length == 0) then
            call set_error('diagnostic sector sanitizes to an empty label: ' // trim(name), &
                           error_code, error_message)
        end if

    end subroutine sanitize_sector_label


    !-------------------------------------------------------------------
    !  routine validate_returned_diagnostics
    !-------------------------------------------------------------------
    subroutine validate_returned_diagnostics(sectors, categories, layered_diagnostics, &
                                             error_code, error_message)

        character(len=*), intent(in) :: sectors(:)
        integer, intent(in) :: categories(:)
        logical, intent(in) :: layered_diagnostics
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: i

        error_code = 0
        error_message = ''
        if (.not. associated(emissions%sector_ordering) .or. &
            emissions%sector_ordering%size() /= size(sectors)) then
            call set_error('MUSICA returned an inconsistent diagnostic sector count', &
                           error_code, error_message)
            return
        end if
        do i = 1, size(sectors)
            if (trim(emissions%sector_ordering%name(i)) /= trim(sectors(i))) then
                call set_error('MUSICA changed diagnostic sector ordering at position ' // &
                               integer_string(i), error_code, error_message)
                return
            end if
        end do
        if (size(categories) == 0) then
            if (associated(emissions%category_ids)) then
                if (size(emissions%category_ids) /= 0) then
                    call set_error('MUSICA returned unexpected diagnostic categories', &
                                   error_code, error_message)
                    return
                end if
            end if
        else
            if (.not. associated(emissions%category_ids) .or. &
                size(emissions%category_ids) /= size(categories)) then
                call set_error('MUSICA returned an inconsistent diagnostic category count', &
                               error_code, error_message)
                return
            end if
            do i = 1, size(categories)
                if (emissions%category_ids(i) /= categories(i)) then
                    call set_error('MUSICA changed diagnostic category ordering at position ' // &
                                   integer_string(i), error_code, error_message)
                    return
                end if
            end do
        end if
        if (emissions%layered_diagnostics .neqv. layered_diagnostics) then
            call set_error('MUSICA changed the layered-diagnostics selection', &
                           error_code, error_message)
        end if

    end subroutine validate_returned_diagnostics

    !-------------------------------------------------------------------
    !  routine miem_init
    !
    !> \brief Construct rank-local MIEM state in owned MPAS cell order
    !-------------------------------------------------------------------
    subroutine miem_init(config_path, number_of_global_cells, number_of_vertical_levels, &
                         selected_global_cell_ids, on_a_sphere, is_periodic, sphere_radius, &
                         area_cell, lat_cell, lon_cell, x_cell, y_cell, z_cell, &
                         error_code, error_message, diagnostic_sectors, &
                         diagnostic_categories, layered_diagnostics, max_diagnostic_fields)

        use mpas_log, only: mpas_log_write
        use mpas_kind_types, only: RKIND

        character(len=*), intent(in) :: config_path
        integer, optional, intent(in) :: number_of_global_cells
        integer, optional, intent(in) :: number_of_vertical_levels
        integer, optional, intent(in) :: selected_global_cell_ids(:)
        logical, optional, intent(in) :: on_a_sphere, is_periodic
        real(kind=RKIND), optional, intent(in) :: sphere_radius
        real(kind=RKIND), optional, intent(in) :: area_cell(:)
        real(kind=RKIND), optional, intent(in) :: lat_cell(:), lon_cell(:)
        real(kind=RKIND), optional, intent(in) :: x_cell(:), y_cell(:), z_cell(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        character(len=*), optional, intent(in) :: diagnostic_sectors(:)
        integer, optional, intent(in) :: diagnostic_categories(:)
        logical, optional, intent(in) :: layered_diagnostics
        integer, optional, intent(in) :: max_diagnostic_fields

        type(error_t) :: error
        character(len=:), allocatable :: name
        integer :: i, j

        call miem_release()
        error_code = 0
        error_message = ''

        if (len_trim(config_path) == 0) then
            call mpas_log_write('[MIEM] Disabled (config_miem_file is empty).')
            return
        end if

        if (.not. present(number_of_global_cells) .or. &
            .not. present(number_of_vertical_levels) .or. &
            .not. present(selected_global_cell_ids) .or. &
            .not. present(on_a_sphere) .or. .not. present(is_periodic) .or. &
            .not. present(sphere_radius) .or. .not. present(area_cell) .or. &
            .not. present(lat_cell) .or. .not. present(lon_cell) .or. &
            .not. present(x_cell) .or. .not. present(y_cell) .or. .not. present(z_cell)) then
            call set_error('enabled selected-cell construction requires complete MPAS mesh metadata', &
                           error_code, error_message)
            return
        end if

        if (number_of_global_cells < 1) then
            call set_error('global cell count must be positive', error_code, error_message)
            return
        end if
        if (number_of_vertical_levels < 1) then
            call set_error('vertical level count must be positive', error_code, error_message)
            return
        end if
        if (size(selected_global_cell_ids) < 1) then
            call set_error('selected-cell construction requires at least one owned cell', &
                           error_code, error_message)
            return
        end if
        if (size(area_cell) < size(selected_global_cell_ids) .or. &
            size(lat_cell) < size(selected_global_cell_ids) .or. &
            size(lon_cell) < size(selected_global_cell_ids) .or. &
            size(x_cell) < size(selected_global_cell_ids) .or. &
            size(y_cell) < size(selected_global_cell_ids) .or. &
            size(z_cell) < size(selected_global_cell_ids)) then
            call set_error('rank-local MPAS mesh arrays are smaller than the owned-cell selection', &
                           error_code, error_message)
            return
        end if

        call validate_selected_ids(selected_global_cell_ids, number_of_global_cells, &
                                   error_code, error_message)
        if (error_code /= 0) return

        if (present(diagnostic_sectors) .or. present(diagnostic_categories) .or. &
            present(layered_diagnostics) .or. present(max_diagnostic_fields)) then
            if (.not. present(diagnostic_sectors) .or. &
                .not. present(diagnostic_categories) .or. &
                .not. present(layered_diagnostics) .or. &
                .not. present(max_diagnostic_fields)) then
                call set_error('diagnostic construction options must be supplied together', &
                               error_code, error_message)
                return
            end if
            call validate_diagnostic_inputs(diagnostic_sectors, diagnostic_categories, &
                                            number_of_vertical_levels, 0, layered_diagnostics, &
                                            max_diagnostic_fields, error_code, error_message)
            if (error_code /= 0) return
        end if

        mechanism => mechanism_t(trim(config_path), error)
        if (.not. error%is_success()) then
            call set_musica_error('could not read emissions configuration', error, &
                                  error_code, error_message)
            call miem_release()
            return
        end if
        if (.not. associated(mechanism)) then
            call set_error('MUSICA returned no mechanism object', error_code, error_message)
            call miem_release()
            return
        end if

        if (present(diagnostic_sectors)) then
            emissions => emissions_t(mechanism, number_of_global_cells, &
                                     number_of_vertical_levels, selected_global_cell_ids, error, &
                                     diagnostic_sectors=diagnostic_sectors, &
                                     diagnostic_categories=diagnostic_categories, &
                                     layered_diagnostics=layered_diagnostics, &
                                     max_diagnostic_fields=max_diagnostic_fields)
        else
            emissions => emissions_t(mechanism, number_of_global_cells, &
                                     number_of_vertical_levels, selected_global_cell_ids, error)
        end if
        if (.not. error%is_success()) then
            call set_musica_error('could not construct emissions object', error, &
                                  error_code, error_message)
            call miem_release()
            return
        end if
        if (.not. associated(emissions)) then
            call set_error('MUSICA returned no emissions object', error_code, error_message)
            call miem_release()
            return
        end if
        if (.not. associated(emissions%species_ordering)) then
            call set_error('MUSICA returned no emissions species ordering', &
                           error_code, error_message)
            call miem_release()
            return
        end if

        if (emissions%global_number_of_cells /= number_of_global_cells .or. &
            emissions%number_of_cells /= size(selected_global_cell_ids) .or. &
            emissions%number_of_vertical_levels /= number_of_vertical_levels) then
            call set_error('MUSICA returned inconsistent selected-grid dimensions', &
                           error_code, error_message)
            call miem_release()
            return
        end if
        if (.not. associated(emissions%selected_global_cell_ids)) then
            call set_error('MUSICA returned no selected global cell IDs', &
                           error_code, error_message)
            call miem_release()
            return
        end if
        if (size(emissions%selected_global_cell_ids) /= size(selected_global_cell_ids)) then
            call set_error('MUSICA returned an inconsistent selected-cell ID count', &
                           error_code, error_message)
            call miem_release()
            return
        end if
        do i = 1, size(selected_global_cell_ids)
            if (int(emissions%selected_global_cell_ids(i), int64) /= &
                int(selected_global_cell_ids(i), int64)) then
                call set_error('MUSICA changed selected-cell ordering at local cell ' // &
                               integer_string(i), error_code, error_message)
                call miem_release()
                return
            end if
        end do

        n_species = emissions%number_of_species
        if (n_species < 1 .or. emissions%species_ordering%size() /= n_species) then
            call set_error('emissions species ordering is empty or inconsistent', &
                           error_code, error_message)
            call miem_release()
            return
        end if

        if (present(diagnostic_sectors)) then
            call validate_diagnostic_inputs(diagnostic_sectors, diagnostic_categories, &
                                            number_of_vertical_levels, n_species, &
                                            layered_diagnostics, max_diagnostic_fields, &
                                            error_code, error_message)
            if (error_code /= 0) then
                call miem_release()
                return
            end if
            call validate_returned_diagnostics(diagnostic_sectors, diagnostic_categories, &
                                               layered_diagnostics, error_code, error_message)
            if (error_code /= 0) then
                call miem_release()
                return
            end if
            diagnostic_layered_output = layered_diagnostics
            allocate(diagnostic_sector_names(size(diagnostic_sectors)))
            allocate(diagnostic_sector_labels(size(diagnostic_sectors)))
            allocate(diagnostic_category_ids(size(diagnostic_categories)))
            diagnostic_sector_names = diagnostic_sectors
            diagnostic_category_ids = diagnostic_categories
            do i = 1, size(diagnostic_sectors)
                call sanitize_sector_label(diagnostic_sectors(i), diagnostic_sector_labels(i), &
                                           error_code, error_message)
                if (error_code /= 0) then
                    call miem_release()
                    return
                end if
            end do
        else
            allocate(diagnostic_sector_names(0), diagnostic_sector_labels(0))
            allocate(diagnostic_category_ids(0))
        end if

        cell_stride = int(emissions%cell_stride, int64)
        species_stride = int(emissions%species_stride, int64)
        layer_cell_stride = int(emissions%layer_cell_stride, int64)
        layer_level_stride = int(emissions%layer_level_stride, int64)
        layer_species_stride = int(emissions%layer_species_stride, int64)
        if (cell_stride < 1_int64 .or. species_stride < 1_int64 .or. &
            layer_cell_stride < 1_int64 .or. layer_level_stride < 1_int64 .or. &
            layer_species_stride < 1_int64) then
            call set_error('MUSICA returned invalid emissions flux strides', &
                           error_code, error_message)
            call miem_release()
            return
        end if

        allocate(species_names(n_species))
        allocate(species_storage_indices(n_species))
        allocate(local_cumulative_mass_kg(n_species))
        species_names = ''
        species_storage_indices = 0
        local_cumulative_mass_kg = 0.0_real64

        do i = 1, n_species
            name = emissions%species_ordering%name(i)
            if (len_trim(name) < 1) then
                call set_error('MUSICA returned an empty emissions species name', &
                               error_code, error_message)
                call miem_release()
                return
            end if
            if (len_trim(name) > MIEM_NAME_LEN) then
                call set_error('emissions species name exceeds 64 characters: ' // trim(name), &
                               error_code, error_message)
                call miem_release()
                return
            end if
            do j = 1, i - 1
                if (trim(species_names(j)) == trim(name)) then
                    call set_error('duplicate emissions species name: ' // trim(name), &
                                   error_code, error_message)
                    call miem_release()
                    return
                end if
            end do

            species_names(i) = trim(name)
            species_storage_indices(i) = &
                emissions%species_ordering%index(trim(name), error)
            if (.not. error%is_success()) then
                call set_musica_error('could not resolve emissions species ' // trim(name), &
                                      error, error_code, error_message)
                call miem_release()
                return
            end if
            if (species_storage_indices(i) < 1) then
                call set_error('invalid storage index for emissions species ' // trim(name), &
                               error_code, error_message)
                call miem_release()
                return
            end if
        end do

        global_ncells = number_of_global_cells
        local_ncells = size(selected_global_cell_ids)
        n_vert_levels = number_of_vertical_levels
        expected_on_a_sphere = on_a_sphere
        expected_is_periodic = is_periodic
        expected_sphere_radius = real(sphere_radius, real64)
        allocate(expected_global_ids(local_ncells), expected_area_cell(local_ncells))
        allocate(expected_lat_cell(local_ncells), expected_lon_cell(local_ncells))
        allocate(expected_x_cell(local_ncells), expected_y_cell(local_ncells))
        allocate(expected_z_cell(local_ncells))
        expected_global_ids = selected_global_cell_ids
        expected_area_cell = real(area_cell(1:local_ncells), real64)
        expected_lat_cell = real(lat_cell(1:local_ncells), real64)
        expected_lon_cell = real(lon_cell(1:local_ncells), real64)
        expected_x_cell = real(x_cell(1:local_ncells), real64)
        expected_y_cell = real(y_cell(1:local_ncells), real64)
        expected_z_cell = real(z_cell(1:local_ncells), real64)
        enabled = .true.

        call mpas_log_write('[MIEM] Initialized selected-cell emissions object: $i of $i cells, $i levels, $i species.', &
                            intArgs=[local_ncells, global_ncells, n_vert_levels, n_species])
        do i = 1, n_species
            call mpas_log_write('[MIEM] Emissions species: ' // trim(species_names(i)))
        end do
        if (size(diagnostic_sector_names) + size(diagnostic_category_ids) > 0) then
            call mpas_log_write('[MIEM] Requested $i sector and $i category diagnostics.', &
                                intArgs=[size(diagnostic_sector_names), size(diagnostic_category_ids)])
            if (diagnostic_layered_output) then
                call mpas_log_write('[MIEM] Requested group and total layered emissions diagnostics.')
            end if
        else if (diagnostic_layered_output) then
            call mpas_log_write('[MIEM] Requested total layered emissions diagnostics only.')
        end if

    end subroutine miem_init


    !-------------------------------------------------------------------
    !  routine miem_get_species
    !
    !> \brief Return a copy of the configured MIEM species names
    !-------------------------------------------------------------------
    subroutine miem_get_species(number_of_species, names)

        integer, intent(out) :: number_of_species
        character(len=MIEM_NAME_LEN), allocatable, intent(out) :: names(:)

        number_of_species = n_species
        allocate(names(n_species))
        if (n_species > 0) names = species_names

    end subroutine miem_get_species


    !-------------------------------------------------------------------
    !  routine miem_get_diagnostic_selection
    !
    !> \brief Return immutable requested group names and output labels
    !-------------------------------------------------------------------
    subroutine miem_get_diagnostic_selection(number_of_sectors, sectors, sector_labels, &
                                             number_of_categories, categories, layered_diagnostics)

        integer, intent(out) :: number_of_sectors, number_of_categories
        character(len=MIEM_NAME_LEN), allocatable, intent(out) :: sectors(:)
        character(len=MIEM_NAME_LEN), allocatable, intent(out) :: sector_labels(:)
        integer, allocatable, intent(out) :: categories(:)
        logical, intent(out) :: layered_diagnostics

        number_of_sectors = 0
        number_of_categories = 0
        if (allocated(diagnostic_sector_names)) number_of_sectors = size(diagnostic_sector_names)
        if (allocated(diagnostic_category_ids)) number_of_categories = size(diagnostic_category_ids)
        allocate(sectors(number_of_sectors), sector_labels(number_of_sectors))
        allocate(categories(number_of_categories))
        if (number_of_sectors > 0) then
            sectors = diagnostic_sector_names
            sector_labels = diagnostic_sector_labels
        end if
        if (number_of_categories > 0) categories = diagnostic_category_ids
        layered_diagnostics = diagnostic_layered_output

    end subroutine miem_get_diagnostic_selection


    !-------------------------------------------------------------------
    !  routine miem_validate_owned_ids
    !
    !> \brief Validate rank-owned global IDs and return their local maximum
    !-------------------------------------------------------------------
    subroutine miem_validate_owned_ids(index_to_cell_id, ncells_solve, local_max_id, &
                                       error_code, error_message)

        integer, intent(in) :: index_to_cell_id(:)
        integer, intent(in) :: ncells_solve
        integer, intent(out) :: local_max_id
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer, allocatable :: sorted_ids(:)
        integer :: i

        error_code = 0
        error_message = ''
        local_max_id = 0

        if (ncells_solve < 0 .or. ncells_solve > size(index_to_cell_id)) then
            call set_error('owned-cell count is outside indexToCellID bounds', &
                           error_code, error_message)
            return
        end if
        if (ncells_solve == 0) return

        allocate(sorted_ids(ncells_solve))
        sorted_ids = index_to_cell_id(1:ncells_solve)
        do i = 1, ncells_solve
            if (sorted_ids(i) < 1) then
                call set_error('owned indexToCellID values must be positive', &
                               error_code, error_message)
                deallocate(sorted_ids)
                return
            end if
        end do

        call heap_sort(sorted_ids)
        do i = 2, ncells_solve
            if (sorted_ids(i) == sorted_ids(i - 1)) then
                call set_error('duplicate owned indexToCellID value: ' // &
                               integer_string(sorted_ids(i)), error_code, error_message)
                deallocate(sorted_ids)
                return
            end if
        end do

        local_max_id = sorted_ids(ncells_solve)
        deallocate(sorted_ids)

    end subroutine miem_validate_owned_ids


    !-------------------------------------------------------------------
    !  routine validate_selected_ids
    !-------------------------------------------------------------------
    subroutine validate_selected_ids(ids, number_of_global_cells, error_code, error_message)

        integer, intent(in) :: ids(:)
        integer, intent(in) :: number_of_global_cells
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: local_max_id

        call miem_validate_owned_ids(ids, size(ids), local_max_id, error_code, error_message)
        if (error_code /= 0) return
        if (local_max_id > number_of_global_cells) then
            call set_error('selected global cell ID is outside 1..global_nCells: ' // &
                           integer_string(local_max_id), error_code, error_message)
        end if

    end subroutine validate_selected_ids


    !-------------------------------------------------------------------
    !  routine miem_run
    !
    !> \brief Run MIEM and return column/layer fluxes in owned-cell order
    !-------------------------------------------------------------------
    subroutine miem_run(epoch_seconds, dt_seconds, local_flux, local_layer_flux, &
                        error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use mpas_dmpar, only: mpas_dmpar_get_time
        real(real64), intent(in) :: epoch_seconds
        real(real64), intent(in) :: dt_seconds
        real(real64), allocatable, intent(out) :: local_flux(:,:)
        real(real64), allocatable, intent(out) :: local_layer_flux(:,:,:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: error
        integer(int64) :: buffer_index, layer_buffer_index
        real(real64) :: column_sum, closure_scale, run_start, run_stop
        integer :: i, k, s

        error_code = 0
        error_message = ''

        if (.not. enabled) then
            allocate(local_flux(0, 0), local_layer_flux(0, 0, 0))
            return
        end if
        if (.not. associated(emissions)) then
            call set_error('enabled adapter has no emissions object', error_code, error_message)
            return
        end if
        if (.not. ieee_is_finite(epoch_seconds)) then
            call set_error('epoch seconds must be finite', error_code, error_message)
            return
        end if
        if (.not. ieee_is_finite(dt_seconds) .or. dt_seconds <= 0.0_real64) then
            call set_error('emissions timestep must be finite and positive', &
                           error_code, error_message)
            return
        end if
        allocate(local_flux(n_species, local_ncells))
        allocate(local_layer_flux(n_species, n_vert_levels, local_ncells))
        local_flux = 0.0_real64
        local_layer_flux = 0.0_real64

        call mpas_dmpar_get_time(run_start)
        call emissions%run(epoch_seconds, dt_seconds, error)
        call mpas_dmpar_get_time(run_stop)
        if (.not. error%is_success()) then
            call set_musica_error('emissions run failed', error, error_code, error_message)
            return
        end if
        if (first_run_seconds < 0.0_real64) then
            first_run_seconds = max(0.0_real64, run_stop - run_start)
        else
            warm_run_seconds = warm_run_seconds + max(0.0_real64, run_stop - run_start)
            warm_run_calls = warm_run_calls + 1
        end if
        if (.not. associated(emissions%surface_flux)) then
            call set_error('MUSICA returned no surface-flux buffer after run', &
                           error_code, error_message)
            return
        end if
        if (.not. associated(emissions%layer_flux)) then
            call set_error('MUSICA returned no layer-flux buffer after run', &
                           error_code, error_message)
            return
        end if

        if (.not. grid_metadata_validated) then
            call validate_grid_metadata(error_code, error_message)
            if (error_code /= 0) return
            grid_metadata_validated = .true.
        end if

        do s = 1, n_species
            do i = 1, local_ncells
                buffer_index = 1_int64 &
                    + int(i - 1, int64) * cell_stride &
                    + int(species_storage_indices(s) - 1, int64) * species_stride
                if (buffer_index < 1_int64 .or. &
                    buffer_index > int(size(emissions%surface_flux), int64)) then
                    call set_error('surface-flux buffer index is out of bounds for species ' // &
                                   trim(species_names(s)), error_code, error_message)
                    return
                end if
                local_flux(s, i) = emissions%surface_flux(buffer_index)
                if (.not. ieee_is_finite(local_flux(s, i)) .or. &
                    (local_flux(s, i) < 0.0_real64 .and. .not. allows_net_flux(s))) then
                    call set_error('invalid surface flux for species ' // &
                                   trim(species_names(s)) // ', global cell ' // &
                                   integer_string(expected_global_ids(i)), error_code, error_message)
                    return
                end if

                do k = 1, n_vert_levels
                    layer_buffer_index = 1_int64 &
                        + int(i - 1, int64) * layer_cell_stride &
                        + int(k - 1, int64) * layer_level_stride &
                        + int(species_storage_indices(s) - 1, int64) * layer_species_stride
                    if (layer_buffer_index < 1_int64 .or. &
                        layer_buffer_index > int(size(emissions%layer_flux), int64)) then
                        call set_error('layer-flux buffer index is out of bounds for species ' // &
                                       trim(species_names(s)), error_code, error_message)
                        return
                    end if
                    local_layer_flux(s, k, i) = emissions%layer_flux(layer_buffer_index)
                    if (.not. ieee_is_finite(local_layer_flux(s, k, i)) .or. &
                        (local_layer_flux(s, k, i) < 0.0_real64 .and. &
                         .not. allows_net_flux(s))) then
                        call set_error('invalid layer flux for species ' // trim(species_names(s)) // &
                                       ', level ' // integer_string(k) // ', global cell ' // &
                                       integer_string(expected_global_ids(i)), error_code, error_message)
                        return
                    end if
                end do

                column_sum = sum(local_layer_flux(s, :, i))
                closure_scale = max(abs(local_flux(s, i)), tiny(1.0_real64))
                if (abs(column_sum - local_flux(s, i)) > &
                    64.0_real64 * epsilon(1.0_real64) * closure_scale) then
                    call set_error('layer flux does not close to column flux for species ' // &
                                   trim(species_names(s)) // ', global cell ' // &
                                   integer_string(expected_global_ids(i)), error_code, error_message)
                    return
                end if
            end do
        end do

    end subroutine miem_run


    !-------------------------------------------------------------------
    !  routine miem_get_diagnostic_fluxes
    !
    !> \brief Copy explicitly requested group buffers into owned-cell order
    !-------------------------------------------------------------------
    subroutine miem_get_diagnostic_fluxes(sector_flux, category_flux, &
                                          sector_layer_flux, category_layer_flux, &
                                          error_code, error_message)

        real(real64), allocatable, intent(out) :: sector_flux(:,:,:)
        real(real64), allocatable, intent(out) :: category_flux(:,:,:)
        real(real64), allocatable, intent(out) :: sector_layer_flux(:,:,:,:)
        real(real64), allocatable, intent(out) :: category_layer_flux(:,:,:,:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: group

        error_code = 0
        error_message = ''
        if (.not. enabled) then
            allocate(sector_flux(0, 0, 0), category_flux(0, 0, 0))
            allocate(sector_layer_flux(0, 0, 0, 0), category_layer_flux(0, 0, 0, 0))
            return
        end if
        if (.not. associated(emissions)) then
            call set_error('diagnostic extraction has no emissions object', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(diagnostic_sector_names) .or. &
            .not. allocated(diagnostic_category_ids)) then
            call set_error('diagnostic selection state is unavailable', error_code, error_message)
            return
        end if
        if (.not. allocated(emissions%sector_flux) .or. &
            size(emissions%sector_flux) /= size(diagnostic_sector_names) .or. &
            .not. allocated(emissions%category_flux) .or. &
            size(emissions%category_flux) /= size(diagnostic_category_ids)) then
            call set_error('MUSICA diagnostic column-buffer counts changed after construction', &
                           error_code, error_message)
            return
        end if

        allocate(sector_flux(n_species, size(diagnostic_sector_names), local_ncells))
        allocate(category_flux(n_species, size(diagnostic_category_ids), local_ncells))
        sector_flux = 0.0_real64
        category_flux = 0.0_real64

        if (diagnostic_layered_output) then
            if (.not. allocated(emissions%sector_layer_flux) .or. &
                size(emissions%sector_layer_flux) /= size(diagnostic_sector_names) .or. &
                .not. allocated(emissions%category_layer_flux) .or. &
                size(emissions%category_layer_flux) /= size(diagnostic_category_ids)) then
                call set_error('MUSICA diagnostic layer-buffer counts changed after construction', &
                               error_code, error_message)
                return
            end if
            allocate(sector_layer_flux(n_species, n_vert_levels, &
                                       size(diagnostic_sector_names), local_ncells))
            allocate(category_layer_flux(n_species, n_vert_levels, &
                                         size(diagnostic_category_ids), local_ncells))
            sector_layer_flux = 0.0_real64
            category_layer_flux = 0.0_real64
        else
            allocate(sector_layer_flux(0, 0, 0, 0), category_layer_flux(0, 0, 0, 0))
        end if

        do group = 1, size(diagnostic_sector_names)
            call copy_column_diagnostic_buffer('sector ' // trim(diagnostic_sector_names(group)), &
                                               emissions%sector_flux(group)%values, &
                                               sector_flux(:, group, :), error_code, error_message)
            if (error_code /= 0) return
            if (diagnostic_layered_output) then
                call copy_layer_diagnostic_buffer('sector ' // trim(diagnostic_sector_names(group)), &
                                                  emissions%sector_layer_flux(group)%values, &
                                                  sector_flux(:, group, :), &
                                                  sector_layer_flux(:, :, group, :), &
                                                  error_code, error_message)
                if (error_code /= 0) return
            end if
        end do
        do group = 1, size(diagnostic_category_ids)
            call copy_column_diagnostic_buffer('category ' // &
                                               integer_string(diagnostic_category_ids(group)), &
                                               emissions%category_flux(group)%values, &
                                               category_flux(:, group, :), error_code, error_message)
            if (error_code /= 0) return
            if (diagnostic_layered_output) then
                call copy_layer_diagnostic_buffer('category ' // &
                                                  integer_string(diagnostic_category_ids(group)), &
                                                  emissions%category_layer_flux(group)%values, &
                                                  category_flux(:, group, :), &
                                                  category_layer_flux(:, :, group, :), &
                                                  error_code, error_message)
                if (error_code /= 0) return
            end if
        end do

    end subroutine miem_get_diagnostic_fluxes


    !-------------------------------------------------------------------
    !  routine copy_column_diagnostic_buffer
    !-------------------------------------------------------------------
    subroutine copy_column_diagnostic_buffer(group_name, buffer, output, &
                                             error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        character(len=*), intent(in) :: group_name
        real(real64), pointer, intent(in) :: buffer(:)
        real(real64), intent(out) :: output(:,:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer(int64) :: buffer_index
        integer :: cell, species

        error_code = 0
        error_message = ''
        if (.not. associated(buffer)) then
            call set_error('MUSICA returned no column buffer for diagnostic ' // &
                           trim(group_name), error_code, error_message)
            return
        end if
        if (size(output, 1) /= n_species .or. size(output, 2) /= local_ncells) then
            call set_error('internal column diagnostic output has inconsistent dimensions', &
                           error_code, error_message)
            return
        end if
        do species = 1, n_species
            do cell = 1, local_ncells
                buffer_index = 1_int64 &
                    + int(cell - 1, int64) * cell_stride &
                    + int(species_storage_indices(species) - 1, int64) * species_stride
                if (buffer_index < 1_int64 .or. buffer_index > int(size(buffer), int64)) then
                    call set_error('column buffer index is out of bounds for diagnostic ' // &
                                   trim(group_name), error_code, error_message)
                    return
                end if
                output(species, cell) = buffer(buffer_index)
                if (.not. ieee_is_finite(output(species, cell)) .or. &
                    (output(species, cell) < 0.0_real64 .and. &
                     .not. allows_net_flux(species))) then
                    call set_error('invalid column flux for diagnostic ' // trim(group_name) // &
                                   ', species ' // trim(species_names(species)) // &
                                   ', global cell ' // integer_string(expected_global_ids(cell)), &
                                   error_code, error_message)
                    return
                end if
            end do
        end do

    end subroutine copy_column_diagnostic_buffer


    !-------------------------------------------------------------------
    !  routine copy_layer_diagnostic_buffer
    !-------------------------------------------------------------------
    subroutine copy_layer_diagnostic_buffer(group_name, buffer, column_flux, output, &
                                            error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        character(len=*), intent(in) :: group_name
        real(real64), pointer, intent(in) :: buffer(:)
        real(real64), intent(in) :: column_flux(:,:)
        real(real64), intent(out) :: output(:,:,:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer(int64) :: buffer_index
        real(real64) :: closure_scale
        integer :: cell, level, species

        error_code = 0
        error_message = ''
        if (.not. associated(buffer)) then
            call set_error('MUSICA returned no layer buffer for diagnostic ' // trim(group_name), &
                           error_code, error_message)
            return
        end if
        if (size(output, 1) /= n_species .or. size(output, 2) /= n_vert_levels .or. &
            size(output, 3) /= local_ncells) then
            call set_error('internal layer diagnostic output has inconsistent dimensions', &
                           error_code, error_message)
            return
        end if
        do species = 1, n_species
            do cell = 1, local_ncells
                do level = 1, n_vert_levels
                    buffer_index = 1_int64 &
                        + int(cell - 1, int64) * layer_cell_stride &
                        + int(level - 1, int64) * layer_level_stride &
                        + int(species_storage_indices(species) - 1, int64) * layer_species_stride
                    if (buffer_index < 1_int64 .or. buffer_index > int(size(buffer), int64)) then
                        call set_error('layer buffer index is out of bounds for diagnostic ' // &
                                       trim(group_name), error_code, error_message)
                        return
                    end if
                    output(species, level, cell) = buffer(buffer_index)
                    if (.not. ieee_is_finite(output(species, level, cell)) .or. &
                        (output(species, level, cell) < 0.0_real64 .and. &
                         .not. allows_net_flux(species))) then
                        call set_error('invalid layer flux for diagnostic ' // trim(group_name) // &
                                       ', species ' // trim(species_names(species)) // &
                                       ', global cell ' // integer_string(expected_global_ids(cell)), &
                                       error_code, error_message)
                        return
                    end if
                end do
                closure_scale = max(abs(column_flux(species, cell)), tiny(1.0_real64))
                if (abs(sum(output(species, :, cell)) - column_flux(species, cell)) > &
                    64.0_real64 * epsilon(1.0_real64) * closure_scale) then
                    call set_error('layer flux does not close for diagnostic ' // trim(group_name) // &
                                   ', species ' // trim(species_names(species)) // &
                                   ', global cell ' // integer_string(expected_global_ids(cell)), &
                                   error_code, error_message)
                    return
                end if
            end do
        end do

    end subroutine copy_layer_diagnostic_buffer


    !-------------------------------------------------------------------
    !  routine validate_grid_metadata
    !
    !> \brief Compare selected inventory identity with rank-local MPAS mesh
    !-------------------------------------------------------------------
    subroutine validate_grid_metadata(error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use mpas_log, only: mpas_log_write

        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=3) :: expected_on_sphere, expected_periodic
        integer :: i, expected_geometry

        error_code = 0
        error_message = ''

        if (.not. associated(emissions)) then
            call set_error('grid metadata validation has no emissions object', &
                           error_code, error_message)
            return
        end if
        if (.not. emissions%grid_metadata%available .or. &
            .not. emissions%grid_metadata%exact_grid) then
            call set_error('selected inventory does not expose exact-grid metadata', &
                           error_code, error_message)
            return
        end if
        if (emissions%grid_metadata%global_number_of_cells /= global_ncells) then
            call set_error('inventory global cell count differs from the MPAS mesh', &
                           error_code, error_message)
            return
        end if

        if (expected_on_a_sphere) then
            expected_geometry = MUSICA_EMISSIONS_GRID_GEOMETRY_SPHERICAL
            expected_on_sphere = 'YES'
        else
            expected_geometry = MUSICA_EMISSIONS_GRID_GEOMETRY_PLANAR
            expected_on_sphere = 'NO '
        end if
        if (expected_is_periodic) then
            expected_periodic = 'YES'
        else
            expected_periodic = 'NO '
        end if

        if (emissions%grid_metadata%geometry /= expected_geometry) then
            call set_error('inventory geometry class differs from the MPAS mesh', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(emissions%grid_metadata%on_a_sphere) .or. &
            trim(emissions%grid_metadata%on_a_sphere) /= trim(expected_on_sphere)) then
            call set_error('inventory on_a_sphere metadata differs from the MPAS mesh', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(emissions%grid_metadata%is_periodic) .or. &
            trim(emissions%grid_metadata%is_periodic) /= trim(expected_periodic)) then
            call set_error('inventory is_periodic metadata differs from the MPAS mesh', &
                           error_code, error_message)
            return
        end if
        if (emissions%grid_metadata%has_sphere_radius) then
            if (.not. ieee_is_finite(emissions%grid_metadata%sphere_radius) .or. &
                .not. scalar_values_match(emissions%grid_metadata%sphere_radius, &
                                          expected_sphere_radius)) then
                call set_error('inventory sphere_radius differs from the MPAS mesh', &
                               error_code, error_message)
                return
            end if
        end if
        if (.not. allocated(emissions%grid_metadata%fingerprint_algorithm) .or. &
            trim(emissions%grid_metadata%fingerprint_algorithm) /= &
            'chempas-mesh-sha256-v1') then
            call set_error('inventory mesh fingerprint algorithm is not chempas-mesh-sha256-v1', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(emissions%grid_metadata%fingerprint) .or. &
            len_trim(emissions%grid_metadata%fingerprint) /= 64) then
            call set_error('inventory mesh fingerprint digest is missing or malformed', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(emissions%grid_metadata%field_manifest) .or. &
            len_trim(emissions%grid_metadata%field_manifest) < 1) then
            call set_error('inventory mesh field manifest is missing', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(emissions%grid_metadata%index_to_cell_id_units) .or. &
            trim(emissions%grid_metadata%index_to_cell_id_units) /= '-') then
            call set_error('inventory indexToCellID units must be -', &
                           error_code, error_message)
            return
        end if
        if (.not. associated(emissions%grid_metadata%index_to_cell_id) .or. &
            size(emissions%grid_metadata%index_to_cell_id) /= local_ncells) then
            call set_error('inventory selected indexToCellID metadata has the wrong size', &
                           error_code, error_message)
            return
        end if
        do i = 1, local_ncells
            if (int(emissions%grid_metadata%index_to_cell_id(i), int64) /= &
                int(expected_global_ids(i), int64)) then
                call set_error('inventory selected indexToCellID differs at local cell ' // &
                               integer_string(i), error_code, error_message)
                return
            end if
        end do

        call compare_grid_field('areaCell', MUSICA_EMISSIONS_GRID_FIELD_AREA_CELL, &
                                emissions%grid_metadata%area_cell, &
                                emissions%grid_metadata%area_cell_units, &
                                expected_area_cell, 'm^2', .true., error_code, error_message)
        if (error_code /= 0) return
        call compare_grid_field('latCell', MUSICA_EMISSIONS_GRID_FIELD_LAT_CELL, &
                                emissions%grid_metadata%lat_cell, &
                                emissions%grid_metadata%lat_cell_units, &
                                expected_lat_cell, 'rad', expected_on_a_sphere, &
                                error_code, error_message)
        if (error_code /= 0) return
        call compare_grid_field('lonCell', MUSICA_EMISSIONS_GRID_FIELD_LON_CELL, &
                                emissions%grid_metadata%lon_cell, &
                                emissions%grid_metadata%lon_cell_units, &
                                expected_lon_cell, 'rad', expected_on_a_sphere, &
                                error_code, error_message)
        if (error_code /= 0) return
        call compare_grid_field('xCell', MUSICA_EMISSIONS_GRID_FIELD_X_CELL, &
                                emissions%grid_metadata%x_cell, &
                                emissions%grid_metadata%x_cell_units, &
                                expected_x_cell, 'm', .not. expected_on_a_sphere, &
                                error_code, error_message)
        if (error_code /= 0) return
        call compare_grid_field('yCell', MUSICA_EMISSIONS_GRID_FIELD_Y_CELL, &
                                emissions%grid_metadata%y_cell, &
                                emissions%grid_metadata%y_cell_units, &
                                expected_y_cell, 'm', .not. expected_on_a_sphere, &
                                error_code, error_message)
        if (error_code /= 0) return
        call compare_grid_field('zCell', MUSICA_EMISSIONS_GRID_FIELD_Z_CELL, &
                                emissions%grid_metadata%z_cell, &
                                emissions%grid_metadata%z_cell_units, &
                                expected_z_cell, 'm', .not. expected_on_a_sphere, &
                                error_code, error_message)
        if (error_code /= 0) return

        call mpas_log_write('[MIEM] Validated selected exact-grid metadata: ' // &
                            trim(emissions%grid_metadata%fingerprint))

    end subroutine validate_grid_metadata


    !-------------------------------------------------------------------
    !  routine compare_grid_field
    !-------------------------------------------------------------------
    subroutine compare_grid_field(name, field, observed, observed_units, expected, &
                                  expected_units, required, error_code, error_message)

        character(len=*), intent(in) :: name, expected_units
        integer, intent(in) :: field
        real(real64), pointer, intent(in) :: observed(:)
        character(len=:), allocatable, intent(in) :: observed_units
        real(real64), intent(in) :: expected(:)
        logical, intent(in) :: required
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        logical :: present_in_inventory

        error_code = 0
        error_message = ''
        present_in_inventory = btest(emissions%grid_metadata%field_mask, field)
        if (.not. present_in_inventory) then
            if (required) then
                call set_error('inventory exact-grid metadata is missing required ' // trim(name), &
                               error_code, error_message)
            end if
            return
        end if
        if (.not. associated(observed) .or. size(observed) /= local_ncells) then
            call set_error('inventory ' // trim(name) // ' metadata has the wrong size', &
                           error_code, error_message)
            return
        end if
        if (.not. allocated(observed_units) .or. trim(observed_units) /= trim(expected_units)) then
            call set_error('inventory ' // trim(name) // ' units differ from the MPAS contract', &
                           error_code, error_message)
            return
        end if
        if (.not. array_values_match(observed, expected)) then
            call set_error('inventory ' // trim(name) // ' differs from the MPAS mesh', &
                           error_code, error_message)
        end if

    end subroutine compare_grid_field


    !-------------------------------------------------------------------
    !  function scalar_values_match
    !-------------------------------------------------------------------
    logical function scalar_values_match(left, right)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(real64), intent(in) :: left, right
        real(real64) :: scale

        if (.not. ieee_is_finite(left) .or. .not. ieee_is_finite(right)) then
            scalar_values_match = .false.
            return
        end if
        scale = max(1.0_real64, abs(left), abs(right))
        scalar_values_match = abs(left - right) <= &
                              64.0_real64 * epsilon(1.0_real64) * scale

    end function scalar_values_match


    !-------------------------------------------------------------------
    !  function array_values_match
    !-------------------------------------------------------------------
    logical function array_values_match(left, right)

        real(real64), intent(in) :: left(:), right(:)
        integer :: i

        array_values_match = size(left) == size(right)
        if (.not. array_values_match) return
        do i = 1, size(left)
            if (.not. scalar_values_match(left(i), right(i))) then
                array_values_match = .false.
                return
            end if
        end do

    end function array_values_match


    !-------------------------------------------------------------------
    !  routine miem_commit_mass
    !
    !> \brief Add one successful step to rank-local emitted-mass totals
    !-------------------------------------------------------------------
    subroutine miem_commit_mass(local_flux, area_cell, dt_seconds, &
                                error_code, error_message)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(real64), intent(in) :: local_flux(:,:)
        real(real64), intent(in) :: area_cell(:)
        real(real64), intent(in) :: dt_seconds
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: i, s

        error_code = 0
        error_message = ''
        if (.not. enabled) return

        if (size(local_flux, 1) /= n_species .or. &
            size(local_flux, 2) /= size(area_cell)) then
            call set_error('mass-accounting array dimensions are inconsistent', &
                           error_code, error_message)
            return
        end if
        if (.not. ieee_is_finite(dt_seconds) .or. dt_seconds <= 0.0_real64) then
            call set_error('mass-accounting timestep must be finite and positive', &
                           error_code, error_message)
            return
        end if

        do i = 1, size(area_cell)
            if (.not. ieee_is_finite(area_cell(i)) .or. area_cell(i) <= 0.0_real64) then
                call set_error('areaCell must be finite and positive', &
                               error_code, error_message)
                return
            end if
            do s = 1, n_species
                if (.not. ieee_is_finite(local_flux(s, i)) .or. &
                    (local_flux(s, i) < 0.0_real64 .and. .not. allows_net_flux(s))) then
                    call set_error('surface flux is invalid during mass accounting', &
                                   error_code, error_message)
                    return
                end if
            end do
        end do

        do s = 1, n_species
            local_cumulative_mass_kg(s) = local_cumulative_mass_kg(s) + &
                sum(local_flux(s, :) * area_cell(:)) * dt_seconds
        end do

    end subroutine miem_commit_mass


    !-------------------------------------------------------------------
    !  routine miem_finalize
    !
    !> \brief Release MIEM before its parsed mechanism dependency
    !-------------------------------------------------------------------
    subroutine miem_finalize(dminfo)

        use mpas_derived_types, only: dm_info
        use mpas_log, only: mpas_log_write
        use mpas_dmpar, only: mpas_dmpar_max_real, mpas_dmpar_sum_real_array
        use mpas_kind_types, only: RKIND

        type(dm_info), intent(in) :: dminfo

        real(kind=RKIND), allocatable :: local_mass(:)
        real(kind=RKIND), allocatable :: global_mass(:)
        real(kind=RKIND) :: local_first, global_first
        real(kind=RKIND) :: local_warm_mean, global_warm_mean
        integer :: s

        if (enabled) then
            allocate(local_mass(n_species), global_mass(n_species))
            local_mass = real(local_cumulative_mass_kg, RKIND)
            global_mass = 0.0_RKIND
            call mpas_dmpar_sum_real_array(dminfo, n_species, local_mass, global_mass)
            local_first = real(max(0.0_real64, first_run_seconds), RKIND)
            if (warm_run_calls > 0) then
                local_warm_mean = real(warm_run_seconds / real(warm_run_calls, real64), RKIND)
            else
                local_warm_mean = 0.0_RKIND
            end if
            call mpas_dmpar_max_real(dminfo, local_first, global_first)
            call mpas_dmpar_max_real(dminfo, local_warm_mean, global_warm_mean)

            if (dminfo%my_proc_id == 0) then
                do s = 1, n_species
                    if (allows_net_flux(s)) then
                        call mpas_log_write('[MIEM] Net applied surface exchange mass for ' // &
                            trim(species_names(s)) // ': $r kg', realArgs=[global_mass(s)])
                    else
                        call mpas_log_write('[MIEM] Total emitted mass for ' // &
                            trim(species_names(s)) // ': $r kg', realArgs=[global_mass(s)])
                    end if
                end do
                call mpas_log_write('[MIEM] First-open sampled step wall time: $r s', &
                                    realArgs=[global_first])
                call mpas_log_write('[MIEM] Warm sampled step mean wall time: $r s over $i calls', &
                                    realArgs=[global_warm_mean], intArgs=[warm_run_calls])
            end if
            deallocate(local_mass, global_mass)
            call mpas_log_write('[MIEM] Finalizing emissions adapter.')
        end if
        call miem_release()

    end subroutine miem_finalize


    !-------------------------------------------------------------------
    !  function miem_is_enabled
    !-------------------------------------------------------------------
    logical function miem_is_enabled()

        miem_is_enabled = enabled

    end function miem_is_enabled


    !-------------------------------------------------------------------
    !  routine miem_release
    !
    !> \brief Reset all adapter state in finalizer-safe dependency order
    !-------------------------------------------------------------------
    subroutine miem_release()

        enabled = .false.
        if (associated(emissions)) then
            deallocate(emissions)
            nullify(emissions)
        end if
        if (associated(mechanism)) then
            deallocate(mechanism)
            nullify(mechanism)
        end if
        if (allocated(species_names)) deallocate(species_names)
        if (allocated(net_flux_species)) deallocate(net_flux_species)
        if (allocated(diagnostic_sector_names)) deallocate(diagnostic_sector_names)
        if (allocated(diagnostic_sector_labels)) deallocate(diagnostic_sector_labels)
        if (allocated(diagnostic_category_ids)) deallocate(diagnostic_category_ids)
        if (allocated(species_storage_indices)) deallocate(species_storage_indices)
        if (allocated(local_cumulative_mass_kg)) deallocate(local_cumulative_mass_kg)
        if (allocated(expected_global_ids)) deallocate(expected_global_ids)
        if (allocated(expected_area_cell)) deallocate(expected_area_cell)
        if (allocated(expected_lat_cell)) deallocate(expected_lat_cell)
        if (allocated(expected_lon_cell)) deallocate(expected_lon_cell)
        if (allocated(expected_x_cell)) deallocate(expected_x_cell)
        if (allocated(expected_y_cell)) deallocate(expected_y_cell)
        if (allocated(expected_z_cell)) deallocate(expected_z_cell)
        grid_metadata_validated = .false.
        diagnostic_layered_output = .false.
        expected_on_a_sphere = .false.
        expected_is_periodic = .false.
        expected_sphere_radius = 0.0_real64
        first_run_seconds = -1.0_real64
        warm_run_seconds = 0.0_real64
        warm_run_calls = 0
        global_ncells = 0
        local_ncells = 0
        n_vert_levels = 0
        n_species = 0
        cell_stride = 0_int64
        species_stride = 0_int64
        layer_cell_stride = 0_int64
        layer_level_stride = 0_int64
        layer_species_stride = 0_int64

    end subroutine miem_release


    !-------------------------------------------------------------------
    !  routine set_musica_error
    !-------------------------------------------------------------------
    subroutine set_musica_error(context, error, error_code, error_message)

        character(len=*), intent(in) :: context
        type(error_t), intent(in) :: error
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        error_code = 1
        error_message = '[MIEM] ' // trim(context) // ': ' // trim(error%message())

    end subroutine set_musica_error


    !-------------------------------------------------------------------
    !  routine set_error
    !-------------------------------------------------------------------
    subroutine set_error(message, error_code, error_message)

        character(len=*), intent(in) :: message
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        error_code = 1
        error_message = '[MIEM] ' // trim(message)

    end subroutine set_error


    !-------------------------------------------------------------------
    !  function integer_string
    !-------------------------------------------------------------------
    function integer_string(value) result(text)

        integer, intent(in) :: value
        character(len=:), allocatable :: text
        character(len=32) :: buffer

        write(buffer, '(i0)') value
        text = trim(buffer)

    end function integer_string


    !-------------------------------------------------------------------
    !  function int64_string
    !-------------------------------------------------------------------
    function int64_string(value) result(text)

        integer(int64), intent(in) :: value
        character(len=:), allocatable :: text
        character(len=32) :: buffer

        write(buffer, '(i0)') value
        text = trim(buffer)

    end function int64_string


    !-------------------------------------------------------------------
    !  routine heap_sort
    !
    !> \brief Sort an integer array without allocating by global ID range
    !-------------------------------------------------------------------
    subroutine heap_sort(values)

        integer, intent(inout) :: values(:)

        integer :: end_pos, root, child, temp

        if (size(values) < 2) return

        do root = size(values) / 2, 1, -1
            call sift_down(values, root, size(values))
        end do
        do end_pos = size(values), 2, -1
            temp = values(1)
            values(1) = values(end_pos)
            values(end_pos) = temp
            call sift_down(values, 1, end_pos - 1)
        end do

    contains

        subroutine sift_down(items, start_pos, final_pos)

            integer, intent(inout) :: items(:)
            integer, intent(in) :: start_pos, final_pos

            integer :: child_pos, root_pos, swap_pos, swap_value

            root_pos = start_pos
            do while (2 * root_pos <= final_pos)
                child_pos = 2 * root_pos
                swap_pos = root_pos
                if (items(swap_pos) < items(child_pos)) swap_pos = child_pos
                if (child_pos + 1 <= final_pos) then
                    if (items(swap_pos) < items(child_pos + 1)) swap_pos = child_pos + 1
                end if
                if (swap_pos == root_pos) return
                swap_value = items(root_pos)
                items(root_pos) = items(swap_pos)
                items(swap_pos) = swap_value
                root_pos = swap_pos
            end do

        end subroutine sift_down

    end subroutine heap_sort

end module mpas_miem
