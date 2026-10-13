! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file.
!
!-----------------------------------------------------------------------
!  Exact-grid prescribed-field provider for CheMPAS-A.
!-----------------------------------------------------------------------
module mpas_prescribed_fields

#ifdef MPAS_USE_MUSICA
    use netcdf
#endif
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: real64
    use mpas_kind_types, only: RKIND

    implicit none
    private

    public :: prescribed_fields_init, prescribed_fields_finalize
    public :: prescribed_fields_is_enabled, prescribed_fields_n_levels
    public :: prescribed_fields_get_structure, prescribed_fields_update_o3
    public :: prescribed_fields_month_bracket

    character(len=*), parameter :: PACKAGE_SCHEMA = 'chempas-prescribed-field-package-v1'
    character(len=*), parameter :: MESH_ALGORITHM = 'chempas-mesh-sha256-v1'
    character(len=*), parameter :: TIME_UNITS = 'days since 2000-01-01 00:00:00'
    integer, parameter :: EXPECTED_MONTHS = 12
    integer, parameter :: EXPECTED_LEVELS = 11
    integer, parameter :: EXPECTED_EDGES = 12
    ! The authoritative x1.40962 coordinates originate as IEEE binary32.
    ! A pole rounded to binary32 can exceed pi/2 by about 4.4e-8 rad; this
    ! range guard allows that representation error while the owned-cell
    ! exact-grid comparison below remains at its existing tight tolerance.
    real(real64), parameter :: COORDINATE_RANGE_TOLERANCE = 5.0e-7_real64
    real(real64), parameter :: COLUMN_REL_TOL = 5.0e-13_real64

#ifdef MPAS_USE_MUSICA
    integer, save :: file_id = -1
    integer, save :: global_cells = 0
    integer, save :: local_cells = 0
    integer, save :: o3_density_var_id = -1
    integer, save :: o3_column_var_id = -1
    integer, allocatable, save :: local_storage_index(:)
    integer, allocatable, save :: storage_to_local(:)
    real(real64), allocatable, save :: upper_edges_m(:)
    real(real64), allocatable, save :: upper_temperature_k(:)
    real(real64), allocatable, save :: upper_air_density_cm3(:)
    real(real64), allocatable, save :: cached_density_left(:,:)
    real(real64), allocatable, save :: cached_density_right(:,:)
    integer, save :: cached_left_month = 0
    integer, save :: cached_right_month = 0
    logical, save :: provider_enabled = .false.
#endif

contains

    !--------------------------------------------------------------------
    !  Initialize and validate one exact-grid monthly package.
    !--------------------------------------------------------------------
    subroutine prescribed_fields_init(filename, n_global_cells, owned_cell_ids, &
                                      owned_area, owned_latitude, owned_longitude, &
                                      error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_global_cells
        integer, intent(in) :: owned_cell_ids(:)
        real(kind=RKIND), intent(in) :: owned_area(:), owned_latitude(:), owned_longitude(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        integer :: status, dim_time, dim_cells, dim_levels, dim_edges, dim_bounds
        integer :: n_time, n_cells_file, n_levels, n_edges, n_bounds
        integer :: id_var, area_var, lat_var, lon_var
        integer :: time_var, bounds_var, month_var, edge_var, layer_height_var
        integer :: layer_bounds_var, temperature_var, air_var
        integer :: i, storage, global_id
        integer, allocatable :: file_ids(:), seen(:), months(:)
        real(real64), allocatable :: file_area(:), file_lat(:), file_lon(:)
        real(real64), allocatable :: times(:), bounds(:,:), layer_column(:,:,:)
        real(real64), allocatable :: layer_density(:,:,:), dz_cm(:)
        real(real64), allocatable :: layer_height(:), layer_bounds(:,:)
        real(real64) :: scale, maximum_relative
        character(len=1024) :: attribute

        call prescribed_fields_finalize()
        error_code = 0
        error_message = ''

        if (len_trim(filename) == 0) then
            error_code = 1
            error_message = '[Prescribed fields] spatial_climatology requires a nonempty file path.'
            return
        end if
        if (n_global_cells <= 0 .or. size(owned_cell_ids) <= 0) then
            error_code = 1
            error_message = '[Prescribed fields] global and owned cell counts must be positive.'
            return
        end if
        if (size(owned_area) /= size(owned_cell_ids) .or. &
            size(owned_latitude) /= size(owned_cell_ids) .or. &
            size(owned_longitude) /= size(owned_cell_ids)) then
            error_code = 1
            error_message = '[Prescribed fields] owned mesh arrays have inconsistent sizes.'
            return
        end if

        status = nf90_open(trim(filename), NF90_NOWRITE, file_id)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('open ' // trim(filename), status, error_code, error_message)
            return
        end if

        call require_dimension('Time', dim_time, n_time, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_dimension('nCells', dim_cells, n_cells_file, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_dimension('nUpperLevels', dim_levels, n_levels, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_dimension('nUpperEdges', dim_edges, n_edges, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_dimension('nBounds', dim_bounds, n_bounds, error_code, error_message)
        if (error_code /= 0) goto 900
        if (n_time /= EXPECTED_MONTHS .or. n_cells_file /= n_global_cells .or. &
            n_levels /= EXPECTED_LEVELS .or. n_edges /= EXPECTED_EDGES .or. n_bounds /= 2) then
            error_code = 1
            error_message = '[Prescribed fields] package dimensions do not match the MVP contract.'
            goto 900
        end if

        call require_global_attribute('chempas_schema_version', PACKAGE_SCHEMA, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_global_attribute('chempas_mesh_fingerprint_algorithm', MESH_ALGORITHM, &
                                      error_code, error_message)
        if (error_code /= 0) goto 900
        call require_global_attribute('Conventions', 'CF-1.10', error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('source_product', 1, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('source_version', 1, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('source_sha256', 64, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('horizontal_remapping', 1, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('vertical_remapping', 1, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('reference_atmosphere_completion', 1, &
                                               error_code, error_message)
        if (error_code /= 0) goto 900
        call require_nonempty_global_attribute('temporal_semantics', 1, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_global_attribute('missing_value_policy', 'fail; no fill', error_code, error_message)
        if (error_code /= 0) goto 900
        call require_global_attribute('negative_value_policy', 'fail; no clipping', error_code, error_message)
        if (error_code /= 0) goto 900
        attribute = ''
        status = nf90_get_att(file_id, NF90_GLOBAL, 'chempas_mesh_fingerprint', attribute)
        if (status /= NF90_NOERR .or. len_trim(attribute) /= 64) then
            error_code = 1
            error_message = '[Prescribed fields] package has no valid mesh fingerprint.'
            goto 900
        end if

        call require_variable('indexToCellID', '1', 1, id_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('areaCell', 'm2', 1, area_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('latCell', 'radian', 1, lat_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('lonCell', 'radian', 1, lon_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('Time', TIME_UNITS, 1, time_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable_attribute(time_var, 'calendar', 'gregorian', error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable_attribute(time_var, 'climatology', 'climatology_bounds', &
                                        error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('climatology_bounds', TIME_UNITS, 2, bounds_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable_attribute(bounds_var, 'calendar', 'gregorian', error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('climatology_month', '1', 1, month_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('upper_edge_height', 'm', 1, edge_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('upper_layer_height', 'm', 1, layer_height_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable_attribute(layer_height_var, 'bounds', 'upper_layer_bounds', &
                                        error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('upper_layer_bounds', 'm', 2, layer_bounds_var, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('temperature_at_upper_edge', 'K', 1, temperature_var, &
                              error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('air_number_density_at_upper_edge', 'molecule cm-3', 1, air_var, &
                              error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('O3_layer_column', 'molecule cm-2', 3, o3_column_var_id, &
                              error_code, error_message)
        if (error_code /= 0) goto 900
        call require_variable('O3_number_density', 'molecule cm-3', 3, o3_density_var_id, &
                              error_code, error_message)
        if (error_code /= 0) goto 900
        call require_o3_dimension_order(o3_column_var_id, error_code, error_message)
        if (error_code /= 0) goto 900
        call require_o3_dimension_order(o3_density_var_id, error_code, error_message)
        if (error_code /= 0) goto 900

        allocate(file_ids(n_global_cells), seen(n_global_cells))
        allocate(file_area(n_global_cells), file_lat(n_global_cells), file_lon(n_global_cells))
        status = nf90_get_var(file_id, id_var, file_ids)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read indexToCellID', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, area_var, file_area)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read areaCell', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, lat_var, file_lat)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read latCell', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, lon_var, file_lon)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read lonCell', status, error_code, error_message)
            goto 900
        end if
        if (.not. all(ieee_is_finite(file_area)) .or. .not. all(ieee_is_finite(file_lat)) .or. &
            .not. all(ieee_is_finite(file_lon)) .or. any(file_area <= 0.0_real64) .or. &
            any(abs(file_lat) > 0.5_real64 * acos(-1.0_real64) + COORDINATE_RANGE_TOLERANCE) .or. &
            any(abs(file_lon) > 2.0_real64 * acos(-1.0_real64) + COORDINATE_RANGE_TOLERANCE)) then
            error_code = 1
            error_message = '[Prescribed fields] package mesh coordinates are invalid.'
            goto 900
        end if
        seen = 0
        do storage = 1, n_global_cells
            global_id = file_ids(storage)
            if (global_id < 1 .or. global_id > n_global_cells) then
                error_code = 1
                error_message = '[Prescribed fields] indexToCellID is outside 1..nCells.'
                goto 900
            end if
            if (seen(global_id) /= 0) then
                error_code = 1
                error_message = '[Prescribed fields] indexToCellID contains a duplicate.'
                goto 900
            end if
            seen(global_id) = storage
        end do

        global_cells = n_global_cells
        local_cells = size(owned_cell_ids)
        allocate(local_storage_index(local_cells), storage_to_local(global_cells))
        storage_to_local = 0
        do i = 1, local_cells
            global_id = owned_cell_ids(i)
            if (global_id < 1 .or. global_id > global_cells) then
                error_code = 1
                error_message = '[Prescribed fields] owned cell ID is outside 1..nCells.'
                goto 900
            end if
            storage = seen(global_id)
            if (storage <= 0 .or. storage_to_local(storage) /= 0) then
                error_code = 1
                error_message = '[Prescribed fields] owned cell IDs are missing or duplicated.'
                goto 900
            end if
            local_storage_index(i) = storage
            storage_to_local(storage) = i
            scale = max(1.0_real64, abs(file_area(storage)), abs(real(owned_area(i), real64)))
            if (abs(file_area(storage) - real(owned_area(i), real64)) > 2.0e-13_real64 * scale .or. &
                abs(file_lat(storage) - real(owned_latitude(i), real64)) > 2.0e-13_real64 .or. &
                abs(file_lon(storage) - real(owned_longitude(i), real64)) > 2.0e-13_real64) then
                error_code = 1
                error_message = '[Prescribed fields] exact-grid mesh metadata differs on an owned cell.'
                goto 900
            end if
        end do

        allocate(times(EXPECTED_MONTHS), bounds(2, EXPECTED_MONTHS), months(EXPECTED_MONTHS))
        status = nf90_get_var(file_id, time_var, times)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read Time', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, bounds_var, bounds)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read climatology_bounds', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, month_var, months)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read climatology_month', status, error_code, error_message)
            goto 900
        end if
        call validate_climatology_coordinates(times, bounds, months, error_code, error_message)
        if (error_code /= 0) goto 900

        allocate(upper_edges_m(EXPECTED_EDGES), upper_temperature_k(EXPECTED_EDGES), &
                 upper_air_density_cm3(EXPECTED_EDGES), layer_height(EXPECTED_LEVELS), &
                 layer_bounds(2, EXPECTED_LEVELS))
        status = nf90_get_var(file_id, edge_var, upper_edges_m)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read upper_edge_height', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, temperature_var, upper_temperature_k)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read temperature_at_upper_edge', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, air_var, upper_air_density_cm3)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read air_number_density_at_upper_edge', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, layer_height_var, layer_height)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read upper_layer_height', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, layer_bounds_var, layer_bounds)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read upper_layer_bounds', status, error_code, error_message)
            goto 900
        end if
        do i = 1, EXPECTED_EDGES
            if (.not. ieee_is_finite(upper_edges_m(i)) .or. &
                .not. ieee_is_finite(upper_temperature_k(i)) .or. &
                .not. ieee_is_finite(upper_air_density_cm3(i)) .or. &
                upper_temperature_k(i) <= 0.0_real64 .or. upper_air_density_cm3(i) <= 0.0_real64) then
                error_code = 1
                error_message = '[Prescribed fields] structural upper-atmosphere fields are invalid.'
                goto 900
            end if
            if (i > 1 .and. upper_edges_m(i) <= upper_edges_m(i - 1)) then
                error_code = 1
                error_message = '[Prescribed fields] upper height edges must be strictly increasing.'
                goto 900
            end if
        end do
        do i = 1, EXPECTED_LEVELS
            if (.not. ieee_is_finite(layer_height(i)) .or. &
                abs(layer_height(i) - 0.5_real64 * (upper_edges_m(i) + upper_edges_m(i + 1))) > &
                    1.0e-10_real64 .or. &
                abs(layer_bounds(1,i) - upper_edges_m(i)) > 1.0e-10_real64 .or. &
                abs(layer_bounds(2,i) - upper_edges_m(i + 1)) > 1.0e-10_real64) then
                error_code = 1
                error_message = '[Prescribed fields] upper-layer coordinates or bounds changed.'
                goto 900
            end if
        end do

        ! Validate one selected slab before promotion. Every newly loaded slab
        ! is checked again, so later corrupt records fail before use.
        allocate(layer_column(EXPECTED_LEVELS, 1, 1), layer_density(EXPECTED_LEVELS, 1, 1), &
                 dz_cm(EXPECTED_LEVELS))
        status = nf90_get_var(file_id, o3_column_var_id, layer_column, &
                              start=(/1, local_storage_index(1), 1/), &
                              count=(/EXPECTED_LEVELS, 1, 1/))
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read O3_layer_column probe', status, error_code, error_message)
            goto 900
        end if
        status = nf90_get_var(file_id, o3_density_var_id, layer_density, &
                              start=(/1, local_storage_index(1), 1/), &
                              count=(/EXPECTED_LEVELS, 1, 1/))
        if (status /= NF90_NOERR) then
            call set_netcdf_error('read O3_number_density probe', status, error_code, error_message)
            goto 900
        end if
        dz_cm = (upper_edges_m(2:) - upper_edges_m(:EXPECTED_LEVELS)) * 100.0_real64
        maximum_relative = maxval(abs(layer_density(:,1,1) * dz_cm - layer_column(:,1,1)) / &
                                  max(layer_column(:,1,1), tiny(1.0_real64)))
        if (.not. all(ieee_is_finite(layer_column)) .or. &
            .not. all(ieee_is_finite(layer_density)) .or. &
            any(layer_column < 0.0_real64) .or. any(layer_density < 0.0_real64) .or. &
            maximum_relative > COLUMN_REL_TOL) then
            error_code = 1
            error_message = '[Prescribed fields] O3 density/column probe failed.'
            goto 900
        end if

        allocate(cached_density_left(EXPECTED_LEVELS, local_cells), &
                 cached_density_right(EXPECTED_LEVELS, local_cells))
        cached_density_left = 0.0_real64
        cached_density_right = 0.0_real64
        cached_left_month = 0
        cached_right_month = 0
        provider_enabled = .true.
        call mpas_log_write('[Prescribed fields] Initialized spatial monthly O3 from ' // trim(filename))
        return

900     continue
        call prescribed_fields_finalize()
#else
        error_code = 1
        error_message = '[Prescribed fields] MUSICA/NetCDF support is not available in this build.'
#endif

    end subroutine prescribed_fields_init


    !--------------------------------------------------------------------
    !  Return the static structural extension arrays.
    !--------------------------------------------------------------------
    subroutine prescribed_fields_get_structure(z_edges_m, temperature_edges_k, &
                                               air_density_edges_cm3, error_code, error_message)
        real(kind=RKIND), allocatable, intent(out) :: z_edges_m(:)
        real(kind=RKIND), allocatable, intent(out) :: temperature_edges_k(:)
        real(kind=RKIND), allocatable, intent(out) :: air_density_edges_cm3(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        error_code = 0
        error_message = ''
#ifdef MPAS_USE_MUSICA
        if (.not. provider_enabled) then
            error_code = 1
            error_message = '[Prescribed fields] provider is not initialized.'
            return
        end if
        allocate(z_edges_m(EXPECTED_EDGES), temperature_edges_k(EXPECTED_EDGES), &
                 air_density_edges_cm3(EXPECTED_EDGES))
        z_edges_m = real(upper_edges_m, RKIND)
        temperature_edges_k = real(upper_temperature_k, RKIND)
        air_density_edges_cm3 = real(upper_air_density_cm3, RKIND)
#else
        error_code = 1
        error_message = '[Prescribed fields] provider is unavailable.'
#endif
    end subroutine prescribed_fields_get_structure


    !--------------------------------------------------------------------
    !  Interpolate the two cached monthly slabs at a Gregorian time.
    !--------------------------------------------------------------------
    subroutine prescribed_fields_update_o3(year, month, day, hour, minute, second, &
                                           o3_density, error_code, error_message)
        integer, intent(in) :: year, month, day, hour, minute, second
        real(kind=RKIND), intent(out) :: o3_density(:,:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        integer :: left_month, right_month
        real(real64) :: weight

        error_code = 0
        error_message = ''
        if (.not. provider_enabled) then
            error_code = 1
            error_message = '[Prescribed fields] update requested before initialization.'
            return
        end if
        if (size(o3_density, 1) /= EXPECTED_LEVELS .or. size(o3_density, 2) /= local_cells) then
            error_code = 1
            error_message = '[Prescribed fields] O3 output array has the wrong shape.'
            return
        end if
        call prescribed_fields_month_bracket(year, month, day, hour, minute, second, &
                                             left_month, right_month, weight, &
                                             error_code, error_message)
        if (error_code /= 0) return
        if (left_month /= cached_left_month .or. right_month /= cached_right_month) then
            call load_month(left_month, cached_density_left, error_code, error_message)
            if (error_code /= 0) return
            call load_month(right_month, cached_density_right, error_code, error_message)
            if (error_code /= 0) return
            cached_left_month = left_month
            cached_right_month = right_month
        end if
        o3_density = real((1.0_real64 - weight) * cached_density_left + &
                          weight * cached_density_right, RKIND)
        if (.not. all(ieee_is_finite(o3_density)) .or. any(o3_density < 0.0_RKIND)) then
            error_code = 1
            error_message = '[Prescribed fields] interpolated O3 is nonfinite or negative.'
            return
        end if
#else
        error_code = 1
        error_message = '[Prescribed fields] provider is unavailable.'
        o3_density = 0.0_RKIND
#endif
    end subroutine prescribed_fields_update_o3


    !--------------------------------------------------------------------
    !  Compute cyclic monthly brackets using actual Gregorian month lengths.
    !--------------------------------------------------------------------
    subroutine prescribed_fields_month_bracket(year, month, day, hour, minute, second, &
                                               left_month, right_month, weight, &
                                               error_code, error_message)
        integer, intent(in) :: year, month, day, hour, minute, second
        integer, intent(out) :: left_month, right_month
        real(real64), intent(out) :: weight
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: days(12), starts(12), i
        real(real64) :: anchors(12), current_seconds, left_anchor, right_anchor

        error_code = 0
        error_message = ''
        left_month = 0
        right_month = 0
        weight = 0.0_real64
        if (year < 1 .or. month < 1 .or. month > 12 .or. &
            hour < 0 .or. hour > 23 .or. minute < 0 .or. minute > 59 .or. &
            second < 0 .or. second > 59) then
            error_code = 1
            error_message = '[Prescribed fields] invalid Gregorian date/time.'
            return
        end if
        do i = 1, 12
            days(i) = days_in_month(year, i)
        end do
        if (day < 1 .or. day > days(month)) then
            error_code = 1
            error_message = '[Prescribed fields] invalid Gregorian day-of-month.'
            return
        end if
        starts(1) = 0
        do i = 2, 12
            starts(i) = starts(i - 1) + days(i - 1)
        end do
        do i = 1, 12
            anchors(i) = (real(starts(i), real64) + 0.5_real64 * real(days(i), real64)) * 86400.0_real64
        end do
        current_seconds = (real(starts(month) + day - 1, real64) * 86400.0_real64) + &
                          real(hour * 3600 + minute * 60 + second, real64)

        if (current_seconds < anchors(1)) then
            left_month = 12
            right_month = 1
            left_anchor = -0.5_real64 * real(days_in_month(year - 1, 12), real64) * 86400.0_real64
            right_anchor = anchors(1)
        else if (current_seconds >= anchors(12)) then
            left_month = 12
            right_month = 1
            left_anchor = anchors(12)
            right_anchor = real(sum(days), real64) * 86400.0_real64 + &
                           0.5_real64 * real(days_in_month(year + 1, 1), real64) * 86400.0_real64
        else
            do i = 1, 11
                if (current_seconds >= anchors(i) .and. current_seconds < anchors(i + 1)) then
                    left_month = i
                    right_month = i + 1
                    left_anchor = anchors(i)
                    right_anchor = anchors(i + 1)
                    exit
                end if
            end do
        end if
        if (left_month == 0 .or. right_anchor <= left_anchor) then
            error_code = 1
            error_message = '[Prescribed fields] could not resolve monthly climatology brackets.'
            return
        end if
        weight = (current_seconds - left_anchor) / (right_anchor - left_anchor)
        if (.not. ieee_is_finite(weight) .or. weight < 0.0_real64 .or. weight > 1.0_real64) then
            error_code = 1
            error_message = '[Prescribed fields] monthly interpolation weight is invalid.'
            return
        end if
    end subroutine prescribed_fields_month_bracket


    logical function prescribed_fields_is_enabled()
#ifdef MPAS_USE_MUSICA
        prescribed_fields_is_enabled = provider_enabled
#else
        prescribed_fields_is_enabled = .false.
#endif
    end function prescribed_fields_is_enabled


    integer function prescribed_fields_n_levels()
#ifdef MPAS_USE_MUSICA
        if (provider_enabled) then
            prescribed_fields_n_levels = EXPECTED_LEVELS
        else
            prescribed_fields_n_levels = 0
        end if
#else
        prescribed_fields_n_levels = 0
#endif
    end function prescribed_fields_n_levels


    subroutine prescribed_fields_finalize()
#ifdef MPAS_USE_MUSICA
        integer :: status
        if (file_id >= 0) status = nf90_close(file_id)
        file_id = -1
        global_cells = 0
        local_cells = 0
        o3_density_var_id = -1
        o3_column_var_id = -1
        cached_left_month = 0
        cached_right_month = 0
        provider_enabled = .false.
        if (allocated(local_storage_index)) deallocate(local_storage_index)
        if (allocated(storage_to_local)) deallocate(storage_to_local)
        if (allocated(upper_edges_m)) deallocate(upper_edges_m)
        if (allocated(upper_temperature_k)) deallocate(upper_temperature_k)
        if (allocated(upper_air_density_cm3)) deallocate(upper_air_density_cm3)
        if (allocated(cached_density_left)) deallocate(cached_density_left)
        if (allocated(cached_density_right)) deallocate(cached_density_right)
#endif
    end subroutine prescribed_fields_finalize


#ifdef MPAS_USE_MUSICA
    subroutine load_month(month, destination, error_code, error_message)
        integer, intent(in) :: month
        real(real64), intent(out) :: destination(:,:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: storage, run_start, run_stop, run_length, local_index, j, status
        real(real64), allocatable :: density(:,:,:), column(:,:,:), dz_cm(:)
        real(real64) :: relative

        error_code = 0
        error_message = ''
        if (month < 1 .or. month > 12) then
            error_code = 1
            error_message = '[Prescribed fields] requested month is outside 1..12.'
            return
        end if
        destination = 0.0_real64
        allocate(dz_cm(EXPECTED_LEVELS))
        dz_cm = (upper_edges_m(2:) - upper_edges_m(:EXPECTED_LEVELS)) * 100.0_real64
        storage = 1
        do while (storage <= global_cells)
            if (storage_to_local(storage) == 0) then
                storage = storage + 1
                cycle
            end if
            run_start = storage
            run_stop = storage
            do while (run_stop < global_cells)
                if (storage_to_local(run_stop + 1) == 0) exit
                run_stop = run_stop + 1
            end do
            run_length = run_stop - run_start + 1
            allocate(density(EXPECTED_LEVELS, run_length, 1), &
                     column(EXPECTED_LEVELS, run_length, 1))
            status = nf90_get_var(file_id, o3_density_var_id, density, &
                                  start=(/1, run_start, month/), &
                                  count=(/EXPECTED_LEVELS, run_length, 1/))
            if (status /= NF90_NOERR) then
                call set_netcdf_error('read selected O3_number_density slab', status, &
                                      error_code, error_message)
                deallocate(density, column, dz_cm)
                return
            end if
            status = nf90_get_var(file_id, o3_column_var_id, column, &
                                  start=(/1, run_start, month/), &
                                  count=(/EXPECTED_LEVELS, run_length, 1/))
            if (status /= NF90_NOERR) then
                call set_netcdf_error('read selected O3_layer_column slab', status, &
                                      error_code, error_message)
                deallocate(density, column, dz_cm)
                return
            end if
            if (.not. all(ieee_is_finite(density)) .or. .not. all(ieee_is_finite(column)) .or. &
                any(density < 0.0_real64) .or. any(column < 0.0_real64)) then
                error_code = 1
                error_message = '[Prescribed fields] selected O3 slab is nonfinite or negative.'
                deallocate(density, column, dz_cm)
                return
            end if
            do j = 1, run_length
                relative = maxval(abs(density(:,j,1) * dz_cm - column(:,j,1)) / &
                                  max(column(:,j,1), tiny(1.0_real64)))
                if (relative > COLUMN_REL_TOL) then
                    error_code = 1
                    error_message = '[Prescribed fields] selected O3 density/column closure failed.'
                    deallocate(density, column, dz_cm)
                    return
                end if
                local_index = storage_to_local(run_start + j - 1)
                destination(:, local_index) = density(:, j, 1)
            end do
            deallocate(density, column)
            storage = run_stop + 1
        end do
        deallocate(dz_cm)
    end subroutine load_month


    subroutine require_dimension(name, dimension_id, length, error_code, error_message)
        character(len=*), intent(in) :: name
        integer, intent(out) :: dimension_id, length, error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: status

        error_code = 0
        error_message = ''
        status = nf90_inq_dimid(file_id, trim(name), dimension_id)
        if (status == NF90_NOERR) status = nf90_inquire_dimension(file_id, dimension_id, len=length)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('require dimension ' // trim(name), status, error_code, error_message)
        end if
    end subroutine require_dimension


    subroutine require_variable(name, units, expected_ndims, variable_id, error_code, error_message)
        character(len=*), intent(in) :: name, units
        integer, intent(in) :: expected_ndims
        integer, intent(out) :: variable_id, error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: status, ndims
        character(len=256) :: observed_units

        error_code = 0
        error_message = ''
        status = nf90_inq_varid(file_id, trim(name), variable_id)
        if (status == NF90_NOERR) status = nf90_inquire_variable(file_id, variable_id, ndims=ndims)
        if (status /= NF90_NOERR) then
            call set_netcdf_error('require variable ' // trim(name), status, error_code, error_message)
            return
        end if
        if (ndims /= expected_ndims) then
            error_code = 1
            error_message = '[Prescribed fields] variable ' // trim(name) // ' has the wrong rank.'
            return
        end if
        observed_units = ''
        status = nf90_get_att(file_id, variable_id, 'units', observed_units)
        if (status /= NF90_NOERR .or. trim(observed_units) /= trim(units)) then
            error_code = 1
            error_message = '[Prescribed fields] variable ' // trim(name) // ' has wrong or missing units.'
            return
        end if
    end subroutine require_variable


    subroutine require_o3_dimension_order(variable_id, error_code, error_message)
        integer, intent(in) :: variable_id
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: status, ndims, dimids(3), i
        character(len=NF90_MAX_NAME) :: names(3)

        error_code = 0
        error_message = ''
        status = nf90_inquire_variable(file_id, variable_id, ndims=ndims, dimids=dimids)
        if (status /= NF90_NOERR .or. ndims /= 3) then
            error_code = 1
            error_message = '[Prescribed fields] O3 variable rank inquiry failed.'
            return
        end if
        do i = 1, 3
            status = nf90_inquire_dimension(file_id, dimids(i), name=names(i))
            if (status /= NF90_NOERR) then
                call set_netcdf_error('inquire O3 dimension', status, error_code, error_message)
                return
            end if
        end do
        if (trim(names(1)) /= 'nUpperLevels' .or. trim(names(2)) /= 'nCells' .or. &
            trim(names(3)) /= 'Time') then
            error_code = 1
            error_message = '[Prescribed fields] O3 variables must expose Fortran dimensions ' // &
                            '(nUpperLevels,nCells,Time).'
        end if
    end subroutine require_o3_dimension_order


    subroutine require_global_attribute(name, expected, error_code, error_message)
        character(len=*), intent(in) :: name, expected
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: status
        character(len=1024) :: observed

        error_code = 0
        error_message = ''
        observed = ''
        status = nf90_get_att(file_id, NF90_GLOBAL, trim(name), observed)
        if (status /= NF90_NOERR .or. trim(observed) /= trim(expected)) then
            error_code = 1
            error_message = '[Prescribed fields] global attribute ' // trim(name) // ' changed.'
        end if
    end subroutine require_global_attribute


    subroutine require_nonempty_global_attribute(name, minimum_length, error_code, error_message)
        character(len=*), intent(in) :: name
        integer, intent(in) :: minimum_length
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: status
        character(len=2048) :: observed

        error_code = 0
        error_message = ''
        observed = ''
        status = nf90_get_att(file_id, NF90_GLOBAL, trim(name), observed)
        if (status /= NF90_NOERR .or. len_trim(observed) < minimum_length) then
            error_code = 1
            error_message = '[Prescribed fields] global attribute ' // trim(name) // &
                            ' is missing or incomplete.'
        end if
    end subroutine require_nonempty_global_attribute


    subroutine require_variable_attribute(variable_id, name, expected, error_code, error_message)
        integer, intent(in) :: variable_id
        character(len=*), intent(in) :: name, expected
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: status
        character(len=1024) :: observed

        error_code = 0
        error_message = ''
        observed = ''
        status = nf90_get_att(file_id, variable_id, trim(name), observed)
        if (status /= NF90_NOERR .or. trim(observed) /= trim(expected)) then
            error_code = 1
            error_message = '[Prescribed fields] variable attribute ' // trim(name) // ' changed.'
        end if
    end subroutine require_variable_attribute


    subroutine validate_climatology_coordinates(times, bounds, months, error_code, error_message)
        real(real64), intent(in) :: times(:), bounds(:,:)
        integer, intent(in) :: months(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        integer :: expected_days(12), i, start_day
        real(real64) :: expected_left, expected_right

        error_code = 0
        error_message = ''
        expected_days = (/31,29,31,30,31,30,31,31,30,31,30,31/)
        start_day = 0
        do i = 1, 12
            expected_left = real(start_day, real64)
            expected_right = real(start_day + expected_days(i), real64)
            if (months(i) /= i .or. abs(bounds(1,i) - expected_left) > 1.0e-12_real64 .or. &
                abs(bounds(2,i) - expected_right) > 1.0e-12_real64 .or. &
                abs(times(i) - 0.5_real64 * (expected_left + expected_right)) > 1.0e-12_real64) then
                error_code = 1
                error_message = '[Prescribed fields] climatology coordinates or anchors changed.'
                return
            end if
            start_day = start_day + expected_days(i)
        end do
    end subroutine validate_climatology_coordinates


    subroutine set_netcdf_error(context, status, error_code, error_message)
        character(len=*), intent(in) :: context
        integer, intent(in) :: status
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message
        error_code = 1
        error_message = '[Prescribed fields] NetCDF ' // trim(context) // ': ' // &
                        trim(nf90_strerror(status))
    end subroutine set_netcdf_error
#endif


    integer function days_in_month(year, month)
        integer, intent(in) :: year, month
        integer :: days(12)
        days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
        if (month < 1 .or. month > 12) then
            days_in_month = 0
        else
            days_in_month = days(month)
            if (month == 2 .and. is_gregorian_leap_year(year)) days_in_month = 29
        end if
    end function days_in_month


    logical function is_gregorian_leap_year(year)
        integer, intent(in) :: year
        is_gregorian_leap_year = mod(year, 4) == 0 .and. &
            (mod(year, 100) /= 0 .or. mod(year, 400) == 0)
    end function is_gregorian_leap_year

end module mpas_prescribed_fields
