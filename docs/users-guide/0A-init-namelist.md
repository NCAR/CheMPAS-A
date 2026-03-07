# Appendix A: Initialization Namelist Options

This appendix summarizes the complete set of namelist options available when running the MPAS non-hydrostatic atmosphere initialization core. The applicability of certain options depends on the type of initial conditions to be created -- idealized or 'real-data' -- and such applicability is identified in the description when it exists.

Date-time strings throughout all MPAS namelists assume a common format. Specifically, time intervals are of the form `'[DDD_]HH:MM:SS[.sss]'`, where `DDD` is an integer number of days with any number of digits, `HH` is a two-digit hour value, `MM` is a two-digit minute value, `SS` is a two-digit second value, and `sss` are fractions of a second with any number of digits; any part of the time interval format in square brackets (`[ ]`) may be omitted, and if days are omitted, `HH` may be either a one- or two-digit hour specification. Time instants (e.g., start time or end time) are of the form `'YYYY-MM-DD[_HH:MM:SS[.sss]]'`, where `YYYY` is an integer year with any number of digits, `MM` is a two-digit month value, `DD` is a two-digit day value, and `HH:MM:SS.sss` is a time with the same format as in a time interval specification. For both time instants and time intervals, a value of `'none'` represents 'no value'.

## A.1 nhyd_model

### `config_init_case` (integer)

| | |
|---|---|
| Units | - |
| Description | Type of initial conditions to create: 1 = Jablonowski & Williamson baroclinic wave (no initial perturbation), 2 = Jablonowski & Williamson baroclinic wave (with initial perturbation), 3 = Jablonowski & Williamson baroclinic wave (with normal-mode perturbation), 4 = squall line, 5 = super-cell, 6 = mountain wave, 7 = real-data initial conditions from, e.g., GFS, 8 = surface field (SST, sea-ice) update file for use with real-data simulations, 9 = lateral boundary conditions update file for use with real-data simulations, 13 = CAM-MPAS 3-d grid with specified topography and zeta levels |
| Possible Values | 1 -- 9, or 13 *(default: 7)* |

### `config_calendar_type` (character)

| | |
|---|---|
| Units | - |
| Description | Simulation calendar type *(hidden by default)* |
| Possible Values | `'gregorian'`, `'gregorian_noleap'` *(default: gregorian)* |

### `config_start_time` (character)

| | |
|---|---|
| Units | - |
| Description | Time to begin processing first-guess data (cases 7, 8, and 9 only) |
| Possible Values | `'YYYY-MM-DD_hh:mm:ss'` *(default: 2010-10-23_00:00:00)* |

### `config_stop_time` (character)

| | |
|---|---|
| Units | - |
| Description | Time to end processing first-guess data (cases 8 and 9 only) |
| Possible Values | `'YYYY-MM-DD_hh:mm:ss'` *(default: 2010-10-23_00:00:00)* |

### `config_theta_adv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Horizontal advection order for theta |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_coef_3rd_order` (real)

| | |
|---|---|
| Units | - |
| Description | Upwinding coefficient in the 3rd order advection scheme |
| Possible Values | 0 <= config_coef_3rd_order <= 1 *(default: 0.25)* |

### `config_num_halos` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of halo layers for fields *(hidden by default)* |
| Possible Values | Integer values, typically 2 or 3; DO NOT CHANGE *(default: 2)* |

### `config_interface_projection` (character)

| | |
|---|---|
| Units | - |
| Description | Projecting layer values to the interface, linear vertical interpolation or integral |
| Possible Values | `'linear_interpolation'` or `'layer_integral'` (layer average value) *(default: linear_interpolation)* |

## A.2 dimensions

### `config_nvertlevels` (integer)

| | |
|---|---|
| Units | - |
| Description | The number of vertical levels to be used in the model |
| Possible Values | Positive integer values *(default: 55)* |

### `config_nsoillevels` (integer)

| | |
|---|---|
| Units | - |
| Description | The number of vertical soil levels needed by LSM in the model (case 7 only) |
| Possible Values | Positive integer values *(default: 4)* |

### `config_nfglevels` (integer)

| | |
|---|---|
| Units | - |
| Description | The number of atmospheric levels (including surface and sea-level) in the first-guess dataset (cases 7 and 9 only) |
| Possible Values | Positive integer values *(default: 38)* |

### `config_nfgsoillevels` (integer)

| | |
|---|---|
| Units | - |
| Description | The number of vertical soil levels in the first-guess dataset (case 7 only) |
| Possible Values | Positive integer values *(default: 4)* |

### `config_gocartlevels` (integer)

| | |
|---|---|
| Units | - |
| Description | The number of vertical GOCART levels in the climatological GOCART data set |
| Possible Values | Positive integer values *(default: 30)* |

### `config_months` (integer)

| | |
|---|---|
| Units | - |
| Description | The number of months in a year *(hidden by default)* |
| Possible Values | Positive integer values *(default: 12)* |

## A.3 data_sources

### `config_geog_data_path` (character)

| | |
|---|---|
| Units | - |
| Description | Path to the WPS static data files (case 7 only) |
| Possible Values | Any valid path *(default: /glade/work/wrfhelp/WPS_GEOG/)* |

### `config_met_prefix` (character)

| | |
|---|---|
| Units | - |
| Description | Filename prefix of ungrib intermediate file to use for initial conditions (cases 7 and 9 only) |
| Possible Values | Any alpha-numeric string *(default: CFSR)* |

### `config_sfc_prefix` (character)

| | |
|---|---|
| Units | - |
| Description | Filename prefix of ungrib intermediate file to use for SST and sea-ice (cases 7 and 8 only) |
| Possible Values | Any alpha-numeric string *(default: SST)* |

### `config_fg_interval` (integer)

| | |
|---|---|
| Units | - |
| Description | Interval between intermediate files (cases 8 and 9 only) |
| Possible Values | `[DDD_]hh:mm:ss` *(default: 86400)* |

### `config_landuse_data` (character)

| | |
|---|---|
| Units | - |
| Description | The land use classification to use (case 7 only) |
| Possible Values | `'USGS'`, `'MODIFIED_IGBP_MODIS_NOAH'`, or `'MODIFIED_IGBP_MODIS_NOAH_15s'` *(default: MODIFIED_IGBP_MODIS_NOAH)* |

### `config_soilcat_data` (character)

| | |
|---|---|
| Units | - |
| Description | The soil category classification to use |
| Possible Values | `'STATSGO'` or `'BNU'` *(default: STATSGO)* |

### `config_topo_data` (character)

| | |
|---|---|
| Units | - |
| Description | The topography dataset to use for the model terrain and for GWDO static fields (case 7 only) |
| Possible Values | `'GTOPO30'` or `'GMTED2010'` *(default: GMTED2010)* |

### `config_vegfrac_data` (character)

| | |
|---|---|
| Units | - |
| Description | The climatological monthly vegetation fraction dataset to use (case 7 only) |
| Possible Values | `'MODIS'` or `'NCEP'` *(default: MODIS)* |

### `config_albedo_data` (character)

| | |
|---|---|
| Units | - |
| Description | The climatological monthly albedo dataset to use (case 7 only) |
| Possible Values | `'MODIS'` or `'NCEP'` *(default: MODIS)* |

### `config_maxsnowalbedo_data` (character)

| | |
|---|---|
| Units | - |
| Description | The maximum snow albedo dataset to use (case 7 only) |
| Possible Values | `'MODIS'` or `'NCEP'` *(default: MODIS)* |

### `config_supersample_factor` (integer)

| | |
|---|---|
| Units | - |
| Description | The supersampling factor to be used for MODIS maximum snow albedo and monthly albedo datasets (case 7 only) |
| Possible Values | Positive integer values *(default: 3)* |

### `config_lu_supersample_factor` (integer)

| | |
|---|---|
| Units | - |
| Description | The supersampling factor to be used for 30s or 15s MODIS land use, or for 30s USGS land use, as selected by `config_landuse_data` (case 7 only) |
| Possible Values | Positive integer values *(default: 1)* |

### `config_30s_supersample_factor` (integer)

| | |
|---|---|
| Units | - |
| Description | The supersampling factor to be used for 30s terrain, soil category, and MODIS FPAR monthly vegetation fraction (case 7 only) |
| Possible Values | Positive integer values *(default: 1)* |

### `config_noahmp_static` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to process, read, and write static fields only required by Noah-MP *(hidden by default)* |
| Possible Values | true or false *(default: true)* |

### `config_use_spechumd` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to use specific-humidity as the first-guess moisture variable. If this option is false, relative humidity will be used. |
| Possible Values | true or false *(default: false)* |

## A.4 vertical_grid

### `config_ztop` (real)

| | |
|---|---|
| Units | m |
| Description | Model top height |
| Possible Values | Positive real values *(default: 30000.0)* |

### `config_nsmterrain` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of smoothing passes to apply to the interpolated terrain field |
| Possible Values | Non-negative integer values *(default: 1)* |

### `config_smooth_surfaces` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to smooth zeta surfaces |
| Possible Values | true or false *(default: true)* |

### `config_dzmin` (real)

| | |
|---|---|
| Units | - |
| Description | Minimum thickness of layers as a fraction of nominal thickness |
| Possible Values | Real values in the interval (0,1) *(default: 0.3)* |

### `config_nsm` (integer)

| | |
|---|---|
| Units | - |
| Description | Maximum number of smoothing passes for coordinate surfaces |
| Possible Values | Positive integer values *(default: 30)* |

### `config_tc_vertical_grid` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to use the vertical layer profile that was developed for use in real-time TC experiments |
| Possible Values | true or false *(default: true)* |

### `config_specified_zeta_levels` (character)

| | |
|---|---|
| Units | - |
| Description | The name of a text file containing a list of vertical coordinate (zeta) values in increasing order at layer interfaces *(hidden by default)* |
| Possible Values | Any valid filename *(default: )* |

### `config_blend_bdy_terrain` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to blend terrain along domain boundaries with first-guess terrain. Only useful for limited-area domains. |
| Possible Values | true or false *(default: false)* |

## A.5 interpolation_control

### `config_extrap_airtemp` (character)

| | |
|---|---|
| Units | - |
| Description | Method of extrapolation of air temperature above/below first-guess levels |
| Possible Values | `'constant'` (last valid value), `'linear'` (linear extrapolation based on last two values), `'lapse-rate'` (0.0065 K/m from last valid value) *(default: lapse-rate)* |

## A.6 preproc_stages

### `config_static_interp` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to interpolate WPS static data (case 7 only) |
| Possible Values | true or false *(default: true)* |

### `config_native_gwd_static` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to recompute subgrid-scale orography statistics directly on the native MPAS mesh (case 7 only) |
| Possible Values | true or false *(default: true)* |

### `config_native_gwd_gsl_static` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to recompute subgrid-scale orography statistics for the GSL drag directly on the native MPAS mesh (case 7 only) |
| Possible Values | true or false *(default: false)* |

### `config_gwd_cell_scaling` (real)

| | |
|---|---|
| Units | - |
| Description | Scaling factor for the effective grid cell diameter used in computation of GWD static fields *(hidden by default)* |
| Possible Values | Positive real values *(default: 1.0)* |

### `config_vertical_grid` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to generate vertical grid |
| Possible Values | true or false *(default: true)* |

### `config_met_interp` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to interpolate first-guess fields from intermediate file |
| Possible Values | true or false *(default: true)* |

### `config_input_sst` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to re-compute SST and sea-ice fields from surface input data set; should be set to `.true.` when running case 8 |
| Possible Values | true or false *(default: false)* |

### `config_frac_seaice` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to switch sea-ice threshold from 0.5 to 0.02 |
| Possible Values | true or false *(default: true)* |

## A.7 physics

### `config_tsk_seaice_threshold` (real)

| | |
|---|---|
| Units | K |
| Description | Surface temperature threshold below which water points are set to sea-ice points *(hidden by default)* |
| Possible Values | Positive real values *(default: 100.)* |

## A.8 io

### `config_pio_num_iotasks` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of tasks to perform file I/O |
| Possible Values | Integer valued, 0 <= config_pio_num_iotasks <= # MPI tasks, 0 indicates all tasks perform I/O *(default: 0)* |

### `config_pio_stride` (integer)

| | |
|---|---|
| Units | - |
| Description | Stride between file I/O tasks |
| Possible Values | Integer valued, <= (# MPI tasks) / config_pio_num_iotasks *(default: 1)* |

## A.9 decomposition

### `config_block_decomp_file_prefix` (character)

| | |
|---|---|
| Units | - |
| Description | Prefix of graph decomposition file, to be suffixed with the MPI task count |
| Possible Values | Any valid filename *(default: x1.40962.graph.info.part.)* |

### `config_number_of_blocks` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of blocks to assign to each MPI task *(hidden by default)* |
| Possible Values | Positive integer values *(default: 0)* |

### `config_explicit_proc_decomp` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to use an explicit mapping of blocks to MPI tasks *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_proc_decomp_file_prefix` (character)

| | |
|---|---|
| Units | - |
| Description | Prefix of block mapping file *(hidden by default)* |
| Possible Values | Any valid filename *(default: graph.info.part.)* |
