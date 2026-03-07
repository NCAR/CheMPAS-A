# Appendix B: Model Namelist Options

This appendix summarizes the complete set of namelist options available when running the MPAS non-hydrostatic atmosphere model. All date-time string specifications are of the form described at the beginning of [Appendix A](0A-init-namelist.md).

## B.1 nhyd_model

### `config_time_integration` (character)

| | |
|---|---|
| Units | - |
| Description | Time integration scheme *(hidden by default)* |
| Possible Values | `'SRK3'` *(default: SRK3)* |

### `config_time_integration_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Order for RK time integration |
| Possible Values | 2 or 3 *(default: 2)* |

### `config_dt` (real)

| | |
|---|---|
| Units | s |
| Description | Model time step, seconds |
| Possible Values | Positive real values *(default: 720.0)* |

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
| Description | Starting time for model simulation |
| Possible Values | `'YYYY-MM-DD_hh:mm:ss'` *(default: 2010-10-23_00:00:00)* |

### `config_stop_time` (character)

| | |
|---|---|
| Units | - |
| Description | Stopping time for model simulation *(hidden by default)* |
| Possible Values | `'YYYY-MM-DD_hh:mm:ss'` *(default: none)* |

### `config_run_duration` (character)

| | |
|---|---|
| Units | - |
| Description | Length of model simulation |
| Possible Values | `[DDD_]hh:mm:ss` *(default: 5_00:00:00)* |

### `config_split_dynamics_transport` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to super-cycle scalar transport |
| Possible Values | Logical values *(default: true)* |

### `config_number_of_sub_steps` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of acoustic steps per full RK step |
| Possible Values | Positive, even integer values, typically 2 or 6 depending on transport splitting *(default: 2)* |

### `config_dynamics_split_steps` (integer)

| | |
|---|---|
| Units | - |
| Description | When config_split_dynamics_transport = T, the number of RK steps per transport step |
| Possible Values | Positive integer values *(default: 3)* |

### `config_h_mom_eddy_visc2` (real)

| | |
|---|---|
| Units | m^2 s^-1 |
| Description | Laplacian eddy viscosity for horizontal diffusion of momentum *(hidden by default)* |
| Possible Values | Positive real values *(default: 0.0)* |

### `config_h_mom_eddy_visc4` (real)

| | |
|---|---|
| Units | m^4 s^-1 |
| Description | Biharmonic eddy hyper-viscosity for horizontal diffusion of momentum *(hidden by default)* |
| Possible Values | Positive real values *(default: 0.0)* |

### `config_v_mom_eddy_visc2` (real)

| | |
|---|---|
| Units | m^2 s^-1 |
| Description | Laplacian eddy viscosity for vertical diffusion of momentum *(hidden by default)* |
| Possible Values | Positive real values *(default: 0.0)* |

### `config_h_theta_eddy_visc2` (real)

| | |
|---|---|
| Units | m^2 s^-1 |
| Description | Laplacian eddy viscosity for horizontal diffusion of theta *(hidden by default)* |
| Possible Values | Positive real values *(default: 0.0)* |

### `config_h_theta_eddy_visc4` (real)

| | |
|---|---|
| Units | m^4 s^-1 |
| Description | Biharmonic eddy hyper-viscosity for horizontal diffusion of theta *(hidden by default)* |
| Possible Values | Positive real values *(default: 0.0)* |

### `config_v_theta_eddy_visc2` (real)

| | |
|---|---|
| Units | m^2 s^-1 |
| Description | Laplacian eddy viscosity for vertical diffusion of theta *(hidden by default)* |
| Possible Values | Positive real values *(default: 0.0)* |

### `config_horiz_mixing` (character)

| | |
|---|---|
| Units | - |
| Description | Formulation of horizontal mixing |
| Possible Values | `'2d_fixed'` or `'2d_smagorinsky'` *(default: 2d_smagorinsky)* |

### `config_len_disp` (real)

| | |
|---|---|
| Units | m |
| Description | Horizontal length scale, used by the Smagorinsky formulation of horizontal diffusion and by 3-d divergence damping *(hidden by default)* |
| Possible Values | Positive real values. A zero value implies that the length scale is prescribed by the nominalMinDc value in the input file. *(default: 0.0)* |

### `config_visc4_2dsmag` (real)

| | |
|---|---|
| Units | - |
| Description | Scaling coefficient of dx^3 to obtain biharmonic diffusion coefficient |
| Possible Values | Non-negative real values *(default: 0.05)* |

### `config_del4u_div_factor` (real)

| | |
|---|---|
| Units | - |
| Description | Scaling factor for the divergent component of the biharmonic u calculation *(hidden by default)* |
| Possible Values | Positive real values *(default: 10.0)* |

### `config_w_adv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Horizontal advection order for w *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_theta_adv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Horizontal advection order for theta *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_scalar_adv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Horizontal advection order for scalars *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_u_vadv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Vertical advection order for normal velocities (u) *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_w_vadv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Vertical advection order for w *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_theta_vadv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Vertical advection order for theta *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_scalar_vadv_order` (integer)

| | |
|---|---|
| Units | - |
| Description | Vertical advection order for scalars *(hidden by default)* |
| Possible Values | 2, 3, or 4 *(default: 3)* |

### `config_scalar_advection` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to advect scalar fields |
| Possible Values | .true. or .false. *(default: true)* |

### `config_positive_definite` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to enable positive-definite advection of scalars *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_monotonic` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to enable monotonic limiter in scalar advection |
| Possible Values | .true. or .false. *(default: true)* |

### `config_coef_3rd_order` (real)

| | |
|---|---|
| Units | - |
| Description | Upwinding coefficient in the 3rd order advection scheme |
| Possible Values | 0 <= config_coef_3rd_order <= 1 *(default: 0.25)* |

### `config_smagorinsky_coef` (real)

| | |
|---|---|
| Units | - |
| Description | Dimensionless empirical parameter relating the strain tensor to the eddy viscosity in the Smagorinsky turbulence model *(hidden by default)* |
| Possible Values | Real values typically in the range 0.1 to 0.4 *(default: 0.125)* |

### `config_mix_full` (logical)

| | |
|---|---|
| Units | - |
| Description | Mix full theta and u fields, or mix perturbation from initial state *(hidden by default)* |
| Possible Values | .true. or .false. *(default: true)* |

### `config_epssm` (real)

| | |
|---|---|
| Units | - |
| Description | Off-centering parameter for the vertically implicit acoustic integration |
| Possible Values | Positive real values *(default: 0.1)* |

### `config_smdiv` (real)

| | |
|---|---|
| Units | - |
| Description | 3-d divergence damping coefficient |
| Possible Values | Positive real values *(default: 0.1)* |

### `config_apvm_upwinding` (real)

| | |
|---|---|
| Units | - |
| Description | Amount of upwinding in APVM *(hidden by default)* |
| Possible Values | 0 <= config_apvm_upwinding <= 1 *(default: 0.5)* |

### `config_h_ScaleWithMesh` (logical)

| | |
|---|---|
| Units | - |
| Description | Scale eddy viscosities with mesh-density function for horizontal diffusion *(hidden by default)* |
| Possible Values | .true. or .false. *(default: true)* |

### `config_num_halos` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of halo layers for fields *(hidden by default)* |
| Possible Values | Integer values, typically 2 or 3; DO NOT CHANGE *(default: 2)* |

### `config_relax_zone_divdamp_coef` (real)

| | |
|---|---|
| Units | - |
| Description | Coefficient for the divergent component of the Laplacian filter of momentum in the relaxation zone *(hidden by default)* |
| Possible Values | Positive real values *(default: 6.0)* |

## B.2 damping

### `config_zd` (real)

| | |
|---|---|
| Units | m |
| Description | Height MSL to begin w-damping profile |
| Possible Values | Positive real values *(default: 22000.0)* |

### `config_xnutr` (real)

| | |
|---|---|
| Units | - |
| Description | Maximum w-damping coefficient at model top |
| Possible Values | 0 <= config_xnutr <= 1 *(default: 0.2)* |

### `config_mpas_cam_coef` (real)

| | |
|---|---|
| Units | - |
| Description | Coefficient for scaling the 2nd-order horizontal mixing in the mpas_cam absorbing layer *(hidden by default)* |
| Possible Values | 0 <= config_mpas_cam_coef <= 1, standard value is 0.2 *(default: 0.0)* |

### `config_number_cam_damping_levels` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of layers in which to apply CAM 2nd-order horizontal filter at top of model; viscosity linearly ramps to zero by layer number from the top *(hidden by default)* |
| Possible Values | Positive integer values *(default: 4)* |

### `config_rayleigh_damp_u` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to apply Rayleigh damping on horizontal velocity in the topmost model levels. The number of levels is specified by `config_number_rayleigh_damp_u_levels`, and the damping timescale is specified by `config_rayleigh_damp_u_timescale_days`. *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_rayleigh_damp_u_timescale_days` (real)

| | |
|---|---|
| Units | days |
| Description | Timescale, in days (86400 s), for the Rayleigh damping on horizontal velocity in the top-most model levels *(hidden by default)* |
| Possible Values | Positive real values *(default: 5.0)* |

### `config_number_rayleigh_damp_u_levels` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of layers in which to apply Rayleigh damping on horizontal velocity at top of model; damping linearly ramps to zero by layer number from the top *(hidden by default)* |
| Possible Values | Positive integer values *(default: 6)* |

## B.3 limited_area

### `config_apply_lbcs` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to apply lateral boundary conditions |
| Possible Values | true or false; this option must be set to true for limited-area simulations and false for global simulations *(default: false)* |

## B.4 io

### `config_restart_timestamp_name` (character)

| | |
|---|---|
| Units | - |
| Description | Filename used to store most recent restart time stamp *(hidden by default)* |
| Possible Values | Any valid filename *(default: restart_timestamp)* |

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

## B.5 decomposition

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

## B.6 restart

### `config_do_restart` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether this run of the model is to restart from a previous restart file or not |
| Possible Values | .true. or .false. *(default: false)* |

### `config_do_DAcycling` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to re-compute coupled fields (theta_m, rho_tilde, rho*u, etc.) from uncoupled fields when restarting the model; used for cycling DA experiments that analyze uncoupled fields in restart files *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

## B.7 printout

### `config_print_global_minmax_vel` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to print the global min/max of horizontal normal velocity and vertical velocity each timestep |
| Possible Values | .true. or .false. *(default: true)* |

### `config_print_detailed_minmax_vel` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to print the global min/max of horizontal normal velocity and vertical velocity each timestep, along with the location in the domain where those extrema occurred |
| Possible Values | .true. or .false. *(default: false)* |

### `config_print_global_minmax_sca` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to print the global min/max of scalar fields each timestep *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

## B.8 IAU

### `config_IAU_option` (character)

| | |
|---|---|
| Units | - |
| Description | Incremental Analysis Update scheme |
| Possible Values | `'off'` or `'on'`; `'off'` turns off IAU, `'on'` uses equal weighting of increments across the IAU window *(default: off)* |

### `config_IAU_window_length_s` (real)

| | |
|---|---|
| Units | s |
| Description | Length of window over which analysis increments are applied |
| Possible Values | Non-negative real values *(default: 21600.)* |

## B.9 assimilation

### `config_jedi_da` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether this run is within the JEDI data assimilation framework; used to add temperature and specific humidity as diagnostics *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

## B.10 development

### `config_halo_exch_method` (character)

| | |
|---|---|
| Units | - |
| Description | Method to use for exchanging halos *(hidden by default)* |
| Possible Values | `'mpas_dmpar'`, `'mpas_halo'` *(default: mpas_halo)* |

## B.11 physics

### `input_soil_data` (character)

| | |
|---|---|
| Units | - |
| Description | Input soil use classification *(hidden by default)* |
| Possible Values | `'STAS'` *(default: STAS)* |

### `input_soil_temperature_lag` (integer)

| | |
|---|---|
| Units | - |
| Description | Days over which the deep soil temperature is computed using skin temperature *(hidden by default)* |
| Possible Values | Positive integers *(default: 140)* |

### `num_soil_layers` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of soil layers in Noah land surface scheme *(hidden by default)* |
| Possible Values | Positive integers. For Noah LSM, must be set to 4. *(default: 4)* |

### `months` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of months per year *(hidden by default)* |
| Possible Values | Positive integer values. DO NOT CHANGE. *(default: 12)* |

### `noznlev` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of prescribed pressure levels for input climatological ozone volume mixing ratios *(hidden by default)* |
| Possible Values | Positive integers *(default: 59)* |

### `naerlev` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of prescribed pressure levels for input climatological aerosol mixing ratio *(hidden by default)* |
| Possible Values | Positive integers *(default: 29)* |

### `camdim1` (integer)

| | |
|---|---|
| Units | - |
| Description | Dimension of CAM radiation absorption save array *(hidden by default)* |
| Possible Values | Positive integers *(default: 4)* |

### `config_frac_seaice` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of fractional sea-ice *(hidden by default)* |
| Possible Values | .true. for sea-ice between 0 or 1; .false. for sea-ice equal to 0 or 1 (flag) *(default: true)* |

### `config_sfc_albedo` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of surface albedo *(hidden by default)* |
| Possible Values | .true. for climatologically varying surface albedo; .false. for fixed input data *(default: true)* |

### `config_sfc_snowalbedo` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of maximum surface albedo for snow *(hidden by default)* |
| Possible Values | .true. for geographical distribution; .false. for fixed input data *(default: true)* |

### `config_sst_update` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of sea-surface temperature |
| Possible Values | .true. for time-varying sea-surface temperatures; .false. otherwise *(default: false)* |

### `config_sstdiurn_update` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of diurnal cycle of sea-surface temperatures |
| Possible Values | .true. for applying a diurnal cycle to sea-surface temperatures; .false. otherwise *(default: false)* |

### `config_deepsoiltemp_update` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of deep soil temperatures |
| Possible Values | .true. for slowly time-varying deep soil temperatures; .false. otherwise *(default: false)* |

### `config_o3climatology` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for configuration of input ozone data in RRTMG long- and short-wave radiation *(hidden by default)* |
| Possible Values | .true. for using monthly-varying ozone data; .false. for using fixed vertical profile *(default: true)* |

### `config_microp_re` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for calculation of the effective radii for cloud water, cloud ice, and snow *(hidden by default)* |
| Possible Values | .true. for calculating effective radii; .false. for using defaults in RRTMG radiation *(default: false)* |

### `config_ysu_pblmix` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for turning on/off top-down, radiation-driven mixing *(hidden by default)* |
| Possible Values | .true. to turn on top-down radiation-driven mixing; .false. otherwise *(default: false)* |

### `config_urban_physics` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for turning on/off the urban physics parameterization *(hidden by default)* |
| Possible Values | .true. to turn on the urban physics; .false. otherwise *(default: false)* |

### `config_n_microp` (integer)

| | |
|---|---|
| Units | - |
| Description | Number of microphysics time-steps per physics time-steps *(hidden by default)* |
| Possible Values | Positive integers *(default: 1)* |

### `config_radtlw_interval` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between calls to parameterization of long-wave radiation |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: 00:30:00)* |

### `config_radtsw_interval` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between calls to parameterization of short-wave radiation |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: 00:30:00)* |

### `config_conv_interval` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between calls to parameterization of convection *(hidden by default)* |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: none)* |

### `config_pbl_interval` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between calls to parameterization of planetary boundary layer *(hidden by default)* |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: none)* |

### `config_camrad_abs_update` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between updates of absorption/emission coefficients in CAM radiation *(hidden by default)* |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: 06:00:00)* |

### `config_greeness_update` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between updates of greeness fraction *(hidden by default)* |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: 24:00:00)* |

### `config_bucket_update` (character)

| | |
|---|---|
| Units | - |
| Description | Time interval between updates of accumulated rain and radiation diagnostics |
| Possible Values | `'DD_HH:MM:SS'` or `'none'` *(default: none)* |

### `config_physics_suite` (character)

| | |
|---|---|
| Units | - |
| Description | Choice of physics suite |
| Possible Values | `'mesoscale_reference'`, `'convection_permitting'`, `'none'` *(default: mesoscale_reference)* |

### `config_microp_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for cloud microphysics schemes *(hidden by default)* |
| Possible Values | `'suite'`, `'mp_wsm6'`, `'mp_thompson'`, `'mp_thompson_aerosols'`, `'mp_kessler'`, `'off'` *(default: suite)* |

### `config_convection_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for convection schemes *(hidden by default)* |
| Possible Values | `'suite'`, `'cu_kain_fritsch'`, `'cu_tiedtke'`, `'cu_ntiedtke'`, `'cu_grell_freitas'`, `'off'` *(default: suite)* |

### `config_lsm_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for land-surface schemes *(hidden by default)* |
| Possible Values | `'suite'`, `'sf_noah'`, `'sf_noahmp'`, `'off'` *(default: suite)* |

### `config_pbl_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for planetary boundary layer schemes *(hidden by default)* |
| Possible Values | `'suite'`, `'bl_ysu'`, `'bl_mynn'`, `'off'` *(default: suite)* |

### `config_gwdo_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration of gravity wave drag over orography *(hidden by default)* |
| Possible Values | `'suite'`, `'bl_ysu_gwdo'`, `'bl_ugwp_gwdo'`, `'off'` *(default: suite)* |

### `config_ngw_scheme` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for non-stationary gravity wave drag scheme *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_knob_ugwp_tauamp` (real)

| | |
|---|---|
| Units | - |
| Description | Non-stationary gravity wave drag absolute momentum flux at launch level *(hidden by default)* |
| Possible Values | Non-negative real values *(default: 0.13e-3)* |

### `config_ugwp_diags` (logical)

| | |
|---|---|
| Units | - |
| Description | Logical for outputting UGWP drag suite diagnostic variables *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_radt_cld_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for calculation of horizontal cloud fraction *(hidden by default)* |
| Possible Values | `'suite'`, `'cld_fraction'`, `'cld_incidence'` *(default: suite)* |

### `config_radt_lw_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for long-wave radiation schemes *(hidden by default)* |
| Possible Values | `'suite'`, `'rrtmg_lw'`, `'cam_lw'`, `'off'` *(default: suite)* |

### `config_radt_sw_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for short-wave radiation schemes *(hidden by default)* |
| Possible Values | `'suite'`, `'rrtmg_sw'`, `'cam_sw'`, `'off'` *(default: suite)* |

### `config_sfclayer_scheme` (character)

| | |
|---|---|
| Units | - |
| Description | Configuration for surface layer scheme *(hidden by default)* |
| Possible Values | `'suite'`, `'sf_monin_obukhov'`, `'sf_mynn'`, `'off'` *(default: suite)* |

### `config_radt_cld_overlap` (character)

| | |
|---|---|
| Units | - |
| Description | Cloud overlapping option in the RRTMG LW and SW radiation schemes *(hidden by default)* |
| Possible Values | `'none'`, `'random'`, `'maximum_random'`, `'maximum'`, `'exponential'`, `'exponential_random'` *(default: maximum_random)* |

### `config_radt_cld_dcorrlen` (character)

| | |
|---|---|
| Units | - |
| Description | Decorrelation length for cloud overlapping option in the RRTMG LW and SW radiation schemes *(hidden by default)* |
| Possible Values | `'constant'`, `'latitude_varying'` *(default: constant)* |

### `config_gfconv_closure_deep` (integer)

| | |
|---|---|
| Units | - |
| Description | Closure option for deep convection in Grell-Freitas convection scheme *(hidden by default)* |
| Possible Values | 0 *(default: 0)* |

### `config_gfconv_closure_shallow` (integer)

| | |
|---|---|
| Units | - |
| Description | Closure option for shallow convection in Grell-Freitas convection scheme *(hidden by default)* |
| Possible Values | 8 *(default: 8)* |

### `config_mynn_tkeadvect` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to do MYNN TKE advection *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_mynn_tkebudget` (integer)

| | |
|---|---|
| Units | - |
| Description | Option to add the MYNN TKE budget terms to output (=1 to enable) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_mynn_closure` (real)

| | |
|---|---|
| Units | - |
| Description | 2.5 level or 3.0 level for the MYNN TKE turbulence closure *(hidden by default)* |
| Possible Values | 2.0, 2.5, 3.0 *(default: 2.5)* |

### `config_mynn_cloudpdf` (integer)

| | |
|---|---|
| Units | - |
| Description | Switch to different cloud PDFs to represent subgrid clouds in the MYNN PBL scheme *(hidden by default)* |
| Possible Values | 0, 1, 2 *(default: 2)* |

### `config_mynn_mixlength` (integer)

| | |
|---|---|
| Units | - |
| Description | Mixing length in the MYNN PBL scheme: 0=Nakanishi and Niino 2009; 1=RAP/HRRR; 2=Ito et al. 2015 *(hidden by default)* |
| Possible Values | 0, 1, 2 *(default: 2)* |

### `config_mynn_stfunc` (integer)

| | |
|---|---|
| Units | - |
| Description | Option to switch flux profile relationship for surface (0:Dyer-Hicks, 1:Cheng-Brustaert) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 1)* |

### `config_mynn_topdown` (integer)

| | |
|---|---|
| Units | - |
| Description | Option to add top-down diffusion driven by cloud-top radiative cooling *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_mynn_scaleaware` (integer)

| | |
|---|---|
| Units | - |
| Description | Option to turn on the scale-aware option *(hidden by default)* |
| Possible Values | 0, 1 *(default: 1)* |

### `config_mynn_dheat_opt` (integer)

| | |
|---|---|
| Units | - |
| Description | Option to activate heating due to dissipation of TKE *(hidden by default)* |
| Possible Values | 0, 1 *(default: 1)* |

### `config_mynn_edmf` (integer)

| | |
|---|---|
| Units | - |
| Description | Activate the EDMF option in the MYNN PBL scheme (=1) or not (=0) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 1)* |

### `config_mynn_edmf_dd` (integer)

| | |
|---|---|
| Units | - |
| Description | Activate the EDMF downdraft option in the MYNN PBL scheme (=1) or not (=0) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_mynn_edmf_mom` (integer)

| | |
|---|---|
| Units | - |
| Description | Activate momentum transport with the EDMF option of the MYNN PBL scheme (=1) or not (=0) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_mynn_edmf_tke` (integer)

| | |
|---|---|
| Units | - |
| Description | Activate TKE transport with the EDMF option of the MYNN PBL scheme (=1) or not (=0) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_mynn_edmf_output` (integer)

| | |
|---|---|
| Units | - |
| Description | Suppress (=0) the allocation of variables needed in the EDMF option of the MYNN PBL scheme *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_mynn_mixscalars` (integer)

| | |
|---|---|
| Units | - |
| Description | Activate mixing of scalars in the MYNN PBL scheme (=1) or not (=0) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 1)* |

### `config_mynn_mixclouds` (integer)

| | |
|---|---|
| Units | - |
| Description | Activate mixing of cloud condensates in the MYNN PBL scheme (=1) *(hidden by default)* |
| Possible Values | 0, 1 *(default: 1)* |

### `config_mynn_mixqt` (integer)

| | |
|---|---|
| Units | - |
| Description | =0: activate mixing of moisture species separately; =1: activate mixing of total water in the MYNN PBL scheme *(hidden by default)* |
| Possible Values | 0, 1 *(default: 0)* |

### `config_bucket_radt` (real)

| | |
|---|---|
| Units | - |
| Description | Threshold above which accumulated radiation diagnostics are reset *(hidden by default)* |
| Possible Values | Positive real values *(default: 1.0e9)* |

### `config_bucket_rainc` (real)

| | |
|---|---|
| Units | - |
| Description | Threshold above which the accumulated convective precipitation is reset *(hidden by default)* |
| Possible Values | Positive real values *(default: 100.0)* |

### `config_bucket_rainnc` (real)

| | |
|---|---|
| Units | - |
| Description | Threshold above which the accumulated grid-scale precipitation is reset *(hidden by default)* |
| Possible Values | Positive real values *(default: 100.0)* |

### `config_oml1d` (logical)

| | |
|---|---|
| Units | - |
| Description | Whether to activate the 1-d ocean mixed-layer model *(hidden by default)* |
| Possible Values | .true. or .false. *(default: false)* |

### `config_oml_hml0` (real)

| | |
|---|---|
| Units | m |
| Description | If greater than zero, provides the constant, initial mixed-layer depth; if zero, initial mixed-layer depth will be taken from the `oml_initial` field *(hidden by default)* |
| Possible Values | Non-negative real values *(default: 30.0)* |

### `config_oml_gamma` (real)

| | |
|---|---|
| Units | K m^-1 |
| Description | Deep water lapse rate in 1-d OML model *(hidden by default)* |
| Possible Values | Real values *(default: 0.14)* |

### `config_oml_relaxation_time` (real)

| | |
|---|---|
| Units | s |
| Description | Relaxation time to initial values in 1-d OML *(hidden by default)* |
| Possible Values | Non-negative real values *(default: 864000.)* |

## B.12 soundings

### `config_sounding_interval` (character)

| | |
|---|---|
| Units | - |
| Description | Interval between writing of soundings. A value of `'none'` disables writing of soundings. |
| Possible Values | `'[DDD_]hh:mm:ss'` or `'none'` *(default: none)* |

## B.13 physics_lsm_noahmp

### `config_noahmp_iopt_dveg` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for dynamic vegetation |
| Possible Values | 1 to 6 *(default: 4)* |

### `config_noahmp_iopt_crs` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for canopy stomatal resistance |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_btr` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for soil moisture in stomatal resistance |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_runsrf` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for surface runoff |
| Possible Values | 0 or 1 *(default: 3)* |

### `config_noahmp_iopt_runsub` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for sub-surface option |
| Possible Values | 0 or 1 *(default: 3)* |

### `config_noahmp_iopt_sfc` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for surface drag |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_frz` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for supercooled water |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_inf` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for frozen soil |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_rad` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for radiative transfer |
| Possible Values | 0 or 1 *(default: 3)* |

### `config_noahmp_iopt_alb` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for snow albedo |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_snf` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for precipitation partition option |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_tksno` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for snow thermal conductivity |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_tbot` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for lower boundary of soil temperature |
| Possible Values | 0 or 1 *(default: 2)* |

### `config_noahmp_iopt_stc` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for soil or snow time scheme |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_gla` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for glacier |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_rsf` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for subsurface resistance |
| Possible Values | 0 or 1 *(default: 4)* |

### `config_noahmp_iopt_soil` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for soil data |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_pedo` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for pedo transfer |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_crop` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for crop |
| Possible Values | 0 or 1 *(default: 0)* |

### `config_noahmp_iopt_irr` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for irrigation |
| Possible Values | 0 or 1 *(default: 0)* |

### `config_noahmp_iopt_irrm` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for irrigation method |
| Possible Values | 0 or 1 *(default: 0)* |

### `config_noahmp_iopt_infdv` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for DVIC infiltration |
| Possible Values | 0 or 1 *(default: 1)* |

### `config_noahmp_iopt_tdrn` (integer)

| | |
|---|---|
| Units | - |
| Description | Option for drainage option |
| Possible Values | 0 or 1 *(default: 0)* |
