# Chapter 7: Running the MPAS Non-hydrostatic Atmosphere Model

Given a CVT mesh, prepared using steps from [Chapter 4](04-preparing-meshes.md) as appropriate, this chapter describes the two main steps to running the MPAS-Atmosphere model: creating initial conditions and running the model itself. This chapter makes use of two MPAS cores, `init_atmosphere` and `atmosphere`, which are, respectively, used for initializing and running the non-hydrostatic atmospheric model. Sections 7.1 and 7.2 of this chapter describe the creation of idealized and real-data initial condition files using the `init_atmosphere` core. Section 7.3 describes the basic procedure of running the model itself.

Each section of this chapter follows a familiar pattern of compiling and executing MPAS model components, albeit using different cores depending on its intended use. The compilation will create either an initialization or a model executable, which are named, respectively, `init_atmosphere_model` and `atmosphere_model`. In general, an executable is run with `mpiexec` or `mpirun`, for example:

```
mpiexec -n 8 atmosphere_model
```

where `8` is the number of MPI tasks to be used. In any case where n > 1, there must exist a corresponding graph decomposition file, e.g., `graph.info.part.8`. For more on graph decomposition, see [Section 4.1](04-preparing-meshes.md#41-graph-partitioning-with-metis).

## 7.1 Creating Idealized ICs

There are several idealized test cases supported within the `init_atmosphere` model initialization core:

1. Jablonowski and Williamson baroclinic wave, no initial perturbation[^1]
2. Jablonowski and Williamson baroclinic wave, with initial perturbation
3. Jablonowski and Williamson baroclinic wave, with normal-mode perturbation
4. Squall line
5. Super-cell
6. Mountain wave

[^1]: Jablonowski, C. and D.L. Williamson, 2006, A baroclinic instability test case for atmospheric model dynamical cores, *QJRMS*, 132, 2943-2975. doi:10.1256/qj.06.12.

Creating idealized initial conditions is fairly straightforward, as no external data are required and the starting date/time is irrelevant to building the initial conditions file (hereafter referred to as `init.nc`) that will be used to run the model.

The following steps summarize the creation of `init.nc`:

- Include a `grid.nc` file, which contains the CVT mesh, in the working directory
- If running with more than one MPI task, include a `graph.info.part.*` file in the working directory ([Section 4.1](04-preparing-meshes.md#41-graph-partitioning-with-metis))
- Compile MPAS with the `init_atmosphere` core specified ([Section 3.3](03-building.md#33-compiling-mpas))
- Edit the `namelist.init_atmosphere` configuration file (described below)
- Edit the `streams.init_atmosphere` I/O configuration file (described below)
- Run `init_atmosphere_model` to create the initial condition file, `init.nc`

When the `init_atmosphere_model` executable is built, a default namelist, `namelist.init_atmosphere`, will have been created. A number of the namelist parameters found in `namelist.init_atmosphere` are irrelevant to creating idealized conditions and can be removed or ignored. The following table outlines the namelist parameters that are required; formal explanations for all namelist parameters can be found in [Appendix A](0A-init-namelist.md).

```
&nhyd_model
config_init_case        = 2                    ! a number between 1 and 6 corresponding to the
                                               ! cases listed at the beginning of this section
config_start_time       = '0000-01-01_00:00:00'  ! the starting time for the simulation
config_theta_adv_order  = 3                    ! advection order for theta
config_coef_3rd_order   = 0.25
/

&dimensions
config_nvertlevels      = 26                   ! the number of vertical levels to be used in the model
/

&decomposition
config_block_decomp_file_prefix = 'graph.info.part.'
                                               ! if running in parallel, needs to match the
                                               ! grid decomposition file prefix
/
```

After editing the `namelist.init_atmosphere` namelist file, the name of the input grid file, as well as the name of the initial condition file to be created, must be set in the XML I/O configuration file, `streams.init_atmosphere`. For a detailed description of the format of the XML I/O configuration file, refer to [Chapter 5](05-configuring-io.md). Specifically, the `filename_template` attribute must be set to the name of the grid file in the `"input"` stream definition, and the `filename_template` attribute must be set to name of the initial condition file to be created in the `"output"` stream definition.

## 7.2 Creating Real-Data ICs

The process of generating real-data initial conditions is similar to that of the idealized case described in the previous section, but is more involved as it requires interpolation of static geographic data (e.g., topography, land cover, soil category, etc.), surface fields such as soil temperature and SST, and the atmospheric initial conditions valid at a specific date and time. The static datasets are the same as those used by the WRF model, and the surface fields and atmospheric initial conditions can be obtained from, e.g., NCEP's GFS data using the WRF Pre-processing System (WPS).

Creating real-data initial conditions requires a single compilation of the `init_atmosphere` core, but the actual generation of the IC files will take place using two separate runs of the `init_atmosphere_model` program, where each of these two runs is described individually in the following sub-sections. While it is possible to condense the two real-data initialization steps into a single run, running each step separately will both improve clarity and save a significant amount of time when generating subsequent initial conditions (i.e., when making initial conditions using the same mesh but different starting times).

The first step, described in Section 7.2.1, is the interpolation of static fields onto the mesh to create a `static.nc` file. This step cannot be run in parallel and generally takes considerable time to complete; however, the fields being static, this step need only be run once for a particular mesh, regardless of the number of initial condition files that are ultimately created from the `static.nc` output file. Section 7.2.2 then describes the processing of the atmospheric initial conditions beginning with the `static.nc` file created in Section 7.2.1.

### 7.2.1 Static Fields

The generation of a `static.nc` file requires a set of static geographic data. A suitable dataset can be obtained from the WRF model's download page: <http://www2.mmm.ucar.edu/wrf/users/download/get_source.html>. These static data files should be downloaded to a directory, which will be specified in the `namelist.init_atmosphere` file prior to running this interpolation step. The result of this run will be the creation of a netCDF file (`static.nc`), which is used in Section 7.2.2 to create dynamic initial conditions. Note that `static.nc` can be generated once and then used repeatedly to generate initial condition files for different start times.

The following steps summarize the creation of `static.nc`:

- Download geographic data from the WRF download page (described above)
- Compile MPAS with the `init_atmosphere` core specified ([Section 3.3](03-building.md#33-compiling-mpas))
- Include a `grid.nc` file in the working directory
- Edit the `namelist.init_atmosphere` configuration file (described below)
- Edit the `streams.init_atmosphere` I/O configuration file (described below)
- Run `init_atmosphere_model` *with only one MPI task specified* to create `static.nc`

Note that it is critical for this step that the initialization core is run serially; afterward, however, the steps described in 7.2.2 may be run with more than one MPI task.

```
&nhyd_model
config_init_case        = 7                    ! must be 7, the real-data initialization case
/

&dimensions
config_nvertlevels      = 1                    ! the following dimensions should be set to 1 now,
config_nsoillevels      = 1                    ! and their values will become significant in §7.2.2
config_nfglevels        = 1
config_nfgsoillevels    = 1
/

&data_sources
config_geog_data_path        = '/path/to/WPS_GEOG/'
                                               ! absolute path to static files obtained from the
                                               ! WRF download page
config_landuse_data          = 'MODIFIED_IGBP_MODIS_NOAH'
                                               ! land use dataset selection
config_topo_data             = 'GMTED2010'     ! terrain dataset selection
config_vegfrac_data          = 'MODIS'         ! monthly vegetation fraction dataset selection
config_albedo_data           = 'MODIS'         ! monthly albedo dataset selection
config_maxsnowalbedo_data    = 'MODIS'         ! maximum snow albedo dataset selection
config_supersample_factor    = 1               ! MODIS supersampling factor, generally only
                                               ! needed for meshes with grid distance < 6 km
/

&preproc_stages
config_static_interp         = true            ! only the static_interp and native_gwd_static
config_native_gwd_static     = true            ! stages should be enabled
config_vertical_grid         = false
config_met_interp            = false
config_input_sst             = false
config_frac_seaice           = false
/
```

After editing the `namelist.init_atmosphere` namelist file, the name of the input CVT grid file, as well as the name of the static file to be created, must be set in the XML I/O configuration file, `streams.init_atmosphere`. Specifically, the `filename_template` attribute must be set to the name of the SCVT grid file in the `"input"` stream definition, and the `filename_template` attribute must be set to name of the static file to be created in the `"output"` stream definition.

### 7.2.2 Vertical Grid Generation and Initial Field Interpolation

The second step in creating a real-data initial conditions file (`init.nc`) is to generate a vertical grid, the parameters of which will be specified in the `namelist.init_atmosphere` file, and to obtain an initial conditions dataset and interpolate it onto the model grid. While initial conditions could ultimately be obtained from many different data sources, here we assume the use of WPS "intermediate" data files obtained from GFS data using the WPS `ungrib` program. Detailed instructions for building and running the WPS, and how to generate intermediate data files from GFS data, can be found in Chapter 3 of the WRF user guide: <http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.1/users_guide_chap3.html>.

The following steps summarize the creation of `init.nc`:

- Include a WPS intermediate data file in the working directory
- Include the `static.nc` file in the working directory (Section 7.2.1)
- If running in parallel, include a `graph.info.part.*` file in the working directory ([Section 4.1](04-preparing-meshes.md#41-graph-partitioning-with-metis))
- Edit the `namelist.init_atmosphere` configuration file (described below)
- Edit the `streams.init_atmosphere` I/O configuration file (described below)
- Run `init_atmosphere_model` to create `init.nc`

```
&nhyd_model
config_init_case        = 7                    ! must be 7
config_start_time       = '2010-10-23_00:00:00'  ! time to process first-guess data
config_theta_adv_order  = 3                    ! advection order for theta
config_coef_3rd_order   = 0.25
/

&dimensions
config_nvertlevels      = 55                   ! number of vertical levels to be used in MPAS
config_nsoillevels      = 4                    ! number of soil layers to be used in MPAS
config_nfglevels        = 38                   ! number of vertical levels in intermediate file
config_nfgsoillevels    = 4                    ! number of soil layers in intermediate file
/

&data_sources
config_met_prefix       = 'FILE'               ! the prefix of the intermediate file to be used
                                               ! for initial conditions
config_use_spechumd     = true                 ! if available, use specific humidity rather than
                                               ! relative humidity
/

&vertical_grid
config_ztop                    = 30000.0       ! model top height (m)
config_nsmterrain              = 1             ! number of smoothing passes for terrain
config_smooth_surfaces         = true          ! whether to smooth zeta surfaces
config_blend_boundary_terrain  = false         ! whether to blend terrain along domain boundaries;
                                               ! only for regional simulations as in Section 8.2
/

&preproc_stages
config_static_interp         = false
config_native_gwd_static     = false
config_vertical_grid         = true            ! only these three stages should be enabled
config_met_interp            = true
config_input_sst             = false
config_frac_seaice           = true
/

&decomposition
config_block_decomp_file_prefix = 'graph.info.part.'
                                               ! if running in parallel, needs to match the
                                               ! grid decomposition file prefix
/
```

After editing the `namelist.init_atmosphere` namelist file, the name of the static file, as well as the name of the initial condition file to be created, must be set in the XML I/O configuration file, `streams.init_atmosphere`. Specifically, the `filename_template` attribute must be set to the name of the static file in the `"input"` stream definition, and the `filename_template` attribute must be set to name of the initial condition file to be created in the `"output"` stream definition.

## 7.3 Running the Model

Having generated the model initial conditions, `init.nc`, as described in either Section 7.1 or 7.2, we have completed the prerequisites to run the model. The only step remaining before running the model itself is the configuration of `namelist.atmosphere`. When the `atmosphere` core is built, a default `namelist.atmosphere` namelist file will be automatically generated; this namelist can serve as a starting point for any modifications made following the steps below. This section will discuss both running the model from a cold start and restarting the model from some point in a previous run.

The following steps summarize running the model:

- Include an initial condition netCDF file (e.g., `init.nc`) in the working directory (Section 7.1, Section 7.2)
- (OPTIONAL) If the SST and sea-ice fields are to be periodically updated, include a surface netCDF file (e.g., `surface.nc`) in the working directory ([Section 8.1](08-model-options.md#81-periodic-sst-and-sea-ice-updates))
- (OPTIONAL) If running a regional simulation, include LBC netCDF files (e.g., `lbc.YYYY-MM-DD_HH.mm.ss.nc`) in the working directory ([Section 8.2](08-model-options.md#82-regional-simulation))
- If running in parallel, include a graph decomposition file in the working directory ([Section 4.1](04-preparing-meshes.md#41-graph-partitioning-with-metis))
- If the MPAS directory has not been cleaned since running initialization, run `make clean` with the `atmosphere` core specified
- Compile MPAS with the `atmosphere` core specified ([Section 3.3](03-building.md#33-compiling-mpas))
- Edit the default `namelist.atmosphere` configuration file (described below)
- Edit the `streams.atmosphere` I/O configuration file (described below)
- Run the `atmosphere_model` executable

Below is a list of variables in `namelist.atmosphere` that pertain to model timestepping, explicit horizontal diffusion, and model restarts. A number of namelist variables are not listed here (specifications for dynamical core configuration, physics parameters, etc.) and [Appendix B](0B-model-namelist.md) should be consulted for the purpose and acceptable values of these parameters.

```
&nhyd_model
config_dt               = 720.0               ! the model timestep; an appropriate value must
                                               ! be chosen relative to the grid cell spacing
config_start_time       = '2010-10-23_00:00:00'  ! the model start time corresponding to init.nc
config_run_duration     = '5_00:00:00'         ! the duration of the model run; for format
                                               ! rules, see Appendix B
config_len_disp         = 120000.0             ! the smallest cell-to-cell distance in the mesh,
                                               ! used for computing a dissipation length scale
/

&limited_area
config_apply_lbcs       = false                ! must be set to true only if running a regional
                                               ! simulation, as in Section 8.2
/

&decomposition
config_block_decomp_file_prefix = 'graph.info.part.'
                                               ! if running in parallel, must match the prefix
                                               ! of the graph decomposition file
/

&restart
config_do_restart       = false                ! if true, will select the appropriate restart.nc
                                               ! file generated from a previous run
/

&physics
config_sst_update       = true                 ! if updating sea-ice and SST with a surface.nc
                                               ! file, set to true, and edit the "surface" stream
                                               ! in the streams.atmosphere file accordingly
config_physics_suite    = 'mesoscale_reference'
/
```

When running the model from a cold start, `config_start_time` should match the time that was used when creating `init.nc`.

Configuration of model input and output is accomplished by editing the `streams.atmosphere` file. The following streams exist by default in the atmosphere core:

| Stream | Description |
|--------|-------------|
| `input` | The stream used to read model initial conditions for cold-start simulations |
| `restart` | The stream used to periodically write restart files during model integration, and to read initial conditions when performing a restart model run |
| `output` | The stream responsible for writing model prognostic and diagnostic fields to history files |
| `diagnostics` | The stream responsible for writing (mostly) 2-d diagnostic fields, typically at higher temporal frequency than the history files |
| `surface` | The stream used to read periodic updates of sea-ice and SST from a surface update file created as described in [Section 8.1](08-model-options.md#81-periodic-sst-and-sea-ice-updates) |
| `iau` | The stream used to read analysis increments for the Incremental Analysis Update (IAU) scheme |
| `lbc_in` | The stream used to read lateral boundary updates for limited-area simulations, described in [Section 8.2](08-model-options.md#82-regional-simulation) |

For more information on the options available in the XML I/O configuration file, users are referred to [Chapter 5](05-configuring-io.md).

During the course of a model run, restart files are created at an interval specified by the `output_interval` attribute in the definition of the `"restart"` stream. Running the model from a restart file is similar to running the model from `init.nc`. The required changes are that `config_do_restart` must be set to `true` and `config_start_time` must correspond to a restart file existing in the working directory.
