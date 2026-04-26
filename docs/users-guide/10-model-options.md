# Chapter 10: Model Options

Beyond the basic process of running a global simulation with standard output files outlined in [Chapter 9](09-running.md), the CheMPAS-A model provides several options that can be described in terms of variations on the basic simulation workflow. In the sections that follow, major model options are described in terms of the deviation from the basic global simulation process.

## 10.1 Periodic SST and Sea-ice Updates

The stand-alone CheMPAS-A model is not coupled to fully prognostic ocean or sea-ice models, and accordingly, the model SST and sea-ice fraction fields will in general not change over the course of a simulation. For simulations shorter than a few days, invariant SST and sea-ice fraction fields will generally not be problematic. However, for longer model simulations, it is generally recommended to periodically update the SST and sea-ice fields from an external file.

The surface data to be used for periodic SST and sea-ice updates could originate from any number of sources, though the most straightforward way to obtain a dataset in a usable format is to process GRIB data (e.g., GFS GRIB data) with the *ungrib* program of the WRF model's pre-processing system (WPS). Detailed instructions for building and running the WPS, and the process of generating intermediate data files from GFS data, can be found in Chapter 3 of the WRF User Guide: <http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.1/users_guide_chap3.html>.

The following steps summarize the generation of an SST and sea-ice update file, `surface.nc`, using the `init_atmosphere_model` program:

- Include surface data intermediate files in the working directory
- Include a `static.nc` file in the working directory ([Section 9.2.1](09-running.md#921-static-fields))
- If running in parallel, include a `graph.info.part.*` in the working directory ([Section 4.1](04-preparing-meshes.md#41-graph-partitioning-with-metis))
- Edit the `namelist.init_atmosphere` configuration file (see below)
- Edit the `streams.init_atmosphere` I/O configuration file (described below)
- Run `init_atmosphere_model` to create `surface.nc`

```
&nhyd_model
config_init_case        = 8                    ! must be 8, the surface field initialization case
config_start_time       = '2010-10-23_00:00:00'  ! time to begin processing surface data
config_stop_time        = '2010-10-30_00:00:00'  ! time to end processing surface data
/

&data_sources
config_sfc_prefix       = 'SST'                ! the prefix of the intermediate data files
                                               ! containing SST and sea-ice
config_fg_interval      = 86400                ! interval between intermediate files to use
                                               ! for SST and sea-ice
/

&preproc_stages
config_static_interp         = false           ! only the input_sst and frac_seaice stages
config_native_gwd_static     = false           ! should be enabled
config_vertical_grid         = false
config_met_interp            = false
config_input_sst             = true
config_frac_seaice           = true
/

&decomposition
config_block_decomp_file_prefix = 'graph.info.part.'
                                               ! if running in parallel, needs to match the
                                               ! grid decomposition file prefix
/
```

After editing the `namelist.init_atmosphere` namelist file, the name of the static file, as well as the name of the surface update file to be created, must be set in the XML I/O configuration file, `streams.init_atmosphere`. Specifically, the `filename_template` attribute must be set to the name of the static file in the `"input"` stream definition, and the `filename_template` attribute must be set to name of the surface update file to be created in the `"surface"` stream definition. *Also, for the "surface" stream, ensure that the `output_interval` attribute is set to the interval at which the surface intermediate files are provided.*

## 10.2 Regional Simulation

New in MPAS v7.0 is the capability to run simulations over regional domains on the surface of the sphere. Setting up and running a limited-area simulation requires as a starting point a limited-area SCVT mesh, described in [Section 4.3](04-preparing-meshes.md#43-creating-limited-area-scvt-meshes). Given a limited-area mesh, the key differences from a global simulation are:

- The blending of the MPAS terrain field with the "first-guess" terrain data along the boundaries of the limited-area domain;
- The generation of a set of files containing lateral boundary conditions (LBCs); and
- The application of LBCs during the model integration.

### Terrain blending

The first of these differences -- the blending of terrain data -- takes place when generating the limited-area initial conditions. Limited-area initial conditions are prepared as in [Section 9.2.2](09-running.md#922-vertical-grid-generation-and-initial-field-interpolation), except that the `config_blend_bdy_terrain` option should be set to `true` in the `namelist.init_atmosphere` file. This option instructs the `init_atmosphere_model` program to perform averaging of the model terrain field from the `static.nc` file with the terrain field from the atmospheric initial conditions dataset along the lateral boundaries of the mesh.

### LBC generation

The second difference -- the generation of LBC files -- requires running the `init_atmosphere_model` program one additional time, with namelist options set as described below.

```
&nhyd_model
config_init_case        = 9                    ! the LBCs processing case
config_start_time       = '2010-10-23_00:00:00'  ! time to begin processing LBC data
config_stop_time        = '2010-10-30_00:00:00'  ! time to end processing LBC data
/

&dimensions
config_nfglevels        = 38                   ! number of vertical levels in intermediate file
/

&data_sources
config_met_prefix       = 'GFS'                ! the prefix of intermediate data files to be
                                               ! used for LBCs
config_fg_interval      = 10800                ! interval between intermediate files
/

&decomposition
config_block_decomp_file_prefix = 'graph.info.part.'
                                               ! if running in parallel, needs to match the
                                               ! grid decomposition file prefix
/
```

When running the LBC processing case, the `output_interval` for the `"lbc"` stream in the `streams.init_atmosphere` file must be set to match the value of `config_fg_interval` in the `namelist.init_atmosphere` file. Additionally, the file to be read by the `"input"` stream must contain vertical grid information; typically, the model initial-conditions file can be used as the source for the `"input"` stream. The end result of the LBC processing case should be a set of netCDF files containing LBCs for the model integration.

### LBC application

The final difference -- application of LBCs during the model integration -- simply requires setting `config_apply_lbcs` to `true` in the model's `namelist.atmosphere` file, as well as setting the `input_interval` for the `"lbc_in"` stream in the `streams.atmosphere` file to match the interval at which the LBC netCDF files were produced.

## 10.3 Separate Stream for Invariant Fields

By default, the CheMPAS-A model reads time-invariant fields (e.g., `latCell`, `lonCell`, `areaCell`, `zgrid`, `zz`, etc.) from the `"input"` and `"restart"` streams (for cold-start and restart runs, respectively), and it writes time-invariant fields to the `"restart"` stream. In the case of large ensembles, the time-invariant fields that are replicated in the restart files for all ensemble members can account for a substantial amount of storage. In principle, since these time-invariant fields do not change in time or across ensemble members, only one copy of these fields needs to be stored.

MPAS-Atmosphere v8.1.0 introduces a capability to omit time-invariant fields from model restart files. When the model restarts, a new `"invariant"` stream may be used to read time-invariant fields from a separate file, and many ensemble members can share this file.

In order to make use of the new `"invariant"` stream, several changes to the standard CheMPAS-A workflow are needed.

### 10.3.1 Preparing an Invariant File

Through the use of the `init_atmosphere_model` program, a file containing all required time-invariant fields must be prepared. Of course, since the model initial conditions (typically referred to as the `init.nc` file) file contains time-invariant fields, the initial conditions file from any ensemble member may be used.

If a file containing purely time-invariant fields is desired, the output from the `init_atmosphere_model` program after the following pre-processing stages have been run will suffice:

- `config_static_interp = true`
- `config_native_gwd_static = true`
- `config_vertical_grid = true`

Note that these pre-processing stages do not need to be run all at once. It is possible, for example, to first produce a `static.nc` file using the first two of these pre-processing stages, and to then produce an invariant file (e.g., `invariant.nc`) by running the vertical grid generation stage using the `static.nc` file as input.

### 10.3.2 Activating the Invariant Stream

When running the model itself (`atmosphere_model`), the use of the new invariant stream may be activated by simply defining the `"invariant"` immutable stream in the `streams.atmosphere` file as follows:

```xml
<immutable_stream name="invariant"
                  type="input"
                  filename_template="invariant.nc"
                  input_interval="initial_only" />
```

In the definition of the `"invariant"` stream, the `filename_template` attribute should be set to the actual name of the invariant file.

By the existence of the `"invariant"` stream in the `streams.atmosphere` file, the model will omit all time-invariant fields from any restart files that are written. When the model restarts, all time-invariant fields will be read from the `"invariant"` stream rather than from the `"restart"` stream.
