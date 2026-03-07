# Chapter 4: Preparing Meshes

This chapter describes the steps used to prepare SCVT meshes for use in MPAS-Atmosphere. For quasi-uniform meshes, very little preparation is actually needed, and generally, one only needs to prepare mesh decomposition files -- files that describe the decomposition of the SCVT mesh across processors -- when running MPAS-Atmosphere using multiple MPI tasks. The procedure for creating these mesh decomposition files is described in the first section.

For variable-resolution SCVT meshes, the area of mesh refinement may be rotated to any part of the sphere using a program, `grid_rotate`, described in the second section. This utility program may be obtained from the MPAS-Atmosphere download page.

When running limited-area simulations, the limited-area mesh may be produced by subsetting the elements in an existing SCVT mesh as described in the third section.

## 4.1 Graph Partitioning with METIS

Before MPAS can be run in parallel, a mesh decomposition file with an appropriate number of partitions (equal to the number of MPI tasks that will be used) is required. A limited number of mesh decomposition files, named `graph.info.part.*`, are provided with each mesh, as is the mesh connectivity file, named `graph.info`. If the number of MPI tasks to be used when running MPAS matches one of the pre-computed decomposition files, then there is no need to run METIS.

In order to create new mesh decomposition files for some particular number of MPI tasks, only the `graph.info` file is required. The currently supported method for partitioning a `graph.info` file uses the METIS software (<http://glaros.dtc.umn.edu/gkhome/views/metis>). The serial graph partitioning program, METIS (rather than ParMETIS or hMETIS) should be sufficient for quickly partitioning any mesh usable by MPAS.

After installing METIS, a `graph.info` file may be partitioned into *N* partitions by running:

```
gpmetis -minconn -contig -niter=200 graph.info N
```

where *N* is the required number of partitions. The resulting file, `graph.info.part.N`, can then be copied into the MPAS run directory before running the model with *N* MPI tasks.

## 4.2 Relocating Refinement Regions on the Sphere

The purpose of the `grid_rotate` program is simply to rotate an MPAS mesh file, moving a refinement region from one geographic location to another, so that the mesh can be re-used for different applications. This utility was developed out of the need to save computational resources, since generating an SCVT -- particularly one with a large number of generating points or a high degree of refinement -- can take considerable time.

To build the `grid_rotate` program, the Makefile should first be edited to set the Fortran compiler to be used; if the netCDF installation pointed to by the `NETCDF` environment variable was built with a separate Fortran interface library, it will also be necessary to add `-lnetcdff` just before `-lnetcdf` in the Makefile. After editing the Makefile, running `make` should result in a `grid_rotate` executable file.

Besides the MPAS grid file to be rotated, `grid_rotate` requires a namelist file, `namelist.input`, which specifies the rotation to be applied to the mesh. The namelist variables are summarized in the table below:

| Variable | Description |
|----------|-------------|
| `config_original_latitude_degrees` | Original latitude of any point on the sphere |
| `config_original_longitude_degrees` | Original longitude of any point on the sphere |
| `config_new_latitude_degrees` | Latitude to which the original point should be shifted |
| `config_new_longitude_degrees` | Longitude to which the original point should be shifted |
| `config_birdseye_rotation_counter_clockwise_degrees` | Rotation about a vector from the sphere center through the original point |

Essentially, one chooses any point on the sphere, decides where that point should be shifted to, and specifies any change to the orientation (i.e., rotation) of the mesh about that point.

Having set the rotation parameters in the `namelist.input` file, the `grid_rotate` program should be run with two command-line options specifying the original grid file name and the name of the rotated grid file to be produced, e.g.,

```
grid_rotate grid.nc grid_SE_Asia_refinement.nc
```

The original grid file will not be altered, and a new, rotated grid file will be created. The NCL script `mesh.ncl` may be used to plot either of the original or rotated grid files after suitable setting the name of the grid file in the script.

> *Note: The grid\_rotate program initializes the new, rotated grid file to a copy of the original grid file. If the original grid file has only read permission (i.e., no write permission), then so will the copy, and consequently, the grid\_rotate program will fail when attempting to update the fields in the copy.*

## 4.3 Creating Limited-Area SCVT Meshes

The process of creating a limited-area (regional) mesh for MPAS-Atmosphere involves selecting any existing mesh -- either a quasi-uniform mesh, or a variable-resolution mesh that has been rotated as in [Section 4.2](#42-relocating-refinement-regions-on-the-sphere) -- describing the geographical region to be extracted from that mesh, and running the limited-area Python program to extract all cells, edges, and vertices in the designated region. The result is a new netCDF mesh file that can be used to make a limited-area simulation as described in [Section 8.2](08-model-options.md#82-regional-simulation).

The limited-area Python program may be obtained from the MPAS-Atmosphere download page. Although the set of required Python packages may change over time, the program currently requires the `numpy` and `netCDF4` packages, in addition to other standard packages.
