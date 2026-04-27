# Chapter 11: Visualization

Since the MPAS input and output files are in netCDF format, a wide variety of software tools may be used to manipulate and visualize fields in these files. As a starting point, upstream MPAS-Atmosphere provides several NCL[^1] scripts for making basic mesh, contour, cell-filled, and cross-section plots. Each of these scripts reads the name of the file from which fields should be plotted from the environment variable `FNAME`, and all but the mesh-plotting script read the time frame to be plotted from the environment variable `T`. To plot a field from the first frame (indexed from 0) of the file `output.2010-10-23_00:00:00.nc`, for example, one would set the following environment variables:

```
setenv FNAME output.2010-10-23_00:00:00.nc
setenv T 0
```

before running one of the scripts. In general, the specific field to be plotted from the netCDF file must be set within a script before running that script.

CheMPAS-A also includes Python-based chemistry plotting scripts in the repository's `scripts/` directory, including `plot_chemistry.py` for ABBA-style tracer output and `plot_chemistry_profiles.py` for Chapman/NOx profile diagnostics. These scripts are documented in the CheMPAS-A visualization guide under `docs/chempas/guides/VISUALIZE.md`.

[^1]: NCAR Command Language; http://ncl.ucar.edu

## 11.1 Meshes

A plot showing just an MPAS SCVT mesh can be produced using the `atm_mesh.ncl` script. This script reads a subset of the mesh description fields in [Appendix C](0C-grid-description.md) and uses this information to draw the SCVT mesh over a color-filled map background. Parameters in the script can be used to control the type of map projection (e.g., orthographic, cylindrical equidistant, etc.), the colors used to fill land and water points, and the widths of lines used for the Voronoi cells.

## 11.2 Horizontal Contour Plots

Contour plots of horizontal fields can be produced with the `atm_contours.ncl` script. The particular field to be plotted is set in the script and can in principle be drawn on any horizontal surface (e.g., a constant pressure surface, a constant height surface, sea-level, etc.) if suitable vertical interpolation code is added to the script. Not shown in the figure are horizontal wind vectors, which can also be added to the plot using example code provided in the script.

## 11.3 Horizontal Cell-Filled Plots

For visualizing horizontal fields on their native SCVT grid, the `atm_cells.ncl` script may be used to produce horizontal cell-filled plots. This script draws each MPAS grid cell as a polygon colored according to the value of the field in that cell; the color scale is automatically chosen based on the range of the field and the default NCL color table, though other color tables can be selected instead.

## 11.4 Vertical Cross-Sections

Vertical cross-sections of fields can be created using the `atm_xsec.ncl` script. Before running this script, a starting point and an ending point for the cross section must be given as latitude-longitude pairs near the top of the script, and the number of points along the cross section should be specified. The script evenly distributes the specified number of points along the shortest great-circle arc from the starting point to the ending point, and for each point, the script uses values from the grid cell containing that point (i.e., a nearest-neighbor interpolation to the horizontal cross-section points is performed); no vertical interpolation is performed, and the thicknesses and vertical heights of cells are all drawn according to the MPAS vertical grid.
