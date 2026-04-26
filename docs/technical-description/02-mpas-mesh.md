# Chapter 2: MPAS mesh

MPAS-A integrates the equations of motion that are spatially discretized employing a centroidal Voronoi mesh for its horizontal discretization with a hybrid vertical coordinate based on geometric height. For applications on the sphere, the horizontal mesh is defined using great-circle arcs as opposed to using a planer representation, thus the distances and areas are consistent in all spherical applications. In this chapter we briefly describe the centroidal Voronoi horizontal mesh and the vertical coordinate, along with the C-grid staggering of the MPAS-A prognostic variables. For detailed information on the centroidal Voronoi horizontal mesh readers are advised to consult Ringler et al. (2008) and references therein.

## 2.1 Horizontal Voronoi Mesh

**[Figure 2.1: A portion of an MPAS-Atmosphere horizontal mesh with C-grid staggered velocities. To be added next session.]**

A schematic of a representative portion of the horizontal MPAS-A centroidal Voronoi mesh is given in Figure 2.1. MPAS-A meshes are mostly comprised of hexagons, but may contain some pentagons and possibly heptagons, and mesh cells are referred to by their central points that are located at the centroid (center of mass) of the cells. In Figure 2.1 the cells (cell centers) are $A$, $B$ and $C$. The dual of the centroidal Voronoi mesh is the Delaunay triangular mesh, and a Delaunay triangle $ABC$ is shown in the figure. There are three properties of this centroidal Voronoi mesh that are critical to the MPAS-A discretization:

(a) The centroidal property: The cell centers are centroids of the Voronoi cell, and any point in the cell is closer to the cell center than it is to any other cell center.

(b) The lines connecting cell centers (edges of the dual Delaunay triangular mesh) are bisected by the cell edges: This follows from (a).

(c) The lines connecting cell centers are perpendicular to cell edges at their intersection point. This also follows from (a).

Specifically, the fact that the horizontal line connecting the cell centers is orthogonal to the cell edge, and that the prognostic orthogonal velocity is coincident with that line at the edge, greatly facilitates the discretization of the momentum equation, the transport of the cell-centered variables, and many of the filters, as described in Chapters 4 and 5.

### 2.1.1 Points on the MPAS horizontal mesh

**[Figure 2.2: Horizontal mesh points and lengths used in MPAS. To be added next session.]**

The locations of the three points that are used to define the mesh are given in Figure 2.2:

1. The cell-center points, which are the centroids of the Voronoi cells. Note that they are also the vertices of the dual Delauney triangular mesh cell (see Figure 2.1).

2. The vertex points, which are the vertices of the Voronoi cells and the cell centers for the dual Delauney triangular mesh cell.

3. The edge points, which are located where the lines connecting the cell centers intersect the cell edges (the lines connecting the cell vertices). The edge points bisect the lines connecting the cell centers. However, the edge points do not necessary bisect the cell edges.

Two horizontal meshes are supported in MPAS-Atmosphere, meshes on Cartesian planes and meshes on the sphere. Points are defined on the sphere by their locations in an $(x, y, z)$ coordinate system with its origin at the center of the sphere. MPAS also keeps the latitude and longitude locations of the points on the sphere. The $(x, y, z)$ locations of the points as a function of the latitude $\phi$ and the longitude $\lambda$ are given by

$$
\begin{aligned}
x &= r\cos(\lambda)\cos(\phi), \\
y &= r\sin(\lambda)\cos(\phi), \\
z &= r\sin(\phi),
\end{aligned}
$$

where $r$ is the radius of the sphere. Conversely,

$$
\phi = \arcsin\!\left(\frac{z}{r}\right), \qquad \lambda = \arctan\!\left(\frac{y}{x}\right)
$$

The Cartesian plane is defined by $(x, y, z = 0)$, and the latitudes and longitudes are set to zero unless an $f$-plane approximation is used.

:::{admonition} MPAS code
:class: note

Information about the MPAS mesh is defined in the MPAS `Registry.xml` file that an be found at `MPAS/src/core_atmosphere/Registry.xml` for the atmospheric model and in `MPAS/src/core_init_atmosphere/Registry.xml` for the atmospheric initialization utility. There are `nCells` cell-center points (and hence cells) in an MPAS mesh, and the locations of these points are defined in the Fortran arrays `xCell(nCells)`, `yCell(nCells)`, `zCell(nCells)` and `latCell(nCells)` and `lonCell(nCells)`. Similarly, there are `nEdges` edges and edge points defined in the arrays `xEdge`, `yEdge`, `zEdge` and `latEdge` and `lonEdge`. There are `nVertices` vertex points (and cells in the dual Delauney triangular mesh) defined in the arrays `xVertex`, `yVertex`, `zVertex`, and `latVertex` and `lonVertex`. The $(x,y,z)$ points on the sphere have units of meters.
:::

:::{admonition} Note — unstructured MPAS mesh
:class: note

The MPAS mesh is unstructured, and the order of the points on the mesh is arbitrary. For example, cell 42 will be located at $(x,y,z)$ points `xCell(42)`, `yCell(42)`, `zCell(42)` and latitude and longitude `latCell(42)` and `lonCell(42)`. Likewise for edges and vertices. Neighbor cells (that share an edge) may have any other index between 1 and `nCells`. How to find neighbors, i.e. the connectivity of the mesh, is discussed in section 2.1.3.
:::

### 2.1.2 Lines and arcs on the MPAS horizontal mesh

There are two classes of lines that are define the MPAS horizontal mesh - the lines connecting the cell centers defining the dual Delauney triangular mesh, and the lines connecting vertices that define the cell edges. These lines are depicted in Figure 2.2, and are defined by the points they connect and their length. The points and their definition in the sphere and plane have been discussed in section 2.1.1. The line lengths on the sphere are the lengths of the great circle arcs on the sphere connecting the points. The line lengths on the plane are the lengths of the straight lines connecting the points.

:::{admonition} MPAS code
:class: note

In `MPAS/src/core_atmosphere/Registry.xml`, the line lengths for lines connecting cell centers across an edge are stored in the MPAS array `dcEdge(nEdges)`, and the cells at the endpoints of the line are stored in the array `cellsOnEdge(2,nEdges)`. The line lengths for lines connecting vertices (edge lengths) are stored in the MPAS array `dvEdge(nEdges)`, and the vertices at the endpoints of the edge are stored in the array `verticesOnEdge(2,nEdges)`. The lengths have units of meters.
:::

### 2.1.3 Connectivity

**Components of a cell**

Cells are composed of a cell center, edges, and vertices. The edges are bounded by vertices. As an illustration, Figure 2.3 shows the cells surrounding cell 42 in the left panel, the edges in the center panel, and the vertices in the right panel. MPAS contains arrays that contain this information.

**[Figure 2.3: Cell, edge and vertex neighbors of a cell. Cell indices are given in black, edge indices are red and vertex indices are given in blue. To be added next session.]**

:::{admonition} MPAS code
:class: note

Following Figure 2.3, the array containing the neighbors of cell 42 has the values `cellsOnCell(1:6,42) = (30, 58, 6, 104, 17, 51)`. By convention the order is counterclockwise around the cell. Likewise, the edges around the cell are given as `edgesOnCell(1:6,42) = (17, 14, 86, 54, 67, 201)`. By convention the edges are taken in the same order and start with the edge shared by the center cell and the first cell in the `cellsOnCell` list, i.e. `cellsOnCell(1,42) = 30`. Lastly, vertices for a cell are given as `verticesOnCell(1:6,42) = (13, 11, 55, 46, 112, 81)` and follow the counterclockwise order convention.
:::

**Edge and vertex neighbors**

Each edge in an MPAS mesh is shared by two cells and connects two vertex points. Each vertex has three edges that meet at the vertex and three cells that share the vertex.

**[Figure 2.4: Neighbors of edges and vertices. Cell indices are given in black, edge indices are red and vertex indices are given in blue. To be added next session.]**

:::{admonition} MPAS code
:class: note

In the example given in the left figure in Figure 2.4, the cells sharing edge 17 are given by `cellsOnEdge(1:2,17) = (30, 42)`. By convention the lower number cell is given first followed by the higher number cell. Similarly, the vertices on the edge are given by `verticesOnEdge(1:2,17) = (11, 13)`. This also follows a convention that will be discussed in section 2.1.4. The right figure in Figure 2.4 shows the neighbor components of vertex 11. The cells sharing vertex 11 are given by `cellsOnVertex(1:3,11) = (30, 58, 42)` and the edges meeting at the vertex are given as `edgesOnVertex(1:3,11) = (17, 44, 122)`.
:::

### 2.1.4 Horizontal velocity on the MPAS horizontal mesh

Two horizontal velocities are defined on the MPAS mesh, the edge-normal velocity, denoted as $u$ in the MPAS code, and the edge-tangential velocity that is denoted as $v$. These velocities are defined following the convention given in Figure 2.5. By convention, a positive edge-normal velocity $u$ at edge $Edge$ points from the cell index `cellsOnEdge(1,Edge)` to cell index `cellsOnEdge(2,Edge)`. Also by convention, `cellsOnEdge(1,Edge)` is the smaller the two cell indices. A positive tangential velocity $v$ points from `verticesOnEdge(1,Edge)` to `verticesOnEdge(2,Edge)`. These velocities follow a right-hand rule - the cross product of the unit vectors at the edge points corresponding to positive $u$ and $v$ point outward from the sphere or Cartesian plane.

**[Figure 2.5: Definition of the prognostic edge-normal horizontal velocity and the diagnostic horizontal edge-tangential velocity using the example mesh given in Figure 2.4. To be added next session.]**

### 2.1.5 Areas on the MPAS horizontal mesh

Three areas are defined on the MPAS mesh and they are depicted in Figure 2.6: The Voronoi mesh cell areas, the areas of the Delauney triangular (dual) mesh, and the areas formed by the intersection of a Voronoi cell and a Delauney triangle, denoted as kite areas. The areas are on the plane for the Cartesian plane configuration of MPAS, and the area on the surface of the sphere for the spherical configuration of MPAS.

**[Figure 2.6: Horizontal mesh areas used in MPAS. To be added next session.]**

:::{admonition} MPAS code
:class: note

The arrays containing the areas for the MPAS mesh are defined in the `Registry.xml` file, `MPAS/src/core_atmosphere/Registry.xml`, and are contained in the arrays `areaCells(nCells)`, `areaTriangle(nVertices)` and `kiteAreasOnVertex(3,nVertices)`. The `kiteAreasOnVertex` follow the same ordering as `cellsOnVertex`, e.g. the kite area `kiteAreasOnVertex(1,Vertex)` will be the kite area associated with cell `cellsOnVertex(1,Vertex)`, and likewise for 2 and 3.
:::

## 2.2 MPAS Vertical Grid

In this section we provide an overview of the computational vertical coordinate used in MPAS and its relation to the physical heights of layers and interfaces in the column, along with a description of the hybrid smoothed 3D mesh that is used in MPAS-Atmosphere.

MPAS-Atmosphere uses a computational vertical coordinate $\zeta$ that represents geometric height in the absence of topography, in contrast to the scaled hydrostatic pressure (mass) vertical coordinate used in WRF (Skamarock et al. 2021a). The Jacobian $\partial\zeta/\partial z$ is used to transform the governing equations cast in terms of $z$ to those using $\zeta$ (see Chapter 3).

The vertical grid in MPAS is structured - the vertical coordinate axis points outward from the center of the earth (or upward from the Cartesian plane), and quantities defined on the vertical axis have neighbors adjacent in the arrays in which they are stored. The MPAS vertical grid defines variables at *layers* and at layer *interfaces*. The interfaces are where the vertical velocities $\Omega$ and $w$, and the height above sea level $z$ (or above the Cartesian plane) are defined. The layers are where all other variables ($u, \rho, \theta, q$, etc) are defined, and by definition the layer heights are halfway between the interface heights. For example, referring to Figure 2.7:

**[Figure 2.7: The MPAS-Atmosphere vertical grid with the locations of the staggered variables. The dashed red lines indicate the center of the layers, and solid lines are interfaces. To be added next session.]**

$$
z_{\mathrm{layer}}(2) = \left[z_{\mathrm{int}}(3) + z_{\mathrm{int}}(2)\right]/2.
$$

### 2.2.1 Specifying the interface heights

The computational vertical coordinate $\zeta$ represents the heights of the layer interfaces and layer midpoints in the absence of topography. As noted earlier, the layer heights are just the average of the bounding interface heights, so $\zeta$ is completely specified by the interface heights. In the following, we describe options to specify the heights $\zeta$ available in the real data configurations; idealized configurations have different (usually hardwired) configurations that will be discussed in Chapter 8.

One option in the MPAS initialization sequence is to specify the interface heights in an ascii file supplied by the user.

:::{admonition} MPAS code
:class: note

User specified levels are enabled by supplying an ascii file with the levels and setting the configuration variable `config_specified_zeta_levels` in the `namelist.init_atmosphere` to the name of the ascii file. If no file is specified (the default configuration) then an analytic specification of the interface level distribution is used. Also see the MPAS Users Guide available at <https://mpas-dev.github.io/atmosphere/atmosphere_download.html>
:::

The MPAS model also provides analytic functions to define the interface layer heights $\zeta$. For a given number of interface levels $n$, where level 1 is the surface, $\zeta(1) = 0$ and $\zeta(n) = z_{\mathrm{top}}$, $\zeta(k)$ is given by

$$
\zeta(k) = z_{\mathrm{top}}\left(\frac{k-1}{n-1}\right)^{1.5} \qquad k = 1\cdots n
$$ (eq:2.1)

:::{admonition} MPAS code
:class: note

Use of {eq}`eq:2.1` is enabled by setting the `namelist.init_atmosphere` variable `config_tc_vertical_grid = false`. The variable $z_\mathrm{top}$ is available through the `namelist.init_atmosphere` variable `config_ztop` and the units are meters.
:::

Another option for distributing the $\zeta$ levels given using the following function form:

For $(k-1)/(n-1) < \zeta_l$

$$
\begin{aligned}
\zeta(k) = z_{\mathrm{top}}\Bigl\{{}&a_s\!\left(\frac{k-1}{n-1}\right) \\
&+ \bigl[3(1-a_t) + 2(a_t-a_s)\zeta_l\bigr]\!\left(\frac{1}{\zeta_l}\frac{k-1}{n-1}\right)^{\!2} \\
&- \bigl[2(1-a_t) + (a_t-a_s)\zeta_l\bigr]\!\left(\frac{1}{\zeta_l}\frac{k-1}{n-1}\right)^{\!3}\Bigr\}.
\end{aligned}
$$ (eq:2.2)

else

$$
\zeta(k) = z_{\mathrm{top}}\left(1 - a_t(1-\zeta_l) + a_t\left[\left(\frac{k-1}{n-1}\right) - \zeta_l\right]\right).
$$ (eq:2.3)

:::{admonition} MPAS code
:class: note

Use of {eq}`eq:2.2` and {eq}`eq:2.3` is enabled by setting the `namelist.init_atmosphere` variable `config_tc_vertical_grid = true`, otherwise {eq}`eq:2.1` is used. The coefficients $(a_s, a_t, \zeta_l) = (0.075, 1.23, 0.31)$ produce a $\zeta$ distribution similar to that given in the default WRF model configuration using 41 interface levels, and these values are used for configurations having less than 55 interface levels. For 55 or greater vertical interface levels the value $a_t = 1.70$ is used and the other coefficients are unchanged. The analytic formula are available for real-data cases and are used in the initialization code `MPAS/src/core_init_atmosphere/mpas_init_atm_cases.F` in subroutine `init_atm_case_gfs`.

We note here that users could replace any of these analytic functions for $\zeta$ with their own function if desired.
:::

### 2.2.2 Hybrid smoothed terrain-following coordinate

Given the 1D computational coordinate $\zeta$, the physical heights $z(\zeta, \vec{x})$ of the interface levels must be set. To accomplish this MPAS-Atmosphere employs a hybrid smoothed terrain-following vertical coordinate described in Klemp (2011, equation 4):

$$
z = \zeta + A\,h_s(\zeta, \vec{x}),
$$ (eq:2.4)

where $z$ is the geometric height, $\zeta$ is the vertical computational coordinate where $0 \le \zeta \le z_{\mathrm{top}}$ for a model top (lid) at height $z_{\mathrm{top}}$, and $\vec{x}$ is a horizontal position vector. $A(\zeta)$ controls the rate at which the coordinate transitions from terrain following at the surface toward constant height surfaces aloft,

$$
\begin{aligned}
A(\zeta) &= \cos^6\!\left(\frac{\pi}{2}\frac{\zeta}{z_H}\right) & \text{for}\quad \zeta &< z_H \\
A(\zeta) &= 0 & \text{for}\quad \zeta &\ge z_H,
\end{aligned}
$$ (eq:2.5)

where $z_H = 30$ km is the default value for MPAS-Atmosphere. $h_s(\zeta, \vec{x})$ represents the terrain influence in the vertical-coordinate, and smoother versions of the terrain are used for increasing $\zeta$. The smoother representations of the terrain are generated using multiple passes of a horizontal Laplacian filter. Denoting the actual terrain elevation as $h$, and the discrete model level $\zeta$ at the surface as model level $n = 1$, we can cast the smoothing operator for increasing level $n$ as

$$
\begin{aligned}
h_s^{(1)} &= h \\
h_s^{(n)} &= \bigl[1 + \beta(\zeta, \vec{x})\,\widetilde{\nabla}^2_\zeta\bigr]^{m}\,h_s^{(n-1)}.
\end{aligned}
$$ (eq:2.6)

The horizontal (on $\zeta$ levels) dimensionless Laplacian is discretized as

$$
\widetilde{\nabla}^2_\zeta\,h_s^{(n-1)} = \sum_{e_c}\frac{d_v}{d_c}\,\delta_e\!\left(h_s^{(n-1)}\right).
$$

In the current MPAS-Atmosphere release, the coefficient $\beta$ multiplying the Laplacian is dependent on its horizontal and vertical position on the computational mesh:

$$
\beta(\zeta, \vec{x}) = \max\!\left[0.01,\ \frac{1}{8}\min\!\left(1.0,\ \frac{3000}{\overline{dc(\vec{x})}}\right)\right] \times \min\!\left[\left(\frac{3\zeta}{z_H}\right)^{2},\ 1.0\right].
$$ (eq:2.7)

This formulation includes a function of the horizontal cell spacing $dc(\vec{x})$ that decreases the smoothing coefficient as the cell spacing increases that is not present in the original Klemp (2011) specification. This factor provides some compensation for the increased smoothing that is inherent in interpolating the actual terrain to a coarser horizontal mesh, and provides more consistency in applying the smoothed coordinate over a variable resolution MPAS mesh. Finally, the number of smoothing passes at each vertical coordinate level is defined as $m(n) = N_s + n$ with the default value of $N_s$ set to 30. This also differs from Klemp (2011) where $m$ was specified as a constant. Since the smaller scale terrain features are removed more efficiently by the Laplacian filter, increasing $m$ with height allows greater smoothing of larger larger scale terrain influences that become more dominate at higher levels.

At each application $i = 1, m$ of the smoother {eq}`eq:2.6` the resulting heights are bounded so as to limit the minimum spacing between levels in a given column. Given a smoothed preliminary value of $h_s(k)$ at level $k$, the value is checked to see if the level spacing is above a minimum value:

$$
\left[z(k) - z(k-1)\right]_{\mathrm{min}} < \gamma\,\left[\zeta(k) - \zeta(k-1)\right].
$$ (eq:2.8)

where the default value for $\gamma = 0.3$. On a given iteration $i$ of the smoother, if {eq}`eq:2.8` is not satisfied the increment to $z$ is not applied and the iteration is halted.

Figure 2.8 shows an example of the coordinate surfaces from a traditional terrain following coordinate where $A(\zeta)$ is a linear function with a value of 1 at the surface and 0 at the model top and $h_s = h$, and the smoothed hybrid terrain following coordinate.

**[Figure 2.8: The terrain following coordinate used in MPAS. The left figure shows the traditional terrain following coordinate surfaces and the right figure shows the smoothed vertical coordinate. To be added next session.]**

:::{admonition} MPAS code
:class: note

For earth applications, the generation of the full 3D mesh and its coordinate surfaces occurs in the initialization subroutine `init_atm_case_gfs` in `MPAS/src/core_init_atmosphere/mpas_init_atm_cases.F`. The height where the coordinate surface transitions from hybrid terrain following to constant height at $z_H$ is set in the initialization code to 30 km.
:::

### 2.2.3 Vertical Interpolation

There is a need to interpolate layer variables to interfaces, and to interpolate interface variables to layers. For the latter, we take the average of the interface variables values from the two surrounding interface to set the layer value given that a layer lies halfway between the interfaces. We have two methods to interpolate layer values to an interface. The first is to do a linear in $\zeta$ interpolation of the layer value to the interface. Taking into account that the layers lie halfway between interfaces, we can write this interpolation as

$$
\begin{aligned}
\phi_{\mathrm{int}}(k) = {}&\frac{\Delta\zeta_w(k-1)}{\Delta\zeta_w(k) + \Delta\zeta_w(k-1)}\,\phi_{\mathrm{layer}}(k) \\
&+ \frac{\Delta\zeta_w(k)}{\Delta\zeta_w(k) + \Delta\zeta_w(k-1)}\,\phi_{\mathrm{layer}}(k-1).
\end{aligned}
$$ (eq:2.9)

The second approach is to vertically integrate the quantity between the layers containing the interface and divide by the height difference between the two layers. We assume the quantity is constant between the interfaces, and in this case the result simply switches the weights in {eq}`eq:2.9` between the two layer values, i.e.

$$
\begin{aligned}
\phi_{\mathrm{int}}(k) = {}&\frac{\Delta\zeta_w(k)}{\Delta\zeta_w(k) + \Delta\zeta_w(k-1)}\,\phi_{\mathrm{layer}}(k) \\
&+ \frac{\Delta\zeta_w(k-1)}{\Delta\zeta_w(k) + \Delta\zeta_w(k-1)}\,\phi_{\mathrm{layer}}(k-1).
\end{aligned}
$$ (eq:2.10)

:::{admonition} MPAS code
:class: note

The weights for the vertical interpolation of values from a layer to an interface are computed and stored during model initialization, where a user can choose between a linear interpolation {eq}`eq:2.9` or the integrated average value {eq}`eq:2.10`. In MPAS the weights are stored in the arrays `fzm` and `fzp`, where the weight `fzm` multiplies $\phi_{\mathrm{layer}}(k)$ in {eq}`eq:2.9` and {eq}`eq:2.10`, and likewise `fzp` multiplies $\phi_{\mathrm{layer}}(k-1)$. In some places in the MPAS code the array `fnm` points to `fzm` and `fnp` points to `fzp`.
:::
