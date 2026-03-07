# Appendix C: Grid Description

This chapter provides a brief introduction to the common types of grids used in the MPAS framework.

## C.1 Horizontal grid

The MPAS grid system requires the definition of seven elements. These seven elements are composed of two types of *cells*, two types of *lines*, and three types of *points*. These elements are depicted in Figure C.1 and defined in Table C.1. These elements can be defined on either the plane or the surface of the sphere. The two types of cells form two meshes, a primal mesh composed of Voronoi regions and a dual mesh composed of Delaunay triangles. Each corner of a primal mesh cell is uniquely associated with the "center" of a dual mesh cell and vice versa. So we define the two meshes as either a primal mesh (composed of cells P_i) or a dual mesh (composed of cells D_v). The center of any primal mesh cell, P_i, is denoted by **x**_i and the center of any dual mesh cell, D_v, is denoted by **x**_v. The boundary of a given primal mesh cell P_i is composed of the set of lines that connect the **x**_v locations of associated dual mesh cells D_v. Similarly, the boundary of a given dual mesh cell D_v is composed of the set of lines that connect the **x**_i locations of the associated primal mesh cells P_i. As shown in Figure C.1, a line segment that connects two primal mesh cell centers is uniquely associated with a line segment that connects two dual mesh cell centers. We assume that these two line segments cross and the point of intersection is labeled as **x**_e. In addition, we assume that these two line segments are orthogonal as indicated in Figure C.1. Each **x**_e is associated with two distances: d_e measures the distance between the primal mesh cells sharing **x**_e and l_e measures the distance between the dual mesh cells sharing **x**_e.

Since the two line segments crossing at **x**_e are orthogonal, these line segments form a convenient local coordinate system for each edge. At each **x**_e location a unit vector **n**_e is defined to be parallel to the line connecting primal mesh cells. A second unit vector **t**_e is defined such that **t**_e = **k** x **n**_e.

Table C.2 provides the names of all *elements* and all *sets of elements* as used in the MPAS framework. Elements appear twice in the table when described in the grid file in more than one way, e.g. points are described with both cartesian and latitude/longitude coordinates. An `ncdump -h` of any MPAS grid, output or restart file will contain all variable names shown in the second column of Table C.2.

In addition to these seven element types, we require the definition of *sets of elements*. In all, eight different types of sets are required and these are defined and explained in Table C.3 and Figure C.2. The notation is always of the form of, for example, i in CE(e), where the LHS indicates the type of element to be gathered (cells) based on the RHS relation to another type of element (edges).

The angle of each edge in an MPAS grid is provided in the variable `angleEdge`. The angle given is the angle between a vector pointing north and a vector pointing in the positive tangential direction of the edge. Referring to Fig. C.3,

```
angleEdge = arcsin ||n_hat x v_hat||
```

where **n** is the unit vector pointing north and **v** is the unit vector pointing from `verticesOnEdge(1,iEdge)` to `verticesOnEdge(2,iEdge)`.

Given a wind vector (u_perp, u_par) defined in terms of components orthogonal to and parallel to the edge, the earth-relative wind (u, v) may be recovered as

```
| u |   | cos(a)  -sin(a) | | u_perp |
|   | = |                   | |        |
| v |   | sin(a)   cos(a) | | u_par  |
```

where a = `angleEdge`.

### Table C.1: Definition of elements used to build the MPAS grid

| Element | Type | Definition |
|---------|------|------------|
| **x**_i | point | location of center of primal-mesh cells |
| **x**_v | point | location of center of dual-mesh cells |
| **x**_e | point | location of edge points where velocity is defined |
| d_e | line segment | distance between neighboring **x**_i locations |
| l_e | line segment | distance between neighboring **x**_v locations |
| P_i | cell | a cell on the primal-mesh |
| D_v | cell | a cell on the dual-mesh |

### Table C.2: Variable names used to describe an MPAS grid

| Element | Name | Size | Comment |
|---------|------|------|---------|
| **x**_i | `{x,y,z}Cell` | nCells | cartesian location of **x**_i |
| **x**_i | `{lon,lat}Cell` | nCells | longitude and latitude of **x**_i |
| **x**_v | `{x,y,z}Vertex` | nVertices | cartesian location of **x**_v |
| **x**_v | `{lon,lat}Vertex` | nVertices | longitude and latitude of **x**_v |
| **x**_e | `{x,y,z}Edge` | nEdges | cartesian location of **x**_e |
| **x**_e | `{lon,lat}Edge` | nEdges | longitude and latitude of **x**_e |
| d_e | `dcEdge` | nEdges | distance between **x**_i locations |
| l_e | `dvEdge` | nEdges | distance between **x**_v locations |
| e in EC(i) | `edgesOnCell` | (nEdgesMax, nCells) | edges that define P_i |
| e in EV(v) | `edgesOnVertex` | (3, nCells) | edges that define D_v |
| i in CE(e) | `cellsOnEdge` | (2, nEdges) | primal-mesh cells that share edge e |
| i in CV(v) | `cellsOnVertex` | (3, nVertices) | primal-mesh cells that define D_v |
| v in VE(e) | `verticesOnEdge` | (2, nEdges) | dual-mesh cells that share edge e |
| v in VI(i) | `verticesOnCell` | (nEdgesMax, nCells) | vertices that define P_i |

### Table C.3: Definition of element groups used to reference connections in the MPAS grid

| Syntax | Output |
|--------|--------|
| e in EC(i) | set of edges that define the boundary of P_i |
| e in EV(v) | set of edges that define the boundary of D_v |
| i in CE(e) | two primal-mesh cells that share edge e |
| i in CV(v) | set of primal-mesh cells that form the vertices of dual mesh cell D_v |
| v in VE(e) | the two dual-mesh cells that share edge e |
| v in VI(i) | the set of dual-mesh cells that form the vertices of primal-mesh cell P_i |
| e in ECP(e) | edges of cell pair meeting at edge e |
| e in EVC(v, i) | edge pair associated with vertex v and mesh cell i |

**Example element group references** (see Figure C.2):

```
e in EC(P1)  = [e1, e2, e3, e4, e5, e6]
e in EV(D1)  = [e1, e6, e7]
i in CE(e1)  = [P1, P2]
i in CV(D1)  = [P1, P2, P3]
v in VE(e1)  = [D1, D2]
v in VC(P1)  = [D1, D2, D3, D4, D4, D5, D6]
e in ECP(e1) = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11]
e in ECV(P1, D1) = [e1, e6]
```

## C.2 Vertical grid

The vertical coordinate in MPAS-Atmosphere is zeta and has units of length, where 0 <= zeta <= z_t and z_t is the height of the model top. The relationship between the vertical coordinate and height in the physical domain is given as

```
z = zeta + A * h_s(x, y, zeta)                    (C.1)
```

where (x, y) denotes a location on the horizontal mesh and zeta is the vertical coordinate (zeta is directed radially outward from the surface of the sphere, or perpendicular to the horizontal (x, y) plane in a Cartesian coordinate MPAS-A configuration). MPAS-A can be configured with the traditional Gal-Chen and Somerville terrain-following coordinate by setting h_s(x, y, zeta) = h(x, y) and A = 1 - zeta/z_t, where h(x, y) is the terrain height. Alternatively, A can be modified to allow a more rapid or less rapid transition to the constant-height upper boundary condition. Additionally, a constant-height coordinate can be specified at some intermediate height below z_t.

The influence of the terrain on any coordinate surface zeta can be influenced by the specification of h_s(x, y, zeta). Specifically, h_s can be set such that h_s(x, y, 0) = h(x, y) (i.e. terrain following at the surface), and progressively filtered fields of h(x, y) can be used at zeta > 0 in h_s(x, y, zeta), such that the small-scale features in the topography are quickly filtered from the coordinate.

On the MPAS-A mesh C-grid staggering, the state variables u, rho, theta and scalars are located halfway between w levels in both physical height and in the coordinate zeta. Variables associated with the coordinate systems used in the MPAS-A solver, and possibly appearing in its input, output or history files, are defined in Table C.4.

### Table C.4: Vertical coordinate variables in MPAS-Atmosphere

*level* is the integer model level (usually specified with index k where k = 1 is the lowest model level and physical height increases with increasing k). Delta denotes a vertical difference between levels, and *cell* is a given mesh cell on the primary mesh.

| Variable | Definition |
|----------|------------|
| `zgrid(level, cell)` | physical height of the w points in meters |
| `zw(level)` | zeta at w levels |
| `zu(level)` | zeta at u levels; zu(k) = [zw(k+1) + zw(k)] / 2 |
| `dzw(level)` | Delta-zeta at u levels; dzw(k) = zw(k+1) - zw(k) |
| `dzu(level)` | Delta-zeta at w levels; dzu(k) = [dzw(k+1) + dzw(k)] / 2 |
| `rdzw(level)` | 1 / dzw |
| `rdzu(level)` | 1 / dzu |
| `zz(level, cell)` | Delta-zeta / Delta-z at u levels; (zw(k+1) - zw(k)) / (zgrid(k+1, cell) - zgrid(k, cell)) |
| `fzm(level)` | weight for linear interpolation to w(k) point for u(k) level variable |
| `fzp(level)` | weight for linear interpolation to w(k) point for u(k-1) level variable |

### Reference

Klemp, J. B. (2011). A Terrain-Following Coordinate with Smoothed Coordinate Surfaces. *Mon. Wea. Rev.*, **139**, 2163-2169. doi:10.1175/MWR-D-10-05046.1
