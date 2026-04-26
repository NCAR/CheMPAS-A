# Chapter 4: Spatial Discretization

In Chapter 3 we examined the time integration of the MPAS-Atmosphere temporally discrete prognostic equations {eq}`eq:3.15`–{eq}`eq:3.19` which include the shallow atmosphere approximation. We also explored the vertical discretization of the acoustic-step equations for $\Omega$, $\Theta_m$ and $\tilde{\rho}$. In this chapter we describe the spatial discretization of the right-hand-side terms in {eq}`eq:3.15`–{eq}`eq:3.19`.

## 4.1 Mass Conservation Equation

The continuous equation describing conservation of dry air mass can be expressed as

$$
\frac{\partial \tilde{\rho}_d}{\partial t} = -\nabla\cdot\mathbf{V}
$$ (eq:4.1)

where $\mathbf{V} = (\mathbf{V}_{\mathbf{H}}, \Omega) = (\tilde{\rho}_d\mathbf{v}_{\mathbf{H}}, \tilde{\rho}_d\omega)$. The finite volume formulation employed in MPAS integrates {eq}`eq:4.1` over a cell volume,

$$
\int_D\!\left[\frac{\partial \tilde{\rho}_d}{\partial t} = -\nabla\cdot\mathbf{V}\right]dV,
$$

and applies the divergence theorem to turn the volume integral into an area integral over the cell sides:

$$
\frac{\partial \overline{(\tilde{\rho}_d)}}{\partial t} = -\,\frac{1}{V}\oint_\Sigma \mathbf{V}\cdot\mathbf{n}\,d\sigma,
$$ (eq:4.2)

where the left hand-side-term is the time rate of change of the cell averaged value of $\tilde{\rho}_d$ and the right-hand-side is an integral of the integral form of the mass-flux divergence, i.e. the mass flux through the cell boundaries. $\mathbf{n}$ is the outward-pointing unit vector normal to a cell edge. Note that {eq}`eq:4.2` is still exact.

The left-hand-side term in {eq}`eq:4.2` is discretized by approximating the mass flux as the instantaneous value of the momentum normal to the cell faces (these are the MPAS prognostic variables $\mathbf{V}_{\mathbf{H}}$ and $\Omega$) multiplied by the area of the cell face divided by the cell volume. Thus the discrete form of the flux divergence for cell $i, k$ denoting the $i$th horizontal cell at the $k$th level:

$$
-\,\frac{1}{V_{k,i}}\!\left[\sum_{n_{e_i}}\bigl(A_{e_i}\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}\bigr) + A_c\bigl(\Omega_{k+1,i} - \Omega_{k,i}\bigr)\right],
$$ (eq:4.3)

where $V_{i,k}$ is the volume of the cell, and $A_{e_i}$ and $\mathbf{V}_{\mathbf{H},\mathbf{e_i},\mathbf{k}}$ are the area of the hexagon face and normal mass flux through cell face $e_i$ of cell $i$ at level $k$, respectively. For the shallow-atmosphere approximation, {eq}`eq:4.3` can be written as

$$
-\,\frac{1}{A_i}\sum_{n_{e_i}}\bigl(L_{e_i}\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}\bigr) + \frac{\Omega_{k+1,i} - \Omega_{k,i}}{\Delta\zeta_w(k)},
$$ (eq:4.4)

where $L_{e_i}$ is the length of the horizontal cell edge $e_i$, and $\Delta\zeta_w$ is cell thickness at level $k$. This is the form used in the MPAS-A code.

:::{admonition} MPAS code
:class: note

The mass flux divergence {eq}`eq:4.4` appears in a number of locations in MPAS-A. For example, it is used to compute $R^t_{\tilde{\rho}_d}$ in {eq}`eq:3.33` (the array `rho_tend` in the code), and this occurs early in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The quantity $1/\Delta\zeta_w$ is stored in the array `rdzw` in the MPAS code. {eq}`eq:4.3` and {eq}`eq:4.4` are identical for the shallow atmosphere approximated equations on the discrete MPAS mesh.
:::

We repeat here that the momentum $\mathbf{V} = (\mathbf{V}_{\mathbf{H}}, \Omega) = (\overline{\tilde{\rho}_d}^{c_e}\mathbf{v}_{\mathbf{H}}, \overline{\tilde{\rho}_d}^{k}\omega)$, where here $\overline{\tilde{\rho}_d}^{c_e}$ denotes the horizontal average of $\tilde{\rho}_d$ to the cell edge and $\overline{\tilde{\rho}_d}^{k}$ is interpolation of $\tilde{\rho}_d$ to the interface level $k$. This represents a simple 2nd order projection of $\tilde{\rho}_d$ to the cell faces, and we use it because the density is a very smoothly varying field while other scalar mass fields are often not smooth. In the following section we describe transport for cell-centered scalars other than the dry-air mass, where a higher order projection is described.

## 4.2 Transport for Cell-Centered Variables

For the flux divergence evaluations in the transport of scalars in {eq}`eq:3.10`, in the potential temperature equation {eq}`eq:3.8`, in the equation for $W$ in {eq}`eq:3.7`, and for the vertical flux divergence needed to horizontal momentum in {eq}`eq:3.6`, we use a higher-order projection of the relevant variable to the cell edges to produce a more-accurate approximation to the flux divergence.

The transport equation for a scalar $\phi_i$ in cell $i$,

$$
\frac{\partial(\tilde{\rho}_d\phi)}{\partial t} = -\nabla\cdot\mathbf{V}\phi,
$$ (eq:4.5)

is discretized using the same finite-volume approach used for the dry-air conservation equation {eq}`eq:4.1`. The discrete flux divergence becomes

$$
-\nabla\cdot\mathbf{V}\phi \approx -\,\frac{1}{A_i}\sum_{n_{e_i}}\left[L_{e_i}\bigl(\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}\,\phi_{e_i}\bigr)\right] + \frac{\Omega_{k+1,i}\overline{\phi}^{\,k+1,i} - \Omega_{k,i}\overline{\phi}^{\,k,i}}{\Delta\zeta_w(k)}.
$$ (eq:4.6)

In the following subsections we describe the horizontal flux divergence operators and the vertical flux divergence operator, i.e. we describe how the cell-edge projections of $\phi_{e_i}$ and $\overline{\phi}^{\,k}$ are computed.

### 4.2.1 Horizontal Flux Divergences on the Voronoi Mesh

The evaluation of horizontal flux divergence on the Voronoi mesh and its performance in the RK3 scheme for transport tests are described in Skamarock and Gassmann (2011). This scheme generalizes the flux formula used in Wicker and Skamarock (2002) and first presented in Hundsdorfer et al. (1995). The mass flux per unit area through cell face $e_i$ is given by

$$
\begin{aligned}
\bigl(\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}\,\phi_{e_i}\bigr) = \bigl(\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}\bigr)\biggl[{}&\tfrac{1}{2}(\phi_{c_0} + \phi_{c_i}) \\
&- \tfrac{\Delta x^2_{e_i}}{12}\!\left\{\left(\tfrac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_i} + \left(\tfrac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_0}\right\} \\
&+ \mathrm{sign}(\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n})\,\tfrac{\beta\,\Delta x^2_{e_i}}{12}\!\left\{\left(\tfrac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_i} \right. \\
&\qquad\left. - \left(\tfrac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_0}\right\}\biggr],
\end{aligned}
$$ (eq:4.7)

where a positive flux velocity flows from cell $c_0$ to cell $c_i$ through edge $e_i$, $\Delta x_{e_i}$ is the distance between the cell centers of cells $c_0$ and $c_i$, and the derivative $\partial x_{e_i}$ is the horizontal derivative in the $\overrightarrow{c_0c_i}$ direction. In {eq}`eq:4.7`, the values at the cell centers $\phi_{c_0}$ and $\phi_{c_i}$ are known, but their second derivatives at the cell centers in the direction $\overrightarrow{c_0c_i}$ are not known. $\beta$ is a coefficient scaling the magnitude of the upwinding in the scheme; $\beta = 0$ results in a 4th-order scheme and $\beta = 1$ is the third-order scheme from Hundsdorfer et al. (1995). As described in Skamarock and Gassmann (2011), we compute a least-squares fit polynomial for $\phi$ using values at the cell center and at the neighboring cell for each of the two cells sharing the cell edge $e_i$, and we compute the second derivative of the polynomials in the direction $\overrightarrow{c_0c_i}$ and use this in {eq}`eq:4.7`. The least-squares polynomial computation is described in Appendix B.

**Horizontal Transport Coefficients**

As noted in Appendix B, the least-squares-polynomial computation (B.2) requires a matrix inversion, but given that the matrix depends only on the geometry of the Voronoi mesh it is computed and stored as part of the initialization process by the MPAS initialization core *init_atmosphere*. The cells that contribute to the projection of a mixing ratio for flux on a cell edge requires values from the two cells sharing the edge and all of the neighbors of these two cells (see Figure B.1 in Appendix B). To streamline the transport calculation in the MPAS integration, we assemble all the coefficients needed to compute the flux {eq}`eq:4.7` as part of the start-up of the MPAS model integration such that the flux computation is a vector-vector multiply given the mixing ratios and previously assembled weights.

:::{admonition} MPAS code
:class: note

The coefficients used in the transport routines are precomputed (and stored) in subroutine `atm_adv_coef_compression` found in `MPAS/src/core_atmosphere/mpas_atm_core.F`.
For a given edge, a list of cells needed in the flux computation is constructed and stored. For example, the list would consist of cells $c_0$ through $c_9$ for edge $e_1$ in Figure B.1. Two set of coefficients are computed and stored. The set corresponding to the 4th-order flux for edge $e_1$

$$
L_{e_1}\!\left[\frac{1}{2}(\phi_{c_0} + \phi_{c_i}) - \Delta x^2_{e_i}\,\frac{1}{12}\!\left\{\left(\frac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_1} + \left(\frac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_0}\right\}\right]
$$

are stored in the array `adv_coef(source cells, edge)`. The set corresponding to the 3rd-order upwind piece,

$$
L_{e_1}\Delta x^2_{e_i}\,\frac{\beta}{12}\!\left\{\left(\frac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_1} - \left(\frac{\partial^2\phi}{\partial x^2_{e_i}}\right)\!\bigg|_{c_0}\right\},
$$ (eq:4.8)

are stored in the `adv_coef_3rd(source cells, edge)`. The source cells for each edge flux computation are stored in the array `advCellsForEdge(cell,iEdge)`. In the MPAS code, the coefficient $\beta$ in {eq}`eq:4.8` is included in a subsequent call to subroutine `atm_couple_coef_3rd_order`.
:::

### 4.2.2 Vertical Fluxes

The vertical transport in {eq}`eq:4.6` requires the computation of the vertical fluxes $\Omega\phi$ located at the $(w, \Omega)$ points of the MPAS mesh. We use a similar operator as used for the horizontal fluxes except that we make use of the structured vertical coordinate. The vertical flux at level $k$ in {eq}`eq:4.6` is

$$
\begin{aligned}
\Omega_{k,i}\,\overline{\phi}^{\,k,i} = \Omega_{k,i}\!\biggl\{{}&\tfrac{1}{2}(\phi_{k,i} + \phi_{k-1,i}) - \tfrac{1}{12}\bigl(\delta^2_\zeta\phi_{k,i} + \delta^2_\zeta\phi_{k-1,i}\bigr) \\
&+ \mathrm{sign}(\Omega)\,\tfrac{\beta}{12}\bigl(\delta^2_\zeta\phi_{k,i} - \delta^2_\zeta\phi_{k-1,i}\bigr)\biggr\},
\end{aligned}
$$

where $\delta^2_\zeta\phi_{k,i} = \phi_{k+1,i} - 2\phi_{k,i} + \phi_{k-1,i}$. In the MPAS code the flux is equivalently written as

$$
\begin{aligned}
(\Omega\phi)_{k,i} = \Omega_{k,i}\!\biggl\{{}&\tfrac{1}{12}\bigl[7(\phi_{k,i} + \phi_{k-1,i}) - (\phi_{k+1,i} + \phi_{k-2,i})\bigr] \\
&+ \mathrm{sign}(\Omega_{k,i})\,\tfrac{\beta}{12}\bigl[(\phi_{k+1,i} - \phi_{k-2,i}) - 3(\phi_{k,i} - \phi_{k-1,i})\bigr]\biggr\}.
\end{aligned}
$$ (eq:4.9)

As with the horizontal fluxes, the parameter $\beta$ controls the level of upwinding, and hence damping, from this component of the transport. This form of the flux divergence originally appeared in Hundsdorfer et al. (1995) and was first implemented in the RK3 transport solver in Wicker and Skamarock (2002). It is also used in the WRF model (Skamarock et al. 2021b). In the referenced descriptions, the coefficient $\beta = 1$. For $\beta = 0$ the flux divergence is 4th-order accurate and neutral in the integration. $\beta > 0$ introduces damping in the flux divergence. In MPAS $\beta = 0.25$ is the default; see (Skamarock and Gassmann 2011) for details.

:::{admonition} MPAS code
:class: note

The vertical flux divergence calculations for scalars occur in the scalar transport routines `atm_advance_scalars` and `atm_advance_scalars_mono` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. Given the horizontal and vertical flux formulae described in sections 4.2.1 and 4.2.2, the updates in the first two substeps of the RK3 integration, {eq}`eq:3.46` and {eq}`eq:3.47`, can be computed, and the fluxes are also needed for the final monotonic update. The vertical flux divergence for potential temperature (the remaining cell-centered prognostic variable) is computed in subroutine `atm_compute_dyn_tend`. The quantity $1/\Delta\zeta_w$ in the flux divergence {eq}`eq:4.6` is stored in the array `rdnw` which is computed in the initialization routines for various cases in `MPAS/src/core_init_atmosphere/mpas_init_atm_cases.F`.
:::

### 4.2.3 Boundary Conditions for the Vertical Fluxes

The vertical flux divergence employs a no-flux condition ($\Omega = 0$) at the surface and at the top of the domain. Referring to the vertical mesh in Figure 2.7, the zero flux condition at the lower boundary is at the $k = 1$ level for $\Omega$. To evaluate the flux divergence for all the layer variables $(u, \theta_m, q)$, vertical fluxes are needed at the interfaces except for the lower and upper boundaries. To compute a vertical flux at interface $k$ using {eq}`eq:4.9`, scalar values at layers $k + 1$, $k$, $k − 1$, and $k − 2$ are needed, i.e. from two adjacent layers on either side of the interface where the flux is computed. For the case of the interfaces closest to the lower boundary (at $k = 2$) and interface closest to the upper boundary (the MPAS-A model top), there are not sufficient layer values to fill out {eq}`eq:4.9`. In this case the value of the fluxed quantity is interpolated to the interface from the two layers. The interpolation is either linear in the computation coordinate $\zeta$ or it is the average value of the integral of the quantity between the layers as described in section 2.2.3.

### 4.2.4 Positive-Definite and Monotonic Flux Renormalization

MPAS uses a shape-preserving (monotonic) flux renormalization adapted from Zalesak (1979). In the time integration scheme {eq}`eq:3.46`–{eq}`eq:3.49`, the flux renormalization is applied on the final update {eq}`eq:3.49`. The update, including the flux renormalization, is accomplished as follows. Beginning with the intermediate solution that has been updated with the physics tendencies, {eq}`eq:3.48`, $q_j^{***} = Q_j^{***}/\tilde{\rho}'_d$:

(1) All the fluxes needed to update the cell values, i.e. the fluxes for both the horizontal faces and vertical faces, are computed and saved.

(2) For each cell, possible cell maximum and minimum values are determined by finding the maximum and minimum values of the scalar $\phi_j^{***}$ from the cell center and all neighbor cells (cells that share an edge).

(3) First-order upwind fluxes are used to update the solution. The upwind update is

$$
\begin{aligned}
\tilde{\rho}^{t+\Delta t}_d\phi^{t*}_j = \tilde{\rho}^t_d\phi^{***}_j - \Delta t\!\biggl[{}&\frac{1}{A_i}\sum_{n_{e_i}}L_{e_i}\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}\,\phi^{***}_{\mathrm{upwind}} \\
&+ \frac{1}{\Delta\zeta_k}\bigl(\Omega_{k+1,i}\phi^{***}_{\mathrm{upwind}} - \Omega_{k,i}\phi^{***}_{\mathrm{upwind}}\bigr)\biggr]
\end{aligned}
$$ (eq:4.10)

(4) Perturbation fluxes are formed by subtracting the upwind fluxes from the full fluxes computed in step (1). For the upwind fluxes, $\phi_{\mathrm{upwind}}$ is the value from the cell upwind of the velocity on the cell face.

(5) A test update is computed using the outgoing perturbation fluxes for each cell. If the value produced is less than the minimum cell value obtained in step (2) a scalar factor (`scale_out(cell)`) is computed that multiplies the fluxes such that the minimum value is obtained. The `scale_out` factor initially set to one before the scaling check, and it is bounded - $0 \le$ `scale_out` $\le 1$.

(6) A test update is computed using the incoming perturbation fluxes for each cell. If the value produced is greater than the maximum cell value obtained in step (2), a scalar factor (`scale_in(cell)`) is computed that multiplies the fluxes such that the maximum value is obtained.

(7) The perturbation fluxes on an edge are scaled (renormalized) with the minimum value of the `scale_out` from the upwind cell and `scale_in` from the downwind cell.

(8) To complete the time step, the values are updated, starting from upwind-updated values $\phi^{t*}_j$ from {eq}`eq:4.10`, with the scaled perturbation fluxes.

:::{admonition} MPAS code
:class: note

The shape-preserving (monotonic) transport scheme is contained in subroutine `atm_advance_scalars_mono` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. This routine is called to execute updates {eq}`eq:3.48` and {eq}`eq:3.49` in the described in section 3.5. The first action in this subroutine is to compute the update {eq}`eq:3.48` by adding the physics tendencies to the scalars to produce the intermediate values $q_j^{***}$. The code then follows steps 1 through 8 which is the update {eq}`eq:3.49` with the flux renormalization. The main loops in the code often include calculations from multiple steps. In contrast to the scalar updates in the first two substeps, {eq}`eq:3.46` and {eq}`eq:3.47`, here each scalar is updated individually given the need to store the fluxes and the scaling factors.
:::

## 4.3 Horizontal Momentum Equation

The horizontal momentum equation {eq}`eq:3.6` in MPAS is cast in vector-invariant form. In the derivation of this vector invariant form of the equation, beginning with the advective form of the horizontal momentum equation [using $u_t$ as opposed to $(\rho u)_t$], the horizontal transport of the horizontal velocity, $(-\mathbf{v}_{\mathbf{H}}\cdot\nabla_H\mathbf{v}_{\mathbf{H}})$, is recast as $-\eta\,\mathbf{k}\times\mathbf{v}_{\mathbf{H}} - \nabla_\zeta K$. Coupling the horizontal velocity with the dry-air density brings in the additional term $\mathbf{v}_{\mathbf{H}}\nabla_\zeta\cdot\mathbf{V}$. Also note that the vector-invariant form does not include the vertical transport term which still appears in {eq}`eq:3.6`.

MPAS uses a C-grid staggering of the horizontal velocity, thus the prognostic velocity is defined in the cell edge and is normal to the cell edge as illustrated in Figure 2.1. In the following sections we review the spatial discretization of the terms in the various forms of the horizontal momentum equation used in MPAS-A.

### 4.3.1 Pressure and Kinetic Energy

There are two places where the pressure gradient term in the horizontal momentum equation needs to be evaluated, first in the RK3 integration using {eq}`eq:3.15`, and also in the acoustic integration using {eq}`eq:3.25`. The kinetic energy gradient is only needed for the RK3 in integration {eq}`eq:3.15`. Consider the velocity $u_{13}$ on edge $ab$ given in Figure 2.1. The line connecting cell centers $A$ and $C$ is perpendicular to the cell edge and is bisected by the cell edge, thus the difference of the cell center values divided by the distance between the cell centers will give a second-order accurate representation of the gradient at the $u_{13}$ edge point. With this gradient formula, and with pressure $\Theta_m$, $\rho$, and the diagnostic quantities $K$ (horizontal kinetic energy), pressure, and exner function defined at the cell centers, we can write the discrete terms as

$$
\begin{aligned}
\text{spatially continuous} &\qquad \text{spatially discrete} \\[4pt]
\frac{\rho_d}{\rho_m}\!\left(\nabla_\zeta\!\left(\frac{p'}{\zeta_z}\right) + g\,z_H\tilde{\rho}'_m\right) &\;\to\; \left.\frac{\rho_d}{\rho_m}\right|_{c_e}\!\left[\delta_{c_e}\!\left(\frac{p'}{\delta_z(\zeta)}\right) + g\,(\delta_{c_e}z)\,\overline{\tilde{\rho}'_m}^{c_e}\right],
\end{aligned}
$$ (eq:4.11)

$$
\tilde{\rho}_d\nabla_\zeta K \;\to\; \overline{\tilde{\rho}_d}^{c_e}\,\delta_{c_e}K,
$$ (eq:4.12)

$$
\begin{aligned}
\frac{\rho_d^t}{\rho_m^t}\bigl[\gamma R_d\pi^t\nabla_\zeta\overline{\Theta''_m}^\tau + g\,z_H\tilde{\rho}''_d\bigr] \;\to\; {}&\left.\frac{\rho_d}{\rho_m}\right|_{c_e}\!\left[\gamma R_d\overline{\pi^t}^{c_e}\delta_{c_e}\!\left(\overline{\Theta''_m}^\tau\right)\right. \\
&\left.\hphantom{\left.\tfrac{\rho_d}{\rho_m}\right|_{c_e}\!\bigl[\,} + g\,(\delta_{c_e}z)\,\overline{\tilde{\rho}_m''}^{c_e}\right].
\end{aligned}
$$ (eq:4.13)

{eq}`eq:4.11` and {eq}`eq:4.12` are found in the RK3 integration {eq}`eq:3.15`, and {eq}`eq:4.13` is found in the acoustic integration {eq}`eq:3.25`. The gradient operator $\delta_{c_e}$ is defined as the difference across the edge $e$ of the cell-centered values from the two cells sharing the edge. Referring to Figure 2.1 and velocity $u_{13}$, where a positive value of $u_{13}$ indicates flow from cell $C$ to cell $A$, the gradient operator is

$$
\delta_{c_e}(\phi) = \frac{\phi_A - \phi_C}{|\overrightarrow{CA}|},
$$ (eq:4.14)

where $|\overrightarrow{CA}|$ is the distance between the cell centers (the great-circle arc distance on the sphere). The averaging operator $\overline{\phi}^{c_e}$ as the arithmetic average of the cell centers sharing the edge $e$:

$$
\overline{\phi}^{c_e} = \frac{\phi_A + \phi_C}{2}.
$$ (eq:4.15)

The density ratio is computed at the edge in MPAS as follows:

$$
\left.\left(\frac{\rho_d}{\rho_m}\right)\right|_{c_e} = \left.\left(\frac{1}{1 + \sum q_j}\right)\right|_{c_e} \;\to\; \overline{\left(\frac{1}{1 + \sum q_j}\right)^{c_e}}.
$$ (eq:4.16)

:::{admonition} MPAS code
:class: note

The horizontal pressure gradient and kinetic energy operators {eq}`eq:4.11` and {eq}`eq:4.12` for the RK3 timestep can be found in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The acoustic step horizontal pressure gradient operator {eq}`eq:4.13` can be found in subroutine `atm_advance_acoustic_step` in the same source file. The coefficients $\rho_d/\rho_m$ for the cell edges {eq}`eq:4.16` are computed in subroutine `atm_compute_moist_coefficients` and the array storing the coefficients is `cqu(levels, edge)`. The cell order in the gradient operator is given by the integer array `cellsOnEdge(2, edges)`, where the positive flow direction across an edge is always from `cellsOnEdge(1, edge)` to `cellsOnEdge(2, edge)` by convention.
:::

**Cell-Centered Kinetic Energy Evaluation**

The discrete kinetic gradient in {eq}`eq:4.12` uses a kinetic energies defined at the cell centers. Ringler et al. (2010) showed that for the shallow water equations on the sphere using an SCVT, an energy conserving form of the cell-center KE used with the vector invariant equations is

$$
KE_{c_e} = \frac{1}{A_c}\sum_{n_{e_i}}\frac{u_e^2\,l_e\,d_e}{4},
$$ (eq:4.17)

where $l_e$ is the length of edge $e$, $d_e$ is the length between the cell centers sharing edge $e$, and $A_c$ is the area of cell $c$.

Testing during the development of MPAS-A revealed that noise would develop in 3D atmospheric solutions in idealized baroclinic wave simulations. Skamarock et al. (2012) demonstrate that an augmented form of {eq}`eq:4.17` would remove this noise. This formulation uses a definition of the kinetic energy at the vertices of a cell

$$
KE_v = \frac{1}{A_v}\sum_{n_{e_v}}\frac{u_e^2\,l_e\,d_e}{4},
$$ (eq:4.18)

where $A_v$ is the area of the dual-mesh triangle containing the vertex, and the summation is over the three edges $n_{e_v}$ meeting at the vertex. A weighted sum of the $KE_v$ is next used to compute a cell-center kinetic energy:

$$
KE_{c_v} = \frac{1}{A_{c_v}}\sum_{n_{v_c}}A_{c_v}KE_v
$$ (eq:4.19)

Using this definition the kinetic energy at the cell center $KE_c$ is computed as the weighted sum of {eq}`eq:4.17` and {eq}`eq:4.19`:

$$
KE_c = \alpha\,KE_{c_e} + (1-\alpha)\,KE_{c_v}.
$$ (eq:4.20)

MPAS-A uses a value $\alpha = 0.375$ in its default configuration.

:::{admonition} MPAS code
:class: note

The kinetic energy is computed each RK3 substep in subroutine `atm_compute_solve_diagnostics`. The coefficient $\alpha$ in {eq}`eq:4.20` is set in the code before {eq}`eq:4.20` is evaluated. The subroutine is found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

### 4.3.2 Nonlinear Coriolis Term

The nonlinear Coriolis term in {eq}`eq:3.6`, $-\eta\,\mathbf{k}\times\mathbf{v}_{\mathbf{H}}$, requires the reconstruction of the velocity tangent to the cell edge at the points where the prognostic normal component of velocity is defined. Additionally, the full vertical vorticity $\eta$ also must be computed. The reconstruction of the tangential velocity at a edge is described in Thuburn et al. (2009), and the reconstruction represented a solution to a longstanding problem centered on using a C-grid staggering of the velocities on icosahedral meshes. Using this reconstruction, Ringler et al. (2010) derives fully nonlinear forms of the nonlinear Coriolis term for a shallow water equations discretization possessing various conservation properties. In the following subsections we describe the Thuburn et al. reconstruction, the full vertical vorticity calculation, and the construction of the nonlinear Coriolis term that we have adopted from Ringler et al.

**Tangential Velocity Reconstruction**

Consider the mesh depicted in Figure 2.1. The edge tangent velocities $v_{13} = u_{13}^\perp$, $v_{14}$ and $v_{15}$ are the normal velocities on the dual-mesh triangular cell (dashed lines) with vertices $A$, $B$, and $C$. The Thuburn et al. (2009) reconstruction is designed so that the divergence on a dual mesh triangular cell is equal to the area-weighted sum of the divergences from the three CVT (primary mesh) cells sharing the vertex or triangular cell. Thuburn et al. determined that the reconstruction of the tangential velocity at an edge satisfying this constraint could be constructed from the weighted sum of the edge-normal velocities from the cells sharing the edge. Returning to Figure 2.1, the tangential velocity $v_{13}$, for example, is composed of the weighted sum of velocities on the edges of of cells A and C:

$$
\begin{aligned}
v_{13} = {}&\bigl(w_{1,13}u_1 + w_{2,13}u_2 + w_{3,13}u_3 + w_{4,13}u_4 + w_{14,13}u_{14}\bigr) \\
&+ \bigl(w_{15,13}u_{15} + w_{9,13}u_9 + w_{10,13}u_{10} + w_{11,13}u_{11} + w_{12,13}u_{12}\bigr).
\end{aligned}
$$ (eq:4.21)

Note that the normal velocity at the edge, $u_{13}$, does not contribute to the tangential velocity given that it is perpendicular to the edge. The first term in parentheses on the right-hand-side of {eq}`eq:4.21` are the contributions from the edges of cell $A$ and the second are from cell $C$.

The general formula for the weights $w_{e,e'}$, where $e$ is an edge on the cells sharing the edge $e'$ where the tangential velocity is reconstructed, is given in Thuburn et al. equation (33):

$$
w_{e,e'}\,t_{e,v} = \left(\sum_v R_{i,v} - 1/2\right) n_{e,i}.
$$ (eq:4.22)

This formula is for tangential velocities defined by the right-hand-rule where $\mathbf{k}\times\mathbf{n}_\mathbf{e}$ gives the direction for the positive tangential velocity ($\mathbf{n}_\mathbf{e}$ is the positive unit vector for the edge normal velocity). In {eq}`eq:4.22` the indicator function $n_{e,i}$ is equal to 1 if the edge normal unit vector (direction of positive flow) points outward from cell $i$ for edge $e$, and it is equal to -1 if the normal points inward. The indicator function $t_{e,v}$ is equal to 1 when the vertex $v$ is on the left end of edge $e$ and is equal to -1 when it is on the right edge. The quantities $R_{i,v}$ in the summations are the area of each *kite* associated with vertex $v$ in the cell $i$ divided by the area of the cell $i$. Consider cell $C$ in Figure 2.1. The kite areas associated with a vertex for a given cell is the area of the dual-mesh triangle centered at that vertex that is contained in a given cell. Thus $R_{C,a} = A_{C,a}/A_C$, $R_{C,j} = A_{C,j}/A_C$, and likewise for the remaining vertices $k$, $l$, $m$, and $b$ for cell $C$. The summation of $R_{i,v}$ over all the vertices for a given cell will equal 1. In {eq}`eq:4.22`, the summation for a given edge $e$ is over the vertices traversed to reach the edge $e$ starting from edge $e'$ where the tangential velocity is located. The direction traversing the cell from $e'$ to $e$ can be clockwise or counterclockwise - both will give the same result.

The reader is advised to consult Thuburn et al. (2009) for a detailed description of the derivation of {eq}`eq:4.22` and the constraints underlying it.

:::{admonition} MPAS code
:class: note

The weights $w_{e,e'}$ are computed as one of the final steps in the MPAS mesh generation; the evaluation of {eq}`eq:4.22` is not part of the MPAS release. We are considering putting the calculation into the release in the future. The reconstruction of the tangential velocity, e.g. {eq}`eq:4.21`, takes place in subroutine `atm_compute_solve_diagnostics` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The weights $w_{e,e'}$ are stored in the array `weightsOnEdge(e,e')`, the source edges $e$ are stored in the integer array `edgesOnEdge(e,e')` and the number of edges contributing to each reconstruction is stored in the array `nEdgesOnEdge(e')`.
:::

**Vertical Vorticity**

The discrete vertical vorticity is evaluated at the vertices of the MPAS CVT mesh. Referring to Figure 2.1, the vertical vorticity at vertex $a$ is computed using the circulation theorem applied to the MPAS dual triangular mesh. The relative vertical vorticity at vertex $a$ is computed as

$$
\zeta_a = \frac{u_{13}|\overrightarrow{CA}|t_{13,a} + u_{14}|\overrightarrow{AB}|t_{14,a} + u_{15}|\overrightarrow{BC}|t_{15,a}}{A_a},
$$ (eq:4.23)

where $A_a$ is the area of the triangle centered at vertex $a$. The indicator function $t_{e,v}$ corrects for the direction of a positive velocity from the first cell center to the second in the circulation calculation. In the example {eq}`eq:4.23`, the indicator functions are all equal to 1. In contrast, the indicator functions would be equal to -1 for these same velocities contributing to the circulation about vertices $b$, $f$ and $j$. Also note that these computations are performed on MPAS horizontal surfaces. Finally, the absolute vorticity at the vertices $\eta_a = \zeta_a + 2\Omega\sin(\phi)$, where $\Omega$ is the angular velocity of the earth's rotation and $\phi$ is the latitude.

:::{admonition} MPAS code
:class: note

The vertical vorticity is computed in subroutine `atm_compute_solve_diagnostics` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The relative vertical vorticity at a vertex $\zeta_v$ is stored in the array `vorticity(level,vertex)`. The indicator function $t_{e,v}$ is stored in the array `edgesOnVertex_sign(edge,vertex)`.
The absolute vorticity at the vertices $\eta_v$ is stored in the array `pv_vertex(level,vertex)`. The name goes back to MPAS-A's initial implementation of a shallow water model.
:::

For the discrete nonlinear Coriolis term, the vorticity is needed on the edges where the horizontal velocities are defined. There are two options for defining this edge vorticity. The first option is to average the vorticities from the vertices at the ends of an edge. For example, for the edge containing $u_{13}$ in Figure 2.1, we would average the vorticities from the vertices $a$ and $b$. The second option is to use an upwind transport estimate of the vorticity (Ringler et al. 2010). The process of producing the upwind estimate begins by defining a cell-center vorticity:

$$
\eta_c = \frac{1}{A_c}\sum_v \eta_v A_{c_v}
$$

The sum is over the vertices $v$ of a cell $c$ where $A_c$ is the area of cell $c$ and $A_{c_v}$ are the kite areas in cell $c$ associated with vertex $v$. The edge absolute vorticity is then updated with an advective-form upwind transport:

$$
\eta_e = \overline{\eta_v}^e - \beta\,\Delta t\,\bigl[u_e\delta_{c_e}\eta_c + v_e\delta_{v_e}\eta_v\bigr],
$$ (eq:4.24)

where $\delta_{c_e}$ is the discrete derivative across the edge between cell centers {eq}`eq:4.14` and $\delta_{v_e}$ is the discrete derivative along the edge between the cell vertices. For example, for the edge containing $u_{13}$ in Figure 2.1, the upwinding is

$$
\eta_{13} = \frac{\eta_a + \eta_b}{2} - \beta\,\Delta t\!\left[u_{13}\frac{\eta_A - \eta_C}{|\overrightarrow{CA}|} + v_{13}\frac{\eta_a - \eta_b}{|\overrightarrow{ba}|}\right].
$$ (eq:4.25)

This upwind estimate has a damping effect on the absolute vorticity and it was first considered by Sadourny and Basdevant (1985) and is called the Anticipated Potential Vorticity Method (AVPM). The damping is most appropriate when the flow is characterized by a downscale enstrophy cascade characteristic of large-scale atmospheric flow and ocean dynamics. See (Ringler et al. 2010) and references therein for further information. The default value for $\beta$ is 0.5 in MPAS-A, although we often configure applications with $\beta = 0$ without significant impact on the solutions. It should be noted that on variable-resolution meshes some edge lengths may be very small, hence the discrete derivative along the edge between the vertices may violate stability constraints for the estimate {eq}`eq:4.24`. In these cases we set $\beta = 0$.

:::{admonition} MPAS code
:class: note

The vertical vorticity at edges are computed in subroutine `atm_compute_solve_diagnostics` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The absolute vertical vorticity is stored in the array `pv_edge(level,edge)`. The upwind estimate is also computed in this routine where first `pv_cell(level,cell)` is computed and then {eq}`eq:4.24` is evaluated.
:::

**MPAS-A Formulation for the Nonlinear Coriolis Term**

(Ringler et al. 2010) derived a discretization of the nonlinear Coriolis term that conserved potential vorticity in a shallow water model using the SCVT mesh on the sphere. We have adapted the shallow-water formulation to a 3D nonhydrostatic equations used in MPAS:

$$
\bigl[\eta\,\mathbf{k}\times\mathbf{V}_{\mathbf{H}}\bigr]_{e'} = \overline{\tilde{\rho}_d}^{e'}\sum_e \frac{1}{2}(\eta_{e'} + \eta_e)\,w_{e,e'}\,u_{e,e'}.
$$ (eq:4.26)

Here the summation is over the edges that contribute to the tangential velocity reconstruction, and the absolute vorticity that multiplies each component of the reconstructed velocity is the average of the absolute vorticities at the edge contributing to the reconstruction and the edge of the reconstruction.

:::{admonition} MPAS code
:class: note

The nonlinear Coriolis term is computed in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The absolute vertical vorticity is computed in subroutine `atm_compute_solve_diagnostics` as noted earlier. The nonlinear Coriolis term is added directly to the tendency for the momentum at the edge. Note that the reconstructed tangential velocity at the edges is not directly used in {eq}`eq:4.26`.
:::

### 4.3.3 Vertical Transport for the Horizontal Velocity

The vertical flux divergence for the horizontal velocity is evaluated at the edges MPAS-A mesh where the edge-normal velocities are defined. The flux divergence is written in a similar way to that for scalars {eq}`eq:4.6`, except the vertical mass flux $\Omega$ is averaged to the cell edges, i.e.

$$
(\tilde{\rho}_d u)^{\mathrm{new}}_k = (\tilde{\rho}_d u)^{\mathrm{old}}_k \cdots - \Delta t\,\frac{\overline{(\overline{\Omega}^e\,\overline{u})}_{k+1} - \overline{(\overline{\Omega}^e\,\overline{u})}_k}{\Delta\zeta_w(k)},
$$ (eq:4.27)

The flux is evaluated as in {eq}`eq:4.9` with $u$ replacing $\phi$ and $\overline{\Omega}^e$ replacing $\Omega$, and of course replacing cell $i$ with edge $e$. As in the scalar transport case, a nonzero positive value of $\beta$ introduces damping into the scheme. The boundary conditions for the vertical flux are handled in the same manner as in the scalar transport described in section 4.2.3.

:::{admonition} MPAS code
:class: note

The vertical flux divergence for the horizontal velocity is computed in subroutine `atm_compute_dyn_tend` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

### 4.3.4 Horizontal Mass Flux Divergence

The horizontal mass flux divergence term in the horizontal momentum equation {eq}`eq:3.6` is discretized as

$$
-\mathbf{v}_H\nabla_\zeta\cdot\mathbf{V} \;\to\; -\mathbf{v}_{H_e}\overline{\frac{1}{A_i}\sum_{n_{e_i}}\bigl(L_{e_i}\mathbf{V}_{H,e_i}\cdot\mathbf{n}\bigr)}^{c_e}.
$$ (eq:4.28)

Here the horizontal mass flux divergence in the two cells sharing the edge is averaged to the edge and multiplied by the edge velocity.

:::{admonition} MPAS code
:class: note

The discrete term {eq}`eq:4.28` is computed in subroutine `atm_compute_dyn_tend` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

### 4.3.5 Vertical Momentum Equation

**Pressure Gradient and Body Force Term**

The vertical discretization of the acoustic perturbation equations is discussed in Appendix A. The RK3 integration requires evaluation of the pressure gradient, body force and transport terms, as given in {eq}`eq:3.31`, for the acoustic integration, and its form is analogous to the perturbation terms in the acoustic-step discretization:

$$
\frac{\rho_d}{\rho_m}\!\left(\frac{\partial p'}{\partial\zeta} + g\,\tilde{\rho}'_m\right) \;\to\; \overline{\left(\frac{\rho_d^t}{\rho_m^t}\right)}^{k'}\!\left(\frac{p(k) - p(k-1)}{\Delta\zeta(k)} + g\,\overline{\tilde{\rho}'_m}^k\right).
$$ (eq:4.29)

The operator $\overline{\phi}^{k}$ is an interpolation operator that maps $\phi$ to the interface level $k$ using either the linear interpolation {eq}`eq:2.9` or a cell average {eq}`eq:2.10`. The operator quantity $\overline{(\rho^t_d/\rho^t_m)}^{k'}$ is handled in an analogous way to how it is discretized in the horizontal momentum equation {eq}`eq:4.16` except we use a vertical average of the summed moist mixing ratios:

$$
\overline{\left(\frac{\rho_d}{\rho_m}\right)}^{k'} = \overline{\left(\frac{1}{1 + \sum q_j}\right)}^{k'} \;\to\; \overline{\left(\frac{1}{1 + \sum q_j}\right)^k},
$$ (eq:4.30)

where the interpolation of the moist variables to the interface levels uses {eq}`eq:2.9` or {eq}`eq:2.10`.

:::{admonition} MPAS code
:class: note

The pressure gradient and body force terms {eq}`eq:4.29` are computed in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The moist coefficient {eq}`eq:4.30` is computed in subroutine `atm_compute_moist_coefficients`.
:::

**Transport for the Vertical Velocity**

The 3D transport of the vertical velocity in {eq}`eq:3.31` requires evaluating the term $(\nabla\cdot\mathbf{V}w)_\zeta$, and this is accomplished in a manner similar to the transport of scalar variables described in section 4.2. The 3D flux divergence for scalars {eq}`eq:4.6` is modified to account for the staggered vertical velocity $w$ and mass flux $\Omega$:

$$
-(\nabla\cdot\mathbf{V}w)_\zeta \approx -\,\frac{1}{A_i}\sum_{n_{e_i}}\!\left[L_{e_i}\bigl(\overline{\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}}^k\,w_{e_i}\bigr)\right] + \frac{\overline{\Omega}^{k,i}\,\overline{w}^{k,i} - \overline{\Omega}^{k-1,i}\,\overline{w}^{k-1,i}}{\Delta\zeta_k}.
$$ (eq:4.31)

Here the terms $\overline{\Omega}^{k,i}$ and $\overline{\Omega}^{k-1,i}$ in the vertical flux divergence denote an average of $\Omega$ to the $k$ and $k - 1$ *layers*, i.e. the layers above and below the interface $k$ (see Figure 2.7). In this case the flux is computed as

$$
\begin{aligned}
\overline{\Omega}^{k,i}\,\overline{w}^{k,i} = \overline{\Omega}^{k,i}\biggl\{&\frac{1}{12}\bigl[7(w_{k+1,i} + w_{k,i}) - (w_{k+2,i} + w_{k-1,i})\bigr] \\
&+ \mathrm{sign}(\overline{\Omega}^{k,i})\,\frac{\beta}{12}\bigl[(w_{k+2,i} - w_{k-1,i}) - 3(w_{k+1,i} - w_{k,i})\bigr]\biggr\}.
\end{aligned}
$$

The evaluation of the horizontal portion of the flux divergence occurs in the exact same manner as that for scalars as described in section 4.2.1 except that the horizontal velocities are vertically interpolated to the interface $k$, as denoted by the notation $\overline{\mathbf{V}_{\mathbf{H},\mathbf{e_i}}\cdot\mathbf{n}}^k$. The vertical interpolation uses either the linear interpolation {eq}`eq:2.9` or the cell integrated average {eq}`eq:2.10`. The boundary conditions for the vertical flux are handled in the same manner as in the scalar transport described in section 4.2.3.

:::{admonition} MPAS code
:class: note

The horizontal and vertical transport for the vertical velocity is computed in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::
