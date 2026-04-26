# Chapter 5: Filters

In this chapter we describe the discretization of filters employed in MPAS. Not covered in this chapter are filters described elsewhere:

- Vertical velocity damping in the gravity-wave absorbing layer: See section 3.4.3.
- Acoustic wave filter: See section 3.4.3.
- Upwinding in the horizontal transport scheme for $\theta$, $w$ and scalars: See section 4.2.1.
- Upwinding in the vertical transport: Scalars and potential temperature in section 4.2.2; $w$ in section 4.3.5; horizontal velocity $\mathbf{v}$ in section 4.3.3.

## 5.1 Horizontal Filters for Cell-Centered Variables

There are two horizontal filters available in MPAS for filtering cell-centered variables, a 2nd-order filter employing a horizontal Laplacian and a 4th-order filter that employs two applications of the Laplacian. Both of these filters are discretized on the $\zeta$ surfaces of the MPAS-Atmosphere mesh, and we use the term "horizontal" to describe the smoothing on these coordinate surfaces. In the current MPAS release the 2nd-order filter can be applied to the potential temperature $\theta$ and the vertical velocity $w$. Horizontal filters are not applied to other scalars (e.g. $q_v$, etc) because monotonic transport supplies any necessary filtering.

### 5.1.1 2nd-Order Horizontal Filter

The continuous form for the second-order horizontal filter is given by the RHS term for the update equation for a scalar $\phi$:

$$
\frac{\partial(\tilde{\rho}_d\phi)}{\partial t} = \cdots + \nabla_\zeta\cdot\tilde{\rho}_d K\,\nabla_\zeta\phi,
$$ (eq:5.1)

where the variable $K$ is an eddy viscosity with units $\mathrm{m}^2\,\mathrm{s}^{-1}$. The discrete form of Laplacian in {eq}`eq:5.1` for $\theta$ (a variable that is not vertically staggered) is

$$
\nabla_\zeta\cdot\tilde{\rho}_d K\,\nabla_\zeta\theta \;\to\; \frac{1}{A_i}\sum_{e_i}L_{e_i}\,\overline{\tilde{\rho}_d}\,\overline{K}\,(\mathbf{n}_{e_i}\cdot\nabla\theta).
$$ (eq:5.2)

The overbar terms in {eq}`eq:5.2` represent an average of a quantity from the cell centers of the two cells sharing the edge $e_i$. The gradient operator $\nabla\theta$ is the gradient of the quantity $\theta$ normal to the edge in the direction outward from the cell ($\mathbf{n}_{e_i}$ points outward away from cell $i$). Compared to the horizontal mass flux used for transport as given in {eq}`eq:4.6`, the mass flux $\mathbf{V}_{\mathbf{H},\mathbf{e_i}}$ is replaced by $\tilde{\rho}_d K$ and the transported scalar $\theta$ replaced by its gradient normal to the cell edge $\mathbf{n}_{e_i}\cdot\nabla\theta$. This turbulent flux divergence is constructed such that a constant value of $\theta$ remains constant because $\nabla\phi \equiv 0$. In addition, the filters conserved the quantity $\tilde{\rho}_d\theta$ locally and hence globally.

The 2nd-order horizontal filter can also be applied to the vertical velocity $w$. In this case the discrete of the turbulent flux divergence {eq}`eq:5.2` is given by

$$
\nabla_\zeta\cdot\tilde{\rho}_d K\,\nabla_\zeta w \;\to\; \frac{1}{A_i}\sum_{e_i}L_{e_i}\,\overline{\tilde{\rho}_d}^{H,z}\,\overline{K}^{H,z}\,(\mathbf{n}_{e_i}\cdot\nabla w),
$$ (eq:5.3)

where the averaging of $\tilde{\rho}_d$ and $K$ is over the horizontal (between the cells sharing the edge) and vertical (between the model layers sharing the interface where $w$ is defined).

### 5.1.2 4th-Order Horizontal Filter

The continuous form for the fourth-order horizontal filter is given by the RHS term for the update equation for a scalar $\phi$:

$$
\frac{\partial(\tilde{\rho}_d\phi)}{\partial t} = \cdots - \nabla\cdot\bigl(\tilde{\rho}_d\nu_4\nabla(\nabla_\zeta\cdot\nabla_\zeta\phi)\bigr),
$$ (eq:5.4)

where the hyperviscosity $\nu_4$ has units $\mathrm{m}^4\,\mathrm{s}^{-1}$. Defining a discrete turbulent flux divergence operator $T_2(\phi)$ as

$$
\nabla_\zeta\cdot\nabla_\zeta\phi \;\to\; \frac{1}{A_i}\sum_{e_i}L_{e_i}\,(\mathbf{n}_{e_i}\cdot\nabla\phi) = T_2(\phi),
$$ (eq:5.5)

we can write the discrete 4th-order filter term as

$$
\nabla\cdot\bigl(\tilde{\rho}_d\nu_4\nabla(\nabla_\zeta\cdot\nabla_\zeta\phi)\bigr) \;\to\; \frac{1}{A_i}\sum_{e_i}L_{e_i}\,\overline{\tilde{\rho}_d}\,\nu_4\,(\mathbf{n}_{e_i}\cdot\nabla T_2(\phi)).
$$ (eq:5.6)

The hyperviscosity $\nu_4$ is constant in the current MPAS implementation, thus there is no averaging operator applied to it. As in the second-order filter, the density at the cell edges is an average of the densities from the cell centers of the two cells sharing the edge for vertically unstaggered variables ($\theta$). The density will also be averaged in $z$ when the filter is applied to $w$.

:::{admonition} MPAS code
:class: note

The second-order and fourth-order horizontal filtering for $\theta$ and $w$ are computed in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The computations occur only once in each dynamics timestep, during the first RK3 substep. When the 4th-order filter is applied, the turbulent flux divergences from the 2nd-order filter application {eq}`eq:5.5` are saved and used in the second application of the turbulent flux divergence operator {eq}`eq:5.6` needed to complete the 4th-order operator.
:::

## 5.2 Horizontal Momentum Equation: Horizontal Filters

We cannot easily cast the horizontal filters for the momentum in terms of a flux divergence. Instead, we employ the vector identity to cast the Laplacian of the horizontal velocity in terms of the gradient of divergence across the cell edge minus the gradient of the relative vertical vorticity along the edge:

$$
\begin{aligned}
\frac{\partial(\tilde{\rho}_du_i)}{\partial t} &= \cdots + \tilde{\rho}_dK\nabla^2 u_i \\
&= \cdots + \tilde{\rho}_dK\!\left(\frac{\partial}{\partial x_i}\nabla_H\cdot\mathbf{v} - \frac{\partial\zeta}{\partial x_j}\right),
\end{aligned}
$$ (eq:5.7)

where the relative vertical vorticity $\zeta = \mathbf{n}\cdot\nabla\times\mathbf{v}$ and $\mathbf{n}$ is a unit vector normal to the $\zeta$ coordinate surface.

The fourth-order filter for the horizontal velocity employs two applications of the Laplacian given on the right-hand-side of {eq}`eq:5.7`. If we define an operator for the Laplacian for the horizontal velocity as

$$
L(\mathbf{v}) = \frac{\partial}{\partial x_i}\nabla_H\cdot\mathbf{v} - \frac{\partial}{\partial x_j}\mathbf{n}\cdot\nabla\times\mathbf{v}
$$ (eq:5.8)

we can express the filter as

$$
\begin{aligned}
\frac{\partial(\tilde{\rho}_du_i)}{\partial t} &= \cdots - \tilde{\rho}_d\nu_4 L^2(\mathbf{v}) \\
&= \cdots - \tilde{\rho}_d\nu_4\!\left(\frac{\partial}{\partial x_i}\nabla_H\cdot L(\mathbf{v}) - \frac{\partial}{\partial x_j}\mathbf{n}\cdot\nabla\times L(\mathbf{v})\right)
\end{aligned}
$$ (eq:5.9)

**Discrete Vorticity-Divergence Form of the Laplacian**

The discrete form of the velocity divergence used in {eq}`eq:5.8` for cell $i$ with edges $e_i$, edge lengths $L_{e_i}$, and area $A_i$, is computed using

$$
\nabla_H\cdot\mathbf{v} \;\to\; \frac{1}{A_i}\sum_{n_{e_i}}(L_{e_i}\mathbf{v}_{H,e_i}\cdot\mathbf{n}) = D_i,
$$ (eq:5.10)

where $L_{e_i}$ is the length of the edge $e_i$ and $A_i$ is the area of cell $i$. Referring to the horizontal MPAS mesh given in Figure 2.1, the gradient of the divergence to update the horizontal velocity $u_{13}$, where positive $u_{13}$ indicates flow from cell $C$ to cell $A$, is given as

$$
\frac{\partial}{\partial x_i}\nabla_H\cdot\mathbf{v} \;\to\; \frac{D_A - D_C}{|\overrightarrow{CA}|},
$$

where $|\overrightarrow{CA}|$ is the distance from cell centers $C$ to $A$.

The discrete form of the relative vertical vorticity used in {eq}`eq:5.8` is described in section 4.3.2 and is computed using {eq}`eq:4.23`. The relative vertical vorticity is defined at the vertices of the MPAS horizontal mesh and it is computed by evaluating the circulation about the dual-mesh triangle containing the vertex and then dividing the circulation by the triangle area. Referring to the horizontal MPAS mesh given in Figure 2.1, the gradient of the vorticity along the edge $ab$ for a positive horizontal velocity $u_{13}$ indicating flow from cell $C$ to cell $A$ is

$$
\frac{\partial\zeta}{\partial x_j} \;\to\; \frac{\zeta_a - \zeta_b}{|\overrightarrow{ab}|},
$$ (eq:5.11)

where $\zeta_a$ and $\zeta_b$ are the relative vertical vorticities defined at the vertices $a$ and $b$ using {eq}`eq:4.23` and $|\overrightarrow{ab}|$ the distance between the vertices $a$ and $b$, i.e. $L_{e_i}$ for the cell edge with vertex end points $a$ and $b$.

**Horizontal Divergence Damping**

In convection-permitting applications (horizontal cell spacing of a few kilometers of less) we have found that it is advantageous to more strongly filter the horizontal divergence relative to the vertical vorticity in the vorticity-divergence form of the 4th-order filter for the horizontal momentum given in {eq}`eq:5.9`. This modification to the horizontal filtering of the horizontal momentum equation is accomplished by introducing a coefficient in the filter term, modifying {eq}`eq:5.9`:

$$
\begin{aligned}
\frac{\partial(\tilde{\rho}_du_i)}{\partial t} &= \cdots - \tilde{\rho}_d\nu_4 L^2(\mathbf{v}) \\
&= \cdots - \tilde{\rho}_d\nu_4\!\left(\beta_d\frac{\partial}{\partial x_i}\nabla_H\cdot L(\mathbf{v}) - \frac{\partial}{\partial x_j}\mathbf{n}\cdot\nabla\times L(\mathbf{v})\right)
\end{aligned}
$$ (eq:5.12)

:::{admonition} MPAS code
:class: note

The second-order and fourth-order horizontal filtering for $u$ are computed in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. As with $\theta$ and $w$ filter applications, the computations occur only once in each dynamics timestep, during the first RK3 substep. When the 4th-order filter is applied, the turbulent flux divergence from the 2nd-order filter application is saved and used in the second application of the turbulent flux divergence operator that completes the 4th-order filter evaluation. The coefficient $\beta_d$ is runtime configurable and is set in the `namelist.atmosphere` as the variable `config_del4u_div_factor`.
:::

## 5.3 Eddy Viscosities and Hyper-Viscosities

There are two options for specifying the 2nd-order eddy viscosities and the 4th-order hyper-viscosities. Fixed, constant values of the viscosities can be specified using the *2d_fixed* option. Fixed viscosity (with units of $\mathrm{m}^2\,\mathrm{s}^{-1}$) and hypervisocity (with units of $\mathrm{m}^4\,\mathrm{s}^{-1}$) values are specified for momentum and theta filtering, and the default values are zero, so they must be set when configuring an MPAS application using the *2d_fixed* option. This option is used for idealized simulations, often when a converged solution is sought. The eddy viscosity and hypervisocity used in the filtering of theta is divided by the Prandtl number which has a default value of 1.

The second option for filtering employs a spatially and temporally varying 2nd-order eddy viscosity that is computed using the horizontal deformation following Smagorinsky (1963).

$$
K = c_s^2\,l^2\!\left[\left(\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}\right)^{\!2} + \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^{\!2}\right]^{1/2}.
$$ (eq:5.13)

The velocities $u$ and $v$ are orthogonal horizontal velocities where $u\times v$ produce a upward pointing vector. $c_s$ is the Smagorinsky coefficient where the default value in MPAS is $c_s = 0.125$. $l$ is a horizontal length scale and it is set to the nominal cell spacing in the highest resolution region of the MPAS Voronoi mesh. This cell spacing specified as part of the MPAS mesh file description. The computation of the deformation used in the eddy viscosity calculation {eq}`eq:5.13` is given in Appendix C.

The 4th-order filter is activated automatically when the Smagorinsky filter is specified. The constant hyperviscosity is given as

$$
\nu_4(m^4 s^{-1}) = c_4 l^3,
$$ (eq:5.14)

where $c_4$ is a user specified coefficient and $l$ is the same length scale used in the Smagorinsky scheme. The default value for the coefficient $c_4 = 0.05$.

The eddy viscosities and hyper-viscosities {eq}`eq:5.13` and {eq}`eq:5.14` are given for the uniform MPAS meshes. For variable-resolution meshes the eddy viscosities and hyper-viscosities are scaled such that they give the same physical viscosities and hyper-viscosities as would be employed at the nominal *local cell spacing*. For the 4th-order hyper-viscosity {eq}`eq:5.14`, the length scale $l$ is replaced by

$$
l = (\rho_{\mathrm{mesh}})^{-1/4},
$$ (eq:5.15)

where $\rho_{\mathrm{mesh}}$ is the mesh density function used to generate the Voronoi mesh (see Section 2). This is equivalent to scaling the hyper-viscosity computed for the highest-resolution part of the mesh by a coefficient proportional to $(l/l_f)^3$, where $l_f$ is the length scale (nominal cell spacing) on the highest resolution of the variable-resolution mesh. We use the mesh density function as opposed to some sampling of the local cell-center spacing because the Voronoi mesh cell spacing can vary locally on an unstructured mesh whereas the density function is smooth by design. This helps minimize any grid imprinting in the filter operators.

:::{admonition} MPAS code
:class: note

The coefficient $c_s$ can be set in the `namelist.atmosphere` file with the variable `config_smagorinsky_coef`. The default value is 0.125. The variable $c_4$ also a `namelist.atmosphere` variable named `config_visc4_2dsmag`. It has a default value of 0.05 and units of m/s.
:::

## 5.4 Offcentering in the Vertical Acoustic Integration

In section 3.4.2 describing the acoustic timestep, the evolution equations for the vertical mass flux $\Omega$ {eq}`eq:3.37`, the coupled potential temperature $\Theta_m$ {eq}`eq:3.38`, and the dry air density $\tilde{\rho}_d$ {eq}`eq:3.39`, with their corresponding right-hand-side terms {eq}`eq:3.40`, {eq}`eq:3.41`, and {eq}`eq:3.42`, represent a semi-implicit integration of vertically-propagating acoustic and gravity waves. In these equations the terms evaluated at the new time level $(\tau + \Delta\tau)$ are multiplied by the factor $(1 + \epsilon)/2$ and the time $\tau$ terms are multiplied with a factor of $(1 - \epsilon)/2$. For $\epsilon > 0$ this will damp acoustic and gravity waves, although given the small acoustic time step it is only the high frequency acoustic waves that experience significant damping.

The off-centering in the semi-implicit time step is perhaps most important for stabilizing the integration when the coordinate surfaces are sloped. As described in Ikawa (1988) and Dudhia (1995), a linear analysis shows that the stability of the scheme is achieved when the coordinate surface slope is less than the off-centering parameter in the vertically semi-implicit solution:

$$
\left.\frac{\partial z}{\partial x}\right|_\zeta < \epsilon,
$$

The default value for $\epsilon = 0.1$ in MPAS-Atmosphere. Higher values of the off-centering parameter $\epsilon$ may be needed for higher resolution applications with significant terrain features.

:::{admonition} MPAS code
:class: note

The value for $\epsilon$ is runtime configurable and is set in `namelist.atmosphere` through the variable `config_epssm`.
:::

## 5.5 Gravity-Wave absorbing Layer

### 5.5.1 Rayleigh Damping for w

Vertical velocity damping in the gravity-wave absorbing layer is described in section 3.4.3.

### 5.5.2 2nd-Order Horizontal Filtering

Another option for filtering vertically-propagating gravity waves is to enable a 2nd-order horizontal filter for the horizontal momentum, vertical momentum, and potential temperature, in the layers close to the model lid. The 2nd-order filter for the potential temperature and vertical velocity follows the formulation described in section 5.1.1 using equations {eq}`eq:5.1` and {eq}`eq:5.2`, and the formulation described in section 5.2 for the horizontal momentum using the vorticity-divergence form of the Laplacian. The eddy viscosity used is constant on each coordinate surface level varying linearly from a value of zero at a user-specified level to a maximum value at the model top. When the 2D Smagorinsky filter is active, the eddy viscosity used is chosen to be the maximum of either the fixed value or the value computed in the Smagorinsky formulation.

:::{admonition} MPAS code
:class: note

The 2nd-order horizontal filter is included as part of the Runge-Kutta time-integration tendencies computed in subroutine `atm_compute_dyn_tend` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The filter is enabled by setting the `namelist.atmosphere` parameter `config_mpas_cam_coef` to a non-zero value between 0 and 1. The default value is zero (off). The 2nd-order filter is applied starting at the integer `namelist.atmosphere` variable `config_number_cam_damping_levels` number of levels from the model top. The maximum value of the eddy viscosity is $8.333\times\Delta x\times$`config_mpas_cam_coef` ($\mathrm{m}^2/\mathrm{s}$). Setting `config_mpas_cam_coef` to 1 recovers the typical maximum value used in the climate configuration, while 0.2 is the typical value used in other applications.
:::

### 5.5.3 Rayleigh Damping of the Horizontal Momentum

While not technically a gravity-wave absorbing filter, a Rayleigh damping filter is available for the horizontal momentum and is applied to a specified number of layers below the model top. For a given time step $\Delta t$ the damping term has the form

$$
(\tilde{\rho}_du)^{t+\Delta t} = (\tilde{\rho}_du)^t - \frac{\Delta t}{\Delta t_d}(\tilde{\rho}_du)^t,
$$

where $\Delta t_d$ is the damping timescale. The damping varies linear between a maximum value at the model top to zero at the layer where the damping begins. The filter is used in longer integrations and helps stabilize the integration when abnormally high horizontal velocities can develop because of the artificial nature of the rigid-lid upper boundary. This filter is also used in numerous hydrostatic models using a constant pressure upper boundary condition. The default damping timescale is long (5 days), but as noted we have found this helps stabilize longer integrations.

:::{admonition} MPAS code
:class: note

The Rayleigh damping of the horizontal momentum is included as part of the Runge-Kutta time-integration tendencies computed in subroutine `atm_compute_dyn_tend` found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The filter is enabled by setting the `namelist.atmosphere` variable `config_rayleigh_damp_u` to true (false is the default). The number of layers over which the damping layer extends is set in the `namelist.atmosphere` variable `config_number_rayleigh_damp_u_levels`, and the damping timescale (days) is set in the `namelist.atmosphere` variable `config_rayleigh_damp_u_timescale_days`.
:::
