# Chapter 7: Regional Configuration

The capability to perform regional simulations with MPAS-Atmosphere in real-data configurations over only a portion of the earth, where lateral boundary conditions are prescribed, was introduced in the 8 June 2019 release of MPAS Version 7.0. The regional capability is described in Skamarock et al. (2018). Here we review the methodology of the regional configuration from this reference, provide updates that have been introduced since the initial release, and provide specific notes on the implementation in the MPAS code.

This section begins with the definition of the regions of the horizontal MPAS-Atmosphere mesh where the lateral boundary conditions are implemented, followed by a description of where the boundary conditions are applied within the time integration.

## 7.1 Specified and Relaxation Boundary Regions

Regional MPAS meshes can be constructed in two ways. The first method is to cut the mesh from a global mesh (or another regional mesh) using tools described in the MPAS User's Guide. Examples of such meshes are given in Figure 7.1 that illustrate the flexibility one has in defining an MPAS regional mesh. The boundary zone cells that are in the relaxation region and the specified region are those that are not deep blue in color, and these cells make up the boundary region where the lateral boundary conditions are applied. The utility that extracts the regional mesh defines the interior cells (deep blue) and then builds the boundary region cells.

**[Figure 7.1: Skamarock et al. (2018) Figure 2 — Regional MPAS meshes extracted from a global MPAS mesh by (a) specifying a circular region and (b) specifying a polygon whose faces are great-circle arcs on the sphere. The regional mesh specifications are given by the black lines, and the blue cells fall in the mesh. The boundary-zone cells appended to the outside of the specified mesh are also shown. To be added next session.]**

A closer view of the boundary region cells is given in Figure 7.2. The interior cells have index number 0 (not shown in the figure) and are deep blue in color, and the boundary region cells are different colors and have indices numbered sequentially beginning from 1, starting with those adjacent to (sharing an edge with) an interior cell to those at the outermost edge of the region (having edges possessing no neighbor). Edges are also labelled where the edge index number is that of the lowest cell number of the cells sharing the edge.

In the C-grid staggering of prognostic variables in MPAS (see Chapter 2 Figure 2.1), the cell edge normal horizontal velocity is prognosed and all other prognostic variables are defined at cell centers. On the outer two cell and edge boundary-zone index values are specified from an external analysis or forecast, usually by spatial and temporal interpolation using MPAS-Atmosphere utilities (see the MPAS User's Guide). This is called the *specified zone*. In Figure 7.2 the edges and cells with indices 6 and 7 would have their values specified. The cells and edges with index values 1 through 5 are updated using the standard solution update in the time integration augmented by filtering a relaxation based on the driving analysis/forecast. Cells and edges 1 through 5 are in this *relaxation zone*.

**[Figure 7.2: Skamarock et al. (2018) Figure 3 — Boundary zone for the MPAS horizontal mesh. The blue cells are the interior mesh, and the labeled cells and edges are in the boundary region. To be added next session.]**

There are two additional relaxation/filter terms applied to cells and edges in the relaxation zone. For a prognostic variable $\phi$:

$$
\frac{\partial\phi}{\partial t} = RHS_\phi + F_1(\phi_{LS} - \phi) - F_2\Delta x^2\nabla^2(\phi_{LS} - \phi).
$$ (eq:7.1)

The variable $\phi_{LS}$ represents the analysis or other forecast driving the lateral boundary conditions. The first additional term on the right-hand-side of {eq}`eq:7.1` is a Rayleigh damping term with

$$
F_1 = \gamma_1(i - 1)/m,
$$

where $i$ is the cell or edge index and $m$ is the maximum index of the relaxation zone (5 for the regional mesh shown in Figure 7.2). The second term is a 2nd-order filter that is applied to the difference between the driving values and the solution in the relaxation region. The coefficient is set as

$$
F_2 = \gamma_2(i - 1)/m,
$$

The default configuration for MPAS sets $\gamma_1 = (0.06\Delta x)^{-1}$ and has units of $\mathrm{s}^{-1}$, and the coefficient $\gamma_2 = (0.3\Delta x)^{-1}$ and also has units of $\mathrm{s}^{-1}$. The coefficients $\gamma_1$ and $\gamma_2$ are set based on the local mesh size $\Delta x$ so that the filters are appropriately set for a variable-resolution mesh. See Skamarock et al. (2018) section 2b for details on the configuration of these coefficients.

## 7.2 Time Integration

The time integration sequence in MPAS-Atmosphere is described in chapter 3 and summarized in Figure 3.1. The additions to the MPAS-A solver to accommodate the regional configuration include code sections where lateral boundary conditions and relaxation zone filters are applied, along with masks and a few places where conditionals are used to change the algorithms for the lateral boundaries.

**[Figure 7.3: Pseudo code describing the locations of lateral boundary condition applications in the MPAS timestep. To be added next session.]**

Figure 7.3 expands Figure 3.1 to indicate the locations in the MPAS time integration where lateral boundary adjustments are applied. In the regional boundary zone section (1) in Figure 7.3, the tendency for the specified zone increments for $u$, $\Theta_m$ and $\tilde{\rho}_d$ are set, and the filter terms are applied.

:::{admonition} MPAS code
:class: note

The time integration code is found in subroutine `atm_srk3` in `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. Boundary zone section (1) is found early in the named `RK3_DYNAMICS` loop and is in an `if (config_apply_lbcs)` conditional code section. It consists of calls to subroutine `atm_bdy_adjust_dynamics_speczone_tend` and subroutine `atm_bdy_adjust_dynamics_relaxzone_tend` that modify the tendencies of $u$ and $\theta_m$ using {eq}`eq:7.1`. Tendencies and states from the driving analysis are obtained through `mpas_atm_get_bdy_tend` and the `mpas_atm_get_bdy_state` generic interface, implemented in `src/core_atmosphere/dynamics/mpas_atm_boundaries.F`.
:::

The second boundary zone section occurs after the acoustic step within a specific Runge-Kutta substep. In section (2), the values for $u$ and $\tilde{\rho}u$ are reset to those of the appropriate time level (e.g. $u^*$, $u^{**}$, or $u^{t+\Delta t}$, see chapter 3 equations {eq}`eq:3.12`–{eq}`eq:3.14` and Figure 3.3) after the full variables are reconstructed after the acoustic steps. It is also in this part of the dynamics timestep that the lateral boundary conditions for $w$ are applied to set the vertical velocity in the specified zone. In the current release of MPAS-Atmosphere, $w$ in the specified zone is set such that the horizontal gradient of $w$ is zero for the flux computation across an edge between a relaxation-zone and specified-zone cell. The current implementation sets the specified-zone cell value of $w$ to that of the nearest relaxation-zone cell value on that level.

:::{admonition} MPAS code
:class: note

The code for resetting the $u$ and $\tilde{\rho}u$ values is found in the named `RK3_DYNAMICS` loop after the acoustic-step loop and after the call to subroutine `atm_recover_large_step_variables`, and is in an `if (config_apply_lbcs)` conditional code section. The $w$ values in the specified zone are set after a call to subroutine `atm_compute_solve_diagnostics` and the values are set in a call to subroutine `atm_zero_gradient_w_bdy` in a separate `if (config_apply_lbcs)` conditional code section immediately before the end of the `RK3_DYNAMICS` loop. This subroutine is in the time integration file `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

The third lateral boundary condition section is in the scalar transport Runge-Kutta loop immediately after the call to advance the scalars over a Runge-Kutta substep. Here the values in the specified zone are set and the filtering in the relaxation zone using {eq}`eq:7.1` is applied if the values from the driving analysis exist. If the values do not exist then they are assumed to be zero and the specified zone values are set to zero with the filters applied using this value.

:::{admonition} MPAS code
:class: note

The lateral boundary condition code for scalar transport, section (3), is in an `if (config_apply_lbcs)` conditional code section within the named loop `RK3_SPLIT_TRANSPORT`. In the current implementation individual scalars are checked to see if they exist, and if so the driving analysis state for them is acquired. The relaxation zone filtering and the specified zone values are set in a call to `atm_bdy_adjust_scalars`. All this code is in the time integration file `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

The final boundary condition section in the MPAS time integration resets values for $\theta_m$, $\Theta_m$ and the scalars after the call to microphysics. The microphysics is applied in all the columns, even in the specified zone, hence the need to reset the values. If the driving analysis does not contain values for a given scalar then that scalar is set to zero.

:::{admonition} MPAS code
:class: note

The lateral boundary condition code for scalar transport, section (4), is in an `if (config_apply_lbcs)` conditional code section and appears immediately after the call to `microphysics_driver`. Specified zone values for $\theta_m$ and $\Theta_m$ are reset in a call to `atm_bdy_reset_speczone_values` and the scalars are reset in a call to subroutine `atm_bdy_set_scalars`. All this code, including the subroutines, is in the time integration file `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

### 7.2.1 Divergence Damping in the Relaxation Zone

For the horizontal momentum, the 2nd-order horizontal filter applied in {eq}`eq:7.1` in the relaxation region takes the form given in {eq}`eq:5.7`. However, we have introduced a coefficient to increase the divergent component of the filtering, similar to additional filtering applied to the horizontal filter for the 4th-order background filter given in section 5.2 equation {eq}`eq:5.12`. The resulting 2nd-order filter applied to the horizontal momentum in the relaxation zone takes the form

$$
\begin{aligned}
\frac{\partial(\tilde{\rho}_du_i)}{\partial t} = {}&\cdots - \tilde{\rho}\frac{\gamma_2(i-1)}{m}\Delta x^2\!\left(\beta_d\frac{\partial}{\partial x_i}\nabla_H\cdot(\mathbf{v}_{LBC} - \mathbf{v})\right. \\
&\left.\hphantom{\cdots - \tilde{\rho}\frac{\gamma_2(i-1)}{m}\Delta x^2\!\bigl(} - \frac{\partial(\zeta_{LBC} - \zeta)}{\partial x_j}\right).
\end{aligned}
$$ (eq:7.2)

:::{admonition} MPAS code
:class: note

The relaxation zone filtering for the horizontal momentum is performed in section (1) of the time integration (see Figure 7.3 and the discussion earlier in this section). The coefficient $\beta_d$ is configurable in `namelist.atmosphere` through the namelist variable `config_relax_zone_divdamp_coef`. Its default value is 6. The filtering is performed in subroutine `atm_bdy_adjust_dynamics_relaxzone_tend` in `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::
