# Chapter 8: Initial Conditions

The method for producing an initial state for MPAS-Atmosphere to integrate involves executing a different MPAS model called *init_atmosphere*. The practical aspects of producing an initial state are covered in the *MPAS-Atmosphere User's Guide*. In this chapter we describe the physical configurations of the different states along with some technical aspects concerning how they are produced. Here are the initial states for the cases that can be produced by the *init_atmosphere* model along with the case numbers used within the *init_atmosphere* model:

**(1, 2, 3):** Three variants of the Jablonowski and Williamson (2006) baroclinic wave case.

**(4, 5):** A 2D $(x, z)$ squall-line and 3D supercell thunderstorm initialization in Cartesian domains.

**(6):** A 2D $(x, z)$ mountain-wave initialization.

**(7):** Initialization for full earth atmosphere simulations on the sphere or regional domains on the sphere.

**(8):** Produce sea-surface temperature (SST) and surface fields at specified times. These are typically used to update the SSTs and surface fields during long integrations of MPAS-A.

**(9):** Produce boundary conditions for regional MPAS-A simulations on the sphere.

**(13):** Produce the 3D MPAS-A mesh used in applications outside of MPAS-A integrations, e.g. in NSF NCAR's Community Earth System Model when it is employing the MPAS-A dynamical core.

:::{admonition} MPAS code
:class: note

The numbers attached to the test cases listed above correspond to the case configuration number in the `namelist.init_atmosphere` configuration file for the *init_atmosphere* core. For example in the full earth atmosphere simulations the configuration variable `config_init_case` = 7.
:::

## 8.1 Initialization Commonalities

MPAS-Atmosphere integrates an equation set that uses perturbation variables for density and potential temperature defined relative to a reference state (see chapter 3). The reference state can be specified as a continuous function or a discrete function. The perturbation state density and potential temperature, along with a specified water vapor profile, satisfy an MPAS-A discrete hydrostatic balance.

### 8.1.1 Reference States

An analytic reference state is used in the Jablonowski and Williamson (2006) baroclinic wave case and in the initialization of the full earth atmosphere, often referred to as a *real-data* initialization. The analytic state is isothermal with a reference temperature $T_o = 250$ Kelvin. The pressure, density and potential temperature are given by

$$
\begin{aligned}
\overline{p} &= p_o\exp\!\left(-\frac{gz}{R_d T_o}\right), \\
\overline{\rho} &= \frac{\overline{p}}{R_d T_o}, \\
\overline{\theta} &= T_o\!\left(\frac{\overline{p}}{p_o}\right)^{\!-R_d/c_p},
\end{aligned}
$$ (eq:8.1)

where $p_o = 1000$ hPa is a reference pressure.

### 8.1.2 Hydrostatic Balance

The hydrostatic relation cast in terms of perturbation thermodynamical variables is given in {eq}`eq:3.22`, and for the case where the only nonzero moist species is water vapor $q_v$ it can be expressed as

$$
\frac{\partial p'}{\partial\zeta} = -g\bigl(\tilde{\rho}'_d + \tilde{\rho}_d q_v\bigr).
$$

The discrete formulation used in the initialization is that used for these terms in the MPAS-A vertical momentum equation given in {eq}`eq:4.29`,

$$
p(k) - p(k-1) = -\Delta\zeta(k)\,\overline{g\bigl(\tilde{\rho}'_d + \tilde{\rho}_d q_v\bigr)}^{k},
$$ (eq:8.2)

where $k$ indicates the layer number, with layer height increasing with increasing $k$. The operator $\overline{()}^{k}$ represents either a vertical interpolation to the midpoint of the interval or the average value computed from integrating the quantity over the interval, both options described in Section 2.2.3 with weights defined in {eq}`eq:2.9` and {eq}`eq:2.10`, respectively.

Another form of the hydrostatic relation used in some of the initializations is cast in terms of the Exner function:

$$
\frac{\partial\Pi}{\partial z} = -\frac{g}{c_p\theta\,\dfrac{1+\dfrac{R_v}{R_d}q_v}{1+q_v}}
$$ (eq:8.3)

where the Exner function is defined as

$$
\Pi = \left(\frac{p}{p_o}\right)^{R_d/c_p}
$$

## 8.2 Initialization for Idealized Cases

The idealized atmosphere cases comprise cases 1 through 6 as described below.

### 8.2.1 Idealized Baroclinic Wave Case

The idealized baroclinic wave case is taken from Jablonowski and Williamson (2006) and the cases 1 through 3 involve three different perturbations imposed on the zonally uniform balanced midlatitude northern and southern hemisphere jets prescribed in the initial state. The balanced flow represents a steady-state solution to the equations of motion. However, any zonally nonuniform perturbations will project onto the baroclinically unstable modes of the jets giving rise to baroclinically unstable waves. The perturbations for the three cases are as follows:

**(1):** No perturbation is applied. In models where the mesh is not zonally uniform (e.g. MPAS), truncation errors arising from the discretization on the horizontal mesh will give rise to baroclinic waves. See Skamarock et al. (2012) figure 8 for an example of an MPAS-Atmosphere result for this case.

**(2):** The perturbation is that given in equation (10) from Jablonowski and Williamson (2006) for the horizontal wind field $u$:

$$
u'(\lambda, \phi) = u_p\exp\!\left\{-\left(\frac{r}{R}\right)^{\!2}\right\},
$$ (eq:8.4)

where $(\lambda, \phi)$ are the longitude and latitude, $u_p = 1\ \mathrm{m\,s^{-1}}$, $r$ is the great-circle arc distance on the sphere from the point $(\lambda_o, \phi_o) = (20^\circ\mathrm{E}, 40^\circ\mathrm{N})$, and $R = a/10$ where $a$ is the earth radius. There is no vertical dependence to this perturbation.

See Skamarock et al. (2012) figure 7 for an example of an MPAS-Atmosphere result for this case.

**(3):** The perturbation is designed to be the starting point for growing the most unstable normal mode for the jet which is zonal wavenumber 9:

$$
u'(\lambda, \phi) = u_p\cos[k_x(\lambda - \lambda_o)]\,\sin^2 2\phi,
$$ (eq:8.5)

where $k_x = 9$ is the zonal wavenumber of the perturbation. There is no vertical dependence to this perturbation. Also, there is not a published solution for this perturbation, but the converged normal-mode solution for wavenumber 9 is shown in Park et al. (2013).

The initialization integrates the hydrostatic balance equation {eq}`eq:8.2` from the surface upward to the lid, with a lower boundary condition of $p = 1000$ hPa. There is a gradual terrain slope from the equator to the pole to satisfy a zero zonal wind at the surface for the given analytic formula for the jet. This case uses the isothermal reference state {eq}`eq:8.1`.

:::{admonition} MPAS code
:class: note

The initialization code for this test case is found in subroutine `init_atm_case_jw` in the source code file `src/core_init_atmosphere/mpas_init_atm_cases.F`. The wind perturbations in equations {eq}`eq:8.4` and {eq}`eq:8.5` are integrated between the latitudes of the edge endpoints and projected onto (dotted with) the normal vector in the zonal direction so as to get the average edge-normal wind component on the edge. Additionally, the horizontal winds of the unperturbed jet are computed so that the horizontal mass divergence in the initial state is minimized.
:::

### 8.2.2 Squall Line and Supercell Case

Both the squall line and supercell thunderstorm case use the same initialization procedure. The lower boundary is a flat surface (there is no terrain) and vertical grid spacing is constant. In these initializations the model top is set to $z = 20$ km. The environment is initialized with the sounding taken from Weisman and Klemp (1982):

$$
\theta(z) = \theta_o + (\theta_{tr} - \theta_o)\!\left(\frac{z}{z_{tr}}\right)^{\!5/4}, \qquad \text{for } z \le z_{tr}
$$ (eq:8.6)

$$
= \theta_{tr}\exp\!\left[\frac{g}{c_p T_{tr}}(z - z_{tr})\right], \qquad \text{for } z > z_{tr}
$$ (eq:8.7)

with the relative humidity defined as

$$
H(z) = 1 - \frac{3}{4}\!\left(\frac{z}{z_{tr}}\right)^{\!5/4}, \qquad \text{for } z \le z_{tr}
$$ (eq:8.8)

$$
= 0.25, \qquad \text{for } z > z_{tr}
$$ (eq:8.9)

where

$$
\begin{aligned}
\theta_o &= 300\ \text{K}, \\
\theta_{tr} &= 343\ \text{K}, \\
T_{tr} &= 213\ \text{K}, \\
z_{tr} &= 12\ \text{km}.
\end{aligned}
$$

For the 2D squall-line simulation and the supercell simulation the environmental wind profile is set in the initialization code. It only varies in $z$ and only the $x$ component of the wind field is nonzero:

$$
\begin{aligned}
u(z) &= u_m\frac{z}{z_{ts}} - u_s \qquad \text{for } z \le z_{ts} \\
&= u_m - u_s, \qquad \text{for } z > z_{ts}
\end{aligned}
$$

where $u_m = 12$ m/s, $u_s = 10$ m/s and $z_{ts} = 2500$ m for the squall line case, and $u_m = 30$ m/s, $u_s = 15$ m/s and $z_{ts} = 5000$ m for the supercell case.

The reference state for these simulations is that given in {eq}`eq:8.6` and {eq}`eq:8.7` for a dry atmosphere $(H(z) = 0)$. The hydrostatic balance for the reference state uses the hydrostatic relation cast in terms of the Exner function {eq}`eq:8.3`. This equation is also used to compute the mixing ratios based on the relative humidities {eq}`eq:8.8` and {eq}`eq:8.9`. The hydrostatic balance for the full state uses an upward vertical integration of {eq}`eq:8.3` to estimate the pressure at the top of the domain, using a reference pressure $p_o = 10^5$ Pa as the boundary condition at the surface, and then uses a downward integration of {eq}`eq:8.2` that is used with the state equation to specify the density. This hydrostatic balance includes a warm low-level temperature perturbation (a warm bubble) to initiate convection. The integrations are applied iteratively to arrive at a converged solution that satisfies both the hydrostatic relation and the state equation. Given that the vapor mixing ratio $q_v$ is computed using the Exner-function-based hydrostatic relation, the relative humidity {eq}`eq:8.8` and {eq}`eq:8.9` will not be exactly satisfied in the MPAS atmosphere computed using the hydrostatic relation {eq}`eq:8.2`. We also note that an upper bound of $q_{v_{\max}} = 14$ g/kg is set in this initialization.

:::{admonition} MPAS code
:class: note

The initialization code for this test case is found in subroutine `init_atm_case_squall_line` in the source code file `src/core_init_atmosphere/mpas_init_atm_cases.F`. The hydrostatic balance loops and iterations do not assume a flat lower boundary, thus many of the calculations are redundant, particularly for the reference state. Also note that the variable `p` refers to the Exner function in the hydrostatic balancing of the reference state and the variable `pp` refers to a perturbation pressure in the balancing of the MPAS-A perturbation state.
:::

### 8.2.3 Mountain Wave Case

As in the squall line and supercell cases, the vertical grid has a constant $\Delta z$ and the model top is set to $z = 21$ km. The reference state potential temperature is $\overline{\theta}(z) = \theta_o = 288$ K. The initial profile for the potential temperature characterized by constant stability:

$$
\theta(z) = \theta_o(1 + N^2 z)
$$

where the stability parameter

$$
N^2 = \frac{g}{\theta}\frac{\partial\theta}{\partial z}.
$$

In the existing initialization code the value for this case is $N^2 = 10^{-4}\ \mathrm{s}^{-2}$.

The surface terrain profile is that described in Schär et al. (2002) and Klemp et al. (2003), and is given by

$$
h(x) = h_m\exp\!\left[-\left(\frac{x - x_c}{x_a}\right)^{\!2}\right]\cos^2\!\left[\pi\frac{x - x_c}{x_a}\right]
$$

where $h_m = 250$ m and $x_a = 4000$ m.

The hydrostatic balance in this case uses the same procedure as used in the squall line and supercell initialization.

:::{admonition} MPAS code
:class: note

The initialization code for this test case is found in subroutine `init_atm_case_mtn_wave` in the source code file `src/core_init_atmosphere/mpas_init_atm_cases.F`. There are a number of options for terrain profiles and stability configurations that can be changed by editing and recompiling the *init_atmosphere* core.
:::

## 8.3 Initialization for Earth-Atmosphere Simulations

The initialization for a real-data full-earth-atmosphere case is much more complex compared to the idealized cases. In the following we describe only the hydrostatic balancing and a few other numerical aspects of the initialization. Other aspects of the initialization scheme are described in the MPAS-Atmosphere User's Guide.

### 8.3.1 Initialization of the 3D State

Here we only discuss the initialization of the atmospheric state on the MPAS mesh. The specification of that mesh is described in Chapter 2. The reference state used in the initialization is the analytic isothermal state given by {eq}`eq:8.1`. The hydrostatic balance uses the perturbation state balance {eq}`eq:8.2`. The surface pressure, and the potential temperature and water vapor mixing ratio in the column, are interpolated to the MPAS-Atmosphere mesh from the analysis. The lowest model level is interpolated from the analysis, and the lowest model level density is computed from the state equation. For each layer above the first layer, the hydrostatic relation {eq}`eq:8.2` along with the state equation are evaluated iteratively until the densities and pressures satisfy both the state equation and hydrostatic relation within a certain tolerance which in the current release requires the change in pressure in each iteration be less than 0.0001 Pa or that the iteration count is less than 30.

:::{admonition} MPAS code
:class: note

The initialization code for this test case is found in subroutine `init_atm_case_gfs` in the source code file `src/core_init_atmosphere/mpas_init_atm_cases.F`. Note that most of the subroutine, after the 3D MPAS mesh is created, involves the reading in and interpolating the analysis to the MPAS mesh.
:::

### 8.3.2 Lateral Boundary Conditions for Regional Simulations

The lateral boundary condition state for regional simulations is computed in the same manner as the full 3D initialization, except that the 3D mesh is taken from the full state initialization.
