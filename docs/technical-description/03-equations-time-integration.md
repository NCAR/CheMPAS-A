# Chapter 3: Equations and Time Integration

The MPAS-A dynamics solver integrates the compressible, nonhydrostatic Euler equations and a brief description of the solver is presented in Skamarock et al. (2012). The equations are cast in flux form using variables that have conservation properties, following the philosophy of Ooyama (1990). The equations are formulated using geometric height as an independent variable. The vertical coordinate is terrain following, using a hybrid $\sigma - z$ formulation that includes a filtering of the terrain effect for higher coordinate surfaces following Klemp (2011) as described in section 2.2.

In this chapter we present the moist flux-form equations cast in coordinates for the sphere, and we describe the time integration scheme used in the MPAS-A solver. The current MPAS-A release use equations employing the shallow-atmosphere approximation, and we focus on these equations in the next two chapters. The perturbation equations we employ in MPAS-A are formulated to facilitate the integration scheme we employ, particularly the acoustic splitting used in the MPAS-A time integration. This motivates our presentation of the equations and time integration together in this chapter.

## 3.1 Continuous Equations

MPAS uses the spherical centroidal Voronoi tesselation (mesh) (SCVT) for horizontal tiling of the sphere which is described in Chapter 2. Contrary to a latitude-longitude grid, this mesh is not based on a global coordinate hence equations are cast in general form. For the vertical coordinate, MPAS employs a height-based terrain-following vertical coordinate $\zeta$ following Klemp (2011) as described in section 2.2. The coupled prognostic variables in this coordinate are defined as

$$
\begin{aligned}
\tilde{\rho}_d &= \rho_d/\zeta_z,
\end{aligned}
$$ (eq:3.1)

$$
\begin{aligned}
\mathbf{V}_H &= \tilde{\rho}_d\,\mathbf{v}_H,
\end{aligned}
$$ (eq:3.2)

$$
\begin{aligned}
W &= \tilde{\rho}_d\,w,
\end{aligned}
$$ (eq:3.3)

$$
\begin{aligned}
\Theta_m &= \tilde{\rho}_d\,\theta\bigl[1 + (R_v/R_d)\,q_v\bigr],
\end{aligned}
$$ (eq:3.4)

$$
\begin{aligned}
Q_j &= \tilde{\rho}_d\,q_j,
\end{aligned}
$$ (eq:3.5)

where the vertical derivative of the computational coordinate $\zeta_z = \partial\zeta/\partial z$, $\rho_d$ is the dry-air density, $\mathbf{v}_H$ and $w$ represent the horizontal and vertical velocities $(u\mathbf{i}, v\mathbf{j}, w\mathbf{k})$ where $\mathbf{i} \times \mathbf{j} = \mathbf{k}$, $\theta$ is the potential temperature, $q_v$ is the water vapor mixing ratio with respect to the dry-air density $\rho_d$, $q_j$ are the mixing ratios of other moisture constituents and other scalars, and $R_d$ and $R_v$ are the gas constants for dry air and water vapor, respectively. The full continuous equations for both deep- and shallow-atmosphere configurations can be cast as

$$
\begin{aligned}
\frac{\partial \mathbf{V}_H}{\partial t} = {}&- \frac{\rho_d}{\rho_m}\!\left[\nabla_\zeta\!\left(\frac{p}{\zeta_z}\right) + \frac{\partial}{\partial\zeta}\!\left(\frac{\zeta_H\,p}{\zeta_z}\right)\right] - \eta\,\mathbf{k}\times\mathbf{V}_H - \tilde{\rho}_d\,\nabla_\zeta K \\
&- \mathbf{v}_H\nabla_\zeta\cdot\mathbf{V} - \frac{\partial\Omega\,\mathbf{v}_H}{\partial\zeta} + \mathbf{F}_{V_H} - \beta_d\,\tilde{\rho}_d\!\left(ew\mathbf{i} + \frac{\mathbf{v}_H\,w}{r}\right),
\end{aligned}
$$ (eq:3.6)

$$
\begin{aligned}
\frac{\partial W}{\partial t} = {}&-\frac{\rho_d}{\rho_m}\!\left[\frac{\partial p}{\partial\zeta} + g\!\left(\frac{r_o}{r}\right)^{\!2}\tilde{\rho}_m\right] - \bigl(\nabla\cdot\mathbf{V}\,w\bigr)_\zeta + F_W \\
&+ \beta_d\,\tilde{\rho}_d\!\left(\frac{\mathbf{v}_H\cdot\mathbf{v}_H}{r} + e\,\mathbf{v}_H\cdot\mathbf{i}\right),
\end{aligned}
$$ (eq:3.7)

$$
\frac{\partial \Theta_m}{\partial t} = -\bigl(\nabla\cdot\mathbf{V}\,\theta_m\bigr)_\zeta + F_{\Theta_m},
$$ (eq:3.8)

$$
\frac{\partial \tilde{\rho}_d}{\partial t} = -\bigl(\nabla\cdot\mathbf{V}\bigr)_\zeta,
$$ (eq:3.9)

$$
\frac{\partial Q_j}{\partial t} = -\bigl(\nabla\cdot\mathbf{V}\,q_j\bigr)_\zeta + F_{Q_j}.
$$ (eq:3.10)

The equation of state,

$$
p = p_o\!\left(\frac{R_d\,\Theta_m}{p_o}\right)^{\!c_p/c_v}
$$ (eq:3.11)

is used to diagnose pressure. The moist density $\rho_m = \rho_d(1 + \sum q_m)$, where $q_m$ are mixing ratios for water vapor and all hydrometeors. $\beta_d = 1$ for the deep-atmosphere equations, and the traditional shallow-atmosphere approximation is recovered when $\beta_d = 0$ and the radial distance from the center of the sphere $r$ is approximated by the constant sphere radius $r_o$. The gravity $g$ is defined by its value at the sphere surface. The horizontal momentum equation {eq}`eq:3.6` is not cast in conservative (flux) form. In {eq}`eq:3.6` the absolute vertical vorticity $\eta = \mathbf{n}\cdot\nabla\times\mathbf{v}_H + f$, where $\mathbf{n}$ is the unit vector normal to a $\zeta$ surface, and $K = |\mathbf{v}_H|^2/2$ is the horizontal kinetic energy. In the curvature and Coriolis terms in {eq}`eq:3.6` and {eq}`eq:3.7`, $f = 2\Omega_e\sin\psi$, $e = 2\Omega_e\cos\psi$, $\psi$ is the latitude, $\Omega_e$ is the angular rotation rate of the Earth. The terms $\mathbf{F}_{V_H}$, $F_W$, $F_{\Theta_m}$ and $F_{Q_j}$ represent sources and sinks from physics, sub-grid models and filters. The mass-flux normal to the horizontal coordinate surface is $\Omega = \mathbf{V}\cdot\nabla\zeta$ and the corresponding normal velocity is $\omega = \dot{\zeta}$. The flux divergence

$$
(\nabla\cdot\mathbf{V}\,\phi)_\zeta = \nabla_\zeta\cdot(\mathbf{V}_{\mathbf{H}}\,\phi) + \frac{\partial(\Omega\,\phi)}{\partial\zeta},
$$

where $\mathbf{V}_{\mathbf{H}}$ and $\mathbf{v}_{\mathbf{H}}$ are the horizontal momentum and velocity, respectively.

The current release of MPAS, Version 8, does not contain the solver for the deep atmosphere (full) equations, rather MPAS-Atmosphere employs the shallow atmosphere approximation in the release. A MPAS-Atmosphere solver for the full equations {eq}`eq:3.6`–{eq}`eq:3.10` has been developed and tested as reported in Skamarock et al. (2021a).

## 3.2 Time Integration Overview

MPAS integrates the equations of motion using a Runge-Kutta (RK3) time integration scheme, described in Wicker and Skamarock (2002), within which a split-explicit time integration technique is employed to integrate acoustic and gravity wave modes (Klemp et al. 2007). For the dry dynamics, prognostic equations {eq}`eq:3.6`–{eq}`eq:3.9` for the density-coupled horizontal velocity, vertical velocity, potential temperature $\Theta_m$, and the dry air density are integrated in the RK3 split-explicit dynamics time step. Other scalars, e.g. water vapor, precipitation species, number concentrations in microphysics, etc, are integrated in a separate Runge-Kutta step after the dynamics are advanced.

**[Figure 3.1: The MPAS time integration methodology. To be added next session.]**

The MPAS time step integration sequence is depicted in Figure 3.1. The main MPAS time step, `config_dt` in the `namelist.atmosphere`, is the time step for the scalar transport which is performed as the final process in the time integration, and this scalar integration appears as the final loop in the pseudo-code in Figure 3.1. The dynamics can be integrated with a smaller time step than the main MPAS time step, and this time step is specified through the `namelist.atmosphere` parameters as `config_dt/config_dynamics_split`, where `config_dynamics_split` is an integer. The loop for the dynamics steps is the first loop in Figure 3.1, and within it are the RK3 substeps within which appear a loop for the acoustic time steps that integrate each RK3 substep.

:::{admonition} MPAS code
:class: note

The time step integration is performed in subroutine `atm_srk3` in `MPAS/src/core_atmosphere/dynamics/MPAS_atm_time_integration.F`.
The dynamics RK3 integration is in the named loop `rk3_dynamics`. The number of dynamics steps taken each transport step is given by the `namelist.atmosphere` variable `config_dynamics_split` and the loop executing these steps is the named loop `dynamics_substeps`. The acoustic steps inside the dynamics integration are controlled by the named loop `dynamics_acoustic_steps`. The scalar RK3 transport is performed in in the named loop `rk3_split_transport`.
:::

In the following sections we will describe the RK3 time integration scheme and the perturbation equations it employs, followed by a description of the splitting of the acoustic and gravity-wave modes (the acoustic steps within the dynamics integration), a summary of the dynamics step, and a summary of the scalar transport time integration.

## 3.3 RK3 Time Integration Scheme

The RK3 time integration scheme, described in Wicker and Skamarock (2002), is comprised of three steps to advance a solution over a time step $\Delta t$. To advance a variable $\phi$ from time $t$ to time $t + \Delta t$ using RK3, where $\phi_t = L(\phi)$,

$$
\phi^{*} = \phi^{t} + \frac{\Delta t}{3}\,L(\phi^{t}),
$$ (eq:3.12)

$$
\phi^{**} = \phi^{t} + \frac{\Delta t}{2}\,L(\phi^{*}),
$$ (eq:3.13)

$$
\phi^{t+\Delta t} = \phi^{t} + \Delta t\,L(\phi^{**}).
$$ (eq:3.14)

The scheme {eq}`eq:3.12`–{eq}`eq:3.14` is third-order accurate in time for linear systems and second-order for nonlinear systems. A variant of this scheme available in MPAS replaces the coefficient $\Delta t/3$ in {eq}`eq:3.12` with $\Delta t/2$, and the resulting scheme is second-order accurate for linear and nonlinear systems. The default MPAS configuration uses the second order variant because it allows for a larger time step as shown in linear analysis results given in Figure 3.2.

:::{admonition} MPAS code
:class: note

The RK3 time-step coefficients are set in two places in the main MPAS time integration driver in `MPAS/src/core_atmosphere/dynamics/MPAS_atm_time_integration.F`, in subroutine `atm_srk3`. The first location is before the main Runge-Kutta loop for integrating the dry dynamics, and the coefficients are used for the named loop `rk3_dynamics`. The second location is immediately before the named loop `rk3_split_transport` where, as the name implies, an RK3 integration of the scalar transport is accomplished.
:::

**[Figure 3.2: Runge-Kutta time integration response. RK3 unstable for $k\Delta t > 1.73$; RK3/2 unstable for $k\Delta t > 2$. To be added next session.]**

### 3.3.1 Perturbation Equations Employed in the RK3 Solver

The shallow atmosphere approximation is recovered by setting $\beta_d = 0$ and $r = r_o$ in {eq}`eq:3.6`–{eq}`eq:3.10`. We describe the discrete solver for the shallow atmosphere equations in the following chapters, so for clarity and ease of reference we will describe the shallow-atmosphere equations and solver in the rest of this note.

Following Klemp et al. (2007), two sets of perturbations are used in the MPAS-A solver. Perturbation variables are defined about a hydrostatically-balanced thermodynamic state that are only a function of the geometric height $z$: $\rho_d = \overline{\rho}_d(z) + \rho'_d$, and $\Theta_m = \overline{\rho}_d(z)\overline{\theta}_m(z) + \Theta'_m$. The reference state $(\overline{\rho}_d$ and $\overline{\Theta}_m)$ satisfies the equation of state {eq}`eq:3.11` that defines the reference pressure $\overline{p}$ and hence a perturbation pressure $p' = p - \overline{p}(z)$. The thermodynamic reference state is for a dry atmosphere $(\overline{q_m} = 0)$. Incorporating these perturbation variables into {eq}`eq:3.6`–{eq}`eq:3.10`, making the shallow atmosphere approximation, and temporally discretizing the equations {eq}`eq:3.6`–{eq}`eq:3.10` yields

$$
\begin{aligned}
\delta_t \mathbf{V}_H = {}&\left[-\frac{\rho_d}{\rho_m}\!\left(\nabla_\zeta\!\left(\frac{p'}{\zeta_z}\right) + g\,z_H\tilde{\rho}'_m\right) - \eta\,\mathbf{k}\times\mathbf{V}_H - \tilde{\rho}_d\,\nabla_\zeta K \right. \\
&\left. - \mathbf{v}_H\nabla_\zeta\cdot\mathbf{V} - \frac{\partial\Omega\,\mathbf{v}_H}{\partial\zeta}\right]_{t_1} + \left[\mathbf{F}_{V_H}\right]_{t_2},
\end{aligned}
$$ (eq:3.15)

$$
\delta_t W = \left[-\frac{\rho_d}{\rho_m}\!\left(\frac{\partial p'}{\partial\zeta} + g\tilde{\rho}'_m\right) - \bigl(\nabla\cdot\mathbf{V}\,w\bigr)_\zeta\right]_{t_1} + \left[F_W\right]_{t_2},
$$ (eq:3.16)

$$
\delta_t \Theta'_m = \left[-\bigl(\nabla\cdot\mathbf{V}\,\theta_m\bigr)_\zeta\right]_{t_1} + \left[F_{\Theta_m}\right]_{t_2},
$$ (eq:3.17)

$$
\delta_t \tilde{\rho}'_d = \left[-\bigl(\nabla\cdot\mathbf{V}\bigr)_\zeta\right]_{t_1},
$$ (eq:3.18)

$$
\delta_t Q_j = \left[-\bigl(\nabla\cdot\mathbf{V}\,q_j\bigr)_\zeta\right]_{t_1} + \left[F_{Q_j}\right]_{t_2},
$$ (eq:3.19)

where the discrete time difference over the timestep $\Delta t$ for an arbitrary variable $\phi$ is

$$
\delta_t\phi = \frac{\phi^{t+\Delta t} - \phi^{t}}{\Delta t},
$$ (eq:3.20)

the bracketed terms $[\ ]_{t_*}$ denote terms evaluated at time $t_*$, and

$$
\tilde{\rho}'_m = \tilde{\rho}'_d + \tilde{\rho}_d\sum q_j.
$$ (eq:3.21)

We have also recast the vertical pressure gradient correction in the horizontal momentum equation by its hydrostatic approximation:

$$
\frac{\partial p'}{\partial\zeta} = -\,g\,\tilde{\rho}'_m,
$$ (eq:3.22)

which has been found to improve the accuracy in computing the horizontal pressure gradient by reducing truncation errors in the numerics (Klemp and Skamarock 2021).

### 3.3.2 Operators in the RK3 Scheme

The operator $L(\phi)$ in {eq}`eq:3.12`–{eq}`eq:3.14` includes all the terms except for the partial derivatives with respect to time in the prognostics equations for the dynamics {eq}`eq:3.15`–{eq}`eq:3.18`, i.e. all the terms to the right-hand-side of the equal signs. However, the dynamics integration using $L(\phi)$ in {eq}`eq:3.12`–{eq}`eq:3.14` involves evaluating only some of the RHS terms in {eq}`eq:3.15`–{eq}`eq:3.18` in each of the RK3 substeps. Specifically, the bracketed terms $[\ ]_{t_1}$ in {eq}`eq:3.15`–{eq}`eq:3.18` are the terms evaluated from the states $\phi^t$, $\phi^*$, and $\phi^{**}$. The other RHS terms in {eq}`eq:3.15`–{eq}`eq:3.18`, given as bracketed terms $[\ ]_{t_2}$, are evaluated at prescribed times during the integration. These other terms include atmospheric physics, whose tendencies may be updated either once per time step or even less frequently (see Chapter 9), and subgrid turbulence models and MPAS model filters, that usually are evaluated at the beginning of a time step at time $t$ and held fixed over the time step. Additionally, some filters are evaluated on the acoustic timestep.

:::{admonition} MPAS code
:class: note

The RK3 tendencies for the dynamics, the bracketed terms $[\ ]_{t_1}$ in {eq}`eq:3.15`–{eq}`eq:3.18`, are computed in a call to subroutine `atm_compute_dyn_tend` from the main MPAS time integration routine `atm_srk3` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
For the $[\ ]_{t_2}$ terms in {eq}`eq:3.15`–{eq}`eq:3.18`, the terms associated with MPAS model filters and subgrid transport, aside from the column physics, are computed in subroutine `atm_compute_dyn_tend` during the first RK3 substep and saved. The tendencies from the column physics are acquired in the call to `subroutine physics_get_tend` in routine `atm_srk3` at the beginning of the main MPAS timestep. The dynamics tendencies, the saved physics tendencies and saved model filter and subgrid transport tendencies are summed in `subroutine atm_compute_dyn_tend` to construct $L(\phi)$ in each RK3 substep {eq}`eq:3.12`–{eq}`eq:3.14`.
:::

## 3.4 Acoustic Solver

Within the RK3 substeps, smaller acoustic steps are taken where only the terms responsible for acoustic modes are integrated while the terms responsible for the slower atmospheric modes are held fixed to their RK3 evaluation. Thus the Runge-Kutta time integration takes place within the acoustic steps where the precomputed RK3 tendencies are added to the prognostic variables. The combined timesteps are illustrated in Figure 3.3 which includes the graphical relation of the acoustic timestep to the RK3 timestep in each RK3 substep. The RK3 substeps $(\Delta t/3, \Delta t/2, \Delta t)$ are accomplished over the acoustic steps within them. In the following section we present the perturbation equations for the acoustic steps along with gravity wave and acoustic filters that are employed in MPAS-A.

**[Figure 3.3: Integration timesteps in the split-explicit RK3 scheme. To be added next session.]**

### 3.4.1 Perturbation Equations Employed in the Acoustic Solver

The fully compressible equations contain acoustic modes that can severely limit the allowable timesteps in an explicit numerical time integration. MPAS employs a split-explicit vertically-implicit time integration scheme to advance the acoustic modes with smaller timesteps. This scheme is described in detail in Klemp et al. (2007). The acoustic step advances the dynamics variables defined as perturbations from their respective values at the beginning of the acoustic time steps at time $t$:

$$
\begin{aligned}
\mathbf{V}''_{\mathbf{h}} &= \mathbf{V}_{\mathbf{h}} - \mathbf{V}_{\mathbf{h}}^{t} \\
W'' &= W - W^{t} \\
\Theta''_{m} &= \Theta'_{m} - \Theta_{m}^{\prime\,t} \\
\tilde{\rho}''_d &= \tilde{\rho}'_d - \tilde{\rho}_d^{\prime\,t}
\end{aligned}
$$ (eq:3.23)

Additionally, we convert the vertical momentum equation for the coupled contravariant vertical velocity $W''$ to an equation for the covariant vertical velocity $\Omega = \Omega^t + \Omega''$ (i.e. the coupled velocity normal to the hybrid coordinate horizontal surfaces), defined as

$$
\Omega = V_h\cdot\nabla\zeta + \zeta_z W,
$$ (eq:3.24)

so that we can use the upper and lower boundary conditions $\Omega = \Omega'' = 0$ at the top and bottom of the domain during the acoustic step integration. The acoustic integration follows Klemp et al. (2007) and employs a forward-backward time integration of the system, first advancing the horizontal momentum $\mathbf{V}_{\mathbf{h}}''$ followed by advancing $\Omega'' = V_h''\cdot\nabla\zeta + \zeta_z W''$, $\Theta''_m$ and $\tilde{\rho}''_d$ with semi-implicit numerics, using the updated values of $\mathbf{V}_{\mathbf{h}}''$. In terms of these double-prime perturbation variables, the acoustic step equations are very similar to equations (13)–(16) in Klemp et al. (2007), except we substitute the equation for $W''$ with that for $\Omega''$ and we employ a slightly different linearization as discussed below.

$$
\delta_\tau \mathbf{V}''_{\mathbf{h}} = -\,\frac{\rho_d^t}{\rho_m^t}\!\left[\gamma\,R_d\,\pi^t\,\nabla_\zeta\overline{\Theta''_m}^\tau + g\,z_H\,\tilde{\rho}''_d\right] + \mathbf{R}^t_{V_h}
$$ (eq:3.25)

$$
\begin{aligned}
\delta_\tau \Omega'' + \zeta_z\frac{\rho_d^t}{\rho_m^t}&\!\left[\gamma\,R_d\,\pi^t\,\partial_\zeta\!\left(\zeta_z\overline{\Theta''_m}^\tau\right) - g\,\tilde{\rho}_m^t\frac{R_d}{c_v}\frac{\overline{\Theta''_m}^\tau}{\Theta^t_m}\right] \\
&+ \zeta_z\,g\,\overline{\tilde{\rho}''_d}^\tau = \mathbf{R}^t_{V_H}\cdot\nabla_h\zeta + R^t_W\,\zeta_z
\end{aligned}
$$ (eq:3.26)

$$
\delta_\tau \Theta''_m + \nabla_\zeta\cdot\!\left(\mathbf{V}''^{\,\tau+\Delta\tau}_{\mathbf{h}}\,\theta^t_m\right) + \partial_\zeta\!\left(\overline{\Omega''^\tau}\,\theta^t_m\right) = R^t_{\Theta_m}
$$ (eq:3.27)

$$
\delta_\tau \tilde{\rho}''_d + \nabla_\zeta\cdot\mathbf{V}''^{\,\tau+\Delta\tau}_{\mathbf{h}} + \partial_\zeta\!\left(\overline{\Omega''^\tau}\right) = R^t_{\tilde{\rho}_d}
$$ (eq:3.28)

The operator

$$
\overline{\phi}^{\,\tau} = \frac{1+\epsilon}{2}\,\phi^{\tau+\Delta\tau} + \frac{1-\epsilon}{2}\,\phi^{\tau}
$$ (eq:3.29)

is a time averaging of $\phi$ and represents the vertically semi-implicit time differencing of the vertically-propagating acoustic modes and buoyancy oscillation modes in the compressible non-hydrostatic system. For the system {eq}`eq:3.28`–{eq}`eq:3.26` the vertically semi-implicit scheme will be second-order accurate and neutral (non-damping) for $\epsilon = 0$ and will be damping for $\epsilon > 0$. This damping mechanism is discussed further in section 5.4. The terms $\mathbf{R}^t_{V_H}$, $R^t_W$, $R^t_{\Theta_m}$ and $R^t_{\tilde{\rho}_d}$ in {eq}`eq:3.25`–{eq}`eq:3.28` are the RHS terms in {eq}`eq:3.15`–{eq}`eq:3.18`, that is

$$
\begin{aligned}
\mathbf{R}^t_{V_H} = {}&\left[-\frac{\rho_d}{\rho_m}\!\left(\nabla_\zeta\!\left(\frac{p'}{\zeta_z}\right) + g\,z_H\tilde{\rho}'_m\right) - \eta\,\mathbf{k}\times\mathbf{V}_H - \tilde{\rho}_d\,\nabla_\zeta K \right. \\
&\left. - \mathbf{v}_H\nabla_\zeta\cdot\mathbf{V} - \frac{\partial\Omega\,\mathbf{v}_H}{\partial\zeta}\right]_{t_1} + \left[\mathbf{F}_{V_H}\right]_{t_2},
\end{aligned}
$$ (eq:3.30)

$$
R^t_W = \left[-\frac{\rho_d}{\rho_m}\!\left(\frac{\partial p'}{\partial\zeta} + g\tilde{\rho}'_m\right) - \bigl(\nabla\cdot\mathbf{V}\,w\bigr)_\zeta\right]_{t_1} + \left[F_W\right]_{t_2},
$$ (eq:3.31)

$$
R^t_{\Theta_m} = \left[-\bigl(\nabla\cdot\mathbf{V}\,\theta_m\bigr)_\zeta\right]_{t_1} + \left[F_{\Theta_m}\right]_{t_2},
$$ (eq:3.32)

$$
R^t_{\tilde{\rho}_d} = \left[-\bigl(\nabla\cdot\mathbf{V}\bigr)_\zeta\right]_{t_1}.
$$ (eq:3.33)

**Derivation of the Acoustic Step Momentum Equations**

These acoustic-step equations are readily obtained (for the most part) by substituting the perturbation form of the variables {eq}`eq:3.23` into {eq}`eq:3.15`–{eq}`eq:3.18` and collecting all of the terms with double-prime variables on the left hand side of the equations {eq}`eq:3.25`–{eq}`eq:3.28` and all of the terms at time $t$ on the right hand side {eq}`eq:3.30`–{eq}`eq:3.33`. In these equations we have recast the horizontal and vertical pressure gradients in terms of the prognostic variables in order to gain computational efficiency in the acoustic integration where we treat terms responsible for the vertical propagation of acoustic modes in an implicit manner. These acoustic equations are constructed beginning from the RK3 equations where the pressure gradients in the momentum equations {eq}`eq:3.15` and {eq}`eq:3.16` are modified such that they are cast in terms of the prognostic variable $\Theta_m$ using the equation of state {eq}`eq:3.11`;

$$
\nabla p = \gamma\,R_d\,\pi\,\nabla\Theta_m,
$$ (eq:3.34)

where $\pi = (p/p_o)^{R_d/c_p}$ is the Exner function and $\gamma = c_p/c_v$ is the ratio of the specific heats. Using this representation, the perturbation form of vertical pressure gradient and buoyancy terms in the vertical momentum equation {eq}`eq:3.16` can be written as

$$
\begin{aligned}
\frac{\partial p'}{\partial\zeta} + g\,\tilde{\rho}'_m
&= \gamma R_d(\overline{\pi}+\pi')\,\frac{\partial}{\partial\zeta}\!\left[\zeta_z(\overline{\Theta}_m + \Theta'_m)\right] + g(\overline{\tilde{\rho}}_d + \tilde{\rho}'_m) \\
&= \gamma R_d\!\left[\pi\,\frac{\partial}{\partial\zeta}(\zeta_z\Theta'_m) + \pi'\,\frac{\partial}{\partial\zeta}(\zeta_z\overline{\Theta}_m)\right] + g\,\tilde{\rho}'_m \\
&= \gamma R_d\,\pi\,\frac{\partial}{\partial\zeta}(\zeta_z\Theta'_m) - g\,\overline{\tilde{\rho}}_d\frac{\pi'}{\overline{\pi}} + g\,\tilde{\rho}'_m.
\end{aligned}
$$ (eq:3.35)

The pressure gradients {eq}`eq:3.34` and {eq}`eq:3.35` correspond those used in the horizontal and vertical momentum equations (11) and (12) in Klemp et al. (2007).

Substitution of the acoustic perturbation variables {eq}`eq:3.23` into the modified RK3 equations, with some additional approximations, leads to the acoustic equations {eq}`eq:3.25`–{eq}`eq:3.28`. For the horizontal pressure gradient cast in terms of $\Theta_m$ {eq}`eq:3.34`, $\pi$ is linearized about its value at time $t$ in {eq}`eq:3.25` for the integration over the acoustic steps. The representation of the vertical pressure gradient in the vertical momentum equation merits some further explanation [the terms within the large square brackets in {eq}`eq:3.26`]. It is obtained by substituting the perturbation representation of $\Theta_m$ and $\pi$ into {eq}`eq:3.35`:

$$
\begin{aligned}
\frac{\partial p'}{\partial\zeta}
&= \gamma R_d(\pi^t + \pi'')\,\frac{\partial}{\partial\zeta}\!\left[\zeta_z(\Theta_m^{\prime\,t} + \Theta''_m)\right] - g\,\tilde{\overline{\rho}}_d\frac{\pi^{\prime\,t} + \pi''}{\overline{\pi}} \\
&= \frac{\partial p'^{\,t}}{\partial\zeta} + \gamma R_d\,\pi^t\,\frac{\partial}{\partial\zeta}(\zeta_z\Theta''_m) + \left[\gamma R_d\,\pi^t\,\frac{\partial}{\partial\zeta}(\zeta_z\Theta_m^t) - g\,\tilde{\overline{\rho}}_d\frac{\pi^t}{\overline{\pi}}\right]\frac{\pi''}{\pi^t} \\
&= \frac{\partial p'^{\,t}}{\partial\zeta} + \gamma R_d\,\pi^t\,\frac{\partial}{\partial\zeta}(\zeta_z\Theta''_m) + \left[\frac{\partial p'^{\,t}}{\partial\zeta} + g\,\tilde{\overline{\rho}}_d\frac{\pi^t}{\overline{\pi}} - g\,\tilde{\overline{\rho}}_d\frac{\pi^t}{\overline{\pi}}\right]\frac{R_d}{c_v}\frac{\Theta''_m}{\Theta^t_m} \\
&= \frac{\partial p'^{\,t}}{\partial\zeta} + \gamma R_d\,\pi^t\,\frac{\partial}{\partial\zeta}(\zeta_z\Theta''_m) - g\,\tilde{\rho}_m^t\,\frac{R_d}{c_v}\,\frac{\Theta''_m}{\Theta^t_m}.
\end{aligned}
$$ (eq:3.36)

In deriving this expression, we have employed an approximation using a perturbation form of the equation of state $\pi''/\pi^t \simeq (R_d/c_v)(\Theta''_m/\Theta^t_m)$. Also, in simplifying the coefficient of the $\Theta''_m/\Theta^t_m$ term we have used a hydrostatic approximation for the vertical pressure gradient $p'^{\,t}_\zeta \simeq -g\tilde{\overline{\rho}}^t_m$, as we did in {eq}`eq:3.22`. With this representation for $\partial p'/\partial\zeta$, the coefficient of the $\Theta''_m/\Theta^t_m$ term in {eq}`eq:3.26` differs slightly from that in equation (14) in Klemp et al. (2007), and represents a somewhat more consistent derivation for the double-prime perturbation variables. We also note that in their current form the coefficients of all of the double-prime terms in {eq}`eq:3.25`–{eq}`eq:3.28` are exactly the same as in the corresponding acoustic equations obtained directly from {eq}`eq:3.6`–{eq}`eq:3.9` using the full variables for $p$, $\tilde{\rho}_d$, and $\Theta_m$ without expressing them as perturbations from a specified reference sounding profile.

### 3.4.2 Acoustic Timestep

We define an acoustic timestep as advancing $\mathbf{V}_{\mathbf{h}}''$, $\Omega''$, $\Theta''_m$, and $\tilde{\rho}''_d$ from time $\tau$ to time $\tau + \Delta\tau$, where $\Delta\tau$ is the acoustic timestep. The acoustic timestep begins with advancing the horizontal momentum equation {eq}`eq:3.25` using the known values of $\Theta''_m$ and $\tilde{\rho}''_d$ at $\tau$. Having the values $\mathbf{V}_{\mathbf{h}}''$ at time $\tau+\Delta\tau$, and using the time-weighting operator {eq}`eq:3.29` and the time derivative operator {eq}`eq:3.20`, we can rearrange the acoustic update equations for $\Omega''$, $\Theta''_m$ and $\tilde{\rho}''_d$, {eq}`eq:3.26`, {eq}`eq:3.27` and {eq}`eq:3.28`, as

$$
\begin{aligned}
\Omega''^{\,\tau+\Delta\tau} + \Delta\tau\,\frac{1+\epsilon}{2}&\!\left\{\zeta_z\,\frac{\rho_d^t}{\rho_m^t}\!\left[\gamma R_d\pi^t\partial_\zeta\!\bigl(\zeta_z\Theta''^{\,\tau+\Delta\tau}_m\bigr) - g\,\tilde{\rho}_m\frac{R_d}{c_v}\frac{\Theta''^{\,\tau+\Delta\tau}_m}{\Theta^t_m}\right]\right. \\
&\left.+ \zeta_z\,g\,\tilde{\rho}''^{\,\tau+\Delta\tau}_d\right\} = F_{\Omega^\tau}
\end{aligned}
$$ (eq:3.37)

$$
\Theta''^{\,\tau+\Delta\tau}_m = -\Delta\tau\,\frac{1+\epsilon}{2}\,\partial_\zeta\!\bigl(\Omega''^{\,\tau+\Delta\tau}\theta^t_m\bigr) + F_{\Theta^\tau_m}
$$ (eq:3.38)

$$
\tilde{\rho}''^{\,\tau+\Delta\tau}_d = -\Delta\tau\,\frac{1+\epsilon}{2}\,\partial_\zeta\!\bigl(\Omega''^{\,\tau+\Delta\tau}\bigr) + F_{\tilde{\rho}^\tau}
$$ (eq:3.39)

where the right-hand-side terms $F_{\Omega^\tau}$, $F_{\Theta^\tau_m}$, and $F_{\tilde{\rho}^\tau}$ are known:

$$
\begin{aligned}
F_{\Omega^\tau} = {}&\Omega''^{\,\tau} - \Delta\tau\,\frac{1-\epsilon}{2}\!\left\{\zeta_z\,\frac{\rho_d^t}{\rho_m^t}\!\left[\gamma R_d\pi^t\partial_\zeta\!\bigl(\zeta_z\Theta''^{\,\tau}_m\bigr) - g\,\tilde{\rho}_m\frac{R_d}{c_v}\frac{\Theta''^{\,\tau}_m}{\Theta^t_m}\right]\right. \\
&\left.\hphantom{{}-\Delta\tau\,\frac{1-\epsilon}{2}\!} + \zeta_z\,g\,\tilde{\rho}''^{\,\tau}_d\right\} + \Delta\tau\!\left[\delta_\tau\mathbf{V}''^{\,\tau+\Delta\tau}_{\mathbf{h}}\cdot\nabla_h\zeta + \zeta_z R^t_W\right]
\end{aligned}
$$ (eq:3.40)

$$
F_{\Theta^\tau_m} = \Theta''^{\,\tau}_m + \Delta\tau\!\left\{-\nabla_\zeta\cdot\!\left(\mathbf{V}''^{\,\tau+\Delta\tau}_{\mathbf{h}}\,\theta^t_m\right) - \frac{1-\epsilon}{2}\,\partial_\zeta\bigl(\Omega''^{\,\tau}\theta^t_m\bigr) + R^t_{\Theta_m}\right\}
$$ (eq:3.41)

$$
F_{\tilde{\rho}^\tau} = \tilde{\rho}''^{\,\tau}_d + \Delta\tau\!\left\{-\nabla_\zeta\cdot\mathbf{V}''^{\,\tau+\Delta\tau}_{\mathbf{h}} - \frac{1-\epsilon}{2}\,\partial_\zeta\bigl(\Omega''^{\,\tau}\bigr) + R^t_{\tilde{\rho}_d}\right\}
$$ (eq:3.42)

**Implicit Acoustic Solver**

Equations {eq}`eq:3.37`–{eq}`eq:3.39` form a coupled system for $(\Omega'', \Theta''_m, \rho''_d)$ at time $\tau + \Delta\tau$. The solution of this system is obtained using {eq}`eq:3.38` and {eq}`eq:3.39` to eliminate $\Theta''^{\,\tau+\Delta\tau}_m$ and $\tilde{\rho}''^{\,\tau+\Delta\tau}_d$ from {eq}`eq:3.37`, resulting in an equation for the single unknown variable $\Omega''^{\,\tau+\Delta\tau}$. The variables $\Theta''_m$ and $\tilde{\rho}''_d$ are vertically staggered relative to $\Omega''$, and interpolation of these variables, along with vertical derivatives, results in an equation for $\Omega''^{\,\tau+\Delta\tau}$ where the new value at a level depends on the new values above and below it thus requiring us to solve a tridiagonal matrix. A description of the formulation of the tridiagonal system and solver in the acoustic step is given in Appendix A. Subsequent to the tridiagonal solution for $\Omega''^{\,\tau+\Delta\tau}$, vertical velocity damping is applied to damp vertically-propagating gravity waves (see section 3.4.3), after which $\Theta''^{\,\tau+\Delta\tau}_m$ and $\tilde{\rho}''^{\,\tau+\Delta\tau}_d$ are computed by back substitution of $\Omega''^{\,\tau+\Delta\tau}$ into {eq}`eq:3.38` and {eq}`eq:3.39`.

:::{admonition} MPAS code
:class: note

The acoustic timestep is accomplished in subroutine `atm_advance_acoustic_step`. The tridiagonal solver coefficients are computed in subroutine `atm_compute_vert_imp_coefs`. The tendency for $\Omega$ involving $R_w$ and $R_{V_H}$, the right-hand-side of {eq}`eq:3.26`, is computed in subroutine `atm_set_smlstep_pert_variables`, while other pieces involving multiplication by $\zeta_z$ are picked up directly in subroutine `atm_compute_dyn_tend`. All of these routines can be found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::

### 3.4.3 Filters in the Acoustic Timestep

There are two filters that are applied on each acoustic timestep - a gravity wave filter following Klemp et al. (2008) that damps the vertical velocity and serves as an absorber to prevent wave reflection at the rigid MPAS model top, and an acoustic mode filter as described by (Klemp et al. 2018).

**Gravity-Wave Absorbing Layer**

The upper boundary of MPAS is a rigid lid ($\Omega = 0$) through which there is no flow. Vertically-propagating gravity waves can reflect off the rigid upper boundary of MPAS. In order to prevent spurious gravity-wave reflection that can contaminate the solution as the reflected waves propagate downward, MPAS employs the filter described in Klemp et al. (2008) where the vertically velocity is damped in an implicit manner as part of each acoustic step.

The implicit damping term is included after solving the tridiagonal equation for $\Omega''$ {eq}`eq:3.37`. Denoting the solution to {eq}`eq:3.37` as $\Omega^{\tau *}$, the final solution for $\Omega''^{\,\tau+\Delta\tau}$ is computed by adding the implicit damping term for $\Omega$, recognizing that the $\Omega^\tau = \Omega^t + \Omega''^{\,\tau}$:

$$
\Omega''^{\,\tau+\Delta\tau} = \Omega''^{\,\tau *} - \Delta\tau\,R_\Omega\!\left[\Omega^t + \Omega''^{\,\tau+\Delta\tau}\right].
$$ (eq:3.43)

{eq}`eq:3.43` is equivalent to equation (9) in Klemp et al. (2008) except cast in terms of $\Omega''$ instead of $w$. $R_\Omega$ is the damping coefficient. It typically has a value of zero (no damping) from the surface through most of the atmosphere and uses the square of a sine function to increase from 0 to a specified positive value over a user-specified range below the model top. {eq}`eq:3.43` is trivially solved:

$$
\Omega''^{\,\tau+\Delta\tau} = \frac{\Omega''^{\,\tau *} - \Delta\tau\,R_\Omega\,\Omega^t}{1 + \Delta\tau\,R_\Omega}.
$$ (eq:3.44)

The damping coefficient $R_\Omega$ is specified as

$$
R_\Omega = \nu\,\sin^2\!\left[\frac{\pi}{2}\frac{(z-z_d)}{(z_t-z_d)}\right] \quad \text{for } z_d \le z \le z_t,\ \text{otherwise } R_\Omega = 0,
$$

where $z$ is the height of an $\Omega$ point, $z_d$ is the height above which the absorbing layer is active, $z_t$ is the height of the model top, and $\nu$ is the damping rate and has units $\mathrm{s}^{-1}$. After the velocity $\Omega''^{\,\tau+\Delta\tau}$ is computed using {eq}`eq:3.44`, $\Theta''^{\,\tau+\Delta\tau}_m$ and $\tilde{\rho}''^{\,\tau+\Delta\tau}_d$ are computed by back substitution of $\Omega''^{\,\tau+\Delta\tau}$ into {eq}`eq:3.38` and {eq}`eq:3.39`.

:::{admonition} MPAS code
:class: note

The gravity wave filter application takes place in subroutine `atm_advance_acoustic_step` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The coefficient $R_\Omega$ is set in subroutine `atm_compute_damping_coefs` in in `MPAS/src/core_atmosphere/mpas_atm_core.F` and the values are stored in the array variable named `dss`. The maximum value of the damping coefficient $R_\Omega$ is given by the `namelist.atmosphere` configuration variable `config_xnutr` and the filtering begins at `namelist.atmosphere` configuration variable height `config_zd` (m).
:::

**Acoustic Filter**

The acoustic filter in MPAS follows (Klemp et al. 2018) and damps 3D divergence to filter acoustic waves. (Klemp et al. 2018) propose that the appropriate divergence to damp is $\nabla\cdot(\mathbf{V}\Theta_m)$ where $\mathbf{V} = (\mathbf{V}_H, \Omega)$. In MPAS, we use the time-rate-of-change of $\Theta''_m$, given in {eq}`eq:3.27`, to define the divergence which includes diabatic heating. After {eq}`eq:3.25`–{eq}`eq:3.28` are advanced along with the gravity-wave filter, the acoustic filter is applied as an additional update to the horizontal momentum $\mathbf{V}_{\mathbf{H}}''^{\,*}$:

$$
\mathbf{V}_{\mathbf{H}}''^{\,\tau+\Delta\tau} = \mathbf{V}_{\mathbf{H}}''^{\,*} + \frac{\gamma_d\,\Delta x}{\Theta^t_m}\,\delta_\tau\Theta''_m,
$$ (eq:3.45)

where $\delta_\tau\Theta''_m = (\Theta''^{\,\tau+\Delta\tau}_m - \Theta''^{\,\tau}_m)/\Delta\tau$ as computed in the acoustic step. The damping coefficient $\gamma_d$ has a default value of 0.1. See (Klemp et al. 2018) for an analysis of this approach, specifically the discussion following their equations (21)–(26).

:::{admonition} MPAS code
:class: note

The acoustic filter application takes place in subroutine `atm_divergence_damping_3d` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`. The damping coefficient is `config_smdiv` in `namelist.atmosphere`.
:::

## 3.5 Scalar Transport Time Integration

The time integration of the scalar transport occurs after the dynamics steps, and the scalars are held fixed during the dynamics steps. As noted earlier, the dynamics timestep may be smaller than the scalar timestep, with the constraint that the scalar timestep be an integer multiple of the dynamics timestep. This is enforced in the specification of the dynamics timestep as being the scalar timestep divided by an integer.

Scalar transport is integrated using the RK3 time integration described in section 3.3 where the flux divergence and physics forcing terms are evaluated as the RK3 operator $L(\phi)$. Scalar transport (the flux divergence term in {eq}`eq:3.19`) is integrated with a shape-preserving (monotonic) scheme that employs limiters on the final RK3 substep to enforce shape preservation. The limiters do not take into account the forcing terms $F_{Q_j}$ in {eq}`eq:3.19`, so we advance the solution using the RK3 scheme as follows:

$$
Q_j^{*} = Q_j^{t} + \frac{\Delta t}{3}\,L(Q_j^{t}); \qquad L(Q_j^{t}) = \bigl[-\nabla\cdot\mathbf{V}\,q_j^{t}\bigr] + F_{Q_j},
$$ (eq:3.46)

$$
Q_j^{**} = Q_j^{t} + \frac{\Delta t}{2}\,L(Q_j^{*}); \qquad L(Q_j^{*}) = \bigl[-\nabla\cdot\mathbf{V}\,q_j^{*}\bigr] + F_{Q_j},
$$ (eq:3.47)

$$
Q_j^{***} = Q_j^{t} + \Delta t\,L(Q_j^{**}); \qquad L(Q_j^{**}) = F_{Q_j},
$$ (eq:3.48)

$$
Q_j^{t+\Delta t} = Q_j^{***} + \Delta t\,L(Q_j^{***}); \qquad L(Q_j^{***}) = \bigl[-\nabla\cdot\mathbf{V}\,q_j^{***}\bigr].
$$ (eq:3.49)

To accommodate the limiters, on the third RK3 substep we first update the scalars with the physics tendency, and then we apply the final transport step starting from the physics-updated values but using the $Q_j^{**}$ values in the flux divergence. If the physics forcing $F_{Q_j} = 0$, this scheme defaults to the {eq}`eq:3.12`–{eq}`eq:3.14` where the final step uses an operator that includes the monotonicity constraint of the flux divergence for scalars and the monotonicity constraints are described in Chapter 4.

:::{admonition} MPAS code
:class: note

There are two routines that perform scalar advection in MPAS. Substeps {eq}`eq:3.46` and {eq}`eq:3.47` are performed in subroutine `atm_advance_scalars`, and substeps {eq}`eq:3.48` and {eq}`eq:3.49` are performed in a single call to subroutine `atm_advance_scalars_mono`. These routines are found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
:::
