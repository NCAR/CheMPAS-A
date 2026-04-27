# Appendix A: Vertically-Implicit Acoustic Solution

The solution to the coupled system of equations {eq}`eq:3.37`, {eq}`eq:3.38` and {eq}`eq:3.39`, repeated here for convenience,

$$
\begin{aligned}
\Omega''^{\,\tau+\Delta\tau} + \Delta\tau\,\frac{1+\epsilon}{2}&\!\left\{\zeta_z\,\frac{\rho_d^t}{\rho_m^t}\!\left[\gamma R_d\pi^t\partial_\zeta\!\bigl(\zeta_z\Theta''^{\,\tau+\Delta\tau}_m\bigr) - g\,\tilde{\rho}_m\frac{R_d}{c_v}\frac{\Theta''^{\,\tau+\Delta\tau}_m}{\Theta^t_m}\right]\right. \\
&\left.+ \zeta_z\,g\,\tilde{\rho}''^{\,\tau+\Delta\tau}_d\right\} = \Delta\tau\,F_{\Omega^\tau},
\end{aligned}
$$

$$
\Theta''^{\,\tau+\Delta\tau}_m = -\Delta\tau\,\frac{1+\epsilon}{2}\,\partial_\zeta\!\bigl(\Omega''^{\,\tau+\Delta\tau}\theta^t_m\bigr) + \Delta\tau\,F_{\Theta^\tau_m},
$$

$$
\tilde{\rho}''^{\,\tau+\Delta\tau}_d = -\Delta\tau\,\frac{1+\epsilon}{2}\,\partial_\zeta\!\bigl(\Omega''^{\,\tau+\Delta\tau}\bigr) + \Delta\tau\,F_{\tilde{\rho}^\tau},
$$

is needed to complete the acoustic timestep. Given vertical C-grid staggering, the elimination of $\Theta''^{\,\tau+\Delta\tau}_m$ and $\tilde{\rho}''^{\,\tau+\Delta\tau}_d$ from {eq}`eq:3.37` using {eq}`eq:3.38` and {eq}`eq:3.39` results in a tridiagonal matrix that must be inverted. In this appendix we describe the setup of the tridiagonal coefficients that define the matrix needed to solve for $\Omega''^{\,\tau+\Delta\tau}$.

The C-grid vertical staggering of the MPAS mesh is depicted in chapter 2 section 2.2 figure 2.7. We denote the heights at which $\tilde{\rho}_d$, $\Theta_m$, $\mathbf{V}_{\mathbf{H}}$ and $Q_j$ are defined as *model levels*, and the heights at which $W$ is defined as *layer interfaces*. The MPAS vertical grid is defined such that the model levels lie exactly 1/2 the distance between the upper and lower interface. Based on this staggering and grid definition, the update equations {eq}`eq:3.37`–{eq}`eq:3.39` are vertically discretized:

$$
\begin{aligned}
&\Omega''^{\,\tau+\Delta\tau}(k) + C_{\Omega\Theta z}(k)\!\left[\zeta_z(k)\cdot\Theta''^{\,\tau+\Delta\tau}_m(k) - \zeta_z(k-1)\cdot\Theta''^{\,\tau+\Delta\tau}_m(k-1)\right] \\
&\qquad - C_{\Omega\Theta}(k)\cdot\Theta''^{\,\tau+\Delta\tau}_m(k) - C_{\Omega\Theta}(k-1)\cdot\Theta''^{\,\tau+\Delta\tau}_m(k-1) \\
&\qquad + C_{\Omega\rho}(k)\cdot\!\left[\tilde{\rho}''^{\,\tau+\Delta\tau}_d(k) + \tilde{\rho}''^{\,\tau+\Delta\tau}_d(k-1)\right] \\
&\qquad = \Delta\tau\,F_{\Omega^\tau}(k),
\end{aligned}
$$ (eq:A.1)

$$
\begin{aligned}
\Theta''^{\,\tau+\Delta\tau}_m(k) = {}&-\frac{1}{\Delta\zeta_w(k)}\!\Bigl[C_\Theta(k+1)\cdot\Omega''^{\,\tau+\Delta\tau}(k+1) \\
&\hphantom{{}-\frac{1}{\Delta\zeta_w(k)}\!\Bigl[\,} - C_\Theta(k)\cdot\Omega''^{\,\tau+\Delta\tau}(k)\Bigr] + \Delta\tau\,F_{\Theta^\tau_m}(k),
\end{aligned}
$$ (eq:A.2)

$$
\tilde{\rho}''^{\,\tau+\Delta\tau}_d(k) = -C_\rho(k)\!\left[\Omega''^{\,\tau+\Delta\tau}(k+1) - \Omega''^{\,\tau+\Delta\tau}(k)\right] + \Delta\tau\,F_{\tilde{\rho}^\tau}(k).
$$ (eq:A.3)

The levels $k$ indicate layer interfaces in {eq}`eq:A.1` and model levels in {eq}`eq:A.2` and {eq}`eq:A.3`. The new coefficients given in blue are defined as

$$
C_{\Omega\Theta z}(k) = \Delta\tau\,\frac{1+\epsilon}{2}\,\overline{\zeta_z}^k\,\frac{c_p}{c_v}\,\overline{\left(\frac{\rho_d^t}{\rho_m^t}\right)}^k\,\overline{\pi^t}^k\,\frac{1}{\Delta\zeta(k)}, \quad \text{interface},
$$ (eq:A.4)

$$
C_{\Omega\Theta}(k) = \Delta\tau\,\frac{1+\epsilon}{2}\,\frac{\zeta_z}{2}\,\frac{R_d}{c_v}\,g\,\tilde{\rho}^t_d\,\frac{1}{\Theta^t_M}, \quad \text{level},
$$ (eq:A.5)

$$
C_{\Omega\rho}(k) = \Delta\tau\,\frac{1+\epsilon}{2}\,\frac{\overline{\zeta_z}^{k}}{2}\,g, \quad \text{interface},
$$ (eq:A.6)

$$
C_\Theta(k) = \Delta\tau\,\frac{1+\epsilon}{2}\,\overline{\Theta^t_m}^k, \quad \text{interface},
$$ (eq:A.7)

$$
C_\rho(k) = \Delta\tau\,\frac{1+\epsilon}{2}\,\frac{1}{\Delta\zeta_w(k)}, \quad \text{level},
$$ (eq:A.8)

where we have also indicated where they are defined (levels or interfaces). The operators $\overline{\phi}^{k}$ interpolate values from model levels to layer interfaces, or from layer interfaces to model levels, using the two bounding levels or interfaces as appropriate. Interpolation from the layer interfaces to levels just averages the bounding values.

:::{admonition} MPAS code
:class: note

The coefficients {eq}`eq:A.4`–{eq}`eq:A.8` are set in subroutine `atm_compute_vert_imp_coefs` in `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
In the MPAS-A code the coefficients are named as follows:

$$
\begin{aligned}
C_{\Omega\Theta z}(k) &\;-\; \texttt{cofwz} \\
C_{\Omega\Theta}(k) &\;-\; \texttt{cofwt} \\
C_{\Omega\rho}(k) &\;-\; \texttt{cofwr} \\
C_\Theta(k) &\;-\; \texttt{coftz} \\
C_\rho(k) &\;-\; \texttt{cofrz}
\end{aligned}
$$

These coefficients are used in subroutine `atm_compute_vert_imp_coefs` to compute more coefficients (see below) for the tridiagonal solver and are used in the acoustic timestep in subroutine `atm_advance_acoustic_step` to evaluate contributions to the time level $\tau$ operators ($F_{\tilde{\rho}^\tau}$, $F_{\Theta^\tau_m}$ and $F_{\Omega^\tau}$) and in the back substitution to recover $\Theta''_m$ and $\tilde{\rho}''$ at $\tau + \Delta\tau$.
:::

Using {eq}`eq:A.2` and {eq}`eq:A.3` to eliminate $\Theta''^{\,\tau+\Delta\tau}_m$ and $\tilde{\rho}''^{\,\tau+\Delta\tau}$ in {eq}`eq:A.1` we arrive at the following equation for $\Omega''^{\,\tau+\Delta\tau}$:

$$
\begin{aligned}
\Omega''^{\,\tau+\Delta\tau}(k-1)\cdot{}&\Bigl[-C_{\Omega\Theta z}(k)\cdot C_\Theta(k-1)\cdot\zeta_z(k-1)/\Delta\zeta(k-1) \\
&\quad + C_{\Omega\rho}(k)\cdot C_\rho(k-1) \\
&\quad - C_{\Omega\Theta}(k-1)\cdot C_\Theta(k-1)/\Delta\zeta(k-1)\Bigr] \\
+\Omega''^{\,\tau+\Delta\tau}(k)\cdot{}&\Bigl[1 + C_{\Omega\Theta z}(k)\bigl(C_\Theta(k)\cdot\zeta_z(k)/\Delta\zeta(k) \\
&\qquad\quad + C_\Theta(k-1)\cdot\zeta_z(k-1)/\Delta\zeta(k-1)\bigr) \\
&\quad - C_\Theta(k)\cdot\bigl(C_{\Omega\Theta}(k)/\Delta\zeta(k) - C_{\Omega\Theta}(k-1)/\Delta\zeta(k-1)\bigr) \\
&\quad + C_{\Omega\rho}(k)\cdot\bigl(C_\rho(k) - C_\rho(k-1)\bigr)\Bigr] \\
+\Omega''^{\,\tau+\Delta\tau}(k+1)\cdot{}&\Bigl[-C_{\Omega\Theta z}(k)\cdot C_\Theta(k+1)\cdot\zeta_z(k)/\Delta\zeta(k) \\
&\quad - C_{\Omega\rho}(k)\cdot C_\rho(k) \\
&\quad + C_{\Omega\Theta}(k)\cdot C_\Theta(k+1)/\Delta\zeta(k)\Bigr] \\
={}&+ \Delta\tau\Bigl[F_{\Omega^\tau}(k) - \bigl(C_{\Omega\Theta z}(k)\cdot\zeta_z(k) - C_{\Omega\Theta}(k)\bigr)F_{\Theta^\tau_m}(k) \\
&\quad + \bigl(C_{\Omega\Theta z}(k)\cdot\zeta_z(k-1) + C_{\Omega\Theta}(k-1)\bigr)F_{\Theta^\tau_m}(k-1) \\
&\quad - C_{\Omega\rho}(k)\bigl(F_{\tilde{\rho}^\tau}(k) + F_{\tilde{\rho}^\tau}(k-1)\bigr)\Bigr].
\end{aligned}
$$ (eq:A.9)

This equation can also be written as

$$
C_{-}\cdot\Omega''^{\,\tau+\Delta\tau}(k-1) + C\cdot\Omega''^{\,\tau+\Delta\tau}(k) + C_{+}\cdot\Omega''^{\,\tau+\Delta\tau}(k+1) = R,
$$ (eq:A.10)

where the three coefficients are the three bracketed terms in {eq}`eq:A.9`. This equation is a tridiagonal system with the diagonal coefficient $C$ and the lower and upper diagonal coefficients $C_{-}$ and $C_{+}$. In MPAS this system is solved with a standard tridiagonal solver with the boundary condition $\Omega = 0$ at the upper and lower boundaries.

:::{admonition} MPAS code
:class: note

The coefficients for the tridiagonal solver, $C$, $C_{-}$ and $C_{+}$, are computed in subroutine `atm_compute_vert_imp_coefs` in `src/core_atmosphere/dynamics/mpas_atm_time_integration.F`.
In the code, $C_{-}$ is stored in the array `a_tri`, $C$ is stored in the array `b_tri`, and $C_{+}$ is stored in the array `c_tri`. The coefficients are used to precompute other coefficients (`alpha_tri` and `gamma_tri`) that are used along with `a_tri` to perform the implicit solve in subroutine `atm_advance_acoustic_step`.
:::
