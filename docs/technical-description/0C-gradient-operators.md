# Appendix C: Gradient Operators

Horizontal gradients on the MPAS mesh are needed in many places in MPAS, for example in the calculation of the horizontal deformation and in the calculation of the horizontal components of the frontogenesis function. In this technical note we have already described horizontal gradients across cell edges (cell-centered quantities differenced between cells sharing the edge), gradients along cell edges (e.g. vorticity differenced between vertices), and the horizontal flux divergence operator applied to momentum, coupled potential temperature, and scalar mass. We describe the finite-volume computation of arbitrary horizontal gradients in the first section in the first section of this appendix and there use in specific operators in the latter sections of this appendix.

## C.1 Horizontal Derivatives

Within the finite-volume formulation of MPAS, we define gradient operators as the cell-average of the gradient. Using Green's theorem,

$$
\iint_R\left(\frac{\partial N}{\partial x} - \frac{\partial M}{\partial y}\right)dx\,dy = \oint_C(M\,dx + N\,dy),
$$

with the components

$$
\iint_R\frac{\partial N}{\partial x}\,dx\,dy = \oint_C N\,dy
$$ (eq:C.1)

and

$$
\iint_R -\frac{\partial M}{\partial y}\,dx\,dy = \oint_C M\,dx,
$$ (eq:C.2)

we transform cell-integrated values of derivatives into line integrals of the differentiated quantities. Specifically, we take the line integral in the counterclockwise direction in this application, and the integral is over a closed curve $C$, completely around the cell. Dividing the discrete evaluation of {eq}`eq:C.1` or {eq}`eq:C.2` by the area recovers the cell-averaged gradient.

To compute the discrete horizontal derivatives, we begin by defining a tangent plane to the cell center and map the vertices of the cell to the tangent plane, as described in figure C.1. We define the plane and the locations of the vertices $(x_i, y_i)$ on the tangent plane for vertices $i = 1,\ldots,n$, where $n$ is the number of vertices (and the number of edges). The vertices are

**[Figure C.1: MPAS tangent-plane schematic. To be added next session.]**

ordered such that they progress counterclockwise around the cell; likewise edges are ordered counterclockwise around the cell with edge $e_i$ spanning vertices $(x_i, y_i)$ and $(x_{i+1}, y_{i+1})$ as depicted in figure C.2. The line integrals take place in the $(x, y)$ coordinate system defined in figure C.2, where the $x$ direction is unconstrained, e.g. we do not require any vertex to lie on the $x$ axis nor do we require any normal velocity to lie on the axis. In the MPAS implementation, we define the positive direction $x$ as tangent to latitude circles at the cell edge points and increasing in the direction of increasing longitude, and the direction $y$ as tangent to longitudes and increasing with increasing latitude, i.e. the local latitude-longitude coordinate. Thus the $x$ derivative is the longitudinal derivative and the $y$ derivative is the latitudinal derivative.

Referring to figure C.2, we define the length increments for the edge from vertex $i$ to vertex $i + 1$ as

$$
\begin{aligned}
\Delta x &= -l_e\sin\theta_e, \\
\Delta y &= l_e\cos\theta_e,
\end{aligned}
$$

where

$$
l_e = \sqrt{(x_{i+1} - x_i)^2 + (y_{i+1} - y_i)^2}
$$

is the edge length on the tangent plane for edge $e$. We can now define the gradient operators for cell-centered quantities as

$$
A_c\phi_x = \sum_{e=1}^{n}l_e\cos\theta_e\,\overline{\phi}^e,
$$ (eq:C.3)

$$
A_c\phi_y = \sum_{e=1}^{n}l_e\sin\theta_e\,\overline{\phi}^e,
$$ (eq:C.4)

where $\overline{\phi}^e$ denotes an average of the cell-centered values of $\phi$ from the cells sharing the edge and $n$ is the number of edges on the cell. Next we define coefficients:

$$
c_{x_{c,e}} = \frac{l_e\cos\theta_e}{A_c}
$$ (eq:C.5)

$$
c_{y_{c,e}} = \frac{l_e\sin\theta_e}{A_c}
$$ (eq:C.6)

**[Figure C.2: Horizontal depiction of the tangent plane and cell-edge velocity fields for an arbitrary MPAS mesh cell. To be added next session.]**

To compute the line integrals {eq}`eq:C.1` and {eq}`eq:C.2` for the horizontal velocities $(u, v)$ in the tangent plane coordinate, we rotate the edge velocities into the coordinate:

$$
\begin{aligned}
u &= u_e\cos\theta_e - v_e\sin\theta_e, \\
v &= u_e\sin\theta_e + v_e\cos\theta_e.
\end{aligned}
$$

The edge velocity $u_e$ is recovered from the prognosed momentum $\tilde{\rho}_du_e$, and the edge tangential velocity is diagnosed as described in section 4.3.2. Note that in MPAS the normal velocity may or may not be positive outward from a given cell, and when pointing inward the signs are shifted given the 180 degree rotation of both $u_e$ and $v_e$. Given these definitions, we can compute the cell-averaged derivates of the velocity components $(u, v)$ on the tangent plane:

$$
\begin{aligned}
A_cu_x &= \sum_{e=1}^{n}l_e\cos\theta_e(u_e\cos\theta_e - v_e\sin\theta_e) \\
&= \sum_{e=1}^{n}l_e(u_e\cos^2\theta_e - v_e\cos\theta_e\sin\theta_e)
\end{aligned}
$$ (eq:C.7)

$$
\begin{aligned}
A_cu_y &= \sum_{e=1}^{n}l_e\sin\theta_e(u_e\cos\theta_e - v_e\sin\theta_e) \\
&= \sum_{e=1}^{n}l_e(u_e\sin\theta_e\cos\theta_e - v_e\sin^2\theta_e)
\end{aligned}
$$ (eq:C.8)

$$
\begin{aligned}
A_cv_x &= \sum_{e=1}^{n}l_e\cos\theta_e(u_e\sin\theta_e + v_e\cos\theta_e) \\
&= \sum_{e=1}^{n}l_e(u_e\cos\theta_e\sin\theta_e + v_e\cos^2\theta_e)
\end{aligned}
$$ (eq:C.9)

$$
\begin{aligned}
A_cv_y &= \sum_{e=1}^{n}l_e\sin\theta_e(u_e\sin\theta_e + v_e\cos\theta_e) \\
&= \sum_{e=1}^{n}l_e(u_e\sin^2\theta_e + v_e\sin\theta_e\cos\theta_e)
\end{aligned}
$$ (eq:C.10)

where $A_c$ is the cell area and $n$ is the number of edges on each cell. There are only three distinct coefficients multiplying velocities $u_e$ and $v_e$ in {eq}`eq:C.7`–{eq}`eq:C.10`: $\cos^2\theta_e$, $\sin^2\theta_e$, and $\cos\theta_e\sin\theta_e$, all three multiplied by $l_e$.

## C.2 Horizontal deformation for the Smagorinsky scheme

The 2D (horizontal) Smagorinsky viscosity calculation uses the horizontal deformation term in the formula for the horizontal eddy viscosity $K$:

$$
K = c_s^2l^2\!\left[\left(\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}\right)^{\!2} + \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^{\!2}\right]^{1/2}.
$$ (eq:C.11)

Defining line-integral coefficients $a_{c,e}$, $b_{c,e}$, and $c_{c,e}$ for cell $c$ and edge $e$ as

$$
\begin{aligned}
a_{c,e} &= \bigl[l_e\cos^2\theta_e\bigr]_c/A_c, \\
b_{c,e} &= \bigl[l_e\sin^2\theta_e\bigr]_c/A_c, \\
c_{c,e} &= \bigl[l_e\sin\theta_e\cos\theta_e\bigr]_c/A_c,
\end{aligned}
$$

we can re-write the derivatives for cell $c$ as

$$
u_x = \sum_{e=1}^{n}(u_e a_{c,e} - v_e c_{c,e}),
$$ (eq:C.12)

$$
u_y = \sum_{e=1}^{n}(u_e c_{c,e} - v_e b_{c,e}),
$$ (eq:C.13)

$$
v_x = \sum_{e=1}^{n}(u_e c_{c,e} + v_e a_{c,e}),
$$ (eq:C.14)

$$
v_y = \sum_{e=1}^{n}(u_e b_{c,e} + v_e c_{c,e}).
$$ (eq:C.15)

In the Smagorinsky eddy viscosity {eq}`eq:C.11`, the terms in the deformation operator can be written in terms of the discrete velocity derivatives {eq}`eq:C.12`–{eq}`eq:C.15`:

$$
\left(\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}\right) \;\to\; \sum_{e=1}^{n}\bigl((a_{c,e} - b_{c,e})u_e - 2c_{c,e}v_e\bigr)
$$ (eq:C.16)

$$
\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right) \;\to\; \sum_{e=1}^{n}\bigl(2c_{c,e}u_e + (a_{c,e} - b_{c,e})v_e\bigr)
$$ (eq:C.17)

In MPAS we pre-compute the coefficients $(a_{c,e} - b_{c,e})$ and $2c_{c,e}$ for all edges on each cell used in computing the horizontal derivatives of the velocities and store them for use during the integration.

:::{admonition} MPAS code
:class: note

The two coefficients in {eq}`eq:C.16` and {eq}`eq:C.17`, $(a_{c,e} - b_{c,e})$ and $2c_{c,e}$, are computed in subroutine `atm_initialize_deformation_weights` in `MPAS/src/core_init_atmosphere/mpas_atm_advection.F` in subroutine `atm_initialize_deformation_weights`. The coefficients are stored in the arrays `defc_a` $(a_{c,e} - b_{c,e})$ and `defc_b` $(2c_{c,e})$. The Smagorinsky viscosity {eq}`eq:C.11` is computed in subroutine `atm_compute_dyn_tend` in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F` as are the 2nd-order filter tendencies for the dynamics variables $u$, $w$, and $\theta_m$ for the RK3 time integration.
:::

## C.3 Frontogenesis Function

The frontogenesis function is given in section 6.2 and we repeat it here:

$$
\frac{1}{2}\left|\frac{\partial\nabla\theta}{\partial t}\right| = \frac{1}{2}\!\left[-(u_x + v_y)\,|\nabla\theta|^2 + 2(u_y + v_x)(\theta_x\theta_y) - (u_x - v_y)(\theta_x^2 - \theta_y^2)\right].
$$ (eq:C.18)

There are three quantities involving the velocity needed to evaluate this function: $(u_x + v_y)$, $(u_y + v_x)$, and $(u_x - v_y)$. The first is the velocity divergence, and this is computed using the standard divergence operator for a cell:

$$
(u_x + v_y) = \nabla_\zeta\cdot\mathbf{v} = \frac{1}{A_i}\sum_{n_{e_i}}(L_{e_i}\mathbf{v}_{H,e_i}\cdot\mathbf{n}),
$$ (eq:C.19)

where $n_e$ are the number of edges for the cell and $L_{e_i}$ is the length of each edge. In MPAS the velocity divergence is computed on the sphere (i.e. $L_{e_i}$ are the spherical arc lengths of the edge as opposed to the length of the edges on the tangent plane in figure C.2.

The second term, $(u_y + v_x)$, is the shearing deformation, and it is computed as in the Smagorinsky deformation using {eq}`eq:C.17`. The third term, $(u_x - v_y)$, is the stretching deformation term and it is computed as in the Smagorinsky deformation using {eq}`eq:C.16`.

Derivatives of the potential temperature $(\theta_x, \theta_y)$ are also needed in the frontogenesis function {eq}`eq:C.18`. These are computed using {eq}`eq:C.3` and {eq}`eq:C.4`.

:::{admonition} MPAS code
:class: note

The coefficients involving hearing deformation and the stretching deformation are computed stored as described in the case of the Smagorinsky eddy viscosity calculation. The coefficients for the scalar gradient calculations, {eq}`eq:C.5` and {eq}`eq:C.6`, are computed in subroutine `atm_initialize_deformation_weights` in `MPAS/src/core_init_atmosphere/mpas_atm_advection.F` in subroutine `atm_initialize_deformation_weights`. The coefficients are stored in the arrays `cell_gradient_coef_x` $(c_{x_{c,e}})$ and `cell_gradient_coef_y` $(c_{y_{c,e}})$. They are not presently used in MPAS-A but are used in the frontogenesis computation in CAM-MPAS.
:::
