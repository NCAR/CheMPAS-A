# Appendix B: Least-Squares Polynomial Computation

**[Figure B.1: A limited region in an MPAS horizontal mesh. Adapted from Skamarock and Gassmann (2011) Figure 1. To be added next session.]**

Consider the Voronoi cells depicted in Figure B.1. To compute the flux edge $e_1$ we fit polynomials to cells $c_1$ and its neighbor cells, i.e. $(c_1, c_0, c_6, c_7, c_8, c_9, c_2)$, and also to $c_0$ and its neighbors, i.e. $(c_0, c_1, c_2, c_3, c_4, c_5, c_6)$. More generally, a polynomial fit is needed for each cell, and its 2nd derivative normal to each edge needs to be computed.

To illustrate how the polynomials are computed, we examine the polynomial computation for the cell $c_0$. We begin by constructing a planar projection of the cell values on the sphere using a tangent plane that is coincident with the sphere at the cell center point $c_0$. We then define a 2 dimensional $(x, y)$ coordinate system on the plane such that the $x$-axis extends through the first cell on the cell neighbor list for $c_0$ (see Figure B.2) and define the $y$-axis by the right-hand rule.

Following Skamarock and Gassmann (2011), the tangent plane is defined such that the angles $\theta_i$ are the spherical angles at $c_0$ between the vectors $\overrightarrow{c_0c_1}$ and $\overrightarrow{c_0c_i}$. The length $|\overrightarrow{c_0c_i}|$ is the great circle arc distance on the sphere between $c_0$ and $c_i$. We fit a least-squares-fit polynomial for a

**[Figure B.2: Tangent plane for cell $c_0$. To be added next session.]**

variable $\phi$ of the form

$$
\phi = a_0 + a_xx + a_yy + a_{xx}x^2 + a_{xy}xy + a_{yy}y^2
$$ (eq:B.1)

using the values from the cell $c_0$ and its neighbors. Defining the cell-centered values as $\mathbf{s} = [\phi_0,\ldots,\phi_6]^{\mathrm{T}}$ and the polynomial coefficients $\mathbf{f} = [a_0,\ldots,a_{yy}]^{\mathrm{T}}$, we can define the least-squares-fit polynomial following Strang (1980):

$$
\mathbf{f} = \bigl(\mathbf{P}^{\mathrm{T}}\mathbf{P}\bigr)^{-1}\mathbf{P}^{\mathrm{T}}\mathbf{s} = \mathbf{B}\mathbf{s},
$$ (eq:B.2)

where the matrix $\mathbf{P}$ is an $m \times n$ matrix defined as

$$
\mathbf{P} = \begin{pmatrix} 1 & x_0 & y_0 & x_0^2 & x_0 y_0 & y_0^2 \\ \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\ 1 & x_m & y_m & x_m^2 & x_m y_m & y_m^2 \end{pmatrix}.
$$

$m$ is equal to the number of cells used in the polynomial fit (e.g. $m = 7$ for the configuration shown in figure B.2, 6 when the central cell is a pentagon, and 8 when the central cell is a heptagon) and $n$ is the number of coefficients in the polynomial {eq}`eq:B.1`, in this case 6.

The matrix $\mathbf{B} = \bigl(\mathbf{P}^{\mathrm{T}}\mathbf{P}\bigr)^{-1}\mathbf{P}^{\mathrm{T}}$ is an $n \times m$ matrix that is used to compute the coefficients of the polynomial given the cell values. In the horizontal flux divergence {eq}`eq:4.7`, we need the second derivatives of the polynomials in the direction normal to each cell face. Employing the chain rule, the 2nd derivative of the polynomial {eq}`eq:B.1` in the $x'$ direction, where $x'$ is a rotation of the coordinate $x$ through an angle $\theta_{e_i}$, is

$$
\frac{\partial^2\phi}{\partial x'^2} = 2\bigl(a_{xx}\cos^2\theta_{e_i} + a_{xy}\cos\theta_{e_i}\sin\theta_{e_i} + a_{yy}\sin^2\theta_{e_i}\bigr),
$$ (eq:B.3)

and we can write this using the $\mathbf{B}$ matrix as

$$
\frac{\partial^2\phi}{\partial x'^2} = \sum_{i=1}^{m}\left[\phi_i\bigl(2\mathbf{B}_{4,i}\cos^2\theta_{e_i} + 2\mathbf{B}_{5,i}\cos\theta_{e_i}\sin\theta_{e_i} + 2\mathbf{B}_{6,i}\sin^2\theta_{e_i}\bigr)\right].
$$ (eq:B.4)

For each cell we compute and store the coefficients for each edge $e$ on the cell (the term in parentheses in {eq}`eq:B.4`) and use these in the construction of the flux {eq}`eq:4.7`.

:::{admonition} MPAS code
:class: note

The computation of the least-squares-fit 2nd derivative weights occurs within the *init_atmosphere* core in `src/core_init_atmosphere/mpas_atm_advection.F`. The 2nd derivative coefficients are stored in the array `deriv_two(source_cell,cells,edges)`, where `edges` are all the edges on the mesh, `cells` are the two cells that share the edge, and `source_cells` are the cells used in the least-squares fit polynomial. These operations are performed in subroutine `atm_initialize_advection_rk`.
:::
