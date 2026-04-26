# Chapter 6: Diagnostics

Description of diagnostic quantities computed by MPAS.

## 6.1 Zonal and Meridional Velocity

Description of radial basis function (RBF) formulation to compute the cell-centered zonal and meridional velocities.

## 6.2 Frontogenesis Function

A frontogenesis function is used in the Community Earth System Model (CESM) within the gravity-wave drag parameterization in the physics of CESM's Community Atmosphere Model (CAM). In continuous form, the frontogenesis function is

$$
\frac{1}{2}\left|\frac{\partial\nabla\theta}{\partial t}\right| = \frac{1}{2}\!\left[-(u_x + v_y)\,|\nabla\theta|^2 + 2(u_y + v_x)(\theta_x\theta_y) - (u_x - v_y)(\theta_x^2 - \theta_y^2)\right]
$$ (eq:6.1)

This frontogenesis function considers only the horizontal component of frontogenesis, and it is based on the formulation described in Richter et al. (2010). The finite-volume computation of components of the function are described in Appendix C.

## 6.3 Other Diagnostics to Follow
