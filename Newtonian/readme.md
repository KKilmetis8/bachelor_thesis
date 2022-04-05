## Newtonian Equations

You gotta walk before you run.
These solve the hydrostatic equilibrium for
newtonian WDs and NSs.

$$
\frac{dP}{dr}=-G\frac{\rho (r)M(r)}{r^2}
$$

$$
\frac{dM}{dr}=4\pi r^2 \rho (r)
$$

I've also included my implementation of a simple RK4 (`wd_rk4`) algorithm which works
well enough, considering it can't change its step.