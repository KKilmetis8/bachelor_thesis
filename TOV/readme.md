## Solving the Tolman-Oppenheimer-Volkov equations

In order to study Neutron Stars we can attempt to find out the Pressure and the Mass at any given radius.

Since NSs are compact objects, they are described by the TOV eqs which are the spherically symmetrical solution to Einstein's field equations. These are:

$$
\frac{dP}{dr} = -\frac{G\rho m}{r^2} \left( 1 +\frac{P}{\rho c^2} \right) \left( 1 + \frac{4\pi Pr^3}{mc^2} \right)  \left(1 -\frac{2Gm}{r c^2} \right)^{-1}
$$

$$
\frac{dM}{dr} = 4\pi r^2 \rho
$$

In order to solve this system of differential equations, we need to supply it with an EoS. An EoS relates *density ($\rho$)* with *pressure ($P$)*. 

`ns_tov.py` solves for a degenarate neutron fermi gas
`ns_tovinter.py` solves for a tabulated EoS (SLy4, contained in `haenselEOS.csv`)

The directory CompOSE_solvers contains the solvers for EoS taken from
CompOSE. The solvers work identically, the sole differentiatior is the input method.