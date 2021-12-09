# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 20:04:54 2021

@author: Konstantinos
"""

# We'll solve the Newtonian hydrostatic equilibrium for a WD

import numpy as np # Arrays
import matplotlib.pyplot as plt # Plotting

from scipy.integrate import solve_ivp # ODE Solver
#%% Using scipy
mu = 2 # WDs are charge neutral so you got half an electron per baryon.

# Defining the ODEs
def dSdx(r, S):
    p, m = S
    G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
    c = 2.997e10 # cgs [cm/s]
    if r==0 or m==0:
        print('zero')
        return [0,0] #dP/dr = 0, dM/dr = 0 at the center
    else:
        rho = (p/K)**(1/gamma)
        oros1 = - G*rho*m/r**2
        oros2 = 1+p/(rho*c**2)
        oros3 = 1+ 4*np.pi*p*r**3/(m*c**2)
        oros4 = 1 - 2*G*m / (r*c**2)
        dpdx = oros1*oros2*oros3/oros4
        dmdx =  4*np.pi *r**2 * rho# dM/dr = 4πr^2 (P/K)^(1/Γ)
        return [dpdx, dmdx]

# Event to stop the integration when the Pressure is low enough.
def starend(t,y):
    return y[0] - 1e15 # Stop integrating when P = 1e15, otherwise the graph looks real ugly.
starend.terminal = True # Stop the integration when starend event is True
starend.direction= -1 # Looks for sign changes + -> -

# Distinguish between rel. and non-rel. 
Central_Density = 1e15
if Central_Density>6e15: #cutoff is at 1e16.
    gamma = 4/3 # Relativistic electrons
    K = 1.2435e15/mu**gamma
    Central_Pressure = K*Central_Density**gamma
else:
    gamma = 5/3 # Non-relativistic electrons
    K = 5.3802e9/mu**gamma # Shapiro, Teukolsky p.28
    Central_Pressure = K*Central_Density**gamma

# Initial Conditions at r=0
S_0 = [Central_Pressure,1]
dt = 1000 # Good enough for a fine grid and quick integration.
sol = solve_ivp(fun=dSdx, t_span=(1,1e10), max_step=dt, y0=S_0, method='RK45', 
            atol=0.01, rtol=0.001,  events=starend, dense_output=True) # actually solves it

# Extract the values from the solver, rescale to the desired units
p_sol = sol.y[0]
m_sol = sol.y[1]/(1.989*1e33) # Solar Masses
x = sol.t*1e-5 # Km

#Plotting
plt.figure(1)
den_text = "{:.1e}".format(Central_Density)
fig = plt.plot(x,p_sol, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.legend(loc='best')
plt.title('Pressure of a Newtonian Neutron Star')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")
plt.yscale('log')

plt.figure(2)
fig = plt.plot(x,m_sol, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.legend(loc='best')
plt.title('Mass of a Newtonian Neutron Star')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\oplus}]$")