# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 13:29:26 2021

@author: Konstantinos
"""
# We'll solve the Newtonian hydrostatic equilibrium for a WD

import numpy as np # Arrays
import matplotlib.pyplot as plt # Plotting

from scipy.integrate import solve_ivp # ODE Solver
#%% Using scipy
mu = 2 # WDs are charge neutral so you got half an electron per baryon.

def dSdx(r, S):
    p, m = S
    G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
    if r==0:
        return [0,0] #dP/dr = 0, dM/dr = 0 at the center
    else:
        dpdx = -G*m* p**(1/gamma)/(r**2 * K**(1/gamma)) # dP/dr  -G M P^(1/Γ) / r^2 K^(1/Γ)
        dmdx =  4* p**(1/gamma)* np.pi *r**2 /K**(1/gamma) # dM/dr = 4πr^2 (P/K)^(1/Γ)
        return [dpdx, dmdx]
    
def starend(t,y):
    return y[0] - 1e17 # Stop integrating when P = 1e17, otherwise the graph looks real ugly.
starend.terminal = True # Stop the integration when starend event is True
starend.direction= -1 # Looks for sign changes + -> -

Central_Density = np.logspace(4,8,5) # Range of different central densities. 
Central_Pressure = np.zeros( len(Central_Density)) # Empty array, same length as densities.


for i in range(len(Central_Density)):
    # This checks whether we have relativistic electrons
    # or not, and adjusts the pressure, K, Γ, accordingly
    if Central_Density[i]>=9e5: #cutoff is at 1e6.
        gamma = 4/3 # Relativistic electrons
        K = 1.2435e15/mu**gamma
        Central_Pressure[i] = K*Central_Density[i]**gamma
    else:
        gamma = 5/3 # Non-relativistic electrons
        K = 1.0036e13/mu**gamma
        Central_Pressure[i] = K*Central_Density[i]**gamma
        
    S_0 = [Central_Pressure[i],0]
    dt = 100000 # Good enough for a fine grid and quick integration.
    sol = solve_ivp(fun=dSdx, t_span=(0,1e10), max_step=dt, y0=S_0, method='RK45', 
                atol=0.01, rtol=0.001, events=starend, dense_output=True)

    # Extract the values from the solver, rescale to the desired units
    p_sol = sol.y[0]
    m_sol = sol.y[1]/(1.989*1e33) # Solar Masses
    x = sol.t*1e-5 # Km

    # Plotting
    colarr=['tab:blue','tab:orange','tab:green','tab:red','tab:purple'] # Plotting Colours
    
    plt.figure(1)
    den_text = "{:.1e}".format(Central_Density[i])
    fig = plt.plot(x,p_sol, color=colarr[i], 
                   label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='best')
    plt.title('Pressure of a Newtonian White Dwarf')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Pressure $[erg/cm^3] $")
    plt.yscale('log')
    
    plt.figure(2)
    fig = plt.plot(x,m_sol, color=colarr[i],
                   label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='best')
    plt.title('Mass of a Newtonian White Dwarf')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Mass $[M_{\oplus}]$")
    

