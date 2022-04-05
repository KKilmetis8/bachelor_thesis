# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 15:23:51 2021

@author: Konstantinos
"""

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
# The RuntimeWarning occurs at the end of every integration, whenever P is negative
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
    return y[0] - 1e17 # Stop integrating when P = 1e15, otherwise the graph looks real ugly.
starend.terminal = True # Stop the integration when starend event is True
starend.direction= -1 # Looks for sign changes + -> -

Central_Density = np.logspace(14,17.5,10) # Range of different central densities. 
Central_Pressure = np.zeros( len(Central_Density)) # Empty array, same length as densities.
M_tot = []
R_tot = []
for i in range(len(Central_Density)):
    # This checks whether we have relativistic neutrons
    # or not, and adjusts the pressure, K, Γ, accordingly
    if Central_Density[i]>=6e15: #cutoff is at 6e15.
        gamma = 4/3 # Relativistic neutrons
        K = 1.2293e15 # Shapiro, Teukolsky p.28
        Central_Pressure[i] = K*Central_Density[i]**gamma
    else:
        gamma = 5/3 # Non-relativistic neutrons
        K = 5.3802e9 # Shapiro, Teukolsky p.28
        Central_Pressure[i] = K*Central_Density[i]**gamma
        
    S_0 = [Central_Pressure[i],0] # Initial values [Pressure, Mass]
    sol = solve_ivp(fun=dSdx, t_span=(0,5e6), y0=S_0, method='RK45', 
                    events=starend, dense_output=True)

    # Extract the values from the solver, rescale to the desired units
    p_sol = sol.y[0]
    m_sol = sol.y[1]/(1.989*1e33) # Solar Masses
    x = sol.t*1e-5 # Km
    M_tot.append(max(m_sol)) # Total mass
    R_tot.append(max(x)) # Radius
    
    plt.figure(1)
    den_text = "{:.1e}".format(Central_Density[i])
    fig = plt.plot(x,p_sol, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='lower left')
    plt.title('Pressure of a Newtonian Neutron Star')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Pressure $[erg/cm^3] $")
    plt.yscale('log')
    plt.savefig('nsn_pr.pdf', format='pdf')
    
    plt.figure(2)
    fig = plt.plot(x,m_sol, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='best')
    plt.title('Mass of a Newtonian Neutron Star')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Mass $[M_{\oplus}]$")
    plt.savefig('nsn_mr.pdf', format='pdf')
#%%
plt.figure(3) 
plt.plot(R_tot, M_tot, 'o')
plt.title('Total Masses of Newtonian NS\'s of varying central densities')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\oplus}]$")
plt.savefig('nsn_max.pdf', format='pdf')