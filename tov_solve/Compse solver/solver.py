# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 18:23:32 2021

@author: Konstantinos

This code solves the TOV eqs using a tabulated EOS from CompOSE.
The CompOSE data needs some clean up, which is detailed below
"""
import numpy as np # Arrays
import matplotlib.pyplot as plt # Plotting
import pandas as pd # CSV import

from scipy.integrate import solve_ivp # ODE Solver
from scipy.interpolate import interp1d # Interpolator

#%% Read the Data
# In order to clean up the data, one needs to delete the first two lines in the .thermo.
# Delete the first two lines in the .nb and make sure the first line is .nb
# This is not a great coding practice.

df = pd.read_csv('eos.thermo',sep='\s+', header=None) # Reads the Compose EoS, P/nb is in MeV
names = [ 'it','inb','iyq','P/nb' , 
         's/nb','mub/mn-1','muq/mn','f/nb*mn-1','e/nb*mn-1','idk','idk',]
df.columns = names # Change the column names

df_nb = pd.read_csv('eos.nb') # The nb is on another file, units [fm^-3]

#%% Get the Pressure in CGS

pres = [ df_nb.loc[i,' nb']*df.loc[i,'P/nb'] # P= nb * P/nb, The name is 'space nb', units are [MeV/fm^3]
        for i in range(df.index[-1]+1) ] # last index is df.index[-1] 

pres_cgs =  [ pres[i] * 1.60218e-6 * 1e39 
             for i in range(len(pres)) ] 
# MeV -> 1.60218e-6 ergs
# 1/fm^3 -> 1e39 1/cm^3
# units are cgs [erg/cm^3]

#%% Get the Density in CGS     
                       
mn = 1.674e-24 # Neutron mass, units are cgs, [g]
den = [ mn*df_nb.loc[i,' nb']
       for i in range(df.index[-1]+1) ] # units are [g/fm^3]

den_cgs =  [ den[i] * 1e39 
             for i in range(len(den)) ]
# 1/fm^3 -> 1e39 1/cm^3
# units are cgs [erg/cm^3]
#%% Interpolate

pressure_inter =  interp1d(den_cgs, pres_cgs, kind='cubic')
den_inter = interp1d(pres_cgs, den_cgs, kind='cubic')
#%% One Neutron Star

# The RuntimeWarning occurs at the end of every integration, whenever P is negative
# it is nothing to worry about, I think.
def dSdx(r, S):
    p, m = S
    G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
    c = 2.997e10 # cgs [cm/s]
    if r==0 or m==0:
        print('the code understood we are at the core')
        return [0,0] # dP/dr = 0, dM/dr = 0 at the center
    else:
        try:
            rho = den_inter(p) # Interpolated Density
            print('Hi P is:',"{:.1e}".format(p),'while ρ is',"{:.1e}".format(rho))
            oros1 = - G*rho*m/r**2
            oros2 = 1+p/(rho*c**2)
            oros3 = 1+ 4*np.pi*p*r**3/(m*c**2)
            oros4 = 1 - 2*G*m / (r*c**2)
            dpdx = oros1*oros2*oros3/oros4
            dmdx =  4*np.pi *r**2 * rho # dM/dr = 4πr^2 (P/K)^(1/Γ)
        except ValueError:
            print('valerr')
            return [0,0]
    return [dpdx, dmdx]
    
def starend(t,y):
    return y[0] - 4e25 # Stop integrating when P = 4e24, otherwise we get out of Interpolation range. (3e34 compiles)
starend.terminal = True # Stop the integration when starend event is True
starend.direction= -1 # Looks for sign changes + -> -

Central_Density = 2e15 # <2e15 to be inside of interpolation range
M_tot = []
R_tot = []
Central_Pressure = pressure_inter(Central_Density)    
S_0 = [Central_Pressure,0.01] # Initial values [Pressure, Mass]
sol = solve_ivp(fun=dSdx, t_span=(0,1e10), max_step=1000 ,y0=S_0, method='RK45', 
                events=starend)

# Extract the values from the solver, rescale to the desired units
p_sol = sol.y[0]
m_sol = sol.y[1]/(1.989*1e33) # Solar Masses
x = sol.t*1e-5 # Km
M_tot.append(max(m_sol)) # Total mass
R_tot.append(max(x)) # Radius

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

#%% EoS Comparison
K = 5.3802e9 #non-rel neutrons.
gamma = 5/3 #non-rel neutrons.
plt.figure(4)
plt.plot(den_cgs, pressure_inter(den_cgs), color='darkorange', label='Compose EoS')
plt.scatter(den_cgs, pres_cgs, color='black')

analytic_pressure =[K*den_cgs[i]**gamma for i in range(len(den_cgs))]
plt.plot(den_cgs, analytic_pressure, color='darkblue', label='non-rel Fermi Gas')
plt.title("Comparison of EoS of non-interacting, non-relativistic neutron Fermi Gas vs Tabulated from Haensel & Douchin ")
plt.legend(loc='best')
plt.xlabel('Density $[g/cm^3]$')
plt.ylabel('Pressure $[erg/cm^3]$')
plt.yscale('log')
plt.xscale('log')
# The cutoff for rel. neutrons is at ρ=6e15. Τhe biggest density in the tabulated EoS is 2e15. This means I don't have to care
# about them.
#%% Many Neutron Stars.

Central_Density = np.logspace(np.log10(2e14),np.log10(2e15),20) # <5e15 to be inside of interpolation range. The Liquid core starts at 1.4e14
M_tot = []
R_tot = []
for i in range(len(Central_Density)):
    Central_Pressure = pressure_inter(Central_Density[i])    
    S_0 = [Central_Pressure,0.01] # Initial values [Pressure, Mass]
    sol = solve_ivp(fun=dSdx, t_span=(0,1e10), max_step=1000 ,y0=S_0, method='RK45', 
                    events=starend)
    
    # Extract the values from the solver, rescale to the desired units
    p_sol = sol.y[0]
    m_sol = sol.y[1]/(1.989*1e33) # Solar Masses
    x = sol.t*1e-5 # Km
    M_tot.append(max(m_sol)) # Total mass
    R_tot.append(max(x)) # Radius
    
    den_text = "{:.1e}".format(Central_Density[i])
    plt.figure(1)
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
    plt.ylabel("Mass $[M_{\odot}]$")
        
plt.figure(3) 
plt.plot(R_tot, M_tot, 'o')
plt.title('Total Masses of NS\'s of varying central densities')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\odot}]$")