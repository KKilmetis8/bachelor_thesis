# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 15:09:29 2021

@author: Konstantinos
"""
import numpy as np # Arrays
import matplotlib.pyplot as plt # Plotting
import pandas as pd # CSV import

from scipy.integrate import solve_ivp # ODE Solver
from scipy.interpolate import InterpolatedUnivariateSpline as us # Interpolator
#%% Interpolate

df = pd.read_csv('haenselEOS.csv',sep=';', header=0) # Reads the Haensel & Douchin '01 EOS + the Haensel & Pichon '94 EOS for the outer crust
pressure_inter =  us(df.loc[:,'rhog/cm3'],df.loc[:,'Perg/cm3'])
den_inter = us(df.loc[:,'Perg/cm3'],df.loc[:,'rhog/cm3'])
#%% One Neutron Star
starno = 0
# The RuntimeWarning occurs at the end of every integration, whenever P is negative
# it is nothing to worry about, I think.
def dSdx(r, S):
    global starno
    p, m = S
    G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
    c = 2.997e10 # cgs [cm/s]
    if r==0 or m==0:
        starno = starno + 1
        print('the code understood we are at the core ' + 'star no. ',starno)
        return [0,0] # dP/dr = 0, dM/dr = 0 at the center
    else:
        try:
            rho = den_inter(p) # Interpolated Density
            oros1 = - G*rho*m/r**2
            oros2 = 1+p/(rho*c**2)
            oros3 = 1+ 4*np.pi*p*r**3/(m*c**2)
            oros4 = 1 - 2*G*m / (r*c**2)
            dpdx = oros1*oros2*oros3/oros4
            dmdx =  4*np.pi *r**2 * rho # dM/dr = 4πr^2 (P/K)^(1/Γ)
        except ValueError:
            return [0,0]
    return [dpdx, dmdx]
    
def starend(t,y):
    return y[0] # Stop integrating when P = 0
starend.terminal = True # Stop the integration when starend event is True
starend.direction= -1 # Looks for sign changes + -> -

Central_Density = 4.4e15 # <5e15 to be inside of interpolation range
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
plt.title('Pressure of a Newtonian Neutron Star, Sly 4 H&D')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")
plt.yscale('log')

plt.figure(2)
fig = plt.plot(x,m_sol, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.legend(loc='best')
plt.title('Mass of a Newtonian Neutron Star, SLy4 H&D')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\oplus}]$")

#%% EoS Comparison
K = 5.3802e9 #non-rel neutrons, cgs units
gamma = 5/3 #non-rel neutrons.
plt.figure(4)
plt.plot(df.loc[:,'rhog/cm3'], pressure_inter(df.loc[:,'rhog/cm3']), color='darkorange', label='Haensel & Douchin')
plt.scatter(df.loc[:,'rhog/cm3'], df.loc[:,'Perg/cm3'], color='black')
plt.plot(df.loc[:,'rhog/cm3'], K*df.loc[:,'rhog/cm3']**gamma, color='darkblue', label='non-rel Fermi Gas')
plt.title("Comparison of EoS of non-interacting, non-relativistic neutron Fermi Gas vs Tabulated from Haensel & Douchin ")
plt.legend(loc='best')
plt.xlabel('Density $[g/cm^3]$')
plt.ylabel('Pressure $[erg/cm^3]$')
plt.yscale('log')
plt.xscale('log')
# The cutoff for rel. neutrons is at ρ=6e15. Τhe biggest density in the tabulated EoS is 4e15. This means I don't have to care
# about them.
#%% Many Neutron Stars.

Central_Density = np.logspace(np.log10(3e14),np.log10(2.8e15),30) # <5e15 to be inside of interpolation range. The Liquid core starts at 1.4e14
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
    plt.title('Pressure of a Neutron Star, SLy4 H&D')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Pressure $[erg/cm^3] $")
    plt.yscale('log')
    plt.ylim(1e26)
    plt.savefig('sly4_pr.pdf', format='pdf')
    
    plt.figure(2)
    fig = plt.plot(x,m_sol, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='best')
    plt.title('Mass of a Neutron Star')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Mass $[M_{\odot}]$")
    plt.savefig('sly4_mr.pdf', format='pdf')
        
plt.figure(3) 
plt.plot(R_tot, M_tot, 'o')
plt.title('Total Masses of NS\'s of varying central densities')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\odot}]$")
plt.savefig('sly4_max.pdf', format='pdf')
#%% Variables for the multi plot
M_hdsly4 = M_tot
R_hdsly4 = R_tot
den_hdsly4 = df['rhog/cm3']
pres_hdsly4 = pressure_inter(den_hdsly4)
text_hdsly4 = "{:.1e}".format(Central_Density[1]) + ' ' + '-' +' ' + "{:.1e}".format(Central_Density[-1])
