# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 13:29:26 2021

@author: Konstantinos
"""
# We'll solve the Newtonian hydrostatic equilibrium for a WD

import numpy as np
import math as m
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
#%% Defining the Constants and ODEs

# The Electron degeneracy pressure, for relativistic electrons
# is given by Pe= (3π^2)^1/3 / 4 * hbar * c *(1/m_bar* μ)^4/3 * ρ^4/3
# P=Kρ^4/3
# where m_bar=1.66057 10^{-24} is the mean rest mass per baryon
# and μ is the mean molecular mass per electron
# these constants are computed to the 1.2435e15 value.

mu = 2
gamma = 4/3
K = 1.2435e15/mu**gamma # cgs units. Masti p.9 & Balberg, Shapiro p.5
G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
c = 2.99792458e10 # cgs [cm s^-1]

# We derive E_0= Mc^2 / (4/3)*π*R^3 = 1e27 erg/cm^3
# For M=0.5 M_sol and R=1e4 Km, which are typical values for a WD

E0 = 1.989e33*c**2 / (4*m.pi/3 * (1e4*1e5)**3) #1e4 Km and 1e5 to convert to cm

# It only now occured to me, that I could make python do
# the calculation, instead of doing it by hand. Well, at least
# I'm sure the value is good. This relates to P0 since, P=K (E/c^2)^Γ

P0 = K * (E0/c**2)**(gamma)

# Using P0, we may rescale our eqs.

M0 = m.sqrt(P0**(4*gamma+3) / (4*m.pi* G**3* K**(4*gamma) ) )
R0 = G*K**(gamma)*M0 / P0**(gamma+1)

def pressa(x,P,m):
    if x==0:
        return 0 # dP/dr=0 at the core. The pressure stays constant.
                 # Maybe this is a bad way to implement this...
    a = - G* K**(gamma)*m/(x**2* P**gamma) 
    return a

def maza(x,P):
    if x==0:
        return 0 # dM/dr=0 at the core. The mass stays constant.
    b = 4*m.pi* K**(gamma)*x**2 /P**gamma 
    return b

#%% Initial Conditions

# In the for loop, we should iterate -and thus, append- until P<0
# For x=0 => M=0 dM/dr=0
# For x=0 => P=P_c dP/dr=0

h=1
P_c=1e14 # P=P0 * Pbar and we use Pbar |
M=[]
P=[]
r=np.arange(0,1e9,10000)
nr=r.size
M.append(0) # No mass at the center
P.append(P_c) # Test a 100 values from 1e-24 to 1e-26


i=0
# RK4 Loop, add the P_c loop.
for i in range(nr-1):
         P1=h* pressa(r[i],P[i],M[i])
         M1=h* maza(r[i],P[i])
         
         P2=h* pressa(r[i]+h/2, P[i] + P1/2, M[i]+ M1/2)
         M2=h* maza(r[i]+h/2, P[i]+ P1/2)
         
         P3=h* pressa(r[i]+h/2, P[i] + P2/2, M[i]+M2/2)
         M3=h* maza(r[i]+h/2, P[i]+P2/2)
         
         P4=h* pressa(r[i]+h, P[i] + P3, M[i]+ M3)
         M4=h* maza(r[i]+h, P[i]+ P3)
         
         P.append(P[i] + 1/6 * (P1+2*P2+2*P3+P4))
         M.append(M[i] + 1/6 * (M1+2*M2+2*M3+M4))

# Proper Units, x=X0*x[]
rp= [j*1e-5 for j in r]      # get Km
Mp= [j/(1.989*1e33) for j in M] # get solar masses


# Plotting
plt.figure(1)
fig = plt.plot(rp,P)
plt.legend(['Pressure'])
plt.title('Pressure of a Newtonian White Dwarf')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")

plt.figure(2)
fig = plt.plot(rp,Mp, color='darkorange')
plt.legend(['Mass'])
plt.title('Mass of a Newtonian White Dwarf')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\oplus}]$")
#%% Using scipy

def dSdx(r, S):
    p, m = S
    mu = 2
    gamma = 4/3
    K = 1.2435e15/mu**gamma
    G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
    # print(r,S,-m*G* K**(gamma)/(r**2 * p**gamma), 4* K**(gamma)* 3.14*r**2 /p**gamma )
    if r==0:
        return [0,0]
    else:
        return [-m*G* K**(gamma)/(r**2 * p**gamma), 4* K**(gamma)* np.pi *r**2 /p**gamma]

S_0 = [1e17,0]
dt=100000
sol = solve_ivp(fun=dSdx, t_span=(0,1e19), y0=S_0, max_step=dt, method='RK45')
p_sol = sol.y[0] # Rescaling goes here
m_sol = sol.y[1]/(1.989*1e33) # Solar Masses
x=sol.t*1e-5 # Km

# Plotting
plt.figure(1)
fig = plt.plot(x,p_sol)
plt.legend(['Pressure'])
plt.title('Pressure of a Newtonian White Dwarf')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")

plt.figure(2)
fig = plt.plot(x,m_sol, color='darkorange')
plt.legend(['Mass'])
plt.title('Mass of a Newtonian White Dwarf')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\oplus}]$")
