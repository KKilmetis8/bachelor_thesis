# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 14:40:55 2021

@author: Konstantinos
"""
import numpy as np # Arrays
import matplotlib.pyplot as plt # Plotting
#%% Defining the Constants and ODEs
G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
mu=2
def pressa(r,P,m,gamma,K):
    if r==0:
        return 0 # dP/dr=0 at the core. The pressure stays constant.
    a = - G* P**(1/gamma)*m/(r**2* K**(1/gamma))  # dP/dr = -G M P^(1/Γ) / r^2 K^(1/Γ)
    return a

def maza(r,P,gamma,K):
    if r==0:
        return 0 # dM/dr=0 at the core. The mass stays constant.
    b = 4*np.pi* r**2* P**(1/gamma)/K**(1/gamma)  # 4πr^2 (P/K)^(1/Γ)
    return b

#%% Solving it
h=100000 #works on 100
M=[]
P=[]
x=[]
r=np.arange(0,9e9,h)
nr=r.size
M.append(0) # No mass at the center
i=0
j=0
# RK4 Loop, add the P_c loop.
Central_Density = np.logspace(4,7,20) # Range of different central densities. 
Central_Pressure = np.zeros( len(Central_Density)) # Empty array, same length as densities.
M_tot = []
R_tot = []
for j in range(len(Central_Density)):
    # This checks whether we have relativistic electrons
    # or not, and adjusts the pressure, K, Γ, accordingly
    if Central_Density[j]>=9e5: #cutoff is at 1e6.
        gamma = 4/3 # Relativistic electrons
        K = 1.2435e15/mu**gamma
        P.append(K*Central_Density[j]**gamma)
    else:
        gamma = 5/3 # Non-relativistic electrons
        K = 1.0036e13/mu**gamma
        P.append(K*Central_Density[j]**gamma)
    for i in range(nr-1):
             P1=h* pressa(r[i],P[i],M[i],gamma,K)
             M1=h* maza(r[i],P[i],gamma,K)
             
             P2=h* pressa(r[i]+h/2, P[i] + P1/2, M[i]+ M1/2,gamma,K)
             M2=h* maza(r[i]+h/2, P[i]+ P1/2,gamma,K)
             
             P3=h* pressa(r[i]+h/2, P[i] + P2/2, M[i]+M2/2,gamma,K)
             M3=h* maza(r[i]+h/2, P[i]+P2/2,gamma,K)
             
             P4=h* pressa(r[i]+h, P[i] + P3, M[i]+ M3,gamma,K)
             M4=h* maza(r[i]+h, P[i]+ P3,gamma,K)
             x.append(r[i])
             if P[i] < 1e15:
                 break
             
             P.append(P[i] + 1/6 * (P1+2*P2+2*P3+P4))
             M.append(M[i] + 1/6 * (M1+2*M2+2*M3+M4))
             
             
    # Proper Units, x=X0*x[]
    rp= [n*1e-5 for n in x]      # get Km
    Mp= [n/(1.989*1e33) for n in M] # get solar masses
    M_tot.append(Mp[-1]) # Total mass
    R_tot.append(rp[-1]) # Radius

    plt.figure(1)
    den_text = "{:.1e}".format(Central_Density[j])
    fig = plt.plot(rp,P, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='best')
    plt.title('Pressure of a Newtonian White Dwarf')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Pressure $[erg/cm^3] $")
    plt.yscale('log')

    plt.figure(2)
    fig = plt.plot(rp,Mp, label='Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
    plt.legend(loc='best')
    plt.title('Mass of a Newtonian White Dwarf')
    plt.xlabel("R $[Km]$")
    plt.ylabel("Mass $[M_{\oplus}]$")
    
    # Reset the lists
    M=[]
    M.append(0) # No mass at the center
    P=[]
    x=[]
plt.figure(3) 
plt.plot(R_tot, M_tot, 'o')
plt.title('Total Masses of WD\'s of varying central densities')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\oplus}]$")