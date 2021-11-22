# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 14:40:55 2021

@author: Konstantinos
"""

#%% Defining the Constants and ODEs

# The Electron degeneracy pressure, for relativistic electrons
# is given by Pe= (3π^2)^1/3 / 4 * hbar * c *(1/m_bar* μ)^4/3 * ρ^4/3
# P=Kρ^4/3
# where m_bar=1.66057 10^{-24} is the mean rest mass per baryon
# and μ is the mean molecular mass per electron
# these constants are computed to the 1.2435e15 value.

mu = 2 # WDs are neutral so you got half an electron per baryon.
gamma = 4/3 # Relativistic electrons. For non-rel it should be 5/3
K = 1.2435e15/mu**gamma # cgs units. Masti p.9 & Balberg, Shapiro p.5
G = 6.67430e-8 # cgs [cm^3 g^-1 s^-2]
c = 2.99792458e10 # cgs [cm s^-1]

def pressa(r,P,m):
    if x==0:
        return 0 # dP/dr=0 at the core. The pressure stays constant.
    a = - G* P**(1/gamma)*m/(r**2* K**(1/gamma))  # dP/dr = -G M P^(1/Γ) / r^2 K^(1/Γ)
    return a

def maza(r,P):
    if x==0:
        return 0 # dM/dr=0 at the core. The mass stays constant.
    b = 4*np.pi* r**2* P**(1/gamma)/K**(1/gamma)  # 4πr^2 (P/K)^(1/Γ)
    return b

#%% Initial Conditions
h=10000 #works on 100
P_c=1e24 # P=P0
M=[]
P=[]
r=np.arange(0,1e9,h)
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
         if P[i]<1:
             print('flag!')
             break

# Proper Units, x=X0*x[]
rp= [j*1e-5 for j in r]      # get Km
Mp= [j/(1.989*1e33) for j in M] # get solar masses


# Plotting
plt.figure(1)
fig = plt.scatter(rp,P)
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