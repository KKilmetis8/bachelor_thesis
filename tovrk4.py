# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 23:57:23 2021

@author: Konstantinos
"""


import numpy as np
import math as m
import matplotlib.pyplot as plt
#%% Defining the Constants and the ODEs
P_0=1.285 # GeV/fm^3 = E_0
M_0=2.837 # solar masses
R_0=8.378 # km
# Assuming a polytropic EoS, where P=K*ρ^{(ν+1)/ν} and according to wikipedia https://en.wikipedia.org/wiki/Polytrope
# apparently 0.5<ν<1 is a reasonable approximation for a NS
# I'll implement ν=1 for starters, but this is easily modifiable
nu=1
# The Polytropic EoS for ν=1 results in E=c^2 sqrt{P}
# I will accept J. Piekarewicz 2017 rescaling, maybe it's wrong though.
def pressa(x,P,mass):
    a = -0.5 * (m.sqrt(P)+P) *(mass+ 3*x**3*P) / (x**2 *(1-mass/x) )
    return a

def maza(x,P):
    b = 3* x**2 *m.sqrt(P)
    return b
#%% Doing the Thing
# There are 2 problems
# a. how does rk4 work for a system of eqs instead of a solitary one -> maybe i got it?
# b. How on earth do I get the E[i] => Assume a polytropic. 

h=0.01
r = np.arange(0,1,h)
r[0]=0.001 #else the first step divides by 0
nr = r.size
M = np.zeros(nr) 
M[0] = 0 # m(0)=0, No mass at the center
P= np.zeros(nr)
P[nr-1] = 0 # P(r=R)=0, No pressure at the edge. nr-1 is the last element of the array
            # technically, this is already set through initializing with zeroes
            # instead of something else, but I wanted it to be visible
P[0]=P_0 #test

#%% RK4
#RK4 Loop
for i in range(nr - 1):
         P1=h* pressa(r[i],P[i],M[i])
         P2=h* pressa(r[i]+h/2, P[i] + P1/2, M[i])
         P3=h* pressa(r[i]+h/2, P[i] + P2/2, M[i])
         P4=h* pressa(r[i]+h, P[i] + P3, M[i])
         P[i+1]=P[i] + 1/6 * (P1+2*P2+2*P3+P4)
         
         M1=h*maza(r[i],P[i])
         M2=h*maza(r[i]+h/2, P[i])
         M3=h*maza(r[i]+h/2, P[i])
         M4=h*maza(r[i]+h/2, P[i])
         M[i+1]=M[i] + 1/6 * (M1+2*M2+2*M3+M4)


r[:]=R_0*r[:]
M[:]=M_0*M[:]
fig = plt.plot(r,P)
fig = plt.plot(r,M)
plt.xlabel("R $[Km]$")
plt.ylabel("P $[GeV/fm^3]$ - M $[M_{\odot}]$")
plt.title("Attempt at solving the TOV eqs with RK4")
plt.legend(['Pressure','Gravitational Mass'])
