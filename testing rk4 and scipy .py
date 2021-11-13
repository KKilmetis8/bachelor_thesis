# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:03:17 2021

@author: Konstantinos
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt

from scipy.integrate import odeint

# We'll solve the Lotka-Volterra eqs in order to compare my RK4 solver with scipy's. 
# X is the Prey like rabbits, or treats. Y is the Predators like this prime example of
# a hunter: https://imgur.com/a/0go8yIE
#%% Scipy

def dSdx(t, S):
    x, y = S
    a = 1 # Constants
    b = 0.4
    c = 0.4
    d = 0.1
    return [a*x-b*x*y, 
            d*x*y-c*y] # Lotka - Volterra 

S_0 = [20,5] # We start with 20 Prey and 5 Predators
t = np.arange(0,119,0.1)
sol =  odeint(dSdx, t=t, y0=S_0, tfirst=True)
x_sol = sol.T[0] # Rescaling goes here
y_sol = sol.T[1] # Solar Masses

plt.plot(t,x_sol,label='Prey')
plt.plot(t,y_sol, label='Predator')
plt.legend(loc = 'best')
plt.title("Lotka - Volterra by Scipy ")
plt.xlabel('Time in, i dont know, years?')
plt.ylabel('Population')
#%% My RK4
def Prey(t,x,y):
    a = 1 # Constants
    b = 0.4
    dxdt = a*x-b*x*y
    return dxdt

def Predator(t,x,y):
    c = 0.4
    d = 0.1
    dydt = d*x*y-c*y
    return dydt

# Initial Conditions
h=0.1
Treats=[]
Luigi=[]
t = np.arange(0,119,0.1)
nt=t.size
Treats.append(20)  # Treats is the prey
Luigi.append(5) # Luigi is, of course, the predator

i=0
# RK4 Loop, add the P_c loop.
for i in range(nt-1):
         T1=h* Prey(t[i],Treats[i],Luigi[i])
         L1=h* Predator(t[i],Treats[i],Luigi[i])
         
         T2=h* Prey(t[i]+h/2, Treats[i] + T1/2, Luigi[i] + L1/2)
         L2=h* Predator(t[i]+h/2, Treats[i]+ T1/2, Luigi[i] + L1/2)
         
         T3=h* Prey(t[i]+h/2, Treats[i] + T2/2, Luigi[i]+L2/2)
         L3=h* Predator(t[i]+h/2, Treats[i] + T2/2, Luigi[i]+L2/2)
         
         T4=h* Prey(t[i]+h, Treats[i] + T3, Luigi[i]+ L3)
         L4=h* Predator(t[i]+h, Treats[i] + T3, Luigi[i]+ L3)
         
         Treats.append(Treats[i] + 1/6 * (T1+2*T2+2*T3+T4))
         Luigi.append(Luigi[i] + 1/6 * (L1+2*L2+2*L3+L4))
         
plt.plot(t,Treats,label='Prey')
plt.plot(t,Luigi, label='Predator')
plt.legend(loc = 'best')
plt.title("Lotka - Volterra by my own RK4 ")
plt.xlabel('Time in, i dont know, years?')
plt.ylabel('Population')
#%% Residuals
Res_prey =  np.absolute(x_sol-Treats)
Res_predator = np.absolute(y_sol-Luigi)
plt.plot(t,Res_prey,label='Prey')
plt.plot(t,Res_predator, label='Predator')
plt.legend(loc = 'best')
plt.title("Lotka - Volterra Residuals")
plt.xlabel('Time in, i dont know, years?')
plt.ylabel('Population')