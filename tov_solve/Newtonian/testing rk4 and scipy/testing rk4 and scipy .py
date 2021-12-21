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
t = np.arange(0,1000,0.1)
sol =  odeint(dSdx, t=t, y0=S_0, tfirst=True)
x_sol = sol.T[0] 
y_sol = sol.T[1] 

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
#plt.plot(t,Res_predator, label='Predator')
plt.legend(loc = 'best')
plt.title("Lotka - Volterra Residuals")
plt.xlabel('Time in, i dont know, years?')
plt.ylabel('Population')

#%% Looking at the peaks

from scipy.signal import argrelextrema # need it to find local maxima
treats_np=np.array(Treats) # Treats is a list, we need it to be an array

my_max = treats_np[argrelextrema( treats_np, np.greater )[0]] # The argrelextrema outputs the indices of the extrema, the x[argrelextrema(x,np.greater)] outputs the values
scipy_max = x_sol[argrelextrema( x_sol, np.greater )[0]]

# Scatter
res_max = abs(my_max - scipy_max) # Let's have a nice coolwarm graph

fig, ax = plt.subplots()
plt.figure(1)
ax.scatter(my_max, scipy_max, s=25, edgecolors= "grey", c=res_max, cmap=plt.cm.coolwarm_r, zorder=10) #_r at the end reverses the thing
plt.title("Lotka - Volterra Prey Peaks")
plt.xlabel('My Prey Peaks')
plt.ylabel('Scipy Prey Peaks')
plt.plot(my_max,my_max,'k-',alpha=0.75, zorder=0)
peak_res = abs(my_max - scipy_max)

# How do peak_diffs evolve with time?
fig2, ax2 = plt.subplots()
plt.figure(2)
fig2 = plt.plot(np.arange(0,len(res_max),1), res_max)
plt.title("Residuals of Prey Peaks of My RK4 & Scipy")
plt.xlabel('Peak Index [$\propto$ to time]')
plt.ylabel('Peak Residuals')

# It is actually periodical, right?
my_peak_diff = []
scipy_peak_diff = []
for i in np.arange(0,len(res_max)-1,1):
    my_peak_diff.append(my_max[i+1] - my_max[i])
    scipy_peak_diff.append(scipy_max[i+1] - scipy_max[i])

plt.figure(3)
plt.plot(np.arange(0,len(res_max)-1,1), my_peak_diff, label='My RK4')
plt.plot(np.arange(0,len(res_max)-1,1), scipy_peak_diff, label='Scipy')
plt.legend(loc = 'best')
plt.title("Peak Differences for each method ")
plt.xlabel('Peak Index [$\propto$ to time]')
plt.ylabel('Peak diff')