# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 20:18:05 2021

@author: Konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt

gamma = 4/3
c = 2.99792458e10 # cgs [cm s^-1]
K = 493483302030614.0
h = 1e25
E = np.arange(1e26,1e28,h)
P = np.zeros(E.size)
for i in range(E.size):
    P[i] = K * (E[i]/c**2)**(gamma)
    
fig = plt.plot(E,P)
plt.title('How $E_0$ impacts $P_0$')
plt.xlabel("Εnergy Density $[erg/cm^3]$")
plt.ylabel("Pressure $[erg/cm^3] $")