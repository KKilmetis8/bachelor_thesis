# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 21:54:24 2021

@author: Konstantinos
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
#%% Simple RK4
# y' = 5x^2 - y / e^{x+y}, for 0<x<1, h=0.1
# y(0)=1
# Theory gives us, y(x+h)= y(x) + 1/6 ( f1+2f2+2f3+f4 )
# Where F1=h*f(x,y)
#       F2=h*f(x + h/2, y +  F1/2)
#       F3=h*f(x + h/2, y + F2/2)
#       F4=h*f(x+h, y+F3)
h=0.1

def f(x,y):
    y1 = (5*x**2 - y) / m.exp(x+y) 
    return y1

x = np.arange(0,1,h)
nx = x.size
#ny = y0.size
y = np.zeros(nx) #here one would write y=np.zeros((ny,nx)) for more dims
y[0] = 1 # y(0)=1

for i in range(nx - 1):
         F1=h* f(x[i],y[i])
         F2=h* f(x[i]+h/2, y[i] + F1/2)
         F3=h* f(x[i]+h/2, y[i] + F2/2)
         F4=h* f(x[i]+h, y[i] + F3)
         y[i+1]=y[i] + 1/6 * (F1+2*F2+2*F3+F4)

fig = plt.plot(x,y)
