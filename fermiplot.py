# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 21:32:17 2021

@author: Konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import math as m

mu=0.5 #eV
kb=8.613e-5 #eV/K

def f(E,T):
    b = (E-mu)/(kb*T)
    Ef=E-mu
    a = 1/(m.exp(b)+1)
    return [a,Ef]

T=np.linspace(50,750,5)
Ex=np.arange(-1,1, 0.001)
Fgraph=np.zeros(len(Ex))
diff=np.zeros(len(Ex))
colarr=['tab:blue','tab:orange','tab:green','tab:red','tab:purple']

for i in range(len(T)):
    for j in range(len(Ex)):
        Fgraph[j] = f(Ex[j], T[i])[0]
        diff[j] = f(Ex[j], T[i])[1]

    ax = plt.subplot(111)
    plt.plot(diff, Fgraph, color=colarr[i], label='T='+str(round(T[i]))+' K')
    plt.grid(which='major',axis='y')
    plt.ylim(-0.1,1.1)
    plt.xlim(-0.35,0.35)
    plt.legend(loc='lower left')
    plt.title('Κατανομή Fermi-Dirac για διάφορες θερμοκρασίες')
    plt.xlabel('$E-E_{Fermi} \hspace{1} [eV]$')
    plt.ylabel('f(E)')
    
    akides_y=np.arange(0,1.25,0.25)
    plt.yticks(ticks=akides_y)
    ind = 2 # use np.where(akides_y==0.5) to find index of value 0.5
    gridlines = ax.yaxis.get_gridlines()
    gridlines[ind].set_color("k")
    gridlines[ind].set_linewidth(1.5)