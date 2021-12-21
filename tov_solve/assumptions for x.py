# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:50:32 2021

@author: Konstantinos
"""
import numpy as np

me = 9.1094e-28
mn = 1.674e-24
mb = 1.66057e-24
hbar = 1.0546e-27
h = 6.6261e-27
c = 2.997e10

def compton(m):
    " Returns the compton wavelength in cgs, given the mass"
    return (h/(m*c))

def den_x(m,x):
    " Denisty of a Deg, Cold Fermi Gas for rel. parameter x"
    return (mn*x**3/(3*np.pi**2 *compton(m)**3 ))

print('WD Density, for x=0.1:', "{:.2E}".format(den_x(me,0.1)))
print('WD Density, for x=1:', "{:.2E}".format(den_x(me,1)))
print('WD Density, for x=10:', "{:.2E}".format(den_x(me,10)))
print('---')
print('NS Density, for x=0.1:', "{:.2E}".format(den_x(mn,0.1)))
print('NS Density, for x=1:', "{:.2E}".format(den_x(mn,1)))
print('NS Density, for x=10:', "{:.2E}".format(den_x(mn,10)))

# Conclusion, the assumptions for x>>1 and x<<1 are real bad.
#%% K

def k_rel(m):
    """
    Returns the proportionality constant between Pressure and Density for a fully deg, relativistic fermi gas
    Needs to be divided by μ^Γ. μ: mean molecular weight per electron, Γ:4/3
    """
    return (hbar*c*3**(1/3) *np.pi**(2/3) /(4*mb**(4/3) ))

print('Electron rel K:', "{:.2E}".format(k_rel(me)))
print('Neutron rel K:', "{:.2E}".format(k_rel(mn)))

def k_nonrel(m):
    """
    Returns the proportionality constant between Pressure and Density for a fully deg, ΝΟΝrelativistic fermi gas
    Needs to be divided by μ^Γ. μ: mean molecular weight per electron, Γ:5/3
    """
    return (3**(2/3)*np.pi**(4/3)*hbar**2/(5*m*mb**(5/3)))

print('Electron nonrel K:', "{:.2E}".format(k_nonrel(me)))
print('Neutron nonrel K:', "{:.2E}".format(k_nonrel(mn)))
#%%