# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 19:21:12 2022

@author: Konstantinos
"""
import matplotlib.pyplot as plt # Plotting
import numpy as np

# For this to work you need
# 1. Run all the solvers for each EoS
# I could automate this, but I don't think that the
# relative time gain is not worth it.
plt.figure(2)
start = 9.4
end = 17.5
plt.plot(R_apr, M_apr, 'o' ,label='APR ' + 'Densities: ' + text_apr, zorder=2)
plt.plot(R_sly4,M_sly4, 's', label='SLy4 '+ 'Densities: '  + text_sly4, zorder=2)
#plt.plot(R_hdsly4,M_hdsly4, 's', label='Haensel & Douchin SLy4 '+ 'Densities: '  + text_hdsly4, zorder=2)
plt.plot(R_ds8, M_ds8, '^', label='DS-8 '+ 'Density range:' + text_ds8, zorder=2)
plt.plot(R_bsk20, M_bsk20,'h', label='BSk20 ' + 'Density range: ' + text_bsk20, zorder=2)

plt.title('NS Mass-Radius Diagram for various EoS')
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\odot}]$")
plt.legend(loc='center right')
plt.xlim(start,end)
plt.grid(True)

#plt.axhline(y=1.97, color='black', linestyle='-')
plt.fill_between( (start,end) , 1.94, 2.01, color='plum' ,zorder=1)
plt.text(x=16, y=1.96, s='PSR J1614-2230', fontweight='bold')
plt.savefig('max_mass.pdf', format='pdf')
#%% EoS Comparison
start2 = 1e15
end2 = 3e15
plt.figure(3)

plt.plot(den_apr, pres_apr, label='APR '+ text_apr)
plt.plot(den_sly4, pres_sly4, label='SLy4 ' + text_sly4)
#plt.scatter(den_hdsly4, pres_hdsly4, label='SLy4' + text_sly4)
plt.plot(den_ds8, pres_ds8, label='DS-8 ' + text_ds8)
plt.plot(den_bsk20, pres_bsk20, label='BSk20 ' + text_ds8)

K = 5.3802e9 #non-rel neutrons.
gamma = 5/3 #non-rel neutrons.
gamma_rel = 4/3 # Relativistic neutrons
K_rel = 4.878e14 # Shapiro, Teukolsky p.28


analytic_pressure = [K*den_apr[i]**gamma for i in range(len(den_apr))]
plt.plot(den_apr, analytic_pressure, color='black', label='non-rel Fermi Gas')
# analytic_pressure_rel = [K_rel*den_cgs[i]**gamma_rel for i in range(len(den_cgs))]

plt.title('EoS Comparison Diagram')
plt.legend(loc='best')
plt.xlabel('Density $[g/cm^3]$')
plt.ylabel('Pressyre $[erg/cm^3]$')
plt.grid(True)
plt.savefig('eoss.pdf', format='pdf')
#%% Sly Only

plt.plot(R_hdsly4, M_hdsly4, 's', label='Haensel & Douchin SLy4 '+ 'Densities: '  + text_hdsly4, zorder=2)
plt.plot(R_sly4, M_sly4, 's', label='SLy4 '+ 'Densities: '  + text_sly4, zorder=2)
plt.xlabel("R $[Km]$")
plt.ylabel("Mass $[M_{\odot}]$")
plt.legend(loc='best')

#%% No CompOSE

plt.plot(R_hdsly4,M_hdsly4, label='Haensel & Douchin SLy4 '+ 'Densities: '  + text_hdsly4, zorder=2)
plt.plot(R_bsk20, M_bsk20, label='BSk20 ' + 'Density range: ' + text_bsk20, zorder=2)
plt.legend(loc='best')
plt.savefig('hdsly4_bsk20.eps', format='eps')

#%% Star 1

# Comparing Pressure of 1.9e15 NS.
plt.figure(4)
plt.plot(r1bsk, p1bsk, label ='BSk20')
plt.plot(r1slyhd, p1slyhd, label='SLy4')
plt.legend(loc='best')
plt.title('Pressure of a NS of Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")
plt.yscale('log')
plt.savefig('P1.eps', format='eps')

# Comparing Density
plt.figure(5)
plt.plot(r1bsk, d1bsk, label ='BSk20')
plt.plot(r1slyhd, d1slyhd, label='SLy4')
plt.legend(loc='best')
plt.title('Density of a NS of Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.xlabel("R $[Km]$")
plt.ylabel("Density $[g/cm^3] $")
plt.yscale('log')
plt.savefig('D1.eps', format='eps')
#%% Star 2

# Comparing Pressure of 6e14 NS.
plt.figure(6)
plt.plot(r2bsk, p2bsk, label ='BSk20')
plt.plot(r2slyhd, p2slyhd, label='SLy4')
plt.legend(loc='best')
plt.title('Pressure of a NS of Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")
plt.yscale('log')
plt.savefig('P2.pdf', format='pdf')

# Comparing Density of 6e14 NS.
plt.figure(7)
plt.plot(r2bsk, d2bsk, label ='BSk20')
plt.plot(r2slyhd, d2slyhd, label='SLy4')
plt.legend(loc='best')
plt.title('Density of a NS of Central Density: ' + den_text + '$\hspace{1} [g/cm^3]$')
plt.xlabel("R $[Km]$")
plt.ylabel("Pressure $[erg/cm^3] $")
plt.yscale('log')
plt.savefig('D2.pdf', format='pdf')

