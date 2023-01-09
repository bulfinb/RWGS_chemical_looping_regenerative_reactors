import cantera as ct
import numpy as np
import os
import matplotlib.pyplot as plt

T_range = np.arange(300, 1273, 1)

print("RWGS data:")


H2O = ct.Solution("RWGS-database.cti")
H2O.X = {"H2O":1 }

H2 = ct.Solution("RWGS-database.cti")
H2.X = {"H2":1 }

CO = ct.Solution("RWGS-database.cti")
CO.X = {"CO":1 }

CO2 = ct.Solution("RWGS-database.cti")
CO2.X = {"CO2":1 }

# make lists for the change in gibbs free energy and the equilibrium constant
dG = []
K = []
found = 0 # keep track of when the dG becomes negative
for T in T_range:
    H2O.TP = T, 100000
    H_H2O = H2O.enthalpy_mole/1000
    S_H2O = H2O.entropy_mole/1000
    H2.TP = T, 100000
    H_H2 = H2.enthalpy_mole/1000
    S_H2 = H2.entropy_mole/1000
    CO.TP = T, 100000
    H_CO = CO.enthalpy_mole/1000
    S_CO = CO.entropy_mole/1000
    CO2.TP = T, 100000
    H_CO2 = CO2.enthalpy_mole/1000
    S_CO2 = CO2.entropy_mole/1000
    dH_r = H_H2O + H_CO - H_H2 - H_CO2
    dS_r = S_H2O + S_CO - S_H2 - S_CO2
    dG_r = dH_r -T*dS_r
    dG.append(dG_r/1000)
    K_r = np.exp(-dG_r/(8.314*T))
    K.append(K_r)
    if found==0 and dG_r < 0:
        print("\Delta G becomes negative at T = ", T)
        found = 1

# plot the results of dG and K vs. T
xaxislabel = '$T\;\; \mathrm{[K]}$'
yaxislabel = "$\Delta G^\circ \; \; \mathrm{ [kJ \, mol^{-1}] }$"
y2axislabel = '$K$ [-]'
labels = [" $\mathrm{CO_2}$", "$\mathrm{H_2}$"]

fig = plt.figure(figsize=(4.8, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)

axr = ax.twinx()
axr.spines['right'].set_color('C0')
axr.yaxis.label.set_color('C0')
axr.tick_params(axis = 'y', which = 'both', colors='C0')
axr.set_ylabel( y2axislabel )
axr.plot(T_range, K, color = 'C0')
axr.set_yscale('log')

ax.axhline(y=0, color = 'black', lw=0.7)
a, = ax.plot(T_range, dG, ls='-',  lw=1.0, color = 'black' , label=labels[0])
ax.set_xlim(300,1273)
ax.set_ylim(-10, 30)

plt.tight_layout()
plt.savefig(os.path.join('plots','ESI_dG' + '.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'ESI_dG' + '.pdf'), bbox_inches='tight')
plt.show()
