import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os

"""In this file we calculate oxygen chemical potential profiles as a function of reaction extent as visualised in 
figure 1 b and e of the manuscript. We also calculate mu_O(delta) for CeO_(2-delta) at 1073 K.
For the co-feed RWGS we visualize the oxygen chemical potential in CO2 and H2O for the two half reactions;
(i) CO2 -> CO + 1/2 O2
(ii) H2 + 1/2 O2 -> H2O
For each half reaction the oxygen chemical potential (or oxygen partial pressure) at a given temperature 
is fixed by the ratio (i) p_CO/p_CO2 and (ii) p_H2/p_H2O. when the oxygen chemical potential of these two half reactions
matches then the reaction has reached equilibrium.
For the countercurrent membrane reactor we calculate chemical potential as a function of exchanged oxygen for both half
reactions in countercurrent streams"""

# set the process condtions
p = 100000 # Pa
T = 1073   # K
nH2 = 1.5  # H2/CO2
R = 8.314  # ideal gas constant

# get oxygen chemical potential at the reference pressure 1 bar
O2 = ct.Solution("RWGS-database.cti")
zeros = np.zeros(len(O2.X))
O2.X = zeros
O2.X = {"O2": 1.0}
O2.TP = T, p
h_O2 = O2.enthalpy_mole/1000    #1/1000 to get J/mol
s_O2 = O2.entropy_mole/1000
mu_O2 = h_O2 - s_O2*T           # O2 chemicalpotential at 1 bar and 1073 K
print("T, p = ", T, p)
print("O2, h, s, mu = ", h_O2, s_O2, mu_O2)


# COFEED figure 2b

# define solution objects for the half reactions
CO2 = ct.Solution('RWGS-database.cti')
H2O = ct.Solution('RWGS-database.cti')

# make lists to store the equilibrium oxygen partial pressure
pO2_CO2 = []
pO2_H2O = []
X_CO2 = np.arange(0, 1.004, 0.004)
for x in X_CO2:
    # iterates through CO2 conversion until partial pressures of the half reactions match
    CO2.X = {'CO': 1, 'O2': 0.5 - x/2}
    H2O.X = {'H2': nH2, 'O2': x/2}
    CO2.TP = T, p
    H2O.TP = T, p
    CO2.equilibrate('TP')
    H2O.equilibrate('TP')
    pO2_1 = CO2.X[CO2.species_index('O2')] * p
    pO2_2 = H2O.X[H2O.species_index('O2')] * p
    pO2_CO2.append(pO2_1)
    pO2_H2O.append(pO2_2)
    if (pO2_2 > pO2_1 or x > 0.999):
        # we slightly overshoot equilibrium
        X_CO2_equilibrium = x
        break

# data for plotting
X = np.arange(0, X_CO2_equilibrium+0.004, 0.004)
mu_O_CO2 = 0.5 * (R * T * np.log(pO2_CO2) + mu_O2) / 1000    # mu_O = 0.5(\mu_O2^\circ + RT ln(p_O2/p^\circ))
mu_O_H2O = 0.5 * (R * T * np.log(pO2_H2O) + mu_O2) / 1000     # mu_O = 0.5(\mu_O2^\circ + RT ln(p_O2/p^\circ))

# plot the results
xaxislabel = '$X\;\; \mathrm{[-]}$'
yaxislabel = "$\mathrm{\mu_{\mathrm{O}}\; \; [kJ \, mol^{-1}]}$"
labels = [" $\mathrm{CO_2}$", "$\mathrm{H_2O}$"]

fig = plt.figure(figsize=(4.0, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
C1 = "black"
C2 = "dimgrey"

a, = ax.plot(X, mu_O_CO2, ls='-',  lw=1.0, color = C1 , label=labels[0])
# add arrows to the lines to indicate progression of the reaction
ax.quiver(X[20], mu_O_CO2[20], 1, -1, color=C1 , width=0.005, headwidth=8, scale =100)
ax.quiver(X[90], mu_O_CO2[90], 1, -0.6, color=C1 , width=0.005, headwidth=10, scale =100)

b, = ax.plot(X, mu_O_H2O, ls='-',  lw=1.0, color = C2, label=labels[1])
# add arrows to the lines to indicate progression of the reaction
ax.quiver(X[20], mu_O_H2O[20], 1, 1, color=C2, width=0.005, headwidth=12, scale =100)
ax.quiver(X[90], mu_O_H2O[90], 1, 0.6, color=C2, width=0.005, headwidth=12, scale =100)

# visualise the transfer of oxygen from CO2 to H2 and from higher to lower chemical potential
a, = ax.plot([X[40], X[40]], [mu_O_CO2[40],mu_O_H2O[40]], ls='--',  lw=1.0, color = 'C3' , label=labels[1])
ax.quiver(X[40], (mu_O_CO2[40]+mu_O_H2O[40])/2.0, 0, -1, color='C3', width=0.005, headwidth=12, scale =100)
ax.text(X[40]-0.05, (mu_O_CO2[40]+mu_O_H2O[40])/2.0 -1, '$\mathrm{O}$', color = 'C3')

# mark the equilibrium
a, = ax.plot([X.max()-0.002, X.max()-0.002], [mu_O_H2O.max(), -300.0], ls='--',  lw=1.0, color = 'grey')
ax.text(X.max()+0.01, -280.0, '$X_\mathrm{eq}$', color = 'grey')

# add some labels to the graph
ax.text(X[10]+0.01, mu_O_CO2[10], labels[0], color = C1)
ax.text(X[10]+0.04, mu_O_H2O[10], labels[1], color = C2)
ax.text(0.65, -220, '$T =$ '+ str(T) + ' $\, \mathrm{K}$' + '\n'
        + '$\mathrm{H_2/CO_2} =$ ' + str(nH2), color ="black")

ax.set_xlim(0,1)
ax.set_ylim(-300,-200)
plt.tight_layout()
plt.savefig(os.path.join('plots','figure_2b.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots','figure_2b.pdf'), bbox_inches='tight')
plt.show()


# COUNTERCURRENT - figure 2e

# make lists to store the equilibrium oxygen partial pressure
pO2_CO2 = []
pO2_H2O = []
X_CO2 = np.arange(0, 1.004, 0.004)
for x in X_CO2:
    # iterates through CO2 conversion until partial pressures of the half reactions match
    CO2.X = {'CO': 1, 'O2': 0.5 - x/2}
    H2O.X = {'H2': nH2, 'O2': x/2}
    CO2.TP = T, p
    H2O.TP = T, p
    CO2.equilibrate('TP')
    H2O.equilibrate('TP')
    pO2_1 = CO2.X[CO2.species_index('O2')] * p
    pO2_2 = H2O.X[H2O.species_index('O2')] * p
    pO2_CO2.append(pO2_1)
    pO2_H2O.insert(0, pO2_2)
    if np.any(np.asarray(pO2_CO2)-np.asarray(pO2_H2O) <= 0) or x > 0.999:
        # check if we reach equilibrium.
        # for the given conditions we reach full conversion
        X_CO2_equilibrium = x
        break

# data for plotting
X = np.arange(0, X_CO2_equilibrium+0.004, 0.004)
mu_O_CO2 = 0.5 * (R * T * np.log(pO2_CO2) + mu_O2) / 1000    # mu_O = 0.5(\mu_O2^\circ + RT ln(p_O2/p^\circ))
mu_O_H2O = 0.5 * (R * T * np.log(pO2_H2O) + mu_O2) / 1000     # mu_O = 0.5(\mu_O2^\circ + RT ln(p_O2/p^\circ))


# plot the results
xaxislabel = '$X\;\; \mathrm{[-]}$'
yaxislabel = "$\mathrm{\mu_{\mathrm{O}}\; \; [kJ \, mol^{-1}]}$"
labels = [" $\mathrm{CO_2}$", "$\mathrm{H_2}$"]

fig = plt.figure(figsize=(4.0, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
# colors
C1 = "black"
C2 = "dimgrey"

# shade between the lines to indicate range of chemical potentials that need to be stored in the packed bed
ax.fill_between(X, mu_O_H2O, mu_O_CO2, color = 'lightgrey', alpha=0.5)

a, = ax.plot(X, mu_O_CO2, ls='-',  lw=1.0, color = C1 , label=labels[0])
# add arrows to the lines to indicate progression of the reaction
ax.quiver(X[20], mu_O_CO2[20], 1, -1, color=C1 , width=0.005, headwidth=8, scale =100)
ax.quiver(X[90], mu_O_CO2[90], 1, -0.6, color=C1 , width=0.005, headwidth=10, scale =100)
ax.quiver(X[160], mu_O_CO2[160], 1, -0.6, color=C1 , width=0.005, headwidth=10, scale =100)
ax.quiver(X[230], mu_O_CO2[230], 1, -1, color=C1 , width=0.005, headwidth=8, scale =100)

b, = ax.plot(X, mu_O_H2O, ls='-',  lw=1.0, color = C2, label=labels[1])
# add arrows to the lines to indicate progression of the reaction
ax.quiver(X[20], mu_O_H2O[20], -1, 0.4, color=C2, width=0.005, headwidth=12, scale =100)
ax.quiver(X[90], mu_O_H2O[90], -1, 0.4, color=C2, width=0.005, headwidth=12, scale =100)
ax.quiver(X[160], mu_O_H2O[160], -1, 0.5, color=C2, width=0.005, headwidth=12, scale =100)
ax.quiver(X[230], mu_O_H2O[230], -1, 1, color=C2, width=0.005, headwidth=8, scale =100)

# visualise the transfer of oxygen from CO2 to H2 and from higher to lower chemical potential
a, = ax.plot([X[40], X[40]], [mu_O_CO2[40],mu_O_H2O[40]], ls='--',  lw=1.0, color = 'C3' , label=labels[1])
ax.quiver(X[40], (mu_O_CO2[40]+mu_O_H2O[40])/2.0, 0, -1, color='C3', width=0.005, headwidth=12, scale =100)
ax.text(X[40]-0.05, (mu_O_CO2[40]+mu_O_H2O[40])/2.0 -1, '$\mathrm{O}$', color = 'C3')

# add labels
ax.text(X[10]+0.01, mu_O_CO2[10], labels[0], color = C1)
ax.text(X[240]-0.09, mu_O_H2O[240], labels[1], color = C2)
ax.text(0.65, -220, '$T =$ '+ str(T) + ' $\, \mathrm{K}$' + '\n'
        + '$\mathrm{H_2/CO_2} =$ ' + str(nH2), color ="black")

ax.set_xlim(0,1)
ax.set_ylim(-300,-200)
plt.tight_layout()
plt.savefig(os.path.join('plots','figure_2e.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots','figure_2e.pdf'), bbox_inches='tight')
plt.show()


# Also calcultae and plot mu_O(delta) for CeO2
# CeO2 Thermodymics https://doi.org/10.1039/C6CP03158G  Eq 22 and 23 and table 1
dH = 430000   # [J mol^-1]
ds_th = 165   # J mol^-1 K^-1
d_m = 0.35    # -
n = 2.31      # -

def pO2(delta, T):
    """equilibrium oxygen partial pressure for CeO_2-delta at temperature T"""
    return 100000 * np.exp(2 * ds_th / R) * np.exp(-2 * dH / (R * T))  * ((d_m - delta) / delta) ** (2 * n)

def mu_O(delta, T):
    """equilibrium mu_O [kJ/mol]"""
    return 0.5 * ( mu_O2 + R * T * np.log( pO2(delta, T) /100000) ) / 1000.0


delta = np.arange(0.00001, 0.34101, 0.001)
mu_O_delta_1073 = mu_O(delta, 1073.0)



# plot mu_O vs delta at 1073 K
xaxislabel = '$\delta\;\; \mathrm{[-]}$'
yaxislabel = "$\mathrm{\mu_{\mathrm{O}}\; \; [kJ \, mol^{-1}]}$"

fig = plt.figure(figsize=(4.0, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
ax.plot(delta, mu_O_delta_1073, ls='-', color = 'black',  lw=1.0, label='1073')
ax.text(0.13, -250, '$T = 1073 \, \mathrm{K}$' + '\n' + '$\mathrm{CeO}_{2-\delta}$ ', color ="black")
ax.set_xlim(-0.001,0.2)
ax.set_ylim(-400,-200)
plt.tight_layout()
plt.savefig(os.path.join('plots','mu_O_vs_delta_CeO2.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots','mu_O_vs_delta_CeO2.pdf'), bbox_inches='tight')
plt.show()
