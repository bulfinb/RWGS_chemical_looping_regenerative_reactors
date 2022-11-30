import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os

"""In this file we calculate the thermodynamic limit on CO2 conversion extent for the reverse 
water-gas shift process performed in a countercurrent oxygen permeable membrane reactor. 
The calculation follows the methodology developed in previous work DOI:10.1039/C8CP07077F"""


def maximum_CO2_conversion(T, p, nH2):
    """Determine the maximum conversion of CO2 for countercurrent flows of CO2 and H2 with
     flow rates of F_CO2 = 1, F_H2 = nH2, at temperature T [K], and presure p [Pa].
     kappa is the exchange coefficient which corresponds to the amount of exchanged O2.
     kappa = 0 means no reaction, kappa = 0.5 complete reaction.
     Paramaters
     -T   (temperature [K])
     -p   (pressure [Pa])
     -nH2 (hydrogen excess n_H2/n_CO2 [-])
     returns
     -X_CO2 (conversion extent of CO2 [-])"""
    # define the two flows as canter solution objects using cantera's gri30 database
    flow1 = ct.Solution('RWGS-database.cti')
    flow2 = ct.Solution('RWGS-database.cti')
    # Make arrays of the pO2(kappa) in each flow (note \mu_O2(p_O2))
    pO2_flow2 = []
    pO2_flow1 = []
    kappa_max = 0
    # complete conversion of CO2 to CO gives kappa = 0.5, kappa in range 0-0.5
    kappa_range = np.arange(0, 0.502, 0.002)
    for kappa in kappa_range:
        flow1.X = {'CO': 1, 'O2': 0.5 - kappa}
        flow2.X = {'H2': nH2, 'O2': kappa}
        flow1.TP = T, p
        flow2.TP = T, p
        flow1.equilibrate('TP')
        flow2.equilibrate('TP')
        pO2_f1 = flow1.X[flow1.species_index('O2')] * p
        pO2_f2 = flow2.X[flow2.species_index('O2')] * p
        # add the new pO2 values to the arrays
        pO2_flow1.append(pO2_f1)
        # for the second flow add pO2 value to the start of array to reverse kappa for countercurrent
        pO2_flow2.insert(0, pO2_f2)
        # if the arrays meet at any point, then the chemical potentials of oxygen are equal at that point
        # and we have reached the maximum exchange extent kappa_max
        if np.any(np.asarray(pO2_flow1) - np.asarray(pO2_flow2) <= 0) or kappa > 0.499:
            X_CO2 = 2 * kappa
            break
    return X_CO2


p = 100000 # Pa

# Make arrays over which to calculate the CO2 conversion
nH2_range = np.arange(1.0,4.0,0.05)                 # [-] H2/CO2 = nH2
T_range = np.arange(573, 1183, 10)                  # [K]
p = 100000                                          # [Pa] or 1 bar

# make an array to store the results
X_CO2_map = np.zeros((len(nH2_range), len(T_range)))
# run the calculations
for i, n_H2 in enumerate(nH2_range):
    print("Progress: ", round(i*100 / len(nH2_range),3), " %")
    for j, T in enumerate(T_range):
        CO2_conversion = maximum_CO2_conversion(T, p, n_H2)
        X_CO2_map[i][j] = CO2_conversion


T_range, nH2_range  = np.meshgrid(T_range, nH2_range)

fig = plt.figure(figsize=(4.2, 3.2), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2/CO_2}$ [-]')
im = ax0.pcolormesh(T_range, nH2_range, X_CO2_map, rasterized=True, vmin =0.1, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
ct = ax0.contour(T_range, nH2_range, X_CO2_map, levels = levels, colors = 'black')
ax0.clabel(ct, fmt='%1.1f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$X_\mathrm{CO_2}$')
#plt.title("")
#plt.plot(xCO2i, xCOi)
plt.tight_layout()
plt.savefig(os.path.join('plots',"figure_1f.png"), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots',"figure_1f.pdf"), bbox_inches='tight')
plt.show()

# Calculation is slow. Save output for faster plot modification
np.savetxt('RWGS_Countercurrent_XCO2_map.csv', X_CO2_map, delimiter=",")