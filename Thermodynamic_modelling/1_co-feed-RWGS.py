import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os

"""In this file we calculate the equilibrium conversion of CO2 in the reverse water-gas shift 
as a function of temperature and excess hydrogen feed ratio H2:CO2. The RWGS-database is a 
modified version of the GRI-MECH database. It was modified by removing unwanted species in 
order to avoid side reactions, which also speeds up the equilibrate function"""




# define a solution object
RWGS = ct.Solution('RWGS-database.cti')

# Make arrays over which to calculate the CO2 conversion
nH2_range = np.arange(1.0,4.0,0.05)                 # [-] H2/CO2 = nH2
T_range = np.arange(573, 1183, 10)                  # [K]
p = 100000                                          # [Pa] or 1 bar

# make an array to store the results
X_CO2_map = np.zeros((len(nH2_range), len(T_range)))

for i, n_H2 in enumerate(nH2_range):
    for j , T in enumerate(T_range):
        RWGS.X = {"H2": n_H2, "CO2": 1.0}
        RWGS.TP = T, p
        RWGS.equilibrate('TP')
        CO2_conversion = RWGS.X[RWGS.species_index('CO')]/(RWGS.X[RWGS.species_index('CO')]+ RWGS.X[RWGS.species_index('CO2')])
        X_CO2_map[i][j] = CO2_conversion

T_range, nH2_range  = np.meshgrid(T_range, nH2_range)

# Plot the results on a heatmap contour plot
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
plt.tight_layout()
plt.savefig(os.path.join('plots',"figure_1c.png"), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots',"figure_1c.pdf"), bbox_inches='tight')
plt.show()