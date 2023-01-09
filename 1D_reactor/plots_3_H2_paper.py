#!/usr/bin/env python
import pandas as pd
import os
import matplotlib.pyplot as plt


filename1 = os.path.join("exported_data","RWGS_demo_CeO2_cycle2_H2-flow_d_vs_x.txt")
filename2 = os.path.join("exported_data","RWGS_demo_CeO2_cycle2_H2-flow_xH2O_vs_x.txt")
filename3 = os.path.join("exported_data","RWGS_demo_CeO2_cycle2_H2-flow_deq-d_vs_x.txt")
#correction = "ExpDat_SA_oxygen_correction.txt"
collumn_names= ['x0', 't = 0 s',  '80', '160', '240', '320', '400', '446.7']
y_cols =  ['t = 0 s',  '80', '160', '240', '320', '400', '446.7']
times = [0, 80, 160, 240, 320, 400, 446.7]

data1 = pd.read_csv(filename1, delim_whitespace = True, names=collumn_names, skiprows =8)
data2 = pd.read_csv(filename2, delim_whitespace = True, names=collumn_names, skiprows =8)
data3 = pd.read_csv(filename3, delim_whitespace = True, names=collumn_names, skiprows =8)


filename4 = os.path.join("exported_data","RWGS_demo_CeO2_cycle2_H2-flow_xH2O_out_vs_t.txt")
data4 = pd.read_csv(filename4, delim_whitespace = True, names=['t','X_H2O'], skiprows =8)

filename5 = os.path.join("exported_data", 'X_H2O_vs_delta.txt')
data5 = pd.read_csv(filename5, delim_whitespace = True, names=['delta','X_H2O'], skiprows =8)

# plot the results of conversion vs. T
xaxislabel = 'Position $x$ [m]'
yaxislabel = 'Non-stoicheometry $\delta$ [-]'
#label_x = [500,500,500,800]
#offsets = [(10,0.07),(10,-0.07),(20,-0.07),(0,-0.07)]
#lines = []

fig = plt.figure(figsize=(4.2, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 0.075)
ax.fill_between(data1['x0'].array, data1['t = 0 s'], data1[y_cols[-1]], color = 'lightgrey', alpha=0.7)
for i, y_col in enumerate(y_cols):
    ax.plot(data1['x0'].array, data1[y_col], lw=1.4, label = y_col)
ax.text(0.33, 0.05, '$\mathrm{H_2}$ flow $\\rightarrow$ ', color = 'dimgrey')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join('plots', 'figure_3b.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_3b.pdf'), bbox_inches='tight')
plt.show()

# plot the results of conversion vs. T
xaxislabel = 'Position $x$ [m]'
yaxislabel = '$x_\mathrm{H_2O}$ [-]'


fig = plt.figure(figsize=(4.2, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.05)
for i, y_col in enumerate(y_cols):
    ax.plot(data2['x0'].array, data2[y_col], lw=1.4, label = y_col)
#plt.legend()
plt.tight_layout()
plt.savefig(os.path.join('plots', 'ESI_mod_cycle2_H2-flow_xH2O_vs_x.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'ESI_mod_cycle2_H2-flow_xH2O_vs_x.pdf'), bbox_inches='tight')
plt.show()

# plot the results of conversion vs. T
xaxislabel = 'Position $x$ [m]'
yaxislabel = '$\delta_\mathrm{eq}(C_\mathrm{H_2O})- \delta$ [-]'


fig = plt.figure(figsize=(4.2, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
ax.set_xlim(0, 1.0)
ax.set_ylim(-0.0005, 0.0022)
for i, y_col in enumerate(y_cols):
    ax.plot(data3['x0'].array, data3[y_col], lw=1.4, label = y_col)
#plt.legend()
ax.hlines(y=0.0, xmin=0, xmax=1.0, linewidth=1.0, color='black', ls = '--')
plt.tight_layout()
plt.savefig(os.path.join('plots', 'ESI_mod_cycle2_H2-flow_deq-d_vs_x.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'ESI_mod_cycle2_H2-flow_deq-d_vs_x.pdf'), bbox_inches='tight')
plt.show()


xaxislabel = 'Time [s]'
yaxislabel = '$x_\mathrm{H_2O, \, out}$ [-]'
#label_x = [500,500,500,800]
#offsets = [(10,0.07),(10,-0.07),(20,-0.07),(0,-0.07)]
#lines = []

fig = plt.figure(figsize=(4.0, 3.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
ax.set_xlim(0, times[-1]+5)
ax.set_ylim(0, 1.05)
ax.text(150, 0.5, 'Outlet $x=1$', color = 'dimgrey')
ax.plot(data4['t'], data4['X_H2O'], color = 'dimgrey', lw=1.4)
#add markers at the time stamps
for i, y_col in enumerate(y_cols):
    data4_cut = data4[data4['t']<times[i]]
    if i == 0:
        ax.plot(0, 1.0, lw=0.0, marker='s', color = 'dimgrey')
    else:
        ax.plot(data4_cut['t'].values[-1], data4_cut['X_H2O'].values[-1],
                 color = 'dimgrey', lw = 0.0, marker = 's')
#plt.legend()
plt.tight_layout()
plt.savefig(os.path.join('plots', 'figure_3c.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_3c.pdf'), bbox_inches='tight')
plt.show()

xaxislabel = '$\delta$ [-]'
yaxislabel = '$x_\mathrm{H_2O}$ [-]'

fig = plt.figure(figsize=(3.2, 2.2), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel)
ax.set_ylabel(yaxislabel)
ax.set_xlim(0, 0.1)
ax.set_ylim(0, 1.05)
ax.text(0.065, 0.85, '$T= 1073 \, \mathrm{K}$', color = 'black')
ax.text(0.025, 0.3, '$x_\mathrm{H_2O, \, eq.}(\delta)$ - CeO$_{2-\delta}$', color = 'dimgrey')
ax.plot(data5['delta'], data5['X_H2O'], color = 'dimgrey', lw=1.4)
#add markers at the time stamps
#plt.legend()
plt.tight_layout()
plt.savefig(os.path.join('plots', 'figure_3a.png'), dpi=400, bbox_inches='tight')
plt.show()

