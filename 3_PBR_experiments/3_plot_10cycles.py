#!/usr/bin/env python
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

"""This file plots the 10 cycle experiment data. The data logged during the experiment was  pre-processed 
(not in this script) to correct the H2 signal for cross sensitivity to CO and CO2. 
A background signal of 0.15% H2 was also removed from the data.
"""

filename = "10_cycle_experimental_dataset.csv"
collumn_names= ['time', 'T', 'p', 'v_H2', 'v_CO2', 'v_Ar', 'x_H2', 'x_CO2', 'x_CO','x_O2', 'x_CH4',
                'Ar MFC sp', 'H2 MFC sp', 'CO2 MFC sp']
# Units
# time  s
# T     C
# p     bar (above atmospheric, i.e. p = p_amb - p_reactor)
# v_i   standard litres/min
# x_i   %
# MFC sp  - MFC set point standard litres/min

#  load data as pandas arrays and make an instance of the Tga_data object for each file
data = pd.read_csv( filename, skiprows=1, sep=',', names=collumn_names)

# Plot the enthalpy and entropy data
fig = plt.figure(figsize=(10,4.3))
ax = plt.subplot()
ax.set_xlabel('Time [s]')
ax.set_ylabel('$x_\mathrm{i}$ [%]')

# make to additional axes
axrr = ax.twinx()
axrr.spines['right'].set_color('C1')
axrr.yaxis.label.set_color('C1')
axrr.tick_params(axis='y', colors='C1')
axrr.spines.right.set_bounds((800.0, 815.0))
axr = ax.twinx()
axr.yaxis.set_ticks_position('left')
axr.set_ylabel('$v_\mathrm{i}$ [l/min]')
axrr.set_ylabel('$T$ [°C]', color ='C1')


axrr.plot(data['time'], data['T'], lw = 1.0, color = 'C1', label = 'packed bed')

ax.plot(data['time'], data['x_H2'], lw = 1.0, color = 'C2', label = '$x_\mathrm{H_2}$')
ax.plot(data['time'], data['x_CO'], lw = 1.0,color = 'C3', label = '$x_\mathrm{CO}$')
ax.plot(data['time'], data['x_CO2'], lw = 1.0, color = 'C0', label = '$x_\mathrm{CO_2}$')

axr.plot(data['time'], data['v_H2'], lw = 1.0, color='black', ls ='-', label = "$v_\mathrm{H_2}$")
axr.plot(data['time'], data['v_CO2'], lw = 1.0, color='dimgrey', ls ='--', label = "$v_\mathrm{CO_2}$")
axr.plot(data['time'], data['v_Ar'], lw = 1.0, color='grey', ls =':',  label = "$v_\mathrm{Ar}$")

# add zero line
ax.plot([0,23800], [0,0], lw = 1.0,color = 'black', label = '__nolegend__')

# set axes limits and tick paramaters
ax.set_xlim(0,23800)
ax.set_ylim(-2.5,5.3)
axr.set_ylim(0,3.6)
axrr.set_ylim(800,870)
ax.set_yticks([0, 1, 2, 3, 4, 5])
axr.set_yticks([0, 0.3, 0.6, 0.9])
axrr.set_yticks([800, 805, 810, 815])
ax.yaxis.set_label_coords(-0.05,0.65)
axr.yaxis.set_label_coords(-0.07,0.13)
axrr.yaxis.set_label_coords(1.06,0.13)

# add legends
ax.legend(loc='upper left', ncol=1)
axr.legend(bbox_to_anchor=(0.02, 0.086), ncol=1)

plt.tight_layout()
plt.savefig(os.path.join('plots', 'SI_PBR_10_cycles' + '.png'), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'SI_PBR_10_cycles' + '.pdf'), bbox_inches='tight')
plt.show()