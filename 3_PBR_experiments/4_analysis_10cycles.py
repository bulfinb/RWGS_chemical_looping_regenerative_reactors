import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cantera as ct
from scipy import integrate
import os

print('Counter Current, 10 cycles')

# RELATIVE ERROR IN THE MASS BALANCE
# 2 % error in bottle concentration. 2 % error in gas analysis signals, and 1 % error in the mass flow.
# The errors in the gas-analysis are used twice in each mass balance.
# This neglects the numerical errors introduced during integration of the signals.
relative_error_MB = np.sqrt(2 * (0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))
print("relative error mass balance:", relative_error_MB)
# relative error in n_f outflows is smaller as we only use on gas signal
relative_error_nf = np.sqrt((0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))

# Bottle concentrations and molar volumes
x_CO2_Bottle = 0.0501           # rest Ar => Accuracy 2%
x_H2_Bottle = 0.0486            # rest Ar => Accuracy 2%
V_0 = 22.414                    # Molar Volume [dm^3/mol] = [l/mol]



# Import Data
filename = "10_cycle_experimental_dataset.csv"
collumn_names= ['time', 'T', 'p', 'v_H2', 'v_CO2', 'v_Ar', 'x_H2', 'x_CO2', 'x_CO','x_O2', 'x_CH4', 'Ar MFC sp', 'H2 MFC sp', 'CO2 MFC sp']
# Units
# time  s
# T     C
# p     bar (above atmospheric, i.e. p = p_amb - p_reactor)
# v_i   standard litres/min
# x_i   %
# MFC sp  - MFC set point standard litres/min

#  load data as pandas arrays and make an instance of the Tga_data object for each file
data = pd.read_csv( filename, skiprows=1, sep=',', names=collumn_names)

# convert measured % to mole fractions
data['x_H2'] = data['x_H2']*0.01
data['x_CO'] = data['x_CO']*0.01
data['x_CO2'] = data['x_CO2']*0.01
data['x_CH4'] = data['x_CH4']*0.01
data['x_O2'] = data['x_O2']*0.01

#Check CH4 and O2 are zero
print("CH4, mean x, max x =", data['x_CH4'].mean(),  data['x_CH4'].max())
print("O2, mean x, max x =", data['x_O2'].mean(),  data['x_O2'].max())

# Calculate molar flows [mol/s]
# inflow
data['F_H2_0'] = (data['v_H2']*x_H2_Bottle)/(V_0*60)
data['F_H2_Ar_0'] = data['v_H2']*(1 - x_H2_Bottle)/V_0/60
data['F_CO2_0'] = data['v_CO2']*x_CO2_Bottle/V_0/60
data['F_CO2_Ar_0'] = data['v_CO2']*(1 - x_CO2_Bottle)/V_0/60
data['F_Ar_0'] = data['v_Ar']/V_0/60
data['F_total_0'] = (data['v_H2'] + data['v_CO2'] + data['v_Ar'])/V_0/60
data['F_Ar_total_0'] = data['F_H2_Ar_0']+data['F_CO2_Ar_0']+data['F_Ar_0']

# outflow, for gas analysis mass balance, H2O has been condensed out which we need to account for
data['F_H2_out'] = data['x_H2']*(data['F_Ar_total_0']+data['F_CO2_0'])/(1-data['x_H2'])
data['F_H2O_out'] = data['F_H2_0'] - data['F_H2_out']
# total outflow minus H2O
data['F_total_out'] = data['F_H2_out'] + data['F_CO2_0']  + data['F_Ar_total_0']
data['F_CO2_out'] = data['x_CO2']*data['F_total_out']
data['F_CO_out'] = data['x_CO']*data['F_total_out']


# Analyse the chemical looping cycles by breaking up the data at the switching points between the steps
# where the Ar flow rate set point was 0.1 slm
timeswitches = []
timeswitches.append(0)
time_skip = 10*60  # skip ten minutes to save time
safety = 60        # Safety after switching

i = 0
while i < len(data['Ar MFC sp']):
    Ar_flow = data['Ar MFC sp'].iloc[i]
    if Ar_flow == 0.1:
        timeswitches.append(data['time'].iloc[i])
        i += time_skip  # skip 10min
    else:
        i += 5  # skip 5sec to next value


print("timeswitches ", timeswitches)
# first cycle is reduction (H2) => first timestamp is end of cycle

# store the data for each cycle steps as pandas dataframes in lists
frames = []
frames_reduction = []
frames_oxidation = []

j = 0
while j < (len(timeswitches)-1):
    # loop to slice up data into sections of oxidation and reduction store in frames
    start_time = timeswitches[j]
    end_time = timeswitches[j+1]
    df = data[data['time'] > start_time]
    df = df[df['time'] < (end_time + safety)]
    frames.append(df)
    j += 1

# use the set flow rates to separate oxidation from reduction
for frame in frames:
    if frame['H2 MFC sp'].mean() > 0.1:
        frames_reduction.append(frame)
    if frame['CO2 MFC sp'].mean() > 0.1:
        frames_oxidation.append(frame)



# Make lists to store the results
cycle_number = [1,2,3,4,5,6,7,8,9,10]
n_CO2_in = []       # [mol]
n_CO_out = []       # [mol]
C_balance = []      # [-]
X_CO2 = []          # [-]
X_CO2_peak = []     # [-]

# integrate over time each oxidation step
for frame in frames_oxidation:
    nCO2_0 = integrate.trapz(frame['F_CO2_0'], frame['time'])       # [mol]
    nCO2_f = integrate.trapz(frame['F_CO2_out'], frame['time'])     # [mol]
    nCO_f = integrate.trapz(frame['F_CO_out'], frame['time'])       # [mol]
    n_CO2_in.append(nCO2_0)                                         # [mol]
    n_CO_out.append(nCO_f)                                          # [mol]
    X_CO2.append(round(nCO_f/(nCO2_f+nCO_f), 5))
    X_CO2_peak.append(frame['x_CO'].max()/0.0502)
    C_balance.append(round((nCO_f+nCO2_f)/(nCO2_0), 5))


print('Results')
print('X_CO2 = ', X_CO2)
print('X_CO2_peak = ', X_CO2_peak)
print('C_balance = ', C_balance)

n_H2_in = []        # [mol]
n_H2_out = []       # [-]
n_H2O_out = []      # [-]
X_H2 = []           # [-]

# integrate over time for each reduction step
for frame in frames_reduction:
    nH2_0 = integrate.trapz(frame['F_H2_0'], frame['time'])
    nH2_f = integrate.trapz(frame['F_H2_out'], frame['time'])
    nH2O_f = nH2_0-nH2_f
    X_H2.append(round(1-nH2_f/nH2_0, 5))
    n_H2O_out.append(nH2O_f)                                            # is equal to extent of reaction
    n_H2_in.append(nH2_0)
    n_H2_out.append(nH2_f)

R_feed = np.asarray(n_CO2_in)/np.asarray(n_H2_in)
print('R_feed =' , R_feed, 0.33/0.5)
# oxygen mass balance
n_H2O_out = np.asarray(n_H2O_out)
n_CO_out = np.asarray(n_CO_out)
O_balance = n_H2O_out/n_CO_out

print('X_H2 = ', X_H2)
print('O_balance = ', O_balance)





fig = plt.figure(figsize=(4.3, 3.8), dpi = 200)
plt.errorbar(cycle_number, C_balance,  yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'C0', label = '$(n_{\mathrm{CO},f}+n_{\mathrm{CO_2},f})/n_\mathrm{CO_2,0}$')
plt.errorbar(cycle_number,O_balance,yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'black', label = '$n_{\mathrm{H_2O},f}/n_{\mathrm{CO},f}$')
plt.plot([-0.5, 10.5],[1.0,1.0], linewidth=1.0, ls = '--', color = 'grey', label = '__nolegend__')
plt.xlim(0,11)
plt.ylabel("Mass balance [-]")
plt.xlabel("Cycle number [-]")
plt.xticks(cycle_number)
plt.legend(loc='upper right', ncol=1)
plt.ylim([0.85, 1.15])
plt.tight_layout()
plt.savefig(os.path.join('plots', 'SI_10cycles_MB' + '.png'), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'SI_10cycles_MB' + '.pdf'), bbox_inches='tight')
plt.show()


fig = plt.figure(figsize=(4.3, 3.5), dpi = 200)
ax = plt.subplot()
ax.set_xlabel('Cycle number [-]')
ax.set_ylabel('Millimoles per cycle [mmol]')
axr = ax.twinx()
axr.spines['right'].set_color('dimgrey')
axr.yaxis.label.set_color('dimgrey')
axr.tick_params(axis='y', colors='dimgrey')
axr.set_ylabel('CO$_2$ conversion [-]')
axr.plot(cycle_number, X_CO2,
             marker = '+', lw = 0.0, color = 'dimgrey', label = '$X_\mathrm{CO_2}$')

ax.errorbar(cycle_number, np.asarray(n_CO_out)*1000,  yerr=relative_error_nf*np.asarray(n_CO_out)*1000,
            elinewidth=1.0, capsize=2, marker = 's', lw = 0.0, color = 'C0', label = '$n_\mathrm{CO}$')
ax.errorbar(np.asarray(cycle_number)+0.1,np.asarray(n_H2O_out)*1000,yerr=relative_error_nf*np.asarray(n_H2O_out)*1000,
            elinewidth=1.0, capsize=2, marker = 'o', lw = 0.0, color = 'C3', label = '$n_\mathrm{H_2 O}$')
ax.set_xticks(cycle_number)
ax.set_xlim(0,11)
axr.set_ylim(0,1.0)
ax.set_ylim(0,10)

fig.legend(loc=(0.6,0.19), ncol=1)
plt.tight_layout()
plt.savefig(os.path.join('plots', 'figure_4' + '.png'), dpi =400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_4' + '.pdf'), bbox_inches='tight')
plt.show()