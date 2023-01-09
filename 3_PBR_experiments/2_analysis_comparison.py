import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cantera as ct
from scipy import integrate
import os

"""In this file we analyse the experimental data from the comparison PBR experiment."""

# Get the theoretical thermodynamic limit of co-feeding conversions first.
T = 790+273  # average reactor temperature
p = 100000
RWGS = ct.Solution('RWGS-database.cti')
RWGS.TP = T, p
RWGS.X = {'H2': 1.0, 'CO2': 0.66}
RWGS.equilibrate('TP')
RWGS_X_CO2 = RWGS.X[RWGS.species_index('CO')]/(RWGS.X[RWGS.species_index('CO')]+ RWGS.X[RWGS.species_index('CO2')])
RWGS_X_H2 = RWGS.X[RWGS.species_index('H2O')]/(RWGS.X[RWGS.species_index('H2O')]+ RWGS.X[RWGS.species_index('H2')])

# RELATIVE ERROR IN THE MASS BALANCE
# 2 % error in bottle concentration. 2 % error in gas analysis signals, and 1 % error in the mass flow.
# The errors in the gas-analysis are used twice in each mass balance.
# This neglects the numerical errors introduced during integration of the signals.
relative_error_MB = np.sqrt(2 * (0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))
print("relative error mass balance:", relative_error_MB)

# Bottle concentrations and molar volumes
x_CO2_Bottle = 0.0501           # rest Ar => Accuracy 2%
x_H2_Bottle = 0.0486            # rest Ar => Accuracy 2%
V_0 = 22.414                    # Molar Volume [dm^3/mol] = [l/mol]


# Import Data
filename = "comparison_experimental_dataset.csv"
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

#Check CH4 and O2 signals are zero
print("CH4, mean x, max x =", data['x_CH4'].mean(),  data['x_CH4'].max())
print("O2, mean x, max x =", data['x_O2'].mean(),  data['x_O2'].max())

# Calculate molar inflow rates [mol/s]
data['F_H2_0'] = data['v_H2']*x_H2_Bottle/V_0/60
data['F_H2_Ar_0'] = data['v_H2']*(1 - x_H2_Bottle)/V_0/60
data['F_CO2_0'] = data['v_CO2']*x_CO2_Bottle/V_0/60
data['F_CO2_Ar_0'] = data['v_CO2']*(1 - x_CO2_Bottle)/V_0/60
data['F_Ar_0'] = data['v_Ar']/V_0/60
data['F_total_0'] = (data['v_H2'] + data['v_CO2'] + data['v_Ar'])/V_0/60
data['F_Ar_total_0'] = data['F_H2_Ar_0']+data['F_CO2_Ar_0']+data['F_Ar_0']

# calculate molar outflow rates
# H2O has been condensed out which reduces the flow rate which needs to be accounted for.
data['F_H2_out'] = data['x_H2']*(data['F_Ar_total_0']+data['F_CO2_0'])/(1-data['x_H2'])
data['F_H2O_out'] = data['F_H2_0'] - data['F_H2_out']
# total outflow minus H2O
data['F_total_out'] = data['F_H2_out'] + data['F_CO2_0']  + data['F_Ar_total_0']
data['F_CO2_out'] = data['x_CO2']*data['F_total_out']
data['F_CO_out'] = data['x_CO']*data['F_total_out']


# separate a 100 second period of co-feeding data to analsyse.
cofeed_df = data[data['time'] > 9600]
cofeed_df = cofeed_df[cofeed_df['time'] < 9700]

# calculate conversion and mass balances during co-feeding by averaging over the 100 seconds
cofeed_df['X_H2'] = 1-cofeed_df['F_H2_out']/cofeed_df['F_H2_0']
X_CO2_cofeed = cofeed_df['x_CO'].mean()/(cofeed_df['x_CO'].mean()+cofeed_df['x_CO2'].mean())
X_H2_cofeed_2 = X_CO2_cofeed*0.33/0.5                       # X_H2 = F_H2O_out/F_H2,0
X_H2_cofeed = cofeed_df['X_H2'].mean()             # X_H2 = X_CO2*F_CO2/F_H2,0
C_balance_cofeed = (cofeed_df['F_CO_out'].mean()+cofeed_df['F_CO2_out'].mean())/cofeed_df['F_CO2_0'].mean()
O_balance_cofeed = cofeed_df['F_H2O_out'].mean()/cofeed_df['F_CO_out'].mean()

# Next analyse the chemical looping cycles by breaking up the data at the switching points between the steps
# where the Ar flow rate set point was 0.1 slm
timeswitches = []
timeswitches.append(0)
time_skip = 10*60  # skip ten minutes to ave time
safety = 60        # Safety after switching

i = 0
while i < len(data['Ar MFC sp']):
    Ar_flow = data['Ar MFC sp'].iloc[i]
    if Ar_flow == 0.1:
        timeswitches.append(data['time'].iloc[i])
        i += time_skip  # skip 10min
    else:
        i += 5  # skip 5sec to next value

timeswitches.pop() # drop the co-feeding section


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

print('Results = first cycle, countercurrent, cocurrent, cofeed')
print('X_CO2 = ', X_CO2, X_CO2_cofeed)
print('X_CO2_peak = ', X_CO2_peak)
print('C_balance = ', C_balance, C_balance_cofeed)

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

print('X_H2 = ', X_H2, X_H2_cofeed, X_H2_cofeed_2)
print('O_balance = ', O_balance, O_balance_cofeed)




width = 0.25
# x_label = ['Cycle 1', 'Cycle 2', 'Cycle 3', 'Cycle 4', 'Cycle 5', 'Cycle 6', 'Cycle 7', 'Cycle 8', 'Cycle 9', 'Cycle 10']
x_label = ['Countercurrent', 'Cocurrent', 'Cofeed']

# Plots of conversions
X_CO2_plot = [X_CO2[1], X_CO2[2], X_CO2_cofeed]
X_H2_plot = [X_H2[1], X_H2[2], X_H2_cofeed_2]


fig = plt.figure(figsize=(4.3, 3.8))


#plt.plot([0.0+width/2, 1.0+width/2, 2.0 + width/2],[R_feed[1], R_feed[2], 0.666], marker = 's', lw = 0.0, color = 'black', label = 'CO$_2$/H$_2$')
plt.errorbar([0.0, 1.0, 2.0],[C_balance[1], C_balance[2], C_balance_cofeed],  yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'C0', label = '$(n_{\mathrm{CO},f}+n_{\mathrm{CO_2},f})/n_\mathrm{CO_2,0}$')
plt.errorbar([0.05, 1.05, 2.05],[O_balance[1], O_balance[2], O_balance_cofeed],yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'black', label = '$n_{\mathrm{H_2O},f}/n_{\mathrm{CO},f}$')

plt.plot([-0.5, 2.5],[1.0,1.0], linewidth=1.0, ls = '--', color = 'grey', label = '__nolegend__')

#plt.xlabel("Cycle")
plt.ylabel("Mass balance [-]")
#plt.title("Change of concentration")

# plt.grid(linestyle='--')
plt.xticks(np.arange(len(X_CO2_plot)) , x_label)
plt.legend(loc='upper left', ncol=1)

plt.ylim([0.85, 1.15])
plt.tight_layout()
plt.savefig(os.path.join('plots', 'SI_MB_CT_CC_cofeed' + '.png'), dpi =400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'SI_MB_CT_CC_cofeed' + '.pdf'), bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(4.5, 3.9))

plt.bar(np.arange(len(X_CO2_plot)), X_CO2_plot, color='C0', width=width, edgecolor='black', label='$X_\\mathrm{CO_2}$')
plt.bar(np.arange(len(X_H2_plot)) + width, X_H2_plot, color='C3', width=width, edgecolor='black', label='$X_\\mathrm{H_2}$')
#plt.plot([0.0+width/2, 1.0+width/2, 2.0 + width/2],[R_feed[1], R_feed[2], 0.666], marker = 's', lw = 0.0, color = 'black', label = 'CO$_2$/H$_2$')
plt.plot([0.0, 1.0],[X_CO2_peak[1], X_CO2_peak[2]], marker = 's', markersize = 4,  lw = 0.0, color = 'dimgrey', label = '$X_\\mathrm{CO_2}$ peak')
plt.plot([2.0, 2.0+width],[RWGS_X_CO2, RWGS_X_H2], marker = '+', markersize = 9,  lw = 0.0, color = 'dimgrey', label = 'cofeed TL')

#plt.xlabel("Cycle")
plt.ylabel("Conversion extent [-]")
#plt.title("Change of concentration")
plt.ylim([0,1.0])
# plt.grid(linestyle='--')
plt.xticks(np.arange(len(X_CO2_plot)) + width / 2, x_label)
plt.legend(loc='upper right',
               ncol=2)
plt.tight_layout()

plt.savefig(os.path.join('plots', 'figure_3b' + '.png'))
plt.savefig(os.path.join('plots', 'figure_3b' + '.pdf'), bbox_inches='tight')
plt.show()