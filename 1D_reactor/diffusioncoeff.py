"""Source for Formula for D_AB and diffusion volumes: "The properties of gases and liquids, volume 5. Mcgraw-hill New York, 2001, chapter 11" Fullers method,
 can be found here: https://www.accessengineeringlibrary.com/content/book/9780070116825/chapter/chapter11#/p2000ad2b99711_9002"""
T = 1073
p = 10**5

diff_vol_H2O = 13.1
diff_vol_H2 = 6.12
diff_vol_CO2 = 26.9
diff_vol_CO = 18.0

M_H2 = 2.0
M_H2O = 18.0
M_CO2 = 44.0
M_CO = 28.0

M_H2O_H2 = 2*1/((1/M_H2O)+(1/M_H2))
M_CO_CO2 = 2*1/((1/M_CO)+(1/M_CO2))

D_H2O_H2 = (0.00143 * T ** 1.75) / ((p / 10 ** 5) * (M_H2O_H2**(1/2)*(diff_vol_H2O**(1/3)+diff_vol_H2**(1/3))**2) )
D_CO_CO2 = (0.00143 * T ** 1.75) / ((p / 10 ** 5) * (M_CO_CO2**(1/2)*(diff_vol_CO**(1/3)+diff_vol_CO2**(1/3))**2) )

print("T, p = ", T, p)
print("D_H2O_H2: [cm^2/s]=", D_H2O_H2)
print("D_CO_CO2: [cm^2/s]=", D_CO_CO2)

