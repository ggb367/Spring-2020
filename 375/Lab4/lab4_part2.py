import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as sgl

#part2
file = pd.read_csv('data/part2_reduced.csv')

V_signal = file["Strain"][230000:270000]

b, a = sgl.bessel(3, 0.01)  # create a bessel lowpass filter
V_denoise = sgl.filtfilt(b, a, V_signal)  #denoise strain readings

GF = 2.05
young_mod = 10153 # ksi ref: https://www.azom.com/article.aspx?ArticleID=6619
poisson_ratio = 0.33
V_ex_strain = 3.3
STRAIN_GAIN = 303
thic = 0.0065
radius = 1.28

V_unstrained = np.average(V_denoise[21050:40000])
V_balanced = np.abs(V_denoise - V_unstrained)
V_balanced[21050:40000891]=0
delta_V = np.multiply(V_balanced, STRAIN_GAIN)
V_r = delta_V/V_ex_strain

strain = np.abs(np.divide(np.multiply(-4,V_r),(GF*(1+np.multiply(2, V_r)))))

press = (np.divide(np.multiply(young_mod*thic, strain), radius*(1-poisson_ratio/2)))

part2_df = pd.DataFrame(data={'Strain': strain, 'Pressure': press})
part2_df.to_csv('data_td/part2_calculated.csv')

plt.plot(strain)
plt.title("")
plt.figure()
plt.plot(press)
plt.title("Pressure of a Diet Pepsi Can")
plt.ylabel("Pressure (psig)")
plt.xlabel("Data")
plt.grid(True)
plt.show()
