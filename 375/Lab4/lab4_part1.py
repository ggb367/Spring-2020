import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as sgl

#part 1
file = pd.read_csv('data/part1_reduced.csv')

pressure = np.array(file["Pressure"])
strain = np.array(file["Strain"])


b, a = sgl.bessel(3, 0.01)  # create a bessel lowpass filter
pressure_denoise = sgl.filtfilt(b, a, pressure)  #denoise pressure readings
strain_denoise = sgl.filtfilt(b, a, strain)
pressure_denoise = np.delete(pressure_denoise, np.s_[0:200])
strain_denoise = np.delete(strain_denoise, np.s_[0:200])
# pressure steps: [0, 5, 11, 16, 20, 25, 32]

V_unstrained = np.average(strain_denoise[5000:35000])

GF = 2.05  # gauge factor
r = 2  # radius
thic = 0.125  # thickness
STRAIN_GAIN = 303/4.5
PRESS_GAIN = 500
V_ex_strain = 3.3
V_ex_press = 2.0
poisson_ratio = (0.3182+0.3487)/2 # average from internet
young_mod = (13.5+21.4)/2  # from internet
# poisson_ratio = 0.37 #
# young_mod = 354e3  # what's up with this?

V_balanced = np.abs(strain_denoise - V_unstrained)
delta_V = np.multiply(V_balanced, STRAIN_GAIN)
V_r = delta_V/V_ex_strain

# so this is a bit too high, dividing it by five makes it more reasonable
press_strain = np.divide(np.multiply(4*young_mod*thic, V_r), GF*(r*(1+poisson_ratio)))

plt.plot(press_strain)
pressure_denoise = pressure_denoise- np.average(pressure_denoise[5000:35000]) #set the press to zero
true_pressure = np.multiply(pressure_denoise, PRESS_GAIN*V_ex_press)
pressure_df = pd.DataFrame(data={'True Pressure': true_pressure, 'Strain Pressure': press_strain})
pressure_df.to_csv('data_td/calculated_press.csv')
plt.plot(true_pressure)
plt.legend(["Strain", "True"])
plt.xlabel("Data")
plt.ylabel("Pressure (psi)")
plt.title("Comparison Between Pressure Transducer and Strain Gauge Pressure")
plt.grid(True)

plt.figure()
percent_error = np.divide(np.abs(press_strain[50000::]-true_pressure[50000::]), true_pressure[50000::])
plt.plot(percent_error)
x_range = np.arange(0, np.size(percent_error), 1)

m, b = np.polyfit(x_range, percent_error, 1)
plt.plot(x_range, np.multiply(m, x_range)+b)
plt.title("Percent Error of True vs Calculated Pressure")
plt.ylabel("Percent Error")
plt.xlabel("Data")
plt.legend(["Error", "Linear Regression"])
plt.grid(True)
plt.show()
