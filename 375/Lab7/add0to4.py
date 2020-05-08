import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt

data = pd.read_csv('Data/Lab7.csv')
freq = np.array(data['freq'])
accel = np.array(data['Accel'])
force = np.array(data['Force'])

slope, intercept, r_value, p_value, std_err = linregress(freq[0:20], accel[0:20])

hz = [0, 1, 2, 3, 4]

temp_accel = np.multiply(slope, hz)+intercept
new_accel = np.concatenate((temp_accel, accel), axis = 0)

slope, intercept, r_value, p_value, std_err = linregress(freq[0:8], force[0:8])

temp_force = np.multiply(slope, hz)+intercept
new_force = np.concatenate((temp_force, force), axis = 0)
new_freq = np.concatenate((hz, freq), axis = 0)

# plt.plot(new_freq, new_accel)
# plt.show()

d = {'freq': new_freq, 'Accel': new_accel, 'Force': new_force}

estimated_data = pd.DataFrame(data=d)

estimated_data.to_csv('Data/Lab7_start0.csv')

print('complete')
