import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# part_1 = pd.read_excel('lab_2_data/Part1.xlsx')
# num_23_ab = pd.read_excel('lab_2_data/2_3_ab.xlsx')
num_23_c = pd.read_excel('lab_2_data/2_3_c.xlsx')
# num_23_d = pd.read_excel('lab_2_data/2_3_d.xlsx')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

# Get time array
time_array = num_23_c["Time"]  # I think this is in MJD
time_array = np.array((time_array-time_array[0])*86400) #Normalization to seconds since epoch

# Get data

thermocouple = np.array(num_23_c["Thermocouple"])
thermistor = np.array(num_23_c["Thermistor"])
r_t = np.divide(np.multiply(1500, thermistor), 5-thermistor)
T = np.power(a+np.multiply(b, np.log(r_t))+np.multiply(c, np.power(np.log(r_t), 3)), -1)

#plot data
plt.figure()
plt.plot(time_array, thermocouple)

plt.figure()
plt.plot(time_array, thermistor)

a = 0.044747
b = -0.01388
c = 0.000357


# # calculate time constant
# val_thermocouple = (thermocouple.max())*(1-1/np.e)
# val_thermistor = thermistor.min()*(1-1/np.e)
#
# value_couple_tc, idx_tc_couple = find_nearest(thermocouple[0:250320], val_thermocouple)
# thermocouple_timeconstant = time_array[idx_tc_couple]
# value_istor_tc, idx_tc_istor = find_nearest(thermistor[0:250320], val_thermistor)
# thermistor_timeconstant = time_array[idx_tc_istor]
# print("thermistor time constant:")
# print(thermistor_timeconstant)
# print("thermocouple time constant:")   # Important to note that this is an estimate as data is discrete and there are inconsitencies
# print(thermocouple_timeconstant)
#
# # Time to steady state
# max_thermocouple = np.argmax(thermocouple)
# min_thermistor = np.argmin(thermistor)
#
# couple_ss = time_array[max_thermocouple]
# istor_ss = time_array[min_thermistor]
# print("Thermocouple time to steady state")
# print(couple_ss)
# print("Thermoistor time to steady state")
# print(istor_ss)

plt.show()
