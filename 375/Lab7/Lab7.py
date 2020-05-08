import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as int

def volt2load(volt): # return the load in kilograms
    gain = 1+100000/2700
    return np.divide((np.multiply(10823.017142857841*gain, volt)-5339.357926769584), 1000)

data_td = pd.read_csv('Data/lab7_TimeDomain.csv')
accel_raw_td = np.array(data_td['Accel'])

data = pd.read_csv('Data/Lab7.csv')
force_raw = np.array(data['Force'])
freq_raw = np.array(data['freq'])
accel_raw = np.array(data['Accel'])


# gain = 1+100000/2700
# accel_gain_td = np.multiply(gain, accel_raw_td)
accel_td = volt2load(accel_raw_td)
force = volt2load(force_raw)
freq = (freq_raw)
accel = volt2load(accel_raw)

accel_td = accel_td-np.average(accel_td) # zero the acceleration

vel = int.cumtrapz(accel_td)
pos = int.cumtrapz(vel)

pos_FreqDomain = np.fft.rfft(pos)

plt.plot(freq, np.abs(pos_FreqDomain[3:-1]))
# plt.plot(pos)
plt.title("Tip Position in Frequency Domain")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Average Position")

plt.figure()
plt.plot(freq, force)
plt.title("Tip Force in Frequency Domain")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Average Force")

plt.figure()
# plt.plot(freq, accel)
plt.plot(freq, np.fft.rfft(accel_td)[4:-1])
plt.title("Tip Acceleration in Frequency Domain")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Average Acceleration")

#
freq_resp = np.divide(pos_FreqDomain[3:-1], force)
plt.figure()
plt.plot(freq, freq_resp)
plt.xlabel("Frequency [Hz]")
plt.title("Tip Displacement Divided by Tip Force")
plt.show()
