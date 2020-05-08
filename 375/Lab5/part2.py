import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#identify the files to be used
file_names = []
for path, subdirs, files in os.walk("Data/Part 2"):
    for name in files:
        file_names.append(os.path.join(path, name))

#get data_td from each file
GF = 2.05
first = True
count = 0
for file in file_names:
    csv = pd.read_csv(file)
    #get values from df
    time = np.array(csv["X_Value"])
    strain = np.abs(np.divide(np.multiply(-4,np.array(csv["Bridge (Collected)"])),(GF*(1+np.multiply(2, np.array(csv["Bridge (Collected)"]))))))*1000
    if first:
        strain_store = np.zeros([np.size(file_names), np.size(strain)])
        first = False
    strain_store[count, :] = strain
    count = count+1
#average the strain
strain_avg = np.mean(strain_store, axis=0)
#perform a fourier transform
strain_transform = np.fft.rfft(strain_avg)  # The natural freq is 11Hz and 36Hz
#plot the transform
plt.plot(np.abs(strain_transform))
plt.xlim(2, 50)
plt.ylim(-0.0005, 1)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (Millistrain)")
plt.title("Fourier Transform of Average Strain Measurment")
plt.minorticks_on()
plt.grid(which='minor', linestyle='--', color='silver', axis='x')
plt.grid(which='major', color='grey')
plt.show()
