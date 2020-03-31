import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#get list of files
file_names = []
for path, subdirs, files in os.walk("Data/unbalanced"):
    for name in files:
        file_names.append(os.path.join(path, name))
volt_average = []
# print(file_names)
for file in file_names:
    csv = pd.read_csv(file)

    #get values from df
    volt_average.append(float(csv["Bridge (Arith. Mean)"][0]))
x = [110, 130, 130, 150, 170, 170, 190, 210, 210, 230, 250, 50, 50, 70 , 90, 90, 0, 0]
# plt.plot(x, volt_average)
strain = []
strain_theory = []
GF = 2.05  # gauge factor
r = 2  # radius
thic = 0.0037211  # thickness
STRAIN_GAIN = 500#/4.5
L = 0.33655
b = 0.0254762
poisson_ratio = (0.3182+0.3487)/2 # average from internet
young_mod = 70000000000#(13.5+21.4)/2  # from internet
I =b*thic**3/12

for avg in volt_average:
    strain.append((np.divide(np.multiply(-4,avg),(3*GF*(1+np.multiply(2, avg)))))*1000)
for mass in x:
    force = mass*9.81/1000
    strain_theory.append(force*L*thic/(2*young_mod*I)*1000)

plt.scatter(x, strain-min(strain), label='Measured')
plt.errorbar(x, strain-min(strain), np.std(strain-min(strain)), fmt='none')
plt.plot(x, strain_theory, label='Theory')
plt.legend(loc=4)
plt.xlabel("Mass (Grams)")
plt.ylabel("Strain (Millistrain)")
plt.title("Measured Strain From Unbalanced Bridge vs Mass ")
plt.grid(True)

plt.show()

