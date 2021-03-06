import space_functions as sf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

elements = [70000, 0.01, 20, 180, 25, 0]

r_0, v_0 = sf.elm2cart(elements, sf.earth.mu)

T0 = 2451545.0*86400
TF = ((2451545.0+365.25)*86400 - T0)
T0 = 0
dT = .5*86400

time_series = np.arange(0, 365.25, 0.5)
print("Starting Propagation")
r_vec, v_vec = sf.orbit_prop_3body_RV(r_0, v_0, T0, TF, dT)
print("Propagation Complete!")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot3D(r_vec[:, 0], r_vec[:, 1], r_vec[:, 2])

elements = np.empty([np.shape(r_vec)[0], 6])

for i in range(np.shape(r_vec)[0]):
    elements[i, :] = sf.cart2elm(r_vec[i, :], v_vec[i, :], sf.earth.mu)

print(np.shape(elements))
fig, ax = plt.subplots(3, sharex=True)
elementals=["a[km]", "e", "i[deg]", "RAAN[deg]", "Arg Periapsis [deg]"]
try:
    for j in range(np.size(ax)):
        ax[j].plot(time_series, elements[:, j])
        ax[j].set_ylabel(elementals[j])
except Exception as e:
    print(Fore.RED+"There was a plotting error!")

plt.suptitle("Orbit Propagation with Solar Third Body Perturbations")
plt.xlabel("Time[Days]")

fig, ax = plt.subplots(2, sharex=True)
elementals=["RAAN[deg]", "Arg Periapsis [deg]"]
try:
    for j in range(np.size(ax)):
        ax[j].plot(time_series, elements[:, j+3])
        ax[j].set_ylabel(elementals[j])
except Exception as e:
    print(Fore.RED+"There was a plotting error!")

plt.suptitle("Orbit Propagation with Solar Third Body Perturbations")
plt.xlabel("Time[Days]")

plt.show()
