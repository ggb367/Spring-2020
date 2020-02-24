import space_functions as sf
import matplotlib.pyplot as plt
import numpy as np
from colorama import Fore


MU = 398600

a = 500000
e = 0.3
i = 50
RAAN = 350
arg_peri = 90
anomaly = 0

r_0, v_0 = sf.elm2cart([a, e, i, RAAN, arg_peri, anomaly], MU)

T0 = 2451545.0*86400
tf = ((2451545.0+365.25)*86400 - T0)
T0 = 0
dt = .5*86400

time_series = np.arange(T0, tf, dt)-T0

r_vec, v_vec = sf.orbit_prop_3body(r_0, v_0, T0, tf, dt)

elements = np.empty([np.shape(r_vec)[0], 6])

for i in range(np.shape(r_vec)[0]):
    elements[i, :] = sf.cart2elm(r_vec[i, :], v_vec[i, :], MU)

print(np.shape(elements))
fig, ax = plt.subplots(5, sharex=True)
elementals=["a[km]", "e[deg]", "i[deg]", "RAAN[deg]", "Arg Periapsis [deg]"]
try:
    for j in range(np.size(ax)):
        ax[j].plot(np.divide(time_series,86400), elements[:, j])
        ax[j].set_ylabel(elementals[j])
except Exception as e:
    print(Fore.RED+"There was a plotting error!")

plt.suptitle("Orbit Propagation with Solar Third Body Perturbations")
plt.xlabel("Time[Days]")

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(r_vec[:, 0], r_vec[:, 1], r_vec[:, 2])

h = np.empty(np.shape(r_vec))
for i in range(len(r_vec)):
    h[i, :] = np.cross(r_vec[i, :], v_vec[i, :])

fig, ax = plt.subplots(2, sharex=True)
ax[0].plot(np.divide(time_series,86400), h[:, 0])
ax[0].set_ylabel("h_i [m^2/s]")
ax[1].plot(np.divide(time_series,86400), h[:, 1])
ax[1].set_ylabel("h_j [m^2/s]")
plt.suptitle("Change in Angular Momentum")
plt.xlabel("Time [Days]")

plt.show()

