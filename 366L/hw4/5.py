import space_functions as sf
from mpl_toolkits import mplot3d
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
tf = (2451545.0+365.25)*86400 - T0
T0 = 0
dt = .5*86400

time_series = np.arange(T0, tf, dt)-T0

r_vec, v_vec = sf.orbit_prop_3body(r_0, v_0, T0, tf, dt)

elements = np.empty([np.shape(r_vec)[0], 6])

for i in range(np.shape(r_vec)[0]):
    elements[i, :] = sf.cart2elm(r_vec[i, :], v_vec[i, :], MU)

print(np.shape(elements))
fig, ax = plt.subplots(6, sharex=True)
try:
    for j in range(np.size(ax)):
        ax[j].plot(time_series, elements[:, j])
    # plt.show()
except Exception as e:
    print(Fore.RED+"There was a plotting error!")

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(r_vec[:, 0], r_vec[:, 1], r_vec[:, 2])

plt.show()

