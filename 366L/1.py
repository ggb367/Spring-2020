import math as m

import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lg

import space_functions as sf

MU = 398600.4415

# E - [a, e, i, RAAN, omega, nu]
elements = [26000, 0.72, 75, 90, 270, 0]  # time since periapsis passage is 0 for #1 so nu is 0

# #1
r_0, v_0 = sf.elm2cart(elements, MU)
print("This is r:")
print (r_0)
print("This is v:")
print (v_0)

#2
time_series = np.arange(0, 86400*3, 20)
nu_0 = 0
r_norm = lg.norm(r_0)
h = np.cross(r_0,v_0)
h_norm = lg.norm(h)
eccent = np.cross(v_0, h)/MU - np.divide(r_0, r_norm)
eccent_norm = lg.norm(eccent)
a = (r_norm*(1+eccent_norm*np.cos(nu_0)))/(1-eccent_norm**2)  # semimajor axis
n = np.sqrt(MU/a**3)  # mean motion

nu, E, M = sf.orbit_prop(time_series, n, eccent, 0)
tf = np.size(nu)-1
print(sf.rad2deg(nu[tf]))
print(sf.rad2deg(m.fmod(M[tf], 360)))

# #3
r_truth = np.empty([np.size(nu), 3])
v_truth = np.empty([np.size(nu), 3])
for i in range(np.size(nu)):
    elements[5] = nu[i]*180/m.pi
    r_truth[i], v_truth[i] = sf.elm2cart(elements, MU)
a = (r_truth[8]*(1+eccent_norm*np.cos(nu[8])))/(1-eccent_norm**2)  # use any r and nu paring to find aTo
print("r at tf:")
print(r_truth[tf])
print("v at tf:")
print(v_truth[tf])

# #4
r_3d, v_3d = sf.orbit_prop_rk(r_0, v_0, 0, 86400*3, 20)
r_rk = np.sqrt(np.power(r_3d[:, 0], 2)+np.power(r_3d[:, 1], 2)+np.power(r_3d[:, 2], 2))
v_rk = np.sqrt(np.power(v_3d[:, 0], 2)+np.power(v_3d[:, 1], 2)+np.power(v_3d[:, 2], 2))

energy = np.divide(np.power(v_rk, 2), 2) - np.divide(MU, r_rk)
# remove last element because of rk does one extra
energy = np.delete(energy, np.size(energy)-1)
# plot the energy
plt.plot(time_series/3600, energy-energy[0])
plt.title("Change in Energy")
plt.xlabel("Time [Hours]")
plt.ylabel("Energy [kJ]")

# #5
error_r = r_truth - r_3d[0:-1, :]
error_v = v_truth - v_3d[0:-1, :]
print("r error at tf:")
print(error_r[tf, :])
print("v error at tf:")
print(error_v[tf, :])
fig, axs = plt.subplots(3, sharex=True)
axs[0].plot(time_series/3600, error_r[:, 0])
axs[0].set_title("Error in r_x")
axs[1].set_ylabel("Error [km]")
axs[1].plot(time_series/3600, error_r[:, 1])
axs[1].set_title("Error in r_y")
axs[2].plot(time_series/3600, error_r[:, 2])
axs[2].set_title("Error in r_z")
plt.xlabel("Time [Hours]")
fig.suptitle("Position Error Between Elements and Numeric Propagation")

fig, axs = plt.subplots(3, sharex=True)
axs[0].plot(time_series/3600, error_v[:, 0])
axs[0].set_title("Error in v_x")
axs[1].plot(time_series/3600, error_v[:, 1])
axs[1].set_title("Error in v_y")
axs[1].set_ylabel("Error [km/s]")
axs[2].plot(time_series/3600, error_v[:, 2])
axs[2].set_title("Error in v_z")
plt.xlabel("Time [Hours]")
fig.suptitle("Velocity Error Between Elements and Numeric Propagation")

plt.show()

