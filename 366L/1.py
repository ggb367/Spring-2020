import space_functions as sf
import math as m
import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt


MU = 398600.4415

# E - [a, e, i, RAAN, omega, nu]
E = [26000, 0.72, 75, 90, 270, 0]  # time since periapsis passage is 0 for #1 so nu is 0

# #1
r_0, v_0 = sf.elm2cart(E, MU)
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

# #3
r_truth = np.multiply((h_norm**2/MU), np.divide(1, (1+np.multiply(eccent_norm, np.cos(nu)))))  # orbit equation
a = (r_truth[8]*(1+eccent_norm*np.cos(nu[8])))/(1-eccent_norm**2)  # use any r and nu paring to find aTo
v_truth = np.sqrt(MU*(np.divide(2, r_truth))-1/a)  # use r and a to find v, from vis-viva
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
r_rk = np.delete(r_rk, np.size(r_rk)-1)
v_rk = np.delete(v_rk, np.size(v_rk)-1)
# plot the energy
plt.plot(time_series/3600, energy-energy[0])

# #5
error_r = r_truth - r_rk
error_v = v_truth - v_rk
plt.figure()
plt.plot(time_series, error_r)

plt.show()

