import space_functions as sf
import math as m
import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt


MU = 398600.4415

# E - [a, e, i, RAAN, omega, nu]
E = [26000, 0.72, 75, 90, -90, 0]  # time since periapsis passage is 0 for #1 so nu is 0

# #1
r_0, v_0 = sf.elm2cart(E, MU)
print("This is r:")
print (r_0)
print("This is v:")
print (v_0)

# #2
# time_series = np.arange(0, 86400*3, 20)
# nu_0 = 0
# r_norm = lg.norm(r)
# h = np.cross(r,v)
# h_norm = lg.norm(h)
# energy = np.cross(v, h)/MU - np.divide(r, r_norm)
# energy_norm = lg.norm(energy)
# a = (r_norm*(1+energy_norm*np.cos(nu_0)))/(1-energy_norm**2)  # semimajor axis
# n = np.sqrt(MU/a**3)  # mean motion
#
# nu, E, M = sf.orbit_prop(time_series, n, energy, 0)
#
# r = np.multiply((h_norm**2/MU), np.divide(1, (1+np.multiply(energy_norm, np.cos(nu)))))  # orbit equation
# a = (r[8]*(1+energy_norm*np.cos(nu[8])))/(1-energy_norm**2)  # use any r and nu paring to find aTo
# v = np.sqrt(MU*(np.divide(2, r))-1/a)  # use r and a to find v, from vis-viva
# KE = 0.5*np.square(v)
# PE = -np.divide(MU, r)
# TE = KE+PE
# dTE = TE - (KE[0]+PE[0])

# #3
# Todo

# #4
r, v = sf.orbit_prop_rk(r_0, v_0,0, 86400*3, 20 )
