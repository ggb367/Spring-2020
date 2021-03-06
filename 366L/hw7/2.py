import numpy as np
import numpy.linalg as lg

RE = 6378.1363
MU = 398600.4415
p_srp = 4.57e-6*1000**2  # m^2/km^2 * W/m^2 = W/km^2
c_r = 1.5
A_m = 0.1/1000**2  # km^2/m^2 *m^2/kg = km^2/kg

r_sun = np.array([2379260, 148079334, -1009936])
r = np.array([2781, 5318, -5629])

a_srp = np.multiply(p_srp*c_r*A_m, np.divide(-1*r_sun+r, lg.norm(r_sun-r)))
print(a_srp)

