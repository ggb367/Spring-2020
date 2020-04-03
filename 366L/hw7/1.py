import numpy as np
import numpy.linalg as lg

RE = 6378.1363
MU = 398600.4415
J2 = 0.0010826267
J3 = -0.0000025327

r = [0.009, 4595.737, 4595.731]
r_norm = lg.norm(r)
a_2 = np.multiply(3*MU*J2*RE**2/(2*r_norm**5), [r[0]*(5*(r[2]/r_norm)**2-1), r[1]*(5*(r[2]/r_norm)**2-1), r[2]*(5*(r[2]/r_norm)**2-3)])
a_3 = np.multiply(-5*J3*MU*RE**3/(2*r_norm**7), [r[0]*(3*r[2]-7*r[2]**3/r_norm**2), r[1]*(3*r[2]-7*r[2]**3/r_norm**2), 6*r[2]**2-7*r[2]**4/r_norm**2-(3/5)*r_norm**2])
a_p = a_2+a_3

print(a_p)
