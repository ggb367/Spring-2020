import numpy as np
import space_functions as sf

t1 = 0
t2 = 30
t3 = 60
az_el1 = [651.343, 45.0, 26.411]
az_el2 = [865.398, 45.715, 17.782]
az_el3 = [1088.356, 46.102, 12.210]

r_1 = sf.razel2SEZ(az_el1)
r_2 = sf.razel2SEZ(az_el2)
r_3 = sf.razel2SEZ(az_el3)
# print(r_1, r_2, r_3)
Q_SEZ_ITRF = np.transpose(sf.ITRF2SEZ(0, 0))
# print(Q_SEZ_ITRF)
# print(Q_SEZ_ITRF)
rho_1 = np.matmul(Q_SEZ_ITRF, r_1)
rho_2 = np.matmul(Q_SEZ_ITRF, r_2)
rho_3 = np.matmul(Q_SEZ_ITRF, r_3)
r1 = np.matmul(sf.R3(-sf.earth.rot_speed*t1), (rho_1+np.array([sf.earth.radius, 0, 0])))
r2 = np.matmul(sf.R3(-sf.earth.rot_speed*t2), (rho_2+np.array([sf.earth.radius, 0, 0])))
r3 = np.matmul(sf.R3(-sf.earth.rot_speed*t3), (rho_3+np.array([sf.earth.radius, 0, 0])))
# r1 = [-1970, 4033, -5413]
# r2 = [-2111, 3421, -5800]
# r3 = [-2228, 2770, -6122]
# print(r1, r2, r3)
r2, v2 = sf.IOD_Gibbs(r1, r2, r3, t1, t2, t3)
print(v2)
