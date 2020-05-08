import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

load_ext_al = pd.read_csv("Data/Al/Al_Load_Extension.csv")
load_ext_poly = pd.read_csv("Data/Poly/Poly_Load_Extension.csv")
# print("loaded")
al_time = np.array(load_ext_al["Time"]) # s
al_load = np.array(load_ext_al["Load"]) # N
al_ext = np.array(load_ext_al["Extension"]) # mm

poly_time = np.array(load_ext_poly["Time"]) # s
poly_load = np.array(load_ext_poly["Load"]) # N
poly_ext = np.array(load_ext_poly["Extension"]) # mm

# al_stress = np.divide(al_load, (1.59e-3)*(25.45e-3-9.5e-3))
# al_stress_max = 503e6
# print(max(al_stress))
# K_al_test = np.divide(al_stress_max, al_stress)
# d_D_al = 9.5/25.45
# K_al_theory = 3-3.14*d_D_al+3.667*d_D_al**2-1.527*d_D_al**3
# print(K_al_theory)
# al_strain = (al_ext)/9.5

poly_stress = np.divide(poly_load, (2.4e-3)*(25.5e-3-9.6e-3))
poly_stress_max = max(poly_stress) #62e6
print(max(poly_stress))
K_poly_test = np.divide(poly_stress_max, poly_stress)
d_D_poly = 9.6/25.50
K_poly_theory = 3-3.14*d_D_poly+3.667*d_D_poly**2-1.527*d_D_poly**3
print(K_poly_theory)
poly_strain = (poly_ext)/9.5

print
plt.figure()
plt.plot(al_ext, al_load)
plt.xlabel("Force, N")
plt.ylabel("Displacement, mm")
plt.title("Aluminum Force vs Displacement")

plt.figure()
plt.plot(poly_ext, poly_load)
plt.xlabel("Force, N")
plt.ylabel("Displacement, mm")
plt.title("Polycarbonate Force vs Displacement")

# plt.figure()
# plt.plot(poly_load, K_poly_test, color='orange')
# plt.axhline(y=K_poly_theory, xmin=poly_time[0], xmax=poly_time[-1])
# plt.ylim([0, 10])
# plt.legend(["Test", "Theory"])
# plt.title("Polycarbonate Stress Concentration Factor vs Load")
# plt.xlabel("Load, N")
# plt.ylabel("K_t")

plt.show()
