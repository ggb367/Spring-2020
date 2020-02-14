import numpy as np
import matplotlib.pyplot as plt

altitude = np.arange(0, 80000, 1000)
density_true = np.empty(np.size(altitude))
i = 0
for alt in altitude:
    if alt <=36089:
        tau = 518.69 - alt*3.5662e-3
        density_true[i] = (tau**4.256)*6.6277e-15
    if alt > 36089 and alt <= 65617:
        pressure = 2678.4*np.e**(alt*-4.8063e-5)
        density_true[i] = pressure*1.4939e-6
    if alt > 65617:
        tau = 389.99+(alt-65617)*(5.4864e-4)
        density_true[i] = (tau**-35.164)*2.2099e87
    i = i+1

fig, ax = plt.subplots()

ax.plot(altitude, density_true, label="True")

rho_s = 2.3769e-3
rho_t = 7.0613e-4
rho_plus = 1.7083e-4
h_t = 36089
h_plus = 65617
density_est = np.empty(np.size(altitude))
i = 0
for alt in altitude:
    if alt <=36089:
        density_est[i] = rho_s*np.e**(-alt/29730)
    if alt > 36089 and alt <= 65617:
        density_est[i] = rho_t*np.e**(-(alt-h_t)/20806)
    if alt > 65617:
        density_est[i] = rho_plus*np.e**(-(alt-h_plus)/20770)
    i = i+1

ax.plot(altitude, density_est, label="Exponential")
ax.legend()
plt.show()

