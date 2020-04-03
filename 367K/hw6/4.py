import matplotlib.pyplot as plt
import numpy as np


def density_calc(alt):
    density_true = 0

    if alt <= 36089:
        tau = 518.69 - alt * 3.5662e-3
        density_true = (tau ** 4.256) * 6.6277e-15
    if 36089 < alt <= 65617:
        pressure = 2678.4 * np.e ** (alt * -4.8063e-5)
        density_true = pressure * 1.4939e-6
    if alt > 65617:
        tau = 389.99 + (alt - 65617) * 5.4864e-4
        density_true = (tau ** -35.164) * 2.2099e87
    return density_true

def temp_calc(alt):
    temp = 0

    if alt <= 36089:
        temp = 518.69 - 3.5662e-3 * alt
    if 36089 < alt <= 65617:
        temp = 389.99
    if alt > 65617:
        temp = 389.99 + 5.4864e-4 * (alt - 65617)
    return temp

P = 0.98
W = 11000
Cl_max = 1.24
M_max = 0.81
R = 1716.49
gamma = 1.4
q_max = 300
S = 232

altitude = np.arange(0, 65000, 500)
vel_max = []
dyn_press_vel = []
lift_vel = []
for alt in altitude:
    vel_max.append(M_max*np.sqrt(gamma*R*temp_calc(alt)))
    dyn_press_vel.append(np.sqrt(2*q_max/density_calc(alt)))
    lift_vel.append(np.sqrt((2*W)/(Cl_max*density_calc(alt)*S)))
    thrust =
    # drag_vel =
plt.plot(vel_max, altitude, label='Max Mach')
plt.plot(dyn_press_vel, altitude, label='Max Dynamic Pressure')
plt.plot(lift_vel, altitude, label='C_L Max')
plt.legend()
plt.xlim(0, 1000)
plt.ylim(0,60000)
plt.show()


