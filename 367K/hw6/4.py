import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve


def density_calc(alt):  # density calculator
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

def temp_calc(alt):  # temperature calculator
    temp = 0

    if alt <= 36089:
        temp = 518.69 - 3.5662e-3 * alt
    if 36089 < alt <= 65617:
        temp = 389.99
    if alt > 65617:
        temp = 389.99 + 5.4864e-4 * (alt - 65617)
    return temp
 #constants
P = 0.98
W = 11000
Cl_max = 1.24
M_max = 0.81
R = 1716.49
gamma = 1.4
q_max = 300
S = 232
C_do = 0.23
# allocation
altitude = np.arange(0, 65000, 500)
vel_max = []
dyn_press_vel = []
lift_vel = []
dens_trop = 7.0613e-4
p_s = 2116.2
mach_interp = [0.8, 0.9]
t_c1 = [1309, 1281]
t_c2 = [1921, 1894]
eta_interp = [0.9, 0.95]
n=0.925
K = 0.073
Thrust_corrected = np.interp(n, eta_interp, [np.interp(M_max, mach_interp, t_c1), np.interp(M_max, mach_interp, t_c2)])
thrust_drag_vel = []
thrust_drag_vel2 = []
a_half_TD = []
b_half_TD = []
alt_plot = []
T_t = 710
for alt in altitude:
    speed_sound = np.sqrt(gamma*R*temp_calc(alt))
    V_max = M_max*speed_sound
    V_cl = np.sqrt((2*W)/(density_calc(alt)*Cl_max*S))
    vel_max.append(V_max)
    dyn_press_vel.append(np.sqrt(2*q_max/density_calc(alt)))
    lift_vel.append(np.sqrt((2*W)/(Cl_max*density_calc(alt)*S)))
    # pressure = 2678.4 * np.e ** (alt * -4.8063e-5)
    # delta = pressure*(1+0.2*M_max**2)**3.5/p_s
    # T = Thrust_corrected * delta
    # Thrust.append(T)
    # def drag(V):
    #     L = .5 * Cl_max * V ** 2 * S * density_calc(alt)
    #     return T- .5 * C_do * density_calc(alt) * S * V ** 2 + K*L**2/(S*density_calc(alt)*V**2)
    # thrust_drag_vel2.append(fsolve(drag, 400))
    # thrust_drag_vel.append(np.sqrt((2*W)/(density_calc(alt)*S))*np.sqrt(K/C_do))
    if alt<36089:
        a = 1.2
    else:
        a = 1.0
    Thrust = 2*T_t*(density_calc(alt)/dens_trop)**a
    drag_min = 2*np.sqrt(C_do*K)*W
    if Thrust >= drag_min:
        a_half = lambda V : 0.5*C_do*density_calc(alt)*S*V**2 + ((2*K*W**2)/(density_calc(alt)*S*V**2)) - Thrust
        a_half_TD.append(float(fsolve(a_half, 10)))
        b_half_TD.append(float(fsolve(a_half, V_cl)))
        alt_plot.append(alt)
a_half_TD[a_half_TD==0] = np.nan
b_half_TD[b_half_TD==0] = np.nan
plt.plot(vel_max, altitude, label='Max Mach')
plt.plot(dyn_press_vel, altitude, label='Max Dynamic Pressure')
plt.plot(lift_vel, altitude, label='C_L Max')
plt.plot(a_half_TD, alt_plot, label='T-D=0')
plt.plot(b_half_TD, alt_plot, label='T-D=0')
plt.hlines(36089, 0, 1001, label='Tropopause')
# plt.plot(thrust_drag_vel2, altitude, label='T-D=0')
plt.legend(loc=2)
plt.xlim(0, 1000)
plt.ylim(0,60000)
plt.title("Flight Envelope of an ISBJ")
plt.xlabel("Velocity [ft/s]")
plt.ylabel("Altitude [ft]")
plt.show()


