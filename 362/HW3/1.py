import numpy as np
from scipy import optimize as opt

GAMMA = 1.4

M1 = float(input("What is M1?\n"))

second_check = True

while(second_check):
    second_var_type = input("What is your second input? Type: 'a' for Weak Turn Angle, 'b' for Strong Turn Angle, 'c' for Wave Angle, 'd' M1n)\n")
    if isinstance(second_var_type, str):
        second_var = float(input("Please enter the value for the second input(use degrees)\n"))
        second_check = False
    else:
        print("Please input a, b, c, or d")

if second_var_type is "a":
    # stregth = "weak"
    turn_angle = np.deg2rad(second_var)
    def theta_beta(beta):
        return 2 * (np.cos(beta)/np.sin(beta)) * ((M1 ** 2 * np.sin(beta) ** 2 - 1) / (M1 ** 2 * (GAMMA + np.cos(2 * beta)) + 2)) - np.tan((turn_angle))
    wave_angle = opt.fsolve(theta_beta, 15*np.pi/180)
    M1n = M1*np.sin(wave_angle)

elif second_var_type is "b":
    # strength = "strong"
    turn_angle = np.deg2rad(second_var)
    def theta_beta(beta):
        return (2*np.cos(beta)/np.sin(beta)) * ((M1 ** 2 * np.sin(beta) ** 2 - 1) / ((M1 ** 2 * (GAMMA + np.cos(2 * beta)) )+ 2)) - np.tan((turn_angle))
    wave_angle = opt.fsolve(theta_beta, 90*np.pi/180)
    M1n = M1 * np.sin(wave_angle)

elif second_var_type is "c":
    wave_angle = np.deg2rad(second_var)
    turn_angle = np.arctan2(2*(np.cos(wave_angle)/np.sin(wave_angle))*(((np.sin(wave_angle)**2)*(M1**2)-1)),((M1**2 * (1.4+np.cos(2*wave_angle)))+2))
    M1n = M1 * np.sin(wave_angle)

elif second_var_type is "d":
    M1n = second_var
    wave_angle = np.arcsin(M1n/M1)
    turn_angle = np.arctan2(2*(np.cos(wave_angle)/np.sin(wave_angle))*(((np.sin(wave_angle)**2)*(M1**2)-1)),((M1**2 * (1.4+np.cos(2*wave_angle)))+2))

M2n = np.sqrt((1+.5*(GAMMA-1)*M1n**2)/(GAMMA*M1n**2-.5*(GAMMA-1)))
rho_ratio = ((GAMMA+1)*M1n**2)/(2+(GAMMA-1)*M1n**2)
press_ratio = 1+(2*GAMMA)/(GAMMA+1)*(M1n**2-1)
enthalpy_ratio = press_ratio*rho_ratio**-1
stagnation_ratio = press_ratio*(enthalpy_ratio**-1)**(GAMMA/(GAMMA-1))
M2 = M2n/np.sin(wave_angle-turn_angle)
print("Turn Angle: ", np.rad2deg(turn_angle), "Wave Angle: ", np.rad2deg(wave_angle))
print("M2n: ", M2n, "M2: ", M2, "M1n: ", M1n)
print("Density Ratio: ", rho_ratio, "Pressure Ratio: ", press_ratio)
print("Stagnation Ratio: ", stagnation_ratio, "Enthalpy Ratio: ", enthalpy_ratio)
