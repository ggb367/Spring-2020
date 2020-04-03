import numpy as np

def temp_calc(alt):
    temp = 0

    if alt <= 36089:
        temp = 518.69 - 3.5662e-3 * alt
    if 36089 < alt <= 65617:
        temp = 389.99
    if alt > 65617:
        temp = 389.99 + 5.4864e-4 * (alt - 65617)
    return temp


ro=0.00073654
CD0=0.022332
K=0.073
CD=0.02911
R=1716.5
alt=35000
W=11000
V=(650/5280)*3600
#Values found from HW 5 Problem 5
n=0.925 #eta
T=1050.8 #thrust

a = np.sqrt(1.4 * R * temp_calc(alt))
Mach = V/a
mach_interp = [0.6, 0.7]
ts=518.69
t_=temp_calc(alt)*(1+0.2*Mach**2)
theta=t_/ts
ma=[0.6, 0.7]
eta=[0.9, 0.95]
C_c = np.interp(n, eta, [np.interp(Mach, ma, [1.184, 1.190]), np.interp(Mach, ma, [1.127, 1.129])])
C=C_c*np.sqrt(theta);
F=V/(C*T)
G=1/(C*T)
print("This is F: ", F)
print("This is G: ", G)
