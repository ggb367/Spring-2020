import numpy as np

h = 30000
M = 0.7
W = 11000
S = 232
temp = 411.704
density = 0.0008892768161021751
a = np.sqrt(1.4*1716.5*temp)
V = M*a
C_L = (2*W)/(density*V**2*S)
viscosity = (2.27e-8*temp**(3/2))/(temp+198.6)

# [W, H, V, B, T, N] # Ignore wing nacelles
IF = [1.2, 1.1, 1.1, 1.2,  1.25, 1.5]
FF = [1+1.6*.09+100*.09**4, 1+1.6*.08+100*.08**4, 1+1.6*.1+100*.1**4, 1+60/(41/5.250)**3+0.0025*(41/5.25), 1+60/(14*1.75)**3+.0025*(14/1.75), 1+0.35/(7.7/2.3)]
S_wet = [344, 108, 75.4, 456, 61.2, 55.6]
l = [7, 3.83, 6.92, 41, 14, 7.7]

# there are two wings,  tip tanks and nacelles
steps = [0, 0, 1, 2, 3, 4, 4, 5, 5]
f = 0
for step in steps:
    Re = density*V*l[step]/viscosity
    C_f = 0.455/(np.log10(Re))**2.58
    CF = (1+0.2*M**2)**-0.467
    f_k = C_f*CF*IF[step]*FF[step]*S_wet[step]
    f = f+f_k
C_D_f = 1.1*f/S
print("This is C_D_f: ", C_D_f)

