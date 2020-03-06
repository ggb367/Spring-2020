import numpy as np
from matplotlib import pyplot as plt

Cdo = 0.023
K = 0.073
L = 11000/32.2
rho_s = 2.3769e-3
rho_t = 7.0613e-4
rho_plus = 1.7083e-4
h_t = 36089
h_plus = 65617
S = 54

velocity = np.arange(100, 900, 20)
altitude = [10000, 25000, 35000]

def density(alt):
    if alt <=36089:
        return rho_s*np.e**(-alt/29730)
    if alt > 36089 and alt <= 65617:
        return rho_t*np.e**(-(alt-h_t)/20806)
    if alt > 65617:
        return rho_plus*np.e**(-(alt-h_plus)/20770)

d = np.empty([np.size(velocity), 3])

d[:, 0] = np.multiply(.5*Cdo*density(altitude[0])*S, np.power(velocity, 2))+np.divide(2*K*L**2, np.multiply(S*density(altitude[0]), np.power(velocity, 2)))
d[:, 1] = np.multiply(.5*Cdo*density(altitude[1])*S , np.power(velocity, 2))+np.divide(2*K*L**2, np.multiply(S*density(altitude[1]), np.power(velocity, 2)))
d[:, 2] = np.multiply(.5*Cdo*density(altitude[2])*S , np.power(velocity, 2))+np.divide(2*K*L**2, np.multiply(S*density(altitude[2]), np.power(velocity, 2)))

plt.plot(velocity, d)

plt.figure()

d_p = np.empty([np.size(velocity), 3])
d_p[:, 0] = np.multiply(.5*Cdo*density(altitude[0])*S, np.power(velocity, 2))
d_p[:, 1] = np.multiply(.5*Cdo*density(altitude[1])*S, np.power(velocity, 2))
d_p[:, 2] = np.multiply(.5*Cdo*density(altitude[2])*S, np.power(velocity, 2))

plt.plot(velocity, d_p)

plt.figure()
d_i = np.empty([np.size(velocity), 3])
d_i[:, 0] = np.divide(2*K*L**2, np.multiply(density(altitude[0])*S, np.power(velocity, 2)))
d_i[:, 1] = np.divide(2*K*L**2, np.multiply(density(altitude[1])*S, np.power(velocity, 2)))
d_i[:, 2] = np.divide(2*K*L**2, np.multiply(density(altitude[2])*S, np.power(velocity, 2)))

plt.plot(velocity, d_i)

plt.show()
