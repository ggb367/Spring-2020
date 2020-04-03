import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def standard_atomosphere(radius):
    if isinstance(radius, np.ndarray) or isinstance(radius, list):
        radius = np.linalg.norm(radius)
    density = pd.read_csv('density.csv')
    base_alt = np.array(density['Base Altitude'])
    scl_hgt = np.array(density['Scale Height'])
    nom_dens = np.array(density['Nominal Density'])
    del density
    idx = ((np.divide(np.abs(base_alt - radius), 10)).astype(int)).argmin()  # return the base altitude index
    density = nom_dens[idx]*np.exp(-(radius-base_alt[idx])/scl_hgt[idx])
    return density

altitude = np.arange(0, 1000, .5)

density = np.array([standard_atomosphere(radius) for radius in altitude])

plt.plot(density, altitude)
plt.xscale('log')
plt.show()
