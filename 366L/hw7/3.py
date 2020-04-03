import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def standard_atomosphere(radius, density):
    if isinstance(radius, np.ndarray) or isinstance(radius, list):
        radius = np.linalg.norm(radius)
    # density = pd.read_csv('density.csv')
    base_alt = np.array(density['Base Altitude'])
    scl_hgt = np.array(density['Scale Height'])
    nom_dens = np.multiply(np.array(density['Nominal Density']), 1000**3)
    del density
    idx = ((np.divide(np.abs(base_alt - radius), 10)).astype(int)).argmin()  # return the base altitude index
    density = nom_dens[idx]*np.exp(-(radius-base_alt[idx])/scl_hgt[idx])
    return density

altitude = np.arange(0, 1000, .5)
url = 'https://raw.githubusercontent.com/ggb367/Spring-2020/master/366L/hw7/density.csv'
density = pd.read_csv(url)
density_plot = np.array([standard_atomosphere(radius, density) for radius in altitude])

plt.plot(density_plot, altitude)
plt.xscale('log')
plt.xlabel("Density [kg/km^3]")
plt.ylabel("Altitude [km]")
plt.title("U.S. Standard Atmospheric Model")
plt.grid(True)
plt.show()
