import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import space_functions as sf


def density(radius):
    if isinstance(radius, np.ndarray) or isinstance(radius, list):
        radius = np.linalg.norm(radius)
    url = 'https://raw.githubusercontent.com/ggb367/Spring-2020/master/366L/hw7/density.csv'
    altitude = radius - sf.earth.radius
    density = pd.read_csv(url)
    base_alt = np.array(density['Base Altitude'])
    scl_hgt = np.array(density['Scale Height'])
    nom_dens = np.multiply(np.array(density['Nominal Density']), 1000**3)
    del density
    idx = ((np.divide(np.abs(base_alt - altitude), 10)).astype(int)).argmin()  # return the base altitude index
    density = nom_dens[idx]*np.exp(-(altitude-base_alt[idx])/scl_hgt[idx])
    return density

altitude = np.arange(0, 1000, .5)

density_plot = np.array([density(alt) for alt in altitude])

plt.plot(density_plot, altitude)
plt.xscale('log')
plt.xlabel("Density [kg/km^3]")
plt.ylabel("Altitude [km]")
plt.title("U.S. Standard Atmospheric Model")
plt.grid(True)
plt.show()
