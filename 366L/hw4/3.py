import numpy as np
from matplotlib import pyplot as plt

class moon:
    mu = 4.902e3
    distance = 384.4e3
class sun:
    mu = 1.237e11
    distance = 149.6e6

r = np.arange(6600, 50000, 10)
a_sun = np.empty(np.shape(r))
a_moon = np.empty(np.shape(r))
i = 0
for dist in r:
    try:
        a_sun[i] = sun.mu*((sun.distance-dist)/((sun.distance-dist)**3) - 1/sun.distance**2)
        a_moon[i] = moon.mu * ((moon.distance - dist) / ((moon.distance-dist) ** 3) - 1/ moon.distance ** 2)
    except RuntimeWarning as e:
        print(i)
        print(e)
    i = i + 1

plt.figure()
plt.plot(r, a_sun)
plt.plot(r, a_moon)
plt.title("Accelerations from the Sun and the Moon on a Satelite")
plt.legend(["Sun", "Moon"])
plt.xlabel("Distance from Earth [km]")
plt.ylabel("Acceleration [m/s]")

plt.show()
