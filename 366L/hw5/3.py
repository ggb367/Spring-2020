import space_functions as sf
import numpy as np

r, v = sf.PlanetRV(2451545)

print(np.matmul(sf.J20002GCRF(), r))