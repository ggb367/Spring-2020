import numpy as np

alt = float(input("What is the altitude in feet?\n"))
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

print("The density is: ", density_true)
