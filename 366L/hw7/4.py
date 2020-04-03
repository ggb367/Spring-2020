import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import space_functions as sf

RE = 6378.1363
MU = 398600.4415
J2 = 0.0010826267
J3 = -0.0000025327

class conditions:
    p_srp = 4.57e-6 * 1000 ** 2
    C_r = 1.5
    A_m = 0.1 / 1000 ** 2
    C_D = 2
    url = 'https://raw.githubusercontent.com/ggb367/Spring-2020/master/366L/hw7/density.csv'
    density_table = pd.read_csv(url)
    epoch = 2458200.5  # days
# epoch = sf.JD2MJD(epoch)

r_0, v_0 = sf.elm2cart([6800, 0.005, 71, 300, 78, 0], MU)
r_vec, v_vec = sf.orbit_prop_all_pert(r_0, v_0, 0, 86400, 30, conditions)
r_bad, v_bad = sf.orbit_prop_rk(r_0, v_0, 0, 86400, 30)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot3D(r_vec[:, 0], r_vec[:, 1], r_vec[:, 2])
ax.plot3D(r_bad[:, 0], r_bad[:, 1], r_bad[:, 2])

print("With Perts:")
print("The Final Position is: ", r_vec[-1, :], "The Final Velocity is: ", v_vec[-1, :])
print("Without Perts:")
print("The Final Position is: ", r_bad[-1, :], "The Final Velocity is: ", v_bad[-1, :])


time_series = np.arange(0, 24, 24/2880)
fig, ax = plt.subplots(3, sharex=True)
ax[0].plot(time_series, r_vec[:, 0]-r_bad[:, 0])
ax[1].plot(time_series, r_vec[:, 1]-r_bad[:, 1])
ax[2].plot(time_series, r_vec[:, 2]-r_bad[:, 2])


plt.show()
