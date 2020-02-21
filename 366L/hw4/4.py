import numpy as np
import space_functions as sf
from matplotlib import pyplot as plt
from skyfield.api import load

ts = load.timescale()

UTC_time = ts.utc(2000,1, 1, 12)

test = sf.sun_pos(2457793.5, AU=True)  # This is in MOD frame
test_GCRF = np.matmul(sf.MOD2GCRF(2457793.5), test)
# test = sf.sun_pos(2451545, AU=True)  # This is in MOD frame
# test_GCRF = np.matmul(sf.MOD2GRCF(2451545), test)
print(test_GCRF)
# print(sf.MOD2GRCF(2457793.5))

#part b
sun_series = np.empty([365,3])
for j in range(np.size(sun_series, axis=0)):
    sun_series[j, :] =  np.matmul(np.transpose(sf.J20002GCRF()), np.matmul(sf.MOD2GCRF(2457793.5+j), sf.sun_pos(2457793.5+j, AU=True)))

# plt.plot(range(365), sun_series)
fig, ax = plt.subplots(3, sharex=True)
ax[0].plot(range(365),sun_series[:, 0])
ax[0].set_title("Sun i pos")
ax[1].plot(range(365),sun_series[:, 1])
ax[1].set_title("Sun j pos")
ax[2].plot(range(365),sun_series[:, 2])
ax[2].set_title("Sun k pos")
plt.xlabel("Time [Days]")
ax[1].set_ylabel("Distance [AU]")
plt.suptitle("Position of Sun over 1 Year")
plt.show()
