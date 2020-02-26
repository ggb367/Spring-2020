import space_functions as sf
import numpy as np
import matplotlib.pyplot as plt

r, v = sf.PlanetRV(2451545)

print(np.matmul(sf.J20002GCRF(), r))

# # Part B
#
# time_series = np.arange(0, 365, 1)
#
# r_mat = np.empty([np.shape(time_series)[0], 3])
# v_mat = np.empty(np.shape(r_mat))
#
# count = 0
# print("Begining position propagation")
# for time_step in time_series:
#     r_mat[count, :], v_mat[count, :] = sf.PlanetRV(2457700+time_step)
#     r_mat[count, :] = np.matmul(sf.J20002GCRF(), r_mat[count, :])
#     r_mat[count, :] = np.divide(r_mat[count, :], 149597870.7)
#     count = count+1
#     print(count)
# print("complete")
# fig, ax = plt.subplots(np.shape(r_mat)[1], sharex=True)
# ax[0].plot(time_series, r_mat[:, 0])
# ax[1].plot(time_series, r_mat[:, 1])
# ax[2].plot(time_series, r_mat[:, 2])
# plt.suptitle("Planet RV Sun Position")
# plt.xlabel("Time [Days]")
# plt.ylabel("Distance [AU]")
#
# sun_series = np.empty([365,3])
# for j in range(np.size(sun_series, axis=0)):
#     sun_series[j, :] =  np.matmul(sf.MOD2GCRF(2457700+j), sf.sun_pos(2457700+j, AU=True))
#
# # fig, ax = plt.subplots(3, sharex=True)
# # ax[0].plot(range(365),sun_series[:, 0])
# # ax[0].set_title("Sun i pos")
# # ax[1].plot(range(365),sun_series[:, 1])
# # ax[1].set_title("Sun j pos")
# # ax[2].plot(range(365),sun_series[:, 2])
# # ax[2].set_title("Sun k pos")
# # plt.xlabel("Time [Days]")
# # ax[1].set_ylabel("Distance [AU]")
# # plt.suptitle("Position of Sun over 1 Year")
#
# # r_mat[:, 0] +
#
# fig, ax = plt.subplots(np.shape(r_mat)[1], sharex=True)
# ax[0].plot(time_series, r_mat[:, 0] + sun_series[:, 0])
# ax[1].plot(time_series, r_mat[:, 1] + sun_series[:, 1])
# ax[2].plot(time_series, r_mat[:, 2] + sun_series[:, 2])
# ax[0].plot(time_series, r_mat[:, 0])
# ax[1].plot(time_series, r_mat[:, 1])
# ax[2].plot(time_series, r_mat[:, 2])
# plt.suptitle("Difference in Planet RV and Sun Algorithm")
# plt.xlabel("Time [Days]")
# plt.ylabel("Distance [AU]")
#
# plt.show()
#
