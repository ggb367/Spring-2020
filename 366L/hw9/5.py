import space_functions as sf
import numpy as np
import matplotlib.pyplot as plt

elements = [7000, 0, 75, 0, 0, 0]

rho_0 = [2.5e-2, 0, 0]
rho_dot_0 = [-1e-3, 0, 0]
# part a
print("5.1: ")
r_t, v_t = sf.elm2cart(elements, sf.earth.mu)

r_t_norm = np.linalg.norm(r_t)
R = np.divide(r_t, r_t_norm)
h = np.cross(r_t, v_t)
W = np.divide(h, np.linalg.norm(h))
S = np.cross(W, R)
Q_ijk = np.transpose(np.squeeze([R, S, W]))

rho_ijk = np.matmul(Q_ijk, rho_0)
rho_dot_ijk = np.matmul(Q_ijk, rho_dot_0)
r_c = r_t+rho_ijk
v_c = v_t+rho_dot_ijk
print(r_c, v_c)
# part b
print("5.2: ")
Ang_vel = h/r_t_norm**2

rho_dot_rel = np.multiply(1000, np.matmul(np.transpose(Q_ijk),rho_dot_0-np.cross(Ang_vel, rho_ijk)))
rho_rel = np.multiply(1000, np.matmul(np.transpose(Q_ijk),rho_0))
print(rho_dot_rel)

# part c

T_t = sf.period(elements[0])
elements_chaser = sf.cart2elm(r_c, v_c, sf.earth.mu)
T_c = sf.period(elements_chaser[0])

r_c_vec, v_c_vec = sf.orbit_prop_rk(r_c, v_c, 0, T_c*3, 20)
r_t_vec, v_t_vec = sf.orbit_prop_rk(r_t, v_t, 0, T_t*3, 20)



rho_truth_ijk = r_c_vec - r_t_vec
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot3D(r_c_vec[:, 0]-r_t_vec[:, 0], r_c_vec[:, 1]-r_t_vec[:, 1], r_c_vec[:, 2]-r_t_vec[:, 2])
# # ax.set_zlim(-2000, 2000)
# # ax.plot3D(r_t_vec[:, 0], r_t_vec[:, 1], r_t_vec[:, 2])
# # ax.plot3D(rho_truth[0, 0], rho_truth[0, 1], rho_truth[0, 2], 'g.')
# plt.figure()

rho_truth = np.empty(np.shape(rho_truth_ijk))

for row in range(np.shape(rho_truth)[0]):
    rho_truth[row, :] = np.multiply(1000, np.matmul(np.transpose(Q_ijk), rho_truth_ijk[row, :]))

plt.plot(rho_truth[:, 1], rho_truth[:, 0])
plt.plot(rho_truth[0, 1], rho_truth[0, 0], 'g.')  # start dot
plt.plot(rho_truth[-1, 1], rho_truth[-1, 0], 'r.')  # end dot
plt.grid(True)
plt.ylabel("Radial - R [m]")
plt.xlabel("Along-Track - S [m]")
plt.title("Position of Chaser in RSW Frame")

# part d
n = sf.mean_motion(T_t)
rho_cw, rho_dot_cw = sf.cw_prop(n, rho_rel, rho_dot_rel, 0, 3*T_t, 20)
# rho_cw, rho_dot_cw = sf.cw_prop(n, np.multiply(rho_0,1000), np.multiply(rho_dot_0,1000), 0, 3*T_t, 20)
error = rho_truth-rho_cw
# plt.figure()
# plt.plot(rho_cw[:, 1], rho_cw[:, 0])
# plt.figure()
# plt.plot(error[:, 1], error[:, 0])
fig, ax = plt.subplots(3, sharex=True)
ax[0].plot(error[:, 0])
ax[0].set_ylabel("Radial - R [m]")
ax[1].plot(error[:, 1])
ax[1].set_ylabel("Along-Track - S [m]")
ax[2].plot(error[:, 2])
ax[2].set_ylabel("Cross-Track - W [m]")
plt.suptitle("Error between CW and Numerical in RSW Frame")
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot3D(rho_truth[:, 0], rho_truth[:, 1], rho_truth[:, 2])
# ax.set_zlim(-2000, 2000)
# ax.plot3D(rho_truth[0, 0], rho_truth[0, 1], rho_truth[0, 2], 'g.')
# ax.plot3D(rho_truth[0, 0], rho_truth[0, 1], rho_truth[0, 2], 'g.')
print(5.5)
print(np.linalg.norm(error[-1, :]), "m")

plt.show()

