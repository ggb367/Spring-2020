import space_functions as sf
import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt

r_earth = [1, 0, 0]
_, v_earth = sf.elm2cart([1, 0, 0, 0, 0, 0], 1)
print(v_earth)
t = []
C3 = []
u = np.arange(5, 170, 5)
for step in u:
    r_mars, v_mars = sf.elm2cart([1.524, 0, 1.85, 0, 0, step], 1)
    mu_sun = 1
    t_m = 1
    a_min, e_min, t_abs, v_0 = sf.Lambert_MinEng(r_earth, r_mars)
    t.append(t_abs)
    v_inf = lg.norm(v_0-v_earth)
    C3.append(v_inf**2)

plt.plot(u, t)
plt.xlabel("u_mars, deg")
plt.ylabel("Time [TU]")
plt.title("Argument of Latitude of Mars vs. Absolute Minimum Time of Flight")
plt.figure()
plt.xlabel("u_mars, deg")
plt.ylabel("C_3 [AU^2/TU^2]")
plt.title("Argument of Latitude of Mars vs. Departure C_3")
plt.plot(u, C3)
# plt.yscale('log')
plt.show()

# TODO: Algorithm 57
# print("Absolute minimum time is %f TU" % t_abs)
