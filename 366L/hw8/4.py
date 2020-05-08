import space_functions as sf
import numpy as np
import matplotlib.pyplot as plt

sun_pos = [149597870, 0, 0]  # km
elements = [10000, 0, 90, 0, 0, 0]
period = sf.period(elements[0])
# propogate orbit
r_0, v_0 = sf.elm2cart(elements, sf.earth.mu)
r_vec, v_vec = sf.orbit_prop_rk(r_0, v_0, 0, period, 3)
#get when it's in the shadow
shadow = np.array([sf.cylindrical_shadow(step, sun_pos) for step in r_vec])
plt.plot(shadow)

in_shadow = np.squeeze(np.where(shadow==0))
# print(in_shadow) #leaves shadow at step 1294, enters at step 2023

leave_shadow = sf.cart2elm(r_vec[2023, :], v_vec[2023, :], sf.earth.mu)
enter_shadow = sf.cart2elm(r_vec[1294, :], v_vec[1294, :], sf.earth.mu)

print("Enters Shadow at u=", enter_shadow[-1], "°")
print("Leaves Shadow at u=", leave_shadow[-1], "°")
plt.show()

