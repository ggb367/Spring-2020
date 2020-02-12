import numpy as np
import space_functions as sf
from skyfield.api import load

ts = load.timescale()

obs_a_time_tt = ts.tt(jd=2458165.4375)
obs_a_time_utc = obs_a_time_tt._utc_float()

r_a = [2030.0, 18638.0, 3707.0]
v_a = [-4.5386, 0.52048, 0.10353]

obs_b_time_utc = 2458165.4375

time_diff = np.abs((obs_a_time_utc - obs_b_time_utc)*86400)  # seconds

r_a_f, v_a_f = sf.orbit_prop_rk(r_a, v_a, 0, time_diff, 0.001)

r_b = [1715.7, 18671.5, -3713.7]
v_b = [-4.54600, 0.446796, 0.088875]


print(r_a_f[-1, :]-r_b)

