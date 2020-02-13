from skyfield.api import load
import numpy as np
import space_functions as sf

x_p = sf.deg2rad(0, seconds=0.41075)
y_p = sf.deg2rad(0, seconds=0.3378)
delta_ut1 = -0.200316
dX = sf.deg2rad(0, seconds=0.000146)
dY = sf.deg2rad(0, seconds=-0.000045)
X = sf.deg2rad(0, seconds=396.8722)
Y = sf.deg2rad(0, seconds=-1.481681)
s = sf.deg2rad(0, seconds=-0.000904)

ts = load.timescale()

time_utc = ts.utc(2020, 2, 14, 21)
time_jd_ut1 = time_utc.ut1
time_jd_tt = time_utc.tt
s_prime = sf.s_prime(time_jd_tt)
theta_ERA = sf.JD2ERA(time_jd_ut1)
print(time_jd_ut1)
print(theta_ERA*180/np.pi)


W = sf.polarMotion(x_p, y_p, s_prime)
print(W)
R = sf.R3(-theta_ERA)
print(R)
PN = sf.precession_nutation(X, Y, s)
print(PN)
Q = np.matmul(PN, np.matmul(R, W))
print(Q)

r_itrf = [-742.845, -5463.244, 3196.066]

r_gcrf_error = np.matmul(R, r_itrf)
print(r_gcrf_error)
r_gcrf_true = np.matmul(Q, r_itrf)
print(r_gcrf_true)
print(r_gcrf_error-r_gcrf_true)
print(np.linalg.norm(r_gcrf_error-r_gcrf_true))