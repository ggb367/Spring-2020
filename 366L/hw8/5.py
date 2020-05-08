import space_functions as sf
import numpy as np
import numpy.linalg as lg

sun_pos = np.array([149597870, 0, 0])  # km
# cases:
# r = np.array([-7000, -6400, 600])
r = np.array([-1000, 4695, 4316])

a = np.arcsin(sf.sun.radius/lg.norm(sun_pos+r))
b = np.arcsin(sf.earth.radius/lg.norm(r))
c = 1.415 #np.arccos(np.dot(r, (sun_pos+r))/(lg.norm(r)*lg.norm(sun_pos+r)))

print(a, b, c)
x = (c**2+a**2-b**2)/(2*c)
print(x)
y = np.sqrt(a**2-x**2)
print(y)
A = a**2*np.arccos(x/a)+b**2*np.arccos((c-x)/b)-c*y
print(A)

if c<np.abs(a-b):
    gamma = 0
elif a+b<=c:
    gamma = 1
else:
    gamma = 1-A/(np.pi*a**2)

print(gamma)
