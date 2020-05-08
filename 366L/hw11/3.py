import space_functions as sf
import numpy as np
import numpy.linalg as lg

r_earth = [1, 0, 0]
r_mars = [0, 1.524, 0]
mu_sun = 1
t_m = 1

r_earth_norm = lg.norm(r_earth)
r_mars_norm = lg.norm(r_mars)

cos_nu = np.dot(r_earth, r_mars)/(r_earth_norm*r_mars_norm)
sin_nu = t_m*np.sqrt(1-cos_nu**2)

c = np.sqrt(r_earth_norm**2+r_mars_norm**2-2*r_earth_norm*r_mars_norm*cos_nu)
s = (r_earth_norm+r_mars_norm+c)/2

t_abs = (1/3)*np.sqrt(2/mu_sun)*(s**1.5-(s-c)**1.5)
print("Absolute minimum time is %f TU" % t_abs)

a_min = s/2
print("a_min is : %f AU" % a_min)

p_min = (r_earth_norm*r_mars_norm/c)*(1-cos_nu)

v_earth = np.multiply(np.sqrt(mu_sun*p_min)/(r_earth_norm*r_mars_norm*sin_nu), r_mars - np.multiply((r_mars_norm/p_min)*(1-cos_nu), r_earth))
print("v_earth: ")
print(v_earth)

B_min = 2*np.arcsin(np.sqrt((s-c)/s))
t_of = np.sqrt(a_min**3/mu_sun)*(np.pi-t_m*(B_min-np.sin(B_min)))

print("Minimum Time of Flight is %f TU" % t_of)
