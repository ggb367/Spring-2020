import space_functions as sf
import numpy as np
import numpy.linalg as lg

mu = 1

r1 = np.array([2, 0, 0])
r2 = np.array([-3, 1, 0])

r1_norm = lg.norm(r1)
r2_norm = lg.norm(r2)

cos_nu = np.dot(r1, r2)/(r1_norm*r2_norm)

c = np.sqrt(r1_norm**2+r2_norm**2-2*r1_norm*r2_norm*cos_nu)
s = (r1_norm+r2_norm+c)/2
#part 1
a_min = s/2

print("a_min is : %f" % a_min)

#part 2
a = 3
f1, f2 = sf.Lambert_Focus_Finder(r1, r2, a)

c = lg.norm(f1)/2
print(c)
v = 2+c
e = c/v
print("The eccentricity is %f" % e)
