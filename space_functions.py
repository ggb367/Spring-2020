import math as m
import warnings

import numpy as np
import numpy.linalg as lg
from scipy.integrate import ode
from scipy.optimize import fsolve
from colorama import Fore


def AU2km(AU):
    return AU*149597870.7

class earth:
    mu = 398600.4418
    semimajor = AU2km(1.000001018)
    p_srp = 4.57e-9
    J2 = 0.10826267e-2
    J3 = -0.2532327e-5
    radius = 6378.1363  # km
    rot_speed = 0.04651  # km/s
    def eccentricty(T_TDB):
        return 0.01670862-0.000042037*T_TDB-0.0000001236*T_TDB**2+0.00000000004*T_TDB**3

    def inclination(T_TDB, deg=False):
        if deg:
            return 0+0.0130546*T_TDB-0.00000931*T_TDB**2-0.000000034*T_TDB**3
        else:
            return np.deg2rad(0+0.0130546*T_TDB-0.00000931*T_TDB**2-0.000000034*T_TDB**3)
    def RAAN(T_TDB, deg=False):
        if deg:
            return 174.873174-0.2410908*T_TDB+0.00004067*T_TDB**2-0.000001327*T_TDB**3
        else:
            return np.deg2rad(174.873174-0.2410908*T_TDB+0.00004067*T_TDB**2-0.000001327*T_TDB**3)
    def ARG_PERIHELION(T_TDB, deg=False):
        if deg:
            return 102.937348+0.3225557*T_TDB+0.00015026*T_TDB**2+0.000000478*T_TDB**3
        else:
            return np.deg2rad(102.937348+0.3225557*T_TDB+0.00015026*T_TDB**2+0.000000478*T_TDB**3)
    def Mean_Long(T_TDB, deg=False):
        if deg:
            return 100.466449+35999.3728519*T_TDB-0.00000568*T_TDB**2
        else:
            return np.deg2rad(100.466449+35999.3728519*T_TDB-0.00000568*T_TDB**2)
    def obliquity(T_TT, deg=False):
        if deg:
            return 23.439279-0.0130102*T_TT-5.086e-8*T_TT**2+5.565e-7*T_TT**3+1.6e-10*T_TT**4+1.21e-11*T_TT**5
        else:
            return  np.deg2rad(23.439279-0.0130102*T_TT-5.086e-8*T_TT**2+5.565e-7*T_TT**3+1.6e-10*T_TT**4+1.21e-11*T_TT**5)


class sun:
    mu = 1.32712440042e11
    radius = 6.957e5

class moon:
    mu = 3903


def orbit_prop(time_series, mean_motion, eccent, time_periapsis):  # propogate an eliptical orbit
    # allocate memory for anomalies
    E = np.empty(np.size(time_series))
    M = np.empty(np.size(time_series))
    nu = np.empty(np.size(time_series))
    eccent_norm = lg.norm(eccent)
    # Helper Functions
    def f(x, m_e):
        return x - eccent_norm * np.sin(x) - m_e
    def df(x):
        return 1 - eccent_norm * np.cos(x)
    # propagate through the time series
    for i in range(np.size(time_series)):
        M[i] = mean_motion*(time_series[i]-time_periapsis)  # mean anomaly for this time step
        if M[i] < m.pi:  # inital guess based on mean anomaly
            guess = M[i]+eccent_norm/2
        else:
            guess = M[i]-eccent_norm/2
        it = 0
        error = 100.0
        while error > 10**-10 and it <= 50:  # newton raphson to find eccentric anomaly
            # try:
                E[i] = guess-f(guess, M[i])/df(guess)
                error = np.abs((E[i]-guess)/E[i])
                guess = E[i]
                it = it+1
            # except(ZeroDivisionError, RuntimeWarning, RuntimeError):
            #     print(E[i])
            #     print("Zero Division error")
        nu[i] = 2*m.atan2(np.sqrt(1+eccent_norm)*np.tan(E[i]/2), np.sqrt(1-eccent_norm))  # find anomaly from eccentric anomaly
    return nu, E, M

# class CW:


def hyper_orbit_prop(time_series, mean_motion, eccent, time_periapsis):  # propogate a hyperbolic orbit
    # allocate memory for anomalies
    F = np.empty(np.size(time_series))
    M = np.empty(np.size(time_series))
    nu = np.empty(np.size(time_series))
    eccent_norm = lg.norm(eccent)
    # Helper Functions
    def f(x, m_h):
        return -x+eccent_norm*np.sinh(x)-m_h
    def df(x):
        return -1+eccent_norm*np.cosh(x)
    # propagate through the time series
    for i in range(np.size(time_series)):
        M[i] = mean_motion*(time_series[i]-time_periapsis)  # mean anomaly for this time step
        if M[i] < m.pi:  # inital guess based on mean anomaly
            guess = M[i]+eccent_norm/2
        else:
            guess = M[i]-eccent_norm/2
        it = 0
        error = 100.0
        while error > 10**-10 and it <= 50:  # newton raphson to find eccentric anomaly
            F[i] = guess-f(guess, M[i])/df(guess)
            error = np.abs((F[i]-guess)/F[i])
            guess = F[i]
            it = it+1
        nu[i] = 2*m.atan2(np.sqrt(eccent_norm+1)*np.tanh(F[i]/2), np.sqrt(eccent_norm-1))  # find anomaly from eccentric anomaly
    return nu, F, M


def cart2elm(r, v, mu, deg=True):  # transform position and velocity to classical orbital elements
    h = np.cross(r, v)
    r_norm = lg.norm(r)
    v_norm = lg.norm(v)
    eccent = np.cross(v, h) / mu - np.divide(r, r_norm)  # eccentricity
    eccent_norm = lg.norm(eccent)
    energy = (v_norm**2)/2 - mu/r_norm
    h_norm = lg.norm(h)
    k = (h_norm ** 2) / (r_norm * mu) - 1
    if energy < 0:
        a = -mu/(2*energy)
    elif -10e-12 < energy < 10e-12:
        a = m.inf
    else:
        a = mu/(2*energy)
    i = np.arccos(np.dot(h, [0, 0, 1])/h_norm)
    n = np.cross([0, 0, 1], h)
    n_norm = lg.norm(n)
    if eccent_norm < 10e-12 or eccent_norm > 10e-12:
        nu = np.arccos(k/eccent_norm)
        if np.dot(r,v)<0:
            nu = 2*m.pi-nu
        RAAN = np.arccos(np.dot(n, [1, 0, 0])/n_norm)
        omega = np.arccos(np.dot(n, eccent)/(eccent_norm*n_norm))
    if eccent_norm < 10e-12 and i < 10e-12:
        RAAN = 0
        omega = 0
        nu = np.arccos(r[1]/r_norm)
        if r[1] < 0:
            nu = 2*m.pi-nu
    elif eccent_norm < 10e-12:
        omega = 0
        RAAN = np.arccos(np.dot(n, [1, 0, 0]) / n_norm)
        nu = np.arccos(np.dot((n/n_norm),r)/r_norm)
        if r[2]< 0:
            nu = 2*m.pi-nu
    elif i < 10e-12:
        RAAN = 0
        omega = np.arccos(np.dot(eccent, [1, 0, 0])/eccent_norm)
        if e[1]< 0:
            omega = 2*m.pi-omega
    if deg:
        nu = 180*nu/m.pi
        i = 180*i/m.pi
        RAAN = 180*RAAN/m.pi
        omega = 180*omega/m.pi
    E = [a, eccent_norm, i, RAAN, omega, nu]
    for element in E:
        if not isinstance(element, float):
            print(E)
            raise TypeError("One of the elements is not a float!")
    return np.array(E)


def elm2cart(E, mu, deg=True):  # transform classical orbital elements to cartesian position and velocity
    # E - [a, e, i, RAAN, omega, nu]
    a = E[0]
    e = E[1]
    if deg:
        i = m.pi * E[2] / 180
        RAAN = m.pi * E[3] / 180
        omega = m.pi * E[4] / 180
        nu = m.pi * E[5] / 180
    else:
        i = E[2]
        RAAN = E[3]
        omega = E[4]
        nu = E[5]
    p = a*(1 - e**2)
    r_pqw = np.array([(p/(1+e*np.cos(nu)))*np.cos(nu), (p/(1+e*np.cos(nu)))*np.sin(nu), 0])
    v_pqw = np.array([np.sqrt(mu/p)*(-np.sin(nu)), np.sqrt(mu/p)*(e+np.cos(nu)), 0])
    # R_3(-RAAN)R_1(-i)R_3(-omega)
    c1 = np.cos(-omega)
    c2 = np.cos(-i)
    c3 = np.cos(-RAAN)
    s1 = np.sin(-omega)
    s2 = np.sin(-i)
    s3 = np.sin(-RAAN)
    q1 = np.array([c1*c3-c2*s1*s3, c3*s1+c1*c2*s3, s3*s2])
    q2 = np.array([-c1*s3-c3*c2*s1, c1*c2*c3-s1*s3, c1*s2])
    q3 = np.array([s1*s2, -c1*s2, c2])
    Q = np.array([q1, q2, q3])
    r = np.matmul(Q, r_pqw)
    v = np.matmul(Q, v_pqw)
    return r, v

def R1(phi): # returns R1 transform matrix
    return np.array([np.array([1, 0, 0 ]), np.array([0, np.cos(phi), np.sin(phi)]), np.array([0, -np.sin(phi), np.cos(phi)])])

def R2(phi): # returns R2 transform matrix
    return np.array([np.array([np.cos(phi), 0, -np.sin(phi)]), np.array([0, 1, 0]), np.array([np.sin(phi), 0, np.cos(phi)])])

def R3(phi):  # returns R3 transform matrix
    return np.array([np.array([np.cos(phi), np.sin(phi), 0]), np.array([-np.sin(phi), np.cos(phi), 0]), np.array([0, 0, 1])])

def deg2rad(degree, minutes=0, seconds=0):  # transform an array of degrees to radians
    if isinstance(degree, int) or isinstance(degree, float):
        return (degree+minutes/60+seconds/3600)*np.pi/180
    elif isinstance(degree, list) or isinstance(degree, np.ndarray):
        output = np.empty(np.size(input))
        for i in range(np.size(input)):
            output[i] = input[i] * np.pi / 180
        return output
    else:
        raise TypeError("degree must be a int, float, list, or ndarray, you used a %s", str(type(degree)))

def orbit_prop_rk(r_0, v_0, T0, tF, dT):  # propogate an orbit about Earth using Runge-Kutta Method
    def two_body_orbit(t, Y, mu):
        dY = np.empty([6, 1])
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]
        r = np.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)
        dY[3] = -mu * Y[0] / r ** 3
        dY[4] = -mu * Y[1] / r ** 3
        dY[5] = -mu * Y[2] / r ** 3
        return dY

    MU = 398600.4415

    def derivFcn(t, y):
        return two_body_orbit(t, y, MU)

    Y_0 = np.concatenate([r_0, v_0], axis=0)
    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y_0, T0)
    output = []
    output.append(np.insert(Y_0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:  # rv.successful() and
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.array(output)
    t = output[:, 0]

    r_vec = np.empty([np.shape(output)[0]-1, 3])
    v_vec = np.empty([np.shape(output)[0]-1, 3])

    for i in range(np.shape(output)[0]-1):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
    return r_vec, v_vec

def CRTBP_prop_rk(r_0, v_0, T0, tF, dT, MU):  # propogate an orbit in the CRTBP frame
    def CRTBP_orbit(t, Y, mu):
        dY = np.empty([6, 1])
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]
        r1 = np.sqrt((Y[0]+mu)**2+Y[1]**2+Y[2]**2)
        r2 = np.sqrt((Y[0]+mu-1)**2+Y[1]**2+Y[2]**2)
        dY[3] = 2*dY[1]+Y[0]-(1-mu)*(Y[0]+mu)/r1**3-mu*(Y[0]+mu-1)/r2**3
        dY[4] = -2*dY[0] + Y[1]-(1-mu)*Y[1]/r1**3-mu*Y[1]/r2**3
        dY[5] = -(1-mu)*Y[2]/r1**3-mu*Y[2]/r2**3
        return dY

    def derivFcn(t, y):
        return CRTBP_orbit(t, y, MU)

    Y_0 = np.concatenate([r_0, v_0], axis=0)
    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y_0, T0)
    output = []
    output.append(np.insert(Y_0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:  # rv.successful() and
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.array(output)
    t = output[:, 0]

    r_vec = np.empty([np.shape(output)[0], 3])
    v_vec = np.empty([np.shape(output)[0], 3])

    for i in range(np.shape(output)[0]):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
    return r_vec, v_vec

def lagrange(mu):  # returns a 2x5 vector of lagrange points given mu
    f = lambda r_x: r_x - (1 - mu) * (r_x + mu) / np.abs(r_x + mu) ** 3 - mu * (r_x - (1 - mu)) / np.abs(
        r_x + mu - 1) ** 3

    r_x = np.array([-1, 0, 1])
    r_0_roots = np.array(fsolve(f, r_x))
    roots_x = np.append(r_0_roots, [.5 - mu, .5 - mu])
    roots_y = np.array([0, 0, 0, np.sqrt(3) / 2, -np.sqrt(3) / 2])
    points = np.column_stack((roots_x, roots_y))
    return points

def rad2deg(x):
    return x*180/m.pi

def Date2JD(year, month, day, hour, minute, second):
    return 367*year-np.floor((7*(year+np.floor((month+9)/12)))/4)+np.floor(275*month/9)+day+1721013.5+(1/24)*(hour+(1/60)*(minute+second/60))

def JD2MJD(JD):
    return JD - 2400000.5

def JD2DOY(JulianDate):
    T_1900 = (JulianDate-2415019.5)/365.25
    year = 1900+np.trunc(T_1900)
    LeapYears = np.trunc((year-1901)*.25)
    Days = (JulianDate-2415019.5)-((year-1900)*365+LeapYears)
    if Days < 1.0:
        year = year-1
        LeapYears = np.trunc((year-1901)*25)
        Days = (JulianDate- 2415019.5)-((year-1900)*365+LeapYears)
    return np.trunc(Days)

def time2radians(hour, minute, seconds):
    return 15*(hour+minute/60+seconds/3600)*np.pi/180

def JD2ERA(JulianDate):
    return np.mod(2*np.pi*(0.779057273264+1.00273781191135448*(JulianDate-2451545)), 2*np.pi)

def polarMotion(x_p, y_p, s_prime):
    return np.matmul(R3(-s_prime), np.matmul(R2(x_p),R1(y_p)))

def JulianCenturies(JulianDate):
    return (JulianDate-2451545)/36525

def MJDCenturies(MJD):
    return (MJD-51544.5)/36525

def s_prime(centuries_tt):
    return deg2rad(0, seconds=-0.000047 * centuries_tt)

def precession_nutation(X, Y, s):
    a = (0.5+0.125*(X**2+Y**2))
    return np.matmul(np.array([[1-a*X**2, -a*X*Y, X], [-a*X*Y, 1-a*Y**2, Y], [-X, -Y, 1-a*(X**2+Y**2)]]), R3(s))

def sun_pos(JulianDate, AU=False):
    JulianDate = JulianDate-2400000.5
    T = MJDCenturies(JulianDate)
    longitude_sun = deg2rad(280.46+36000.771*T)
    M = deg2rad(357.52772333 + 35999.0534*T)
    longitude_ecliptic = longitude_sun + np.deg2rad(1.914666471*np.sin(M)+0.019994643*np.sin(2*M))
    radius = 1.000140612-0.016708617*np.cos(M) - 0.000139589*np.cos(2*M)
    obliquity = deg2rad(23.439291-0.0130042*T)
    # print([T, longitude_sun, M, longitude_ecliptic, radius, ecliptic])
    if AU:
        return np.array([radius*np.cos(longitude_ecliptic), radius*np.cos(obliquity)*np.sin(longitude_ecliptic), radius*np.sin(obliquity)*np.sin(longitude_ecliptic)])
    else:
        return np.multiply(np.array([radius*np.cos(longitude_ecliptic), radius*np.cos(obliquity)*np.sin(longitude_ecliptic), radius*np.sin(obliquity)*np.sin(longitude_ecliptic)]),149597870)

def MOD2GCRF(Julian_Date):
    Julian_Date = Julian_Date-2400000.5
    JCTT = MJDCenturies(Julian_Date)
    zeta = deg2rad(0, seconds=2306.2181*JCTT+0.30188*JCTT**2+0.017998*JCTT**3)
    theta = deg2rad(0, seconds=2004.3109*JCTT-0.42665*JCTT**2-0.041833*JCTT**3)
    z = deg2rad(0, seconds=2306.2181*JCTT+1.09468*JCTT**2+0.018203*JCTT**3)
    return np.matmul(R3(zeta),np.matmul(R2(-theta), R3(z)))

def J20002GCRF():
    delta = deg2rad(0, seconds=0.0146)
    zeta = deg2rad(0, seconds=-0.16617)
    eta = deg2rad(0, seconds=-0.0068192)
    return np.matmul(R3(-delta),np.matmul(R2(-zeta), R1(eta)))

def orbit_prop_3body(r_0, v_0, T0, tF, dT):

    def three_body_orbit(t, Y, mu):
        dY = np.empty([6, 1])
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]
        r = lg.norm(Y[0:3])
        t = t/86400+2451545
        sun_range = np.matmul(MOD2GCRF(t),sun_pos(t))
        sat2sun = sun_range - Y[0:3]
        sat2sun_norm = lg.norm(sat2sun)
        sun_range_norm = lg.norm(sun_range)
        dY[3] = (-mu * Y[0] / r ** 3) + sun.mu*(sat2sun[0]/((sat2sun_norm)**3)-sun_range[0]/sun_range_norm**3)
        dY[4] = (-mu * Y[1] / r ** 3) + sun.mu*(sat2sun[1]/((sat2sun_norm)**3)-sun_range[1]/sun_range_norm**3)
        dY[5] = (-mu * Y[2] / r ** 3) + sun.mu*(sat2sun[2]/((sat2sun_norm)**3)-sun_range[2]/sun_range_norm**3)
        return dY


    def derivFcn(t, y):
        return three_body_orbit(t, y, earth.mu)

    Y_0 = np.concatenate([r_0, v_0], axis=0)
    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y_0, T0)
    output = []
    output.append(np.insert(Y_0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:  # rv.successful() and
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    if not rv.successful() and rv.t<tF:
        warnings.warn("Runge Kutta Failed!", RuntimeWarning)
    #  Convert the output a numpy array for later use
    output = np.array(output)
    t = output[:, 0]

    r_vec = np.empty([np.shape(output)[0]-1, 3])
    v_vec = np.empty([np.shape(output)[0]-1, 3])

    for i in range(np.shape(output)[0]-1):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
    return r_vec, v_vec

def KepEqtnE(M, e):
    if -np.pi < M < 0 or M>np.pi:
        E = M-e
    else:
        E = M+e
    E_old = E
    count = 0
    while(count<10e4):
        E = E_old+(M-E_old+e*np.sin(E_old))/(1-e*np.cos(E_old))
        count = count+1
        if (abs(E - E_old) < 10e-6):
            break
        E_old = E
    return  E

def PlanetRV(JD_TDB, MJD=False):
    if not MJD:
        JD_TDB = JD_TDB - 2400000.5
    T_TDB = MJDCenturies(JD_TDB)
    M = earth.Mean_Long(T_TDB) - earth.ARG_PERIHELION(T_TDB)
    arg_periapsis = (earth.ARG_PERIHELION(T_TDB) - earth.RAAN(T_TDB))
    eccentric_anomaly = KepEqtnE(M, earth.eccentricty(T_TDB))
    # elements - a e i RAAN arg peri nu
    nu = 2 * m.atan2(np.sqrt(1 + earth.eccentricty(T_TDB)) * np.tan(eccentric_anomaly / 2), np.sqrt(1 - earth.eccentricty(T_TDB)))
    r, v = elm2cart([earth.semimajor, earth.eccentricty(T_TDB), earth.inclination(T_TDB), earth.RAAN(T_TDB), arg_periapsis, nu], sun.mu, deg=False)
    r = np.matmul(R1(-earth.obliquity(T_TDB)), r)
    v = np.matmul(R1(-earth.obliquity(T_TDB)), v)
    return  r, v

def orbit_prop_3body_RV(r_0, v_0, T0, tF, dT):

    def three_body_orbit(t, Y, mu):
        dY = np.empty([6, 1])
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]
        r = lg.norm(Y[0:3])
        t = t/86400+2451545
        sun_range, _ = PlanetRV(t)
        sun_range = np.matmul(J20002GCRF(), sun_range)
        sun_range = np.multiply(sun_range, -1)
        sat2sun = sun_range - Y[0:3]
        sat2sun_norm = lg.norm(sat2sun)
        sun_range_norm = lg.norm(sun_range)
        dY[3] = (-mu * Y[0] / r ** 3) + sun.mu*(sat2sun[0]/((sat2sun_norm)**3)-sun_range[0]/sun_range_norm**3)
        dY[4] = (-mu * Y[1] / r ** 3) + sun.mu*(sat2sun[1]/((sat2sun_norm)**3)-sun_range[1]/sun_range_norm**3)
        dY[5] = (-mu * Y[2] / r ** 3) + sun.mu*(sat2sun[2]/((sat2sun_norm)**3)-sun_range[2]/sun_range_norm**3)
        return dY


    def derivFcn(t, y):
        return three_body_orbit(t, y, earth.mu)

    Y_0 = np.concatenate([r_0, v_0], axis=0)
    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y_0, T0)
    output = []
    output.append(np.insert(Y_0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:  # rv.successful() and
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    if not rv.successful() and rv.t<tF:
        warnings.warn("Runge Kutta Failed!", RuntimeWarning)
    #  Convert the output a numpy array for later use
    output = np.array(output)
    t = output[:, 0]

    r_vec = np.empty([np.shape(output)[0]-1, 3])
    v_vec = np.empty([np.shape(output)[0]-1, 3])

    for i in range(np.shape(output)[0]-1):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
    return r_vec, v_vec


def J2J3_Pert(r):
    r_norm = lg.norm(r)
    a_2 = np.multiply(3 * earth.mu * earth.J2 * earth.radius ** 2 / (2 * r_norm ** 5),
                      [r[0] * (5 * (r[2] / r_norm) ** 2 - 1), r[1] * (5 * (r[2] / r_norm) ** 2 - 1),
                       r[2] * (5 * (r[2] / r_norm) ** 2 - 3)])
    a_3 = np.multiply(-5 * earth.J3 * earth.mu * earth.radius ** 3 / (2 * r_norm ** 7),
                      [r[0] * (3 * r[2] - 7 * r[2] ** 3 / r_norm ** 2), r[1] * (3 * r[2] - 7 * r[2] ** 3 / r_norm ** 2),
                       6 * r[2] ** 2 - 7 * r[2] ** 4 / r_norm ** 2 - (3 / 5) * r_norm ** 2])
    a_p = a_2 + a_3
    return np.array(a_p)

def SRP_Pert(r, r_sun, C_r, A_m):
    return np.array(np.multiply(earth.p_srp * C_r * A_m, np.divide(-1 * r_sun + r, lg.norm(r_sun - r))))

def drag_pert(r, v, density_table, C_D, A_m):
    if isinstance(r, np.ndarray) or isinstance(r, list):
        radius = np.linalg.norm(r)

    base_alt = np.array(density_table['Base Altitude'])
    scl_hgt = np.array(density_table['Scale Height'])
    nom_dens = np.multiply(np.array(density_table['Nominal Density']), 1000**3)
    del density_table
    idx = ((np.divide(np.abs(base_alt - radius), 10)).astype(int)).argmin()  # return the base altitude index
    density = nom_dens[idx]*np.exp(-(radius-base_alt[idx])/scl_hgt[idx])
    return np.array(np.multiply(-.5*C_D*A_m*density*lg.norm(v), v))
def sun_3body_pert(t, r):
    t = t / 86400 + 2451545
    sun_range, _ = PlanetRV(t)
    sun_range = np.matmul(J20002GCRF(), sun_range)
    sun_range = np.multiply(sun_range, -1)
    sat2sun = sun_range - r
    sat2sun_norm = lg.norm(sat2sun)
    sun_range_norm = lg.norm(sun_range)
    return sat2sun_norm, sat2sun, sun_range_norm, sun_range

def orbit_prop_all_pert(r_0, v_0, T0, tF, dT, conds):

    def three_body_orbit(t, Y, mu):
        dY = np.empty([6, 1])
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]
        r = lg.norm(Y[0:3])
        t = t/86400+conds.epoch
        sat2sun_norm, sat2sun, sun_range_norm, sun_range = sun_3body_pert(t, Y[0:3])
        a_d = drag_pert(Y[0:3], [Y[3:6]], conds.density_table, conds.C_D, conds.A_m)
        a_j = J2J3_Pert(Y[0:3])
        a_srp = SRP_Pert(Y[0:3], sun_range, conds.C_r, conds.A_m)
        a_other = np.squeeze(a_d+a_srp+a_j, axis=0)
        dY[3] = (-mu * Y[0] / r ** 3) + sun.mu*(sat2sun[0]/((sat2sun_norm)**3)-sun_range[0]/sun_range_norm**3) + a_other[0]
        dY[4] = (-mu * Y[1] / r ** 3) + sun.mu*(sat2sun[1]/((sat2sun_norm)**3)-sun_range[1]/sun_range_norm**3) + a_other[1]
        dY[5] = (-mu * Y[2] / r ** 3) + sun.mu*(sat2sun[2]/((sat2sun_norm)**3)-sun_range[2]/sun_range_norm**3) + a_other[2]
        return dY


    def derivFcn(t, y):
        return three_body_orbit(t, y, earth.mu)

    Y_0 = np.concatenate([r_0, v_0], axis=0)
    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y_0, T0)
    output = []
    output.append(np.insert(Y_0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:  # rv.successful() and
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    if not rv.successful() and rv.t<tF:
        warnings.warn("Runge Kutta Failed!", RuntimeWarning)
    #  Convert the output a numpy array for later use
    output = np.array(output)
    t = output[:, 0]

    r_vec = np.empty([np.shape(output)[0]-1, 3])
    v_vec = np.empty([np.shape(output)[0]-1, 3])

    for i in range(np.shape(output)[0]-1):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
    return r_vec, v_vec

def cylindrical_shadow(r_sc, r_sun):
    if np.dot(r_sc, r_sun/lg.norm(r_sun))< -np.sqrt(lg.norm(r_sc) ** 2 - earth.radius ** 2):
        return 0
    else:
        return 1

def period(a):
    return 2 * np.pi * np.sqrt(a ** 3 / earth.mu)

def up_shadow(r, sun_pos):
    a = np.arcsin(sun.radius / lg.norm(sun_pos + r))
    b = np.arcsin(earth.radius / lg.norm(r))
    c = np.arccos(np.dot(r, (sun_pos + r)) / (lg.norm(r) * lg.norm(sun_pos + r)))
    if c < np.abs(a - b):
        return 0
    elif a + b <= c:
        return 1
    else:
        x = (c ** 2 + a ** 2 - b ** 2) / (2 * c)
        y = np.sqrt(a ** 2 - x ** 2)
        A = a ** 2 * np.arccos(x / a) + b ** 2 * np.arccos((c - x) / b) - c * y
        return 1 - A / (np.pi * a ** 2)

def density(radius):
    if isinstance(radius, np.ndarray) or isinstance(radius, list):
        radius = np.linalg.norm(radius)
    url = 'https://raw.githubusercontent.com/ggb367/Spring-2020/master/366L/hw7/density.csv'
    altitude = radius - sf.earth.radius
    density = pd.read_csv(url)
    base_alt = np.array(density['Base Altitude'])
    scl_hgt = np.array(density['Scale Height'])
    nom_dens = np.multiply(np.array(density['Nominal Density']), 1000**3)
    del density
    idx = ((np.divide(np.abs(base_alt - altitude), 10)).astype(int)).argmin()  # return the base altitude index
    density = nom_dens[idx]*np.exp(-(altitude-base_alt[idx])/scl_hgt[idx])
    return density

def cw_prop(n_T, rho_rel_0, rho_dot_rel_0, T_0, T_F, dT):
    if not ((isinstance(rho_rel_0, np.ndarray) or isinstance(rho_rel_0, list)) and (isinstance(rho_dot_rel_0, np.ndarray) or isinstance(rho_dot_rel_0, list))):
        raise TypeError("Rho is of type ", type(rho_rel_0), "and Rho_dot is of type ", type(rho_dot_rel_0), "but they need to be a list or a numpy.ndarray!")

    x_0 = rho_rel_0[0]
    y_0 = rho_rel_0[1]
    z_0 = rho_rel_0[2]
    xd_0 = rho_dot_rel_0[0]
    yd_0 = rho_dot_rel_0[1]
    zd_0 = rho_dot_rel_0[2]

    def x(t):
        nt = n_T*t
        return (4-3*np.cos(nt))*x_0+(np.sin(nt)/n_T)*xd_0+(2/n_T)*(1-np.cos(nt))*yd_0
    def y(t):
        nt = n_T*t
        return 6*(np.sin(nt)-nt)*x_0+y_0+(2/n_T)*(np.cos(nt)-1)*xd_0+(1/n_T)*(4*np.sin(nt)-3*nt)*yd_0
    def z(t):
        nt = n_T*t
        return np.cos(nt)*z_0+(np.sin(nt)/n_T)*zd_0
    def xd(t):
        nt = n_T*t
        return 3*n_T*np.sin(nt)*x_0+np.cos(nt)*xd_0+2*np.sin(nt)*yd_0
    def yd(t):
        nt = n_T*t
        return 6*n_T*(np.cos(nt)-1)*x_0-2*np.sin(nt)*xd_0+(4*np.cos(nt)-3)*yd_0
    def zd(t):
        nt = n_T*t
        return -n_T*np.sin(nt)*z_0+np.cos(nt)*zd_0

    time_series = np.arange(T_0, T_F, dT)
    rho_rel = np.zeros([np.size(time_series), 3])
    rho_dot_rel = np.zeros([np.size(time_series), 3])

    count = 0
    # rho_rel[0, :] = rho_rel_0
    # rho_dot_rel[0, :] = rho_dot_rel_0
    for time_step in time_series:
        rho_rel[count, :] = np.array([x(time_step), y(time_step), z(time_step)])
        rho_dot_rel[count, :] = np.array([xd(time_step), yd(time_step), zd(time_step)])
        count = count+1

    return rho_rel, rho_dot_rel

def mean_motion(period, deg=False):
    if deg:
        return 360/period
    else:
        return 2*np.pi/period

def IOD_Gibbs(r1, r2, r3, t1, t2, t3):
    def Gibbs(r1, r2, r3):
        Z_12 = np.cross(r1, r2)
        Z_23 = np.cross(r2, r3)
        Z_31 = np.cross(r3, r1)
        r1_norm = lg.norm(r1)
        r2_norm = lg.norm(r2)
        r3_norm = lg.norm(r3)
        alpha_cop = np.pi/2-np.arccos(np.dot(Z_23/lg.norm(Z_23), np.divide(r1, r1_norm)))
        if alpha_cop >5*np.pi/180:
            print(Fore.RED+"Alpha_cop is: %f degrees" % (alpha_cop*180/np.pi))
            warnings.warn("Vectors may not be coplanar and IOD may breakdown")
            print(Fore.RESET)
        N = np.multiply(r1_norm, Z_23)+np.multiply(r2_norm, Z_31)+np.multiply(r3_norm, Z_12)
        D = Z_12+Z_23+Z_31
        S = np.multiply(r2_norm-r3_norm, r1)+np.multiply(r3_norm-r1_norm, r2)+np.multiply(r1_norm-r2_norm, r3)
        B = np.cross(D, r2)
        N_norm = lg.norm(N)
        D_norm = lg.norm(D)
        h = np.sqrt((N_norm*earth.mu)/D_norm)
        L_g = np.sqrt(earth.mu/(N_norm*D_norm))
        v2 = np.multiply(L_g/r2_norm, B) + np.multiply(L_g, S)
        return r2, v2

    def Herrick_Gibbs(r1, r2, r3, t1, t2, t3):
        r1_norm = lg.norm(r1)
        Z_23 = np.cross(r2, r3)
        alpha_cop = np.pi/2-np.arccos(np.dot(Z_23, np.divide(r1, r1_norm)))
        if alpha_cop>5*np.pi/180:
            warnings.warn("Vectors may not be coplanar and IOD may breakdown")
        delta_t31 = t3-t1
        delta_t32 = t3-t2
        delta_t21 = t2-t1
        r2_norm = lg.norm(r2)
        r3_norm = lg.norm(r3)
        x = -delta_t32*(1/(delta_t21*delta_t31)+earth.mu/(12*r1_norm**3))
        y = (delta_t32-delta_t21)*(1/(delta_t21*delta_t32)+earth.mu/(12*r2_norm**3))
        z = delta_t21*(1/(delta_t32*delta_t31)+earth.mu/(12*r3_norm**3))
        v2 = np.multiply(x, r1)+np.multiply(y, r2)+np.multiply(z, r3)
        return r2, v2

    case1 = np.dot(np.divide(r1, lg.norm(r1)), np.divide(r2, lg.norm(r2)))
    alpha_12 = np.arccos(case1)
    if r1[2]<0:
        alpha_12 = 2*np.pi-alpha_12
    case2 = np.dot(np.divide(r2, lg.norm(r2)), np.divide(r3, lg.norm(r3)))
    alpha_23 = np.arccos(case2)
    if r2[2]<0:
        alpha_23 = 2*np.pi-alpha_23
    if alpha_12<np.pi/180 or alpha_23<np.pi/180:
        print("Using Herrick-Gibbs because alpha_12 is %f, and alpha_23 is %f" % (alpha_12*180/np.pi, alpha_23*180/np.pi))
        return Herrick_Gibbs(r1, r2, r3, t1, t2, t3)
    elif alpha_12>5*np.pi/180 or alpha_23>5*np.pi/180:
        print("Using Gibbs because alpha_12 is %f, and alpha_23 is %f"%(alpha_12*180/np.pi, alpha_23*180/np.pi))
        return Gibbs(r1, r2, r3)
    else:
        print(Fore.RED+"We are in uncertain waters, alpha_12 is: %f and alpha_23 is %f"%(alpha_12, alpha_23))
        choosing = True
        while(choosing):
            choice = input("If you want to use Herrick-Gibbs, input 0, if you would like to use Gibbs, input 1")
            if choice:
                choosing = False
                return Gibbs(r1, r2, r3)
            elif not choice:
                choosing = False
                return Herrick_Gibbs(r1, r2, r3, t1, t2, t3)
            else:
                print("that was an invalid option")

def razel2SEZ(params , deg=True):
    range = params[0]
    az = params[1]
    el = params[2]
    if deg:
        az = np.deg2rad(az)
        el = np.deg2rad(el)
    r = np.array([-range*np.cos(el)*np.cos(az), range*np.cos(el)*np.sin(az), range*np.sin(el)])
    return r

def ITRF2SEZ(phi_gd, lamda):
    return np.matmul(R2(np.pi/2-phi_gd), R3(lamda))

def SEZ2ITRF(r_SEZ, R_ITRF):  # dumb, for use only in hw10
    return np.matmul(np.transpose(ITRF2SEZ(0, 0)), r_SEZ)+R_ITRF

def Lambert_Gauss(r_0, r, dt, t_m):
    r_0_norm = lg.norm(r_0)
    r_norm = lg.norm(r)
    cos_nu = np.dot(r_0, r)/(r_0_norm*r_norm)
    nu = np.arccos(np.dot(r_0, r)/(r_0_norm*r_norm))
    if r_0[2]<0:
        nu = 2*pi - nu

    t = (r_0_norm+r_norm)/(4*np.sqrt(r_0_norm*r_norm)*np.cos(nu/2))-1/2
    m = 1*dt**2/(2*np.sqrt(r_0_norm*r_norm)*np.cos(nu/2))**3
    y = 1
    error = 1
    count = 0
    while error > 10e-6:
        x1 = m/y**2 - t
        x2 = (4/3)*(1+6*x1/5+(6*8*x1**2)/(5*7)+(6*8*10*x1**3)/(5*7*9))
        this_y = 1+x2*(t+x1)
        error = np.abs(this_y-y)/y
        y = this_y
        count = count +1
        if count >10e5:
            warnings.warn("The method failed to converge")
            break
    cos_E = 1-2*x1
    p = (r_0_norm*r_norm*(1-cos_nu))/(r_0_norm+r_norm-2*np.sqrt(r_0_norm*r_norm)*np.cos(nu/2)*cos_E)
    f = 1-(r_norm/p)*(1-cos_nu)
    g = (r_norm*r_0_norm*np.sin(nu))/np.sqrt(1*p)
    # f_dot = np.sqrt(1/p)*np.tan(nu/2)*((1-np.cos(nu))/p - 1/r_norm-1/r_0_norm)
    g_dot = 1-(r_0_norm/p)*(1-np.cos(nu))
    v_0 = np.divide(r-np.multiply(f, r_0), g)
    v = np.divide(np.multiply(g_dot, r)-r_0, g)
    return v_0, v

def Lambert_MinEng(r_0, r):
    r_0_norm = lg.norm(r_0)
    r_norm = lg.norm(r)
    cos_nu = np.dot(r_0, r)/(r_0_norm*r_norm)
    sin_nu = np.sqrt(1-cos_nu**2)

    c = np.sqrt(r_0_norm**2+r_norm**2-2*r_0_norm*r_norm*cos_nu)
    s = (r_0_norm+r_norm+c)/2

    t_abs = (1/3)*np.sqrt(2/1)*(s**1.5-(s-c)**1.5)
    a_min = s/2
    p_min = (r_0_norm*r_norm/c)*(1-cos_nu)
    e_min = np.sqrt(1-(2*p_min)/s)

    v_0 = np.multiply(np.sqrt(1*p_min)/(r_0_norm*r_norm*sin_nu), (r-np.multiply((1-r_norm/p_min*(1-cos_nu)), r_0)))
    return a_min, e_min, t_abs ,v_0

def Lambert_Focus_Finder(r1 ,r2, a):
    d = np.abs(lg.norm(r1-r2))
    r1_norm = lg.norm(r1)
    r2_norm = lg.norm(r2)
    R = 2*a-r1_norm
    r = 2*a-r2_norm
    x = (d**2-r**2+R**2)/(2*d)
    a = (1-d)*np.sqrt(4*d**2*R**2-(d**2-r**2+R**2)**2)
    y = a/2
    # origin assumed at r1 but the origin is actually at the first focus
    f1 = np.array([x, y, 0])-r1
    f2 = np.array([x, -1*y, 0])-r1
    return f1, f2

