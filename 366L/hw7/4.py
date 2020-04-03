import numpy as np

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


