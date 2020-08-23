from typing import NoReturn
import numpy
import math
from scipy.integrate import odeint
from scipy.linalg import solve_lyapunov


def lyapunov(n, rhs_ext_fcn, fcn_integrator, tstart, stept, tend, ystart):
    """
    fcn_integrator: callable
    """
    n1 = n
    n2 = n1*(n1+1)
    # number of steps
    nit = round((tend - tstart) / stept)
    # memory allocation
    y = numpy.zeros(n2)
    cum = numpy.zeros(n1)
    y0 = numpy.zeros(n2)
    gsc = numpy.zeros(n1)
    znorm = numpy.zeros(n1)
    # result
    Lexp = []
    Texp = []
    # initial values
    y[0:n] = ystart[:]
    for i in range(n1):
        y[n1*(i+1)+i] = 1
    t = tstart
    for ITERLYAP in range(nit):  # main loop
        # solution of extended ode system
        tt = t+stept
        Y = fcn_integrator(rhs_ext_fcn, y, [t, t+stept])
        t += stept
        y = Y[-1]
        #
        for i in range(1, n1+1):
            for j in range(1, n1+1):
                y0[n1*i+j-1] = y[n1*j+i-1]
        #
        znorm[0] = 0
        for j in range(n1):
            # print(n1*(j+1))
            znorm[0] = znorm[0] + y0[n1 * (j+1)]**2
        #
        znorm[0] = math.sqrt(znorm[0])
        for j in range(n1):
            y0[n1*(j+1)] /= znorm[0]
        #
        for j in range(1, n1):
            #
            for k in range(j):
                gsc[k] = 0
                for jj in range(n1):
                    gsc[k] = gsc[k] + y0[n1*(jj+1)+j]*y0[n1*(jj+1)+k]

            #
            for k in range(n1):
                for jj in range(j):
                    y0[n1*(k+1)+j] = y0[n1*(k+1)+j]-gsc[jj]*y0[n1*(k+1)+jj]

            #
            znorm[j] = 0
            for k in range(n1):
                znorm[j] += y0[n1*(k+1)+j]**2

            znorm[j] = math.sqrt(znorm[j])
            for k in range(n1):
                y0[n1*(k+1)+j] = y0[n1*(k+1)+j]/znorm[j]
            pass
        # updata runing vector magnitudes
        for k in range(n1):
            cum[k] += math.log(znorm[k])
        # normalize exponent
        lp = [0] * n1
        for k in range(n1):
            lp[k] = cum[k]/(t-tstart)
        Lexp.append(lp)
        Texp.append(t)
        for i in range(1, n1+1):
            for j in range(1, n1+1):
                y[n1*j+i-1] = y0[n1*i+j-1]
    return Lexp, Texp


def lorenz_jaco(variables, t):
    """Lorenz, Edward N.. 
    "Deterministic Nonperiodic Flow." 
    Journal of the Atmospheric Sciences 20.2 (1963): 130-141.
    DOI: 10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2
    """
    a = 10
    b = 30
    c = 8/3
    x, y, z = variables[0], variables[1], variables[2]
    Y = [[variables[3], variables[6], variables[9]],
         [variables[4], variables[7], variables[10]],
         [variables[5], variables[8], variables[11]]]
    Y = numpy.array(Y)

    dx = numpy.zeros(12)
    dx[0] = -a * (x - y)
    dx[1] = b * x - x * z - y
    dx[2] = -c * z + x * y

    jaco = [[-a, a, 0], [b-z, -1, -x], [y, x, -c]]
    jaco = numpy.array(jaco)
    temp = numpy.matmul(jaco, Y)
    dx[3:12] = temp.reshape((9,), order="F")
    return dx


if __name__ == "__main__":
    res,_ = lyapunov(3, lorenz_jaco, odeint, 0, 0.01, 100, [1, 1, 1])
