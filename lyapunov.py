import numpy
import math
def lyapunov(n, rhs_ext_fcn, fcn_integrator, tstart, stept, tend, ystart, ioutp):
    """
    fcn_integrator: callable
    """
    n1 = n
    n2 = n1*(n1+1)
    # number of steps
    nit = round((tend - tstart) / stept)
    # memory allocation
    y=numpy.zeros(n2)
    cum=numpy.zeros(n1)
    y0=y
    gsc=cum
    znorm=cum
    # initial values
    y[:] = ystart[:]
    for i in range(n1):
        y[(n1+1)*i]=1
    t=tstart
    # main loop
    for ITERLYAP in range(nit):
        # solution of extended ode system
        [T, Y] = fcn_integrator(rhs_ext_fcn, [t, t+stept], y)
        t += stept
        y = Y[size(Y,1)] # todo
        #
        for i in range(n1):
            for j in range(n1):
                y0[n1*i+j] = y[n1*j+1]
        #
        znorm[0] = 0 # todo
        for j in range(n1):
            znorm[0] = znorm[0] + y0[n1 * j + 1]**2
        
        
        
        #
        znorm[0] = math.sqrt(znorm[0])
        for j in range(n1):
            y0[n1*j] /= znorm[0] 
        #
        for i in range(1, n1):
            

