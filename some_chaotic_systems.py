import numpy as np
"""8 chaotic systems\n
include:\n
    lorenz\n
    chua\n
    chen\n
    lu\n
    rucklidge\n
    sprott_1\n
    sprott_2\n
    shimizu_morioka\n
"""


def lorenz(variables, t):
    """Lorenz, Edward N.. 
    "Deterministic Nonperiodic Flow." 
    Journal of the Atmospheric Sciences 20.2 (1963): 130-141.
    DOI: 10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2
    """
    a = 10
    b = 30
    c = 8/3
    x, y, z = variables
    dx = -a * (x - y)
    dy = b * x - x * z - y
    dz = -c * z + x * y
    return np.array([dx, dy, dz])


def chua(variables, t):
    """Chua, Leon O., M. Komuro, and T. Matsumoto. 
    "The double scroll family." 
    IEEE Transactions on Circuits and Systems 33.11 (1986): 1072-1118.
    DOI: 10.1109/TCS.1986.1085869
    """
    def f(x):
        m0 = -1/7
        m1 = 2/7
        return m1*x+0.5*(m0-m1)(abs(x+1)-abs(x-1))
    x, y, z = variables
    dx = 10 * (y - f(x))
    dy = x - y + z
    dz = -15 * y
    return np.array([dx, dy, dz])


def chen(variables, t):
    """Chen, Guanrong, and Tetsushi Ueta. 
    "Yet another chaotic attractor." 
    International Journal of Bifurcation and Chaos 9.7 (1999): 1465-1466.
    DOI: 10.1142/S0218127499001024
    """
    x, y, z = variables
    a = 35
    b = 3
    c = 28
    dx = -a*x+a*y
    dy = (c-a)*x+c*y-x*z
    dz = x*y-b*z
    return np.array([dx, dy, dz])


def lu(variables, t):
    """Lu, Jinhu, and Guanrong Chen. 
    "A NEW CHAOTIC ATTRACTOR COINED." 
    International Journal of Bifurcation and Chaos 12.3 (2002): 659-661.
    DOI: 10.1142/S0218127402004620
    """
    x, y, z = variables
    a = 36
    b = 3
    c = 20
    dx = -a*x+a*y
    dy = c*y-x*z
    dz = x*y-b*z
    return np.array([dx, dy, dz])


def rucklidge(variables, t):
    """Rucklidge, A. M.. 
    "Chaos in models of double convection." 
    Journal of Fluid Mechanics 237.-1 (1992): 209-229.
    DOI: 10.1017/S0022112092003392
    """
    x, y, z = variables
    a = 2
    b = 7.7
    dx = -a*x+b*y-y*z
    dy = x
    dz = y*y-z
    return np.array([dx, dy, dz])


def sprott_1(variables, t):
    """
    Sprott, J. C.. 
        "SOME SIMPLE CHAOTIC FLOWS." 
        Physical Review E 50.2 (1994).
        DOI:10.1103/PhysRevE.50.R647
    Sprott, J. C.. 
        "Simple chaotic systems and circuits." 
        American Journal of Physics 68.8 (2000): 758-763.
        DOI: 10.1119/1.19538
    Sprott, J. C.. 
        "A new class of chaotic circuit." 
        Physics Letters A 266.1 (2000): 19-23.
        DOI: 10.1016/S0375-9601(00)00026-8
    """
    x, y, z = variables
    dx = y*z
    dy = x-y
    dz = 1-x*x
    return np.array([dx, dy, dz])


def sprott_2(variables, t):
    """
    Sprott, J. C.. 
        "SOME SIMPLE CHAOTIC FLOWS." 
        Physical Review E 50.2 (1994).
        DOI:10.1103/PhysRevE.50.R647
    Sprott, J. C.. 
        "Simple chaotic systems and circuits." 
        American Journal of Physics 68.8 (2000): 758-763.
        DOI: 10.1119/1.19538
    Sprott, J. C.. 
        "A new class of chaotic circuit." 
        Physics Letters A 266.1 (2000): 19-23.
        DOI: 10.1016/S0375-9601(00)00026-8
    """
    x, y, z = variables
    dx = -x+y
    dy = -x+y*z
    dz = 1-y*y
    return np.array([dx, dy, dz])


def shimizu_morioka(variables, t):
    """Shimizu, T., and N. Morioka. 
    "On the bifurcation of a symmetric limit cycle to an asymmetric one in a simple model." 
    Physics Letters A (1980): 201-204.
    DOI: 10.1016/0375-9601(80)90466-1
    """
    x, y, z = variables
    a = 0.75
    b = 0.45
    dx = y
    dy = (1-z)*x-a*y
    dz = x*x-b*z
    return np.array([dx, dy, dz])
