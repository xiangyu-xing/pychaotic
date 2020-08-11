import numpy as np
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def dmove(Point, t, sets):
    """
    point：present location index
    sets：super parameters
    """
    p, r, b = sets
    x, y, z = Point
    return np.array([p * (y - x), x * (r - z), x * y - b * z])


t = np.arange(0, 30, 0.001)
P1 = odeint(dmove, (0., 1., 0.), t, args=([10., 28., 3.],))  #
## (0.,1.,0.) is the initial point; args is the set for super parameters
P2 = odeint(dmove, (0., 1.01, 0.), t, args=([10., 28., 3.],))
## slightly change the initial point from y = 1.0 to be y = 1.01

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(P1[:, 0], P1[:, 1], P1[:, 2])
ax.plot(P2[:, 0], P2[:, 1], P2[:, 2])
plt.show()
