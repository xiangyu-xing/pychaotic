import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


g = 9.8
l = 1
def diff2(d_list, t):
	omega, theta = d_list
	return np.array([-g/l*theta, omega])
t = np.linspace(0, 20, 2000)
result = odeint(diff2, [0, 35/180*np.pi], t)
# 结果是一个两列的矩阵， odeint中第二个是初始单摆角度35度
plt.plot(t, result[:, 0])  # 输出omega随时变化曲线
plt.plot(t, result[:, 1])  # 输出theta随时变化曲线，即方程解
plt.show()


t= [[1,2,3],[4,5,6],[7,8,9]]
print(np.array(t).reshape((1,9),order = "F"))
