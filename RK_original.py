import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m = 1
g = 1.622 
h = 10 
v0 = 10

def func(x, t):
    return [x[1],-g]

def main():
	x0 = [h, v0]
	t = np.arange(0, 10, 0.1)

	x = odeint(func, x0, t)


	#print(x)
	#fig = plt.figure()
	plt.plot(t,x[:,0], label="Position")
	plt.plot(t,x[:,1], label="Velocity")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
