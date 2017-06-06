import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Mb = 0.7  #Mass of Body [Kg]
Mp = 0.01 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
k1 = 1.4   #Bias Spring Coefficient [N/mm]
l0 = 100  #Natural Length of Bias Spring [mm]
k2 = 10  #Spring Coefficient of the ground [N/mm]
c2 = 10  #Damping Coefficient [Ns/mm]
Z20 = 0.0  #Initial Position of Pad [mm]
DZ = 50    #Initial Deflextion of Bias Spring [mm]

def func1(x, t):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-(c2/Mp)*x[3]-g]

def func2(x, t):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-g]

def main():
	x0 = [l0-DZ, 0, 0, 0]
	t = np.arange(0, 10, 0.01)

	x = odeint(func1, x0, t)
	Z_cg = Mb/(Mb+Mp)*x[:,0]
	V_cg = Mb/(Mb+Mp)*x[:,1] 

	print(x)
	print(np.max(V_cg))
	print(np.argmax(V_cg))
	print(x[10,:])
	index = np.where((x[:,3]>=-0.01)&(x[:,3]<=-0.001))
	print(index)
	#fig = plt.figure()
	plt.plot(t,x[:,0], label="Position:Body")
	plt.plot(t,x[:,1], label="Velocity:Body")
	plt.plot(t,x[:,2], label="Position:Pad")
	plt.plot(t,x[:,3], label="Velocity:Pad")
	plt.plot(t,Z_cg, label="Position:COG")
	plt.plot(t,V_cg, label="Velocity:COG")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
