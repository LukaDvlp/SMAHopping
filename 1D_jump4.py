import numpy as np
import math
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
h = 0.1   #Interval of RK
d = 1    #Diameter of SMA wire [mm]
D = 6.8  #Diameter of SMA coil [mm]
n = 10   #Number of coil spring
rho = 6.5*10-3 #Density of SMA wire [g/mm^3]
r = math.sqrt(math.pi*d*D*n/4) #radius of pseudo sphere of SMA[mm]
A = math.pi*r**2 #Cross section of the pseudo sphere of SMA[mm^2]
AA = 4*math.pi*r**2 #Surface are of SMA pseudo sphere [mm^2]
Theta_L = 45  #Lattitude [deg]
m = rho*(math.pi*(0.5*d))**2*D*n/1000 #Mass of a SMA spring [g]
c = 440*10**-3  #Specific heat capacity of SMA spring [J/gK]
Sc = 1366*10**-6 #Solar constant [W/mm^2]
delta = 5.67*10**-2 #Stefan-Boltzman constant [W/mm^2K^4]
a = 0.07 #Albedo constant
Tg = 50 + 273 #Temperature of the ground[K]
epsilon = 1.0 #Emissivity of the ground

def func1(x):
    return np.array([x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-(c2/Mp)*x[3]-g])

def func2(x, t):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-g]
	
def func2(x, t):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-g]

#def motion_test(x):
	#return np.array([x[2],x[1],x[0]])

def motion_test(x):
	return np.array([x[1],-g])

def ThermalEq(x):
	q = A*math.cos(math.radians(Theta_L))*Sc/(m*c) + a*A*math.cos(math.radians(Theta_L))*Sc/(m*c) +epsilon*delta*A/(m*c)*(Tg**4-x[0]**4)/2 - delta*AA/(m*c)*(x[0]**4-4**4)
	return q

def RK(x):
	k1 = func1(x)
	k2 = func1(x+0.5*h*k1)
	k3 = func1(x+0.5*h*k2)
	k4 = func1(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			Step = RK(XX[n])
			S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
			#S = np.array(Step)
			XX = np.append(XX, S, axis=0)
			t = t+h
			n = n+1
	
	return np.matrix(XX)

def main():
	t_s = 0.0
	t_f = 10.0 
	t = 0
	n = 0
	m = 1
	v0 = 10
	H = 1
	X0 = [l0-DZ,0,0,0]
	Time = []
	print(X0)
	print("This is for test")
#	print(RK(X0))
#	print(RK(X0)[1])
	XX = Cal_Mtlx(X0, t_s, t_f, 4)
	print(XX)

	print(n)
	print(t)
	print(Time)
	T = np.arange(0, t_f+2*h, h)
	print(T)


	plt.plot(T,XX[:,0], label="Positio")
	plt.plot(T,XX[:,1], label="Velocity")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
