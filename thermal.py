import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

d = 1    #Diameter of SMA wire [mm]
D = 6.8  #Diameter of SMA coil [mm]  
n = 10   #Number of coil spring 
rho = 8.5
r = math.sqrt(math.pi*d*D*n/4) #radius of pseudo sphere of SMA[mm]
A = math.pi*r**2 #Cross section of the pseudo sphere of SMA[mm^2]
AA = 4*math.pi*r**2 #Surface are of SMA pseudo sphere [mm^2]
Theta_L = 45  #Lattitude [deg]
m = rho*(math.pi*(0.5*d))**2*D*n/1000 #Mass of a SMA spring [g]
c = 440 #Specific heat capacity of SMA spring
Sc = 1366*10**-6 #Solar constant [W/mm^2]
delta = 5.67*10**-2 #Stefan-Boltzman constant [W/mm^2K^4]
a = 0.07 #Albedo constant
Tg = 100 + 273 #Temperature of the ground[K]
epsilon = 0.6 #Emissivity of the ground

def func1(x, t):
	return [A*math.cos(math.radians(Theta_L))*Sc/(m*c) + a*A*math.cos(math.radians(Theta_L))*Sc/(m*c) +epsilon*delta*A/(m*c)*(Tg**4-x[0]**4)/2 - delta*AA/(m*c)*(x[0]**4-4**4)] 

def func2(x, t):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-g]

def main():
	x0 = [80+273]
	t = np.arange(0, 10, 0.1)

	x = odeint(func1, x0, t)

	print(x)
	#fig = plt.figure()
	plt.plot(t,x[:,0], label="Temperature of SMA")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
