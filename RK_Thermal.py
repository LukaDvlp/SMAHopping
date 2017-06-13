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
#rho = 6.5*10-3 #Density of SMA wire [g/mm^3]
rho = 1
r = math.sqrt(math.pi*d*D*n/4) #radius of pseudo sphere of SMA[mm]
#A = math.pi*r**2 #Cross section of the pseudo sphere of SMA[mm^2]
A = 1
#AA = 4*math.pi*r**2 #Surface are of SMA pseudo sphere [mm^2]
AA = 1
Theta_L = 45  #Lattitude [deg]
#m = rho*(math.pi*(0.5*d))**2*D*n/1000 #Mass of a SMA spring [g]
m = 1
#c = 440*10**-3  #Specific heat capacity of SMA spring [J/gK]
c = 500 
#Sc = 1366*10**-6 #Solar constant [W/mm^2]
Sc = 1
#delta = 5.67*10**-2 #Stefan-Boltzman constant [W/mm^2K^4]
delta = 1
#a = 0.07 #Albedo constant
a = 1
#Tg = 50 + 273 #Temperature of the ground[K]
Tg = 1
epsilon = 1.0 #Emissivity of the ground

def func1(x):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-(c2/Mp)*x[3]-g]

def func2(x, t):
    return [x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-g]
	
#def motion_test(x):
	#return np.array([x[2],x[1],x[0]])

def motion_test(x):
	return np.array([x[1],-g])

def ThermalEq(x):
	value1 = A*math.cos(math.radians(Theta_L))*Sc/(m*c)
	print("value1")
	print(value1)
	value2 = a*A*math.cos(math.radians(Theta_L))*Sc/(m*c)
	print("value2")
	print(value2)
	value3 = epsilon*delta*A*(Tg**4-x[0]**4)/(2*m*c)
	print("value3")
	print(value3)
	value4 = delta*AA*(x[0]**4-4**4)/(m*c)
	print("value4")
	print(value4)
	q = value1 + value2 + value3 - value4 
	return q
'''	
def RK(x):
	k1 = motion_test(x)
	k2 = motion_test(x+0.5*h*k1)
	k3 = motion_test(x+0.5*h*k2)
	k4 = motion_test(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	print(x_)
	return x_
'''
def RK(x):
	k1 = ThermalEq(x)
	k2 = ThermalEq(x+0.5*h*k1)
	k3 = ThermalEq(x+0.5*h*k2)
	k4 = ThermalEq(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
#'''
def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			Step = RK(XX[n])
			print(n)
			S = np.array([[Step[0]]])
			#S = np.array(Step)
			XX = np.append(XX, S, axis=0)
			t = t+h
			n = n+1
	
	return (XX)

def main():
	t_s = 0.0
	t_f = 10.0 
	t = 0
	n = 0
	m = 1
	v0 = 10
	H = 1
	#X0 = [H,v0]
	X0 = [1]
	#print(X0)
	#print("This is for test")
#	print(RK(X0))
#	print(RK(X0)[1])
	XX = Cal_Mtlx(X0, t_s, t_f, 1)
	#Time = np.empty((0), float)
	#Time = np.append(Time, np.array([[xH,v0]]), axis=0)
	Time = []
#	XX = np.empty((0,2), float)
#	XX = np.append(XX, np.array([[H,v0]]), axis=0)
#	XX = np.append(XX, np.array([[4,5,6]]), axis=0)
	print(XX)
	#print(XX[1])
	#XX = ([0,0,0,0])
	#XX.append([1,1,1,1])
	#print(XX)
	#print(XX(2,1))
	#x0 = [l0-DZ, 0, 0, 0]
	#XX = np.zeros((t_s/h, 4))
	#print(XX)
	#XX(0,:) = [0,0,0,0] 
	#B = np.array([1,1,1,1])
#	while (t < t_f):
	#	k1 = motion_test(XX[n])
	#	k1 = np.array([3.0,2.0,1.0])
#		print(k1)
	#	bb = 3.2*k1
	#	k2 = motion_test(XX[n]+0.5*h*k1)
	#	k3 = motion_test(XX[n]+0.5*h*k2)
	#	k4 = motion_test(XX[n]+h*k3)
	#	Step = XX[n] + (h/6)*(k1+2*k2+2*k3+k4)
#		print(Step)
#		Step = RK(XX[n])
#		S = np.array([[Step[0],Step[1]]])
#		A = np.array([[n,n,n]])
#		XX = np.append(XX, S, axis=0)
		#X(n+1,:) = XX(n,:) + B 
#		Time = Time + [t]
#		t = t+h
#		n = n+1

	#print(n)
	#print(t)
#	print(XX)
#	print(Time)
	T = np.arange(0, t_f+2*h, h)
	print(T)

	#x = odeint(func1, x0, t)
	#X2 = motion_test(x0)
	#Z_cg = Mb/(Mb+Mp)*x[:,0]
	#V_cg = Mb/(Mb+Mp)*x[:,1] 

#	print("This is for test")
#	print(XX[n])
#	print("Call RK function")
#	RK(XX[n])
	#print(X2)
	#print(np.max(V_cg))
	#print(np.argmax(V_cg))
	#print(x[10,:])
	#index = np.where((x[:,3]>=-0.01)&(x[:,3]<=-0.001))
	#print(index)
	#fig = plt.figure()
	#plt.plot(t,x[:,0], label="Position:Body")
	#plt.plot(t,x[:,1], label="Velocity:Body")
	#plt.plot(t,x[:,2], label="Position:Pad")
	#plt.plot(t,x[:,3], label="Velocity:Pad")
	#plt.plot(t,Z_cg, label="Position:COG")
	#plt.plot(t,V_cg, label="Velocity:COG")
	plt.plot(T,XX[:,0], label="Positio")
#	plt.plot(T,XX[:,1], label="Velocity")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
