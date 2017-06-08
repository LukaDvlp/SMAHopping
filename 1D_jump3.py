import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Mb = 0.7  #Mass of Body [Kg]
Mp = 0.01 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
k1 = 1.4   #Bias Spring Coefficient [N/mm]
l0 = 100  #Natural Length of Bias Spring [mm]
k2 = 100  #Spring Coefficient of the ground [N/mm]
c2 = 100  #Damping Coefficient [Ns/mm]
Z20 = 0.0  #Initial Position of Pad [mm]
DZ = 50    #Initial Deflextion of Bias Spring [mm]

def func1(x):
    return np.array([x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-(c2/Mp)*x[3]-g])

def func2(x):
    return np.array([x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-g])
	
def func3(x):
    return np.array([x[1], (k1/Mb)*(l0-(x[0]-x[2]))-g, x[3], (k1/Mp)*(l0-(x[0]-x[2])\
	)-(k2/Mp)*(x[2]-Z20)-g])

def test_func(x):
	return np.array([x[3],x[2],x[1],x[0]])

def motion_test(x):
	return np.array([x[2],x[1],x[0]])

def motion_test(x):
	return np.array([x[1],-g])

def main():
	h = 0.1
	t_s = 0.0
	t_f = 1.0 
	t = t_s 
	n = 0
	m = 1
	v0 = 10
	H = 1
	#Time = np.empty((0), float)
	#Time = np.append(Time, np.array([[xH,v0]]), axis=0)
	Time = []
	XX = np.empty((0,4), float)
	XX = np.append(XX, np.array([[l0-DZ,0,0,0]]), axis=0)
#	XX = np.append(XX, np.array([[4,5,6]]), axis=0)
	#print(XX)
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
	while (t < t_f):
		k1 = test_func(XX[n])
	#	k1 = np.array([3.0,2.0,1.0])
#		print(k1)
	#	bb = 3.2*k1
		k2 = test_func(XX[n]+0.5*h*k1)
		k3 = test_func(XX[n]+0.5*h*k2)
		k4 = test_func(XX[n]+h*k3)
		Step = XX[n] + (h/6)*(k1+2*k2+2*k3+k4)
#		print(Step)
		S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
#		A = np.array([[n,n,n]])
		XX = np.append(XX, S, axis=0)
		#X(n+1,:) = XX(n,:) + B 
#		Time = Time + [t]
		t = t+h
		n = n+1

	print(n)
	print(t)
	print(XX)
	print(Time)
	T = np.arange(0, t_f+2*h, h)
	print(T)

	#x = odeint(func1, x0, t)
	#X2 = motion_test(x0)
	#Z_cg = Mb/(Mb+Mp)*x[:,0]
	#V_cg = Mb/(Mb+Mp)*x[:,1] 

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
	plt.plot(T,XX[:,1], label="Velocity")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
