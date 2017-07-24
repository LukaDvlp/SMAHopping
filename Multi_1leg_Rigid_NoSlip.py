import numpy as np
import math
import matplotlib.pyplot as plt

m1 = 0.65  #Mass of Body [Kg]
m2 = 0.05 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
k1 = 1400   #Bias Spring Coefficient [N/m]
c1 = 1  #Damping Coefficient of Spring [Ns/mm] 
l0 = 0.1  #Natural Length of Bias Spring [m]
k2 = 1000  #Spring Coefficient of the ground [N/mm]
c2 = 100 #Damping Coefficient [Ns/mm]
Z20 = 0.0  #Initial Position of Pad [mm]
DZ = 0.05    #Initial Deflextion of Bias Spring [mm]
h = 0.0001   #Interval of RK
d = 1    #Diameter of SMA wire [mm]
D = 6.8  #Diameter of SMA coil [mm]
n = 10   #Number of coil spring
del_AE = 0.04301 

def func1(x, l):
	#l = np.sqrt((x[0]-x[2])**2 + (x[1]-x[3])**2)
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l)
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	term3 = 0
	term4 = 0 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])

def func2(x, l):
	#l = np.sqrt((x[0]-x[2])**2 + (x[1]-x[3])**2)
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l)
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	term3 = -f*np.cos(th) / m2 
	term4 = -f*np.sin(th) / m2 -g 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])


def func_test(x):
	term1 = -k1/Mb*x[0] - c1/Mb*x[1] + k1/Mb*l0
	print"term1 = {0}".format(term1)
    #return np.array([x[1],term1])
	return np.array([x[1],term1])

def test_func2(x):
	term1 = k1/M1*(l0-(x[0]-x[2])) - c1/M1*(x[1]-x[3]) - g
	print"term1 = {0}".format(term1)
	term2 = -k1/M2*(l0-(x[0]-x[2])) + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - c2/M2*x[3] - g
	print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 
	
def test_func3(x):
	term1 = k1/M1*(l0-(x[0]-x[2])) - c1/M1*(x[1]-x[3]) - g
	print"term1 = {0}".format(term1)
	term2 = -k1/M2*(l0-(x[0]-x[2])) + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - g
	print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def test_func4(x):
	term1 = k1/M1*(l0-(x[0]-x[2])) - c1/M1*(x[1]-x[3]) -g
	print"term1 = {0}".format(term1)
	term2 = -k1/M2*(l0-(x[0]-x[2])) + c1/M2*(x[1]-x[3]) -g
	print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def motion_test(x):
	return np.array([x[1],-g])

def RK(x, f, l):
	k1 = f(x, l)
	k2 = f(x+0.5*h*k1, l)
	k3 = f(x+0.5*h*k2, l)
	k4 = f(x+h*k3, l)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3],\
	                              X0[4],X0[5],X0[6],X0[7]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			l = np.sqrt((XX[n,0]-XX[n,4])**2 + (XX[n,2]-XX[n,6])**2)
			#if XX[n,2]<0 and XX[n,3]<0:
			if l0 - l > 0:
				Step = RK(XX[n], func1, l)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func1")
				#print"term1 = {0}".format(term1)
				#print"term2 = {0}".format(term2)
				t = t+h
				n = n+1
			else:	
				break

	while(t<t_f):
			l = np.sqrt((XX[n,0]-XX[n,4])**2 + (XX[n,2]-XX[n,6])**2)
			if XX[n,2] >= 0 and XX[n,6] >= 0:
				Step = RK(XX[n], func2, l)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
			               Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func2")
				t = t+h
				n = n+1
			else:
				break

	return np.matrix(XX)

def main():
	t_s = 0.0
	t_f = 3.0 
	th0 = np.pi/4
	X0 = [del_AE*np.cos(th0),0,del_AE*np.sin(th0),0,0,0,0,0]
	print(X0)
#	print(RK(X0))
#	print(RK(X0)[1])
	XX = Cal_Mtlx(X0, t_s, t_f, 8)
	print(XX)

	#print(Time)
	#T = np.arange(0, t_f+2*h, h)
#	print(T)
	print("Size of XX")
	row, col = XX.shape
	print(row, col)
	T = np.arange(0, row*h, h)
	print("Length of T={0}".format(len(T)))

	plt.plot(T,XX[:,0], label="x1")
	plt.plot(T,XX[:,2], label="z1")
	plt.plot(T,XX[:,4], label="x2")
	plt.plot(T,XX[:,6], label="z2")
	#plt.plot(T,XX[:,1], label="Velocity of X1")
	#plt.plot(T,XX[:,2], label="Position of X2")
	#plt.plot(T,XX[:,3], label="Velocity of X2")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.title('1-Dimensional Hopping of Spring Only Rover')
	plt.xlabel('Time[s]')
	plt.ylabel('Hopping Height[m]')
	plt.xlim([0,T[row-1]])
	#plt.xlim([0,2.60])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()