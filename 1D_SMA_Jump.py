import numpy as np
import math
import matplotlib.pyplot as plt
import Simu0619
import NoSMAJump 

M1 = 0.65  #Mass of Body [Kg]
M2 = 0.05 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
k1 = 1400   #Bias Spring Coefficient [N/m]
c1 = 1  #Damping Coefficient of Spring [Ns/m] 
l0 = 0.1  #Natural Length of Bias Spring [m]
k2 = 1000  #Spring Coefficient of the ground [N/m]
c2 = 100 #Damping Coefficient [Ns/m]
Z20 = 0.0  #Initial Position of Pad [m]
DZ = 0.05    #Initial Deflextion of Bias Spring [m]
h = 0.00002   #Interval of RK
A = 1/Simu0619.A
print(Simu0619.n)
del_AE = Simu0619.delta_AE/1000
del_s = Simu0619.delta_s/1000
F_s = Simu0619.F_s

def test_func2(x):
	if (x[0]-x[2])-del_AE < del_s:
		term1 = (k1*(l0-(x[0]-x[2]))-A*((x[0]-x[2])-del_AE))/M1 - c1/M1*(x[1]-x[3]) - g
		print"term1 = {0}".format(term1)
		print"SMA resistance = {0}".format(A*(x[0]-x[2]-del_AE))
		term2 = -(k1*(l0-(x[0]-x[2]))-A*((x[0]-x[2])-del_AE))/M2 + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - c2/M2*x[3] - g
		print"term2 = {0}".format(term2)
	else:
		term1 = (k1*(l0-(x[0]-x[2]))-F_s)/M1 - c1/M1*(x[1]-x[3]) - g
		print"term1 = {0}".format(term1)
		print"SMA resistance = {0}".format(F_s)
		term2 = -(k1*(l0-(x[0]-x[2]))-F_s)/M2 + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - c2/M2*x[3] - g
		print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 
	
def test_func3(x):
	if (x[0]-x[2])-del_AE < del_s:
		term1 = (k1*(l0-(x[0]-x[2]))-A*((x[0]-x[2])-del_AE))/M1 - c1/M1*(x[1]-x[3]) - g
		print"term1 = {0}".format(term1)
		print"SMA resistance = {0}".format(A*(x[0]-x[2]-del_AE))
		term2 = -(k1*(l0-(x[0]-x[2]))-A*((x[0]-x[2])-del_AE))/M2 + c1/M2*(x[1]-x[3]) - k2/M2*x[2] -  g
		print"term2 = {0}".format(term2)
	else:
		term1 = (k1*(l0-(x[0]-x[2]))-F_s)/M1 - c1/M1*(x[1]-x[3]) - g
		print"term1 = {0}".format(term1)
		print"SMA resistance = {0}".format(F_s)
		term2 = -(k1*(l0-(x[0]-x[2]))-F_s)/M2 + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - g
		print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def test_func4(x):
	if (x[0]-x[2])-del_AE < del_s:
		term1 = (k1*(l0-(x[0]-x[2]))-A*((x[0]-x[2])-del_AE))/M1 - c1/M1*(x[1]-x[3]) - g
		print"term1 = {0}".format(term1)
		print"SMA resistance = {0}".format(A*(x[0]-x[2]-del_AE))
		term2 = -(k1*(l0-(x[0]-x[2]))-A*((x[0]-x[2])-del_AE))/M2 + c1/M2*(x[1]-x[3])- g
		print"term2 = {0}".format(term2)
	else:
		term1 = (k1*(l0-(x[0]-x[2]))-F_s)/M1 - c1/M1*(x[1]-x[3]) - g
		print"term1 = {0}".format(term1)

		print"SMA resistance = {0}".format(F_s)
		term2 = -(k1*(l0-(x[0]-x[2]))-F_s)/M2 + c1/M2*(x[1]-x[3]) - g
		print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def NoSMA_func2(x):
	term1 = k1/M1*(l0-(x[0]-x[2])) - c1/M1*(x[1]-x[3]) - g
	print"term1 = {0}".format(term1)
	term2 = -k1/M2*(l0-(x[0]-x[2])) + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - c2/M2*x[3] - g
	print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def NoSMA_func3(x):
	term1 = k1/M1*(l0-(x[0]-x[2])) - c1/M1*(x[1]-x[3]) - g
	print"term1 = {0}".format(term1)
	term2 = -k1/M2*(l0-(x[0]-x[2])) + c1/M2*(x[1]-x[3]) - k2/M2*x[2] - g
	print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def NoSMA_func4(x):
	term1 = k1/M1*(l0-(x[0]-x[2])) - c1/M1*(x[1]-x[3]) -g
	print"term1 = {0}".format(term1)
	term2 = -k1/M2*(l0-(x[0]-x[2])) + c1/M2*(x[1]-x[3]) -g
	print"term2 = {0}".format(term2)
	return np.array([x[1], term1, x[3], term2]) 

def motion_test(x):
	return np.array([x[1],-g])

def RK(x, f):
	k1 = f(x)
	k2 = f(x+0.5*h*k1)
	k3 = f(x+0.5*h*k2)
	k4 = f(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3]]]), axis=0)
	#XX = np.append(XX, np.array([[X0[0],X0[1]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			if XX[n,2]<0 and XX[n,3]<0:
			#if 1>0:
				Step = RK(XX[n], test_func2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				#S = np.array([[Step[0],Step[1]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func1")
				#print"term1 = {0}".format(term1)
				#print"term2 = {0}".format(term2)
				t = t+h
				n = n+1
			elif XX[n,2]<0 and XX[n,3]>=0:
				Step = RK(XX[n], test_func3)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				#S = np.array([[Step[0],Step[1]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func2")
				#print"term1 = {0}".format(term1)
				#print"term2 = {0}".format(term2)
				t = t+h
				n = n+1
			else:	
				Step = RK(XX[n], test_func4)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				#S = np.array([[Step[0],Step[1]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func3")
				t = t+h
				n = n+1

	return np.matrix(XX)
	
def Cal_Mtl_NoSMA(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3]]]), axis=0)
	#XX = np.append(XX, np.array([[X0[0],X0[1]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			if XX[n,2]<0 and XX[n,3]<0:
			#if 1>0:
				Step = RK(XX[n], NoSMA_func2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				#S = np.array([[Step[0],Step[1]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func1")
				#print"term1 = {0}".format(term1)
				#print"term2 = {0}".format(term2)
				t = t+h
				n = n+1
			elif XX[n,2]<0 and XX[n,3]>=0:
				Step = RK(XX[n], NoSMA_func3)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				#S = np.array([[Step[0],Step[1]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func2")
				#print"term1 = {0}".format(term1)
				#print"term2 = {0}".format(term2)
				t = t+h
				n = n+1
			else:	
				Step = RK(XX[n], NoSMA_func4)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				#S = np.array([[Step[0],Step[1]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func3")
				t = t+h
				n = n+1

	return np.matrix(XX)

def main():
	t_s = 0.0
	t_f = 3.0 
	v0 = 10
	H = 1
	X0 = [del_AE,0,0.0,0]
	#X0 = [l0-0.05,0]
	#Time = []
	print(X0)
#	print(RK(X0))
#	print(RK(X0)[1])
	XX = Cal_Mtlx(X0, t_s, t_f, 4)
	YY = Cal_Mtl_NoSMA(X0, t_s, t_f, 4)
	print(XX)
	#print(Time)
	#T = np.arange(0, t_f+2*h, h)
	T = np.arange(0, t_f+h, h)
#	print(T)
	print("Length of T={0}".format(len(T)))
	print("Size of XX")
	row, col = XX.shape
	print(row, col)

	plt.plot(T,XX[:,0], label="SMA Rover")
	plt.plot(T,YY[:,0], label="Spring Only", color='black', linestyle=":")
	#plt.plot(T,XX[:,2], label="Height of Rover's Pad")
	#plt.plot(T,XX[:,2], label="Position of X2")
	#plt.plot(T,XX[:,3], label="Velocity of X2")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.title('1-Dimensional Hopping of SMA Rover')
	plt.xlabel('Time[s]')
	plt.ylabel('Hopping Height[m]')
	plt.xlim([0,2.16])
	plt.legend()
	plt.show()
	#plt.plot(T,XX[:,1], label="Height of Rover's Body")
	#plt.plot(T,XX[:,3], label="Velocity of Rover's Pad")
	#plt.title('1-Dimensional Hopping of SMA Rover')
	#plt.xlabel('Time[s]')
	#plt.ylabel('Velocity[m/s]')
	#plt.xlim([0,2.55])
	#plt.legend()
	#plt.show()
if __name__ == '__main__':
		main()
