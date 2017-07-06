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
h = 0.002   #Interval of RK
A = 1/Simu0619.A
print(Simu0619.n)
del_AE = Simu0619.delta_AE/1000
del_s = Simu0619.delta_s/1000
F_s = Simu0619.F_s
theta = 45 #degree
theta = math.radians(theta) #Convert to radian
mu_d = 0.8
I = 0.1 
M_G = M1+M2
r2 = 1

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

def TwoD_func1_NoSlip(x):
	term1 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.cos(theta)
	term3 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	term4 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta) + c1/M2*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) -k2/M2*x[6] -c2/M2*x[7] - g
	if x[6]<=0:
			term2 = 0
	else:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)
	print"term1 = {0}".format(term1)
	print"term2 = {0}".format(term2)
	print"term3 = {0}".format(term3)
	print"term4 = {0}".format(term4)
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7], term4])
	
def TwoD_func1_Slip(x):
	term1 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.cos(theta)
	term3 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	term4 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta) + c1/M2*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) -k2/M2*x[6] -c2/M2*x[7] - g
	if x[6]<=0:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta) + mu_d*(-k2/M2*x[6] -c2/M2*x[7])

	else:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)
	print"term1 = {0}".format(term1)
	print"term2 = {0}".format(term2)
	print"term3 = {0}".format(term3)
	print"term4 = {0}".format(term4)
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7], term4])
	
def TwoD_func2_NoSlip(x):
	term1 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.cos(theta)
	term3 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	term4 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta) + c1/M2*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) -k2/M2*x[6] - g
	if x[6]<=0:
			term2 = 0
	else:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)
	print"term1 = {0}".format(term1)
	print"term2 = {0}".format(term2)
	print"term3 = {0}".format(term3)
	print"term4 = {0}".format(term4)
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7], term4])
	
def TwoD_func2_Slip(x):
	term1 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.cos(theta)
	term3 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	term4 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta) + c1/M2*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) -k2/M2*x[6] - g
	if x[6]<=0:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta) + mu_d*(-k2/M2*x[6] -c2/M2*x[7])

	else:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)
	print"term1 = {0}".format(term1)
	print"term2 = {0}".format(term2)
	print"term3 = {0}".format(term3)
	print"term4 = {0}".format(term4)
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7], term4])

def TwoD_func3_NoSlip(x):
	term1 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.cos(theta)
	term3 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	term4 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta) + c1/M2*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	if x[6]<=0:
			term2 = 0
	else:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)
	print"term1 = {0}".format(term1)
	print"term2 = {0}".format(term2)
	print"term3 = {0}".format(term3)
	print"term4 = {0}".format(term4)
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7], term4])

def TwoD_func3_Slip(x):
	term1 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.cos(theta)
	term3 = k1/M1*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta)	- c1/M1*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	term4 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.sin(theta) + c1/M2*math.sqrt((x[1]-x[3])**2+(x[5]-x[7])**2)*math.sin(theta) - g
	if x[6]<=0:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta) + mu_d*(-k2/M2*x[6] -c2/M2*x[7])

	else:
			term2 = -k1/M2*(l0-math.sqrt((x[0]-x[2])**2+(x[4]-x[6])**2))*math.cos(theta)
	print"term1 = {0}".format(term1)
	print"term2 = {0}".format(term2)
	print"term3 = {0}".format(term3)
	print"term4 = {0}".format(term4)
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7], term4])

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

def MotionOfCOG(x, F1, F2):	
	term1 = F2/M_G
	term2 = F1/M_G -g
	term3 = M2*r2*F1*math.cos(x[4])/I + M2*r2*F2*math.sin(x[4])/I 
	print(term3)
	return np.array([x[1], term1, x[3], term2, term3])

def RK(x, f, F1, F2):
	k1 = f(x, F1, F2)
	k2 = f(x+0.5*h*k1, F1, F2)
	k3 = f(x+0.5*h*k2, F1, F2)
	k4 = f(x+h*k3, F1, F2)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3],X0[4],X0[5],X0[6],X0[7]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			if XX[n,6]<0 and XX[n,7]<0:
				Step = RK(XX[n], TwoD_func1_Slip)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func1")
				t = t+h
				n = n+1
			elif XX[n,6]<0 and XX[n,7]>=0:
				Step = RK(XX[n], TwoD_func2_Slip)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func2")
				t = t+h
				n = n+1
			else:	
				Step = RK(XX[n], TwoD_func3_Slip)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func3")
				t = t+h
				n = n+1

	return np.matrix(XX)
	
def Cal_Mtl_NoSMA(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3]]]), axis=0)
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

def Cal_Test(X0, t_s, t_f, t_i, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3], X0[4]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			if t<t_i:
				F1=40
				F2=20
				Step = RK(XX[n], MotionOfCOG, F1, F2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3], Step[4]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				t = t+h
				n = n+1
			else:
				F1=0
				F2=0
				Step = RK(XX[n], MotionOfCOG, F1, F2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3], Step[4]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				t = t+h
				n = n+1

	return np.matrix(XX)

def main():
	t_s = 0.0
	t_f = 10.0 
	t_i = 0.2
	v0 = 10
	H = 1
	theta_g = math.pi/4
	X0 = [math.cos(theta_g),0,math.sin(theta_g),0,theta_g]
	#X0 = [l0-0.05,0]
	print(X0)
	XX = Cal_Test(X0, t_s, t_f, t_i, 5)
	print(XX)
	T = np.arange(0, t_f+h, h)
	print("Length of T={0}".format(len(T)))
	print("Size of XX")
	row, col = XX.shape
	print(row, col)

	plt.plot(T,XX[:,0], label="Position X")
	plt.plot(T,XX[:,2], label="Velocity X")
	plt.plot(T,XX[:,2], label="Position Z")
	plt.plot(T,XX[:,3], label="Velocity Z")
	plt.plot(T,XX[:,4], label="Angle")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.title('2-Dimensional Behavior of the COG')
	plt.xlabel('Time[s]')
	plt.ylabel('Values')
	#plt.xlim([0,2.16])
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
