import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

h = 0.1   #Interval of RK
d = 0.001    #Diameter of SMA wire [m]
D = 0.0068  #Diameter of SMA coil [m]
n = 10   #Number of coil spring
number = 2 #number of used SMA
rho = 6.5*10**3 #Density of SMA wire [kg/m^3]
r = math.sqrt(math.pi*d*D*n/4) #radius of pseudo sphere of SMA[mm]
print("radius of pseudo sphere of SMA[mm]")
print(r)
A = math.pi*r**2 #Cross section of the pseudo sphere of SMA[mm^2]
print("Cross section of SMA pseudo sphere [mm^2]")
print(A)
#A = 1
AA = 4*math.pi*r**2 #Surface are of SMA pseudo sphere [mm^2]
print("Surface area of SMA pseudo sphere [mm^2]")
print(AA)
#AA = 1
Theta_L = 45  #Lattitude [deg]
m = number*rho*(math.pi*(0.5*d))**2*D*n #Mass of a SMA spring [kg]
c = 440  #Specific heat capacity of SMA spring [J/kg K]
Sc = 1366 #Solar constant [W/m^2]
delta = 5.67*10**-8 #Stefan-Boltzman constant [W/m^2K^4]
a = 0.07 #Albedo constant
epsilon = 0.95 #Emissivity of the ground
deltaB = 1
Sb = 0.03*0.03  #Area of heat sink [m^2] 
theta = 45 #angle of hopping leg 
Fb_g = 0.5*(1+math.cos(math.radians(theta)))
Fb_s = 0.5*(1-math.cos(math.radians(theta)))
alpha = 1.0 #Emissivity of heat sink
R = 1 #Resistance of SMA wire

def ThermalEq1(x): #where Tg >= Tb
	value1 = A*math.cos(math.radians(Theta_L))*Sc/(m*c)
	print("value1")
	print(value1)
	value2 = a*A*math.cos(math.radians(Theta_L))*Sc/(m*c)
	print("value2")
	print(value2)
	value3 = epsilon*delta*A*(Tg**4-x[0]**4)/(2*m*c)
	print"(Tg/TB)**4 = {0}" 
	print(Tg/TB)
	print"(Tb/TB)**4 = {0}".format((x[0]/TB)**4)
	print("value3")
	print(value3)
	value4 = delta*AA*(x[0]**4-4**4)/(m*c)
	print("value4")
	print(value4)
	q = value1 + value2 + value3 - value4 
#	q = value1 + value2 
	return q

def ThermalEq2(x):   #where Tb > Tg
	value1 = A*math.cos(math.radians(Theta_L))*Sc/(m*c)
	print("value1")
	print(value1)
	value2 = a*A*math.cos(math.radians(Theta_L))*Sc/(m*c)
	print("value2")
	print(value2)
	value3 = delta*A*((x[0])**2 + (Tg)**2)*((x[0])**2 - (Tg)**2)/(2*m*c)
	print("value3")
	print(value3)
	value4 = delta*AA*((x[0]/TB)**4-(4/TB)**4)/(m*c)
	#value4 = delta*AA*((x[0]/TB)**4-(4/TB)**4)/(m*c)
	print("value4")
	print(value4)
#	q = value1 + value2 - value3 - value4 
	q = value1 - value3 - value4
	return q

def T_Eq1(x, Tg, I):   #Where Tg>Tb, with plate
	value1 = Sb*Fb_g*epsilon*delta*(Tg**4-x[0]**4)/(2*m*c)
	value2 = Sb*Fb_s*alpha*delta*(x[0]**4-4**4)/(m*c)
	value3 = R*I**2/(m*c) 
	print("value1={0}".format(value1))
	print("value2={0}".format(value2))
	print("value3={0}".format(value3))
	q = value1-value2+value3
	return q

def T_Eq2(x, Tg, I):   #Where Tb>Tg, with plate
	value1 = Sb*Fb_g*alpha*delta*(x[0]**4-Tg**4)/(m*c)
	value2 = Sb*Fb_s*alpha*delta*(x[0]**4-4**4)/(m*c)
	value3 = R*I**2/(m*c) 
	print("value1={0}".format(value1))
	print("value2={0}".format(value2))
	print("value3={0}".format(value3))
	q = -value1-value2+value3
	return q

def RK(x,f):  
	k1 = f(x)
	k2 = f(x+0.5*h*k1)
	k3 = f(x+0.5*h*k2)
	k4 = f(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_

def RK1(x, Tg, I): #Where Tg >= Tb 
	k1 = T_Eq1(x, Tg, I)
	k2 = T_Eq1(x+0.5*h*k1, Tg, I)
	k3 = T_Eq1(x+0.5*h*k2, Tg, I)
	k4 = T_Eq1(x+h*k3,Tg, I)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_

def RK2(x, Tg, I): #Where Tb > Tg 
	k1 = T_Eq2(x,Tg, I)
	k2 = T_Eq2(x+0.5*h*k1,Tg, I)
	k3 = T_Eq2(x+0.5*h*k2, Tg, I)
	k4 = T_Eq2(x+h*k3, Tg, I)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_

def Cal_Mtlx(X0, t_s, t_f, l, Tg):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			if XX[n] <= Tg:
					if XX[n] < 80+273:
							I = 2
							Step = RK1(XX[n], Tg, I)		
							print(n)
							S = np.array([[Step[0]]])
							XX = np.append(XX, S, axis=0)
							t = t+h
							n = n+1
					else:	
							I = 0
							Step = RK1(XX[n], Tg, I)		
							print(n)
							S = np.array([[Step[0]]])
							XX = np.append(XX, S, axis=0)
							t = t+h
							n = n+1
			else:
					if XX[n] < 80+273:
							I =1 
							Step = RK2(XX[n], Tg, I)		
							print(n)
							S = np.array([[Step[0]]])
							XX = np.append(XX, S, axis=0)
							t = t+h
							n = n+1
					else:	
							I =0 
							Step = RK2(XX[n], Tg, I)		
							print(n)
							S = np.array([[Step[0]]])
							XX = np.append(XX, S, axis=0)
							t = t+h
							n = n+1

	return (XX)

def main():
	t_s = 0.0
	t_f = 2000.0 
	t = 0
	n = 0
	X0 = [20+273]
	print(X0)
#	Tg = 120+273
#	XX1 = Cal_Mtlx(X0, t_s, t_f, 1, Tg)
#	Tg = 80+273
#	XX2 = Cal_Mtlx(X0, t_s, t_f, 1, Tg)
#	Tg = 50+273
#	XX3 = Cal_Mtlx(X0, t_s, t_f, 1, Tg)
#	Tg = 20+273
#	XX4 = Cal_Mtlx(X0, t_s, t_f, 1, Tg)
	Tg = 0+273
	XX5 = Cal_Mtlx(X0, t_s, t_f, 1, Tg)
	Time = []
	print(XX5)
#	print(XX2)
	T = np.arange(0, t_f+2*h, h)
	#T = np.arange(0, t_f+h, h)
	#T = 0.01*T
	print(T)
	print("mass of SMA={0}".format(m))
	print("Length of T = {0}".format(len(T)))
	print("size of XX5")
	rows1, cols1 = XX5.shape
	#print(rows1,cols1)
	#print("size of XX2")
	#rows2, cols2 = XX2.shape
	#print(rows2,cols2)
#	plt.plot(T,XX1[:,0], label="Tg:393[K]")
#	plt.plot(T,XX2[:,0], label="Tg:353[K]")
#	plt.plot(T,XX3[:,0], label="Tg:323[K]")
#	plt.plot(T,XX4[:,0], label="Tg:293[K]")
	plt.plot(T,XX5[:,0], label="Tg:273[K]")
	plt.xlabel('Time [s]')
	plt.ylabel('Temperature [s]')
	plt.title('Thermal change of SMA')
	plt.xlim([0,t_f])
	plt.legend()
	plt.show()

if __name__ == '__main__':
		main()
