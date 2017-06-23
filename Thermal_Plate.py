import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

h = 0.1   #Interval of RK
#d = 1    #Diameter of SMA wire [mm]
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
#m = 0.003 
#c = 440  #Specific heat capacity of SMA spring [J/kg K]
c = 440  #Specific heat capacity of SMA spring [J/kg K]
Sc = 1366 #Solar constant [W/m^2]
#Sc = 1
delta = 5.67*10**-8 #Stefan-Boltzman constant [W/m^2K^4]
#delta = 1
a = 0.07 #Albedo constant
#a = 1
Tg = 80 + 273.0 #Temperature of the ground[K]
#Tg = 1
epsilon = 0.6 #Emissivity of the ground
#ThetaB = 1
#epsilonB = 1
deltaB = 1
Sb = 0.02*0.02  #Area of heat sink [m^2] 
theta = 30 #angle of hopping leg 
Fb_g = 0.5*(1+math.cos(math.radians(theta)))
Fb_s = 0.5*(1-math.cos(math.radians(theta)))
alpha = 1.0 #Emissivity of heat sink


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
#	value4 = delta*A*((x[0]/TB)**4-(4/TB)**4)/(m*c)
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

def T_Eq1(x):   #Where Tg>Tb, with plate
	value1 = Sb*Fb_g*epsilon*delta*(Tg**4-x[0]**4)/(2*m*c)
	value2 = Sb*Fb_s*alpha*delta*(x[0]**4-4**4)/(m*c)
	print("value1={0}".format(value1))
	print("value2={0}".format(value2))
	q = value1-value2
	return q

def T_Eq2(x):   #Where Tg>Tb, with plate
	value1 = Sb*Fb_g*alpha*delta*(x[0]**4-Tg**4)/(m*c)
	value2 = Sb*Fb_s*alpha*delta*(x[0]**4-4**4)/(m*c)
	print("value1={0}".format(value1))
	print("value2={0}".format(value2))
	q = value1-value2

def RK(x,f):  
	k1 = f(x)
	k2 = f(x+0.5*h*k1)
	k3 = f(x+0.5*h*k2)
	k4 = f(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_

def RK1(x): #Where Tg >= Tb 
	k1 = ThermalEq1(x)
	k2 = ThermalEq1(x+0.5*h*k1)
	k3 = ThermalEq1(x+0.5*h*k2)
	k4 = ThermalEq1(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
#'''
def RK2(x): #Where Tb > Tg 
	k1 = ThermalEq2(x)
	k2 = ThermalEq2(x+0.5*h*k1)
	k3 = ThermalEq2(x+0.5*h*k2)
	k4 = ThermalEq2(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_

def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			if XX[n] <= Tg:
				Step = RK(XX[n], T_Eq1)
				print(n)
				S = np.array([[Step[0]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				t = t+h
				n = n+1
			else:
				Step = RK(XX[n], T_Eq2)
				print(n)
				S = np.array([[Step[0]]])
				#S = np.array(Step)
				XX = np.append(XX, S, axis=0)
				t = t+h
				n = n+1

	return (XX)

def main():
	t_s = 0.0
	t_f = 6000.0 
	t = 0
	n = 0
	X0 = [80+273]
	print(X0)
	XX = Cal_Mtlx(X0, t_s, t_f, 1)
	Time = []
	print(XX)
	T = np.arange(0, t_f+2*h, h)
	#T = np.arange(0, t_f+h, h)
	#T = 0.01*T
	print(T)
	print("mass of SMA={0}".format(m))
	print("Length of T = {0}".format(len(T)))
	print("size of XX")
	rows, cols = XX.shape
	print(rows,cols)
	plt.plot(T,XX[:,0], label="Temperature of SMA")
	plt.legend()
	plt.show()

if __name__ == '__main__':
		main()
