#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from numpy import sin,cos
import matplotlib.pyplot as plt
import math

m1 = 0.7  #Mass of Body [Kg]
M2 = 0.05 #Mass of Pad [Kg]
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

def func1(x,l,th):
	term1 = k1*(l0-l)*cos(th)/m1
	term2 = k1*(l0-l)*sin(th)/m1 - g
	return np.array([x[1], term1, x[3],term2]) 

def func2(x,l,th):
	term1 = 0
	term2 = -g
	#print("term2 = {0}".format(term2))
	return np.array([x[1], term1, x[3],term2]) 

def RK(x,l,th,f):
	k1 = f(x,l,th)
	k2 = f(x+0.5*h*k1,l,th)
	k3 = f(x+0.5*h*k2,l,th)
	k4 = f(x+h*k3,l,th)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, s,th):
	XX = np.empty((0,s), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3]]]), axis=0)
	#XX = np.append(XX, np.array([[X0[0],X0[1]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			l = np.sqrt(XX[n,0]**2 + XX[n,2]**2)
			if l0-l > 0:
				Step = RK(XX[n],l,th,func1)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func1")
				t = t+h
				n = n+1
				k = t
			elif XX[n,2] < 0:
				break
			else:	
				Step = RK(XX[n],l,th,func2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func2")
				t = t+h
				n = n+1

	return [np.matrix(XX),k]

def main():
	t_s = 0.0
	t_f = 3.0 
	dz = 0.05
	th = np.pi/4
	x0 = (l0-dz)*cos(th)
	y0 = (l0-dz)*sin(th)
	#X0 = [del_AE,0,0.0,0]
	X0 = [x0,0,y0,0]
	XX = Cal_Mtlx(X0, t_s, t_f, 4,th)[0]
	k = Cal_Mtlx(X0, t_s, t_f, 4,th)[1]
	print(XX)
	print("Launch Time was t={0}[s]".format(k))

	print("Size of XX")
	row, col = XX.shape
	print(row, col)
#	T = np.arange(0, row, h)
	T = np.arange(0, row*h, h)
	print("Length of T={0}".format(len(T)))

	plt.plot(T,XX[:,0], label="X")
	plt.plot(T,XX[:,2], label="Z")
	#plt.plot(T,XX[:,2], label="Position of X2")
	#plt.plot(T,XX[:,3], label="Velocity of X2")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	plt.title('2-Dimensional Hopping of Spring Only Rover')
	plt.xlabel('Time[s]')
	plt.ylabel('Hopping distance[m]')
	plt.xlim([0,T[row-1]])
	plt.legend()
	plt.show()
if __name__ == '__main__':
		main()
