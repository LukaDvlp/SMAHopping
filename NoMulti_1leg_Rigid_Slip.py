#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from numpy import sin,cos
import matplotlib.pyplot as plt
import math

m1 = 0.65  #Mass of Body [Kg]
m2 = 0.05 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
#g = 9.8
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
r1 = m2*l0/(m1+m2)
r2 = m1*l0/(m1+m2)
#I = m1*r1**2 + m2*r2**2
I = m1*m2/(m1+m2)*l0**2
m = m1+m2

def func1(x,l,th):
	f = k1*(l0-l) 
	term1 = f*cos(th)/m
	term2 = f*sin(th)/m - g
	term3 = -l*f*sin(th)/sin(ita)*sin(ita-th)/I
	print("term3 = {0}".format(term3))
	return np.array([x[1], term1, x[3],term2, x[5],term3]) 

def func2(x,l,th):
	term1 = 0
	term2 = -g
	term3 = 0
	#print("term2 = {0}".format(term2))
	return np.array([x[1], term1, x[3],term2, x[5],term3]) 

def trans(XX,XA_):
	R = np.matrix([[cos(XX[4]),-sin(XX[4])],[sin(XX[4]),cos(XX[4])]])
	return R*XA_ + np.matrix([[XX[0]],[XX[2]]])

def RK(x,l,th,f):
	k1 = f(x,l,th)
	k2 = f(x+0.5*h*k1,l,th)
	k3 = f(x+0.5*h*k2,l,th)
	k4 = f(x+h*k3,l,th)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, s):
	t = t_s
	n = 0
	XX = np.empty((0,s), float)
	XA_ = np.matrix([[r1],[0]])
	XA = np.empty((0,2),float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3],X0[4],X0[5]]])\
	              , axis=0)
	#XX = np.append(XX, np.array([[X0[0],X0[1]]]), axis=0)
	XA = trans(XX[n],XA_)
	while(t<t_f):
			l = np.sqrt(XX[n,0]**2 + XX[n,2]**2)
			th = XX[n,4]
			if l0-l > 0:
				Step = RK(XX[n],l,th,func1)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5]]])
				XX = np.append(XX, S, axis=0)
				SA = trans(XX[n+1],XA_)
				XA = np.append(XA,SA,axis=1)
				print(n)
				print("Using func1")
				t = t+h
				n = n+1
				k = t
			elif XX[n,2] < 0:
				break
			else:	
				Step = RK(XX[n],l,th,func2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5]]])
				XX = np.append(XX, S, axis=0)
				SA = trans(XX[n+1],XA_)
				XA = np.append(XA,SA,axis=1)
				print(n)
				print("Using func2")
				t = t+h
				n = n+1

	return [np.matrix(XX),np.matrix(XA),k]

def main():
	t_s = 0.0
	t_f = 10.0 
	dz = 0.05
	th0 = np.pi/4
	mu = 0.3
	global ita
	ita  = np.arctan(1/mu)
	x0 = (l0-dz)*cos(th0)
	y0 = (l0-dz)*sin(th0)
	#X0 = [del_AE,0,0.0,0]
	X0 = [x0,0,y0,0,th0,0]
	XX,XA,k = Cal_Mtlx(X0, t_s, t_f, 6)
	#k = Cal_Mtlx(X0, t_s, t_f, 6,th)[1]
	XA = XA.T
	print(XX)
	print("Launch Time was t={0}[s]".format(k))

	print("Size of XX")
	row, col = XX.shape
	print(row, col)
	print("Size of XA")
	row_a, col_a = XA.shape
	print(row_a, col_a)
#	T = np.arange(0, row, h)
	T = np.arange(0, row_a*h, h)
	print("Length of T={0}".format(len(T)))

	plt.figure()
	plt.plot(T,XA[:,0], label="X")
	plt.plot(T,XA[:,1], label="Z")
	plt.figure()
	plt.plot(T,XX[:,4], label="theta")
	plt.plot(T,XX[:,5], label="omega")
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
