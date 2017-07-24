#!/usr/bin/env python
# -*- coding: utf-8 -*-
#=============================================================
#DPend2: Dynamic analysis program for Double Rigid Pendulum
#=============================================================
import numpy as np
from numpy import sin,cos
import matplotlib.pyplot as plt

#Setting of Parameters
l0 = 0.1 
Df = 0.05
m1,m2,g,=0.65,0.05,1.62
r1,r2 = 0.05,0.05
#I1,I2 = 2/5*m1*r1**3, 2/5*m2*r2**3
I1,I2 = 0,0
k1,k2 = 1400,1000

#Initial State of the machanism
th10,th20 = -np.pi/4, -np.pi/4
x0 = np.array([[(l0-Df)*cos(th10),0,(l0-Df)*sin(th10),0,th10,0, \
               0,0,0,0, th20,0]])

#Time Counter
Ts,Te,h=0,5,0.0001

def Eq_Motion(x,Qa,Qc):
	term1 = (Qa[0,0]+Qc[0,0])/m1
	term2 = (Qa[0,1]+Qc[0,1])/m1
	term3 = (Qa[0,2]+Qc[0,2])/(I1)
	term4 = (Qa[0,3]+Qc[0,3])/m2
	term5 = (Qa[0,4]+Qc[0,4])/m2
	term6 = (Qa[0,5]+Qc[0,5])/(I2)
	return np.array([x[1],term1,x[3],term2,x[5],term3,x[7],term4,x[9],term5,x[11],term6])

def RK(x,f,Qa,Qc):
	k1 = f(x, Qa, Qc)
	k2 = f(x+0.5*h*k1, Qa, Qc)
	k3 = f(x+0.5*h*k2, Qa, Qc)
	k4 = f(x+h*k3, Qa, Qc)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	return x_

def trans(XX,XA_):
	R = np.matrix([[cos(XX[10]),-sin(XX[10])],[sin(XX[10]),cos(XX[10])]])
	return R*XA_ + np.matrix([[XX[6]],[XX[8]]])

n = 0
t = Ts
XA_ = np.matrix([[1],[0]])
XX = np.empty((0,12),float)
XA = np.empty((0,2),float)
XX = np.append(XX,x0, axis=0)
XA = trans(XX[n],XA_)
while(t<Te):
	#Fk = k1*(l0-np.sqrt((XX[n,0]-XX[n,6])**2 + (XX[n,2]-XX[n,8])**2))
	l = np.sqrt((XX[n,0]-XX[n,6])**2 + (XX[n,2]-XX[n,8])**2)
	f = k1*(l0-l)
	sth1,cth1=sin(XX[n,4]),cos(XX[n,4])
	sth2,cth2=sin(XX[n,10]),cos(XX[n,10])
	#Fkx = Fk*sth1
	#Fkx = Fk*sth2
	#Fky = Fk*cth1
	#Fky = Fk*cth2
	A = np.array([[m1,0,0,0,0,0,-sth1],[0,m1,0,0,0,0,cth1],\
	              [0,0,I1,0,0,0,-sth1*(XX[n,2]-XX[n,8])-cth1*(XX[n,0]-XX[n,6])],[0,0,0,m2,0,0,sth1], \
				  [0,0,0,0,m2,0,-cth1],[0,0,0,0,0,I2,0], \
				  [-sth1,cth1,-sth1*(XX[n,2]-XX[n,8])-cth1*(XX[n,0]-XX[n,6]),sth1,-cth1,0,0]]) 
	
	#The right side
	b = np.array([f/l*(XX[n,0]-XX[n,6]),f/l*(XX[n,2]-XX[n,8]),0,\
	             -f/l*(XX[n,0]-XX[n,6]),-f/l*(XX[n,2]-XX[n,8]),0, \
				  XX[n,5]*cth1*(XX[n,5]*(XX[n,2]-XX[n,8])+2*(XX[n,1]-XX[n,7])) \
				  - XX[n,5]*sth1*(XX[n,5]*(XX[n,0]-XX[n,6])-2*(XX[n,3]-XX[n,9]))])
	
	#Solution
	x = np.linalg.solve(A,b)
	
	Phi_q = np.array([[-sth1,cth1,-sth1*(XX[n,2]-XX[n,8])-cth1*(XX[n,0]-XX[n,6]),sth1,-cth1,0]])
	Qa = np.array([[f/l*(XX[n,0]-XX[n,6]),f/l*(XX[n,2]-XX[n,8]),0,\
	               -f/l*(XX[n,0]-XX[n,6]),-f/l*(XX[n,2]-XX[n,8]),0]])
	lam = np.matrix([[x[6]]]).T
	Q = -Phi_q.T*lam
	Qc = Q.T
	
	Step = RK(XX[n],Eq_Motion,Qa,Qc)
	S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],\
	               Step[6],Step[7],Step[8],Step[9],Step[10],Step[11]]])
	XX = np.append(XX,S,axis=0)
	SA = trans(XX[n+1],XA_)
	XA = np.append(XA,SA,axis=1)
	print(n)
	t = t+h
	n= n+1

XA = XA.T
#T = np.arange(0, Te+2*h, h)
T = np.arange(0, Te+h, h)
print("Length of T = {0}".format(len(T)))
print("size of XA")
rows1, cols1 = XA.shape
print(rows1,cols1)
print XA
plt.plot(T, XX[:,0], label="x1")
plt.plot(T, XX[:,2], label="z1")
plt.plot(T, XX[:,4], label="theta1")
plt.plot(T, XX[:,6], label="x2")
plt.plot(T, XX[:,8], label="z2")
#plt.plot(T, XX[:,10], label="theta2")
plt.legend()
plt.show()
