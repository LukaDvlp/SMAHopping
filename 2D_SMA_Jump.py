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
h = 0.0002   #Interval of RK
A = 1/Simu0619.A
print(Simu0619.n)
del_AE = Simu0619.delta_AE/1000
del_s = Simu0619.delta_s/1000
F_s = Simu0619.F_s
theta = 45 #degree
theta = math.radians(theta) #Convert to radian
mu_d = 0.8
I = 2*M1*M2*l0/(M1+M2) 
M_G = M1+M2
r2 = M1*l0/(M1+M2) 

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
	term3 = r2*F1*math.cos(x[4])/I + r2*F2*math.sin(x[4])/I 
	print(term3)
	return np.array([x[1], term1, x[3], term2, x[5], term3])

def Eq_Motion1(x):
	R = np.matrix([[np.sin(x[8]), np.cos(x[8])],[-np.cos(x[8]), np.sin(x[8])]])
	R_ = np.linalg.inv(R)
	l2 = M1/(M1+M2) *(x[2]-x[6])
	r2 = R_*np.matrix([[x[4]],[x[6]]])
	v2 = R_*np.matrix([[x[5]],[x[7]]])
	F1_ = k1*(l0-(x[2]-x[6])) - c1*(x[3]-x[7])
	f = np.matrix([[0], [M1*g]])
	FG1_ = R*f
	FG2_ = R*np.matrix([[0],[M2*g]])
	F3_ = R*np.matrix([[0],[-k2*r2[1,0]-c2*v2[1,0]]]) 
	F2_ = mu_d*F3_
	F2 = R_*F2_
	F3 = R_*F3_
	term_r = l2*F2[0,0]*np.cos(x[8]) + l2*F3[1,0]*np.sin(x[8]) 
	Fcent_x1 = M1*x[0]*(x[9])**2
	Fcent_z1 = M1*x[2]*(x[9])**2 
	Fcent_x2 = M2*x[4]*(x[9])**2
	Fcent_z2 = M2*x[6]*(x[9])**2
	Fcor_x1 = 2*M1*x[3]*x[9]
	Fcor_z1 = -2*M1*x[1]*x[9]
	Fcor_x2 = 2*M2*x[7]*x[9]
	Fcor_z2 = -2*M2*x[5]*x[9]
	Feul_x1 = M1*x[2]*term_r
	Feul_z1 = -M1*x[0]*term_r
	Feul_x2 = M2*x[6]*term_r
	Feul_z2= -M2*x[4]*term_r
	
	return np.array([x[1], FG1_[0,0]/M1+Fcent_x1/M1+Fcor_x1/M1+Feul_x1/M1, x[3], F1_/M1 +Fcent_z1/M1+Fcor_z1/M1+Feul_z1/M1+ FG1_[1,0]/M1, x[5], F2_[0,0]/M2 + F3_[0,0]/M2 + FG2_[0,0]/M2+Fcent_x2/M2+Fcor_x2/M2+Feul_x2/M2, x[7], F2_[1,0]/M2 + F3_[1,0]/M2 + FG2_[1,0]/M2 - F1_/M2+Fcent_z2/M2+Fcor_z2/M2+Feul_z2/M2, x[9], term_r ])

def Eq_Motion2(x):
	R = np.matrix([[np.sin(x[8]), np.cos(x[8])],[-np.cos(x[8]), np.sin(x[8])]])
	R_ = np.linalg.inv(R)
	l2 = M1/(M1+M2) *(x[2]-x[6])
	r2 = R_*np.matrix([[x[4]],[x[6]]])
	v2 = R_*np.matrix([[x[5]],[x[7]]])
	F1_ = k1*(l0-(x[2]-x[6])) - c1*(x[3]-x[7])
	FG1_ = R*np.matrix([[0],[M1*g]])
	FG2_ = R*np.matrix([[0],[M2*g]])
	F3_ = R*np.matrix([[0],[-k2*r2[1,0]]]) 
	F2_ = mu_d*F3_
	F2 = R_*F2_
	F3 = R_*F3_
	term_r = l2*F2[0,0]*np.cos(x[8]) + l2*F3[1,0]*np.sin(x[8]) 
	
	Fcent_x1 = M1*x[0]*(x[9])**2
	Fcent_z1 = M1*x[2]*(x[9])**2 
	Fcent_x2 = M2*x[4]*(x[9])**2
	Fcent_z2 = M2*x[6]*(x[9])**2
	Fcor_x1 = 2*M1*x[3]*x[9]
	Fcor_z1 = -2*M1*x[1]*x[9]
	Fcor_x2 = 2*M2*x[7]*x[9]
	Fcor_z2 = -2*M2*x[5]*x[9]
	Feul_x1 = M1*x[2]*term_r
	Feul_z1 = -M1*x[0]*term_r
	Feul_x2 = M2*x[6]*term_r
	Feul_z2= -M2*x[4]*term_r
	
	return np.array([x[1], FG1_[0,0]/M1+Fcent_x1/M1+Fcor_x1/M1+Feul_x1/M1, x[3], F1_/M1 +Fcent_z1/M1+Fcor_z1/M1+Feul_z1/M1+ FG1_[1,0]/M1, x[5], F2_[0,0]/M2 + F3_[0,0]/M2 + FG2_[0,0]/M2+Fcent_x2/M2+Fcor_x2/M2+Feul_x2/M2, x[7], F2_[1,0]/M2 + F3_[1,0]/M2 + FG2_[1,0]/M2 - F1_/M2+Fcent_z2/M2+Fcor_z2/M2+Feul_z2/M2, x[9], term_r ])
#	return np.array([x[1], FG1_[0,0]/M1, x[3], F1_/M1 + FG1_[1,0]/M1, x[5], F2_[0,0]/M2 + F3_[0,0]/M2 + FG2_[0,0]/M2, x[7], F2_[1,0]/M2 + F3_[1,0]/M2 + FG2_[1,0]/M2 - F1_/M2, x[9], term_r ])
	

def Eq_Motion3(x):
	R = np.matrix([[np.sin(x[8]), np.cos(x[8])],[-np.cos(x[8]), np.sin(x[8])]])
	R_ = np.linalg.inv(R)
	l2 = M1/(M1+M2) *(x[2]-x[6])
	r2 = R_*np.matrix([[x[4]],[x[6]]])
	v2 = R_*np.matrix([[x[5]],[x[7]]])
	F1_ = k1*(l0-(x[2]-x[6])) - c1*(x[3]-x[7])
	FG1_ = R*np.matrix([[0],[M1*g]])
	FG2_ = R*np.matrix([[0],[M2*g]])
	F3_ = R*np.matrix([[0],[0]]) 
	F2_ = mu_d*F3_
	F2 = R_*F2_
	F3 = R_*F3_
	term_r = l2*F2[0,0]*np.cos(x[8]) + l2*F3[1,0]*np.sin(x[8]) 
	Fcent_x1 = M1*x[0]*(x[9])**2
	Fcent_z1 = M1*x[2]*(x[9])**2 
	Fcent_x2 = M2*x[4]*(x[9])**2
	Fcent_z2 = M2*x[6]*(x[9])**2
	Fcor_x1 = 2*M1*x[3]*x[9]
	Fcor_z1 = -2*M1*x[1]*x[9]
	Fcor_x2 = 2*M2*x[7]*x[9]
	Fcor_z2 = -2*M2*x[5]*x[9]
	Feul_x1 = M1*x[2]*term_r
	Feul_z1 = -M1*x[0]*term_r
	Feul_x2 = M2*x[6]*term_r
	Feul_z2= -M2*x[4]*term_r
	
	return np.array([x[1], FG1_[0,0]/M1+Fcent_x1/M1+Fcor_x1/M1+Feul_x1/M1, x[3], F1_/M1 +Fcent_z1/M1+Fcor_z1/M1+Feul_z1/M1+ FG1_[1,0]/M1, x[5], F2_[0,0]/M2 + F3_[0,0]/M2 + FG2_[0,0]/M2+Fcent_x2/M2+Fcor_x2/M2+Feul_x2/M2, x[7], F2_[1,0]/M2 + F3_[1,0]/M2 + FG2_[1,0]/M2 - F1_/M2+Fcent_z2/M2+Fcor_z2/M2+Feul_z2/M2, x[9], term_r ])

#	return np.array([x[1], FG1_[0,0]/M1, x[3], F1_/M1 + FG1_[1,0]/M1, x[5], F2_[0,0]/M2 + F3_[0,0]/M2 + FG2_[0,0]/M2, x[7], F2_[1,0]/M2 + F3_[1,0]/M2 + FG2_[1,0]/M2 - F1_/M2, x[9], term_r ])


def RK(x, f):
	k1 = f(x)
	k2 = f(x+0.5*h*k1)
	k3 = f(x+0.5*h*k2)
	k4 = f(x+h*k3)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0,Y0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3],X0[4],X0[5],X0[6],X0[7], X0[8], X0[9]]]), axis=0)
	YY = np.empty((0,l-2), float)
	YY = np.append(YY, np.array([[Y0[0],Y0[1],Y0[2],Y0[3],Y0[4],Y0[5],Y0[6],Y0[7]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
		R = np.matrix([[np.sin(XX[n,8]), np.cos(XX[n,8])],[-np.cos(XX[n,8]), np.sin(XX[n,8])]])
		R_ = np.linalg.inv(R)
		r2 = R_*np.matrix([[XX[n,4]],[XX[n,6]]])
		v2 = R_*np.matrix([[XX[n,5]],[XX[n,7]]])
		if r2[1,0]<0 and v2[1,0]<0:
				Step = RK(XX[n], Eq_Motion1)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],Step[6],Step[7],Step[8],Step[9]]])
				Sxx1 = R_*np.matrix([[S[0,0]],[S[0,2]]])
				Sxx2 = R_*np.matrix([[S[0,4]],[S[0,6]]])
				Svv1 = R_*np.matrix([[S[0,1]],[S[0,3]]])
				Svv2 = R_*np.matrix([[S[0,5]],[S[0,7]]])
				Sy = np.array([[Sxx1[0,0],Svv1[0,0],Sxx1[1,0],Svv1[1,0],Sxx2[0,0],Svv2[0,0],Sxx2[1,0],Svv2[1,0]]])
				XX = np.append(XX, S, axis=0)
				YY = np.append(YY, Sy, axis=0)
				print(n)
				print("Using func1")
				t = t+h
				n = n+1
		elif r2[1,0]<0 and v2[1,0]>=0:
				Step = RK(XX[n], Eq_Motion2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],Step[6],Step[7],Step[8],Step[9]]])
				Sxx1 = R_*np.matrix([[S[0,0]],[S[0,2]]])
				Sxx2 = R_*np.matrix([[S[0,4]],[S[0,6]]])
				Svv1 = R_*np.matrix([[S[0,1]],[S[0,3]]])
				Svv2 = R_*np.matrix([[S[0,5]],[S[0,7]]])
				Sy = np.array([[Sxx1[0,0],Svv1[0,0],Sxx1[1,0],Svv1[1,0],Sxx2[0,0],Svv2[0,0],Sxx2[1,0],Svv2[1,0]]])
				XX = np.append(XX, S, axis=0)
				YY = np.append(YY, Sy, axis=0)

				print(n)
				print("Using func2")
				t = t+h
				n = n+1
		else:	
				Step = RK(XX[n], Eq_Motion3)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],Step[4],Step[5],Step[6],Step[7],Step[8],Step[9]]])
				Sxx1 = R_*np.matrix([[S[0,0]],[S[0,2]]])
				Sxx2 = R_*np.matrix([[S[0,4]],[S[0,6]]])
				Svv1 = R_*np.matrix([[S[0,1]],[S[0,3]]])
				Svv2 = R_*np.matrix([[S[0,5]],[S[0,7]]])
				Sy = np.array([[Sxx1[0,0],Svv1[0,0],Sxx1[1,0],Svv1[1,0],Sxx2[0,0],Svv2[0,0],Sxx2[1,0],Svv2[1,0]]])
				XX = np.append(XX, S, axis=0)
				YY = np.append(YY, Sy, axis=0)

				print(n)
				print("Using func3")
				t = t+h
				n = n+1

	XX = np.matrix(XX)
	YY = np.matrix(YY)
	return (XX, YY)
	

def main():
	t_s = 0.0
	t_f = 5.0 
	t_i = 0.1

	theta_g = math.pi/4
	R = np.matrix([[np.sin(theta_g), np.cos(theta_g)],[-np.cos(theta_g), np.sin(theta_g)]])
	R_ = np.linalg.inv(R)
	X0 = [0, 0, del_AE, 0, 0, 0, 0, 0, theta_g, 0]
	yy1 = R_*np.matrix([[X0[0]],[X0[3]]])
	Y0 = [yy1[0,0], 0, yy1[1,0], 0, 0, 0, 0, 0]
	print(X0)
	XX, YY = Cal_Mtlx(X0,Y0, t_s, t_f, 10)
	print(XX)
	T = np.arange(0, t_f+h, h)
	#T = np.arange(0, t_f+2*h, h)
	print("Length of T={0}".format(len(T)))
	print("Size of XX")
	row, col = XX.shape
	print(row, col)

	print("Length of T={0}".format(len(T)))
	print("Size of YY")
	row, col = YY.shape
	print(row, col)

	plt.plot(T,XX[:,0], label="Position X1\'")
	#plt.plot(T,XX[:,1], label="Velocity X1\'")
	plt.plot(T,XX[:,2], label="Position Z1\'")
	#plt.plot(T,XX[:,3], label="Velocity Z1\'")
	plt.plot(T,XX[:,4], label="Position X2\'")
	#plt.plot(T,XX[:,5], label="Velocity X2\'")
	plt.plot(T,XX[:,6], label="Position Z2\'")
	#plt.plot(T,XX[:,7], label="Velocity Z2\'")
	plt.plot(T,XX[:,8], label="Angle")
	plt.plot(T,XX[:,9], label="Angular Velocity")
	plt.title('2-Dimensional Behavior on the non-inertial flame')
	plt.xlabel('Time[s]')
	plt.ylabel('Values')
	#plt.xlim([0,2.16])
	plt.legend()
	plt.show()
	
	plt.plot(T,YY[:,0], label="Position X1")
	#plt.plot(T,YY[:,1], label="Velocity X1")
	plt.plot(T,YY[:,2], label="Position Z1")
	#plt.plot(T,YY[:,3], label="Velocity Z1")
	plt.plot(T,YY[:,4], label="Position X2")
	#plt.plot(T,YY[:,5], label="Velocity X2")
	plt.plot(T,YY[:,6], label="Position Z2")
	#plt.plot(T,YY[:,7], label="Velocity Z2")
	plt.plot(T,XX[:,8], label="Angle")
	plt.plot(T,XX[:,9], label="Angular Velocity")
	plt.title('2-Dimensional Behavior on the non-inertial flame')
	plt.xlabel('Time[s]')
	plt.ylabel('Values')
	#plt.xlim([0,2.16])
	plt.legend()
	plt.show()

if __name__ == '__main__':
		main()
