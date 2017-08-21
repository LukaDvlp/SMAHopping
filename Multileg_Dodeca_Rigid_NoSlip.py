import numpy as np
from numpy import sin, cos, pi
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Simu0619

m1 = 0.6 + 0.05*11  #Mass of Body [Kg]
m2 = 0.05 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
k1 = 1200   #Bias Spring Coefficient [N/m]
c1 = 1  #Damping Coefficient of Spring [Ns/mm] 
l0 = 0.1  #Natural Length of Bias Spring [m]
Z20 = 0.0  #Initial Position of Pad [mm]
DZ = 0.05    #Initial Deflextion of Bias Spring [mm]
h = 0.0001   #Interval of RK
d = 1    #Diameter of SMA wire [mm]
D = 6.8  #Diameter of SMA coil [mm]
n = 10   #Number of coil spring
del_AE = Simu0619.delta_AE / 1000
del_s = Simu0619.delta_s / 1000
del_ME = Simu0619.delta_ME / 1000
F_s = Simu0619.F_s 
A = 1/Simu0619.A * 1000 
L = 0.1
rad32 = math.radians(32)
rad64 = math.radians(64)
rad112 = math.radians(112)
rad212 = math.radians(212)
rad244 = math.radians(244)
rad302 = math.radians(302) 
rad116 = math.radians(116)
R32 = np.matrix([[cos(rad32), -sin(rad32)],[sin(rad32), cos(rad32)]])
R64 = np.matrix([[cos(rad64), -sin(rad64)],[sin(rad64), cos(rad64)]])
R112 = np.matrix([[cos(rad112), -sin(rad112)],[sin(rad112), cos(rad112)]])
R212 = np.matrix([[cos(rad212), -sin(rad212)],[sin(rad212), cos(rad212)]])
R244 = np.matrix([[cos(rad244), -sin(rad244)],[sin(rad244), cos(rad244)]])
R302 = np.matrix([[cos(rad302), -sin(rad302)],[sin(rad302), cos(rad302)]])
th0 = np.arctan(L*sin(rad116)/(del_AE+L*cos(rad116))) 

def func1(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - A*(l-del_AE)
	print f
	print A*(l-del_AE)
	#f = k1 * (l0 -l) - c1 * dl 
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	term3 = 0
	term4 = 0 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])

def func2(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - F_s
	print f
	print F_s
	#f = k1 * (l0 -l) - c1 * dl
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	term3 = 0
	term4 = 0 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])


def func3(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - F_s
	print f
	print F_s 
	#f = k1 * (l0 -l) - c1 * dl 
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	term3 = -f*np.cos(th) / m2 
	term4 = -f*np.sin(th) / m2 -g 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])


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

	YY = np.empty((0,14), float)
	P1 = np.array([[del_AE*cos(th0)],[del_AE*sin(th0)]])
	P3 = 1/del_AE * P1 * L + P1
	P4 = 1/del_AE * R32 * P1 * L*sin(rad32) + P1
	P5 = 1/del_AE * R64 * P1 * L + P1
	P6 = 1/del_AE * R112 * P1 * L*cos(rad32) + P1
	P7 = 1/del_AE * R212 * P1 * L*sin(rad32) + P1
	P8 = 1/del_AE * R244 * P1 * L + P1
	P9 = 1/del_AE * R302 * P1 * L*cos(rad32) + P1
	YY = np.append(YY, np.array([[P3[0][0],P3[1][0],P4[0][0],P4[1][0],P5[0][0],P5[1][0]\
	              ,P6[0][0],P6[1][0],P7[0][0],P7[1][0],P8[0][0],P8[1][0],P9[0][0],P9[1][0]]]), axis=0)

	t = t_s
	n = 0
	while(t<t_f):
			l = np.sqrt((XX[n,0]-XX[n,4])**2 + (XX[n,2]-XX[n,6])**2)
			if l - del_AE < del_s:
				Step = RK(XX[n], func1, l)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				P1 = np.array([[Step[0]],[Step[2]]]) 
				P2 = np.array([[Step[4]], [Step[6]]]) 
				P3 = 1/del_AE * P1 * L + P1
				P4 = 1/del_AE * R32 * P1 * L*sin(rad32) + P1
				P5 = 1/del_AE * R64 * P1 * L + P1
				P6 = 1/del_AE * R112 * P1 * L*cos(rad32) + P1
				P7 = 1/del_AE * R212 * P1 * L*sin(rad32) + P1
				P8 = 1/del_AE * R244 * P1 * L + P1
				P9 = 1/del_AE * R302 * P1 * L*cos(rad32) + P1

				Sy = np.array([[P3[0][0],P3[1][0],P4[0][0],P4[1][0],P5[0][0],P5[1][0],P6[0][0],P6[1][0]\
				              ,P7[0][0],P7[1][0],P8[0][0],P8[1][0],P9[0][0],P9[1][0]]])
				YY = np.append(YY, Sy, axis=0)
				print(n)
				print("Using func1")
				t = t+h
				n = n+1
			elif l - del_AE >= del_s and l -del_AE < del_ME:
				Step = RK(XX[n], func2, l)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				P1 = np.array([[Step[0]],[Step[2]]]) 
				P2 = np.array([[Step[4]], [Step[6]]]) 
				P3 = 1/del_AE * P1 * L + P1
				P4 = 1/del_AE * R32 * P1 * L*sin(rad32) + P1
				P5 = 1/del_AE * R64 * P1 * L + P1
				P6 = 1/del_AE * R112 * P1 * L*cos(rad32) + P1
				P7 = 1/del_AE * R212 * P1 * L*sin(rad32) + P1
				P8 = 1/del_AE * R244 * P1 * L + P1
				P9 = 1/del_AE * R302 * P1 * L*cos(rad32) + P1

				Sy = np.array([[P3[0][0],P3[1][0],P4[0][0],P4[1][0],P5[0][0],P5[1][0],P6[0][0],P6[1][0]\
				              ,P7[0][0],P7[1][0],P8[0][0],P8[1][0],P9[0][0],P9[1][0]]])
				YY = np.append(YY, Sy, axis=0)
				print(n)
				print("Using func2")
				t = t+h
				n = n+1
			else:	
				break

	while(t<t_f):
			l = np.sqrt((XX[n,0]-XX[n,4])**2 + (XX[n,2]-XX[n,6])**2)
			if XX[n,2] >= 0 and XX[n,6] >= 0 and YY[n,1]>=0 and YY[n,3]>=0 and YY[n,5]>=0:
				Step = RK(XX[n], func3, l)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
			               Step[4],Step[5],Step[6],Step[7]]])
				XX = np.append(XX, S, axis=0)
				P1 = np.array([[Step[0]],[Step[2]]]) 
				P2 = np.array([[Step[4]], [Step[6]]]) 
				P3 = 1/del_AE * P1 * L + P1
				P4 = 1/del_AE * R32 * P1 * L*sin(rad32) + P1
				P5 = 1/del_AE * R64 * P1 * L + P1
				P6 = 1/del_AE * R112 * P1 * L*cos(rad32) + P1
				P7 = 1/del_AE * R212 * P1 * L*sin(rad32) + P1
				P8 = 1/del_AE * R244 * P1 * L + P1
				P9 = 1/del_AE * R302 * P1 * L*cos(rad32) + P1

				Sy = np.array([[P3[0][0],P3[1][0],P4[0][0],P4[1][0],P5[0][0],P5[1][0],P6[0][0],P6[1][0]\
				              ,P7[0][0],P7[1][0],P8[0][0],P8[1][0],P9[0][0],P9[1][0]]])
				YY = np.append(YY, Sy, axis=0)
				print(n)
				print("Using func3")
				t = t+h
				n = n+1
			else:
				break

	return [np.matrix(XX), np.matrix(YY)]

def main():
	t_s = 0.0
	t_f = 8.0 
	#th0 = np.pi/4
	X0 = [del_AE*np.cos(th0),0,del_AE*np.sin(th0),0,0,0,0,0]
	print(X0)
	XX,YY = Cal_Mtlx(X0, t_s, t_f, 8)
	print(XX)

	#print(Time)
	#T = np.arange(0, t_f+2*h, h)
	print("Size of XX")
	row, col = XX.shape
	print(row, col)
	T = np.arange(0, row*h, h)
	print("Length of T={0}".format(len(T)))
	print("del_AE={0}\n del_s={1}\n, del_ME={2}\n, F_s={3}\n, A={4}" .format(del_AE,del_s,del_ME, F_s, A))

	LimMax = 2.2
	LimMaxY = 1.2
	plt.figure()
	plt.title('Hopping Distance X,Z w.r.t Time')
	plt.xlabel('Time[s]')
	plt.ylabel('Hopping Height[m]')
	plt.plot(T,XX[:,0], label="Body_x")
	plt.plot(T,XX[:,2], label="Body_z")
	plt.plot(T,XX[:,4], label="Pad_x")
	plt.plot(T,XX[:,6], label="Pad_z")
	plt.xlim([0,T[row-1]])
	plt.legend()
	plt.savefig("cube_Rigid_NoSlip_T-XZ.png")
	
	plt.figure()
	plt.plot(XX[:,0],XX[:,2], label="Body")
	plt.title('Hopping Trajectory Z w.r.t X')
	plt.xlabel('Hopping Distance [m]')
	plt.ylabel('Hopping Height [m]')
	plt.xlim([0,XX[row-1,0]])
	plt.ylim([0,LimMaxY])
	plt.legend()
	plt.savefig("Cube_Rigid_NoSlip_XZ.png")
	#plt.show()

	fig = plt.figure()
	#ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.05, XX[row-1,0]+0.1),\
	#                                              ylim=(-0.05, XX[row-1,0]+0.1))
	ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.05,LimMax),\
	                                              ylim=(-0.05, LimMax))
	ax.set_xlabel=('Hello')
	ax.grid()
	
	line1, = ax.plot([], [], '-o', lw=3, color='blue')
	line2, = ax.plot([], [], '-o', lw=3, color='blue')
	line3, = ax.plot([], [], '-o', lw=3, color='blue')
	line4, = ax.plot([], [], '-o', lw=3, color='blue')
	time_template = 'time = %.1fs'
	time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
	
	
	def init():
	    line1.set_data([], [])
	    line2.set_data([], [])
	    line3.set_data([], [])
	    line4.set_data([], [])
	    time_text.set_text('')
	    return line1, line2, line3,line4, time_text
	
	
	TM_thin = 30	
	def animate(i):
	    thisx1 = [XX[i*TM_thin,0], XX[i*TM_thin,4]]
	    thisy1 = [XX[i*TM_thin,2], XX[i*TM_thin,6]]
	    thisx2 = [XX[i*TM_thin,0], YY[i*TM_thin,0]]
	    thisy2 = [XX[i*TM_thin,2], YY[i*TM_thin,1]]
	    thisx3 = [XX[i*TM_thin,0], YY[i*TM_thin,2]]
	    thisy3 = [XX[i*TM_thin,2], YY[i*TM_thin,3]]
	    thisx4 = [XX[i*TM_thin,0], YY[i*TM_thin,4]]
	    thisy4 = [XX[i*TM_thin,2], YY[i*TM_thin,5]]
	
	    #line.set_data(thisx, thisy)
	    line1.set_data(thisx1, thisy1)
	    line2.set_data(thisx2, thisy2)
	    line3.set_data(thisx3, thisy3)
	    line4.set_data(thisx4, thisy4)
	    time_text.set_text(time_template % (i*30*h))
	    return line1,line2,line3,line4, time_text
	
	ani = animation.FuncAnimation(fig, animate, np.arange(1, len(T)/30),
	                              interval=1, blit=False, init_func=init)
	
	#ani.save("Cube_Rigid_NoSlip.gif", writer = 'imagemagick')
	plt.show()
if __name__ == '__main__':
		main()
