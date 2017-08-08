import numpy as np
from numpy import sin, cos, pi
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Simu0619

m1 = 0.65  #Mass of Body [Kg]
m2 = 0.05 #Mass of Pad [Kg]
m3 = 0.05 #Mass of Pad [Kg]
g = 1.622  #Gravitational Acceleration [m/s^2] 
k1 = 1060   #Bias Spring Coefficient [N/m]
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
del_AE = Simu0619.delta_AE / 1000
del_s = Simu0619.delta_s / 1000
del_ME = Simu0619.delta_ME / 1000
F_s = Simu0619.F_s 
A = 1/Simu0619.A * 1000 

def func1(x, l1, l2):
	dl1 = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l1
	dl2 = ((x[8]-x[0])*(x[9]-x[1])+(x[10]-x[2])*(x[11]-x[3]))/l2
	th1 = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	th2 = np.arctan((x[10]-x[2])/(x[8]-x[0]))
	f1 = k1 * (l0 -l1) - c1 * dl1 - A*(l1-del_AE)
	f2 = k1 * (l0 -l2) - c1 * dl2 - A*(l1-del_AE)
	#print f
	#print A*(l-del_AE)
	term1 = (f1*np.cos(th1) - f2*cos(th2))/ m1
	term2 = (f1*np.sin(th1) - f2*sin(th2)) / m1 - g
	term3 = 0
	term4 = 0 
	term5 = f2*cos(th2) / m3 
	term6 = f2*sin(th2) / m3 -g 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4, x[9],term5,x[11],term6])

def func2(x, l1, l2):
	dl1 = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l1
	dl2 = ((x[8]-x[0])*(x[9]-x[1])+(x[10]-x[2])*(x[11]-x[3]))/l2
	th1 = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	th2 = np.arctan((x[10]-x[2])/(x[8]-x[0]))
	f1 = k1 * (l0 -l1) - c1 * dl1 - F_s
	f2 = k1 * (l0 -l2) - c1 * dl2 - F_s 
	#print f
	print F_s
	#f = k1 * (l0 -l) - c1 * dl
	#term1 = f*np.cos(th) / m1
	term1 = (f1*np.cos(th1) - f2*cos(th2))/ m1
	#term2 = f*np.sin(th) / m1 - g
	term2 = (f1*np.sin(th1) - f2*sin(th2)) / m1 - g
	term3 = 0
	term4 = 0 
	term5 = f2*cos(th2) / m3 
	term6 = f2*sin(th2) / m3 -g 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4,x[9],term5,x[11],term6])


def func3(x, l1, l2):
	dl1 = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l1
	dl2 = ((x[8]-x[0])*(x[9]-x[1])+(x[10]-x[2])*(x[11]-x[3]))/l2
	th1 = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	th2 = np.arctan((x[10]-x[2])/(x[8]-x[0]))
	f1 = k1 * (l0 -l1) - c1 * dl1 - F_s
	f2 = k1 * (l0 -l2) - c1 * dl2 - F_s 
	#print f
	print F_s 
	#f = k1 * (l0 -l) - c1 * dl 
	#term1 = f*np.cos(th) / m1
	term1 = (f1*np.cos(th1) - f2*cos(th2))/ m1
	#term2 = f*np.sin(th) / m1 - g
	term2 = (f1*np.sin(th1) - f2*sin(th2)) / m1 - g
	term3 = -f1*np.cos(th1) / m2 
	term4 = -f1*np.sin(th1) / m2 -g 
	term5 = f2*cos(th2) / m3 
	term6 = f2*sin(th2) / m3 -g 
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4,x[9],term5,x[11],term6])


def RK(x, f, l1, l2):
	k1 = f(x, l1, l2)
	k2 = f(x+0.5*h*k1, l1, l2)
	k3 = f(x+0.5*h*k2, l1, l2)
	k4 = f(x+h*k3, l1, l2)
	x_ = x + (h/6)*(k1+2*k2+2*k3+k4)
	#print(x_)
	return x_
	
def Cal_Mtlx(X0, t_s, t_f, l):
	XX = np.empty((0,l), float)
	XX = np.append(XX, np.array([[X0[0],X0[1],X0[2],X0[3],\
	                              X0[4],X0[5],X0[6],X0[7],X0[8],X0[9],X0[10],X0[11]]]), axis=0)
	t = t_s
	n = 0
	while(t<t_f):
			l1 = np.sqrt((XX[n,0]-XX[n,4])**2 + (XX[n,2]-XX[n,6])**2)
			l2 = np.sqrt((XX[n,8]-XX[n,0])**2 + (XX[n,10]-XX[n,2])**2)
			if l1 - del_AE < del_s:
				Step = RK(XX[n], func1, l1, l2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5],Step[6],Step[7],Step[8],Step[9],Step[10],Step[11]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func1")
				t = t+h
				n = n+1
			elif l1 - del_AE >= del_s and l1 -del_AE < del_ME:
				Step = RK(XX[n], func2, l1, l2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
				               Step[4],Step[5],Step[6],Step[7],Step[8],Step[9],Step[10],Step[11]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func2")
				t = t+h
				n = n+1
			else:	
				break

	while(t<t_f):
			l1 = np.sqrt((XX[n,0]-XX[n,4])**2 + (XX[n,2]-XX[n,6])**2)
			l2 = np.sqrt((XX[n,8]-XX[n,0])**2 + (XX[n,10]-XX[n,2])**2)
			if XX[n,2] >= 0 and XX[n,6] >= 0:
				Step = RK(XX[n], func3, l1, l2)
				S = np.array([[Step[0],Step[1],Step[2],Step[3],\
			               Step[4],Step[5],Step[6],Step[7],Step[8],Step[9],Step[10],Step[11]]])
				XX = np.append(XX, S, axis=0)
				print(n)
				print("Using func3")
				t = t+h
				n = n+1
			else:
				break

	return np.matrix(XX)

def main():
	t_s = 0.0
	t_f = 3.0 
	th0 = np.pi/4
	X0 = [del_AE*np.cos(th0),0,del_AE*np.sin(th0),0,0,0,0,0,\
	      del_AE*np.cos(th0)+l0*cos(th0),0,del_AE*sin(th0)+l0*sin(th0), 0]
	print(X0)
	XX = Cal_Mtlx(X0, t_s, t_f, 12)
	print(XX)

	#print(Time)
	#T = np.arange(0, t_f+2*h, h)
	print("Size of XX")
	row, col = XX.shape
	print(row, col)
	T = np.arange(0, row*h, h)
	print("Length of T={0}".format(len(T)))
	print("del_AE={0}\n del_s={1}\n, del_ME={2}\n, F_s={3}\n, A={4}" .format(del_AE,del_s,del_ME, F_s, A))

	plt.figure()
	plt.title('Hopping Distance X,Z w.r.t Time')
	plt.xlabel('Time[s]')
	plt.ylabel('Hopping Height[m]')
	plt.plot(T,XX[:,0], label="x1")
	plt.plot(T,XX[:,2], label="z1")
	plt.plot(T,XX[:,4], label="x2")
	plt.plot(T,XX[:,6], label="z2")
	plt.plot(T,XX[:,8], label="x3")
	plt.plot(T,XX[:,10], label="z3")
	plt.xlim([0,T[row-1]])
	plt.legend()
	#plt.savefig("Rigid_NoSlip_T-XZ.png")
	
	plt.figure()
	plt.plot(XX[:,0],XX[:,2], label="Body")
	plt.title('Hopping Trajectory Z w.r.t X')
	plt.xlabel('Distance X[m]')
	plt.ylabel('Distance Z[m]')
	plt.xlim([0,XX[row-1,0]])
	plt.legend()
	plt.savefig("Rigid_NoSlip_XZ.png")
	#plt.show()

	fig = plt.figure()
	ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.05, XX[row-1,0]+0.1),\
	                                              ylim=(-0.05, XX[row-1,0]+0.1))
	ax.grid()
	
	line, = ax.plot([], [], '-o', lw=3)
	time_template = 'time = %.1fs'
	time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
	
	
	def init():
	    line.set_data([], [])
	    time_text.set_text('')
	    return line, time_text
	
	
	def animate(i):
	    thisx = [XX[i*20,0], XX[i*20,4], XX[i*20,8]]
	    thisy = [XX[i*20,2], XX[i*20,6], XX[i*20,10]]
	
	    line.set_data(thisx, thisy)
	    time_text.set_text(time_template % (i*20*h))
	    return line, time_text
	
	ani = animation.FuncAnimation(fig, animate, np.arange(1, len(T)/20),
	                              interval=1, blit=False, init_func=init)
	
	#ani.save("Rigid_NoSlip.gif", writer = 'imagemagick')
	plt.show()
if __name__ == '__main__':
		main()
