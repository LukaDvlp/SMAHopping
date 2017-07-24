import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Simu0619

m1 = 0.65  #Mass of Body [Kg]
m2 = 0.05 #Mass of Pad [Kg]
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
mu = 0.9
mu_d = 0.8
del_AE = Simu0619.delta_AE / 1000
del_s = Simu0619.delta_s / 1000
del_ME = Simu0619.delta_ME / 1000
F_s = Simu0619.F_s 
A = 1/Simu0619.A * 1000 


def func1(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - A*(l-del_AE)
	#f = k1 * (l0 -l) - 0 * dl
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	if f*np.cos(th) < mu*np.sin(th):
		term3 = 0
		term4 = (-f*np.sin(th) - m2*g - k2*x[6] -c2*x[7]) / m2
	else:
		term3 = mu_d*f*np.sin(th) - f*np.cos(th)
		term4 = (-f*np.sin(th) - m2*g - k2*x[6] -c2*x[7]) / m2
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])

def func2(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - F_s 
	#f = k1 * (l0 -l) - 0 * dl
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	if f*np.cos(th) < mu*np.sin(th):
		term3 = 0
		term4 = (-f*np.sin(th) - m2*g - k2*x[6] -c2*x[7]) / m2
	else:
		term3 = mu_d*f*np.sin(th) - f*np.cos(th)
		term4 = (-f*np.sin(th) - m2*g - k2*x[6] -c2*x[7]) / m2
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])

def func3(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - A*(l-del_AE)
	#f = k1 * (l0 -l) - 0 * dl
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	if f*np.cos(th) < mu*np.sin(th):
		term3 = 0
		term4 = (-f*np.sin(th) - m2*g - k2*x[6]) / m2
	else:
		term3 = mu_d*f*np.sin(th) - f*np.cos(th)
		term4 = (-f*np.sin(th) - m2*g - k2*x[6]) / m2
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])
	
def func4(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
	f = k1 * (l0 -l) - c1 * dl - F_s 
	#f = k1 * (l0 -l) - 0 * dl
	term1 = f*np.cos(th) / m1
	term2 = f*np.sin(th) / m1 - g
	if f*np.cos(th) < mu*np.sin(th):
		term3 = 0
		term4 = (-f*np.sin(th) - m2*g - k2*x[6]) / m2
	else:
		term3 = mu_d*f*np.sin(th) - f*np.cos(th)
		term4 = (-f*np.sin(th) - m2*g - k2*x[6]) / m2
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])

def func5(x, l):
	#l = np.sqrt((x[0]-x[2])**2 + (x[1]-x[3])**2)
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	f = k1 * (l0 -l) - c1 * dl - A*(l-del_AE)
	#f = k1 * (l0 -l) - 0 * dl
	if x[0]-x[4]>=0 and x[2]-x[6]>=0:
		th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
		term1 = f*np.cos(th) / m1
		term2 = f*np.sin(th) / m1 - g
		term3 = -f*np.cos(th) /m2  
		term4 = (-f*np.sin(th) - m2*g) / m2
	elif x[0]-x[4]>=0 and x[2]-x[6]<0:
		th = np.arctan((x[6]-x[2])/(x[0]-x[4]))
		term1 = f*np.cos(th) / m1
		term2 = -f*np.sin(th) / m1 - g
		term3 = -f*np.cos(th) /m2  
		term4 = (f*np.sin(th) - m2*g) / m2
	elif x[0]-x[4]<0 and x[2]-x[6]<0:
		th = np.arctan((x[6]-x[2])/(x[4]-x[0]))
		term1 = -f*np.cos(th) / m1
		term2 = -f*np.sin(th) / m1 - g
		term3 = f*np.cos(th) /m2  
		term4 = (f*np.sin(th) - m2*g) / m2
	else:
		th = np.arctan((x[2]-x[6])/(x[4]-x[0]))
		term1 = -f*np.cos(th) / m1
		term2 = f*np.sin(th) / m1 - g
		term3 = f*np.cos(th) /m2  
		term4 = (-f*np.sin(th) - m2*g) / m2
	return np.array([x[1], term1, x[3], term2, x[5], term3, x[7],term4])

def func6(x, l):
	dl = ((x[0]-x[4])*(x[1]-x[5])+(x[2]-x[6])*(x[3]-x[7]))/l
	f = k1 * (l0 -l) - c1 * dl - F_s 
	#f = k1 * (l0 -l) - 0 * dl
	if x[0]-x[4]>=0 and x[2]-x[6]>=0:
		th = np.arctan((x[2]-x[6])/(x[0]-x[4]))
		term1 = f*np.cos(th) / m1
		term2 = f*np.sin(th) / m1 - g
		term3 = -f*np.cos(th) /m2  
		term4 = (-f*np.sin(th) - m2*g) / m2
	elif x[0]-x[4]>=0 and x[2]-x[6]<0:
		th = np.arctan((x[6]-x[2])/(x[0]-x[4]))
		term1 = f*np.cos(th) / m1
		term2 = -f*np.sin(th) / m1 - g
		term3 = -f*np.cos(th) /m2  
		term4 = (f*np.sin(th) - m2*g) / m2
	elif x[0]-x[4]<0 and x[2]-x[6]<0:
		th = np.arctan((x[6]-x[2])/(x[4]-x[0]))
		term1 = -f*np.cos(th) / m1
		term2 = -f*np.sin(th) / m1 - g
		term3 = f*np.cos(th) /m2  
		term4 = (f*np.sin(th) - m2*g) / m2
	else:
		th = np.arctan((x[2]-x[6])/(x[4]-x[0]))
		term1 = -f*np.cos(th) / m1
		term2 = f*np.sin(th) / m1 - g
		term3 = f*np.cos(th) /m2  
		term4 = (-f*np.sin(th) - m2*g) / m2
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
	t = t_s
	n = 0
	while(t<t_f):
			l = np.sqrt(np.absolute(XX[n,0]-XX[n,4])**2 + np.absolute(XX[n,2]-XX[n,6])**2)
			if l - del_AE < del_s:
				if XX[n,6]< 0 and XX[n,7]<0:
					Step = RK(XX[n], func1, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func1")
					t = t+h
					n = n+1
				elif XX[n,6]<0 and XX[n,7]>=0:	
					Step = RK(XX[n], func3, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func2")
					t = t+h
					n = n+1
				elif XX[n,2]<0 or XX[n,6]<0:
					break
				else:
					Step = RK(XX[n], func5, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func3")
					t = t+h
					n = n+1
			elif l - del_AE >= del_s and l -del_AE < del_ME:
				if XX[n,6]< 0 and XX[n,7]<0:
					Step = RK(XX[n], func2, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func2")
					t = t+h
					n = n+1
				elif XX[n,6]<0 and XX[n,7]>=0:	
					Step = RK(XX[n], func4, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func4")
					t = t+h
					n = n+1
				elif XX[n,2]<0 or XX[n,6]<0:
					break
				else:
					Step = RK(XX[n], func6, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func6")
					t = t+h
					n = n+1
			else:
				if XX[n,6]< 0 and XX[n,7]<0:
					Step = RK(XX[n], func2, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func1")
					#print"term1 = {0}".format(term1)
					#print"term2 = {0}".format(term2)
					t = t+h
					n = n+1
				elif XX[n,6]<0 and XX[n,7]>=0:	
					Step = RK(XX[n], func4, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func2")
					t = t+h
					n = n+1
				elif XX[n,2]<0 or XX[n,6]<0:
					break
				else:
					Step = RK(XX[n], func6, l)
					S = np.array([[Step[0],Step[1],Step[2],Step[3],\
					               Step[4],Step[5],Step[6],Step[7]]])
					XX = np.append(XX, S, axis=0)
					print(n)
					print("Using func3")
					t = t+h
					n = n+1
	return np.matrix(XX)

def main():
	t_s = 0.0
	t_f = 3.0 
	th0 = np.pi/4
	X0 = [del_AE*np.cos(th0),0,del_AE*np.sin(th0),0,0,0,0,0]
	print(X0)
	XX = Cal_Mtlx(X0, t_s, t_f, 8)
	print(XX)

	print("Size of XX")
	row, col = XX.shape
	print(row, col)
	T = np.arange(0, row*h, h)
	print("Length of T={0}".format(len(T)))

	plt.figure()
	plt.plot(T,XX[:,0], label="x1")
	plt.plot(T,XX[:,2], label="z1")
	plt.plot(T,XX[:,4], label="x2")
	plt.plot(T,XX[:,6], label="z2")
	plt.xlim([0,T[row-1]])
	plt.legend()
	
	plt.figure()
	plt.plot(XX[:,0],XX[:,2], label="m1")
	plt.title('1-Dimensional Hopping of Spring Only Rover')
	plt.xlabel('Time[s]')
	plt.ylabel('Hopping Height[m]')
	plt.xlim([0,XX[row-1,0]])
	plt.legend()
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
	
	TM_thin = 20	
	def animate(i):
	    thisx = [XX[i*TM_thin,0], XX[i*TM_thin,4]]
	    thisy = [XX[i*TM_thin,2], XX[i*TM_thin,6]]
	
	    line.set_data(thisx, thisy)
	    time_text.set_text(time_template % (i*TM_thin*h))
	    return line, time_text
	
	ani = animation.FuncAnimation(fig, animate, np.arange(1, len(T)/TM_thin),
	                              interval=1, blit=False, init_func=init)
	
	#ani.save("double_pendulum.gif", writer = 'imagemagick')
	plt.show()
if __name__ == '__main__':
		main()
