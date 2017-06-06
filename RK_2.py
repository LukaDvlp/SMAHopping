import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m = 1
g = 1.622 
H = 10 
v0 = 10
h = 0.01
t_s = 0.0
t_f = 10.0
n =0
XX = np.zeros([2,t_f/h])
A = np.matrix([[0,1],[0,0]])
b = np.matrix([[0],[-g]])

def func(x, t):
	return[x[1], -g]

def main():
	n = 0
	t = t_s
	XX[:,0] = [H, v0]
	#print XX
	while t<t_f:
		k1 = A*XX[:,n]+b	
		k2 = A*(XX[:,n]+0.5*h*k1) +b
		k3 = A*(XX[:,n]+0.5*h*k2) +b
		k4 = A*(XX[:,n]+h*k3) +b
		XX[:,n+1] = XX[:,n] + (h/6)*(k1+k2+k3+k4)
		n = n+1
		t = t+h 
	#print(x)
	#fig = plt.figure()
	#plt.plot(t,x[:,0], label="Position")
	#plt.plot(t,x[:,1], label="Velocity")
	#ax = fig.gca(projection='3d')
	#ax.plot(v[:, 0], v[:, 1], v[:, 2])
	#plt.legend()
	#plt.show()

if __name__ == '__main__':
		main()
