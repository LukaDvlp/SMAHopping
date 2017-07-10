import numpy as np

def func(x):
	return np.array([[2*x], [x[0], x[1]]])

def main():
	x = [1,1]
	y = [6,7,8,9,10]
	theta = np.pi/4 
	#ZZ = np.array([])
	#ZZ = np.append(ZZ, func(x)[1], axis=0)
	X = np.matrix([x]).T
	#X = X.T
	ZZ = np.matrix([[np.sin(theta), np.cos(theta)], [-np.cos(theta), np.sin(theta)]])
	print(ZZ)
	Z1 = ZZ*np.matrix([[1],[2]]) 
	print(Z1)
	xx = Z1[0,0]
	print(xx)
	yy = Z1[1,0]
	print(yy)
	XX = np.matrix([[1,2],[3,4]])
	ZZ1_ = np.linalg.inv(XX)
	print(ZZ1_)

	M1 = 1
	g = 1
	f = np.matrix([[0], [M1*g]])
	print(f)
	
	l = np.array([1,2])
	g = np.matrix([[l[0]], [l[1]]]).T
	print(g)
	
	cc = ZZ1_*f
	print(cc)

if __name__ == '__main__':
		main()
