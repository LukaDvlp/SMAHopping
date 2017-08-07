import numpy as np
from numpy import sin, cos, pi

def test():
	x=1
	k=2
	return [x,k]

def main():
	a = 9
	b = np.sqrt(a)
	print b
	XX = np.matrix([[1,2,3],[4,5,6],[7,8,9]])
	print XX[1,1]

	M = test()[1]
	N = test()[0]
	print M,N

	print pi
	print np.pi
if __name__ == '__main__':
		main()
