import numpy as np

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
if __name__ == '__main__':
		main()
