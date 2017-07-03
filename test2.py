import numpy as np

def func(x):
	return np.array([[2*x], [x[0], x[1]]])

def main():
	x = [1,2,3,4,5]
	ZZ = np.array([])
	ZZ = np.append(ZZ, func(x)[1], axis=0)
	print(ZZ)
	Z1 = 3*ZZ
	print(Z1)


if __name__ == '__main__':
		main()
