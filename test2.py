import numpy as np

def func(x):
	return np.array([[2*x], [x[0], x[1]]])

def main():
	x = [1,2,3,4,5]
	ZZ = np.array([])
	y, ZZ = func(x)[0], np.a
	ZZ = np.append(ZZ, func(x)[1], axis=0)
	print(y)
	print(ZZ)


if __name__ == '__main__':
		main()
