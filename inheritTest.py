import numpy as np

def get_max(Result):
	index = Result[:,6].argmax()
	print Result[index,:]

def getResult(Result): #modify from here
	get_max(Result)
	print index

#arr = [[1,2,3,4,5,6,7],[10,20,30,40,50,60,70],[100,200,300,400,500,600,700]]
#Arr = np.array(arr)

Result = np.empty((0),float)
k=0
while (k<10):
		Result = np.append(Result, np.array([k]), axis=0)
		k += 1

index = Result.argmax()
print index
add = np.array([100])
Result = np.append(Result, add, axis=0)
print Result
