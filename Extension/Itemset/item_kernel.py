import numpy as np

datasets=["dna", "mushrooms", "a9a"]


def readX(path):
	X=[]
	with open(path) as f:
		while True:
			line = f.readline().strip()
			if not line:
				break
			X.append(set(line.split(" ")[1:]))
	return X


def jaccard(X):
	n=len(X)
	K=np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			K[i,j]=float(len(X[i]&X[j])) / len(X[i]|X[j])
	return K

for data in datasets:
	X=readX("../ItemData/"+data)
	K=jaccard(X)
	np.savetxt(data, K, fmt="%.8e")
