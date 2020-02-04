import numpy as np
from strkernel.mismatch_kernel import MismatchKernel

datasets=["promoters", "splice"]

for data in datasets:
	yX=np.loadtxt("../SequenceData/"+data,dtype=int)
	X=yX[:,1:]-1
	for k in [3,4,5]:
		mismatch_kernel = MismatchKernel(l=np.unique(X).shape[0], k=k, m=1).get_kernel(X).kernel
		np.savetxt(data+str(k), mismatch_kernel,fmt="%.8e")
