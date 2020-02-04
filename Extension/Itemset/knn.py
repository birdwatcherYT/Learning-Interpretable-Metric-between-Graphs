import numpy as np
import subprocess
import sys

def getVal(pattern, test):
	value = []
	for i in range(test):
		result = subprocess.getoutput('grep ' + str(i+1) + '/' + dataset + ' -e ' + pattern)
		if len(result)==0:
			return -1, -1
		result = np.fromstring(result, sep=',')
		result = np.reshape(result, (-1,2))
		best = np.nanargmax(result[:,0])
		value.append(result[best,1])
	return np.mean(value), np.std(value)

argvs = sys.argv

if len(argvs)!=3:
	print("python knn.py dataset test")
	exit()

dataset = argvs[1]
test = int(argvs[2])

print(dataset)

print("----- Fmacro -----")

mean, std = getVal( '"-nn: KernelValid-Fmacro " | sed "s/.* KernelValid-Fmacro = //" | sed "s/, KernelTest-Fmacro = /,/" | tr "\n" ","', test)
print('kernel mean = {:.5f}, std = {:.5f}'.format(mean, std))

mean, std = getVal( '"-nn: Valid-Fmacro " | sed "s/.* Valid-Fmacro = //" | sed "s/, Test-Fmacro = /,/" | tr "\n" ","', test)
print('proposed diag mean = {:.5f}, std = {:.5f}'.format(mean, std))

mean, std = getVal( '"-nn: FullValid-Fmacro " | sed "s/.* FullValid-Fmacro = //" | sed "s/, FullTest-Fmacro = /,/" | tr "\n" ","', test)
print('proposed full mean = {:.5f}, std = {:.5f}'.format(mean, std))

print("----- Fmicro -----")

mean, std = getVal( '"-nn: KernelValid-Fmicro " | sed "s/.* KernelValid-Fmicro = //" | sed "s/, KernelTest-Fmicro = /,/" | tr "\n" ","', test)
print('kernel mean = {:.5f}, std = {:.5f}'.format(mean, std))

mean, std = getVal( '"-nn: Valid-Fmicro " | sed "s/.* Valid-Fmicro = //" | sed "s/, Test-Fmicro = /,/" | tr "\n" ","', test)
print('proposed diag mean = {:.5f}, std = {:.5f}'.format(mean, std))

mean, std = getVal( '"-nn: FullValid-Fmicro " | sed "s/.* FullValid-Fmicro = //" | sed "s/, FullTest-Fmicro = /,/" | tr "\n" ","', test)
print('proposed full mean = {:.5f}, std = {:.5f}'.format(mean, std))

