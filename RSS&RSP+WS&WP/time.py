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
		value.append(np.sum(result))
	return np.array(value)

argvs = sys.argv

if len(argvs)!=3:
	print("python time.py dataset test")
	exit()

dataset = argvs[1]
test = int(argvs[2])

print(dataset)
print("----- Time -----")

all = getVal(r'"time = " | sed "s/..*time.*//" | grep "time = " | sed "s/time = //" | tr "\n" ","', test)
traverse = getVal(r'":time = " | sed "s/.*:time = //" | tr "\n" ","', test)

print('traverse mean = {:.1f}, std = {:.1f}'.format(np.mean(traverse), np.std(traverse)))
print('solve mean = {:.1f}, std = {:.1f}'.format(np.mean(all-traverse), np.std(all-traverse)))
print('all mean = {:.1f}, std = {:.1f}'.format(np.mean(all), np.std(all)))
