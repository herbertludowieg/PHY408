#!/usr/bin/env python
#
# this is to be used much like with the *-t1.py script except that it will work
# for the spin echo rather than the magnetization values
# will only use the values that are above a certain value that must be set or
# it will take all of the inputted data
#
import sys
import numpy as np
def main():
	data_raw = []
	scaling = []
	#print '#######################'
	fn = open(sys.argv[1],'r')
	for i in fn.readlines():
		#print i[0:3]
		if i[0:8] == 'Vertical':
			all_scale = i.split(',,')
			a,scaling_str = all_scale[0].split(',')
			try:
				scaling.append(float(scaling_str))
			except:
				pass
			#print scaling
		if i[0:3] == ',,,':
			d = i[3:-2].split(',')
#			print i
#			print d
			x = d[0]
			y = d[1]
			#y = d.split(',')
			x = x.strip()
			y = y.strip()
#			print x,y
			data_raw.append(float(y))
#	print max(data_raw)*scaling[0]
	minpoint = 0
	maxpoint = 0
	ave = 0
	for i in range(len(data_raw)):
		if data_raw[i] >= 0.33:
			minpoint = i
		if minpoint != 0 and data_raw[i] <= 0.33:
			maxpoint = i
			break 
#	print max(data_raw)#*scaling[0]
	j=0
	while (data_raw[j] < 0.30):
		ave = ave + data_raw[j]
		j += 1
	if j == 0:
		j = 100
		while (abs(data_raw[j]) < 0.30):
			ave = ave + data_raw[j]
			j += 1

	ave = ave/j
	print (max(data_raw))*scaling[0]
			
main()
