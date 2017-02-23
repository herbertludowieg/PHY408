#!/usr/bin/env python
#
# program takes the file that is inputted either from the command line or by
# the shell script and selects data according to the maximum value that is set
# for the data
# will then print the result on the command line to view or save to a file
# will also calculate the average offset within the data and proceed to
# subtract the average offset from the data values that are of interest
# scales the data as it prints
# will only work when there is a negative number contained in the data so that 
# it can select and splice the correct data
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
			#print d
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
		if data_raw[i] <= 0.0:
			minpoint = i
		if minpoint != 0 and data_raw[i] <= 0.4:
			maxpoint = i
	j=0
	while (data_raw[j] > 0.20):
		ave = ave + data_raw[j]
		j += 1
	ave = ave/j
	print (max(data_raw[minpoint:maxpoint])-ave)*scaling[0]
			
main()
