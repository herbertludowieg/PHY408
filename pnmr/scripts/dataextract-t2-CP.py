#!/usr/bin/env python
#
# this script is used when reading one file with many different maxima
# it will only get the maximum value when the values of the data are above
# a certain threshold
# be mindful of the limits that are set
#
import sys
import numpy as np
def main():
	data_raw = []
	scaling = []
	fn = open(sys.argv[1],'r')
	for i in fn.readlines():
		if i[0:8] == 'Vertical':
			all_scale = i.split(',,')
			a,scaling_str = all_scale[0].split(',')
			try:
				scaling.append(float(scaling_str))
			except:
				pass
		if i[0:3] == ',,,':
			d = i[3:-2].split(',')
			x = d[0]
			y = d[1]
			x = x.strip()
			y = y.strip()
			data_raw.append((float(x)*1000,float(y)))
	maxvalues = []
	for i in range(1,len(data_raw)-1):
		if data_raw[i][1] >= data_raw[i-1][1] and data_raw[i][1] >= 0.32:
			if data_raw[i][1] > data_raw[i+1][1]:
#				print data_raw[i][0],data_raw[i][1]*scaling[0]
				maxvalues.append((data_raw[i][0],data_raw[i][1]))
	toremove = []
	for i in range(len(maxvalues)-1):
#		print len(maxvalues),i
		try:
			if abs(maxvalues[i][0] - maxvalues[i+1][0]) < 0.15:
				maxvalues.remove((maxvalues[i+1][0],maxvalues[i+1][1]))
		except:
			break
	j=0
	ave = 0
	if data_raw[0][1] <= 0.3:
		j = 50
	while (data_raw[j][1] <= 0.31):
		ave = ave + data_raw[j][1]
		j+=1
	ave = ave/j
	maxvalues.remove(maxvalues[0])
	for i in maxvalues:
		print i[0],(i[1]-ave)*scaling[0]
main()
