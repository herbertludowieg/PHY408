#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import numpy as np
import plotly.tools as tls
def main():
	if len(sys.argv) == 2:
		fn = open(sys.argv[1],'r')
		magnetdatatran = []
		for i in fn.readlines():
			d = i.split('	')
			d[-1] = d[-1].strip('\n')
			magnetdatatran.append((d))
		magnetdata = np.zeros((11,13))
		for i in range(0,11):
			for j in range(0,13):
				magnetdata[-i-1][j] = float(magnetdatatran[j][i])
	else:
		print 'ERROR make sure you input data'
		sys.exit()
	
	
	plt.imshow(magnetdata,extent=[-3,3,10,20],interpolation='nearest')
	plt.colorbar()
	plt.xlabel('X-axis movement (-3,3)')
	plt.ylabel('Y-axis movement (10,20)')

	plt.show()
	
main()
