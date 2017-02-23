#!/usr/bin/env python
#
# calculates the exponential decay function for a data set
# will be able to be used with only one data set rather than two as was the
# case from before
# will then find the y-intercept and calculate the T2 time constant
#
import matplotlib.pyplot as plt
import sys
import numpy as np
import math as mt
from scipy.optimize import curve_fit
def func(x,a,k,b):
	return a * (np.exp(-k*x)) + b
def main():
	fn = open(sys.argv[1],'r')
	x = []
	y = []
	for i in fn.readlines():
		d = i.split()
		x.append(float(d[0]))
		y.append(float(d[1]))
	print len(x), len(y)
	p0 = (1,0.5e-2,1)
	param, pcov = curve_fit(func, x[1:], y[1:], p0)
	a, k, b = param
	x2 = np.linspace(0,0.03,250)
	y2 = func(x2,a,k,b)
	plt.plot(x,y,'rx',x2,y2,'-b')
	plt.grid()
	plt.show()

main()
