#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
import numpy as np
import math as mt
from scipy.optimize import curve_fit
def func(x,a,k,b):
	return a * (np.exp(-k*x)) + b
def main():
	time = open(sys.argv[1],'r')
	maxima = open(sys.argv[2],'r')
	x = []
	y = []
	for i in maxima.readlines():
		y.append(float(i))
	for n in time.readlines():
		x.append(float(n))
#	print len(x), len(y)
#	for j in range(len(x)-1):
#		if y[j] <= 0.5:
#			if y[j] <= y[j+1]:
#				if y[j+1] <= y[j+2]:
#					startpoint = j
#					break
#	lnx = np.zeros(len(x[startpoint:]))
#	lny = np.zeros(len(y[startpoint:]))
#	for i in range(len(lnx)):
#		lnx[i] = 1/(x[i+startpoint])
#		lny[i] = mt.log(y[i+startpoint])
#	xx = np.zeros(len(x[:-4]))
#	yy = np.zeros(len(y[:-4]))
#	xx [:] = x[:-4]
#	yy [:] = y[:-4]
	p0 = (1.,0.5e-2,1)
	param, pcov = curve_fit(func, x, y, p0, method='lm')
	a,k,b = param
	print 'Best Fit Equation:\nY = Ae^(-kx)+B'
	print 'A = '+str(a)
	print 'B = '+str(b)
	print 'k = '+str(k)
#	T1 = mt.log(b/a)/-k
#	print 'Spin-Lattice relaxation time = '+str(T1)
	x2 = np.linspace(0,100,250)
	y2 = func(x2,a,k,b)
#	raw = plt.figure(1)
#	plt.plot(x,y,'rx')
	plt.plot(x,y,'rx',x2,y2,'-b')
#	plt.ylim([-0.5,1.5])
#	plt.title('$M_z$ vs. delay time for mineral water')
#	plt.xlabel('Delay Time (ms)')
#	plt.ylabel('Magnetization (V)')
#	plt.text(75,1.0, \
#		'Best fit equation:\nY = $-A * e^{-k*x} + B$\nA = '+str(a)+ \
#		'\nB = '+str(b)+'\nk = '+str(k))
	plt.grid()
	plt.show()
	raw_input
main()
