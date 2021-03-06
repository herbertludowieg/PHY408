#!/usr/bin/env python
#
# finds the exponential fit function for the inputted data
# must input two seperate data files
# calculates a rising exponential function and then finds the x-intercept
# for the purpose that this program is used for this finds the T1 relaxation
# time for a PNMR experiment
#
import matplotlib.pyplot as plt
import sys
import numpy as np
import math as mt
from scipy.optimize import curve_fit
def sig_fig(x,y):
	return round(x,-int(mt.floor(mt.log10(abs(x))))), \
				round(y,-int(mt.floor(mt.log10(abs(x)))))
def func(x,a,k,b):
	return a * (-np.exp(-k*x)) + b
def main():
	time = open(sys.argv[1],'r')
	maxima = open(sys.argv[2],'r')
	x = []
	y = []
	for i in maxima.readlines():
		y.append(float(i))
	for n in time.readlines():
		x.append(float(n))
	for j in range(len(x)-1):
#		print y[j-1],y[j],y[j+1]
		if y[j] <= y[j+1] and y[j] < 0.5:
			if y[j] < y[j-1]:
				for l in range(len(y[:j+1])):
					y[l] = -1*y[l]
				break
#	print y
#	xx = np.zeros(len(x[startpoint:]))
#	yy = np.zeros(len(y[startpoint:]))
#	xx [:] = x[startpoint:]
#	yy [:] = y[startpoint:]
	p0 = (1.,0.5e-5,0.947)
	param, pcov = curve_fit(func, x, y, p0, method='lm')
	a,k,b = param
	sigma_a,sigma_k,sigma_b = np.sqrt(np.diag(pcov))
	round_a = sig_fig(sigma_a,a)
	round_b = sig_fig(sigma_b,b)
	round_k = sig_fig(sigma_k,k)
	print 'Best Fit Equation:\nY = -Ae^(-kx)+B'
	print 'A = '+str(round_a[1])+' +/- '+str(round_a[0])
	print 'B = '+str(round_b[1])+' +/- '+str(round_b[0])
	print 'k = '+str(round_k[1])+' +/- '+str(round_k[0])
	T1 = mt.log(a/b)/k
	sigma_T1 = np.sqrt((sigma_a/(k*a))**2+(sigma_b/(k*b))**2+ \
				(((mt.log(b)-mt.log(a))/(k*k))*sigma_k)**2)
	round_T1 = sig_fig(sigma_T1,T1)
	print 'Spin-Lattice relaxation time = '+str(round_T1[1])+' +/- '+ \
							str(round_T1[0])
	x2 = np.linspace(0,x[-1],10000)
	y2 = func(x2,a,k,b)
#	raw = plt.figure(1)
	plt.plot(x,y,'rx',x2,y2,'-b')
#	plt.ylim([-0.5,1.5])
	plt.xlim([0,25])
	plt.title(r'$M_z$ vs. $\tau$')
	plt.xlabel(r'$\tau$ (ms)')
	plt.ylabel('Magnetization (V)')
	plt.text(1000,-0.5, \
		'Best fit equation:\nY = $-A * e^{-k*x} + B$\nA = '+ \
		str(round_a[1])+' +/- '+str(round_a[0])+'\nB = '+ \
		str(round_b[1])+' +/- '+str(round_b[0])+'\nk = '+ \
		str(round_k[1])+' +/- '+str(round_k[0]))
	plt.grid()
	plt.show()
#	raw_input
main()
