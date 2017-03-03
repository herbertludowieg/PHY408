#!/usr/bin/env python
#
# finds the exponential decay function for the inputted data
# this program will then calculate the y-intercept to find the T2 relaxation
# time for a sample
#
import matplotlib.pyplot as plt
import sys
import numpy as np
import math as mt
from scipy.optimize import curve_fit
def sig_fig(x,y):
	return round(x,-int(mt.floor(mt.log10(abs(x))))), \
				round(y,-int(mt.floor(mt.log10(abs(x)))))

def func(x,a,k):
	return a * (np.exp(-k*x))
def main():
	if len(sys.argv) == 3:
		time = open(sys.argv[1],'r')
		maxima = open(sys.argv[2],'r')
		x = []
		y = []
		for i in maxima.readlines():
			y.append(float(i))
		for n in time.readlines():
			x.append(float(n))
	elif len(sys.argv) == 2:
		fn = open(sys.argv[1],'r')
		x = []
		y = []
		for i in fn.readlines():
			d = i.split()
			x.append(float(d[0]))
			y.append(float(d[1]))
	else:
		print 'ERROR cannot parse arguments\n'+ \
			'USAGE: (Two single column files)'+sys.argv[0]+ \
						' [file1] [file2]\n'+ \
			'USAGE: (One double column file) '+sys.argv[0]+ \
							' [file1]'
		sys.exit()
#	p0 = (1.,0.5e-2,1)
	p0 = (1.,0.5e-2)
	param, pcov = curve_fit(func, x, y, p0, method='lm')
#	a,k,b = param
	a,k = param
#	sigma_a,sigma_k,sigma_b = np.sqrt(np.diag(pcov))
	sigma_a,sigma_k = np.sqrt(np.diag(pcov))

	round_a = sig_fig(sigma_a,a)
#	round_b = sig_fig(sigma_b,b)
	round_k = sig_fig(sigma_k,k)
	print 'Best Fit Equation:\nY = Ae^(-kx)'
	print 'A = '+str(round_a[1])+' +/- '+str(round_a[0])
#	print 'B = '+str(round_b[1])+' +/- '+str(round_b[0])
	print 'k = '+str(round_k[1])+' +/- '+str(round_k[0])
	t2 = 1/k
	sigma_t2 = sigma_k/(k*k)
	round_t2 = sig_fig(sigma_t2,t2)
#	T2_all = np.zeros(len(x))
#	sigma_T2_all = np.zeros(len(x))
#	sigma_T2 = 0
#	T2 = 0
#	for i in range(len(x)):
#		T2_all[i] = x[i]*(mt.log(a+b)-mt.log(y[i]))
#		T2 += T2_all[i]
#		sigma_T2_all[i] = (x[i]/(a+b)) * np.sqrt(sigma_a*sigma_a + \
#								sigma_b*sigma_b)
#		sigma_T2 += sigma_T2_all[i]
#	T2 = T2/len(x)
#	sigma_T2 = sigma_T2/len(x)
#	round_T2 = sig_fig(sigma_T2,T2)
	print 'Spin-Spin relaxation time = '+str(round_t2[1])+' +/- '+ \
							str(round_t2[0])
	x2 = np.linspace(0,x[-1],250)
	y2 = func(x2,a,k)
#	raw = plt.figure(1)
#	plt.plot(x,y,'rx')
	plt.plot(x,y,'rx',x2,y2,'-b')
#	plt.ylim([-0,1.2])
	plt.title(r'$M_z$ vs. $2\tau$')
	plt.xlabel(r'$2\tau$ (ms)')
	plt.ylabel('Magnetization (V)')
	plt.text(150,0.8, \
		'Best fit equation:\nY = $-A * e^{-k*x}$\nA = '+ \
		str(round_a[1])+' +/- '+str(round_a[0])+'\nk = '+ \
#		str(round_b[1])+' +/- '+str(round_b[0])+'\nk = '+ \
		str(round_k[1])+' +/- '+str(round_k[0]))
	plt.grid()
	plt.show()
main()
