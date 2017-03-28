#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from math import floor,log10
import sys
def func(x,a,b):
	return a*x+b
def sig_fig(x,y):
	return round(x,-int(floor(log10(abs(x))))), \
			round(y,-int(floor(log10(abs(x)))))
def main():
	if len(sys.argv) != 2:
		print 'ERROR cannot parse commands\n'+ \
			'USAGE: '+sys.argv[0]+' [B field file]'
		sys.exit()
	fn = open(sys.argv[1],'r')
	x = []
	y = []
	for i in fn.readlines():
		d = i.split()
		x.append(float(d[0]))
		y.append(float(d[1]))
	p0 = (1,1)
	param,pcov = curve_fit(func,x,y,p0)
	a,b = param

	sigma_a,sigma_b = np.sqrt(np.diag(pcov))

	round_a = sig_fig(sigma_a,a)
	round_b = sig_fig(sigma_b,b)
	print 'Parameters:\n'+ \
		'a = '+str(round_a[1])+' +/- '+str(round_a[0])+ \
		'\nb = '+str(round_b[1])+' +/- '+str(round_b[0])
	x2 = np.linspace(-0.67,0.67,1000)
	y2 = func(x2,a,b)
	plt.plot(x,y,'rx',x2,y2,'-b')
	plt.show()
main()
