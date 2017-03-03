#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit

def sig_fig(x,y):
	return round(x,-int(np.floor(np.log10(abs(x))))), \
				round(y,-int(np.floor(np.log10(abs(x)))))

def func(x,a,b):
	return a*x+b
def main():
	if len(sys.argv) > 2:
		fn1 = open(sys.argv[1],'r')
		fn2 = open(sys.argv[2],'r')
	else:
		print 'Cannot parse arguments\n'+ \
			'Please input the arguments as\n'+ \
			sys.argv[0]+' [T1 Molar Data] [T2 Molar Data]'
		sys.exit()
	t1 = []
	mol = []
	sigma_t1 = []
	for i in fn1.readlines():
		d = i.split()
		t1.append(1/float(d[1]))
		mol.append(float(d[0]))
		sigma_t1.append(float(d[2]))
	t2 = []
	sigma_t2 = []
	for i in fn2.readlines():
		d = i.split()
		t2.append(1/float(d[1]))
		sigma_t2.append(float(d[2]))
	p0 = (1,0)
	param1,pcov1 = curve_fit(func,mol,t1,p0)
	a1,b1 = param1
	param2,pcov2 = curve_fit(func,mol,t2,p0)
	a2,b2 = param2
	mol1 = np.linspace(0,mol[-1],100)
	t11 = func(mol1,a1,b1)
	mol2 = np.linspace(0,mol[-1],100)
	t22 = func(mol2,a2,b2)
	sigma_a1,sigma_b1 = np.sqrt(np.diag(pcov1))
	sigma_a2,sigma_b2 = np.sqrt(np.diag(pcov2))
	round_a1 = sig_fig(sigma_a1,a1)
	round_b1 = sig_fig(sigma_b1,b1)
	round_a2 = sig_fig(sigma_a2,a2)
	round_b2 = sig_fig(sigma_b2,b2)
	t1_plot = plt.figure(1)
	plt.plot(mol,t1,'ro',mol1,t11,'-b')
#	plt.ylim([0,1.2])
	plt.title('$1/T_1$ $vs.$ $M$')
	plt.xlabel('$Molarity$')
	plt.ylabel('$T_1 (ms)$')
	plt.text(0.6,0.1, \
		'Best fit equation:\nY = $a * x + b$\na = '+ \
		str(round_a1[1])+' +/- '+str(round_a1[0])+'\nb = '+ \
		str(round_b1[1])+' +/- '+str(round_b1[0]))
	plt.grid()
	t2_plot = plt.figure(2)
	plt.plot(mol,t2,'ro',mol2,t22,'-b')
#	plt.errorbar(t2,mol,xerr=sigma_t2)
#	plt.ylim([0,1.2])
	plt.title('$1/T_2$ $vs.$ $M$')
	plt.xlabel('$Molarity$')
	plt.ylabel('$T_2 (ms)$')
	plt.text(0.6,0.1, \
		'Best fit equation:\nY = $a * x + b$\na = '+ \
		str(round_a2[1])+' +/- '+str(round_a2[0])+'\nb = '+ \
		str(round_b2[1])+' +/- '+str(round_b2[0]))
	plt.grid()
	t1_plot.show()
	t2_plot.show()
	raw_input()
main()
