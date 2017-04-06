#!/usr/bin/env python
import sys
from scipy.optimize import curve_fit
from numpy import log10,floor,sqrt,diag,zeros,linspace,log,pi
import matplotlib.pyplot as plt

def sig_fig(x,y):
	return round(x,-int(floor(log10(abs(x))))), \
			round(y,-int(floor(log10(abs(x)))))

def func(x,a,b):
	return a * x + b

#def carriers():

def colors(i):
	if i == 0: return 'rx'
	elif i == 1: return 'bx'
	elif i == 2: return 'gx'
	elif i == 3: return 'cx'
	elif i == 4: return 'mx'
	elif i == 5: return 'yx'
	elif i == 6: return 'kx'
	elif i == 7: return 'ro'
	elif i == 8: return 'bo'
	elif i == 9: return 'go'
	else: return 'ko'

def fit_colors(i):
	if i == 0: return '-r'
	elif i == 1: return '-b'
	elif i == 2: return '-g'
	elif i == 3: return '-c'
	elif i == 4: return '-m'
	elif i == 5: return '-y'
	elif i == 6: return '-k'
	elif i == 7: return '-.r'
	elif i == 8: return '-.b'
	elif i == 9: return '-.g'
	else: return '-.k'

def room_temp():
	fn = open(sys.argv[-1],'r')
	x = []
	y = []
	for i in fn.readlines():
		d = i.split()
		x.append(float(d[0]))
		y.append(float(d[1]))
	p0 = (1,1)
	param,pcov = curve_fit(func,x,y,p0)
	a,b = param
	sigma_a,sigma_b = sqrt(diag(pcov))
	round_a = sig_fig(sigma_a,a)
	round_b = sig_fig(sigma_b,b)
	print 'Room Temperature Hall Voltage Parameters:\n'+ \
		'a = '+str(round_a[1])+' +/- '+str(round_a[0])+ \
		'\nb = '+str(round_b[1])+' +/- '+str(round_b[0])+ \
		'\n================================================='
	x2 = linspace(-0.67,0.67,200)
	y2 = func(x2,a,b)
	rt_plot,rt = plt.subplots(1)
	rt.plot(x,y,'rx',x2,y2,'-b')
	rt.set_title(r'Room Temperature Hall Voltage vs. $\vec B$')
	rt.set_xlabel('Magnetic Field (T)')
	rt.set_ylabel('Hall Voltage (mV)')
	rt_plot.show()

def vanderpauw():
	fn = open(sys.argv[-2],'r')
	raw_data = []
	temperature = []
	for i in fn.readlines():
		if i[0] == "#": 
			continue
		elif i[0] == "S" or i[0] == "s":
			a,b,scale = i.split()
			scale = float(scale)
			continue
		elif i[0] == "C" or i[0] == "c":
			a,b,current = i.split()
			current = float(current)
			continue
		d = i.split(";")
		temperature.append(float(d[0]))
		raw_data.append((float(d[1]),float(d[2])))
	ln2 = log(2.)
	resistivity = zeros(len(raw_data))
	print 'Temperature and Resistivity results'
	for i in range(len(raw_data)):
		#print raw_data[i][0]
		r_val = (raw_data[i][0]-raw_data[i][1])/(raw_data[i][0]+raw_data[i][1])
		f = 1 - (r_val**2)*(ln2/2.) - (r_val**4)*((ln2**2/4)-(ln2**3/12))
		resistivity[i] = (pi*T/ln2) * \
					(((raw_data[i][0]+raw_data[i][1])*scale)/(2*current)) * f
		print temperature[i],resistivity[i]
	print '================================================='
	res_plot,res = plt.subplots(1)
	res.plot(temperature,resistivity,'bx')
	res.set_title("Resistivity vs. Temperature")
	res.set_xlabel("Temperature (K)")
	res.set_ylabel("Resistivity ($\Omega \cdot m$)")
	res_plot.show()
	return temperature[:-2],resistivity[:-2]

def carriers(all_val,temp,resistivity):
	#all_val = [bfa,bfa_sigma,aea,aea_sigma]
	bf_n = zeros(len(all_val[0]))
	ae_n = zeros(len(all_val[0]))
	bf_mob_d = zeros(len(all_val[0]))
	ae_mob_d = zeros(len(all_val[0]))
	#print len(bf_n),len(temp)
	for i in range(len(all_val[0])):
		bf_n[i] = CURRENT/(1.6e-19*all_val[0][i]*T)
		#print all_val[0][i],all_val[2][i]
		ae_n[i] = CURRENT/(1.6e-19*all_val[2][i]*T)
		bf_mob_d[i] = (all_val[0][i]*T) / (CURRENT*resistivity[i])
		ae_mob_d[i] = (all_val[2][i]*T) / (CURRENT*resistivity[i])
		#print temp[i],bf_mob_d[i]
	#bf_car_plot,bf_car = plt.subplots(1)
	#ae_car_plot,ae_car = plt.subplots(1)
	#bf_mob_plot,bf_mob = plt.subplots(1)
	#ae_mob_plot,ae_mob = plt.subplots(1)
	bf_fig,(bf_car,bf_mob) = plt.subplots(nrows = 2)
	ae_fig,(ae_car,ae_mob) = plt.subplots(nrows = 2)
	bf_car.plot(temp,bf_n)
	ae_car.plot(temp,ae_n)
	bf_car.set_title('Carrier density vs. Temperature\nConfiguration 1,2')
	bf_car.set_xlabel('Temperature (K)')
	bf_car.set_ylabel('Carrier density ($m^{-3}$)')
	ae_car.set_title('Carrier density vs. Temperature\nConfiguration 3,4')
	ae_car.set_xlabel('Temperature (K)')
	ae_car.set_ylabel('Carrier density ($m^{-3}$)')
	bf_mob.plot(temp,bf_mob_d)
	ae_mob.plot(temp,ae_mob_d)
	bf_mob.set_title('Carrier Mobility vs. Temperature\nConfiguration 1,2')
	bf_mob.set_xlabel('Temperature (K)')
	bf_mob.set_ylabel('Carrier Mobility ($\Omega \cdot C / m^{2}$)')
	ae_mob.set_title('Carrier Mobility vs. Temperature\nConfiguration 3,4')
	ae_mob.set_xlabel('Temperature (K)')
	ae_mob.set_ylabel('Carrier Mobility ($\Omega \cdot C / m^{2}$)')
	#bf_car_plot.show()
	#ae_car_plot.show()
	#bf_mob_plot.show()
	#ae_mob_plot.show()
	ae_fig.subplots_adjust(hspace = 0.5)
	ae_fig.show()
	bf_fig.show()

# This function will use all of the input files and extract the magnetic
# field for all of them along with the respective voltages.
# This script is designed for two differen input voltages from the files.
# Will apply a linear fit with scipy.optimize.curve_fit for the linear
# function y = ax+b solving for a and b along with their uncertainties.
# Later it will plot them by use of a for loop and with the colors
# and fit_colors functions above it will select a color for each base on the
# for loop index.
def hall_voltage():
	current = 0
	scale = 0
	bfield = []
	bf_voltage = []
	ae_voltage = []
	for j in range(len(sys.argv[1:-2])):
		fn = open(sys.argv[j+1],'r')
		bfield.append([])
		bf_voltage.append([])
		ae_voltage.append([])
		for i in fn.readlines():
			if i[0] == "#":
				continue
			elif i[0] == "C" or i[0] == "c":
				a,b,current = i.split()
				continue
			elif i[0] == "S" or i[0] == "s":
				a,b,scale = i.split()
				continue
			d = i.split(";")
			bfield[j].append(float(d[0]))
			bf_voltage[j].append(float(d[1]))
			ae_voltage[j].append(float(d[2]))
	p0 = (1,1)
	# bf?/ae? store the whole value
	# where round_* stores only the rounded value for output
	bfa = 		zeros(len(bfield))
	bfa_sigma = 	zeros(len(bfield))
	round_bfa =	zeros((len(bfield),2))
	bfb = 		zeros(len(bfield))
	bfb_sigma = 	zeros(len(bfield))
	round_bfb = 	zeros((len(bfield),2))
	aea = 		zeros(len(bfield))
	aea_sigma = 	zeros(len(bfield))
	round_aea = 	zeros((len(bfield),2))
	aeb = 		zeros(len(bfield))
	aeb_sigma = 	zeros(len(bfield))
	round_aeb = 	zeros((len(bfield),2))
	bf_main, bf = plt.subplots(1)
	ae_main, ae = plt.subplots(1)
	print 'Hall Voltage / Magnetic field fits'
	for i in range(len(bfield)):
		param,pcov = curve_fit(func,bfield[i],bf_voltage[i],p0)
		bfa[i],bfb[i] = param
		bfa_sigma[i],bfb_sigma[i] = sqrt(diag(pcov))
		param,pcov = curve_fit(func,bfield[i],ae_voltage[i],p0)
		aea[i],aeb[i] = param
		aea_sigma[i],aeb_sigma[i] = sqrt(diag(pcov))
		bfx = linspace(-0.6678,0.6678,500)
		bfy = func(bfx,bfa[i],bfb[i])
		aex = linspace(-0.6678,0.6678,500)
		aey = func(aex,aea[i],aeb[i])
		bf.plot(bfield[i],bf_voltage[i],colors(i),label=sys.argv[i+1])
		bf.plot(bfx,bfy,fit_colors(i),label=sys.argv[i+1]+' Best-Fit')
		ae.plot(bfield[i],ae_voltage[i],colors(i),label=sys.argv[i+1])
		ae.plot(aex,aey,fit_colors(i),label=sys.argv[i+1]+' Best-Fit')
		round_bfa[i] = sig_fig(bfa_sigma[i],bfa[i])
		round_bfb[i] = sig_fig(bfb_sigma[i],bfb[i])
		round_aea[i] = sig_fig(aea_sigma[i],aea[i])
		round_aeb[i] = sig_fig(aeb_sigma[i],aeb[i])
		print "Best-Fit line for "+sys.argv[i+1]+"\nBF configuration y = "+ \
			str(bfa[i])+" + "+str(bfb[i])+"\nAE configuration y = "+ \
			str(aea[i])+" + "+str(aeb[i])
	print '==============================================='
	bf.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
	bf.set_xlabel('Magnetic Field (T)')
	bf.set_ylabel('Hall Voltage (mV)')
	bf.set_title(r'$\vec B$ vs. $V_{Hall}$'+'\nConfiguration 1,2')
	#bf.set_xlim([-0.67,0.67])
	ae.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
	ae.set_xlabel('Magnetic Field (T)')
	ae.set_ylabel('Hall Voltage (mV)')
	ae.set_title(r'$\vec B$ vs. $V_{Hall}$'+'\nConfiguration 3,4')
	bf_main.show()
	ae_main.show()
	return bfa,bfa_sigma,aea,aea_sigma

# This function serves to administer the flow of control from one function to
# the next.
# Also, makes it look neater when all of the code to perform one function is
# contained inside of one single function rather than all the code to perform
# all functions is in one function.
# raw_input() at the end of the program serves the purpose of pausing the 
# program before all of the plots created are deleted and the program terminates
def main():
	if len(sys.argv) < 2:
		print "ERROR Cannot parse command.\n"+ \
				"USAGE: "+sys.argv[0]+" [temp voltages] [vanderpauw voltages] [Hall voltage RT]\n"+ \
				"Recommended to make data files as data* to make\n"+ \
				"it easier to input the data. Make the last two\n"+ \
				"the files to find resistivity and the room\ntemperature Hall Voltage."
		sys.exit()
	room_temp()
	temp,resistivity = vanderpauw()
	all_val = hall_voltage()
	moar_stuff = carriers(all_val,temp,resistivity)
	raw_input()
CURRENT = 1.e-6
SCALE = 1.e-3
T = 300.e-6
main()
