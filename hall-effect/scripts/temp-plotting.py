#!/usr/bin/env python
import sys
from scipy.optimize import curve_fit
from numpy import log10,floor,sqrt,diag,zeros,linspace,log,pi
import matplotlib.pyplot as plt

# will reduce the significant digits of the inputs by rounding to the first
# significant digit of the uncertainty
def sig_fig(x,y):
	return round(x,-int(floor(log10(abs(x))))), \
			round(y,-int(floor(log10(abs(x)))))

# simple function for the linear fitting with curve_fit
def func(x,a,b):
	return a * x + b

def func_0(x,a):
	return a*x

# function that makes a switch statement in essence and returns line
# parameters for plotting
def colors(i):
	if i == 0: return (1,0,0)
	elif i == 1: return (0,0,1)
	elif i == 2: return (0.3,0.5,0.3)
	elif i == 3: return (0.3,1,0.9)
	elif i == 4: return (1,0.2,0.9)
	elif i == 5: return (1,1,0)
	elif i == 6: return (0.5,1,0)
	elif i == 7: return (0.5,0,1)
	elif i == 8: return (0.8,0.8,0.5)
	elif i == 9: return (0.8,0.2,0.5)
	else: return (0.5,0.5,0.5)

def magvscurrent():
	x = [0,5,10,15,20,25,30,35]
	y = [0,0.1056,0.2125,0.3179,0.4265,0.5258,0.6058,0.6678]
	param,pcov = curve_fit(func_0,x,y)
	a = param
	a_sigma = sqrt(pcov)
	xx = linspace(0,35,50)
	yy = func_0(xx,a)
	round_a = sig_fig(a_sigma,a)
	#round_b = sig_fig(b_sigma,b)
	print "Magnetic Filed vs Current Fit\ny = ax + b\nParameters: \n"+ \
		"a = "+str(round_a[1])+" +/- "+str(round_a[0])
		#"\nb = "+str(round_b[1])+" +/- "+str(round_b[0])
	print "==============================================="
	mag_plot,mag = plt.subplots(1)
	mag.plot(x,y,'ro',xx,yy,'-r')
	mag.set_title(r'$\vec B$ vs. Current')
	mag.set_xlabel('Current (A)')
	mag.set_ylabel('Magnetic Field (T)')
	mag_plot.show()

# function that plots the room temperature hall voltage as a function of
# magnetic field the slope of theis tells the type of doping for the semi
# conductor.
# the data for this must be in the last argument on the command line
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
	rt.plot(x,y,'ro',x2,y2,'-r')
	rt.set_title(r'Room Temperature Hall Voltage vs. $\vec B$')
	rt.set_xlabel('Magnetic Field (T)')
	rt.set_ylabel('Hall Voltage (mV)')
	rt_plot.show()

# function that will take the data from the second to last input file
# and find the resistivity using the van der pauw method.
# the file must be second to last
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
		r_val = (raw_data[i][0]-raw_data[i][1])/ \
				(raw_data[i][0]+raw_data[i][1])
		f = 1 - (r_val**2)*(ln2/2.) - (r_val**4)*((ln2**2/4)-(ln2**3/12))
		resistivity[i] = (pi*T/ln2) * \
					(((raw_data[i][0]+raw_data[i][1])*scale)/(2*current)) * f
		print temperature[i],resistivity[i]
	print '================================================='
	res_plot,res = plt.subplots(1)
	res.plot(temperature,resistivity,'bo')
	res.set_title("Resistivity vs. Temperature")
	res.set_xlabel("Temperature (K)")
	res.set_ylabel("Resistivity ($\Omega \cdot m$)")
	res_plot.show()
	return temperature[:-2],resistivity[:-2]

# finds the carrier mobility and carrier density.
# receives input from hall voltage
def carriers(all_val,temp,resistivity):
	#all_val = [tot_a,tot_a_sigma]
	bf_n = zeros(len(all_val[0]))
	bf_mob_d = zeros(len(all_val[0]))
	for i in range(len(all_val[0])):
		bf_n[i] = CURRENT/(-1.6e-19*all_val[0][i]*SCALE*T)
		bf_mob_d[i] = (all_val[0][i]*SCALE*T) / (CURRENT*resistivity[i])
	p0 = (0.1,0.001)
	# find the best fit line for bf_n
	param,pcov = curve_fit(func,temp,bf_n,(1e20,2.8e20))
	bf_n_a,bf_n_b = param
	bf_n_a_sigma,bf_n_b_sigma = sqrt(diag(pcov))
	bf_nx = linspace(80,260,250)
	bf_ny = func(bf_nx,bf_n_a,bf_n_b)
	# find the best fit line for bf_mob
	param,pcov = curve_fit(func,temp,bf_mob_d,p0)
	bf_mob_a,bf_mob_b = param
	bf_mob_a_sigma,bf_mob_b_sigma = sqrt(diag(pcov))
	bf_moby = func(bf_nx,bf_mob_a,bf_mob_b)
	# plot absolutely everything
	bf_fig,(bf_car,bf_mob) = plt.subplots(nrows = 2)
	bf_car.plot(temp,bf_n,'ro',bf_nx,bf_ny,'-r')
	bf_car.set_title('Carrier density vs. Temperature')
	bf_car.set_xlabel('Temperature (K)')
	bf_car.set_ylabel('Carrier density ($m^{-3}$)')
	bf_mob.plot(temp,bf_mob_d,'ro',bf_nx,bf_moby,'-r')
	bf_mob.set_title('Carrier Mobility vs. Temperature')
	bf_mob.set_xlabel('Temperature (K)')
	bf_mob.set_ylabel('Carrier Mobility ($\Omega \cdot C / m^{2}$)')
	bf_fig.subplots_adjust(hspace = 0.5)
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
	total_voltages = []
	tot_a = zeros(len(bfield))
	tot_a_sigma = zeros(len(bfield))
	round_tot_a = zeros((len(bfield),2))
	tot_b = zeros(len(bfield))
	tot_b_sigma = zeros(len(bfield))
	round_tot_b = zeros((len(bfield),2))	
	total_volt_plot,total_volt = plt.subplots(1)
	# due to our experimental set-up the ae_voltages list must be inverted
	# to get a similar trend in both the ae_voltages and bf_voltages
	# this is then passed through the for loop to average them out 
	# and later plot the average
	for i in range(len(bfield)):
		for j in range(len(bfield[0])/2):
			temp = 0
			temp = ae_voltage[i][-j-1]
			ae_voltage[i][-j-1] = ae_voltage[i][j]
			ae_voltage[i][j] = temp
	print 'Hall Voltage / Magnetic field fits'
	for i in range(len(bfield)):
		total_voltages.append([])
		for j in range(len(bfield[0])):
			total_voltages[i].append((bf_voltage[i][j]+ae_voltage[i][j])/2.)
	for i in range(len(bfield)):
		# find best fit for average
		param,pcov = curve_fit(func,bfield[i],total_voltages[i],p0)
		tot_a[i],tot_b[i] = param
		tot_a_sigma[i],tot_b_sigma[i] = sqrt(diag(pcov))
		# create best fit line using parameters
		totx = linspace(-0.6678,0.6678,500)
		toty = func(totx,tot_a[i],tot_b[i])
		# plot all of the data seperately according to temperature
		# use colors function to define colors
		total_volt.plot(bfield[i],total_voltages[i],'o',color=colors(i),
						label=sys.argv[i+1][-4:])
		total_volt.plot(totx,toty,'-',color=colors(i),
						label=sys.argv[i+1][-4:]+' Fit')
		# feed parameters into sig_fig function to round to correct sig figs
		round_tot_a[i] = sig_fig(tot_a_sigma[i],tot_a[i])
		round_tot_b[i] = sig_fig(tot_b_sigma[i],tot_b[i])
		# print parameters
		print "Best-Fit line for "+sys.argv[i+1][-4:]+"\nAverage y = ax + b"+ \
			"\nParameters: \na = "+str(round_tot_a[i][1])+" +/- "+ \
			str(round_tot_a[i][0])+"\nb =  "+str(round_tot_b[i][1])+" +/- "+ \
			str(round_tot_b[i][0])
	print '=================================================='
	# plot properties of the Average hall voltage plot
	total_volt.set_title(r'Average Hall Voltage vs. $\vec B$')
	total_volt.set_xlabel('Magnetic Field (T)')
	total_volt.set_ylabel('Hall Voltage (mV)')
	# puts the legend outside the plot area
	total_volt.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
	# configures the whitespace in the plot
	total_volt_plot.subplots_adjust(right=0.75,top=0.92,left=0.09,bottom=0.08)
	total_volt_plot.show()
	return tot_a,tot_a_sigma

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
	magvscurrent()
	room_temp()
	temp,resistivity = vanderpauw()
	all_val = hall_voltage()
	moar_stuff = carriers(all_val,temp,resistivity)
	raw_input()
CURRENT = 1.e-6
SCALE = 1.e-3
T = 300.e-6
main()
