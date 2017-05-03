#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from _sig_fig import *

def light_func(x,a,b,k):
  return a*np.exp(-k*(x-2e16))+b

def dens_func(x,a,k):
  return a*np.exp(k*(x-100))

def period_func(x,a,b,k):
  return a*np.exp(-k*(x-1))+b

# Finds exponential fit to a*e^(-k*(x-z))+b
def exponential_2(x,y,p0,func,sigma=0):
  #print x,y,p0
  if sigma == 0:
    param,pcov = curve_fit(func,x,y,p0)
  else:
    param,pcov = curve_fit(func,x,y,p0,sigma=sigma)
  print "==========BEGIN FIT DATA========================"
  print "Initial guess parameters:"
  print "a = "+str(p0[0])
  print "k = "+str(p0[1])
  print "Calculated parameters with curve_fit function:"
  print "Format (a,k)"
  print param[0],param[1]
  print "Covariance matrix from curve_fit function:"
  print pcov
  a,k = param
  sigma_a,sigma_k = np.sqrt(abs(np.diag(pcov)))
  print "Sigma values:"
  print "Format (a,k)"
  print sigma_a,sigma_k
  print "Percentages:"
  print "Format (a,k)"
  print sigma_a/a*100,sigma_k/k*100
  print "==========END FIT DATA=========================="
  round_a = sig_fig(sigma_a,a)
  round_k = sig_fig(sigma_k,k)
  #round_z = sig_fig(sigma_z,z)
  return a,k,round_a,round_k

# Finds exponential fit to a*e^(-k*(x-z))+b
def exponential_3(x,y,p0,func,sigma=0):
  #print x,y,p0
  if sigma == 0:
    param,pcov = curve_fit(func,x,y,p0)
  else:
    param,pcov = curve_fit(func,x,y,p0,sigma=sigma)
  print "==========BEGIN FIT DATA========================"
  print "Initial guess parameters:"
  print "a = "+str(p0[0])
  print "b = "+str(p0[1])
  print "k = "+str(p0[2])
  print "Calculated parameters with curve_fit function:"
  print "Format (a,b,k)"
  print param[0],param[1],param[2]
  print "Covariance matrix from curve_fit function:"
  print pcov
  a,b,k = param
  sigma_a,sigma_b,sigma_k = np.sqrt(abs(np.diag(pcov)))
  print "Sigma values:"
  print "Format (a,b,k)"
  print sigma_a,sigma_b,sigma_k
  print "Percentages:"
  print "Format (a,b,k)"
  print sigma_a/a*100,sigma_b/b*100,sigma_k/k*100
  print "==========END FIT DATA=========================="
  round_a = sig_fig(sigma_a,a)
  round_b = sig_fig(sigma_b,b)
  round_k = sig_fig(sigma_k,k)
  #round_z = sig_fig(sigma_z,z)
  return a,b,k,sigma_a,sigma_b,sigma_k,round_a,round_b,round_k

def ringing_vs_rfamp():
  period_fn = open("../data-iv/period-data",'r')
  rf_fn = open("../data-iv/rf-amplitudes",'r')
  rb87 = []
  rb85 = []
  rf = []
  bool_rb87 = 0
  bool_rb85 = 0
  for i in period_fn.readlines():
    #print i[:4]
    if i[:4] == "Rb87":
      bool_rb87 = 1
      bool_rb85 = 0
      continue
    elif i[:4] == "Rb85":
      bool_rb85 = 1
      bool_rb87 = 0
      continue
    #print bool_rb87
    if bool_rb87:
      rb87.append(float(i)*1e6)
    elif bool_rb85:
      rb85.append(float(i)*1e6)
  #print rb87,rb85
  for i in rf_fn.readlines():
    rf.append(float(i))
  print "*************************************************"
  print "Fit for Rb 85 period of ringing as a function of\nRF Amplitude"
  round_a_85 = np.zeros(2)
  round_k_85 = np.zeros(2)
  round_b_85 = np.zeros(2)
  a_85,b_85,k_85,n,m,l,round_a_85,round_b_85,round_k_85 = \
                                exponential_3(rf,rb85,(900,300,1),period_func)
  print "a = "+str(round_a_85[1])+" +/- "+str(round_a_85[0])
  print "b = "+str(round_b_85[1])+" +/- "+str(round_b_85[0])
  print "k = "+str(round_k_85[1])+" +/- "+str(round_k_85[0])
  print "z = 1.0"
  print "*************************************************"
  print "*************************************************"
  print "Fit for Rb 87 period of ringing as a function of\nRF Amplitude"
  round_a_87 = np.zeros(2)
  round_k_87 = np.zeros(2)
  round_b_87 = np.zeros(2)
  a_87,b_87,k_87,n,m,l,round_a_87,round_b_87,round_k_87 = \
                                exponential_3(rf,rb87,(550,200,1),period_func)
  print "a = "+str(round_a_87[1])+" +/- "+str(round_a_87[0])
  print "b = "+str(round_b_87[1])+" +/- "+str(round_b_87[0])
  print "k = "+str(round_k_87[1])+" +/- "+str(round_k_87[0])
  print "z = 1.0"
  print "*************************************************"

  xx = np.linspace(rf[0],rf[-1],1000)
  #y_85 = period_func(xx,900,1,300)
  #y_87 = period_func(xx,550,1,200)
  y_85 = period_func(xx,a_85,b_85,k_85)
  y_87 = period_func(xx,a_87,b_87,k_87)
  ax,bx = plt.subplots(1)
  bx.plot(rf,rb85,'ro',label="Rb 85")
  bx.plot(xx,y_85,'r-',label="Fit Rb 85")
  bx.plot(rf,rb87,'bo',label="Rb 87")
  bx.plot(xx,y_87,'b-',label="Fit Rb 87")
  bx.legend()
  bx.set_title("Period of ringing vs. RF Amplitude")
  bx.set_xlabel("RF Amplitude ($V$)")
  bx.set_ylabel("Period of ringing ($\mu s$)")
  ax.show()

def temperature_vs_density(density):
  x = np.zeros(len(density))
  y = np.zeros(len(density))
  for i in range(len(density)):
    x[i] = density[i][0]
    y[i] = density[i][1]
  xx = np.linspace(x[0],x[-1],100)
  n,m,l = 1e18,8.5e-2,7e15
  yy = dens_func(xx,n,m)
  print "*************************************************"
  print "Fit for density as a function of temperature"
  #a,b,k,round_a,round_b,round_k = exponential_3(x,y,(n,l,m),dens_func)
  a,k,round_a,round_k = exponential_2(x,y,(n,m),dens_func)
  yyy = dens_func(xx,a,k)
  print "\nFit parameters with uncertainties:"
  print "a = "+str(round_a[1])+" +/- "+str(round_a[0])
  print "k = "+str(round_k[1])+" +/- "+str(round_k[0])
  print "z = 100"
  print "*************************************************"
  ax,bx = plt.subplots(1)
  bx.plot(x,y,'ro')
  #bx.plot(xx,yy,'b-')
  bx.plot(xx,yyy,'r-')
  bx.text(10,1.75e19, \
          "Best-Fit Equation:\n$Y = A*e^{k*(x-100)}$\nA = "+ \
          str(round_a[1])+" +/- "+str(round_a[0])+"\nk = "+ \
          str(round_k[1])+" +/- "+str(round_k[0]), \
          color='black',bbox=dict(facecolor='none',edgecolor='black'))
  bx.set_xlabel("Temperature ($^oC$)")
  bx.set_ylabel("Density ($m^{-3}$)")
  bx.set_title("Temperature vs. Density")
  ax.show()

def density_vs_light(density):
  data_file = open("../i/data","r")
  data = []
  x = []
  y = []
  for i in data_file.readlines():
    if i[0] == "#":
      continue
    d=i.split(";")
    x.append(float(d[0]))
    y.append(float(d[1]))
  for i in range(len(x)):
    data.append([])
    data[i].append(x[i])
    data[i].append(y[i])
  density_slim = np.zeros((len(data),len(data[0])))
  for i in range(len(density)):
    #print i,len(data)
    if i >= len(data):
      break
    elif density[i][0] != data[0][0]:
      continue
    for j in range(i,len(data)+i):
      density_slim[j-i] = density[j]
      #print density_slim[j-i],j,i
    break
  #print density_slim
  round_a = np.zeros(2) 
  round_b = np.zeros(2)
  round_k = np.zeros(2)
  sigma = np.zeros(3)
  x = np.zeros(len(data))
  y = np.zeros(len(data))
  for i in range(len(data)):
    x[i] = density_slim[i][1]
    y[i] = data[i][1]
  xx = np.linspace(x[0],x[-1],1000)
  yy = light_func(xx,10,1,6e-18)
  #yy = func(xx,12.8356290,1.25988601,6.08534284e-18,-4.93643723e16)
  #print yy
  print "*************************************************"
  print "Fit for light intensity as a function of density"
  a,b,k,sigma[0],sigma[1],sigma[2],round_a,round_b,round_k = \
                                exponential_3(x,y,(10,1,6e-18),light_func,0.002)
  xsarea = k/L
  sigma_xsarea = sigma[2]/L
  round_xsarea = sig_fig(sigma_xsarea,xsarea)
  print "\nFit parameters with uncertainties:"
  print "a = "+str(round_a[1])+" +/- "+str(round_a[0])
  print "b = "+str(round_b[1])+" +/- "+str(round_b[0])
  print "k = "+str(round_k[1])+" +/- "+str(round_k[0])
  print "z = 2e16"
  print "\nCross sectional area:"
  print "sigma = "+str(round_xsarea[1])+" +/- "+str(round_xsarea[0])
  print "*************************************************"
  yyy = light_func(xx,a,b,k)
  ax,bx = plt.subplots(1)
  bx.plot(x,y,'ro',label="Data")
  bx.plot(xx,yyy,'r-',label='Fit')
  #bx.plot(xx,yy,'b-',label='Guess Fit')
  bx.set_xlabel("Density ($m^{-3}$)")
  bx.set_ylabel("Light Intensity ($V$)")
  bx.set_title("Density vs. Light Intensity")
  #bx.legend()
  bx.text(2.7e18,8.0,\
         "Best-Fit Equation:\n$Y = A*e^{-k*(x-2e16)}+B$\nA = "+ \
         str(round_a[1])+" +/- "+str(round_a[0])+"\nB = "+ \
         str(round_b[1])+" +/- "+str(round_b[0])+"\nk = "+ \
         str(round_k[1])+" +/- "+str(round_k[0]), \
         color='black',bbox=dict(facecolor='none',edgecolor='black'))
  #bx.set_ylim([0,12])
  #bx.set_xlim([4.0e18])
  ax.show()

def magnetic_field(N,I,R):
  return (8*mu0*N*I)/(R*np.sqrt(125))*1e3

def current_vs_field(magnets):
  x = np.linspace(0,3,100)
  y_vert = magnetic_field(magnets['Vertical'][1],x,magnets['Vertical'][0])
  y_horiz = magnetic_field(magnets['Horizontal'][1],x,magnets['Horizontal'][0]) 
  y_sweep = magnetic_field(magnets['Sweep'][1],x,magnets['Sweep'][0])
  ax,bx = plt.subplots(1)
  bx.plot(x,y_vert,'r-',label='Vertical Coils')
  bx.plot(x,y_horiz,'b-',label='Horizontal Coils')
  bx.plot(x,y_sweep,'c-',label='Sweep Coils')
  bx.legend(loc='upper left')
  bx.set_xlabel("Current ($A$)")
  bx.set_ylabel("Magnetic Field ($mT$)")
  bx.set_title("Cuurent vs. Magnetic field")
  ax.show()

def main():
  dens_file = open("density-vs-temp.dat",'r')
  density = []
  x = []
  y = []
  for i in dens_file.readlines():
    if i[0] == "#":
      continue
    d = i.split(";")
    x.append(float(d[0]))
    y.append(float(d[2]))
  for i in range(len(x)):
    density.append([])
    density[i].append(x[i])
    density[i].append(y[i])
  magnets = {}
  mag_file = open("magnet-prop.dat",'r')
  for i in mag_file.readlines():
    if i[0] == '#':
      continue
    d = i.split(';')
    magnets[d[0]] = [float(d[1])*1e-2,float(d[2]),float(d[3]),float(d[4])]
  loop = 1
  while (loop):
    print "/////////////////////////////////////////////////////"
    print "Which would you like to print?"
    print "Enter the name in parentheses for the specific plot."
    print "density versus light intensity? (den v light)"
    print "temperature versus density? (temp v den)"
    print "current versus magnetic field? (curr v field)"
    print "period of ringing versus rf amplitude? (per v rf)"
    print "all? (all)"
    print "Enter quit or hit enter key to exit the program."
    which = raw_input("Enter name here:\n")
    if which == 'all':
      density_vs_light(density)
      temperature_vs_density(density)
      current_vs_field(magnets)
      ringing_vs_rfamp()
    elif which == "den v light":
      density_vs_light(density)
    elif which == "temp v den":
      temperature_vs_density(density)
    elif which == "curr v field":
      current_vs_field(magnets)
    elif which == "per v rf":
      ringing_vs_rfamp()
    elif which == "quit" or which == "":
      break
    else:
      print "\nERROR Can not parse command."
      print "Please try again.\n"

L = 0.033
mu0 = 1.2566370614e-6
main()
