#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from _sig_fig import *

def light_func(x,a,b,k):
  return a*np.exp(-k*(x-4e16))+b

def dens_func(x,a,b,k):
  return a*np.exp(k*(x-87))+b

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
  print "==========END FIT DATA=========================="
  round_a = sig_fig(sigma_a,a)
  round_b = sig_fig(sigma_b,b)
  round_k = sig_fig(sigma_k,k)
  #round_z = sig_fig(sigma_z,z)
  return a,b,k,round_a,round_b,round_k

def temperature_vs_density(density):
  x = np.zeros(len(density))
  y = np.zeros(len(density))
  for i in range(len(density)):
    x[i] = density[i][0]
    y[i] = density[i][1]
  xx = np.linspace(x[0],x[-1],100)
  n,m,l = 1e18,8.5e-2,7e15
  yy = dens_func(xx,n,l,m)
  a,b,k,round_a,round_b,round_k = exponential_3(x,y,(n,l,m),dens_func)
  yyy = dens_func(xx,a,b,k)
  ax,bx = plt.subplots(1)
  bx.plot(x,y,'ro')
  bx.plot(xx,yy,'b-')
  bx.plot(xx,yyy,'r-')
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
  round_z = np.zeros(2)
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
  a,b,k,round_a,round_b,round_k = exponential_3(x,y,(10,1,6e-18),light_func,0.002)
  print "\nFit parameters with uncertainties:"
  print "a = "+str(round_a[1])+" +/- "+str(round_a[0])
  print "b = "+str(round_b[1])+" +/- "+str(round_b[0])
  print "k = "+str(round_k[1])+" +/- "+str(round_k[0])
  print "z = 4e16"
  print "*************************************************"
  yyy = light_func(xx,a,b,k)
  ax,bx = plt.subplots(1)
  bx.plot(x,y,'ro',label="Data")
  bx.plot(xx,yyy,'r-',label='Fit')
  bx.plot(xx,yy,'b-',label='Guess Fit')
  bx.set_xlabel("Density ($m^{-3}$)")
  bx.set_ylabel("Light Intensity ($V$)")
  bx.set_title("Density vs. Light Intensity")
  bx.legend()
  #bx.set_ylim([0,12])
  #bx.set_xlim([4.0e18])
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
  density_vs_light(density)
  temperature_vs_density(density)
  raw_input()

main()
