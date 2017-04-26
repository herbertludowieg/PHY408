#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from _sig_fig import *

def func(x,a,b,k,z):
  return a*np.exp(-k*(x-z))+b

# Finds exponential fit to a*e^(-k*(x-z))+b
def exponential_4(x,y,p0,sigma=0):
  #print x,y,p0
  if sigma == 0:
    param,pcov = curve_fit(func,x,y,p0)
  else:
    param,pcov = curve_fit(func,x,y,p0,sigma=sigma)
  print "Initial guess parameters:"
  print "a = "+str(p0[0])
  print "b = "+str(p0[1])
  print "k = "+str(p0[2])
  print "z = "+str(p0[3])
  print "Calculated parameters with curve_fit function:"
  print "Format (a,b,k,z)"
  print param
  print "Covariance matrix from curve_fit function:"
  print pcov
  a,b,k,z = param
  sigma_a,sigma_b,sigma_k,sigma_z = np.sqrt(abs(np.diag(pcov)))
  round_a = sig_fig(sigma_a,a)
  round_b = sig_fig(sigma_b,b)
  round_k = sig_fig(sigma_k,k)
  round_z = sig_fig(sigma_z,z)
  return a,b,k,z,round_a,round_b,round_k,round_z

def open_file(fn):
  x = []
  y = []
  for i in fn.readlines():
    if i == "#":
      continue
    

def temp(density):
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
  xx = np.linspace(-4.0e18,4.0e18,1000)
  yy = func(xx,10,1,6e-18,2e16)
  #yy = func(xx,12.8356290,1.25988601,6.08534284e-18,-4.93643723e16)
  #print yy
  a,b,k,z,round_a,round_b,round_k,round_z = exponential_4(x,y,(10,1,6e-18,2e16),0.002)
  #print round_a,round_b,round_k,round_z
  yyy = func(xx,a,b,k,z)
  ax,bx = plt.subplots(1)
  bx.plot(x,y,'ro',xx,yyy,'r-',xx,yy,'b-')
  bx.set_ylim([0,12])
  bx.set_xlim([-4.0e18,4.0e18])
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
  temp(density)
  raw_input()

main()
