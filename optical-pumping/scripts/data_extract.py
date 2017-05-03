#!/usr/bin/env python
import numpy as np
import sys

def main():
  fn = open(sys.argv[1])
  start = 0
  raw = []
  for i in fn.readlines():
    d = i.split(',')
    if d[0] == "Vertical Scale":
      volt_scale = float(d[1])
    elif d[0] == "Horizontal Scale":
      time_scale = float(d[1])
    d[-3] = d[-3].strip()
    d[-2] = d[-2].strip()
    x = float(d[-3])
    y = float(d[-2])
    if y <= 0:
      start = 1
    if start:
      raw.append((x,y))
  time = []
  signal = []
  j = 0
  for i in range(len(raw)):
    #print raw[i][1]
    if raw[i][1] > 0. and raw[i-1][1] <= 0.:
      time.append([])
      signal.append([])
      while ( raw[i][1] > 0. and i < len(raw) ):
        #time.append([])
        #signal.append([])
        time[j].append(raw[i][0])
        signal[j].append(raw[i][1])
        i += 1
      j += 1
    else:
      continue
  #print len(time)
  max_sig = []
  for i in range(len(time)):
    if len(time[i]) < 10:
      continue
    #print time[i]
    for j in range(int(len(time[i])*0.45),int(len(time[i])*0.55)):
      if signal[i][j] >= signal[i][j-1] and signal[i][j] > signal[i][j-2]:
        if signal[i][j] >= signal[i][j+1] and signal[i][j] >= signal[i][j+2]:
           #and signal[i][j] > signal[i][j+3]:
          max_sig.append((time[i][j],signal[i][j]))
  for n in range(len(time)):
    for i in range(len(time)):
      if len(time[i]) < 10:
        continue
      j = 1
      while ( j < len(max_sig)):
        print time[i][0],time[i][-1],max_sig[j][0],max_sig[j-1][0]
        if max_sig[j][0] < time[i][-1] and max_sig[j][0] > time[i][0]:
          print 1
          if max_sig[j-1][0] < time[i][-1] and max_sig[j-1][0] > time[i][0]:
            print 2
            if max_sig[j][1] <= max_sig[j-1][1]:
              max_sig.remove(max_sig[j])
            else:
              max_sig.remove(max_sig[j-1])
        j += 1
  print max_sig #len(max_sig)
  freq = 0.
  ave = 0.
  for i in range(len(max_sig)-1):
    ave += abs(max_sig[i][0]-max_sig[i+1][0])
  freq = ave / (len(max_sig)-1)
  print freq

main()
