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
  for i in range(len(raw)):
    

main()
