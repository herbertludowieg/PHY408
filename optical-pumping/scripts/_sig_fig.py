from numpy import log10,floor

def sig_fig(x,y):
  return round(x,-int(floor(log10(abs(x))))), \
         round(y,-int(floor(log10(abs(x)))))
