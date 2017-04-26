from scipy.optimize import curve_fit
from numpy import sqrt, diag
import _sig_fig

# Finds exponential fit to a*e^(-k*x)+b
def exponential_3(x,y,p0):
  param,pcov = curve_fit(func,x,y,p0)
  a,b,k = param
  sigma_a,sigma_b,sigma_k = sqrt(diag(pcov))
  round_a = _sig_fig.sig_fig(sigma_a,a)
  round_b = _sig_fig.sig_fig(sigma_b,b)
  round_k = _sig_fig.sig_fig(sigma_k,k)
  return a,b,k,round_a,round_b,round_k
