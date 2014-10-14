# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 10:23:06 2014

@author: christoforos
"""

import numpy as np
from scipy.optimize import leastsq

################ fit_leastsq( x, y, dy, fit_func, residual_func, param_guess ) ################
def fit_leastsq( x, y, dy, fit_func, residual_func, param_guess ):
  #Ensure that the data is stored in np.arrays:
  #can sanity check with curve_fit
  x  = np.array(x)
  y  = np.array(y)
  dy = np.array(dy)
  
  fitResult = leastsq( residual_func, param_guess, args=(x,y,dy,fit_func), full_output =1 )
  popt, pcov, infodict, mesg, ier = fitResult
  
  chiSq_red = (residual_func(popt, x, y, dy, fit_func)**2).sum()/(len(y)-len(param_guess))
  
  perr = [np.sqrt(pcov[i][i]) for i in range(len(pcov)) ]
  
  return popt, perr, chiSq_red, pcov
  
############################## residuals_weighted( p, x, y, dy ) ##############################
def residuals_weighted( p, x, y, dy, fit_func ):
  #Ensure that the data is stored in np.arrays:
  x  = np.array(x)
  y  = np.array(y)
  dy = np.array(dy)

  y_fit = fit_func(x, *p) # '*' is for unpacking the array of parameters
  err = (y - y_fit)/dy
  
  return err
