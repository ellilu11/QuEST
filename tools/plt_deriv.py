#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
ndots = 1
beta = '1.0' #'1.79e-04'
dir = '../outsr/testderiv/'

derivfile0 = dir+'fld_'+str(ndots)+'dots.dat'
derivfile1 = dir+'fldanl_'+str(ndots)+'dots.dat'

def read_deriv(filename,dt) :
  deriv_raw = np.loadtxt(filename)
  nderivs = 4

  ntimes = deriv_raw.shape[0]
  ndots = int(deriv_raw.shape[1]/nderivs)

  deriv = np.zeros( (ntimes, ndots, nderivs) )

  for j in range(nderivs) :
    idx = [i for i in range(j*ndots,(j+1)*ndots, 1)]
    deriv[:,:,j] = deriv_raw[:, idx]

  deriv_mean = np.mean( deriv, axis=(1) )
  return deriv, deriv_mean

def l2_relerror(data0, data1) :
  ntimes = data0.shape[0];

  sqerr = 0
  sqsum = 0

  for time in range(ntimes) :
    sqerr += (data0[time] - data1[time])**2
    sqsum += data0[time]**2

  return np.sqrt(sqerr / sqsum)

def main() :

  ti, tf = 0, 10
 
  dt0 = 1e-3
  tincr = 1
  dt = dt0 * tincr
  tdata = np.linspace( ti, tf, (tf-ti)/dt )
	 
  deriv, deriv_mean = read_deriv(derivfile0,dt0)
  deriv_anl, deriv_anl_mean = read_deriv(derivfile1,dt0)

  nderivs = 4
  # Plot
  fig = plt.figure()
  get_mean = 0
  if get_mean == 1 :
    for i in range(nderivs) :
      plt.subplot(1,5,i+1)
      plt.plot( tdata, deriv_mean[:,i] )

  elif get_mean == 0 :
    dot = 0
    for i in range(nderivs) :
      plt.subplot(1,5,i+1)
#      plt.plot( tdata, deriv[:,dot,i] - deriv_anl[:,dot,i] ) 
#      plt.plot( tdata, deriv[:,dot,i] )
      plt.plot( tdata, deriv[:,dot,i],
                tdata, deriv_anl[:,dot,i] )
#    plt.xlim( (0, 10) )
#    plt.semilogy()
      print( l2_relerror(deriv_anl[:,dot,i], deriv[:,dot,i]) )

  plt.legend(['Interp','Analytic'])
  plt.show()

  # Rel l2 errors

if __name__ == '__main__':
  main()
