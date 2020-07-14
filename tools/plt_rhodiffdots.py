#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
ndots = 4
ndots1 = 6
dir = '../outsr/beta1.79e-04/'

rhofile = dir+'rho_'+str(ndots)+'dots_pc_dt5e-3.dat'
rhofile1 = dir+'rho_'+str(ndots1)+'dots_pc_dt5e-3.dat'
#derivfile = dir+'rhoanl_'+str(ndots)+'dots.dat'

def read_rho(filename) :
  rho_raw = np.loadtxt(filename)

  ntimes = rho_raw.shape[0]
  ndots = int(rho_raw.shape[1]/3)

  rho = np.zeros( (ntimes, ndots, 3) )

  for j in range(3) :
    idx = [i for i in range(j, rho_raw.shape[1], 3)]
    rho[:,:,j] = rho_raw[:, idx] 

  rho_mean = np.mean( rho, axis=(1) )
  return rho, rho_mean

def abs_rho(rho) :

  ntimes = rho.shape[0]
  ndots = rho.shape[1]

  rho_abs = np.zeros( (ntimes, ndots) )

  for time in range(ntimes) :
    for dot in range(ndots) :
      rho_abs[time,dot] = np.sqrt( rho[time,dot,1]**2 + rho[time,dot,2]**2 )

  rho_abs_mean = np.mean( rho_abs, axis=(1) )

  return rho_abs, rho_abs_mean


def read_deriv(filename) :
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

  ti, tf = 0, 50000
 
  dt0 = 5e-3
  tincr = 10
  dt = dt0 * tincr
  tdata = np.linspace( ti, tf, (tf-ti)/dt )
	 
  rho, rho_mean = read_rho(rhofile)
  rho1, rho1_mean = read_rho(rhofile1)
#  deriv, deriv_mean = read_deriv(derivfile)# analytic derivatives!

  rho_abs, rho_abs_mean = abs_rho(rho)
  rho1_abs, rho1_abs_mean = abs_rho(rho1)
#
  # Plot
  fig = plt.figure()
  for i in range(3) :
    plt.subplot(1,4,i+1)
    plt.plot( tdata, rho_mean[:,i],
              tdata, rho1_mean[:,i] )
    plt.xlim( (0,50000) )
  plt.subplot(1,4,4)
  plt.plot( tdata, rho_abs_mean,
            tdata, rho1_abs_mean )
  plt.xlim( (0,50000) )

#    plt.semilogy()
#    print( l2_relerror(rho[:,dot,0], deriv[:,dot,deriv_order]) )

  plt.legend(['N=4','N=6'])
  plt.show()

if __name__ == '__main__':
  main()
