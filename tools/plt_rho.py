#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
ndots = 2
dir = '../build/out/beta1.79e-04/'
#dir = '../build/out/beta0/'

rhofile = dir+'rho_'+str(ndots)+'dots_newhistory.dat'

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

def fft_rho(rho_data, ti, tf, dt) :

  ntimes = rho_data.shape[0]

  ti_idx = int(ti / dt)
  tf_idx = int(tf / dt) - 1

#  rho_fft = np.zeros( (tf_idx - ti_idx, ndots) )

  rho_fft = np.fft.fft( rho_data[ti_idx:tf_idx] )

  return rho_fft


def main() :

  ti, tf = 0, 20000
 
  dt0 = 5e-3
  tincr = 50
  dt = dt0 * tincr
  tdata = np.linspace( ti, tf, (tf-ti)/dt )
  ntimes = tdata.shape[0]

  rho = np.zeros( (ntimes, ndots, 4) )
  rho_mean = np.zeros( (ntimes, 4) )

  rho[:,:,:3], rho_mean[:,:3] = read_rho(rhofile)
  rho[:,:,3], rho_mean[:,3] = abs_rho(rho[:,:,:])
#
  # Plot
  fig = plt.figure()
  get_mean = 0

  if get_mean == 1 :
    for i in range(4) :
      plt.subplot(1,4,i+1)
      plt.plot( tdata, rho_mean[:,i] )

  elif get_mean == 0 :
    dot = 0
    for i in range(4) :
      plt.subplot(1,4,i+1)
      plt.plot( tdata, rho[:,dot,i] )
#    plt.semilogy()
  else :
    ti_fft = 50
    tdata_fft = np.linspace( 0, 1.0 / dt, (tf-ti_fft)/dt+1 )
    rhofft = abs( fft_rho( rho_mean[:,1], ti_fft, tf, dt ) )
    plt.plot( tdata_fft[:len(tdata_fft)-2], rhofft ) 
#    plt.semilogx()
    plt.xlim( (0, 2e-2) )
    plt.xlabel('f (1/ps)')
 

  plt.show()

if __name__ == '__main__':
  main()
