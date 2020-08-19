#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
ndots = 2

#dir = '../build/out/beta1.79e-04/'
dir = '../build/out/beta0/'
#dir = '../build/out/beta1.79e-04_realfld/'

rhofile0 = dir+'rho_'+str(ndots)+'dots_dt5e-3.dat'
rhofile1 = dir+'rho_'+str(ndots)+'dots_dt5e-3_newhistory.dat'
#rhofile0 = dir+'rho_'+str(ndots)+'dots_fixed_dt5e-5.dat'
#rhofile1 = dir+'rho_'+str(ndots)+'dots_dt5e-3.dat'

rhofile = np.array([rhofile0, rhofile1])
nfiles = rhofile.shape[0]

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

def fft_rho(rho_data, ti, tf, dt) :

  ntimes = rho_data.shape[0]

  ti_idx = int(ti / dt)
  tf_idx = int(tf / dt) - 1

#  rho_fft = np.zeros( (tf_idx - ti_idx, ndots) )

  rho_fft = np.fft.fft( rho_data[ti_idx:tf_idx] )

  return rho_fft

def main() :

  ti, tf = 0, 10
 
  dt0 = 5e-3
  tincr = 1
  dt = dt0 * tincr
  tdata = np.linspace( ti, tf, (tf-ti)/dt )
  ntimes = tdata.shape[0]

  rho = np.zeros( (ntimes, ndots, 4, nfiles) )
  rho_mean = np.zeros( (ntimes, 4, nfiles) )

  for i in range(nfiles) :
    rho[:,:,:3,i], rho_mean[:,:3,i] = read_rho(rhofile[i])
    rho[:,:,3,i], rho_mean[:,3,i] = abs_rho(rho[:,:,:,i])

  rho_mean_err = abs( rho_mean[:,:,1] - rho_mean[:,:,0] ) / 1.0 #abs (rho_mean[:,:,0])
  rho_err = abs( rho[:,:,:,1] - rho[:,:,:,0] ) / 1.0 #abs (rho[:,:,:,0])

  # Plot
  fig = plt.figure()
  get_mean = 1

  if get_mean == 1 :
    for i in range(4) :
      plt.subplot(2,4,i+1)
      plt.plot( tdata, rho_mean[:,i,:] )
  #    plt.xlim( (0,10) )
  
      plt.subplot(2,4,i+5)
      plt.plot( tdata, rho_mean_err[:,i] )
  #    plt.semilogy()
  #    plt.xlim( (0,10) )

  elif get_mean == 0 :
    dot = 0
    for i in range(4) :
      plt.subplot(2,4,i+1)
      plt.plot( tdata, rho[:,dot,i,:] )
      plt.legend(['Fix (dt=5e-5)', 'Rot (dt=5e-3)'])
      #plt.xlim( (0,10000) )

      plt.subplot(2,4,i+5)
      #plt.plot( tdata, rho[:,dot+1,i,:] )
      plt.plot( tdata, rho_err[:,dot,i] )
      #plt.semilogy()
      plt.legend(['Rel err'])

  else :
    ti_fft = 50
    tdata_fft = np.linspace( 0, 1.0 / dt, (tf-ti_fft)/dt+1 )
    rhofft = abs( fft_rho( rho_mean[:,1,2], ti_fft, tf, dt ) )
    plt.plot( tdata_fft[:len(tdata_fft)-2], rhofft ) 
    plt.xlabel('f (1/ps)')
  #  plt.legend(['Rot','Fix'])
  plt.show()
 
if __name__ == '__main__':
  main()