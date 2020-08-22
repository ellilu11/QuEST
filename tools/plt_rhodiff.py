#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
ndots = 2
nvars = 4

dir = '../build/out/beta1.79e-04/'
#dir = '../build/out/beta0/'

rhofile0 = dir+'rho_'+str(ndots)+'dots_dt5e-3.dat'
rhofile1 = dir+'rho_'+str(ndots)+'dots_dt5e-3_newhistory.dat'

rhofile = np.array([rhofile0, rhofile1])
nfiles = rhofile.shape[0]

def read_rho(filename) :
  rho_raw = np.loadtxt(filename)

  ntimes = rho_raw.shape[0]
  ndots = int(rho_raw.shape[1]/nvars)

  rho = np.zeros( (ntimes, ndots, nvars) )

  for j in range(nvars) :
    idx = [i for i in range(j, rho_raw.shape[1], nvars)]
    rho[:,:,j] = rho_raw[:, idx] 

  rho_mean = np.mean( rho, axis=(1) )
  return rho, rho_mean

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

def main() :

  ti, tf = 0, 20000
 
  dt0 = 5e-3
  tincr = 50
  dt = dt0 * tincr
  tdata = np.linspace( ti, tf, (tf-ti)/dt )
  ntimes = tdata.shape[0]

  rho = np.zeros( (ntimes, ndots, nvars, nfiles) )
  rho_mean = np.zeros( (ntimes, nvars, nfiles) )

  for i in range(nfiles) :
    rho[:,:,:,i], rho_mean[:,:,i] = read_rho(rhofile[i])

  rho_mean_err = abs( rho_mean[:,:,1] - rho_mean[:,:,0] ) / 1.0 #abs (rho_mean[:,:,0])
  rho_err = abs( rho[:,:,:,1] - rho[:,:,:,0] ) / 1.0 #abs (rho[:,:,:,0])

  # Plot
  fig = plt.figure()
  get_mean = 1

  if get_mean == 1 :
    for i in range(4) :
      plt.subplot(2,4,i+1)
      plt.plot( tdata, rho_mean[:,i,:] )
  
      plt.subplot(2,4,i+5)
      plt.plot( tdata, rho_mean_err[:,i] )
  #    plt.semilogy()

  else :
    dot = 0
    for i in range(4) :
      plt.subplot(2,4,i+1)
      plt.plot( tdata, rho[:,dot,i,:] )
      plt.legend(['Fix (dt=5e-5)', 'Rot (dt=5e-3)'])

      plt.subplot(2,4,i+5)
      plt.plot( tdata, rho_err[:,dot,i] )
      #plt.semilogy()
      plt.legend(['Rel err'])

  #  plt.legend(['Rot','Fix'])
  plt.show()
 
if __name__ == '__main__':
  main()
