#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Times New Roman'

# dirstr = '../outsr/beta1.79e-04/'
dirstr = '../outsr/APS/'

linestyle_str = [ 'solid', 'dotted', 'dashed', 'dashdot' ]

def read_fld(filename,ntrials) :
  fld_raw = np.loadtxt(filename)

  ntimes = int(fld_raw.shape[0]/ntrials)
  ndots = fld_raw.shape[1]

  fld = np.zeros( (ntimes, ndots, ntrials) )

  for k in range(ntrials) :
    idx_row = [i for i in range(ntimes*k, ntimes*(k+1), 1)]
    fld[:,:,k] = fld_raw[idx_row,:]
    
  fld_mean = np.mean( fld, axis=(1,2) )
  return fld, fld_mean


def read_rho(filename,ntrials) :
  rho_raw = np.loadtxt(filename)

  ntimes = int(rho_raw.shape[0]/ntrials)
  ndots = int(rho_raw.shape[1]/3)

  rho = np.zeros( (ntimes, ndots, ntrials, 3) )

  rho_over_trials = np.zeros( (rho_raw.shape[0], ndots, 3) )

  for j in range(3) :
    idx_col = [i for i in range(j, rho_raw.shape[1], 3)]
    rho_over_trials[:,:,j] = rho_raw[:,idx_col]
    for k in range(ntrials) :
      idx_row = [i for i in range(ntimes*k, ntimes*(k+1), 1)]
      rho[:,:,k,j] = rho_over_trials[idx_row,:,j]
    
  rho_mean = np.mean( rho, axis=(1,2) )
  return rho, rho_mean

def main() :

  dt0 = 5e-3
  tincr = 100
  dt = dt0 * tincr
  ti, tf = 10, 10000
  tdata = np.linspace( ti, tf, (tf-ti)/dt )
	
  ndots = [10,20,40,80]
#  ndots = [10,20]
  ntrials = [80,40,40,20]
  nndots = len(ndots)

  fig = plt.figure(dpi=150)
  ax = fig.add_subplot(1,1,1)

  pltflag = 1

  for idx in range(nndots) :

    if pltflag :
      rhofile = dirstr+'rho_'+str(ndots[idx])+'dots.dat';
      rho, rho_mean = read_rho(rhofile, ntrials[idx])

      plt.plot( tdata, 1.0 - 2.0 * rho_mean[int(ti/dt):int(tf/dt),0],
                linestyle=linestyle_str[idx],
                label='N='+str(ndots[idx]) )
      plt.ylim( (-1, 1) )
      plt.ylabel('w', fontsize=24)
    else :
      fldfile = dirstr+'fld_'+str(ndots[idx])+'dots.dat';
      fld, fld_mean = read_fld(fldfile, ntrials[idx])

      plt.plot( tdata, fld_mean[int(ti/dt):int(tf/dt)],
                linestyle=linestyle_str[idx],
                label='N='+str(ndots[idx]) )
      plt.ylim( (0, 10000) )
      plt.ylabel('|E| (meV / $\mu$m, e = 1 )', fontsize=24)
 
	# plt.semilogy()
  ax.tick_params(labelsize=24)
  plt.xlabel('Time (ps)', fontsize=24)
  plt.legend( prop={'size': 24}, loc='upper right' )
  plt.show()

if __name__ == '__main__':
  main()
