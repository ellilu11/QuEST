#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
ndots = 2
beta = '1.79e-04'
dir = '../outsr/beta'+beta+'/'

rhofile_pc = dir+'rho_'+str(ndots)+'dots_pc_dt1e-3.dat'
rhofile_nt0 = dir+'rho_'+str(ndots)+'dots_nt_dt1e-3.dat'
#rhofile_nt1 = dir+'rho_'+str(ndots)+'dots_nt_dt1e-5_nint.dat'
#rhofile_nt2 = dir+'rho_'+str(ndots)+'dots_nt_dt5e-6.dat'

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

def fft_rho(rho_data, ti, tf, dt) :

  ntimes = rho_data.shape[0]

  ti_idx = int(ti / dt)
  tf_idx = int(tf / dt) - 1

#  rho_fft = np.zeros( (tf_idx - ti_idx, ndots) )

  rho_fft = np.fft.fft( rho_data[ti_idx:tf_idx] )

  return rho_fft

def main() :

  ti, tf = 0, 10000
 
  dt0_pc = 1e-3
  tincr_pc = 1
  dt_pc = dt0_pc * tincr_pc
  tdata = np.linspace( ti, tf, (tf-ti)/dt_pc )
	 
  rho_pc, rho_mean_pc = read_rho(rhofile_pc)
  rho_nt0, rho_mean_nt0 = read_rho(rhofile_nt0)
#  rho_nt1, rho_mean_nt1 = read_rho(rhofile_nt1)
#  rho_nt2, rho_mean_nt2 = read_rho(rhofile_nt2)
 
  rho_abs_pc, rho_abs_mean_pc = abs_rho(rho_pc)
  rho_abs_nt0, rho_abs_mean_nt0 = abs_rho(rho_nt0)
#  rho_abs_nt1, rho_abs_mean_nt1 = abs_rho(rho_nt1)
#  rho_abs_nt2, rho_abs_mean_nt2 = abs_rho(rho_nt2)

  # Plot
  fig = plt.figure()
  dot = 0
  get_mean = 0
  if get_mean == 1 :
    for i in range(3) :
      plt.subplot(1,4,i+1)
      plt.plot( tdata, rho_mean_nt[:,i], 
                tdata, rho_mean_pc[:,i] )
    plt.subplot(1,4,4)
    plt.ylim( (0, 1) )
    plt.plot( tdata, rho_abs_mean_nt,
              tdata, rho_abs_mean_pc )
    plt.ylim( (0, 1) )
  
  elif get_mean == 0 :
    for i in range(3) :
      plt.subplot(1,4,i+1)
#      plt.plot( tdata, rho_nt[:,dot,i] - rho_pc[:,dot,i] )
      plt.plot( tdata, rho_pc[:,dot,i], 
                tdata, rho_nt0[:,dot,i] )
#                tdata, rho_nt1[:,dot,i] )
#                tdata, rho_nt2[:,dot,i] )
      plt.xlim( (0, tf) )
      plt.ylim( (0, 1) )
    plt.subplot(1,4,4)
    plt.plot( tdata, rho_abs_pc[:,dot],
              tdata, rho_abs_nt0[:,dot] )
#              tdata, rho_abs_nt1[:,dot] )
#              tdata, rho_abs_nt2[:,dot] )
    plt.xlim( (0, tf) )
    plt.ylim( (0, 1) )
  else :
    ti_fft = 15
    tdata_fft = np.linspace( 0, 1.0 / dt, (tf-ti_fft)/dt+1 )
    rhofft_nt = abs( fft_rho( rho_nt[:,dot,2], ti_fft, tf, dt ) )
    rhofft_pc = abs( fft_rho( rho_pc[:,dot,2], ti_fft, tf, dt ) )
    plt.plot( tdata_fft, rhofft_nt, 
              tdata_fft, rhofft_pc )
    plt.xlim( (0, 1) )
    plt.xlabel('f (1/ps)')
#  plt.ylim( (-1, 1) )

	# plt.semilogy()
  plt.legend(['P-C (rot)', 'Newton (Rot, 1e-4)','Newton (Rot, 1e-5)', 'Newton (Rot, 5e-6)'])
  plt.show()

if __name__ == '__main__':
  main()
