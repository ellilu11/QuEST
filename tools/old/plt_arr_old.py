import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Times New Roman'

#Parameters
nparts = [10,20,40,80]
ntrials = [80,40,39,40]

#nparts = [10]
#ntrials = [80]

nnparts = len(nparts)
tincr0 = 1000 # timestep incr in QuEST
dt0 = 5e-3 # timestep in QuEST
dt = dt0*tincr0 # timestep here
hbar = 0.65821193

betastr=['1.79e-04']
# dirstr = '../outsr/'
dirstr = '../outsr/beta'+betastr[0]+'/'
# t1_1e-4/'

ti = 0
tf = 10000

getfld = 0
getrho = 1

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

if getfld :
    for idx in range(nnparts) :

        fldRaw = np.genfromtxt( \
            dirstr+'fld_'+str(nparts[idx])+'dots.dat',delimiter=' ', dtype=None )

        nsteps = int(np.size(fldRaw,0) / ntrials[idx])
        tmax = nsteps*dt
        tdata = np.linspace( 0, tmax, nsteps )
        tdatasub = np.linspace( ti, tf, (tf-ti)/dt )

        fldAbs = np.zeros( (ntrials[idx], nsteps, nparts[idx]) )

        for trial in range(ntrials[idx]) :
            for step in range(nsteps) :
                if nparts[idx] > 1 :
                    fldAbs[trial,step,:] = fldRaw[trial*nsteps+step,:]
                else :
                    fldAbs[trial,step] = fldRaw[trial*nsteps+step]

        fldAbsMean = np.mean( fldAbs, axis=(0,2) )

        plt.plot( tdatasub, \
              fldAbsMean[int(ti/dt):int(tf/dt)], \
              label='N='+str(nparts[idx]) )

elif getrho :
     for idx in range(nnparts) :

        rhoRaw = np.genfromtxt( \
            dirstr+'rho_'+str(nparts[idx])+'dots.dat',delimiter=' ', dtype=None )

        nsteps = int(np.size(rhoRaw,0) / ntrials[idx])
        tmax = nsteps*dt
        tdata = np.linspace( 0, tmax, nsteps )
        tdatasub = np.linspace( ti, tf, (tf-ti)/dt )

        rho00 = np.zeros( (ntrials[idx], nsteps, nparts[idx]) )

        for trial in range(ntrials[idx]) :
            for step in range(nsteps) :
                for part in range(nparts[idx]) :
                    if nparts[idx] > 1 :
                        rho00[trial,step,part] = rhoRaw[trial*nsteps+step,part*3]
                    else :
                        rho00[trial,step] = rhoRaw[trial*nsteps+step]

        rho00Mean = np.mean( rho00, axis=(0,2) )

        plt.plot( tdatasub, \
                  1.0 - 2.0 * rho00Mean[int(ti/dt):int(tf/dt)], \
                  label='N='+str(nparts[idx]) )

plt.xlabel('Time (ps)', fontsize=32)
ax.tick_params(labelsize=24)
if getfld :
    plt.ylabel('|E| (meV / $\mu$m, e = 1 )', fontsize=32)
#    plt.ylim(0, 1e4)
else :
    plt.ylabel('w', fontsize=32)
    plt.ylim(-1, 1)
plt.legend( prop={'size': 32} )
plt.show()


