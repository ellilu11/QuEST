import numpy as np
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd

#from scipy.optimize import curve_fit

# Parameters
nparts = 2
ntrials = 1
dt0 = 1e-5
tincr0 = 100
dt = dt0 * tincr0
hbar = 0.65821193

betastr=['1.79e-04']
# dirstr = '../outsr/'
dirstr = '../outsr/'

ti = 0
tf = 50

getfld = 0
getrho = 1

# Get fld
#fldRaw = pd.read_csv( dirstr+'fld_'+str(nparts)+'dots.dat' )
# print(fldRaw)
if getfld :
    fldRaw = np.genfromtxt( \
        dirstr+'fld_'+str(nparts)+'dots.dat',delimiter=' ', dtype=None )

    nsteps = int(np.size(fldRaw,0) / ntrials)
    tmax = nsteps*dt
    tdata = np.linspace( 0, tmax, nsteps )

    fldAbs = np.zeros( (ntrials, nsteps, nparts) )

    for trial in range(ntrials) :
        for step in range(nsteps) :
            if nparts > 1 :
                fldAbs[trial,step,:] = fldRaw[trial*nsteps+step,:]
            else :
                fldAbs[trial,step] = fldRaw[trial*nsteps+step]

    fldAbsMean = np.mean( fldAbs, axis=(0,2) )

# Get rho
if getrho :
    rhoRaw = np.genfromtxt( \
        dirstr+'rho_'+str(nparts)+'dots.dat',delimiter=' ', dtype=None )
    nsteps = int(np.size(rhoRaw,0) / ntrials)
    tmax = nsteps*dt
    tdata = np.linspace( 0, tmax, nsteps )
    rho = np.zeros( (ntrials, nsteps, nparts, 3) )

    for trial in range(ntrials) :
        for step in range(nsteps) :
            for part in range(nparts) :
                for i in range(3) :
                    if nparts > 1 :
                        rho[trial,step,part,i] = rhoRaw[trial*nsteps+step,part*3+i]
                    else :
                        rho[trial,step,0,i] = rhoRaw[trial*nsteps+step,i]

    rhoMean = np.mean( rho, axis=(0) )

# Plot
fig = plt.figure()
tdatasub = np.linspace( ti, tf, (tf-ti)/dt )

dot = 0

if getrho :
    plt.plot( tdatasub, rhoMean[int(ti/dt):int(tf/dt),dot] )

    plt.title('w')
    plt.xlabel('Time (ps)')
    plt.ylabel('w')
else :
    plt.plot( tdatasub, fldAbsMean[int(ti/dt):int(tf/dt)] )
    plt.title('|E|')
    plt.xlabel('Time (ps)')
    plt.ylabel('|E| (meV / um, e = 1 )')

#plt.xlim( (0, 2) )
#plt.legend()
plt.show()


