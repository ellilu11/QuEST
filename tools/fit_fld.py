import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Times New Roman'

#rc('text', usetex=True)
#plt.style.use('ggplot')

#Parameters
#nparts = [10,20,40]
#ntrials = [80,40,40]
nparts = [2]
ntrials = [1]

nnparts = len(nparts)
tincr0 = 1 # timestep incr in QuEST
tincr = 1 # timestep incr here
dt0 = 1e-4 # timestep in QuEST
dt = dt0*tincr0 # timestep here
hbar = 0.65821193

betastr=['1.79e-04']
# dirstr = '../outsr/'
dirstr = '../outsr/beta'+betastr[0]+'/'
# t1_1e-4/'

ti = 0
tf = 10

getfld = 1
getrho = 0
getfit = 0

def fitfunc( t, a, b, c, t0 ):
    # return a * t + b
    # return a / np.cosh( b * ( t - t0 ) ) - 1
    return a * ( 1 - 1.0/2.0 * ( b * ( t - t0 ) )**2 + 5.0/24.0 * ( b * ( t - t0 ) )**4 ) + c

fitfile = open('fit.dat','a')
# fit_ti = [2000,1000,750,600]
fit_ti = [ti,ti,ti,ti]
fit_tf = tf

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
                        rho00[trial,step,:] = rhoRaw[trial*nsteps+step,part*3]
                    else :
                        rho00[trial,step] = rhoRaw[trial*nsteps+step]

        rho00Mean = np.mean( rho00, axis=(0,2) )

        plt.plot( tdatasub, \
                  1.0 - 2.0 * rho00Mean[int(ti/dt):int(tf/dt)], \
                  label='N='+str(nparts[idx]) )
        #Fit
        if getfit :
            # Extract values from tdata & fldAbsMean modulo incr
            nsteps2 = int(nsteps/tincr)
            tdataFit = np.linspace( 0, tmax, nsteps2 )
            dataFit = np.zeros( nsteps2 )
            for step in range(nsteps2) :
                dataFit[step] = rho00Mean[step*tincr]
            
            popt, pcov = curve_fit( fitfunc, tdataFit[int(fit_ti[idx]/dt/tincr):int(fit_tf/dt/tincr)], \
                                             dataFit[int(fit_ti[idx]/dt/tincr):int(fit_tf/dt/tincr)], maxfev=1000000 )
            perr = np.sqrt(np.diag(pcov))

            #plt.plot( tdataFit[int(fit_ti[idx]/dt/tincr):int(tf/dt/tincr)], \
            #    dataFit[int(fit_ti[idx]/dt/tincr):int(tf/dt/tincr)], \
            #    label='N='+str(nparts[idx]) )
 
            plt.plot( tdataFit, \
                  fitfunc(tdataFit, *popt), '--', \
                  label = 'A: %5.3f, 1/tau: %5.3f, t0: %5.3f' \
                    % ( popt[0], popt[1], popt[3] ) )
            print(nparts[idx], popt[0], popt[1], popt[3], perr[0], perr[1], perr[3], file=fitfile)


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


