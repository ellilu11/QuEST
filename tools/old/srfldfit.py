import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Parameters
nparts = 25
nsteps = 2000
ntrials = 10
dt = 0.5e-2
hbar = 0.65821193

dirstr = '../outsr/cube_1um/'+str(nparts)+'dots/'
minbeta = 2
maxbeta = 7
betastrall = ['0','0.1','0.2','0.4','0.8','1.6','3.2']
betastr = betastrall[minbeta:maxbeta]
nbetas = len(betastr)

getnint = 0
getfit = 1

# Prepare fit
def fldfunc( t, a, b, td ):
    return a / np.cosh( b * ( t - td ) )

tiall = [1.00, 1.50, 1.20, 1.10, 1.05, 1.00, 0.95] # 25, 50 dots
# tiall = [1.00, 1.20, 1.15, 1.10, 1.00, 0.95, 0.90] # 100 dots
#tiall = [1.00, 1.10, 1.10, 1.00, 1.00, 0.95, 0.90] # 200 dots
# tiall = [1.00, 1.10, 1.05, 1.00, 1.00, 0.95, 0.90] # 400 dots

ti = tiall[minbeta:maxbeta]

tfall = [1.00, 2.50, 2.00, 1.50, 1.30, 1.20, 1.15] # 25, 50 dots
tfall = [1.00, 2.00, 1.50, 1.40, 1.25, 1.20, 1.10] # 100, 200 dots
tfall = [1.00, 1.50, 1.40, 1.30, 1.25, 1.20, 1.10] # 400 dots

# tfall = [1.00, 1.50, 1.35, 1.20, 1.10, 1.10] # 50, 100 dots 
# tfall = [1.00, 1.20, 1.20, 1.10, 1.10, 1.05] # 200 dots 
tf = tfall[minbeta:maxbeta]

tdata = np.linspace( 0, nsteps*dt, nsteps )

fitfile = open('fit.dat','a')

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

for beta in range(nbetas) :

    # Get interacting fld
    #fld = np.loadtxt('../outtest/fld_int_beta10.dat',delimiter=',')
    fld = np.genfromtxt( \
        dirstr+'fld_int_beta'+betastr[beta]+'.dat',delimiter=',', dtype=None, usecols=range(0,nparts) )
    fldAbs = np.zeros( (ntrials, nsteps, nparts) )

    for trial in range(ntrials) :
        for step in range(nsteps) :
            fldAbs[trial,step,:] = fld[trial*nsteps+step,:]

    fldAbsMean = np.mean( fldAbs, axis=(0,2) )

    # Get non-interacting fld (optional)
    if getnint :
        fldNint = np.genfromtxt( \
            dirstr+'fld_nint_beta'+betastr[beta]+'.dat',delimiter=',', dtype=None, usecols=range(0,nparts) )
        fldNintAbs = np.zeros( (ntrials, nsteps, nparts) )

        for trial in range(nsteps) :
            for step in range(nparts) :
                fldNintAbs[trial,step,:] = fldNint[trial*nsteps+step,:]

        fldNintAbsMean = np.mean( fldNintAbs, axis=(0,2) )
    
    # Fit & plot
    if getfit :
        stepmin = int(ti[beta]/dt)
        stepargmax = stepmin + np.argmax( fldAbsMean[stepmin:int(tf[beta]/dt)] )
        step = stepargmax

        while fldAbsMean[step] >= fldAbsMean[step+1] and step < nsteps-1 :
            step = step + 1
        stepmax = np.min( [step, int(tf[beta]/dt)] );

        print(stepmin*dt,stepmax*dt)

        tdatasub = tdata[stepmin:stepmax]
        popt, pcov = curve_fit( fldfunc, tdatasub, fldAbsMean[stepmin:stepmax], maxfev = 5000 )
        perr = np.sqrt(np.diag(pcov))

        stepmax = int(tf[0]/dt)
        plt.plot( tdata[:stepmax], fldAbsMean[:stepmax] )
        plt.plot( tdatasub, fldfunc(tdatasub, *popt), '--', \
            label = 'beta: %5.3f, A: %5.3f, 1/tau: %5.3f, td: %5.3f' \
                % ( float(betastr[beta]), popt[0], popt[1], popt[2] ) )

        print(nparts, betastr[beta], popt[0], popt[1], popt[2], perr[0], perr[1], perr[2], file=fitfile)
    else : 
        # plt.plot( tdata, fldAbsMean, tdata, fldNintAbsMean, 'r-' )
        plt.plot( tdata, fldAbsMean, label = 'beta: %5.3f' % ( float(betastr[beta]) ) )
        # plt.xlim( (0, 2) )
 
plt.figure(1)
plt.xlabel('Time (ps)')
plt.ylabel('|E| (meV / um, e = 1 )')
plt.ylim(0, 2e4)
plt.legend()
plt.show()

fitfile.close()

