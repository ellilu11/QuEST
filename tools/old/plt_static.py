import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
#from scipy.optimize import curve_fit

# Parameters
nparts = 2
nsteps = 10240
dt = 1e-3
tmult = 1
hbar = 0.65821193

dirstr = '../outstatic/'
#dirstr = '../outsr/beta'+betastr[0]+'/'
tdata = np.linspace( 0, tmult*nsteps*dt, nsteps )

# Get interacting fld
fldRaw = np.genfromtxt( \
    dirstr+'fld_dir.dat',delimiter=' ', dtype=None, usecols=range(0,3*nparts) )

fld = np.zeros( (nsteps, 3, nparts) )
fldAbs = np.zeros( (nsteps, nparts) )

if nparts > 1:
    for step in range(nsteps) :
        for part in range(nparts) :
            fld[step,0,part] = fldRaw[step,3*part]
            fld[step,1,part] = fldRaw[step,3*part+1]
            fld[step,2,part] = fldRaw[step,3*part+2]
            fldAbs[step,part] = LA.norm(fld[step,:,part])

# fldAbsMean = np.mean( fldAbs, axis=(0,2) )

# Get analytic fld
fldRaw = np.genfromtxt( \
    dirstr+'fld_anl.dat',delimiter=' ', dtype=None, usecols=range(0,3*nparts) )

fldAnl = np.zeros( (nsteps, 3, nparts) )
fldAbsAnl = np.zeros( (nsteps, nparts) )

if nparts > 1:
    for step in range(nsteps) :
        for part in range(nparts) :
            fldAnl[step,0,part] = fldRaw[step,3*part]
            fldAnl[step,1,part] = fldRaw[step,3*part+1]
            fldAnl[step,2,part] = fldRaw[step,3*part+2]
            fldAbsAnl[step,part] = LA.norm(fldAnl[step,:,part])

# Plot
# plt.plot( tdata, fldAbs[:,1], tdata, fldAbsAnl[:,1])

fig = plt.figure()
titlestr = ['E_x','E_y','E_z']
part = 1

for i in range(1,7) :
    plt.subplot(2, 3, i)
 
    if i < 4 :
        plt.plot( tdata, fld[:,i-1,part], tdata, fldAnl[:,i-1,part])
        plt.title(titlestr[i-1])
    else :
        plt.plot( tdata, fldAnl[:,i-4,part] - fld[:,i-4,part])
        plt.title(titlestr[i-4]+' (Anl - Num)')
 
# plt.xlim( (0, 2) )
# plt.ylim(0, 2e4)

# plt.xlabel('Time (ps)')
# plt.ylabel('|E| (meV / um, e = 1 )')
plt.legend()
plt.show()


