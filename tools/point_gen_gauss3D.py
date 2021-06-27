#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
dir = '/mnt/home/luelliot/QuEST/build/dots/'
filestr = dir+'dots'+sys.argv[1]+'.cfg'
dotfile = open(filestr, 'w')

def main() :

  omega = 4823.67
  T1 = 20
  T2 = 2*T1
  dipx = 0.002536
  dipy = 0
  dipz = 0
  
  c0 = 299.792458
  waveleng = 2.0*np.pi*c0/omega
  mu = [0,0,0]
  cov = [[0.5*waveleng, 0, 0], [0, 0.5*waveleng, 0], [0, 0, 0.5*waveleng]]

  ndots = 200
  X, Y, Z = np.random.multivariate_normal(mu, cov, ndots).T
  R = np.sqrt(X**2+Y**2+Z**2)
  idxsort = np.argsort(R)
  #print(idxsort)
 
  XYZ = np.vstack((X,Y,Z)).T
  XYZ = XYZ[idxsort]
  #print(XY)

  #plt.plot(X, Y)

  for dot in range(ndots) :
    print( XYZ[dot,0], XYZ[dot,1], XYZ[dot,2], omega, T1, T2, dipx, dipy, dipz, file=dotfile )

  #plt.show()

if __name__ == '__main__':
  main()
