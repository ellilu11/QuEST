#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
dir = '/mnt/home/luelliot/QuEST/build/dots/'
filestr = dir+'dots'+sys.argv[2]+'.cfg'
dotfile = open(filestr, 'w')

def main() :

  omega = 4823.67
  T1 = 8.31
  T2 = 2*T1
  dipx = 0.002536 * 1.0
  dipy = 0
  dipz = 0
  
  c0 = 299.792458
  n = 1.54 
  c = c0/n
  waveleng = 2.0*np.pi*c/omega
  sigma = 1.0*waveleng
  sigsq = sigma*sigma
  mu = [0,0,0]
  cov = [[sigsq, 0, 0], [0, sigsq, 0], [0, 0, sigsq]]

  ndots = int(sys.argv[1])
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
