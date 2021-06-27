#!/usr/bin/env python3
import sys
import numpy as np

# Parameters
dir = '/mnt/home/luelliot/QuEST/build/dots/'

dotfile = open(dir+'obss_sphsurf.cfg','a')

def main() :

  num_src = 5
  src_dz = 5.0e-3
  r0 = 1000 # (num_src+2) * src_dz / 2.0
 
  nth = 50
  nph = 100
  nobs = nph*nth

  dph = 2*np.pi/nph
  dlph = r0 * dph
 
  [x0, w0] = np.polynomial.legendre.leggauss(nth)

  area_sum = 0
  #circ_sum = 0
 
  for iph in range(nph) :
    semicirc_sum = 0
    # circ_sum += dlph
    for ith in range(nth) :
      phi = iph*dph
      theta = np.pi/2 * (x0[ith] + 1)

      # Assume dipole along x
      x = r0 * np.cos( theta )
      y = r0 * np.sin( theta ) * np.cos( phi )
      z = r0 * np.sin( theta ) * np.sin( phi )
     
      dlth = np.pi / 2 * w0[ith] * r0
      weight = dlth * dlph * np.sin(theta)

      dotdata = np.array([x, y, z, weight, 0, 0, 0, 0, 0])
      np.savetxt(dotfile, dotdata, fmt='%1.9f', newline=' ')
      dotfile.write("\n")

      semicirc_sum += dlth
      area_sum += weight

  #circ = 2*np.pi*r0
  #print("Sum of line elements: "+str(circ_sum))
  #print("Actual circumference: "+str(circ))

  semicirc = np.pi*r0
  print("Sum of line elements: "+str(semicirc_sum))
  print("Actual semicircumference: "+str(semicirc))

  area = 4.0*np.pi*r0**2
  print("Sum of area elements: "+str(area_sum))
  print("Actual surface area: "+str(area))

if __name__ == '__main__':
  main()
