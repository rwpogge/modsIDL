#!/usr/bin/env python
# PROGRAM: SURFACEFIT.PY
# BY: Scott Adams - sadams@astronomy.ohio-state.edu
# PURPOSE: Given .csv file with randomly sampled x & y coordinates 
#		and correspoending wavelengths return array of 
#		interpolated wavelengths for each X, Y within sampled
#		boundaries
# DATE: 25 Sep 2013
#======================================================================

import asciitable
import numpy as np
from scipy import interpolate

data = asciitable.read('coords.dat')

tck = interpolate.bisplrep(data['XCOORDS'],data['YCOORDS'],data['WAVELENGTH'])
x = np.linspace(int(min(data['XCOORDS'])),int(max(data['XCOORDS'])),int(max(data['XCOORDS']))-int(min(data['XCOORDS']))+1)
y = np.linspace(int(min(data['YCOORDS'])),int(max(data['YCOORDS'])),int(max(data['YCOORDS']))-int(min(data['YCOORDS']))+1)
#x = np.linspace(0,8192-1,8192)
#y = np.linspace(0,3088-1,3088)
fit = interpolate.bisplev(x,y,tck)
#print np.shape(fit)
#print np.shape(x)
#import matplotlib.pyplot as mpl
#fig = mpl.figure()
#mpl.imshow(fit)
#mpl.contour(y,x,fit)
#mpl.colorbar()
#mpl.show()
np.savetxt('fit.dat',fit,delimiter="\t")
