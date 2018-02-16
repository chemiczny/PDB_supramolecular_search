# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 15:51:05 2018
@author: michal
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt("logs/anionPiLigandFullExtracted.log", usecols = (5, 6), skiprows = 1)
distances = data[:, 0]
angles = data[:, 1]

#for i in range(len(angles)):
#    if angles[i] > 90:
#        angles[i] =180 - angles[i]

plt.figure()
plt.plot(angles, distances, '.')
plt.xlabel("angle")
plt.ylabel("distance")

plt.figure()
n, bins, patches = plt.hist(angles, 100,  facecolor='g')
plt.xlabel('Angle')
plt.title('Histogram of angles')
plt.grid(True)
plt.show()

plt.figure()
n, bins, patches = plt.hist(distances, 100, facecolor='g')
plt.xlabel('Distance')
plt.title('Histogram of distances')
plt.grid(True)
plt.show()

#bins_no = 30
#hist, xedges, yedges = np.histogram2d(angles, distances, bins = bins_no )
#
#x,y,z = [], [], []
#
#for i in range(bins_no):
#      
#    for j in range(bins_no):
#        x.append( xedges[i])  
#        y.append(yedges[j])
#        z.append(hist[i][j])
#        
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
#ax.plot_trisurf(x, y, z,cmap=plt.cm.CMRmap )
#
#plt.show()