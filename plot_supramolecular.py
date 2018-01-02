# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 15:51:05 2018

@author: michal
"""

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("logs/anionPiLigand.log", skiprows=1, usecols = (5, 6))
distances = data[:, 0]
angles = data[:, 1]

for i in range(len(angles)):
    if angles[i] > 90:
        angles[i] =180 - angles[i]

plt.figure()
plt.plot(angles, distances, '*')
plt.xlabel("angle")
plt.ylabel("distance")

plt.figure()
n, bins, patches = plt.hist(angles, 50,  facecolor='g')
plt.xlabel('Angle')
plt.title('Histogram of angles')
plt.grid(True)
plt.show()

plt.figure()
n, bins, patches = plt.hist(distances, 50, facecolor='g')
plt.xlabel('Distance')
plt.title('Histogram of distances')
plt.grid(True)
plt.show()
