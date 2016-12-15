#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3
#Script to plot the output data from seabird_sim.py
#Takes folder containing position data file as a command line argument.

from sys import argv,exit
import numpy as np
import matplotlib.pyplot as plt
import os

positiondata = np.genfromtxt(os.path.join(argv[1],"positiondata.csv"),delimiter=",")

#Import ground map
earth = np.genfromtxt("earth1440x720.CSV",delimiter=",")
lattice_shape = np.shape(earth)
for i in range(0,lattice_shape[0]):
    for j in range(0,lattice_shape[1]):
        if earth[i,j] == 99999:
            earth[i,j] = 0
        else:
            earth[i,j] = 1
d_latlong = 180/lattice_shape[0]


#Create plot of bird path on lattice
fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
ax3.imshow(earth,cmap="Greens")
ax3.scatter(positiondata[:,2],positiondata[:,1],s=0.5)
ax3.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
ax3.set_xlim([0,lattice_shape[1]])
ax3.set_ylim([lattice_shape[0],0])
fig2.savefig(os.path.join(argv[1],"lattice.png"),format='png')


#Create map of bird path on basemap
from mpl_toolkits.basemap import Basemap
map = Basemap(projection="robin",lon_0=0)
map.fillcontinents(color='coral',lake_color='aqua')
map.drawmapboundary(fill_color='aqua')
lats = (lattice_shape[0]/2-positiondata[:,1]-0.5)*d_latlong
lons = (positiondata[:,2]+0.5-lattice_shape[1]/2)*d_latlong
x,y=map(lons,lats)
map.scatter(x,y,color="black",s=0.5)
plt.savefig(os.path.join(argv[1],"map.png"),format='png')
