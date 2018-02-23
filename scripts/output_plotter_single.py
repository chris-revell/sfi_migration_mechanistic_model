#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

# Christopher Revell, University of Cambridge, 2016

#Script to plot the output data from seabird_sim.py
#Takes folder containing position data file as a command line argument.

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import os

#Command line input gives path to directory containing data
positiondata = np.genfromtxt(os.path.join(argv[1],"latlongdata.csv"),delimiter=",")
lats = positiondata[:,1]
lons = positiondata[:,2]
colors = []
for i in range(np.shape(positiondata)[0]):
    if positiondata[i,0] < 720:
        colors.append("black")
    elif positiondata[i,0] < 1440:
        colors.append("red")
    elif positiondata[i,0] < 2160:
        colors.append("green")
    else:
        colors.append("blue")
#Create map of bird path on basemap
from mpl_toolkits.basemap import Basemap
#map = Basemap(projection="robin",lon_0=0)
map = Basemap(projection="cyl",lon_0=0)
map.fillcontinents(color='coral',lake_color='aqua')
map.drawmapboundary(fill_color='aqua')
x,y=map(lons,lats)
map.scatter(x,y,color=colors,s=0.1)

plt.savefig(os.path.join(argv[1],"map.png"),format='png',dpi=600,bbox_inches="tight")
