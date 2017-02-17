#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

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

#Create map of bird path on basemap
from mpl_toolkits.basemap import Basemap
map = Basemap(projection="robin",lon_0=0)
map.fillcontinents(color='coral',lake_color='aqua')
map.drawmapboundary(fill_color='aqua')
x,y=map(lons,lats)
map.scatter(x,y,color="black",s=0.1)

plt.savefig(os.path.join(argv[1],"map.png"),format='png',dpi=600)
