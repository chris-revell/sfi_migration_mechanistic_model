#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

#Script to plot the output data from seabird_sim.py
#Takes folder containing position data file as a command line argument.

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import os

#Command line input gives path to directory containing datafolders
parameterpath = argv[1]
#Find all datafolders in parameter directory
runs = [f for f in os.listdir(parameterpath) if os.path.isdir(os.path.join(parameterpath,f))]

#Create map of bird path on basemap
from mpl_toolkits.basemap import Basemap
map = Basemap(projection="robin",lon_0=0)
map.fillcontinents(color='coral',lake_color='aqua')
map.drawmapboundary(fill_color='aqua')

for i in runs:
    positiondata = np.genfromtxt(os.path.join(parameterpath,i,"latlongdata.csv"),delimiter=",")

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

    lats = positiondata[:,1]
    lons = positiondata[:,2]
    x,y=map(lons,lats)
    map.scatter(x,y,color=colors,s=0.1)

plt.savefig(os.path.join(parameterpath,parameterpath.split("/")[-1]+"_map.png"),format='png',dpi=600)
