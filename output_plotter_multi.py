#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

# Christopher Revell, University of Cambridge, 2017

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

parametercomponents = parameterpath.split("_")
if "lat-52.3_lon169.1" in parameterpath:
    location = "Campbell Island"
elif "lat-51.0_lon-61.1" in parameterpath:
    location = "Falkland Islands"
elif "lat-49.4_lon70.0" in parameterpath:
    location = "Iles Kerguelen"
elif "lat-51.3_lon-75.2" in parameterpath:
    location = "Islas Diego de Almagro"
elif "lat-55.4_lon-69.3" in parameterpath:
    location = "Islas Ildefonso"
else:
    location = "South Georgia"

#Create map of bird path on basemap
from mpl_toolkits.basemap import Basemap
#map = Basemap(projection="robin",lon_0=0)
map = Basemap(projection="cyl",lon_0=0)
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
    map.scatter(x,y,alpha=0.02,color="black",s=0.1)
    #map.scatter(x,y,alpha=0.02,color=colors,s=0.1)

plt.title(location+", a="+parametercomponents[-3][1:]+", kT="+parametercomponents[-2][2:])
plt.savefig(os.path.join(parameterpath,parameterpath.split("/")[-1]+"_map2.png"),format='png',dpi=600,bbox_inches="tight")
