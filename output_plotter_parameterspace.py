#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

# Christopher Revell, University of Cambridge, 2016

#Script to plot the output data from seabird_sim.py
#Takes folder containing all datafolders for a given location, and uses data to plot the full parameter space.

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import Basemap

#Command line input gives path to directory containing datafolders
datapath = argv[1]

#Find all datafolders in breeding location directory
parametersets = [f for f in os.listdir(datapath) if os.path.isdir(os.path.join(datapath,f))]

kT_order  = []
a_order   = []
kT_values = []
a_values  = []

# Find all values of a and kT used in this parameter space
for i in parametersets:
    a_str  = i.split("_")[2]
    kT_str = i.split("_")[3]
    a  = float(a_str[1:])
    kT = float(kT_str[2:])
    kT_order.append(kT)
    a_order.append(a)
    if a in a_values:
        pass
    else:
        a_values.append(a)
    if kT in kT_values:
        pass
    else:
        kT_values.append(kT)
kT_values.reverse()

fig1 = plt.figure()
subplotnumbera = len(kT_values)
subplotnumberb = len(a_values)

for n,parameterpair in enumerate(parametersets):

    #Find run folders
    runs = [f for f in os.listdir(os.path.join(datapath,parameterpair)) if os.path.isdir(os.path.join(datapath,parameterpair,f))]

    #Create map of bird path on basemap as subplot
    subplotnumberc = a_values.index(a_order[n])+kT_values.index(kT_order[n])*len(a_values)
    print(subplotnumbera,subplotnumberb,subplotnumberc+1)
    ax1 = fig1.add_subplot(subplotnumbera,subplotnumberb,subplotnumberc+1)
    ax1.set_title("a="+str(a_order[n])+" kT="+str(kT_order[n]))
    map = Basemap(projection="cyl",lon_0=0)
    map.fillcontinents(color='coral',lake_color='aqua')
    map.drawmapboundary(fill_color='aqua')

    for i in runs:
        positiondata = np.genfromtxt(os.path.join(datapath,parameterpair,i,"latlongdata.csv"),delimiter=",")
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


fig1.set_tight_layout(True)
#fig1.subplots_adjust(wspace=0.05, hspace=0.05)
fig1.savefig(os.path.join(datapath,datapath.split("/")[-1]+"_map.png"),format='png',dpi=1200,bbox_inches="tight")
