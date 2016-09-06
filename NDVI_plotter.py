import matplotlib.pyplot as pyplot
import numpy as np
import os
from sys import argv

#Folder path for datafile given at command line
importfolderpath = argv[1]

datafiles = [f for f in os.listdir(importfolderpath) if os.path.isfile(os.path.join(importfolderpath, f)) and f[-1]=='t']
print(datafiles)

for NDVIfile in datafiles:
    NDVI_array = np.genfromtxt(os.path.join(importfolderpath,NDVIfile),delimiter=" ", skip_header=1, usecols = range(1,2001))
    pyplot.imshow(NDVI_array,cmap="Greens")
    outfilename = os.path.join(importfolderpath,NDVIfile[0:-4]+".png")
    pyplot.savefig(outfilename)
