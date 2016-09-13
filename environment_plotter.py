import matplotlib.pyplot as pyplot
import numpy as np
import os
from sys import argv

#Folder path for datafile given at command line
importfolderpath = argv[1]

datafiles = [f for f in os.listdir(importfolderpath) if os.path.isfile(os.path.join(importfolderpath, f)) and f[-4:]=='.CSV']
print(datafiles)

i=0

for environmentfile in datafiles:
    environment_array = np.genfromtxt(os.path.join(importfolderpath,environmentfile),delimiter=",",dtype='float')
    environment_shape = np.shape(environment_array)
    for i in range(0,environment_shape[0]):
        for j in range(0,environment_shape[1]):
            if environment_array[i,j] > 1:
                environment_array[i,j] = 0.0
            else:
                pass
    i=i+1
    pyplot.figure(i)
    pyplot.imshow(environment_array,cmap="BuGn")
    outfilename = os.path.join(importfolderpath,environmentfile[0:-4]+".png")
    pyplot.savefig(outfilename)
    i=i+1
    pyplot.figure(2)
    flat_env = environment_array.flatten()
    pyplot.hist(flat_env,100)
    outfilename = os.path.join(importfolderpath,environmentfile[0:-4]+"_hist.png")
    pyplot.savefig(outfilename)
