
from sys import argv,exit
import numpy as np
from random import random
from math import exp, acos, sin, cos, pi, radians
import os
import time
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

initialbird_position = (900,720)
kT = 1

#Import chloro data
resources = np.genfromtxt(argv[1],delimiter=",")
lattice_shape = np.shape(resources)

d_latlong = 180/lattice_shape[0]

#Threshold chloro data
resources_filtered = np.zeros(lattice_shape)
count = 0
for i in range(0,lattice_shape[0]):
    for j in range(0,lattice_shape[1]):
        if resources[i,j] > 9:
            resources_filtered[i,j] = resources[i,j]
            count = count + 1

currentposition = initialbird_position

output_data_store = np.zeros((100,2))


#Define function to calculate the real distance between two given lattice points
def realdistance(a,b):
    latlonga   = ((lattice_shape[0]/2-a[0]-0.5)*d_latlong, (a[1]+0.5-lattice_shape[1]/2)*d_latlong) #(latitude,longitude) position of lattice point a
    latlongb   = ((lattice_shape[0]/2-b[0]-0.5)*d_latlong, (b[1]+0.5-lattice_shape[1]/2)*d_latlong) #(latitude,longitude) position of lattice point b
    delta_long = latlongb[1]-latlonga[1]
    term1      = sin(radians(latlonga[0]))*sin(radians(latlongb[0]))
    term2      = cos(radians(latlonga[0]))*cos(radians(latlongb[0]))*cos(radians(delta_long))
    #When the bird doesn't move, and latlonga = latlongb, small rounding errors can lead to taking the arccos of a number a tiny bit higher than 1, eg 1.0000000000000002, so it's safer to set the distance equal to 0 in this case rather than doing the full calculation
    if latlonga == latlongb:
        dist = 0.0
    else:
        try:
            dist = 6371*acos(term1+term2) # 6371 is the radius of the earth in km (assuming spherical)
        except Exception as e:
            print(latlonga,latlongb,term1,term2)
    return dist



for t in range(0,100):
    print(t)

    #Calculate potentials in new states
    #Convert to Boltzmann factors
    possible_state_boltzmann_factors = np.zeros((3,3))
    for i in range(-1,2):
        for j in range(-1,2):
            if i == 0 and j == 0:
                pass
            else:
                state_potential = 0
                state_index     = ((currentposition[0]+i),(currentposition[1]+j))
                for k in range(0,lattice_shape[0]):
                    for l in range(0,lattice_shape[1]):
                        if 99999 > resources_filtered[k,l] > 0:
                            state_potential = state_potential + resources_filtered[k,l]/realdistance(state_index,currentposition)
            possible_state_boltzmann_factors[i+1,j+1] = exp(-state_potential/kT)

    #Update position
    #Sum Boltzmann factors for possible states
    boltzmann_sum = 0
    for i in range(0,3):
        for j in range(0,3):
            boltzmann_sum = boltzmann_sum + possible_state_boltzmann_factors[i,j]

    #Use a random number generator and probabilities defined by Boltzmann factors
    #to decide which lattice point the bird moves to next
    probability_sum = 0
    random_number = boltzmann_sum*random()
    for i in range(0,3):
        for j in range(0,3):
            probability_sum = probability_sum + possible_state_boltzmann_factors[i,j]
            if random_number < probability_sum:
                currentposition = tuple(map(sum, zip(currentposition, (i-1,j-1))))
                break
            else:
                pass

    output_data_store[t,0] = (lattice_shape[0]/2-currentposition[0]-0.5)*d_latlong
    output_data_store[t,1] = (currentposition[1]+0.5-lattice_shape[1]/2)*d_latlong


np.savetxt("testdata.txt",output_data_store)

map = Basemap(projection="hammer",lon_0=0)
map.fillcontinents()
lats = output_data_store[:,0]
lons = output_data_store[:,1]
x,y=map(lats,lons)
map.plot(x,y)
plt.savefig("test.pdf")
