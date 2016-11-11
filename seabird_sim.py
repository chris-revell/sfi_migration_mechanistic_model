
from sys import argv,exit
import numpy as np
from random import random
from math import exp, acos, sin, cos, pi, radians, sqrt
import os
import time
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

datafiles = [os.path.join(argv[1],f) for f in os.listdir(argv[1]) if f[-4:].lower()==".csv"]

currentposition = (500,720)
kT = 100
t_max=200
A = 1

#Import ground map
earth = np.genfromtxt(datafiles[0],delimiter=",")
earth_shape = np.shape(earth)
for i in range(0,earth_shape[0]):
    for j in range(0,earth_shape[1]):
        if earth[i,j] == 99999:
            earth[i,j] = 0
        else:
            earth[i,j] = 1


#Import chloro data
resources = np.genfromtxt(datafiles[1],delimiter=",")
resources_shape = np.shape(resources)

#Import wind data
wind_merid = np.genfromtxt(datafiles[2],delimiter=",") #North to south wind speed
wind_zonal = np.genfromtxt(datafiles[3],delimiter=",") #West to east wind speed

d_latlong = 180/resources_shape[0]

#Threshold chloro data
resources_filtered = np.zeros(resources_shape)
for i in range(0,resources_shape[0]):
    for j in range(0,resources_shape[1]):
        if 99999 > resources[i,j] > 5:
            resources_filtered[i,j] = resources[i,j]

output_data_store = np.zeros((t_max,2))

#Define function to calculate the real distance between two given lattice points
def realdistance(a,b):
    latlonga   = ((resources_shape[0]/2-a[0]-0.5)*d_latlong, (a[1]+0.5-resources_shape[1]/2)*d_latlong) #(latitude,longitude) position of lattice point a
    latlongb   = ((resources_shape[0]/2-b[0]-0.5)*d_latlong, (b[1]+0.5-resources_shape[1]/2)*d_latlong) #(latitude,longitude) position of lattice point b
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


output_data_store[0,0] = currentposition[0]
output_data_store[0,1] = currentposition[1]


for t in range(1,t_max):
    print(t)
    #Calculate potentials in new states
    #Convert to Boltzmann factors
    possible_state_boltzmann_factors = np.zeros((3,3))
    for i in range(-1,2):
        for j in range(-1,2):
            state_index     = ((currentposition[0]+i),(currentposition[1]+j))
            if i == 0 and j == 0:
                pass
            elif earth[state_index] == 1:
                pass
            else:
                state_potential = 0
                for k in range(0,resources_shape[0]):
                    for l in range(0,resources_shape[1]):
                        if earth[k,l] == 0 and resources_filtered[k,l] > 0:
                            state_potential = state_potential + resources_filtered[k,l]/realdistance(state_index,currentposition)


                wind_vector = np.array([wind_merid[currentposition],wind_zonal[currentposition]]) #In form [y,x] for ease of translation to np arrays.
                wind_magnitude = sqrt(np.dot(wind_vector,wind_vector))
                displacement_vector = np.array([i,j])

                state_potential = state_potential + A*wind_magnitude*np.dot(wind_vector,displacement_vector)

                possible_state_boltzmann_factors[i+1,j+1] = exp(state_potential/kT)

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

    output_data_store[t,0] = currentposition[0]
    output_data_store[t,1] = currentposition[1]


#Output data

if os.path.exists("../output_data"):
    pass
else:
    os.mkdir("../output_data")
run_folder = os.path.join("../output_data/",time.strftime("%y%m%d%H%M")+"_A"+str(A))
os.mkdir(run_folder)


np.savetxt(os.path.join(run_folder,"positiondata.txt"),output_data_store)

map = Basemap(projection="hammer",lon_0=0)
map.fillcontinents()
lats = (resources_shape[0]/2-output_data_store[:,0]-0.5)*d_latlong
lons = (output_data_store[:,1]+0.5-resources_shape[1]/2)*d_latlong
x,y=map(lats,lons)
map.plot(x,y)
plt.savefig(os.path.join(run_folder,"map.pdf"))

fig2 = plt.figure()
ax2 = fig2.add_subplot(211)
cax = ax2.imshow(resources_filtered,cmap="GnBu")
cbar = fig2.colorbar(cax)
ax2.plot(output_data_store[:,1],output_data_store[:,0])
ax2.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='off')
ax2.set_xrange([0,resources_shape[1]])
ax2.set_yrange([resources_shape[0],0])
ax2.set_yrange
fig2.savefig(os.path.join(run_folder,"lattice_chloro.pdf"))

fig3 = plt.figure()
ax3 = fig3.add_subplot(212)
ax3.imshow(earth,cmap="Greys")
ax3.plot(output_data_store[:,1],output_data_store[:,0])
ax3.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')
fig3.savefig(os.path.join(run_folder,"lattice_earth.pdf"))
