#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3
from sys import argv,exit
import numpy as np
from random import random
from math import exp, acos, sin, cos, pi, radians, sqrt
import os
import time
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

chloro_datafiles = [os.path.join("data_chloro",f) for f in os.listdir("data_chloro") if f[-4:].lower()==".csv"]
wind_merid_datafiles = [os.path.join("data_wind",f) for f in os.listdir("data_wind") if f[-4:].lower()==".csv" and f[0]=="m"]
wind_zonal_datafiles = [os.path.join("data_wind",f) for f in os.listdir("data_wind") if f[-4:].lower()==".csv" and f[0]=="z"]

initialposition = (576,573) #Position of South Georgia and the Sandwich Islands
a = float(argv[1])
kT = float(argv[2])

if len(argv) > 3:
    start_month = argv[3]
else:
    start_month = 1
if len(argv) > 4:
    end_month = argv[4]
else:
    end_month = 12

t_max=(end_month-start_month+1)*30*24   #Number of hours in months specified

bird_speed = 60

#Import ground map
earth = np.genfromtxt("earth1440x720.CSV",delimiter=",")
resources_shape = np.shape(earth)
for i in range(0,resources_shape[0]):
    for j in range(0,resources_shape[1]):
        if earth[i,j] == 99999:
            earth[i,j] = 0
        else:
            earth[i,j] = 1


d_latlong = 180/resources_shape[0]

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
    elif term1+term2 > 1:
        dist = 6371*acos(1)
    elif term1+term2 < -1:
        dist = 6371*acos(-1)
    else:
        dist = 6371*acos(term1+term2) # 6371 is the radius of the earth in km (assuming spherical)
    return dist

t=0
currentposition = initialposition
#Create folder
if os.path.exists("../output_data"):
    pass
else:
    os.mkdir("../output_data")
run_folder = os.path.join("../output_data/",time.strftime("%y%m%d%H%M")+"_a"+str(a))
os.mkdir(run_folder)
#Save run parameters
parameterfile = open(os.path.join(run_folder,"parameters.txt"),'w')
parameterfile.write("a = "+str(a)+"\nkt = "+str(kT)+"\nt_max = "+str(t_max)+"\ninitialposition = "+str(initialposition))
parameterfile.write("\nstart_month = "+str(start_month)+"\nend_month = "+str(end_month)+"\nbird_speed = "+str(bird_speed))
parameterfile.close()
output_data_file = open(os.path.join(run_folder,"positiondata.csv"),'w')
output_data_file.write(str(t)+","+str(currentposition[0])+","+str(currentposition[1])+"\n")

t=0.00001
dt = 0.001
previous_position = initialposition
while t < t_max:

    if t%720 < dt:
        #Refresh data arrays every 720 hours (~1 month)
        #Import chloro data
        resources = np.genfromtxt(chloro_datafiles.pop(0),delimiter=",")
        resources_shape = np.shape(resources)

        #Threshold chloro data
        resources_filtered = np.zeros(resources_shape)
        for i in range(0,resources_shape[0]):
            for j in range(0,resources_shape[1]):
                if 99999 > resources[i,j] > 5:
                    resources_filtered[i,j] = resources[i,j]

        #Import wind data
        wind_merid = np.genfromtxt(wind_merid_datafiles.pop(0),delimiter=",") #North to south wind speed
        wind_zonal = np.genfromtxt(wind_zonal_datafiles.pop(0),delimiter=",") #West to east wind speed

        refresh_timer = 0


    #Calculate potentials in new possible states and convert to Boltzmann factors
    possible_state_boltzmann_factors = np.zeros((3,3))
    for i in range(-1,2):
        for j in range(-1,2):
            state_index = ((currentposition[0]+i),(currentposition[1]+j)%resources_shape[1])
            if i == 0 and j == 0:
                pass
            elif earth[state_index] == 1:
                pass
            else:
                state_potential = 0
                for k in range(0,resources_shape[0]):
                    for l in range(0,resources_shape[1]):
                        if earth[k,l] == 0 and resources_filtered[k,l] > 0 and (k,l) != state_index:
                            state_potential = state_potential + resources_filtered[k,l]/realdistance((k,l),state_index)

                wind_vector = np.array([wind_merid[currentposition],wind_zonal[currentposition]]) #In form [y,x] for ease of translation to np arrays.
                wind_magnitude = sqrt(np.dot(wind_vector,wind_vector))
                displacement_vector = np.array([i,j])
                displacement_vector_magnitude = sqrt(np.dot(displacement_vector,displacement_vector))
                state_potential = state_potential + a*wind_magnitude*np.dot(wind_vector,displacement_vector)/displacement_vector_magnitude

                possible_state_boltzmann_factors[i+1,j+1] = exp(state_potential/kT)

    #Update position
    #Sum Boltzmann factors for possible states
    boltzmann_sum = 0
    for i in range(0,3):
        for j in range(0,3):
            boltzmann_sum = boltzmann_sum + possible_state_boltzmann_factors[i,j]
    #Use a random number generator and probabilities defined by Boltzmann factors to decide which lattice point the bird moves to next
    probability_sum = 0
    random_number = boltzmann_sum*random()
    #Use the moved variable to ensure that the bird can only move once, otherwise the loop below cannot be exited cleanly and the bird will be moved multiple times.
    moved = 0
    for i in range(0,3):
        for j in range(0,3):
            probability_sum = probability_sum + possible_state_boltzmann_factors[i,j]
            if moved == 0:
                if random_number < probability_sum:
                    currentposition = (currentposition[0]+i-1,(currentposition[1]+j-1)%resources_shape[1]) # Use of mod % allows birds to move off one side of the grid and appear at the other. Ignore north and south poles for now because birds should never reach this point.
                    moved = 1
                else:
                    pass
            else:
                pass

    dt = realdistance(currentposition,previous_position)
    t = t + dt/bird_speed


    print(t,currentposition)
    output_data_file.write(str(t)+","+str(currentposition[0])+","+str(currentposition[1])+"\n")



"""
#Save position data
np.savetxt(os.path.join(run_folder,"positiondata.txt"),output_data_store)

#Create map of bird path on basemap
map = Basemap(projection="robin",lon_0=0)
map.fillcontinents(color='coral',lake_color='aqua')
map.drawmapboundary(fill_color='aqua')
lats = (resources_shape[0]/2-output_data_store[:,0]-0.5)*d_latlong
lons = (output_data_store[:,1]+0.5-resources_shape[1]/2)*d_latlong
x,y=map(lons,lats)
map.plot(x,y)
plt.savefig(os.path.join(run_folder,"map.pdf"))

#Create plot of bird path on lattice
fig2 = plt.figure()
ax2 = fig2.add_subplot(211)
cax = ax2.imshow(resources_filtered,cmap="viridis")
cbar = fig2.colorbar(cax)
ax2.plot(output_data_store[:,1],output_data_store[:,0])
ax2.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
ax2.set_xlim([0,resources_shape[1]])
ax2.set_ylim([resources_shape[0],0])
ax3 = fig2.add_subplot(212)
ax3.imshow(earth,cmap="Greys")
ax3.plot(output_data_store[:,1],output_data_store[:,0])
ax3.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
ax3.set_xlim([0,resources_shape[1]])
ax3.set_ylim([resources_shape[0],0])
fig2.savefig(os.path.join(run_folder,"lattice.pdf"))
"""
