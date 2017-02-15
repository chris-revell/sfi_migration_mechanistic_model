#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3
from sys import argv,exit
import numpy as np
from random import random
from math import sqrt
import os
import time
import fortran_subroutines

#Find all datafiles
chloro_datafiles = [os.path.join("data_chloro",f) for f in os.listdir("data_chloro") if f[-4:].lower()==".csv"]
wind_merid_datafiles = [os.path.join("data_wind",f) for f in os.listdir("data_wind") if f[-4:].lower()==".txt" and "merid" in f]
wind_zonal_datafiles = [os.path.join("data_wind",f) for f in os.listdir("data_wind") if f[-4:].lower()==".txt" and "zonal" in f]

#Initial conditions
initial_lat = float(argv[1])
initial_lon = float(argv[2])
a           = float(argv[3]) # Relative contribution of wind potential to total potential where contribution of chlorophyll is 1
kT          = float(argv[4]) # Effective temperature of bird - high value increases randomness in path, lower value collapses to lowest energy state.
start_month = int(argv[5])   # Start month defines the end of the breeding season when birds leave the breeding colony, whose position is defined by initial_lat and initial_lon

bird_speed  = 60

if len(argv) <= 6:
    t_max = 30*24*12 #Full year
else:
    end_month   = int(argv[6])              # End month specified if we want to stop simulation before a full year.
    if end_month < start_month:
        end_month = end_month+12
    t_max=(end_month-start_month+1)*30*24   #Number of hours in months specified

#Function to rotate list of datafiles to begin at start_month
def rotate(l, n):
    return l[n-1:] + l[:n-1]

#Rotate datafile lists to begin with start_month
chloro_datafiles = rotate(chloro_datafiles,start_month)
wind_merid_datafiles = rotate(wind_merid_datafiles,start_month)
wind_zonal_datafiles = rotate(wind_zonal_datafiles,start_month)

#Import earth map - earth[i,j]=1 for dry land, 0 for ocean. Used to restrict movement to ocean and mask chlorophyll from inland lakes.
earth = np.asfortranarray(np.genfromtxt("earth.txt",delimiter=" "))
resources_shape = np.shape(earth)

d_latlong = 180/resources_shape[0] #Change in latitude and longitude angle per lattice point.

#Define function to convert lattice positions to latitude and longitude values.
def xytolatlong(xy):
    lat = (resources_shape[0]/2.0-xy[0]-0.5)*d_latlong
    lon = (xy[1]+0.5-resources_shape[1]/2.0)*d_latlong
    return (lat,lon)

#Initialise system time and position.
t=0.00001   #This needs to be close to but not equal to 0 for import loop later
dt = 0.001  #As for t, initial value of dt is just an abritrary small value and will be recalculated during simulations.
initialposition = (int(resources_shape[0]/2-initial_lat/d_latlong),int(resources_shape[1]/2+initial_lon/d_latlong))
currentposition = initialposition

#Apply a small random jiggle to the initial position
#jiggle = (int(2*random())-1,int(2*random())-1)
#currentposition = (initialposition[0]+int((-1)**jiggle[0]),initialposition[1]+int((-1)**jiggle[1]))

#Test that initial position is not trapped on dry land
onearth = 1
for i in range(3):
    for j in range(3):
        if i == 1 and j ==1:
            pass
        else:
            if earth[initialposition[0]-1+i,initialposition[1]-1+j] == 0:
                onearth = 0
            else:
                pass
if onearth == 1:
    exit("Error: Initial position is on dry land with no surrounding water - pick a different initial position\n")


#Create general data folder if it does not already exist
if os.path.exists("../output_data"):
    pass
else:
    os.mkdir("../output_data")
#Create folder for this particular run. If this set of parameters has been run before, program loops over run number until if finds a value that has not yet been used.
run_number = 0
run_folder = "../output_data/lat{:03.1f}_lon{:04.1f}_a{}_kT{}_m{:02d}-{:02d}_run{:02d}".format(initial_lat,initial_lon,argv[3],argv[4],start_month,end_month,run_number)
while os.path.exists(run_folder):
    run_number = run_number + 1
    run_folder = "../output_data/lat{:03.1f}_lon{:04.1f}_a{}_kT{}_m{:02d}-{:02d}_run{:02d}".format(initial_lat,initial_lon,argv[3],argv[4],start_month,end_month,run_number)
os.mkdir(run_folder)

#Save run parameters
parameterfile = open(os.path.join(run_folder,"parameters.txt"),'w')
parameterfile.write("a  = "+str(a)+"\n")
parameterfile.write("kt = "+str(kT)+"\n")
parameterfile.write("initialposition = "+argv[1]+", "+argv[2])
parameterfile.write("\nstart_month = "+str(start_month)+"\nend_month = "+str(end_month)+"\nbird_speed = "+str(bird_speed)+"\n")
parameterfile.write("t_max = "+str(t_max)+"\n")
parameterfile.close()
#Save positions for t=0
output_latlong_file = open(os.path.join(run_folder,"latlongdata.csv"),'w')
currentlatlon = xytolatlong(currentposition)
output_latlong_file.write(str(t)+","+str(currentlatlon[0])+","+str(currentlatlon[1])+"\n")

while t < t_max:

    previousposition = currentposition

    if t%720 < dt:
        #Refresh data arrays every 720 hours (~1 month)
        #Import chloro data
        chloro_filename = chloro_datafiles.pop(0)
        print(chloro_filename)
        resources = np.genfromtxt(chloro_filename,delimiter=",")
        resources_shape = np.shape(resources)

        #Use ocean-earth map to exclude inland lakes from chlorophyll data
        resources_filtered = np.asfortranarray(np.zeros(resources_shape))
        for i in range(0,resources_shape[0]):
            for j in range(0,resources_shape[1]):
                if earth[i,j] == 0:
                    resources_filtered[i,j] = resources[i,j]
                else:
                    pass

        #Import wind data
        merid_filename = wind_merid_datafiles.pop(0)
        print(merid_filename)
        wind_merid = np.asfortranarray(np.genfromtxt(merid_filename)) #North to south wind speed
        zonal_filename = wind_zonal_datafiles.pop(0)
        print(zonal_filename)
        wind_zonal = np.asfortranarray(np.genfromtxt(zonal_filename)) #West to east wind speed

    #Calculate potentials in new possible states and convert to Boltzmann factors
    possible_state_boltzmann_factors = np.asfortranarray(np.zeros((3,3)))

    #Call boltzmanncalc subroutine from fortran_subroutines library to calculate boltzmann factors for possible states
    fortran_subroutines.boltzmanncalc(possible_state_boltzmann_factors,currentposition[0],currentposition[1],initialposition[0],initialposition[1],earth,wind_merid[currentposition],wind_zonal[currentposition],resources_filtered,a,0.0,0.0,kT,t)
    
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

    #Update time
    wind_vector = np.array([wind_merid[currentposition],wind_zonal[currentposition]]) #In form [y,x] for ease of translation to np arrays.
    #Calculate dx vector on earth, not on euclidean lattice
    theta = fortran_subroutines.realdistance(currentposition[0],currentposition[1],previousposition[0],currentposition[1])
    zeta = fortran_subroutines.realdistance(currentposition[0],currentposition[1],currentposition[0],previousposition[1])
    dx = np.array([theta,zeta])

    #Calculate speed including bird speed and wind component, then update time for moving between lattice points
    speed = bird_speed + np.dot(dx,wind_vector)/sqrt(np.dot(dx,dx))
    dt = fortran_subroutines.realdistance(currentposition[0],currentposition[1],previousposition[0],previousposition[1])/speed
    t = t + dt

    #Output data
    #print(t,currentposition)
    currentlatlon = xytolatlong(currentposition)
    output_latlong_file.write(str(t)+","+str(currentlatlon[0])+","+str(currentlatlon[1])+"\n")
