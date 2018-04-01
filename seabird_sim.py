#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

# Christopher Revell, University of Cambridge, 2016

from sys import argv,exit
import numpy as np
from random import random
from math import sqrt,pi,sin,radians
import os
import time
import seabirdsubroutines

#%%
#Initial conditions
initial_lat = float(argv[1])
initial_lon = float(argv[2])
dest_lat    = float(argv[3])
dest_lon    = float(argv[4])
a           = float(argv[5]) # Relative contribution of potential components
kT          = float(argv[6]) # Effective temperature of bird - high value increases randomness in path, lower value collapses to lowest energy state.
bird_speed  = 60
b=1

elevation = np.genfromtxt("input_data/elevation.txt")

t_max=6*30*24   #Number of hours in months specified

"""
#Function to rotate list of datafiles to begin at start_month
def rotate(l, n):
    return l[n-1:] + l[:n-1]

#Rotate datafile lists to begin with start_month
chloro_datafiles = rotate(chloro_datafiles,start_month)
wind_merid_datafiles = rotate(wind_merid_datafiles,start_month)
wind_zonal_datafiles = rotate(wind_zonal_datafiles,start_month)
"""
#Import earth map - earth[i,j]=1 for dry land, 0 for ocean. Used to restrict movement to ocean and mask chlorophyll from inland lakes.
earth = np.asfortranarray(np.genfromtxt("input_data/earth.txt",delimiter=" "))
resources_shape = np.shape(earth)
for i in range(resources_shape[0]):
    for j in range(resources_shape[1]):
        if earth[i,j] < 0.0001:
            earth[i,j] = 1
        else:
            earth[i,j] = 0
d_latlong = 180/resources_shape[0] #Change in latitude and longitude angle per lattice point.

#Define function to convert lattice positions to latitude and longitude values.
def xytolatlong(xy):
    lat = (resources_shape[0]/2.0-xy[0]-0.5)*d_latlong
    lon = (xy[1]+0.5-resources_shape[1]/2.0)*d_latlong
    return (lat,lon)

#Initialise system time and position.
t = 0.00001   #This needs to be close to but not equal to 0 for import loop later
dt = 0.001  #As for t, initial value of dt is just an abritrary small value and will be recalculated during simulations.
initialposition = (int(resources_shape[0]/2-initial_lat/d_latlong),int(resources_shape[1]/2+initial_lon/d_latlong))
destinationposition = (int(resources_shape[0]/2-dest_lat/d_latlong),int(resources_shape[1]/2+dest_lon/d_latlong))
currentposition = initialposition

"""
#Test that the initial position is not trapped on dry land
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
"""

#Create general data folder if it does not already exist
if os.path.exists("output_data"):
    pass
else:
    os.mkdir("output_data")
#Create folder for this particular run. If this set of parameters has been run before, program loops over run number until if finds a value that has not yet been used.
run_number = 0
run_folder = "output_data/lat{:04.1f}_lon{:04.1f}_lat{:04.1f}_lon{:04.1f}_a{:05.4f}_kT{:03.2f}_run{:02d}".format(initial_lat,initial_lon,dest_lat,dest_lon,a,kT,run_number)
while os.path.exists(run_folder):
    run_number = run_number + 1
    run_folder = "output_data/lat{:04.1f}_lon{:04.1f}_lat{:04.1f}_lon{:04.1f}_a{:05.4f}_kT{:03.2f}_run{:02d}".format(initial_lat,initial_lon,dest_lat,dest_lon,a,kT,run_number)
os.mkdir(run_folder)

#Save run parameters
parameterfile = open(os.path.join(run_folder,"parameters.txt"),'w')
parameterfile.write("initialposition = {},{}\ndestination = {},{}\na  = {}\nkt = {}\nbird_speed = {:04.2f}\n""t_max = {:05.2f}\n".format(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],bird_speed,t_max))
parameterfile.close()
#Save positions for t=0
output_latlong_file = open(os.path.join(run_folder,"latlongdata.csv"),'w')
output_latlong_file.write("{:08.4f}, {:08.4f}, {:08.4f}\n".format(0.0,initial_lat,initial_lon))

while t < t_max:

    previousposition = currentposition

    #Calculate potentials in new possible states and convert to Boltzmann factors
    possible_state_boltzmann_factors = np.asfortranarray(np.zeros((3,3)))

    elevationslice = np.asfortranarray(elevation[currentposition[0]-1:currentposition[0]+2,currentposition[1]-1:currentposition[1]+2])
    print(elevationslice)
    #Call boltzmanncalc subroutine from seabirdsubroutines library to calculate boltzmann factors for possible states
    seabirdsubroutines.boltzmanncalc(possible_state_boltzmann_factors,currentposition[0],currentposition[1],destinationposition[0],destinationposition[1],elevationslice,earth,a,b,kT)

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

    """
    #Update time
    wind_vector = np.array([wind_merid[currentposition],wind_zonal[currentposition]]) #In form [y,x] for ease of translation to np arrays.
    #Calculate dx vector on earth, not on euclidean lattice
    theta = seabirdsubroutines.realdistance(currentposition[0],currentposition[1],previousposition[0],currentposition[1])
    zeta = seabirdsubroutines.realdistance(currentposition[0],currentposition[1],currentposition[0],previousposition[1])
    dx = np.array([theta,zeta])

    #Calculate speed including bird speed and wind component, then update time for moving between lattice points
    speed = bird_speed + np.dot(dx,wind_vector)/sqrt(np.dot(dx,dx))
    dt = seabirdsubroutines.realdistance(currentposition[0],currentposition[1],previousposition[0],previousposition[1])/speed
    """
    t = t + 6

    #Output data
    print(t,currentposition)
    currentlatlon = xytolatlong(currentposition)
    output_latlong_file.write(str(t)+","+str(currentlatlon[0])+","+str(currentlatlon[1])+"\n")
