#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

# Christopher Revell, University of Cambridge, 2016

from sys import argv,exit
import numpy as np
from random import random
from math import sqrt,pi,sin,radians
import os
import time
import seabirdsubroutines

#Find all datafiles
chloro_datafiles = ["data_chloro/chloromean01.csv", "data_chloro/chloromean02.csv", "data_chloro/chloromean03.csv", "data_chloro/chloromean04.csv", "data_chloro/chloromean05.csv", "data_chloro/chloromean06.csv", "data_chloro/chloromean07.csv", "data_chloro/chloromean08.csv", "data_chloro/chloromean09.csv", "data_chloro/chloromean10.csv", "data_chloro/chloromean11.csv", "data_chloro/chloromean12.csv"]
wind_merid_datafiles = ["data_wind/merid_mean01.txt", "data_wind/merid_mean02.txt", "data_wind/merid_mean03.txt", "data_wind/merid_mean04.txt", "data_wind/merid_mean05.txt", "data_wind/merid_mean06.txt", "data_wind/merid_mean07.txt", "data_wind/merid_mean08.txt", "data_wind/merid_mean09.txt", "data_wind/merid_mean10.txt", "data_wind/merid_mean11.txt", "data_wind/merid_mean12.txt"]
wind_zonal_datafiles = ["data_wind/zonal_mean01.txt", "data_wind/zonal_mean02.txt", "data_wind/zonal_mean03.txt", "data_wind/zonal_mean04.txt", "data_wind/zonal_mean05.txt", "data_wind/zonal_mean06.txt", "data_wind/zonal_mean07.txt", "data_wind/zonal_mean08.txt", "data_wind/zonal_mean09.txt", "data_wind/zonal_mean10.txt", "data_wind/zonal_mean11.txt", "data_wind/zonal_mean12.txt"]

class bird:
    def __init__(self,position):
        self.position = position
    velocity = [0,0]

#Initial conditions
initial_lat = float(argv[1])
initial_lon = float(argv[2])
a           = float(argv[3]) # Relative contribution of wind potential to total potential where contribution of chlorophyll is 1
kT          = float(argv[4]) # Effective temperature of bird - high value increases randomness in path, lower value collapses to lowest energy state.
start_month = int(argv[5])   # Start month defines the end of the breeding season when birds leave the breeding colony, whose position is defined by initial_lat and initial_lon

bird_speed  = 60

migratingbird = bird([initial_lat,initial_lon])


if len(argv) <= 6:
    t_max = 30*24*12 #Full year
    end_month = start_month-1
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
#
##Define function to convert lattice positions to latitude and longitude values.
def xytolatlong(xy):
    lat = (resources_shape[0]/2.0-xy[0]-0.5)*d_latlong
    lon = (xy[1]+0.5-resources_shape[1]/2.0)*d_latlong
    return (lat,lon)

#Initialise system time and position.
#t=0.00001   #This needs to be close to but not equal to 0 for import loop later
dt = 1  #As for t, initial value of dt is just an abritrary small value and will be recalculated during simulations.
#initialposition = (int(resources_shape[0]/2-initial_lat/d_latlong),int(resources_shape[1]/2+initial_lon/d_latlong))
#currentposition = initialposition

#Test that the initial position is not trapped on dry land
#onearth = 1
#for i in range(3):
#    for j in range(3):
#        if i == 1 and j ==1:
#            pass
#        else:
#            if earth[initialposition[0]-1+i,initialposition[1]-1+j] == 0:
#                onearth = 0
#            else:
#                pass
#if onearth == 1:
#    exit("Error: Initial position is on dry land with no surrounding water - pick a different initial position\n")


#Create general data folder if it does not already exist
if os.path.exists("output_data"):
    pass
else:
    os.mkdir("output_data")
#Create folder for this particular run. If this set of parameters has been run before, program loops over run number until if finds a value that has not yet been used.
run_number = 0
run_folder = "output_data/lat{:03.1f}_lon{:04.1f}_a{:05.4f}_kT{:03.2f}_m{:02d}-{:02d}_run{:02d}".format(initial_lat,initial_lon,a,kT,start_month,end_month,run_number)
while os.path.exists(run_folder):
    run_number = run_number + 1
    run_folder = "output_data/lat{:03.1f}_lon{:04.1f}_a{}_kT{}_m{:02d}-{:02d}_run{:02d}".format(initial_lat,initial_lon,argv[3],argv[4],start_month,end_month,run_number)
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
currentlatlon = xytolatlong(migratingbird.position)
output_latlong_file.write(str(0.0)+","+str(currentlatlon[0])+","+str(currentlatlon[1])+"\n")

t=0
while t < t_max:

    currentposition = migratingbird.position

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
#            lat1 = radians(xytolatlong([i,1])[0]-d_latlong/2)
#            lat2 = lat1 + radians(d_latlong)
#            area = (6371.0**2)*abs(sin(lat1)-sin(lat2))*abs(radians(d_latlong))
#            print(i,area)
            for j in range(0,resources_shape[1]):
                if earth[i,j] == 0:
                    resources_filtered[i,j] = resources[i,j]#*area
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
    #possible_state_boltzmann_factors = np.asfortranarray(np.zeros((3,3)))

    wind_merid_current = wind_merid[int(resources_shape[0]/2-migratingbird.position[0]/d_latlong),int(resources_shape[1]/2+migratingbird.position[1]/d_latlong)]
    wind_merid_current = wind_merid[int(resources_shape[0]/2-migratingbird.position[0]/d_latlong),int(resources_shape[1]/2+migratingbird.position[1]/d_latlong)]

    #Call boltzmanncalc subroutine from seabird_subroutines library to calculate boltzmann factors for possible states
    force = np.asfortranarray(np.zeros((2),dtype='float32'))
    seabirdsubroutines.seabird_subroutines.forcecalc(force,currentposition[0],currentposition[1],resources_filtered)#wind_merid_current,wind_zonal_current,resources_filtered,a)


    #Update time. Calculate speed including bird speed and wind component, then update time for moving between lattice points
    #Calculate dx vector on earth, not on euclidean lattice
    #wind_vector = np.array([wind_merid[currentposition],wind_zonal[currentposition]]) #In form [y,x] for ease of translation to np arrays.
    #theta = seabirdsubroutines.seabird_subroutines.realdistance(currentposition[0],currentposition[1],previousposition[0],currentposition[1])
    #zeta = seabirdsubroutines.seabird_subroutines.realdistance(currentposition[0],currentposition[1],currentposition[0],previousposition[1])
    #dx = np.array([theta,zeta])
    #speed = bird_speed + np.dot(dx,wind_vector)/sqrt(np.dot(dx,dx))

    migratingbird.position = migratingbird.position + force*dt
    t = t + dt

    #Output data
    print(t,migratingbird.position[0],migratingbird.position[1])
    output_latlong_file.write('{:08.4f}, {:08.4f}, {:08.4f}\n'.format(t,migratingbird.position[0],migratingbird.position[1]))
