#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

# Christopher Revell, University of Cambridge, 2016

from sys import argv,exit
import numpy as np
from random import random
from math import sqrt,pi,sin,radians,cos
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

#Define functions to convert lattice positions to latitude and longitude values and vice versa
def xytolatlong(xy):
    lat = (resources_shape[0]/2.0-xy[0]-0.5)*d_latlong
    lon = (xy[1]+0.5-resources_shape[1]/2.0)*d_latlong
    return (lat,lon)
def latlongtoxy(lat,lon):
    x = int(lon/d_latlong)+resources_shape[1]/2.0
    y = resources_shape[0]/2.0-int(lat/d_latlong)
    return (int(y),int(x))

dt = 1

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
parameterfile.write("a  = {}\n""kt = {}\n""initialposition = {}, {}\nstart_month = {}\nend_month = {}\nbird_speed = {}\n""t_max = {}\n".format(str(a),str(kT),argv[1],argv[2],str(start_month),str(end_month),str(bird_speed),str(t_max)))
parameterfile.close()
#Save positions for t=0
output_latlong_file = open(os.path.join(run_folder,"latlongdata.csv"),'w')
currentlatlon = xytolatlong(migratingbird.position)
output_latlong_file.write("{:08.4f}, {:08.4f}, {:08.4f}\n".format(0.0,currentlatlon[0],currentlatlon[1]))

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

    #Calculate force
    force = np.asfortranarray(np.zeros((2),dtype='float32'))
    #Call fortran subroutine for resources component
    seabirdsubroutines.seabird_subroutines.forcecalc(force,currentposition[0],currentposition[1],resources_filtered,earth)#wind_merid_current,wind_zonal_current,resources_filtered,a)
    #Calculate wind component
    currentwindlatticeposition = latlongtoxy(currentposition[0],currentposition[1])
    wind_merid_current = wind_merid[currentwindlatticeposition]
    wind_zonal_current = wind_zonal[currentwindlatticeposition]
    windforce = np.asfortranarray(np.array([wind_merid_current,wind_zonal_current]))
    #Calculate stochastic component
    stochasticvals = [random(),2*pi*random()]
    stochasticforce = np.asfortranarray(np.array([stochasticvals[0]*sin(stochasticvals[1]),stochasticvals[0]*cos(stochasticvals[1])]))
    #Sum components
    force = force + a*windforce + kT*stochasticforce

    #Treat bird as particle moving with overdamped Langevin dynamics
    migratingbird.position = migratingbird.position + force*dt
    t = t + dt

    #Output data
    print(t,migratingbird.position[0],migratingbird.position[1])
    output_latlong_file.write('{:08.4f}, {:08.4f}, {:08.4f}\n'.format(t,migratingbird.position[0],migratingbird.position[1]))
