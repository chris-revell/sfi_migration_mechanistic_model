
from sys import argv,exit
import numpy as np
import random
from math import exp, acos, sin, cos, pi
import os
import time
import matplotlib.pyplot as pyplot

if len(argv) != 3:
    exit("Please provide 2 arguments: folder path of environment data and the value of kT")

#Folder path for datafile given at command line
importfolderpath = argv[1]

#Set system parameters
#A       = float(argv[2])# Prefactor for breeding site gravitational attraction. Given at command line. ~10^5
kT      = float(argv[2])# Measure of bird temperature or "restlessness". Given at command line. ~10^3
n_runs  = 5             # Number of runs with this set of parameters. Program will produce an average and standard deviation over all runs.
n_output= 1000          # Number of data outputs to file
#Define position of breeding ground and initial position of bird
#breeding_position    = (279,1147)
initialbird_position = (1500,0)
bird_speed           = 100    #Average flight speed of birds in km/h

#Add some parameters of the chlorophyll concentration lattice
d_latlong         = 2*pi/3600 #Change in latitude and longitude values between neighbouring lattice points. Same value for latitude and longitude.
origin_position   = (0,-pi)   #latitude,longitude values for lattice point (0,0)
radius_earth      = 6371      #Radius of the earth in km (assumed spherical for simplicity)

#From folder path provided at command line, find list of files to import chloro data from.
#Each file corresponds to half a month. "isfile" checks that we find only files, not directories.
#f[-1]=='t' excludes anything but .txt files (specifically to exclude .ds_store file on Mac)
#Note that list will be ordered alphabetically, so alphabetical order of filenames must match temporal order of months.
datafiles = [f for f in os.listdir(importfolderpath) if os.path.isfile(os.path.join(importfolderpath, f)) and f[-4:]=='.CSV']
datafiles.sort()
print(datafiles)

#If each datafile is half a month, this corresponds to roughly 15 days per file, which is 360 hours.
#Thus set the total run time in hours from the total number of input files minus 1, multiplied by 360.
#Minus 1 required because the system runs on the interpolated spaces between files.
t_max = 360*(len(datafiles)-1)  # Total time in hours for simulation
#Set output interval
output_interval = t_max/n_output #output_interval ~ 1 hour

#Create array to store data from all runs
output_data_store = np.empty((n_output,(n_runs*3)),dtype=float)

#Create date and time labelled folder to store the data from this run
#First ensure the output_data folder exists
if os.path.exists("../output_data"):
    pass
else:
    os.mkdir("../output_data")
run_folder = '../output_data/'+time.strftime("%y%m%d%H%M")+'_'+'kT'+argv[2]
os.mkdir(run_folder)

#Set interval for importing new datafiles. Note -1 because the system ends once interpolation reaches the state of the final file.
update_interval = t_max/(len(datafiles)-1)

#Use numpy.genfromtxt to import matrix in .txt file into array.
#Data is imported in float format.
#Higher chlorophyll concentration gives higher values in data. We want the potential to be lower at greener areas, so we multiply every value by -1 when using as a potential.
#Import initial chlorophyll file
initialfilename = os.path.join(importfolderpath,datafiles[0])
chloro_next = np.genfromtxt(initialfilename, dtype=float, delimiter=',')

#Find the shape of the environment array
arrayshape = np.shape(chloro_next)

#Use the shape of the imported chlorophyll concentration array to define the shape of other system arrays.
boltzmann_factors   = np.zeros(arrayshape,dtype=float)
chloro_gradient     = np.empty(arrayshape,dtype=float)
chloro_interpolated = np.empty(arrayshape,dtype=float)
time_updated        = np.zeros(arrayshape,dtype=float)

#Imported data ranges from 0 to 1 apart from regions in which it was not possible to collect data, such as on land.
#Such areas are given values of 999.0, but for the purposes of a potential we need to change these values to 0
for x in range(0,arrayshape[0]):
    for y in range(0,arrayshape[1]):
        if chloro_next[x,y] > 1.0:
            chloro_next[x,y] = 0.0
        else:
            pass

#Define list of possible lattice points that the bird can be in at the next timestep.
possible_states = []

#Define subroutine to find the set of possible lattice points that the bird can be in at the next timestep.
def find_possible_states():
    global possible_states

    #Start by removing all elements in previous possible states list in case the number of elements in the new list is smaller
    del possible_states[:]

    #Identify current neighbouring elements
    #Boundary conditions allow birds to move off one edge of the array and onto the opposite edge.
    for x in range(-1,2):
        for y in range(-1,2):
            possible_states.append((bird_position[0]+x,bird_position[1]+y))

#Define function to calculate the real distance between two given lattice points
def realdistance(a,b):
    latlonga       = (pi*(origin_position[0]-a[0]*d_latlong)/180, pi*(origin_position[1]+a[1]*d_latlong)/180) #(latitude,longitude) position of lattice point a, in radians
    latlongb       = (pi*(origin_position[0]-b[0]*d_latlong)/180, pi*(origin_position[1]+b[1]*d_latlong)/180) #(latitude,longitude) position of lattice point b, in radians
    delta_long     = latlongb[1]-latlonga[1]
    term1          = sin(latlonga[0])*sin(latlongb[0])
    term2          = cos(latlonga[0])*cos(latlongb[0])*cos(delta_long)
    #When the bird doesn't move, and latlonga = latlongb, small rounding errors can lead to taking the arccos of a number a tiny bit higher than 1, eg 1.0000000000000002, so it's safer to set the distance equal to 0 in this case rather than doing the full calculation
    if latlonga == latlongb:
        dist = 0.0
    else:
        dist = radius_earth*acos(term1+term2)
    return dist

#Subroutine to update the boltzmann factors for the possible states in the next timestep after the chloro values of these states have been interpolated.
def boltzmann_update(possible_states):
    global boltzmann_factors
    global chloro_interpolated
    for element in possible_states:
        #Define potential at (x,y) from chloro data
        #Note that removing factor of -1 in exponential is equivalent to multiplying chloro values by -1 to form a potential
        boltzmann_factors[element] = exp(chloro_interpolated[element]/kT)

#Define subroutine to update the position of the bird
def system_update(t):
    global boltzmann_factors
    global bird_position
    global possible_states
    global output_data_store

    #Sum Boltzmann factors for possible states
    boltzmann_sum = 0
    for z in possible_states:
        boltzmann_sum = boltzmann_sum+boltzmann_factors[z]

    #Use a random number generator and probabilities defined by Boltzmann factors
    #to decide which lattice point the bird moves to next
    probability_sum = 0
    random_number = boltzmann_sum*random.random()
    for p in possible_states:
        probability_sum = probability_sum + boltzmann_factors[p[0],p[1]]
        if random_number < probability_sum:
            bird_position = [p[0],p[1]]
            break
        else:
            pass

#Subroutine to import a new chloro file and redefine the previous "next" file as the new "current" file.
#Calculates the gradient array between the two files.
def importnext(t,d):
    global chloro_interpolated
    global chloro_gradient
    global chloro_next
    global time_updated
    #Set the current chloro dataset to the previous "next" dataset
    chloro_interpolated = chloro_next
    #Use datafiles list to find the filename of the next file to import chloro data from
    filenameatinterval = os.path.join(importfolderpath,datafiles[d])
    #Use this filename to import data from file into array
    chloro_next = np.genfromtxt(filenameatinterval, dtype=float, delimiter=',')
    for x in range(0,arrayshape[0]):
        for y in range(0,arrayshape[1]):
            if chloro_next[x,y] > 1.0:
                chloro_next[x,y] = 0.0

    #Calculate gradient for interpolation by broadcasting subtraction and then division by interval over all components in data arrays.
    chloro_gradient = (chloro_next-chloro_interpolated)/update_interval
    print (t, filenameatinterval)
    #Change the value of the time the element was last updated for every element of the array.
    time_updated.fill(t)

#Function to interpolate chloro data between two files
def interpolate(possible_states,t):
    global chloro_gradient
    global chloro_interpolated
    global time_updated
    #Loop over all elements in array
    #Use difference between current time t and the last time an element was updated to determine how many multiples of the gradient should be added to interpolate that element.
    for element in possible_states:
        chloro_interpolated[element] = chloro_interpolated[element] + (t-time_updated[element])*chloro_gradient[x,y] #(t-time_updated[element]) gives the time that has elapsed since that component of the array was last interpolated, and therefore the magnitude of the prefactor required when adding some multiple of the gradient
        #Element updated, so set the time_updated to the current time t.
        time_updated[element] = t

for i in range (0,n_runs):
    print('run '+str(i))
    #Reset system for each new run
    bird_position      = initialbird_position
    prev_bird_position = bird_position
    chloro_next = np.genfromtxt(initialfilename, dtype=float, delimiter=',')
    for x in range(0,arrayshape[0]):
        for y in range(0,arrayshape[1]):
            if chloro_next[x,y] > 1.0:
                chloro_next[x,y] = 0.0

    #Find possible states for first run of system
    find_possible_states()
    boltzmann_update(possible_states)
    #Initialise time and counter for update interval
    t = 0.0
    counter = 0
    #Iteratively update system for full time period
    while t<t_max:
        #For every import interval, import a new file into chloro_next and redefine chloro_interpolated to hold the old values of chloro_next
        if (t>=counter*update_interval):
            counter = counter + 1
            importnext(t,counter)
        #Update system state according to current interpolated chloro values and corresponding Boltzmann factors.
        prev_bird_position = bird_position #Store the previous bird position before updating
        system_update(t)
        timetesttuple = divmod(t,output_interval)
        if int(timetesttuple[1]) == 0:
            #Output data to storage array at every output interval
            output_data_store[int(timetesttuple[0]),i*3]   = t
            output_data_store[int(timetesttuple[0]),i*3+1] = bird_position[0]
            output_data_store[int(timetesttuple[0]),i*3+2] = bird_position[1]

        distance_travelled = realdistance(bird_position,prev_bird_position)
        if distance_travelled == 0:
            t=t+0.5
        else:
            t = t + distance_travelled/bird_speed #Time taken for the bird to travel this distance between lattice points.
        #Find possible states for next run of system
        find_possible_states()
        #Update the interpolated chloro array
        interpolate(possible_states,t)
        #Update Boltzmann factors according to the new interpolated chloro values
        boltzmann_update(possible_states)


#Write stored data array to file
np.savetxt(run_folder+'/bird_positions.txt', output_data_store, delimiter='  ')
outfile2 = open(run_folder+'/distance_average.txt','w')

#Calculate centre of mass position for all runs at each time point
COM_array = np.zeros((n_output,2))
for i in range(0,n_runs):
    COM_array[:,0] = COM_array[:,0] + output_data_store[:,3*i+1]
    COM_array[:,1] = COM_array[:,1] + output_data_store[:,3*i+2]
COM_array = COM_array/n_runs
np.savetxt(run_folder+'/COM_path.txt', COM_array, delimiter='  ')

#Plot data with pyplot
#Plot paths
pyplot.figure(2)
for i in range(0,n_runs):
    pyplot.plot(output_data_store[:,3*i+2], output_data_store[:,3*i+1])
#pyplot.plot(breeding_position[1],breeding_position[0],'kx',markersize=12)
pyplot.plot(initialbird_position[1],initialbird_position[0],'kx',markersize=12)
pyplot.axis([-arrayshape[1]/2,arrayshape[1]/2,arrayshape[0],0])
pyplot.xlabel('Longitude')
pyplot.ylabel('Latitude')
pyplot.title('Simulated paths of birds')
pyplot.savefig(os.path.join(run_folder,'path.pdf'))

#Plot centre of mass path (centre of mass of positions of all runs at each timepoint)
pyplot.figure(3)
pyplot.plot(COM_array[:,1], COM_array[:,0])
#pyplot.plot(breeding_position[1],breeding_position[0],'kx',markersize=12)
pyplot.plot(initialbird_position[1],initialbird_position[0],'kx',markersize=12)
pyplot.axis([-arrayshape[1]/2,arrayshape[1]/2,arrayshape[0],0])
pyplot.xlabel('Longitude')
pyplot.ylabel('Latitude')
pyplot.title('Centre of mass at each timepoint of several simulation runs')
pyplot.savefig(os.path.join(run_folder,'COM_path.pdf'))
