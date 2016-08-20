
from sys import argv
import numpy as np
import random
from math import exp
import os
import time
import matplotlib.pyplot as pyplot

#Folder path for datafile given at command line
importfolderpath = argv[1]

#Set system parameters
A     = int(argv[2]) # Prefactor for breeding site gravitational attraction. Given at command line. ~10^5
kT    = int(argv[3]) # Measure of goose temperature or "restlessness". Given at command line. ~10^3
n_runs = 5 # Number of runs with this set of parameters. Program will produce an average and standard deviation over all runs.
#Define position of breeding ground and initial position of goose
breeding_position = (279,1147)
goose_position    = (495,560)

#Create array to store data from all runs
output_data_store = np.empty((t_max,(n_runs*3+1)))
for i in range(0,t_max):
    output_data_store[i,0] = i

#From folder path provided at command line, find list of files to import NDVI data from.
#Each file corresponds to half a month. "isfile" checks that we find only files, not directories.
#f[-1]=='t' excludes anything but .txt files (specifically to exclude .ds_store file on Mac)
#Note that list will be ordered alphabetically, so alphabetical order of filenames must match temporal order of months.
datafiles = [f for f in os.listdir(importfolderpath) if os.path.isfile(os.path.join(importfolderpath, f)) and f[-1]=='t']
datafiles.sort()
print(datafiles)

#If each datafile is half a month, this corresponds to roughly 15 days per file, which is 360 hours.
#Thus set the total run time in hours from the total number of input files minus 1, multiplied by 360.
#Minus 1 required becauase the system runs on the interpolated spaces between files.
t_max = 360*(len(datafiles)-1)  # Total time in hours for simulation
#Set output interval
output_interval = t_max/1000 #output_interval ~ 1 hour

#Create date and time labelled folder to store the data from this run
#First ensure the output_data folder exists
if os.path.exists("../output_data"):
    pass
else:
    os.mkdir("../output_data")
run_folder = '../output_data/A'+argv[2]+'kT'+argv[3]+'_'+time.strftime("%y%m%d%H%M")
os.mkdir(run_folder)

#Write initial conditions to file
#Wintering position and breeding position
winterbreedingpositionfile = open(run_folder+'/winterbreedingposition.txt','w')
winterbreedingpositionfile.write(str(goose_position[0])+' '+str(goose_position[1])+'\n'+str(breeding_position[0])+' '+str(breeding_position[1]))
winterbreedingpositionfile.close()
#Write simulation parameters to data file
parameterfile = open(run_folder+'/parameters.txt','w')
parameterfile.write("A  "+str(A)+"\nkT  "+str(kT)+"\nt_max  "+str(t_max)+'\nDatafiles:')
parameterfile.write("  ".join(datafiles))

#Set interval for importing new datafiles. Note -1 because the system ends once interpolation reaches the state of the final file.
update_interval = t_max/(len(datafiles)-1)

#Use numpy.genfromtxt to import matrix in .txt file into array, skipping first row and column
#Data is imported in string format and contains "NA" values for sea areas.
#Also in data file, greener points have higher NDVI values.
#We want the potential to be lower at greener areas, so we multiply every value by -1 when using as a potential
infile = open(importfolderpath+'/'+datafiles[0],'r')
#Number of columns and rows in data file
nrows=0
ncols=0
for line in infile:
    nrows = nrows+1
    if nrows == 2:
        ncols = len(line.split(' '))
infile.close()

#Use number of columns and rows to define the shape arrays needed later, which must have the same shape as the NDVI data array.
#These are the boltzmann_factors array, array of distances from breeding ground, NDVI_interpolated array and NDVI_gradient array for interpolating.
boltzmann_factors = np.zeros((nrows-1,ncols-1))
r_i_array         = np.zeros((nrows-1,ncols-1))
NDVI_gradient     = np.empty((nrows-1,ncols-1),dtype=float)
NDVI_interpolated = np.empty((nrows-1,ncols-1),dtype=float)
NDVI_next         = np.empty((nrows-1,ncols-1),dtype=int)
time_updated      = np.zeros((nrows-1,ncols-1),dtype=int)

#Fill r_i_array with distances from breeding site.
for x in range(0,nrows-1):
    for y in range(0,ncols-1):
        #Define distance of (x,y) from breeding position
        r_i_vector = [(breeding_position[0]-x),(breeding_position[1]-y)]
        r_i        = (np.vdot(r_i_vector,r_i_vector))**0.5+0.1
        r_i_array[x,y] = r_i

#Import initial NDVI file
initialfilename = importfolderpath+'/'+datafiles.pop(0)
NDVI_import = np.genfromtxt(initialfilename, dtype=str, skip_header=1, usecols=range(1,ncols), delimiter=' ')
#Data imported as strings. Convert to interger format, setting "NA" values to unique intger label 12345
for x in range(0,nrows-1):
    for y in range(0,ncols-1):
        if NDVI_import[x,y] == "NA":
            NDVI_next[x,y] = 12345
        else:
            NDVI_next[x,y] = int(NDVI_import[x,y])

#Define array of possible lattice points that the bird can be in at the next timestep.
#Generally has 9 components, but can have fewer at the edges of the system.
possible_states = []

#Define subroutine to find the set of possible lattice points that the bird can be in at the next timestep.
def find_possible_states():
    global possible_states

    #Start by removing all elements in previous possible states list in case the number of elements in the new list is smaller
    del possible_states[:]

    #Identify current neighbouring elements
    #Do not include values outside the bounds of the array
    for x in range(-1,2):
        if 0 <= (goose_position[0]+x) < nrows-1:
            for y in range(-1,2):
                if 0 <= (goose_position[1]+y) < ncols-1 :
                    possible_states.append((goose_position[0]+x,goose_position[1]+y))
                else:
                    pass
        else:
            pass

#Define functional form of breeding site gravity potential
def breeding_gravity(radius):
    return (A/radius)

#Subroutine to update the boltzmann factors for the possible states in the next timestep after the NDVI values of these states have been interpolated.
def boltzmann_update(possible_states):
    global boltzmann_factors
    global NDVI_interpolated
    for element in possible_states:
        if NDVI_interpolated[element] == 12345:
            #12345 corresponds to "NA" value in original NDVI data. These points correspond to sea. We exclude birds from these points by setting the boltzmann factor to 0.
            boltzmann_factors[x][y] = 0
        else:
            #Define potential at (x,y) from NDVI data and "breeding location gravity"
            #Note that factor of -1 in exponential is included in "potential" value
            potential = (float(NDVI_interpolated[element]) + breeding_gravity(r_i_array[element]))/kT
            boltzmann_factors[element] = exp(potential)

#Define subroutine to update the position of the goose
def system_update(t,n):
    global boltzmann_factors
    global goose_position
    global possible_states
    global output_data_store

    #Sum Boltzmann factors for possible states
    boltzmann_sum = 0
    for z in possible_states:
        boltzmann_sum = boltzmann_sum+boltzmann_factors[z]

    #Use a random number generator and probabilities defined by Boltzmann factors
    #to decide which lattice point the goose moves to next
    probability_sum = 0
    random_number = boltzmann_sum*random.random()
    for p in possible_states:
        probability_sum = probability_sum + boltzmann_factors[p[0],p[1]]
        if random_number < probability_sum:
            goose_position = [p[0],p[1]]
            break
        else:
            pass

#Subroutine to import a new NDVI file and redefine the previous "next" file as the new "current" file.
#Calculates the gradient array between the two files.
def importnext(t):
    global NDVI_interpolated
    global NDVI_gradient
    global NDVI_next
    global time_updated
    #Set the current NDVI dataset to the previous "next" dataset
    NDVI_interpolated = NDVI_next
    #Use datafiles list to find the filename of the next file to import NDVI data from
    filenameatinterval = importfolderpath+'/'+datafiles.pop(0)
    #Use this filename to import data from file into array
    NDVI_import = np.genfromtxt(filenameatinterval, dtype=str, skip_header=1, usecols=range(1,ncols), delimiter=' ')
    #NDVI data imported as strings. Need to convert to integers. "NA" values set to integer value 12345, which is beyond the bounds of NDVI data and therefore a unique integer label.
    for x in range(0,nrows-1):
        for y in range(0,ncols-1):
            if NDVI_import[x,y] == "NA":
                NDVI_next[x,y] = 12345
            else:
                NDVI_next[x,y] = int(NDVI_import[x,y])

    #Calculate gradient for interpolation by broadcasting subtraction and then division by interval over all components in data arrays.
    NDVI_gradient = (NDVI_next-NDVI_interpolated)/update_interval
    print (t, filenameatinterval)
    #Change the value of the time the element was last updated for every element of the array.
    time_updated.fill(t-1)

#Function to interpolate NDVI data between two files
def interpolate(possible_states,t):
    global NDVI_gradient
    global NDVI_interpolated
    global time_updated
    #Loop over all elements in array
    #Use difference between current time t and the last time an element was updated to determine how many multiples of the gradient should be added to interpolate that element.
    for element in possible_states:
        NDVI_interpolated[element] = NDVI_interpolated[element] + (t-time_updated[element])*NDVI_gradient[x,y] #(t-time_updated[element]) gives the time that has elapsed since that component of the array was last interpolated, and therefore the magnitude of the prefactor required when adding some multiple of the gradient
        #Element updated, so set the time_updated to the current time t.
        time_updated[element] = t

for i in range (0,n_runs):
    print('run '+str(i))
    #Reset system for each new run
    goose_position    = (495,560)
    datafiles = [f for f in os.listdir(importfolderpath) if os.path.isfile(os.path.join(importfolderpath, f)) and f[-1]=='t']
    datafiles.sort()
    NDVI_import = np.genfromtxt(initialfilename, dtype=str, skip_header=1, usecols=range(1,ncols), delimiter=' ')
    #Data imported as strings. Convert to interger format, setting "NA" values to unique intger label 12345
    for x in range(0,nrows-1):
        for y in range(0,ncols-1):
            if NDVI_import[x,y] == "NA":
                NDVI_next[x,y] = 12345
            else:
                NDVI_next[x,y] = int(NDVI_import[x,y])
    #Find possible states for first run of system


    find_possible_states()
    #Loop over timesteps
    boltzmann_update(possible_states)
    while t<t_max:

        #For every import interval, import a new file into NDVI_next and redefine NDVI_interpolated to hold the old values of NDVI_next
        if (int(t%update_interval) == 0):
            importnext(t)

        #Update system state according to current interpolated NDVI values and corresponding BOltzmann factors.
        system_update(t,i)

        #Find possible states for next run of system
        find_possible_states()

        #Update the interpolated NDVI array
        interpolate(possible_states,t)

        #Update Boltzmann factors according to the new interpolated NDVI values
        boltzmann_update(possible_states)

        t=

#Write stored data array to file
np.savetxt(run_folder+'/goose_positions.txt', output_data_store, delimiter='  ')
outfile2 = open(run_folder+'/distance_average.txt','w')

#Calculate mean and standard deviation at all timepoints.
mean_list    = np.zeros(t_max)
std_dev_list = np.zeros(t_max)
time_list    = np.zeros(t_max)
for i in range (0,t_max):
    mean    = 0
    std_dev = 0
    for j in range (0,n_runs):
        mean    = mean + output_data_store[i,3*j+3]
        std_dev = std_dev + output_data_store[i,3*j+3]**2
    mean    = mean/n_runs
    std_dev = (std_dev/n_runs - mean**2)**0.5
    mean_list[i]    = mean
    std_dev_list[i] = std_dev
    time_list[i]    = i
    #Write data to file
    outfile2.write(str(i)+'  '+str(mean)+'  '+str(std_dev)+'\n')

#Calculate centre of mass position for all runs at each time point
COM_array = np.zeros((t_max,2))
for i in range(0,n_runs):
    COM_array[:,0] = COM_array[:,0] + output_data_store[:,3*i+1]
    COM_array[:,1] = COM_array[:,1] + output_data_store[:,3*i+2]
COM_array = COM_array/n_runs

#Plot data with pyplot
#Plot distance against time
pyplot.figure(1)
pyplot.plot(mean_list)
pyplot.fill_between(time_list,(mean_list+std_dev_list),(mean_list-std_dev_list), alpha=0.5)
pyplot.axis([0,t_max,0,700])
pyplot.xlabel('Time')
pyplot.ylabel('Distance from breeding ground')
pyplot.title('Distance from breeding ground against time')
#Need to add title and monthly dates
pyplot.savefig(os.path.join(run_folder,'distance.pdf'))

#Plot paths
pyplot.figure(2)
for i in range(0,n_runs):
    pyplot.plot(output_data_store[:,3*i+2], output_data_store[:,3*i+1])
pyplot.axis([0,2000,700,0])
pyplot.xlabel('Longitude')
pyplot.ylabel('Latitude')
pyplot.title('Simulated paths of geese')
pyplot.savefig(os.path.join(run_folder,'path.pdf'))

#Plot centre of mass path (centre of mass of positions of all runs at each timepoint)
pyplot.figure(3)
pyplot.plot(COM_array[:,1], COM_array[:,0])
pyplot.axis([0,2000,700,0])
pyplot.xlabel('Longitude')
pyplot.ylabel('Latitude')
pyplot.title('Centre of mass at each timepoint of several simulation runs')
pyplot.savefig(os.path.join(run_folder,'COM_path.pdf'))
