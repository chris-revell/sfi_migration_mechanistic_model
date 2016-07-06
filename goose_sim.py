
from sys import argv
import numpy as np
import random
from math import exp
import matplotlib.pyplot as plt

infilename = argv[1]

#File path for datafile given at command line
#Use numpy.genfromtxt to import matrix in .txt file into array, skipping first row and column
#Data is imported in string format and contains "NA" values for sea areas.
#Also in data file, greener points have higher NDVI values.
#We want the potential to be lower at greener areas, so we multiply every value by -1 when using as a potential
infile = open(infilename,'r')

#Number of columns and rows in data file
nrows=0
ncols=0
for line in infile:
    nrows = nrows+1
    if nrows == 2:
        ncols = len(line.split(' '))
infile.close()

#Import data matrix
NDVI_import = np.genfromtxt(infilename, dtype=str, skip_header=1, usecols=range(1,ncols), delimiter=' ')

#Convert to number of columns and rows in array.
nrows=nrows-1
ncols=ncols-1

plt.axis([0, ncols, nrows, 0])
plt.ion()

#Create array to store the distance from the breeding ground at each lattice point.
r_i_array = np.zeros(NDVI_import.shape)

#Define position of breeding ground and initial position of goose
breeding_position = (279,1147) #(0,ncols-1)

#x = "NA"
#while x == "NA":
#    coordinates = (int(random.random()*699),int(random.random()*1000))
#    x = NDVI_import[coordinates]
#    goose_position = coordinates
goose_position = (495,560)

#Also define array to hold Boltzmann factors with same dimensions as data array
boltzmann_factors = np.zeros(NDVI_import.shape)
#Calculate values for all elements of boltzmann factor array.
#If the corresponding element in the data array is "NA", the Boltzmann factor is set to 0.
for x in range(0,nrows):
    for y in range(0,ncols):
        if NDVI_import[x,y] == "NA":
            boltzmann_factors[x][y] = 0
        else:
            #Define distance of (x,y) from breeding position - put this into a function?
            r_i_vector = [(breeding_position[0]-x),(breeding_position[1]-y)]
            r_i        = (np.vdot(r_i_vector,r_i_vector))**0.5+0.1
            r_i_array[x,y] = r_i
            #Define potential at (x,y) from NDVI data and "breeding location gravity"
            #Note that factor of -1 in exponential is included in "potential" value
            potential = (float(NDVI_import[x,y])/100 + 10/r_i)/1000
            boltzmann_factors[x,y] = exp(potential)

#We have now defined the Boltzmann factor array that will be used to determine probabilities as the goose moves through the lattice.

#Open file into which results are printed
outfile = open(infilename[0:-4]+'_goose_positions.txt','w')

#Loop over timesteps
for t in range (0,1000000):
    #Define possible_states array to hold the neighbouring lattice points that a goose can move into. Generally 9, fewer at system edges.
    possible_states = []
    #Identify current neighbouring elements
    #Do not include values outside the bounds of the array
    for x in range(-1,2):
        if 0 <= (goose_position[0]+x) < nrows:
            for y in range(-1,2):
                if 0 <= (goose_position[1]+y) < ncols :
                    possible_states.append((goose_position[0]+x,goose_position[1]+y))
                else:
                    pass
        else:
            pass

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

    #Now that we have updated the goose position, write it to output file
    output = str(t)+'  '+str(goose_position[0])+'  '+str(goose_position[1])+'  '+str(r_i_array[goose_position[0],goose_position[1]])+'\n'
    outfile.write(output)
    if t%1000 == 0:
        plt.scatter(goose_position[0], goose_position[1])
        plt.pause(0.)

#while True:
#    plt.pause(0.000005)

#print (boltzmann_factors[goose_position[0]-1,goose_position[1]-1],boltzmann_factors[goose_position[0]-1,goose_position[1]],boltzmann_factors[goose_position[0]-1,goose_position[1]+1])
#print (boltzmann_factors[goose_position[0],goose_position[1]-1],boltzmann_factors[goose_position[0],goose_position[1]],boltzmann_factors[goose_position[0],goose_position[1]+1])
#print (boltzmann_factors[goose_position[0]+1,goose_position[1]-1],boltzmann_factors[goose_position[0]+1,goose_position[1]],boltzmann_factors[goose_position[0]+1,goose_position[1]+1])
#print (NDVI_import[goose_position[0]-1,goose_position[1]-1],NDVI_import[goose_position[0]-1,goose_position[1]],NDVI_import[goose_position[0]-1,goose_position[1]+1])
#print (NDVI_import[goose_position[0],goose_position[1]-1],NDVI_import[goose_position[0],goose_position[1]],NDVI_import[goose_position[0],goose_position[1]+1])
#print (NDVI_import[goose_position[0]+1,goose_position[1]-1],NDVI_import[goose_position[0]+1,goose_position[1]],NDVI_import[goose_position[0]+1,goose_position[1]+1])
#final_outfile = open('finalgooseposition.txt','w')
#final_outfile.write(str(goose_position[0])+'  '+str(goose_position[1]))

#outfile2 = open('NDVI_import.txt','w+b')
#outfile3 = open('boltzmann_factors.txt','w+b')
#np.savetxt(outfile2,NDVI_import)
#np.savetxt(outfile3,boltzmann_factors)
#outfile.close()
#print(ncols,nrows)
