
import numpy as np
import random
from math import exp

#NDVI_import = np.array([])

#NDVI_file = open('filename','r')

#for line in NDVI_file:
#    line_stripped = line.rstrip()
#    line_array    = line_stripped.split(" ")
    #Add line to data array
#    NDVI_import = np.vstack([NDVI_import,line_array])

#Should now have NDVI values imported into array NDVI_import


outfile = open('test_data.txt','w')
outfile2 = open('potentials.txt','w+b')

NDVI_import = np.zeros((10,10))
boltzmann_factors = np.zeros((10,10))

breeding_position = [9,9]
goose_position = [0,0]

for x in range (0,10):
    for y in range (0,10):
        r_i_vector = [(breeding_position[0]-x),(breeding_position[1]-y)]
        r_i        = (np.vdot(r_i_vector,r_i_vector))**0.5+0.1
        NDVI_import[x,y] = random.random() - 0.1/r_i
        boltzmann_factors[x,y] = exp(-1.0*NDVI_import[x,y])

#Now have an array of random potential values tilted by attractor at (9,9)

for t in range (0,200):
    possible_states = []
    for x in range(-1,2):
        if 0 <= (goose_position[0]+x) <= 9:
            #Do not include neighbours outside bounds or array
            for y in range(-1,2):
                if 0 <= (goose_position[1]+y) <= 9:
                    possible_states.append((goose_position[0]+x,goose_position[1]+y))
                else:
                    pass
        else:
            pass

    boltzmann_sum = 0
    for z in possible_states:
        boltzmann_sum = boltzmann_sum+boltzmann_factors[z]

    probability_sum = 0
    random_number = boltzmann_sum*random.random()
    for p in possible_states:
        probability_sum = probability_sum + boltzmann_factors[p[0],p[1]]
        if random_number < probability_sum:
            goose_position = [p[0],p[1]]
            break
        else:
            pass

    output = str(t)+'  '+str(goose_position[0])+'  '+str(goose_position[1])+'\n'
#    outfile.write(output)
#    print(goose_position)

final_outfile = open('finalgooseposition.txt','w')
final_outfile.write(str(goose_position[0])+'  '+str(goose_position[1]))

np.savetxt(outfile2,NDVI_import)
np.savetxt(outfile2,boltzmann_factors)
