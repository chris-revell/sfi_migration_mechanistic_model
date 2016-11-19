
This is a Python3 program to run a mechanistic simulation of bird migration, devised by Marius Somveille and Christopher Revell at the Santa Fe Institute in June 2016.

The program imports lattices of chlorophyll density data and wind data. The chlorophyll data is filtered such that low density regions can be ignored. 

A bird is placed within the lattice and allowed to move to one of the eight surrounding lattice points. Each of these lattice points has a potential associated with it. This potential is defined by a combination of the chlorophyll and wind data. Each non-zero point in the cholorophyll concentration data produces a gravitational 1/r potential and all of these potentials contribute to the total value at each possible lattice point. The potential in each lattice is also weighted by the wind at the current location. The adjustment to the potential of each surrounding lattice point is given by |v|(u.v) where v is the wind vector at the current location and u is the unit vector from the current location towards the new lattice point in question. The relative contributions of the wind and chlorophyll components to the potential are set by parameter a. Once potentials, E, have been calculated for each of the 8 possible states, the probability that the bird will move into each state is defined by the Boltzmann factor of the state, e^(-E/kT) where kT is another parameter of the model, representing a form of energy or restlessness in the birds. Lattice points on solid ground are excluded from all calculations. 

This program is designed to be run in python3. Running in python2 will cause bugs. Matplotlib, mpl_toolkits.basemap and numpy must be installed for the program to run. 

Run the program by typing the following at the command line, where path is the location of the input data:

python3 seabird_sim.py <path>

Link to Overleaf LaTeX document for this project:
https://www.overleaf.com/5621453ctqmfc#/18190037/
