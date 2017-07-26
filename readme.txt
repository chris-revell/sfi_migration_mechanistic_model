
This is a Python3 program to run a mechanistic simulation of bird migration, devised by Christopher Revell and Marius Somveille at the Santa Fe Institute in June 2016.

This work is published in Revell, C. and Somveille, M. (2017). A Physics-Inspired Mechanistic Model of Migratory Movement Patterns in Birds. Scientific Reports.

The program imports lattices of chlorophyll density data and wind data from folders data_chloro and data_wind. Every 720 hours (30 x 24 hours) new data files are imported so that the data changes in accordance with passing seasons. A lattice of ocean map is imported from earth.txt to mask chlorophyll data from inland lakes.

A bird is placed within the lattice and allowed to move to one of the eight surrounding lattice points. Each of these lattice points has a potential associated with it. This potential is defined by a combination of the chlorophyll and wind data. Once potentials, E, have been calculated for each of the 8 possible states, the probability that the bird will move into each state is defined by the Boltzmann factor of the state, e^(-E/kT) where kT is another parameter of the model. Lattice points on solid ground are excluded from all calculations.

This program is designed to be run in python3. Running in python2 will cause bugs. NumPy must be installed. In order to run the accompanying data processing scripts, Matplotlib and mpl_toolkits.basemap must also be installed.

Before running the program for the first time, it is necessary to compile the Fortran subroutine. This requires f2py and is achieved by entering the following at the command line:
f2py -c -m seabird_subroutines seabird_subroutines.f95

Run the program by typing the following at the command line:

python3 seabird_sim.py (arguments)
or,
./seabird_sim.py (arguments)

Where the arguments are a list in the following order:
initial_lat - Latitude of breeding location
initial_lon - Longitude of breeding location
a           - Relative contribution of wind potential to total potential where contribution of chlorophyll is 1
kT          - Effective temperature of bird - high value increases randomness in path, lower value collapses to lowest energy state.
start_month - Start month (range 1-12) defines the end of the breeding season when birds leave the breeding colony, whose position is defined by initial_lat and initial_lon
end_month   - (Optional) specifies when to stop simulation if a full year is not required
run         - Run label (integer) for when multiple runs with the same set of parameters are performed

The program will store data in a new folder labelled with the run parameters. If it does not already exist, the program will also create folder output_data in which to store each individual run data folder.
