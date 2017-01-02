
This is a Python3 program to run a mechanistic simulation of bird migration, devised by Marius Somveille and Christopher Revell at the Santa Fe Institute in June 2016.

The program imports lattices of chlorophyll density data and wind data. The chlorophyll data is filtered such that low density regions can be ignored.
Every 720 hours (30 x 24 hours) new datafiles are imported so that the data changes in accordance with passing seasons.

A bird is placed within the lattice and allowed to move to one of the eight surrounding lattice points. Each of these lattice points has a potential associated with it. This potential is defined by a combination of the chlorophyll and wind data. Each non-zero point in the cholorophyll concentration data produces a gravitational 1/r potential and all of these potentials contribute to the total value at each possible lattice point. The potential in each lattice is also weighted by the wind at the current location. The adjustment to the potential of each surrounding lattice point is given by |v|(u.v) where v is the wind vector at the current location and u is the unit vector from the current location towards the new lattice point in question. The relative contributions of the wind and chlorophyll components to the potential are set by parameter a. Once potentials, E, have been calculated for each of the 8 possible states, the probability that the bird will move into each state is defined by the Boltzmann factor of the state, e^(-E/kT) where kT is another parameter of the model, representing a form of energy or restlessness in the birds. Lattice points on solid ground are excluded from all calculations.

This program is designed to be run in python3. Running in python2 will cause bugs. Matplotlib, mpl_toolkits.basemap and numpy must be installed for the program to run.

Run the program by typing the following at the command line:

python3 seabird_sim.py (arguments)
or,
./seabird_sim.py (arguments)

Where the arguments are a list in the following order:
initial_lat - Latitude of breeding location
initial_lon - Longitude of breeding location
a           - Relative contribution of wind potential to total potential where contribution of chlorophyll is 1
b           - Relative contribution of breeding ground attraction potential to total potential where contribution of chlorophyll is 1
c           - Exponential factor for variation of breeding ground attraction with time
kT          - Effective temperature of bird - high value increases randomness in path, lower value collapses to lowest energy state.
start_month - Start month (range 1-12) defines the end of the breeding season when birds leave the breeding colony, whose position is defined by initial_lat and initial_lon
end_month   - (Optional) specifies when to stop simulation if a full year is not required

Link to Overleaf LaTeX document for this project:
https://www.overleaf.com/5621453ctqmfc#/18190037/

Things left to do:
  - Implement continuous time and bird speed rather than discrete time. *DONE*
  - Implement updating of environment data arrays with passage of time, eg. July data replaced with August data after ~720 hours. *DONE*
  - Restrict distance over which chlorophyll can exert attraction?
  - Exclude very high chlorophyll values? (often due to lakes etc?) *DONE but need to justify thresholding bounds*
  - Introduce penalty for change of direction? *Probably no need*
  - Identify how to collect data.
    > Perform many runs and create heat maps of bird distribution at different time points?
    > Marius to identify relevant measurements from system that would be useful to ornithologists?
  - Change bird speed in accordance with wind speed? *DONE*
  - Improve run time with fortran subroutines? *DONE*

Compile fortran subroutines with:
f2py -c -m fortran_subroutines fortran_subroutines.f95

Overleaf document for this project:
https://www.overleaf.com/7609498qptrwdvytdkn
