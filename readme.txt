
This is a Python3 program to run a mechanistic simulation of bird migration, devised by Marius Somveille and Christopher Revell at the Santa Fe Institute in June 2016.

The program takes as input a set of NDVI data files. The folder containing these files should be passed to the program at the command line. The files should be .txt format, have alphabetical order identical to temporal order, and should be the only .txt files in the directory.

The model has 3 parameters: the number of time steps and factors A and kT. These 3 parameters are defined within the first few lines of the code.  Value of t_max currently set to 300000. Values of A and kT given at command line. Appropriate values of A are of the order 10^7; appropriate values of kT are of the order 10^2-10^3.

A starting position for the bird, and the location of the breeding site should also be defined within the code. Currently these are set as observed from tracking data of White-fronted Geese.

The simulation outputs data in a directory labelled by date and time of the run. Within this directory can be found a file containing the simulation parameters, a file containing the starting position and breeding location of the geese, and a file giving a time series of 2D goose positions and distance from breeding ground. The program also produces a set of gnuplot commands that will draw the data.

This program is designed to be run in python3. In python2, a bug arises that causes .pop(0) to be applied to an empty list. This is because python3 uses float division by default, whereas python2 uses integer division.

Run the program by using the following command line command:

python3 goose_sim.py <data/folder/path> <value_of_A> <value_of_kT>

Link to Overleaf LaTeX document for this project:
https://www.overleaf.com/5621453ctqmfc#/18190037/

Bash command for multiple runs:
for i in {1..10}; do for j in {1..10}; do python3 goose_sim.py NDVI_data/spring $(($i*100000)) $(($j*250)); done; done
