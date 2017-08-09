#######################################################################
#
#   Ethasham Ahmed, December 2016
#
#   This is a simple molecular dynamics program that 
#   simulates the collison of two clusters. 
#
#   Argon has been chosen arbitrarily as the particle
#   within the clusters. 
#
#   The program uses Lennard Jones potential. 
#
#   Velocity verlet algorithm is used to solve newtons 
#   equation of motion  
#
#   Equipartition theorem is used to calculate temperature
#
#   Spatial dimentions (DIM) is set to 3 by default
#
#   Time step (dt) is set to 1.7e-4 by default
#
#   Within the code, calculations are done in reduced units
#   which means mass, sigma and epsiolon are 1
#
#   However, by the end of the code, all values are 
#   converted back to SI units. 
#
#   Initial positions of the particles are determined 
#   from two input files
#
#   The program writes to a .xyz file every 10 steps by defult.
#######################################################################
