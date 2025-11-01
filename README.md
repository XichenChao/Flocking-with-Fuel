# Flocking-with-Fuel
This is the C++ code which generates the simulations shown in the paper https://arxiv.org/abs/2507.14294. Using this code is welcome, and users should cite the paper https://arxiv.org/abs/2507.14294 (or its published version if applicable).

## Library
Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page) is required to run the code. Details can be found in CMakeLists.txt.

## Coding Platform
We suggest using Clion (https://www.jetbrains.com/clion/) as a simple way to test the code. Each single simulation could take around 3 hours, depending on the machine & chips. Running on supercomputing cloud may be faster. This time estimate is based on simulations with $N=1000$ particles.

## Visualisation
The simulation output is of the form '.dump'. Users can use the open visualization tool OVITO (https://www.ovito.org) to open and visualise the simulation results. 

## Parameters
