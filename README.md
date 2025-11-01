# Flocking-with-Fuel
This is the C++ code which generates the flocking simulations with fuel consumption. Starting with random positions and orientations, particles aggregate and flock when consuming the fuel. After running out fuel the system becomes passive and particles disperse.

## Library
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) is required to run the code. Details can be found in CMakeLists.txt.

## Coding Platform
We suggest using [CLion](https://www.jetbrains.com/clion/) as a simple way to test the code. Each single simulation could take around 3 hours, depending on the machine & chips. Running on supercomputing cloud may be faster. This time estimate is based on simulations with $N=1000$ particles.

## Visualisation
The simulation output is of the form '.dump'. Users can use the open visualization tool [OVITO](https://www.ovito.org) to open and visualise the simulation results. 

## Parameters
The parameters set in the main.cpp will produce a simulation showing the transition: disorder (with fuel) - order (less fuel) - disorder (fuel runs out). Changing $\zeta_\theta = \frac{5}{64}, J_0=\frac{7}{2048}$ and varying the noise amplitude $\Theta$, defined by $\Theta \equiv q_r = q_t = q_n$ will generate a set of simulations, where the mean order parameter $v_a$ decays as the noise amplitude $\Theta$ increases.
