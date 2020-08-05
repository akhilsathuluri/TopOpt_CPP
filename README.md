# TopOpt_CPP
This repository is a single core implementation of the 88 line topology optimisation code[1] in C++ using the Eigen library. 

## Usage
You can set up the problem variables like discretisation, constraints, forces etc. in the main.cpp file and run the same to start the optimisation problem. 
The resulting structure is displayed at the end of each iteration.

The solution of the MBB problem with 60x30 discretisation is shown in Fig.1.

![Fig.1: Solution for MBB problem 60x30](sol_60_30.gif)

## References:
[1] Andreassen, Erik, et al. "Efficient topology optimization in MATLAB using 88 lines of code." Structural and Multidisciplinary Optimization 43.1 (2011): 1-16.
