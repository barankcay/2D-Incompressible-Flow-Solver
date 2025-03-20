# 2D-Incompressible-Flow-Solver
Flow solver based on Mark Owkes paper, with C++. </br>


References</br>
-Mark Owkes, "A guide to writing your first CFD solver", April 2024. </br>
-Benjamin Seibold, "A compact and fast Matlab code solving the incompressible
Navier-Stokes equations on rectangular domains", March 2008.



# 2D-Incompressible-Flow-Solver

This project implements a 2D incompressible flow solver based on computational fluid dynamics (CFD) principles. The solver is inspired by Mark Owkes' paper and includes multiple implementations for solving channel flow and lid-driven cavity (LDC) problems using C++.

## Features

- **Channel Flow Solver**: Simulates fluid flow in a channel.
- **Lid-Driven Cavity Solver**: Simulates fluid flow in a square cavity driven by the motion of one of its walls.
- **Stable and Parallel Implementations**: Includes stable and parallelized versions of the solvers for improved performance.
- **Visualization Support**: Outputs results in `.vtk` format for visualization in tools like ParaView.

## Code Overview

### [01_solvers/v1_ChannelFlow_stable.cpp](01_solvers/v1_ChannelFlow_stable.cpp)
- Implements a stable solver for channel flow problems.
- Outputs results in `.vtk` format for visualization.

### [01_solvers/v1_LDC_stable.cpp](01_solvers/v1_LDC_stable.cpp)
- Implements a stable solver for the lid-driven cavity problem.
- Focuses on numerical stability and accuracy.

### [01_solvers/v1_LDC_stable_paraview.cpp](01_solvers/v1_LDC_stable_paraview.cpp)
- Extends the stable LDC solver with enhanced support for ParaView visualization.
- Outputs detailed `.vtk` files for post-processing.

### [01_solvers/v1_LDC.cpp](01_solvers/v1_LDC.cpp)
- A basic implementation of the lid-driven cavity solver.
- Serves as a starting point for understanding the problem setup and solution.

## How to Run

1. **Compile the Code**: Use a C++ compiler (e.g., `g++`) to compile the desired solver file. For example:
   ```sh
   g++ -o ldc_solver 01_solvers/v1_LDC_stable.cpp