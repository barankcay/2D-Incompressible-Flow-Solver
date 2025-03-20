# 2D-Incompressible-Flow-Solver

This project implements a 2D incompressible flow solver based on computational fluid dynamics (CFD) principles. The solver is inspired by Mark Owkes' paper and includes multiple implementations for solving channel flow and lid-driven cavity (LDC) problems using C++.

## Features

- **Channel Flow Solver**: Simulates fluid flow in a channel.
- **Lid-Driven Cavity Solver**: Simulates fluid flow in a square cavity driven by the motion of one of its walls.
- **Stable Implementations**: Includes different stability controls such as Peclet & Courant number and Courant number alone.
- **Visualization Support**: Outputs results in `.vtk` format for visualization in tools like ParaView.

## Code Overview

### [01_solvers/v1_ChannelFlow_CrPe.cpp](01_solvers/v1_ChannelFlow_CrPe.cpp)
- Implements a stable solver for channel flow problems.
- Includes Courant and Peclet number control.
- Outputs results in `.vtk` format for visualization.

### [01_solvers/v1_LDC_CrPe.cpp](01_solvers/v1_LDC_CrPe.cpp)
- Implements a stable solver for the lid-driven cavity problem.
- Includes Courant and Peclet number control.
- Focuses on numerical stability and accuracy.

### [01_solvers/v1_LDC_CrPe_paraview.cpp](01_solvers/v1_LDC_CrPe_paraview.cpp)
- Extends the stable LDC solver with enhanced support for ParaView 
visualization.
- Includes Courant and Peclet number control.
- Outputs detailed `.vtk` files for post-processing.

### [01_solvers/v1_LDC_Cr.cpp](01_solvers/v1_LDC_Cr.cpp)
- A basic implementation of the lid-driven cavity solver.
- Includes Courant and Peclet number control only.
- Serves as a starting point for understanding the problem setup and solution.

## How to Run

1. **Compile the Code**: Use a C++ compiler (e.g., `g++`) to compile the desired solver file. For example:
   ```sh
   g++  01_solvers/v1_LDC_CrPe_paraview.cpp -o IncompSolver.exe -O3 -ffast-math

2. **Run the Executable**: Execute the compiled binary to run the solver.
    ```sh
    .\IncompSolver.exe

3. **Visualize**: Open the generate `.vtk` file in ParaView for visualization

## References
-Mark Owkes, "A guide to writing your first CFD solver", April 2024.
-U. Ghia, K.N. Ghia, and C.T. Shin, "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method," Journal of Computational Physics, vol. 48, no. 3, pp. 387–411, 1982, doi: 10.1016/0021-9991(82)90058-4
