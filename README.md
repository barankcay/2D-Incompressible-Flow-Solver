# 2D-Incompressible-Flow-Solver

This project implements a 2D incompressible flow solver based on computational fluid dynamics (CFD) principles. The solver is inspired by Mark Owkes' paper and includes multiple implementations for solving channel flow and lid-driven cavity (LDC) problems using C++.

## Features

- **Channel Flow Solver**: Simulates fluid flow in a channel.
- **Lid-Driven Cavity Solver**: Simulates fluid flow in a square cavity driven by the motion of one of its walls.
- **Stable Implementations**: Includes different stability controls such as Peclet & Courant number and Courant number alone.
- **Visualization Support**: Outputs results in `.vtk` format for visualization in tools like ParaView.

## Code Overview

The solver is implemented in C++ and consists of the following key components:

- **Boundary Conditions**: Functions to apply boundary conditions for velocity and pressure fields.
- **Fractional Step Method**: The solver uses a fractional step method to decouple the velocity and pressure calculations. This involves:
  1. **Predictor Step**: Computes intermediate velocity fields (`u*`, `v*`) without considering pressure.
  2. **Pressure Correction**: Solves the Poisson equation to compute the pressure correction.
  3. **Corrector Step**: Updates the velocity fields using the corrected pressure to ensure incompressibility.
- **Dynamic Time Step Adjustment**: Ensures stability using Courant and Peclet numbers.
- **Convergence Check**: Monitors residuals to determine when the solution has converged.
- **VTK Output**: Exports results in `.vtk` format for visualization in ParaView.

The main program initializes parameters, runs the simulation loop, and outputs results for post-processing.

## How to Compile

1. **Install a C++ Compiler**:
   - Ensure you have a C++ compiler installed, such as `g++` (part of GCC).

2. **Navigate to the Solver Directory**:
   - Open a terminal or command prompt and navigate to the directory containing `v1_Final.cpp`.

3. **Compile the Code**:
   - Use the following command to compile the code with optimizations:
     ```bash
     g++ -std=c++17 .\v1_Final.cpp -o son.exe -O3 -ffast-math
     ```
   - The `-O3` flag enables high-level optimizations, and `-ffast-math` allows faster floating-point calculations.

4. **Run the Executable**:
   - Execute the compiled program:
     ```bash
     .\son.exe
     ```

5. **Answer the Prompts**:
   - The program will ask for input parameters. Answer the questions based on the options provided:
     - **Reynolds Number**: Enter a numerical value (e.g., `100`).
     - **Time Step Calculation**: Choose `(Y/N)` to calculate the time step size based on Courant and Peclet numbers. If `Y`, both Cr and Pe is included. If `N`, then only Cr is included.
     - **Write VTK File**: Choose `(Y/N)` to enable or disable `.vtk` file output.
     - **Lid-Driven Cavity Case**: Choose `(Y/N)` to specify the type of simulation.
       - If `N`, the case becomes a channel flow case and ensure to provide additional inputs for outlet pressure and inlet velocities.

6. **Output**:
   - The program will run the simulation and generate results, including a `.vtk` file for visualization.

7. **Visualize Results**:
   - Open the `.vtk` file in ParaView or another visualization tool to analyze the simulation output.

## References
- Mark Owkes, "A guide to writing your first CFD solver", April 2024.
- U. Ghia, K.N. Ghia, and C.T. Shin, "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method," Journal of Computational Physics, vol. 48, no. 3, pp. 387–411, 1982, doi: 10.1016/0021-9991(82)90058-4
