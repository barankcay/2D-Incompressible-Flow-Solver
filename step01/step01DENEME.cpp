/* This is STEP 01 of a series of incompresible Navier-Stokes solvers.

   It solves 2D, incompressible Navier-Stokes equations using the projection (fractional step) method.
   It uses the finite difference method to simulate the lid driven cavity flow.
   The mesh is Cartesian, uniform & staggered.
   Horizontal velocity components (u) are placed to the left of the pressure nodes.
   Vertical velocity components (v) are placed to the left of the pressure nodes.
   
   Boundary conditions are implemented using ghost nodes.
   The left, right and bottom boundaries are no slip walls.
   The top boundary is a moving wall with velocity u = 1, v = 0.
   All boundaries have zero pressure gradient.
   
   Time discretization is done using the first-order, explicit Euler method.
   The spatial discretization is done using second-order central differencing.
   The Poisson equation for pressure is solved using the Gauss-Seidel method.

   Average changes of u, v and p unknowns are calculated to check for Gauss-Seidel's convergence,
   and als to check steady state convergence.
  
   To compile: g++ fileName.cpp -o fileName.exe -O3 -ffast-math
   
   To compile with real-time plotting: 
                1 - Uncomment lines 44,45 and REALTIME PLOTTING section in the code.
                2 - Install matplotlibcpp library and its dependencies.
                To compile with Python 3.12 and NumPy support:
                        g++ fileName.cpp -o fileName.exe -O3 -ffast-math 
                        -I C:\Python\Python312\include -I C:\Python\Python312\Lib\site-packages\numpy\core\include 
                        -L C:\Python\Python312\libs -lpython312

   Authors: Ibrahim Baran Kucukcay
            Dr. Cuneyt Sert
            Dept. of Mechanical Eng., Middle East Technical Uni., Ankara, Turkey
*/


#include <iostream>
#include <vector>
#include <fstream>  // Include for file handling
#include <iomanip>  // Include for fixed precision formatting
#include <cmath>
#include <thread>
#include <chrono>
using namespace std;

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

// Add these variables after your other variable declarations
vector<double> time_plot;
vector<double> velocity_plot;
void calculateGhostNodeValues(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall,
                              double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny);

int main()
{
    auto start = std::chrono::steady_clock::now();

    //////////////////////////////////////////////////
    // CREATION OF FIELD PARAMETERS
    //////////////////////////////////////////////////
    vector<vector<double>> u;      // u velocity field
    vector<vector<double>> uStar;  // intermediate u velocity field
    vector<vector<double>> v;      // v velocity field
    vector<vector<double>> vStar;  // intermediate v velocity field
    vector<vector<double>> p;      // pressure field
    vector<vector<double>> pPrev;  // previous pressure field for convergence check
    vector<vector<double>> pOld;   // old pressure field for Gauss-Seidel method

    vector<vector<double>> uPrev; // previous u velocity field for convergence check
    vector<vector<double>> vPrev; // previous v velocity field for convergence check

    double uAdvection; // Total advection term for u. v * du/dy + u * du/dx
    double vAdvection; // Total advection term for v. u * dv/dx + v * dv/dy

    double uDiffuse;   // d²u/dx² + d²u/dy² (diffusion in X and Y for u)
    double vDiffuse;   // d²v/dx² + d²v/dy² (diffusion in X and Y for v)

    double velocityStarGrad; // Gradient of (uStar, vStar) vector used in the pressure Poisson equation. du*/dx + dv*/dy
    

    
    //////////////////////////////////////////////////
    // PROBLEM INPUTS
    //////////////////////////////////////////////////
    double Re = 3200;
    double density = 1.0;
    double kinematicViscosity = 1.0 / Re;
    double dynamicViscosity = kinematicViscosity * density;
    

    
    //////////////////////////////////////////////////
    // BOUNDARY CONDITIONS
    //////////////////////////////////////////////////
    double uTopWall    = 1.0;   // x direction velocity at the top wall
    double uBottomWall = 0.0;   // x direction velocity at the bottom wall
    double uLeftWall   = 0.0;   // etc.
    double uRightWall  = 0.0; 
    double vTopWall    = 0.0; 
    double vLeftWall   = 0.0;
    double vRightWall  = 0.0;
    double vBottomWall = 0.0;   
    

    
    //////////////////////////////////////////////////
    // GRID PARAMETERS
    //////////////////////////////////////////////////
    double lengthX = 1;  // Length of the domain in the x direction (not including the ghost nodes)
    double lengthY = 1;  // Length of the domain in the y direction (not including the ghost nodes)
    int Nx = 182;        // Number of nodes in the x direction (including ghost nodes)
                         // Note: Specify an even value for parctical reasons.
                         // TODO: Check whether an even value is specified or not.
    if (Nx % 2 != 0) {
        cout << "Nx should be an even number. Please change it." << endl;
        return -1; 
    }
    double h = lengthX / (Nx - 2);   // Uniform grid spacing (h = dx = dy)
    int Ny = (lengthY / h) + 2;      // +2 for ghost nodes
    // Note: If the domain is square, Nx and Ny will be equal
    //       For non square domains, this approach still gives uniform grid spacing in both x and y directions.
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // ITERATION NUMBERS and CONVERGENCE TOLERANCES
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    double maxGSiter = 50;     // Max. number of iterations for the Gauss-Seidel method
    double GSerror;     // Error for the Gauss-Seidel method
    double GStolerance = 1e-5; // Tolerance for the Gauss-Seidel method
    double gsPPEdiff; // Difference in pressure for the Gauss-Seidel method
    int iteration; // Iteration counter for the Gauss-Seidel method

                               // TODO: Change its name to maxGSiter.
    double pChangeLim = 1e-4;  // Convergence tolerance for pressure
    double uChangeLim = 1e-4;  // Convergence tolerance for u velocity
    double vChangeLim = 1e-4;  // Convergence tolerance for v velocity
    double aveChangeU;
    double aveChangeV;
    double aveChangeP;
    double uCenter;

    double errorP; // Error for pressure convergence check
    double errorU; // Error for u velocity convergence check
    double errorV; // Error for v velocity convergence check
    // The average change of u, v and p unknowns are calculated as the sum of the absolute differences between the
    // current and the previous values of all nodes divided by the number of nodes. When this value is less than the
    // specified tolerances, the Gauss-Seidel method (or the overall solution) is considered to be converged.    
    
    
    //////////////////////////////////////////////////////////
    // TIME PARAMETERS OF THE SIMULATION
    //////////////////////////////////////////////////////////
    // Time starts at 0 and ends at endTime.
    // endTime should be set to a large number since we are solving a steady problem and have
    // a convergence criteria. It can be changed to a specific time if needed.
    double startTime = 0;
    double endTime   = 10000;
    
    
    
    //////////////////////////////////////////////////
    // OUTPUT CONTROL PARAMETERS
    //////////////////////////////////////////////////
    double periodOfOutput = 100; // Time period for outputting average change and screen output
    fstream U_outputFile;        // Output file vertical centerline u velocity
    fstream P_outputFile;        // Output file vertical centerline pressure.
    fstream averageChangeFile;   // Output file for average changes of unknonws
    fstream vtkFile;             // Output file in VTK format
    
    averageChangeFile.open("01_average_change.txt", std::ios::out);
    averageChangeFile << "Max. number of Gauss-Seidel iterations: " << maxGSiter << "\n"
                      << "Pressure change limit: "   << pChangeLim << "\n"
                      << "U velocity change limit: " << uChangeLim << "\n"
                      << "V velocity change limit: " << vChangeLim << "\n"
                      << "# Time  U_change  V_change  P_change  uMid\n"; // Header for average change file

    

    //////////////////////////////////////////////////////////
    // TIME STEP SIZE CALCULATION
    //////////////////////////////////////////////////////////
    // CRITERIA 1 ===> h^2 * h^2 / (2 * dynamicViscosity * (h^2 + h^2))
    // CRITERIA 2 ===> 2 * dynamicViscosity / uTopWall^2
    // The time step size is calculated using the minimum of the two criteria
    // Lid driven cavity provblem's lid speed (uTopWall) is used in the second criteria for the maximum reference speed.
    double timeStepSize = min(h*h / (4*dynamicViscosity), 2 * dynamicViscosity / (uTopWall*uTopWall));
    
    //////////////////////////////////////////////////////////
    // ALLOCATING MEMORY FOR THE UNKNOWNS AND INITIALIZATION
    //////////////////////////////////////////////////////////
    cout << "timeStepSize: " << timeStepSize << endl;

    // Initialize u, v and p fields as vectors of vectors of appropriate sizes.
    // Note that the grid is staggered with ghost nodes.
    u.resize(Nx+1, vector<double>(Ny));
    uStar.resize(Nx+1, vector<double>(Ny));
    vStar.resize(Nx, vector<double>(Ny+1));
    v.resize(Nx, vector<double>(Ny+1));
    p.resize(Nx, vector<double>(Ny));

    // TODO: Change the formatting of all the for's, if's while's, etc. in the following way so that the "{" is in the same line as for, if, while, etc. 
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny + 1; j++) {
            v[i][j] = 0;
            vStar[i][j] = 0;
        }
    }
    
    for (int i = 0; i < Nx + 1; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = 0;
            uStar[i][j] = 0;
        }
    }
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            p[i][j] = 0;
        }
    }
    calculateGhostNodeValues(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    calculateGhostNodeValues(2, v, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    calculateGhostNodeValues(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
 

    
    ///////////////////////////////////////////
    // MAIN TIME LOOP
    ///////////////////////////////////////////
    int n = 0; // Time step counter
    for (double t = startTime; t <= endTime; t = t + timeStepSize) {
        
        // PREVIOUS FIELDS ARE STORED FOR CONVERGENCE CHECK
        pPrev = p;
        uPrev = u;
        vPrev = v;

        ///////////////////////////////////////////
        // PREDICTOR STEP
        ///////////////////////////////////////////
        // Predictor step for the u velocity field
        for (int i = 1; i < Nx; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                uAdvection = u[i][j] * (u[i+1][j] - u[i-1][j]) / (2 * h) +
                             0.25 * (v[i-1][j] + v[i][j] + v[i-1][j+1] + v[i][j + 1]) * (u[i][j + 1] - u[i][j - 1]) / (2 * h);
                uDiffuse = (u[i+1][j] + u[i][j - 1] - 4 * u[i][j] + u[i-1][j] + u[i][j + 1]) / (h * h); // d²u/dx² + d²u/dy² (diffusion in X for u)
                uStar[i][j] = u[i][j] + timeStepSize * (kinematicViscosity * uDiffuse - (uAdvection));
            }
        }
        
        // Predictor step for the v velocity field
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny; j++) {
                vAdvection = v[i][j] * (v[i][j + 1] - v[i][j - 1]) / (2 * h) +
                             0.25 * (u[i][j] + u[i][j-1] + u[i+1][j] + u[i+1][j-1]) * (v[i+1][j] - v[i-1][j]) / (2 * h);
                vDiffuse = (v[i+1][j] + v[i][j - 1] - 4 * v[i][j] + v[i-1][j] + v[i][j + 1]) / (h * h); // d²v/dx² + d²v/dy² (diffusion in Y for v)
                vStar[i][j] = v[i][j] + timeStepSize * (kinematicViscosity * (vDiffuse) - (vAdvection));
            }
        }

        calculateGhostNodeValues(1, uStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        calculateGhostNodeValues(2, vStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        
        
        ///////////////////////////////////////////
        // POISSON EQUATION SOLVER
        ///////////////////////////////////////////

        iteration = 0; // Reset iteration counter for the Gauss-Seidel method
        GSerror = 100; // Initialize GSerror to a large value to enter the loop
        while (GSerror>GStolerance&& iteration < maxGSiter) {
            // store old pressures
            pOld = p;


            // Gauss–Seidel update
            for (int i = 1; i < Nx - 1; i++) {
                for (int j = 1; j < Ny - 1; j++) {
                    velocityStarGrad = (uStar[i+1][j] - uStar[i][j] + vStar[i][j+1] - vStar[i][j]) / h;
                    p[i][j] = 0.25 * (p[i+1][j] + p[i-1][j] + p[i][j + 1] + p[i][j - 1]) -
                              (h * h * 0.25 * density / timeStepSize) * (velocityStarGrad); // Update pressure using the Poisson equation
                }
            }


            // // compute L2 norm of change
            GSerror = 0.0;
            for (int i = 1; i < Nx-1; i++) {
                for (int j = 1; j < Ny-1; j++) {
                    gsPPEdiff = abs(p[i][j] - pOld[i][j]);
                    if (gsPPEdiff > GSerror) {
                        GSerror = gsPPEdiff;
                    }
                }
            }

            iteration++;
        }
        ///////////////////////////////////////////
        // ANCHORING THE PRESSURE FIELD
        ///////////////////////////////////////////
        // Modify the pressure field so that the lower left-most pressure value always remains at zero. This is a good practice
        // for the lid driven cavity problem which has no pressure related boundary conditions, therefore has a floating pressure
        // field unless this is done. 
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                p[i][j] = p[i][j] - p[Nx-1][Ny-1];
            }
        }
        calculateGhostNodeValues(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);

        
        ///////////////////////////////////////////
        // CORRECTOR STEP
        ///////////////////////////////////////////
        for (int i = 1; i < Nx; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                u[i][j] = uStar[i][j] - (timeStepSize / density) * ((p[i][j] - p[i-1][j]) / h);
            }
        }
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny; j++) {
                v[i][j] = vStar[i][j] - (timeStepSize / density) * ((p[i][j] - p[i][j-1]) / h);
            }
        }
        calculateGhostNodeValues(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        calculateGhostNodeValues(2, v, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        

        ///////////////////////////////////////////
        // CONVERGENCE CHECK
        ///////////////////////////////////////////
        // Calculate the average change in u, v and p from the previous time step to this one
        aveChangeU = 0.0;
        aveChangeV = 0.0;
        aveChangeP = 0.0;
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                errorP = abs(p[i][j] - pPrev[i][j])/abs(pPrev[i][j]+1e-10);
                if (errorP > aveChangeP) {
                    aveChangeP = errorP;
                }
            }
        }
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                errorU = abs(u[i][j] - uPrev[i][j])/abs(uPrev[i][j]+1e-10);
                if (errorU > aveChangeU) {
                    aveChangeU = errorU;
                }
            }
        }
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                errorV = abs(v[i][j] - vPrev[i][j])/abs(vPrev[i][j]+1e-10);
                if (errorV > aveChangeV) {
                    aveChangeV = errorV;
                }
            }
        }
        // aveChangeU = aveChangeU / (Nx * Ny);
        // aveChangeV = aveChangeV / (Nx * Ny);
        // aveChangeP = aveChangeP / (Nx * Ny);
        uCenter=0.5*(u[(Nx) / 2][(Ny - 2) / 2] + u[(Nx) / 2][(Ny) / 2]);

        // Output the average change values and center u velocity to the console and to the average change file
        if (remainder(n, periodOfOutput) == 0) {
            cout << "Time: " << std::fixed << t;

            cout << std::scientific << "  Uchange: "     << aveChangeU
                                    << "  Vchange: "     << aveChangeV
                                    << "  PressChange: " << aveChangeP;

            // Print Center U velocity in normal (default) notation
            cout << std::fixed << "  Center U velocity: " << uCenter << endl;

            // Write to file (unchanged)
            averageChangeFile << std::fixed << std::setprecision(9) << t << " " << aveChangeU << " " << aveChangeV << " " << aveChangeP << " " << uCenter<< endl;
            
            
            ///////////////////////////////////////////
            // REAL TIME PLOTTING
            ///////////////////////////////////////////
            time_plot.push_back(t);
            velocity_plot.push_back(uCenter);
            plt::clf();
            plt::plot(time_plot, velocity_plot, "b-");
            plt::xlabel("Time");
            plt::ylabel("Center U Velocity");
            plt::title("Real-time Center Velocity");
            plt::grid(true);
            plt::pause(0.001);
            plt::save("real_time_plot.png");

        
        }


        
        // Check for convergence
        // If the average change in u, v and p is less than the specified tolerances, the simulation is considered converged.
        // If the simulation is converged, break the loop and output the final time step size.
        if (aveChangeU < uChangeLim && aveChangeV < vChangeLim && aveChangeP < pChangeLim) {
            cout << "Converged at time: " << t << endl;
            break;
        }

        n++;
    }

    cout << timeStepSize << endl;

    std::cout << "\nEnd of the main function is reached. Stopping.\n\n";
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;
    
    
    
    // TODO: Move all the following file writing code to a new function.
    
    ////////////////////////////////////////////////////////////////
    // WRITE X VELOCITY ON THE VERTICAL CENTERLINE TO A FILE
    ////////////////////////////////////////////////////////////////
    U_outputFile.open("02_U_output.txt", std::ios::out);
    U_outputFile << "nx     = " << Nx << "\n";
    U_outputFile << "ny     = " << Ny << "\n";
    U_outputFile << "dt     = " << timeStepSize << "\n";
    U_outputFile << "Re     = " << Re << "\n";
    U_outputFile << "n      = " << n << "\n";
    U_outputFile << "t      = " << n * timeStepSize << "\n";
    U_outputFile << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms.\n";
    U_outputFile << "y            u \n";
    for (int j = Ny-2; j >= 1; j--) {
        U_outputFile << std::fixed << std::setprecision(7) << (j - 1) * h + h / 2 << "    "
                     << std::fixed << std::setprecision(7) << u[Nx / 2][j] << "\n";
    }

    
    
    ////////////////////////////////////////////////////////////////
    // WRITE PRESSURE ON THE VERTICAL CENTERLINE TO A FILE
    ////////////////////////////////////////////////////////////////
    // This is the average of the pressure at the nodes on the left and the right of the vertical
    // centerline because there is no pressure node directly at the center vertical line

    P_outputFile.open("03_P_output.txt", std::ios::out);
    P_outputFile << "# nx = " << Nx -2 << "\n";
    P_outputFile << "# ny = " << Ny -2 << "\n";
    P_outputFile << "# dt = " << timeStepSize << "\n";
    P_outputFile << "# Re = " << Re << "\n";
    P_outputFile << "# n = " << n << "\n";
    P_outputFile << "# t = " << n * timeStepSize << "\n";
    P_outputFile << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms.\n";
    P_outputFile << "y            p \n";
    for (int j = Ny-2; j >= 1; j--) {
        P_outputFile << std::fixed << std::setprecision(7) << (j - 1) * h + h / 2 << "    "
                     << std::fixed << std::setprecision(7) << 0.5 * (p[(Nx / 2) - 1][j] + p[(Nx / 2)][j]) << "\n";
    }


    ////////////////////////////////////////////////////////////////
    // WRITE THE CURRENT SOLUTION AS A VTK FILE FOR VISUALIZATION
    ////////////////////////////////////////////////////////////////
    vtkFile.open("04_velocityField.txt", std::ios::out);
    vtkFile << "# vtk DataFile Version 2.0\n";
    vtkFile << "Lid Driven Cavity Flow\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << Nx-1 << " " << Ny-1 << " 1\n";
    vtkFile << "POINTS " << (Nx-1) * (Ny-1) << " float\n";
    for (int j = 0; j < Ny - 1; j++) {
        for (int i = 0; i < Nx - 1; i++) {
            vtkFile << std::fixed << std::setprecision(7) << i * h << " "
                    << std::fixed << std::setprecision(7) << j * h << " "
                    << "0.0\n";
        }
    }
    vtkFile << "POINT_DATA " << (Nx-1) * (Ny-1) << "\n";
    vtkFile << "VECTORS velocity float\n";
    for (int j = 0; j < Ny - 1; j++) {
        for (int i = 0; i < Nx - 1; i++) {
            vtkFile << std::fixed << std::setprecision(7) << 0.5*(u[i+1][j+1]+u[i+1][j]) << " "
                    << std::fixed << std::setprecision(7) << 0.5*(v[i][j+1]+v[i+1][j+1]) << " "
                    << "0.0\n";
        }
    }

    U_outputFile.close();
    P_outputFile.close();
    averageChangeFile.close();
    vtkFile.close();
    


    return 0;
}  // End of the main function



// TODO: I do not like this function. b and M variables are hard to follow. It only works for the LDC problem I guess, but it is still too complicated.
//       We have to think about it and change it.
//       We may have boundary types such as wall, inlet, outlet, etc. and associate the boundaries of our problem to these
//       types and the code should automatically take care of ghost nodes. For Step 01, we need to implement only the wall type.
void calculateGhostNodeValues(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny)
/* This is a function to set ghost node values based on specified boundary conditions.
   b is the flow field type. 0: pressure, 1: u velocity, 2: velocity
*/
{
    if (b == 0) {
        for (int i = 1; i < Nx-1; i++) {
            M[i][0] = M[i][1];            // Bottom wall, dp/dx = 0
            M[i][Ny - 1] = M[i][Ny - 2];  // Top wall, dp/dx = 0
        }
        for (int j = 1; j < Nx-1; j++) {
            M[0][j] = M[1][j];            // Left wall, dp/dy = 0
            M[Nx - 1][j] = M[Nx - 2][j];  // Right wall, dp/dy = 0
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                                // Top left corner
        M[0][Ny - 1] = 0.5 * (M[0][Ny - 2] + M[1][Ny - 1]);                 // Top right corner
        M[Nx - 1][0] = 0.5 * (M[Nx - 2][0] + M[Nx - 1][1]);                 // Bottom left corner
        M[Nx - 1][Ny - 1] = 0.5 * (M[Nx - 1][Ny - 2] + M[Nx - 2][Ny - 1]);  // Bottom right corner
    }
    else if (b == 1) {
        for (int i = 1; i < Nx-1; i++) {
            M[i][0] = 2 * uBottomWall - M[i][1];
            M[i][Ny - 1] = 2*uTopWall - M[i][Ny-2];  // Right wall, u = 0
        }
        for (int j = 1; j <= Ny-1; j++) {
            M[0][j] = 2 * uLeftWall - M[1][j];  // Top wall, ghost node
            M[1][j] = uLeftWall;
            M[Nx-1][j] = uRightWall;
            M[Nx][j] = 2 * uRightWall - M[Nx - 2][j];  // Bottom wall, ghost node
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                               // Top left corner
        M[0][Ny] = 0.5 * (M[0][Ny - 1] + M[1][Ny]);                        // Top right corner
        M[Nx - 1][0] = 0.5 * (M[Nx - 2][0] + M[Nx - 1][1]);                // Bottom left corner
        M[Nx - 1][Ny - 1] = 0.5 * (M[Nx - 1][Ny - 2] + M[Nx - 2][Ny - 1]); // Bottom right corner
    }
    else if (b == 2) {
        for (int j = 1; j <= Ny-1; j++) {
            M[0][j] = 2 * vLeftWall - M[1][j];         // Top wall, v = 0
            M[Nx - 1][j] = 2*vRightWall - M[Nx-2][j];  // Bottom wall, ghost node, v = 0
        }
        
        for (int i = 1; i < Nx-1; i++) {
            M[i][1] = vBottomWall;
            M[i][0] = 2 * vBottomWall - M[i][2];  // Left wall, ghost node
            
            M[i][Ny - 1] = vTopWall;              // Right wall, ghost node
            M[i][Ny] = 2*vTopWall - M[i][Ny-2];
        }
        
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                  // Top left corner
        M[0][Ny] = 0.5 * (M[0][Ny - 1] + M[1][Ny]);           // Top right corner
        M[Nx-1][0] = 0.5 * (M[Nx - 2][0] + M[Nx-1][1]);       // Bottom left corner
        M[Nx-1][Ny] = 0.5 * (M[Nx-2][Ny] + M[Nx - 1][Ny-1]);  // Bottom right corner
    }
}  // End of function calculateGhostNodeValues
