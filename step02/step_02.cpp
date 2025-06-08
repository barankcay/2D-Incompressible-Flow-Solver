// g++ fdm_stag_uniform.cpp -o fdm_stag_uniform.exe -O3 -ffast-math

// Ibrahim Baran KUCUKCAY
//  Departmant of Mechanical Engineering, Middle East Technical University, Ankara, Turkey

#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting
#include <cmath>
#include <thread>
#include <chrono>
using namespace std;
// Case Summary:
//  This code is a 2D staggered grid finite volume method for simulating lid driven cavity flow.
//  The solution method of Navier Stokes equations is fractional step or predictor corrector method.
//  The code is based on uniform grid spacing in both x and y directions.
//  X direction velocities, U's, are placed at the left face of the cells.
//  Y direction velocities, V's, are placed to the bottom face of the cells.

// The time discretization is done using explicit Euler method, which is first order accurate.
// The spatial discretization is done using central difference scheme, which is second order accurate.
// The Poisson equation for pressure is solved using Gauss-Seidel method.

// Average change of the parameters is calculated to check for Gauss-Seidel convergence and also to check for convergence of the simulation.
// The average change limits for convergence are set for pressure, u velocity and v velocity.

// The boundary conditions are implemented using ghost cells.
// The west, east and south walls are no slip walls.
// The north wall is a moving wall with velocity U = 1, V = 0.
// All walls have zero pressure gradient.

//!!!!!!!!!!!!!!!!!!
// Why indexing is j,i instead of i,j?
//  This is because matrix indexing in C++ is done in row-major order.
//  Our 2D domain is represented as matrix with rows and columns.
//  The first index represents the row (y direction) and the second index represents the column (x direction).
//  To move in y direction, we change the index j. And to move in y direction for matrix, we change the first index.
//  That makes our first index j and second index i.

//-------------------------------------------------------------------------//

////// This is a function to set ghost cell values based on boundary conditions
////// int b is the flow field type:
////// 0 - pressure, 1 - u velocity, 2 - v velocity
////// For lid driven cavity, ALL OF THE DIRICHLET BOUNDARY CONDITIONS ARE VELOCITIES
////// So, specific input of velocity values are required for function call
////// Pressure boundary conditions are all Neumann, which means there is not a specific pressure value at the boundaries
void calculateGhostCellValues(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny);
int main()
{
    auto start = std::chrono::steady_clock::now();

    //////////////////////////////////////////////////
    ////////// CREATION OF FIELD PARAMETERS //////////
    //////////////////////////////////////////////////
    vector<vector<double>> u;     // u velocity field
    vector<vector<double>> uStar; // intermediate u velocity field
    vector<vector<double>> v;     // v velocity field
    vector<vector<double>> vStar; // intermediate v velocity field
    vector<vector<double>> p;     // pressure field
    vector<vector<double>> pPrev; // previous pressure field for convergence check
    vector<vector<double>> pOld;  // old pressure field for Gauss-Seidel method

    vector<vector<double>> uPrev; // previous u velocity field for convergence check
    vector<vector<double>> vPrev; // previous v velocity field for convergence check

    double un;
    double vn;
    double us;
    double vs;
    double ue;
    double ve;
    double uw;
    double vw;

    double uAdvection; // Total advection term for u. v * du/dy + u * du/dx
    double vAdvection; // Total advection term for v. u * dv/dx + v * dv/dy

    double uDiffuse; // d²u/dx² + d²u/dy² (diffusion in X and Y for u)
    double vDiffuse; // d²v/dx² + d²v/dy² (diffusion in X and Y for v)

    double velocityStarGrad; // Gradient of uStar and vStar in X and Y directions for pressure Poisson equation. du*/dx + dv*/dy

    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// CHARACTERISTICS OF THE FLOW ///////////
    //////////////////////////////////////////////////
    double Re = 1000;
    double density = 1.0;
    double kinematicViscosity = 1.0 / Re;
    double dynamicViscosity = kinematicViscosity * density;

    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// BOUNDARY CONDITIONS //////////////////
    //////////////////////////////////////////////////
    double uTopWall = 1.0;    // x direction velocity at the top wall
    double uBottomWall = 0.0; // x direction velocity at the bottom wall
    double uLeftWall = 0.0;   // x direction velocity at the left wall
    double uRightWall = 0.0;  // x direction velocity at the right wall
    double vTopWall = 0.0;    // y direction velocity at the top wall
    double vBottomWall = 0.0; // y direction velocity at the bottom wall
    double vLeftWall = 0.0;   // y direction velocity at the left wall
    double vRightWall = 0.0;  // y direction velocity at the right wall
    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// GRID PARAMETERS ///////////////////////
    //////////////////////////////////////////////////
    // lengthx and lengthy are the lengths of the domain in the x and y directions, not including ghost cells
    double lengthX = 1; // Length of the domain in the x direction
    double lengthY = 1; // Length of the domain in the y direction
    // Nx and Ny are the number of cells in the x and y directions, including ghost cells
    int Nx = 130;
    double h = lengthX / (Nx - 2);
    int Ny = (lengthY / h) + 2; // +2 for ghost cells
    //!!!!!! Since uniform grid spacing is used, this way of calculating Ny is still valid.
    //!!!!!! If the domain is square, Nx and Ny will be equal
    //!!!!!! For non square domains, this approach still gives uniform grid spacing in both x and y directions.
    //--------------END-OF-THE-SECTION------------------------//

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////// AVERAGE CHANGE LIMITS OF PARAMETERS OF THE SIMULATION AND PRESSURE POISSON EQ  /////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    double gsChangeLim = 1e-3; // Gauss-Seidel convergence criteria for pressure Poisson equation
    double pChangeLim = 1e-9; // Average change limit for pressure
    double uChangeLim = 1e-9; // Average change limit for u velocity
    double vChangeLim = 1e-9; // Average change limit for v velocity

    double gsChange; // Average change in pressure for Gauss-Seidel method
    // The average change is calculated as the sum of the absolute differences between the current and previous values divided by the number of cells
    // When this value is less than the specified limit, the Gauss-Seidel method is considered converged.

    int iteration;                // Counter for the number of iterations in the Gauss-Seidel method
    int gsIterationLimit = 10000; // Maximum number of iterations for the Gauss-Seidel method
    // If the Gauss-Seidel method does not converge within this number of iterations, the simulation will stop.
    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// OUTPUT CONTROL PARAMETERS /////////////
    //////////////////////////////////////////////////
    double periodOfOutput = 100; // Time period for outputting average change and screen output
    fstream U_outputFileStag;        // Output file vertical centerline u velocity
    fstream P_outputFileStag;        // Output file vertical centerline pressure.
    fstream averageChangeFile;   // Output file for average change values
    averageChangeFile.open("01_average_change.txt", std::ios::out);
    averageChangeFile <<  "Gauss-Seidel change limit: " << gsChangeLim << "\n"
                      << "Pressure change limit: " << pChangeLim << "\n"
                      << "U velocity change limit: " << uChangeLim << "\n"
                      << "V velocity change limit: " << vChangeLim << "\n"
                      << "# Time  U_change  V_change  P_change  uMid\n"; // Header for average change file
    //////////// END OF THE OUTPUT CONTRO PARAMETERS SECTION ////////////

    //////////////////////////////////////////////////////////
    //////////// TIME PARAMETERS OF THE SIMULATION////////////
    //////////////////////////////////////////////////////////
    // Time starts at 0 and ends at endTime
    double startTime = 0;
    // EndTime should be set to a large number since we have a convergence criteria
    // This can be changed to a specific time if needed
    double endTime = 10000;
    //////////// END OF THE TIME CONTROLS SECTION ////////////

    //////////////////////////////////////////////////////////
    //////////// TIME STEP SIZE CALCULATION///////////////////
    //////////////////////////////////////////////////////////
    // CRITERIA 1 ===> h^2 * h^2 / (2 * dynamicViscosity * (h^2 + h^2))
    // CRITERIA 2 ===> 2 * dynamicViscosity / (uTopWall^2 + vLeftWall^2)
    // The time step size is calculated using the minimum of the two criteria
    // Lid driven cavity is constructed with 3 walls.
    // So, velocity boundary conditions are applicable only to the top wall.
    // For time step size calculation, we consider the top wall velocities only.
    double timeStepSize = min((h * h) / (4 * dynamicViscosity), 2 * dynamicViscosity / (uTopWall * uTopWall + vTopWall * vTopWall));
    //////////// END OF THE TIME STEP SIZE CALCULATION ////////////
    //////////////////////////////////////////////////////////////////////////////
    //////////// INITIALIZATION and ASSIGNING DIMENSION OF THE FIELDS ////////////
    //////////////////////////////////////////////////////////////////////////////
    // Initialize the fields and assign number of cells in each direction
    // To represent ghost cells, adding one extra cell in each direction is needed.
    // u and v fields are have one extra set (not cell) in the x and y directions respectively.
    // reason is to place velocity field around domain other than placing on the boundary.
    // without this approach, for example, far right u velocity point would be on the boundary.
    u.resize(Nx + 1, vector<double>(Ny));
    uStar.resize(Nx + 1, vector<double>(Ny));
    vStar.resize(Nx, vector<double>(Ny + 1));
    v.resize(Nx, vector<double>(Ny+1));
    p.resize(Nx, vector<double>(Ny));

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {

            v[i][j] = 0;

            vStar[i][j] = 0;
        }
    }
    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            u[i][j] = 0;

            uStar[i][j] = 0;
        }
    }
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {

            p[i][j] = 0;
        }
    }
    //////////// END OF THE INITIALIZATION ////////////

    //////////// BOUNDARY CONDITIONS////////////
    calculateGhostCellValues(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    calculateGhostCellValues(2, v, vLeftWall, vRightWall, uTopWall, uBottomWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    calculateGhostCellValues(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    //////////// END OF THE INITIAL BOUNDARY CONDITION ASSIGN SECTION ////////////

    int n = 0; // counter for the time steps
    for (double t = startTime; t <= endTime; t = t + timeStepSize)
    {

        //////////// PREVIOUS FIELDS ARE STORED FOR CONVERGENCE CHECK ///////////
        pPrev = p;
        uPrev = u;
        vPrev = v;
        //////// PREDICTOR STEP //////
        // Predictor step for the u velocity field
        for (int j = 1; j < Ny - 1; j++)
        {
            for (int i = 1; i < Nx; i++)
            {
                vn = (v[i-1][j+1] + v[i][j+1]) / 2; // v velocity field
                un = (u[i][j] + u[i][j+1]) / 2;
                ue = (u[i][j] + u[i+1][j]) / 2; // u velocity field
                uw = (u[i-1][j] + u[i][j]) / 2; // u velocity field
                vs = (v[i-1][j] + v[i][j]) / 2;
                us = (u[i][j-1] + u[i][j]) / 2;

                uAdvection = h * (vn * un - vs * us + ue * ue - uw * uw);                                              // Total advection term for u. v * du/dy + u * du/dx
                uDiffuse = kinematicViscosity * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4 * u[i][j]); // diffusion term
                uStar[i][j] = u[i][j] + (timeStepSize / (h * h)) * (-uAdvection + uDiffuse);                           // u star velocity field
            }
        }
        // Predictor step for the v velocity field
        for (int j = 1; j < Ny; j++)
        {
            for (int i = 1; i < Nx - 1; i++)
            {

                vn = (v[i][j] + v[i][j+1]) / 2;
                uw = (u[i][j-1] + u[i][j]) / 2;
                vw = (v[i-1][j] + v[i][j]) / 2;
                vs = (v[i][j-1] + v[i][j]) / 2;
                ue = (u[i+1][j-1] + u[i+1][j]) / 2;
                ve = (v[i][j] + v[i+1][j]) / 2;
                vAdvection = h * (vn * vn - vs * vs + ue * ve - uw * vw);                                              // Total advection term for v. u * dv/dx + v * dv/dy
                vDiffuse = kinematicViscosity * (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] - 4 * v[i][j]); // diffusion term
                vStar[i][j] = v[i][j] + (timeStepSize / (h * h)) * (-vAdvection + vDiffuse);                           // v velocity field
            }
        }
        calculateGhostCellValues(1, uStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        calculateGhostCellValues(2, vStar, vLeftWall, vRightWall, uTopWall, uBottomWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        
        
        
        ////////// POISSON EQUATION SOLVER //////////
        gsChange = 1.0; // This is set to 1 to enter the while loop
        iteration = 0;  // Counter for the number of iterations in the Gauss-Seidel method
        while (gsChange > gsChangeLim)
        {
            pOld = p; // Store the old pressure values for convergence check
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    velocityStarGrad = (density / timeStepSize) * ((uStar[i+1][j] - uStar[i][j] + vStar[i][j+1] - vStar[i][j]) / h);
                    p[i][j] = ((p[i-1][j] + p[i+1][j] + p[i][j+1] + p[i][j-1]) - velocityStarGrad * h * h) * 0.25; // Pressure Poisson equation
                }
            }

            gsChange = 0.0; // Reset residual for each iteration
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    // Calculate the total change in pressure for convergence check
                    gsChange = gsChange + abs(p[i][j] - pOld[i][j]);
                }
            }
            gsChange = gsChange / (Nx * Ny); // Average change in pressure

            if (iteration > gsIterationLimit)
            {
                cout << "Gauss-Seidel method did not converge within the maximum number of iterations." << endl;
                cout << "Exiting the simulation." << endl;
                exit(1); // Exit the program if Gauss-Seidel method does not converge
            }
            iteration = iteration + 1; // Increment the iteration counter
        }
        ////// END OF POISSON EQUATION SOLVER ///////

        /////////// ANCHORING THE PRESSURE FIELD ////////////
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                p[i][j] = p[i][j] - p[Ny-2][1];
            }
        }
        //////////// END OF ANCHORING THE PRESSURE FIELD ////////////

        //////////// GHOST CELL VALUES ARE CALCULATED FOR THE PRESSURE FIELD ////////////
        calculateGhostCellValues(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);



        ////////////////////////////////////////
        //////////// CORRECTOR STEP ////////////
        ////////////////////////////////////////
        for (int i = 1; i < Nx; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                u[i][j] = uStar[i][j] - (timeStepSize / (h)) * ((p[i][j] - p[i-1][j]) / density); // Corrector step for u velocity field
            }
        }
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny; j++)
            {
                v[i][j] = vStar[i][j] - (timeStepSize / (h)) * ((p[i][j] - p[i][j-1]) / density); // Corrector step for v velocity field
            }
        }
        calculateGhostCellValues(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        calculateGhostCellValues(2, v, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        
        
        
        //////////// CONVERGENCE CHECK ////////////
        // Calculate the average change in u, v and p fields to check for convergence
        double aveChangeU = 0.0;
        double aveChangeV = 0.0;
        double aveChangeP = 0.0;
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                aveChangeP = aveChangeP + abs(p[i][j] - pPrev[i][j]);
            }
        }
        for (int i = 1; i < Nx; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                aveChangeU = aveChangeU + abs(u[i][j] - uPrev[i][j]);
            }
        }
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny; j++)
            {
                aveChangeV = aveChangeV + abs(v[i][j] - vPrev[i][j]);
            }
        }        
        aveChangeU = aveChangeU / (Nx * Ny);
        aveChangeV = aveChangeV / (Nx * Ny);
        aveChangeP = aveChangeP / (Nx * Ny);
        //////////// END OF CONVERGENCE CHECK ////////////


        // Output the average change values and center u velocity to the console and average change file
        if (remainder(n, periodOfOutput) == 0)
        {
            // Print Time in normal (default) notation
            cout << "Time: " << std::fixed << t;

            // Print changes in scientific notation
            cout << std::scientific << "  Uchange: " << aveChangeU
                 << "  Vchange: " << aveChangeV
                 << "  PressChange: " << aveChangeP;

            // Print Center U velocity in normal (default) notation
            cout << std::fixed << "  Center U velocity: " << u[(Ny - 1) / 2][(Nx - 1) / 2] << endl;

            // Write to file (unchanged)
            averageChangeFile << std::fixed << t << " " << log(aveChangeU) << " " << log(aveChangeV) << " " << log(aveChangeP) << " " << u[(Ny - 1) / 2][(Nx - 1) / 2] << endl;
        } // Check for convergence.
        // If the average change in u, v and p is less than the specified limits, the simulation is considered converged.
        // If the simulation is converged, break the loop and output the final time step size.
        if (aveChangeU < uChangeLim && aveChangeV < vChangeLim && aveChangeP < pChangeLim)
        {
            cout << "Converged at time: " << t << endl;
            break;
        }
        n++;
    }

    cout << timeStepSize << endl;

    std::cout << "\nEnd of the main function is reached. Stopping.\n\n";
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;
    //////////////////////////////////////////////////////////////////////////////
    //////////// OUTPUTTING VERTICAL X VELOCITY TO OUTPUT FILE ///////////////////
    //////////////////////////////////////////////////////////////////////////////
    U_outputFileStag.open("02_U_output.txt", std::ios::out);
    U_outputFileStag<< "nx     = " << Nx << "\n";
    U_outputFileStag<< "ny     = " << Ny << "\n";
    U_outputFileStag<< "dt     = " << timeStepSize << "\n";
    U_outputFileStag<< "Re     = " << Re << "\n";
    U_outputFileStag<< "n      = " << n << "\n";
    U_outputFileStag<< "t      = " << n * timeStepSize << "\n";
    U_outputFileStag<< "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms.\n";
    U_outputFileStag<< "y            u \n";
    for (int j = Ny-2; j >= 1; j--)
    {
        U_outputFileStag<< std::fixed << std::setprecision(7) << (j - 1) * h + h / 2 << "    "
                     << std::fixed << std::setprecision(7) << u[Nx / 2][j] << "\n";
    }

    //////////////////////////////////////////////////////////////////////////////
    //////////// OUTPUTTING VERTICAL PRESSURE TO OUTPUT FILE /////////////////////
    //////////// This is the average of the pressure at the left and right nodes//
    //////////// because there is no pressure node at the center vertical line of the domain //
    //////////////////////////////////////////////////////////////////////////////
    P_outputFileStag.open("03_P_output.txt", std::ios::out);
    P_outputFileStag<< "# nx = " << Nx << "\n";
    P_outputFileStag<< "# ny = " << Ny << "\n";
    P_outputFileStag<< "# dt = " << timeStepSize << "\n";
    P_outputFileStag<< "# Re = " << Re << "\n";
    P_outputFileStag<< "# n = " << n << "\n";
    P_outputFileStag<< "# t = " << n * timeStepSize << "\n";
    P_outputFileStag<< "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms.\n";
    P_outputFileStag<< "y            p \n";
    for (int j = Ny-2; j >= 1; j--)
    {
        P_outputFileStag<< std::fixed << std::setprecision(7) << (j - 1) * h + h / 2 << "    "
                     << std::fixed << std::setprecision(7) << 0.5 * (p[(Nx / 2) - 1][j] + p[(Nx / 2) + 1][j]) << "\n";
    }

    U_outputFileStag.close();
    P_outputFileStag.close();
    averageChangeFile.close();

    return 0;
}

void calculateGhostCellValues(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Ny, int Nx)
{
    if (b == 0)
    {
        for (int i = 1; i < Nx-1; i++)
        {
            M[i][0] = M[i][1];           // Bottom wall, dp/dx = 0
            M[i][Ny - 1] = M[i][Ny - 2]; // Top wall, dp/dx = 0
        }
        for (int j = 1; j < Nx-1; j++)
        {
            M[0][j] = M[1][j];           // Left wall, dp/dy = 0
            M[Nx - 1][j] = M[Nx - 2][j]; // Right wall, dp/dy = 0
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                               // Top left corner
        M[0][Ny - 1] = 0.5 * (M[0][Ny - 2] + M[1][Ny - 1]);                // Top right corner
        M[Nx - 1][0] = 0.5 * (M[Nx - 2][0] + M[Nx - 1][1]);                // Bottom left corner
        M[Nx - 1][Ny - 1] = 0.5 * (M[Nx - 1][Ny - 2] + M[Nx - 2][Ny - 1]); // Bottom right corner
    }
    else if (b == 1)
    {
        
        for (int i = 1; i < Nx-1; i++)
        {
        
            M[i][0] = 2 * uBottomWall - M[i][1];
            M[i][Ny - 1] =2*uTopWall- M[i][Ny-2]; // Right wall, u = 0

        }
        for (int j = 1; j <= Ny-1; j++)
        {
            M[0][j] = 2 * uLeftWall - M[1][j];              // Top wall, ghost cell
            M[1][j]=uLeftWall;
            M[Nx-1][j]=uRightWall;
            M[Nx][j] = 2 * uRightWall - M[Nx - 2][j]; // Bottom wall, ghost cell
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                               // Top left corner
        M[0][Ny-1] = 0.5 * (M[0][Ny - 2] + M[1][Ny-1]);                        // Top right corner
        M[Nx][0] = 0.5 * (M[Nx - 1][0] + M[Nx][1]);                // Bottom left corner
        M[Nx][Ny - 1] = 0.5 * (M[Nx - 1][Ny - 1] + M[Nx][Ny - 2]); // Bottom right corner
    }
    else if (b == 2)
    {
        for (int j = 1; j <= Ny-1; j++)
        {
            
            M[0][j] = 2 * vLeftWall - M[1][j]; // Top wall, v = 0

            M[Nx - 1][j] = 2*vRightWall-M[Nx-2][j]; // Bottom wall, ghost cell, v = 0

        }
        
        for (int i = 1; i < Nx-1; i++)
        {
            
            M[i][1]=vBottomWall;
            M[i][0] = 2 * vBottomWall - M[i][2];            // Left wall, ghost cell
            
            M[i][Ny - 1] = vTopWall; // Right wall, ghost cell
            M[i][Ny]=2*vTopWall-M[i][Ny-2];
        }
        
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                // Top left corner
        M[0][Ny] = 0.5 * (M[0][Ny - 1] + M[1][Ny]); // Top right corner
        M[Nx-1][0] = 0.5 * (M[Nx - 2][0] + M[Nx-1][1]);         // Bottom left corner
        M[Nx-1][Ny] = 0.5 * (M[Nx-2][Ny] + M[Nx - 1][Ny-1]);  // Bottom right corner
    }
}