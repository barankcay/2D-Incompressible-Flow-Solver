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

//Case Summary:
// This code is a 2D staggered grid finite difference method for simulating lid driven cavity flow.
// The solution method of Navier Stokes equations is done using fractional step or predictor corrector method.
// The code is based on uniform grid spacing in both x and y directions.
// X direction velocities, U's, are placed to the left of the pressure nodes.
// Y direction velocities, V's, are placed to the bottom of the pressure nodes.

// The time discretization is done using explicit Euler method, which is first order accurate.
// The spatial discretization is done using central difference scheme, which is second order accurate.
// The Poisson equation for pressure is solved using Gauss-Seidel method.

// Average change of the parameters is calculated to check for Gauss-Seidel convergence and also to check for convergence of the simulation.
// The average change limits for convergence are set for pressure, u velocity and v velocity.

// The boundary conditions are implemented using ghost nodes.
// The west, east and south walls are no slip walls.
// The north wall is a moving wall with velocity U = 1, V = 0.
// All walls have zero pressure gradient.


//-------------------------------------------------------------------------//

////// Function to set ghost node values based on boundary conditions
////// int b is the flow field type:
////// 0 - pressure, 1 - u velocity, 2 - v velocity
////// For lid driven cavity, ALL OF THE DIRICHLET BOUNDARY CONDITIONS ARE VELOCITIES
////// So, input of velocity values are required for function call
////// Pressure boundary conditions are all Neumann, which means there is not a specific pressure value at the boundaries
void ghostNodeValues(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny);

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
    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// CHARACTERISTICS OF THE FLOW ///////////
    //////////////////////////////////////////////////
    double Re                   = 5000;
    double density              = 1.0;
    double kinematicViscosity   = 1.0 / Re;
    double dynamicViscosity     = kinematicViscosity * density;

    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// BOUNDARY CONDITIONS //////////////////
    //////////////////////////////////////////////////
    double uTopWall     = 1.0; // x direction velocity at the top wall
    double uBottomWall  = 0.0; // x direction velocity at the bottom wall
    double uLeftWall    = 0.0; // x direction velocity at the left wall
    double uRightWall   = 0.0; // x direction velocity at the right wall
    double vTopWall     = 0.0; // y direction velocity at the top wall
    double vBottomWall  = 0.0; // y direction velocity at the bottom wall
    double vLeftWall    = 0.0; // y direction velocity at the left wall
    double vRightWall   = 0.0; // y direction velocity at the right wall
    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// GRID PARAMETERS ///////////////////////
    //////////////////////////////////////////////////
    // lengthx and lengthy are the lengths of the domain in the x and y directions, not including ghost nodes
    double lengthX  = 1; // Length of the domain in the x direction
    double lengthY  = 1; // Length of the domain in the y direction
    // Nx and Ny are the number of cells in the x and y directions, including ghost cells
    int Nx          = 172; 
    double h        = lengthX / (Nx - 2);
    int Ny          = (lengthY / h) + 2; // +2 for ghost cells
    //!!!!!! Since uniform grid spacing is used, this way of calculating Ny is still valid.
    //!!!!!! If the domain is square, Nx and Ny will be equal
    //!!!!!! For non square domains, this approach still gives uniform grid spacing in both x and y directions.
    //--------------END-OF-THE-SECTION------------------------//

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////// AVERAGE CHANGE LIMITS OF PARAMETERS OF THE SIMULATION AND PRESSURE POISSON EQ  /////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    double GSchangeLim  = 1e-3; // Gauss-Seidel convergence criteria for pressure Poisson equation
    double pChangeLim   = 1e-9; // Average change limit for pressure
    double uChangeLim   = 1e-9; // Average change limit for u velocity
    double vChangeLim   = 1e-9; // Average change limit for v velocity
    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////
    ////////// OUTPUT CONTROL PARAMETERS /////////////
    //////////////////////////////////////////////////
    double periodOfOutput = 100; // Time period for outputting average change and screen output
    fstream U_outputFile;        // Output file vertical centerline u velocity
    fstream P_outputFile;        // Output file vertical centerline pressure.
    fstream averageChangeFile;   // Output file for average change values
    averageChangeFile.open("01_average_change.txt", std::ios::out);
    averageChangeFile << "Gauss-Seidel change limit: " << GSchangeLim << "\n";
    averageChangeFile << "Pressure change limit: " << pChangeLim << "\n";
    averageChangeFile << "U velocity change limit: " << uChangeLim << "\n";
    averageChangeFile << "V velocity change limit: " << vChangeLim << "\n";
    averageChangeFile << "# Time U_change V_change P_change uMid\n"; // Header for average change file
    //--------------END-OF-THE-SECTION------------------------//

    //////////////////////////////////////////////////////////
    //////////// TIME PARAMETERS OF THE SIMULATION////////////
    //////////////////////////////////////////////////////////
    // Time starts at 0 and ends at endTime
    double startTime    = 0; 
    // EndTime should be set to a large number since we have a convergence criteria
    // This can be changed to a specific time if needed
    double endTime      = 10000;
    //--------------END-OF-THE-SECTION------------------------//



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
    //--------------END-OF-THE-SECTION------------------------//



    //////////////////////////////////////////////////////////////////////////////
    //////////// INITIALIZATION and ASSIGNING DIMENSION OF THE FIELDS ////////////
    //////////////////////////////////////////////////////////////////////////////
    // Initialize the fields and assign number of nodes in each direction
    // adding one extra node in each direction. u and v fields are have one extra node in the x and y directions respectively.
    // reason is to place velocity field around domain other than placing on the boundary.
    // without this approach, for example, far right u velocity node would be on the boundary.
    u.resize(Ny, vector<double>(Nx + 1));
    uStar.resize(Ny, vector<double>(Nx + 1));
    vStar.resize(Ny + 1, vector<double>(Nx));
    v.resize(Ny + 1, vector<double>(Nx));
    p.resize(Ny, vector<double>(Nx));

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {

            v[j][i] = 0;

            vStar[j][i] = 0;
        }
    }
    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            u[j][i] = 0;

            uStar[j][i] = 0;
        }
    }
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {

            p[j][i] = 0;
        }
    }
    //////////// END OF INITIALIZATION ////////////

    
    //////////// BOUNDARY CONDITIONS////////////
    ghostNodeValues(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    ghostNodeValues(2, v, vLeftWall, vRightWall, uTopWall, uBottomWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    ghostNodeValues(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    //////////// END OF BOUNDARY CONDITIONS////////////

    int n = 0; // counter for the time steps
    for (double t = startTime; t <= endTime; t = t + timeStepSize)
    {

        //////////// PREVIOUS FIELDS ARE STORED FOR CONVERGENCE CHECK ///////////
        vector<vector<double>> pPrev = p;
        vector<vector<double>> uPrev = u;
        vector<vector<double>> vPrev = v;

        //////// PREDICTOR STEP //////
        for (int i = 1; i < Nx; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                double uAdvectX     = u[j][i] * (u[j][i + 1] - u[j][i - 1]) / (2 * h);                                                        // u * du/dx (advection in X for u)
                double uAdvectY     = 0.25 * (v[j][i - 1] + v[j][i] + v[j + 1][i - 1] + v[j + 1][i]) * (u[j - 1][i] - u[j + 1][i]) / (2 * h); // v * du/dy (advection in Y for u)
                double uDiffuseX    = (u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / (h * h);                                                   // d²u/dx² (diffusion in X for u)
                double uDiffuseY    = (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / (h * h);                                                   // d²u/dy² (diffusion in Y for u)
                uStar[j][i]         = u[j][i] + timeStepSize * (kinematicViscosity * (uDiffuseX + uDiffuseY) - (uAdvectX + uAdvectY));
            }
        }

        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny; j++)
            {
                double vAdvectX     = 0.25 * (u[j - 1][i] + u[j][i] + u[j - 1][i + 1] + u[j][i + 1]) * (v[j][i + 1] - v[j][i - 1]) / (2 * h); // u * dv/dx (advection in X for v)
                double vAdvectY     = v[j][i] * (v[j - 1][i] - v[j + 1][i]) / (2 * h);                                                        // v * dv/dy (advection in Y for v)
                double vDiffuseX    = (v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / (h * h);                                                   // d²v/dx² (diffusion in X for v)
                double vDiffuseY    = (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / (h * h);                                                   // d²v/dy² (diffusion in Y for v)
                vStar[j][i]         = v[j][i] + timeStepSize * (kinematicViscosity * (vDiffuseX + vDiffuseY) - (vAdvectX + vAdvectY));
            }
        }

        ghostNodeValues(1, uStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        ghostNodeValues(2, vStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);

        ///////// END OF PREDICTOR STEP /////////


        ////////// POISSON EQUATION SOLVER //////////
        double GSchange = 1.0; // This is set to 1 to enter the while loop
        int iteration   = 0;
        while (GSchange > GSchangeLim)
        {
            vector<vector<double>> pOld = p; // Store the old pressure values for convergence check
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    double uStarGradX   = (uStar[j][i + 1] - uStar[j][i]) / (h); // Gradient of uStar in X direction
                    double vStarGradY   = (vStar[j][i] - vStar[j + 1][i]) / (h); // Gradient of vStar in Y direction
                    p[j][i]             = 0.25 * (p[j][i + 1] + p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) - (h * h * 0.25 * density / timeStepSize) * (uStarGradX + vStarGradY); // Update pressure using the Poisson equation
                }
            }
            GSchange = 0.0; // Reset residual for each iteration
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    // Calculate the total change in pressure for convergence check
                    GSchange = GSchange + abs(p[j][i] - pOld[j][i]);
                }
            }
            GSchange = GSchange / (Nx * Ny); // Average change in pressure
            iteration++;
            ghostNodeValues(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        }
        ////// END OF POISSON EQUATION SOLVER ///////

        //////////// CORRECTOR STEP ////////////
        for (int i = 1; i < Nx; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                u[j][i] = uStar[j][i] - (timeStepSize / density) * ((p[j][i] - p[j][i - 1]) / (h));
            }
        }
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny; j++)
            {
                v[j][i] = vStar[j][i] - (timeStepSize / density) * ((p[j - 1][i] - p[j][i]) / (h));
            }
        }
        ghostNodeValues(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        ghostNodeValues(2, v, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        //////////// FINAL VELOCITY AND PRESSURE VALUES ARE OBTAINED. END OF CORRECTOR STEP ////////////

        //////////// CONVERGENCE CHECK ////////////
        // Calculate the average change in u, v and p fields to check for convergence
        double aveChangeU = 0.0;
        double aveChangeV = 0.0;
        double aveChangeP = 0.0;
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                aveChangeU = aveChangeU + abs(u[j][i] - uPrev[j][i]);
                aveChangeV = aveChangeV + abs(v[j][i] - vPrev[j][i]);
                aveChangeP = aveChangeP + abs(p[j][i] - pPrev[j][i]);
            }
        }
        aveChangeU = aveChangeU / (Nx * Ny);
        aveChangeV = aveChangeV / (Nx * Ny);
        aveChangeP = aveChangeP / (Nx * Ny);
        //////////// END OF CONVERGENCE CHECK ////////////
        
        // Output the average change values and center u velocity to the console and average change file
        if (remainder(n, periodOfOutput) == 0)
        {
            cout << "Time:" << t << "   Uchange: " << aveChangeU << " Vchange: " << aveChangeV << " PressChange: " << aveChangeP << endl;
            cout << "Center U velocity: " << u[(Ny - 1) / 2][(Nx - 1) / 2] << endl;
            averageChangeFile << t << " " << aveChangeU << " " << aveChangeV << " " << aveChangeP << " " << u[(Ny - 1) / 2][(Nx - 1) / 2] << endl;
        }

        // Check for convergence.
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
    U_outputFile.open("02_U_output.txt", std::ios::out);
    U_outputFile << "nx     = " << Nx << "\n";
    U_outputFile << "ny     = " << Ny << "\n";
    U_outputFile << "dt     = " << timeStepSize << "\n";
    U_outputFile << "Re     = " << Re << "\n";
    U_outputFile << "n      = " << n << "\n";
    U_outputFile << "t      = " << n * timeStepSize << "\n";
    U_outputFile << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms.\n";
    U_outputFile << "y            u \n";
    for (int j = 1; j < Ny - 1; j++)
    {
        U_outputFile << std::fixed << std::setprecision(7) << lengthY - (j - 1) * h - h / 2 << "    "
                     << std::fixed << std::setprecision(7) << u[j][Nx / 2] << "\n";
    }



    //////////////////////////////////////////////////////////////////////////////
    //////////// OUTPUTTING VERTICAL PRESSURE TO OUTPUT FILE /////////////////////
    //////////// This is the average of the pressure at the left and right nodes//
    //////////// because there is no pressure node at the center vertical line of the domain //
    //////////////////////////////////////////////////////////////////////////////
    P_outputFile.open("03_P_output.txt", std::ios::out);
    P_outputFile << "# nx = " << Nx << "\n";
    P_outputFile << "# ny = " << Ny << "\n";
    P_outputFile << "# dt = " << timeStepSize << "\n";
    P_outputFile << "# Re = " << Re << "\n";
    P_outputFile << "# n = " << n << "\n";
    P_outputFile << "# t = " << n * timeStepSize << "\n";
    P_outputFile << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms.\n";
    P_outputFile << "y            p \n";
    for (int j = 1; j < Ny - 1; j++)
    {
        P_outputFile << std::fixed << std::setprecision(7) << lengthY - (j - 1) * h - h / 2 << "    "
                     << std::fixed << std::setprecision(7) << 0.5 * (p[j][(Nx / 2) - 1] + p[j][(Nx / 2) + 1]) << "\n";
    }

    U_outputFile.close();
    P_outputFile.close();
    averageChangeFile.close();

    return 0;
}

void ghostNodeValues(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny)
{

    if (b == 0)
    {
        for (int j = Ny - 2; j >= 1; j--)
        {
            M[j][0] = M[j][1];           // Left wall, dp/dx = 0
            M[j][Nx - 1] = M[j][Nx - 2]; // Right wall, dp/dx = 0
        }
        for (int i = Nx - 2; i >= 1; i--)
        {
            M[0][i] = M[1][i];           // Bottom wall, dp/dy = 0
            M[Ny - 1][i] = M[Ny - 2][i]; // Top wall, dp/dy = 0
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                               // Top left corner
        M[0][Nx - 1] = 0.5 * (M[0][Nx - 2] + M[1][Nx - 1]);                // Top right corner
        M[Ny - 1][0] = 0.5 * (M[Ny - 2][0] + M[Ny - 1][1]);                // Bottom left corner
        M[Ny - 1][Nx - 1] = 0.5 * (M[Ny - 1][Nx - 2] + M[Ny - 2][Nx - 1]); // Bottom right corner
    }
    else if (b == 1)
    {
        for (int j = Ny - 2; j >= 1; j--)
        {
            M[j][1] = uLeftWall; // Left wall, u = 0
            M[j][0] = 2 * M[j][1] - M[j][2];
            M[j][Nx - 1] = uRightWall; // Right wall, u = 0
            M[j][Nx] = 2 * M[j][Nx - 1] - M[j][Nx - 2];
        }
        for (int i = Nx - 1; i >= 1; i--)
        {
            M[0][i] = 2 * uTopWall - M[1][i];              // Top wall, ghost cell
            M[Ny - 1][i] = 2 * uBottomWall - M[Ny - 2][i]; // Bottom wall, ghost cell
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                               // Top left corner
        M[0][Nx] = 0.5 * (M[0][Nx - 1] + M[1][Nx]);                        // Top right corner
        M[Ny - 1][0] = 0.5 * (M[Ny - 2][0] + M[Ny - 1][1]);                // Bottom left corner
        M[Ny - 1][Nx - 1] = 0.5 * (M[Ny - 1][Nx - 2] + M[Ny - 2][Nx - 1]); // Bottom right corner
    }
    else if (b == 2)
    {

        for (int i = Nx - 2; i >= 1; i--)
        {
            M[1][i] = vTopWall;
            M[0][i] = 2 * M[1][i] - M[2][i]; // Top wall, v = 0

            M[Ny - 1][i] = vBottomWall; // Bottom wall, ghost cell, v = 0
            M[Ny][i] = 2 * M[Ny - 1][i] - M[Ny - 2][i];
        }
        for (int j = Ny - 1; j >= 1; j--)
        {
            M[j][0] = 2 * vLeftWall - M[j][1];            // Left wall, ghost cell
            M[j][Nx - 1] = 2 * vRightWall - M[j][Nx - 2]; // Right wall, ghost cell
        }
        M[0][0] = 0.5 * (M[0][1] + M[1][0]);                // Top left corner
        M[0][Nx - 1] = 0.5 * (M[0][Nx - 2] + M[1][Nx - 1]); // Top right corner
        M[Ny][0] = 0.5 * (M[Ny - 1][0] + M[Ny][1]);         // Bottom left corner
        M[Ny][Nx] = 0.5 * (M[Ny][Nx - 1] + M[Ny - 1][Nx]);  // Bottom right corner
    }
}
