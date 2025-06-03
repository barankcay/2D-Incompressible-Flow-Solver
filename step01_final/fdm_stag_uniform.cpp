#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting
#include <cmath>
#include <thread>
#include <chrono>
using namespace std;

//////////// THIS IS A 2D STAGGERED GRID FINITE DIFFERENCE METHOD FOR SIMULATING LID DRIVEN CAVITY FLOW ////////////
//////////// X DIRECTION VELOCITIES, U's, ARE PLACED TO THE LEFT OF THE PRESSURE NODES ////////////
//////////// Y DIRECTION VELOCITIES, V's, ARE PLACED TO THE BOTTOM OF THE PRESSURE NODES ////////////
//////////// GRID IS CONSTRUCTED USING UNIFORM GRID SPACING IN BOTH X AND Y DIRECTIONS ////////////

//////////// TIME DISCRETIZATION IS DONE USING EULER EXPLICIT METHOD AND IS FIRST ORDER ACCURATE ////////////
//////////// SPATIAL DISCRETIZATION IS DONE USING CENTRAL DIFFERENCE SCHEME AND IS SECOND ORDER ACCURATE ////////////

//////////// SOLUTION METHOD OF NAVIER STOKES EQUATIONS IS DONE USING PREDICTOR-CORRECTOR METHOD ////////////
//////////// GAUSS SEIDEL METHOD IS USED TO SOLVE THE POISSON EQUATION FOR PRESSURE ////////////

//////////// BOUNDARY CONDITIONS ARE IMPLEMENTED USING GHOST CELLS ////////////
//////////// WEST, EAST AND SOUTH WALLS ARE NO SLIP WALLS ////////////
//////////// NORTH WALL IS A MOVING WALL WITH VELOCITY U = 1, V = 0 ////////////
//////////// ALL WALLS HAVE ZERO PRESSURE GRADIENT ////////////

void setBoundaryConditions(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny);

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
    double Re                   = 1000;
    double density              = 1.0;
    double kinematicViscosity   = 1.0 / Re;
    double dynamicViscosity     = kinematicViscosity * density;

    //--------------END-OF-THE-SECTION------------------------//


    //////////////////////////////////////////////////
    ////////// BOUNDARY CONDITIONS //////////////////
    //////////////////////////////////////////////////
    double uTopWall     = 1.0; // Velocity at the top wall
    double uBottomWall  = 0.0;
    double uLeftWall    = 0.0;
    double uRightWall   = 0.0;
    double vTopWall     = 0.0;
    double vBottomWall  = 0.0;
    double vLeftWall    = 0.0;
    double vRightWall   = 0.0;
    //--------------END-OF-THE-SECTION------------------------//


    //////////////////////////////////////////////////
    ////////// GRID PARAMETERS ///////////////////////
    //////////////////////////////////////////////////
    // lengthx and lengthy are the lengths of the domain in the x and y directions, not including ghost nodes
    double lengthX      = 1;
    double lengthY      = 1;
    // Nx and Ny are the number of cells in the x and y directions, including ghost cells
    int Nx              = 20;
    double h            = lengthX / (Nx - 2);
    int Ny              = (lengthY / h) + 2; // +2 for ghost cells
    //!!!!!! Since uniform grid spacing is used, this way of calculating Ny is still valid.
    //!!!!!! If the domain is square, Nx and Ny will be equal
    //!!!!!! For non square domains, this approach still gives uniform grid spacing in both x and y directions.
    //--------------END-OF-THE-SECTION------------------------//
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////// AVERAGE CHANGE LIMITS OF PARAMETERS OF THE SIMULATION AND PRESSURE POISSON EQ  /////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    double GSchangeAve     = 1e-3;
    double pChangeAve      = 1e-9;
    double uChangeAve      = 1e-9;
    double vChangeAve      = 1e-9;
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
    double timeStepSize = min((h*h/4*dynamicViscosity), 2 * dynamicViscosity / (uTopWall*uTopWall + vTopWall*vTopWall));
    //--------------END-OF-THE-SECTION------------------------//



    //////////// INITIALIZATION and ASSIGNING DIMENSION OF THE FIELDS ////////////
    // Initialize the fields with the number of cells in the x and y directions
    u.resize(Ny, vector<double>(Nx+1));
    uStar.resize(Ny, vector<double>(Nx+1));
    vStar.resize(Ny+1, vector<double>(Nx));
    v.resize(Ny+1, vector<double>(Nx));
    p.resize(Ny, vector<double>(Nx));

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            u[j][i] = 0;
            v[j][i] = 0;
            p[j][i] = 0;
            uStar[j][i] = 0;
            vStar[j][i] = 0;
        }
    }
    //////////// END OF INITIALIZATION ////////////

    //////////// BOUNDARY CONDITIONS////////////
    setBoundaryConditions(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    setBoundaryConditions(2, v, vLeftWall, vRightWall, uTopWall, uBottomWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    setBoundaryConditions(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
    //////////// END OF BOUNDARY CONDITIONS////////////

    int n = 0; // counter for the time steps
    for (double t = startTime; t <= endTime; t = t + timeStepSize)
    {
        //////////// PREVIOUS FIELDS ARE STORED FOR CONVERGENCE CHECK ///////////
        vector<vector<double>> pPrev = p;
        vector<vector<double>> uPrev = u;
        vector<vector<double>> vPrev = v;

        //////// PREDICTOR STEP //////
        for (int i = 1; i < Nx ; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                double uAdvectX = u[j][i] * (u[j][i + 1] - u[j][i - 1]) / (2 * h);                                                        // u * du/dx (advection in X for u)
                double uAdvectY = 0.25 * (v[j][i - 1] + v[j][i] + v[j + 1][i - 1] + v[j + 1][i]) * (u[j - 1][i] - u[j + 1][i]) / (2 * h); // v * du/dy (advection in Y for u)
                double uDiffuseX = (u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / (h * h);                                                   // d²u/dx² (diffusion in X for u)
                double uDiffuseY = (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / (h * h);                                                   // d²u/dy² (diffusion in Y for u)
                uStar[j][i] = u[j][i] + timeStepSize * (kinematicViscosity * (uDiffuseX + uDiffuseY) - (uAdvectX + uAdvectY));
            }
        }

        for (int i = 1; i < Nx-1; i++)
        {
            for (int j = 1; j < Ny; j++)
            {
                double vAdvectX = 0.25 * (u[j -1][i] + u[j][i] + u[j -1][i + 1] + u[j][i + 1]) * (v[j][i + 1] - v[j][i - 1]) / (2 * h); // u * dv/dx (advection in X for v)
                double vAdvectY = v[j][i] * (v[j - 1][i] - v[j + 1][i]) / (2 * h);                                                        // v * dv/dy (advection in Y for v)
                double vDiffuseX = (v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / (h * h);                                                   // d²v/dx² (diffusion in X for v)
                double vDiffuseY = (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / (h * h);                                                   // d²v/dy² (diffusion in Y for v)
                vStar[j][i] = v[j][i] + timeStepSize * (kinematicViscosity * (vDiffuseX + vDiffuseY) - (vAdvectX + vAdvectY));
            }
            
        }
        
        setBoundaryConditions(1, uStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        setBoundaryConditions(2, vStar, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);

        ///////// END OF PREDICTOR STEP /////////
        ////////// POISSON EQUATION SOLVER //////////
        double residual = 1.0;
        int iteration = 0;
        while (residual > GSchangeAve)
        {
            vector<vector<double>> pOld = p; // Store the old pressure values for convergence check
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    double laplacianP = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) / (h*h); // Laplacian of pressure
                    double uStarGradX = (uStar[j][i + 1] - uStar[j][i]) / (h); // Gradient of uStar in X direction
                    double vStarGradY = (vStar[j][i] - vStar[j+1][i]) / (h); // Gradient of vStar in Y direction
                    p[j][i] =((h*h)/4)*(laplacianP - (density / timeStepSize) * (uStarGradX + vStarGradY)); // Update pressure using the Poisson equation
                }
            }
            residual = 0.0; // Reset residual for each iteration
            for (int i = 1; i < Nx-1; i++)
            {
                for (int j = 1; j < Ny-1; j++)
                {
                    // Calculate the residual
                    residual =residual+ abs(p[j][i] - pOld[j][i]);
                }
                
            }
            residual =residual/(Nx * Ny);
            iteration++;
            setBoundaryConditions(0, p, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);

        }
        ////// END OF POISSON EQUATION SOLVER ///////

        cout << "Number of iterations for Poisson equation: " << iteration << endl;

        //////////// CORRECTOR STEP ////////////
        for (int i = 1; i < Nx ; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                u[j][i] = uStar[j][i] - (timeStepSize / density) * ((p[j][i] - p[j][i - 1]) / (h));
                
            }
        }
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny ; j++)
            {
                v[j][i] = vStar[j][i] - (timeStepSize / density) * ((p[j][i] - p[j - 1][i]) / (h));
            }
        }
        setBoundaryConditions(1, u, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        setBoundaryConditions(2, v, uTopWall, uBottomWall, uLeftWall, uRightWall, vLeftWall, vRightWall, vTopWall, vBottomWall, Nx, Ny);
        //////////// FINAL VELOCITY AND PRESSURE VALUES ARE OBTAINED. END OF CORRECTOR STEP ////////////
        double residualU = 0.0;
        double residualV = 0.0;
        double residualP = 0.0;
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                residualU += abs(u[j][i] - uPrev[j][i]);
                residualV += abs(v[j][i] - vPrev[j][i]);
                residualP += abs(p[j][i] - pPrev[j][i]);
            }
        }
        residualU /= (Nx * Ny);
        residualV /= (Nx * Ny);
        residualP /= (Nx * Ny);

        // cout << "Time: " << t << " U velocity: " << residualU << " V velocity: " << residualV << " pressure: " << residualP << " n: " << n << endl;
        if (residualU < uChangeAve && residualV < vChangeAve && residualP < pChangeAve)
        {
            cout << "Converged at time: " << t << endl;
            break;
        }
        n++;
    }
    for (int j = 0; j < Ny; j++)
    {

        cout << (u[j][(Nx) / 2]) << endl;
    }
    cout << endl;

    cout << timeStepSize << endl;

    std::cout << "\nEnd of the main function is reached. Stopping.\n\n";
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;

    return 0;
}

void setBoundaryConditions(int b, vector<vector<double>> &M, double uTopWall, double uBottomWall, double uLeftWall, double uRightWall, double vLeftWall, double vRightWall, double vTopWall, double vBottomWall, int Nx, int Ny)
{

    if (b == 0)
    {
        for (int j = Ny - 1; j >= 0; j--)
        {
            M[j][0] = M[j][1];           // Left wall, dp/dx = 0
            M[j][Nx - 1] = M[j][Nx - 2]; // Right wall, dp/dx = 0
        }
        for (int i = Nx - 1; i >= 0; i--)
        {
            M[0][i] = M[1][i];           // Bottom wall, dp/dy = 0
            M[Ny - 1][i] = M[Ny - 2][i]; // Top wall, dp/dy = 0
        }
    }
    else if (b == 1)
    {
        for (int j = Ny - 1; j >= 0; j--)
        {

            M[j][1] = uLeftWall;       // Left wall, u = 0
            M[j][Nx - 1] = uRightWall; // Right wall, u = 0
        }
        for (int i = Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * uTopWall - M[1][i];              // Top wall, ghost cell
            M[Ny - 1][i] = 2 * uBottomWall - M[Ny - 2][i]; // Bottom wall, ghost cell
        }
    }
    else if (b == 2)
    {

        for (int i = Nx - 1; i >= 0; i--)
        {
            M[0][i] = vTopWall; // Top wall, v = 0

            M[Ny - 2][i] = vBottomWall; // Bottom wall, ghost cell, v = 0
        }
        for (int j = Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * vLeftWall - M[j][1];            // Left wall, ghost cell
            M[j][Nx - 1] = 2 * vRightWall - M[j][Nx - 2]; // Right wall, ghost cell
        }
    }

    M[0][0] = 0.5 * (M[0][1] + M[1][0]);                               // Top left corner
    M[0][Nx - 1] = 0.5 * (M[0][Nx - 2] + M[1][Nx - 1]);                // Top right corner
    M[Ny - 1][0] = 0.5 * (M[Ny - 2][0] + M[Ny - 1][1]);                // Bottom left corner
    M[Ny - 1][Nx - 1] = 0.5 * (M[Ny - 1][Nx - 2] + M[Ny - 2][Nx - 1]); // Bottom right corner
}
