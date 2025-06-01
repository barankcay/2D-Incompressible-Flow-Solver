#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting
#include <cmath>
#include <thread>
#include <chrono>
using namespace std;

struct solverConfig
{
    int Nx;                    // Number of cells in the x-direction, including ghost cells
    int Ny;                    // Number of cells in the y-direction, including ghost cells
    double lengthX;            // Length of the domain in the x-direction, not including ghost cells
    double lengthY;            // Length of the domain in the y-direction, not including ghost cells
    double hx;                 // Grid spacing in the x-direction
    double hy;                 // Grid spacing in the y-direction
    double startTime;          // Start time of the simulation
    double endTime;            // End time of the simulation. Since we have a convergence criteria, this is set to a large number
    double timeStepSize;       // Time step size. Calculated based on the Courant number and peclet number
    double kinematicViscosity; // Kinematic viscosity of the fluid. Used to calculate the Reynolds number
    double dynamicViscosity;   // Dynamic viscosity of the fluid. Dynamic viscosity = density * kinematic viscosity
    double density;            // Density of the fluid.

    double Re;

    double uTopWall;    // u Velocity at the top wall
    double uBottomWall; // u Velocity at the bottom wall
    double uLeftWall;   // v Velocity at the left wall
    double uRightWall;  // v Velocity at the right wall

    double vLeftWall;   // v Velocity at the left wall
    double vRightWall;  // v Velocity at the right wall
    double vTopWall;    // u Velocity at the top wall
    double vBottomWall; // u Velocity at the bottom wall

    double pressureOutlet; // Pressure at the outlet
    double pressureInlet;  // Pressure at the inlet

    double courantNumber; // Courant number used to calculate the time step size

    double poissonTolerance; // Tolerance for the Poisson equation solver
    double timeToleranceU;   // Tolerance for the time loop convergence
    double timeToleranceV;   // Tolerance for the time loop convergence
    double timeToleranceP;   // Tolerance for the time loop convergence

    // non-uniform mesh properties

    double firstCellHeightY; // height of the first cell for the non-uniform mesh
    double firstCellHeightX; // height of the first cell for the non-uniform mesh

    double smallestCellHeightY; // height of the first cell for the non-uniform mesh
    double smallestCellHeightX; // height of the first cell for the non-uniform mesh

    double maxUvelocity;
    double maxVvelocity;

    double inflationLayerThicknessN; // boundary layer thickness for the non-uniform mesh
    double inflationLayerThicknessS; // boundary layer thickness for the non-uniform mesh
    double inflationLayerThicknessE; // boundary layer thickness for the non-uniform mesh
    double inflationLayerThicknessW; // boundary layer thickness for the non-uniform mesh

    double growthRateE; // growth rate for the non-uniform mesh
    double growthRateW; // growth rate for the non-uniform mesh
    double growthRateN; // growth rate for the non-uniform mesh
    double growthRateS; // growth rate for the non-uniform mesh

    int numberOfMeshLayersE; // number of mesh layers for the non-uniform mesh
    int numberOfMeshLayersW; // number of mesh layers for the non-uniform mesh
    int numberOfMeshLayersN; // number of mesh layers for the non-uniform mesh
    int numberOfMeshLayersS; // number of mesh layers for the non-uniform mesh
};

struct fields
{
    vector<double> ux; // x coordinate of u velocity field
    vector<double> uy; // y coordinate of u velocity field

    vector<double> vx; // x coordinate of v velocity field
    vector<double> vy; // y coordinate of v velocity field

    vector<double> px; // x coordinate of pressure field
    vector<double> py; // y coordinate of pressure field

    vector<vector<double>> areaUn;
    vector<vector<double>> areaUs;
    vector<vector<double>> areaUe;
    vector<vector<double>> areaUw;

    vector<vector<double>> areaVn;
    vector<vector<double>> areaVs;
    vector<vector<double>> areaVe;
    vector<vector<double>> areaVw;

    vector<vector<double>> areaPn;
    vector<vector<double>> areaPs;
    vector<vector<double>> areaPe;
    vector<vector<double>> areaPw;

    vector<vector<double>> volumeU; // volume of u control volume
    vector<vector<double>> volumeV; // volume of v control volume

    vector<vector<double>> u;     // u velocity field
    vector<vector<double>> uStar; // intermediate u velocity field
    vector<vector<double>> uNew;  // next time step u velocity field
    vector<vector<double>> v;     // v velocity field
    vector<vector<double>> vStar; // intermediate v velocity field
    vector<vector<double>> vNew;  // next time step v velocity field
    vector<vector<double>> p;     // pressure field

    fields(int Nx, int Ny)
    {
        ux.resize(Nx - 2);
        uy.resize(Ny - 2);

        vx.resize(Nx - 2);
        vy.resize(Ny - 2);

        px.resize(Nx);
        py.resize(Ny);

        u.resize(Ny, vector<double>(Nx));
        v.resize(Ny, vector<double>(Nx));
        p.resize(Ny, vector<double>(Nx));

        uStar.resize(Ny, vector<double>(Nx));
        vStar.resize(Ny, vector<double>(Nx));
        uNew.resize(Ny, vector<double>(Nx));
        vNew.resize(Ny, vector<double>(Nx));

        areaUn.resize(Ny - 2, vector<double>(Nx - 2));
        areaUs.resize(Ny - 2, vector<double>(Nx - 2));
        areaUe.resize(Ny - 2, vector<double>(Nx - 2));
        areaUw.resize(Ny - 2, vector<double>(Nx - 2));

        areaVn.resize(Ny - 2, vector<double>(Nx - 2));
        areaVs.resize(Ny - 2, vector<double>(Nx - 2));
        areaVe.resize(Ny - 2, vector<double>(Nx - 2));
        areaVw.resize(Ny - 2, vector<double>(Nx - 2));

        areaPn.resize(Ny - 2, vector<double>(Nx - 2));
        areaPs.resize(Ny - 2, vector<double>(Nx - 2));
        areaPe.resize(Ny - 2, vector<double>(Nx - 2));
        areaPw.resize(Ny - 2, vector<double>(Nx - 2));

        volumeU.resize(Ny - 2, vector<double>(Nx - 2));
        volumeV.resize(Ny - 2, vector<double>(Nx - 2));
    }
};
void setBoundaryConditionsLDC(int b, vector<vector<double>> &M, solverConfig &cfg)
{

    if (b == 0)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = M[j][1];                   // Left wall, dp/dx = 0
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2]; // Right wall, dp/dx = 0
        }
        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = M[1][i];                   // Bottom wall, dp/dy = 0
            M[cfg.Ny - 1][i] = M[cfg.Ny - 2][i]; // Top wall, dp/dy = 0
        }
    }
    else if (b == 1)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 0;          // Left wall,ghost cell,  u = 0
            M[j][1] = 0;          // Left wall, u = 0
            M[j][cfg.Nx - 1] = 0; // Right wall, u = 0
        }
        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * cfg.uTopWall - M[1][i];                      // Top wall, ghost cell
            M[cfg.Ny - 1][i] = 2 * cfg.uBottomWall - M[cfg.Ny - 2][i]; // Bottom wall, ghost cell
        }
    }
    else if (b == 2)
    {

        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 0;          // Top wall, v = 0
            M[cfg.Ny - 1][i] = 0; // Bottom wall, v = 0
            M[cfg.Ny - 2][i] = 0; // Bottom wall, ghost cell, v = 0
        }
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * cfg.vLeftWall - M[j][1];                    // Left wall, ghost cell
            M[j][cfg.Nx - 1] = 2 * cfg.vRightWall - M[j][cfg.Nx - 2]; // Right wall, ghost cell
        }
    }
    M[0][0] = 0.5 * (M[0][1] + M[1][0]);                                                       // Top left corner
    M[0][cfg.Nx - 1] = 0.5 * (M[0][cfg.Nx - 2] + M[1][cfg.Nx - 1]);                            // Top right corner
    M[cfg.Ny - 1][0] = 0.5 * (M[cfg.Ny - 2][0] + M[cfg.Ny - 1][1]);                            // Bottom left corner
    M[cfg.Ny - 1][cfg.Nx - 1] = 0.5 * (M[cfg.Ny - 1][cfg.Nx - 2] + M[cfg.Ny - 2][cfg.Nx - 1]); // Bottom right corner
}

void setBoundaryConditionsCF_pressureDriven(int b, vector<vector<double>> &M, solverConfig &cfg)
{
    if (b == 0)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * cfg.pressureInlet - M[j][1]; // pressure driven inlet

            M[j][cfg.Nx - 1] = 2 * cfg.pressureOutlet - M[j][cfg.Nx - 2]; // Right wall, p=0
        }
        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = M[1][i];                   // Bottom wall, dp/dy = 0
            M[cfg.Ny - 1][i] = M[cfg.Ny - 2][i]; // Top wall, dp/dy = 0
        }
    }
    else if (b == 1)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][1] = M[j][2]; // Left wall, u = 0

            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2]; // Right wall, du/dx = 0
        }
        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * cfg.uTopWall - M[1][i];                      // Top wall, ghost cell
            M[cfg.Ny - 1][i] = 2 * cfg.uBottomWall - M[cfg.Ny - 2][i]; // Bottom wall, ghost cell
        }
    }
    else if (b == 2)
    {

        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = cfg.vTopWall;                                    // Top wall, v = 0
            M[cfg.Ny - 1][i] = 2 * cfg.vBottomWall - M[cfg.Ny - 3][i]; // Bottom wall, v = 0
            M[cfg.Ny - 2][i] = cfg.vBottomWall;                        // Bottom wall, ghost cell, v = 0
        }
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = M[j][1];
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2]; // Right wall, ghost cell
        }
    }

    M[0][0] = 0.5 * (M[0][1] + M[1][0]);                                                       // Top left corner
    M[0][cfg.Nx - 1] = 0.5 * (M[0][cfg.Nx - 2] + M[1][cfg.Nx - 1]);                            // Top right corner
    M[cfg.Ny - 1][0] = 0.5 * (M[cfg.Ny - 2][0] + M[cfg.Ny - 1][1]);                            // Bottom left corner
    M[cfg.Ny - 1][cfg.Nx - 1] = 0.5 * (M[cfg.Ny - 1][cfg.Nx - 2] + M[cfg.Ny - 2][cfg.Nx - 1]); // Bottom right corner
}

void setBoundaryConditionsCF_velocityDriven(int b, vector<vector<double>> &M, solverConfig &cfg)
{
    if (b == 0)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = M[j][1]; // Left wall, dp/dx = 0

            M[j][cfg.Nx - 1] = 2 * cfg.pressureOutlet - M[j][cfg.Nx - 2]; // Right wall, p=0
            // M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2];
        }
        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = M[1][i];                   // Bottom wall, dp/dy = 0
            M[cfg.Ny - 1][i] = M[cfg.Ny - 2][i]; // Top wall, dp/dy = 0
        }
    }
    else if (b == 1)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * cfg.uLeftWall - M[j][2]; // Left wall,ghost cell,  u = 0
            M[j][1] = cfg.uLeftWall;               // Left wall, u = 0

            // M[j][0] = M[j][1];
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2]; // Right wall, du/dx = 0
        }
        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * cfg.uTopWall - M[1][i];                      // Top wall, ghost cell
            M[cfg.Ny - 1][i] = 2 * cfg.uBottomWall - M[cfg.Ny - 2][i]; // Bottom wall, ghost cell
        }
    }
    else if (b == 2)
    {

        for (int i = cfg.Nx - 1; i >= 0; i--)
        {
            M[0][i] = cfg.vTopWall;                                    // Top wall, v = 0
            M[cfg.Ny - 1][i] = 2 * cfg.vBottomWall - M[cfg.Ny - 3][i]; // Bottom wall, v = 0
            M[cfg.Ny - 2][i] = cfg.vBottomWall;                        // Bottom wall, ghost cell, v = 0
        }
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            // M[j][0] = 2 * cfg.vLeftWall - M[j][1]; // Left wall, ghost cell
            M[j][0] = 2 * cfg.vLeftWall - M[j][1];
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2]; // Right wall, ghost cell
        }
    }

    M[0][0] = 0.5 * (M[0][1] + M[1][0]);                                                       // Top left corner
    M[0][cfg.Nx - 1] = 0.5 * (M[0][cfg.Nx - 2] + M[1][cfg.Nx - 1]);                            // Top right corner
    M[cfg.Ny - 1][0] = 0.5 * (M[cfg.Ny - 2][0] + M[cfg.Ny - 1][1]);                            // Bottom left corner
    M[cfg.Ny - 1][cfg.Nx - 1] = 0.5 * (M[cfg.Ny - 1][cfg.Nx - 2] + M[cfg.Ny - 2][cfg.Nx - 1]); // Bottom right corner
}
// function to initialize the fields
void initialization(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx; i++)
    {
        for (int j = 0; j < cfg.Ny; j++)
        {
            f.u[j][i] = 0;
            f.v[j][i] = 0;
            f.p[j][i] = 0;
            f.uStar[j][i] = 0;
            f.vStar[j][i] = 0;
            f.uNew[j][i] = 0;
            f.vNew[j][i] = 0;
        }
    }
}

void nonUniformMesh(solverConfig &cfg)
{
    if (cfg.numberOfMeshLayersN + cfg.numberOfMeshLayersS > cfg.Ny - 2)
    {
        cout << "Number of mesh layers is greater than the number of cells in Y direction. Please reduce the number of mesh layers." << endl;
        exit(1);
    }
    double geomSeriesN = ((1 - pow(cfg.growthRateN, cfg.numberOfMeshLayersN)) / (1 - cfg.growthRateN));

    if (cfg.numberOfMeshLayersN == 0)
    {
        geomSeriesN = 0; // height of the first cell for the non-uniform mesh
    }
    double geomSeriesS = ((1 - pow(cfg.growthRateS, cfg.numberOfMeshLayersS)) / (1 - cfg.growthRateS));

    if (cfg.numberOfMeshLayersS == 0)
    {
        geomSeriesS = 0; // height of the first cell for the non-uniform mesh
    }

    if (cfg.numberOfMeshLayersE + cfg.numberOfMeshLayersW > cfg.Nx - 2)
    {
        cout << "Number of mesh layers is greater than the number of cells in X direction. Please reduce the number of mesh layers." << endl;
        exit(1);
    }
    double geomSeriesE = ((1 - pow(cfg.growthRateE, cfg.numberOfMeshLayersE)) / (1 - cfg.growthRateE));

    if (cfg.numberOfMeshLayersE == 0)
    {
        geomSeriesE = 0; // height of the first cell for the non-uniform mesh
    }
    double geomSeriesW = ((1 - pow(cfg.growthRateW, cfg.numberOfMeshLayersW)) / (1 - cfg.growthRateW));

    if (cfg.numberOfMeshLayersW == 0)
    {
        geomSeriesW = 0; // height of the first cell for the non-uniform mesh
    }
    double Gx = geomSeriesE + geomSeriesW + (cfg.Nx - 2 - cfg.numberOfMeshLayersE - cfg.numberOfMeshLayersW) * 0.5 * (pow(cfg.growthRateE, cfg.numberOfMeshLayersE) + pow(cfg.growthRateW, cfg.numberOfMeshLayersW)); // Geometric series for the non-uniform mesh
    double Gy = geomSeriesS + geomSeriesN + (cfg.Ny - 2 - cfg.numberOfMeshLayersN - cfg.numberOfMeshLayersS) * 0.5 * (pow(cfg.growthRateN, cfg.numberOfMeshLayersN) + pow(cfg.growthRateS, cfg.numberOfMeshLayersS)); // Geometric series for the non-uniform mesh
    cfg.firstCellHeightY = cfg.lengthY / Gy;
    cfg.firstCellHeightX = cfg.lengthX / Gx; // height of the first cell for the non-uniform mesh

    cfg.inflationLayerThicknessN = 0;
    for (int i = 0; i < cfg.numberOfMeshLayersN; i++)
    {
        cfg.inflationLayerThicknessN = cfg.firstCellHeightY + cfg.inflationLayerThicknessN * cfg.growthRateN; // boundary layer thickness
    }
    cfg.inflationLayerThicknessS = 0;
    for (int i = 0; i < cfg.numberOfMeshLayersS; i++)
    {
        cfg.inflationLayerThicknessS = cfg.firstCellHeightY + cfg.inflationLayerThicknessS * cfg.growthRateS; // boundary layer thickness
    }

    if (cfg.inflationLayerThicknessN + cfg.inflationLayerThicknessS > cfg.lengthY + 1e-10)
    {
        cout << "Inflation layer thickness in y is " << cfg.inflationLayerThicknessN + cfg.inflationLayerThicknessS << " and greater than the length of the domain. Please reduce the number of mesh layers or the growth rate." << endl;
        exit(1);
    }

    cfg.inflationLayerThicknessE = 0;
    for (int i = 0; i < cfg.numberOfMeshLayersE; i++)
    {
        cfg.inflationLayerThicknessE = cfg.firstCellHeightX + cfg.inflationLayerThicknessE * cfg.growthRateE; // boundary layer thickness
    }

    cfg.inflationLayerThicknessW = 0;
    for (int i = 0; i < cfg.numberOfMeshLayersW; i++)
    {
        cfg.inflationLayerThicknessW = cfg.firstCellHeightX + cfg.inflationLayerThicknessW * cfg.growthRateW; // boundary layer thickness
    }
    if (cfg.inflationLayerThicknessE + cfg.inflationLayerThicknessW > cfg.lengthX + 1e-10)
    {
        cout << "Inflation layer thickness in x is " << cfg.inflationLayerThicknessE + cfg.inflationLayerThicknessW << " and greater than the length of the domain. Please reduce the number of mesh layers or the growth rate." << endl;
        exit(1);
    }
    cfg.hy = (cfg.lengthY - cfg.inflationLayerThicknessN - cfg.inflationLayerThicknessS) / (cfg.Ny - 2 - cfg.numberOfMeshLayersN - cfg.numberOfMeshLayersS); // Grid spacing in the y-direction
    cfg.hx = (cfg.lengthX - cfg.inflationLayerThicknessE - cfg.inflationLayerThicknessW) / (cfg.Nx - 2 - cfg.numberOfMeshLayersE - cfg.numberOfMeshLayersW); // Grid spacing in the x-direction
}

// function to create the coordinates of the nodes of cells
void createCoordinates(fields &f, solverConfig &cfg)
{
    f.vy.insert(f.vy.begin(), cfg.lengthY); // y coordinate of v velocity field
    f.ux.push_back(cfg.lengthX);            // x coordinate of u velocity field
    f.ux[0] = 0;
    f.vy.back() = 0; // y coordinate of v velocity field

    for (int j = 0; j < cfg.numberOfMeshLayersN; j++)
    {
        f.vy[j + 1] = cfg.lengthY - (cfg.firstCellHeightY + (cfg.lengthY - f.vy[j]) * cfg.growthRateN);
        f.uy[j] = 0.5 * (f.vy[j + 1] + f.vy[j]);
        // y coordinate of u velocity field
    }
    for (int i = 0; i < cfg.numberOfMeshLayersS; i++)
    {
        f.vy[cfg.Ny - 3 - i] = cfg.firstCellHeightY + (f.vy[cfg.Ny - 2 - i]) * cfg.growthRateS;
        f.uy[cfg.Ny - 3 - i] = 0.5 * (f.vy[cfg.Ny - 3 - i] + f.vy[cfg.Ny - 2 - i]);
    }
    for (int i = 0; i < cfg.numberOfMeshLayersW; i++)
    {
        f.ux[i + 1] = cfg.firstCellHeightX + (f.ux[i]) * cfg.growthRateW; // x coordinate of u velocity field
        f.vx[i] = 0.5 * (f.ux[i + 1] + f.ux[i]);                          // x coordinate of v velocity field
    }

    for (int i = 0; i < cfg.numberOfMeshLayersE; i++)
    {
        f.ux[cfg.Nx - 3 - i] = cfg.lengthX - (cfg.firstCellHeightX + (cfg.lengthX - f.ux[cfg.Nx - 2 - i]) * cfg.growthRateE); // x coordinate of u velocity field
        f.vx[cfg.Nx - 3 - i] = 0.5 * (f.ux[cfg.Nx - 3 - i] + f.ux[cfg.Nx - 2 - i]);                                           // y coordinate of u velocity field
    }

    cout << "Grid spacing in the y-direction: " << cfg.hy << endl;
    for (int j = cfg.numberOfMeshLayersN; j < cfg.Ny - 2 - cfg.numberOfMeshLayersS; j++)
    {
        f.vy[j + 1] = f.vy[j] - cfg.hy;          // y coordinate of v velocity field
        f.uy[j] = 0.5 * (f.vy[j + 1] + f.vy[j]); // y coordinate of u velocity field
    }
    for (int i = cfg.numberOfMeshLayersW; i < cfg.Nx - 2 - cfg.numberOfMeshLayersE; i++)
    {
        f.ux[i + 1] = f.ux[i] + cfg.hx;          // x coordinate of u velocity field
        f.vx[i] = 0.5 * (f.ux[i + 1] + f.ux[i]); // x coordinate of v velocity field
    }

    f.ux.insert(f.ux.begin(), -f.ux[1]); // x coordinate of u velocity field
    f.uy.push_back(-f.uy[cfg.Ny - 3]);   // y coordinate of u velocity field
    f.uy.insert(f.uy.begin(), f.vy[0] + (f.vy[0] - f.uy[0]));

    f.vx.push_back(f.ux[cfg.Nx - 1] + (f.ux[cfg.Nx - 1] - f.vx[cfg.Nx - 3])); // x coordinate of v velocity field
    f.vx.insert(f.vx.begin(), -f.vx[0]);

    f.vy.push_back(-f.vy[cfg.Ny - 3]);

    f.px = f.vx; // x coordinate of pressure field
    f.py = f.uy; // y coordinate of pressure field
}

void smallestCellHeight(fields &f, solverConfig &cfg)
{
    cfg.smallestCellHeightY = f.vy[0] - f.vy[1];
    cfg.smallestCellHeightX = f.ux[2] - f.ux[1];

    for (int i = 2; i < cfg.Nx - 1; i++)
    {
        if (f.ux[i + 1] - f.ux[i] < cfg.smallestCellHeightX)
        {
            cfg.smallestCellHeightX = f.ux[i + 1] - f.ux[i];
        }
    }
    for (int j = 1; j < cfg.Ny - 2; j++)
    {
        if (f.vy[j] - f.vy[j + 1] < cfg.smallestCellHeightY)
        {
            cfg.smallestCellHeightY = f.vy[j] - f.vy[j + 1];
        }
    }
}

void maximumVelocity(fields &f, solverConfig &cfg)
{
    cfg.maxUvelocity = f.u[0][0];
    cfg.maxVvelocity = f.v[0][0];

    for (int i = 0; i < cfg.Nx; i++)
    {
        for (int j = 0; j < cfg.Ny; j++)
        {
            if (f.u[j][i] > cfg.maxUvelocity)
            {
                cfg.maxUvelocity = f.u[j][i];
            }
            if (f.v[j][i] > cfg.maxVvelocity)
            {
                cfg.maxVvelocity = f.v[j][i];
            }
        }
    }
}
void createAreaU(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.areaUn[j][i] = 1 * (f.vx[i + 1] - f.vx[i]); // area of the cell in the x-direction
            f.areaUs[j][i] = 1 * (f.vx[i + 1] - f.vx[i]); // area of the cell in the y-direction
            f.areaUe[j][i] = 1 * (f.vy[j] - f.vy[j + 1]); // area of the cell in the x-direction
            f.areaUw[j][i] = 1 * (f.vy[j] - f.vy[j + 1]);
        }
    }
}

void createAreaV(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.areaVn[j][i] = 1 * (f.ux[i + 2] - f.ux[i + 1]); // area of the cell in the x-direction
            f.areaVs[j][i] = 1 * (f.ux[i + 2] - f.ux[i + 1]); // area of the cell in the y-direction
            f.areaVe[j][i] = 1 * (f.uy[j + 1] - f.uy[j + 2]); // area of the cell in the x-direction
            f.areaVw[j][i] = 1 * (f.uy[j + 1] - f.uy[j + 2]);
        }
    }
}
void createAreaP(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.areaPn[j][i] = 1 * (f.ux[i + 2] - f.ux[i + 1]); // area of the cell in the x-direction
            f.areaPs[j][i] = 1 * (f.ux[i + 2] - f.ux[i + 1]); // area of the cell in the y-direction
            f.areaPe[j][i] = 1 * (f.vy[j] - f.vy[j + 1]);     // area of the cell in the x-direction
            f.areaPw[j][i] = 1 * (f.vy[j] - f.vy[j + 1]);     // area of the cell in the y-direction
        }
    }
}

void createVolumeU(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.volumeU[j][i] = 1 * f.areaUn[j][i] * f.areaUe[j][i]; // volume of the cell in the x-direction
        }
    }
}

void createVolumeV(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.volumeV[j][i] = 1 * f.areaVn[j][i] * f.areaVe[j][i]; // volume of the cell in the x-direction
        }
    }
}

void predictor(fields &f, solverConfig &cfg)
{
    // Predictor step for the u velocity field
    for (int j = 1; j < cfg.Ny - 1; j++)
    {
        for (int i = 1; i < cfg.Nx - 1; i++)
        {
            double vn = f.v[j - 1][i - 1] + (f.v[j - 1][i] - f.v[j - 1][i - 1]) * (f.ux[i] - f.vx[i - 1]) / (f.vx[i] - f.vx[i - 1]); // v velocity field
            double un = f.u[j][i] + (f.u[j - 1][i] - f.u[j][i]) * (f.vy[j - 1] - f.uy[j]) / (f.uy[j - 1] - f.uy[j]);
            double ue = f.u[j][i] + (f.u[j][i + 1] - f.u[j][i]) * (f.vx[i] - f.ux[i]) / (f.ux[i + 1] - f.ux[i]);             // u velocity field
            double uw = f.u[j][i - 1] + (f.u[j][i] - f.u[j][i - 1]) * (f.vx[i - 1] - f.ux[i - 1]) / (f.ux[i] - f.ux[i - 1]); // u velocity field
            double vs = f.v[j][i - 1] + (f.v[j][i] - f.v[j][i - 1]) * (f.ux[i] - f.vx[i - 1]) / (f.vx[i] - f.vx[i - 1]);
            double us = f.u[j + 1][i] + (f.u[j][i] - f.u[j + 1][i]) * (f.vy[j] - f.uy[j + 1]) / (f.uy[j] - f.uy[j + 1]);
            double An = f.areaUn[j - 1][i - 1];
            double As = f.areaUs[j - 1][i - 1];
            double Ae = f.areaUe[j - 1][i - 1];
            double Aw = f.areaUw[j - 1][i - 1];
            double advection = vn * un * An - vs * us * As + ue * ue * Ae - uw * uw * Aw;
            double diffusion = cfg.kinematicViscosity * ((Ae * (f.u[j][i + 1] - f.u[j][i]) / (f.ux[i + 1] - f.ux[i])) - (Aw * (f.u[j][i] - f.u[j][i - 1]) / (f.ux[i] - f.ux[i - 1])) + (An * (f.u[j - 1][i] - f.u[j][i]) / (f.uy[j - 1] - f.uy[j])) - (As * (f.u[j][i] - f.u[j + 1][i]) / (f.uy[j] - f.uy[j + 1]))); // diffusion term
            f.uStar[j][i] = f.u[j][i] + (cfg.timeStepSize / f.volumeU[j - 1][i - 1]) * (-advection + diffusion);                                                                                                                                                                                                     // u velocity field
        }
    }
    // Predictor step for the v velocity field
    for (int j = 1; j < cfg.Ny - 1; j++)
    {
        for (int i = 1; i < cfg.Nx - 1; i++)
        {
            double vn = f.v[j][i] + (f.v[j - 1][i] - f.v[j][i]) * ((f.uy[j] - f.vy[j]) / (f.vy[j - 1] - f.vy[j]));
            double uw = f.u[j + 1][i] + (f.u[j][i] - f.u[j + 1][i]) * ((f.vy[j] - f.uy[j + 1]) / (f.uy[j] - f.uy[j + 1]));
            double vw = f.v[j][i - 1] + (f.v[j][i] - f.v[j][i - 1]) * ((f.ux[i] - f.vx[i - 1]) / (f.vx[i] - f.vx[i - 1]));
            double vs = f.v[j + 1][i] + (f.v[j][i] - f.v[j + 1][i]) * ((f.uy[j + 1] - f.vy[j + 1]) / (f.vy[j] - f.vy[j + 1]));
            double ue = f.u[j + 1][i + 1] + (f.u[j][i + 1] - f.u[j + 1][i + 1]) * ((f.vy[j] - f.uy[j + 1]) / (f.uy[j] - f.uy[j + 1]));
            double ve = f.v[j][i] + (f.v[j][i + 1] - f.v[j][i]) * ((f.ux[i + 1] - f.vx[i]) / (f.vx[i + 1] - f.vx[i]));
            double An = f.areaVn[j - 1][i - 1];
            double As = f.areaVs[j - 1][i - 1];
            double Ae = f.areaVe[j - 1][i - 1];
            double Aw = f.areaVw[j - 1][i - 1];
            double advection = vn * vn * An - vs * vs * As + ue * ve * Ae - uw * vw * Aw;
            double diffusion = cfg.kinematicViscosity * ((Ae * (f.v[j][i + 1] - f.v[j][i]) / (f.vx[i + 1] - f.vx[i])) - (Aw * (f.v[j][i] - f.v[j][i - 1]) / (f.vx[i] - f.vx[i - 1])) + (An * (f.v[j - 1][i] - f.v[j][i]) / (f.vy[j - 1] - f.vy[j])) - (As * (f.v[j][i] - f.v[j + 1][i]) / (f.vy[j] - f.vy[j + 1])));
            f.vStar[j][i] = f.v[j][i] + (cfg.timeStepSize / f.volumeV[j - 1][i - 1]) * (-advection + diffusion);
        }
    }
    setBoundaryConditionsLDC(1, f.uStar, cfg);
    setBoundaryConditionsLDC(2, f.vStar, cfg);
}

void pressurePoisson(fields &f, solverConfig &cfg)
{
    double residual = 1.0;
    int iteration = 0;
    while (residual > cfg.poissonTolerance)
    {
        residual = 0;
        for (int i = 1; i < cfg.Nx - 1; i++)
        {
            for (int j = 1; j < cfg.Ny - 1; j++)
            {
                double pOld = f.p[j][i];

                double delX1 = f.px[i - 1] - f.px[i];
                double delX2 = f.px[i + 1] - f.px[i];

                double delY1 = f.py[j + 1] - f.py[j];
                double delY2 = f.py[j - 1] - f.py[j];
                double velocStarGrad = (1 / cfg.timeStepSize) * (((f.uStar[j][i + 1] - f.uStar[j][i]) / (f.ux[i + 1] - f.ux[i])) + ((f.vStar[j - 1][i] - f.vStar[j][i]) / (f.vy[j - 1] - f.vy[j])));

                double A = delX1 * delX2 * delX2 - delX2 * delX1 * delX1;

                double B = delY1 * delY2 * delY2 - delY2 * delY1 * delY1;

                double pCoeff = (2 * (delX1 - delX2) / A) + (2 * (delY1 - delY2) / B);
                f.p[cfg.Ny - 2][1] = 0; //
                f.p[j][i] = ((2 * f.p[j][i + 1] * delX1 / A) - (2 * f.p[j][i - 1] * delX2 / A) + (2 * f.p[j - 1][i] * delY1 / B) - (2 * f.p[j + 1][i] * delY2 / B) - velocStarGrad) / pCoeff;
                f.p[cfg.Ny - 2][1] = 0; //
                residual += abs(f.p[j][i] - pOld);
            }
        }
        setBoundaryConditionsLDC(0, f.p, cfg);
        residual /= ((cfg.Nx) * (cfg.Ny));
        iteration++;
    }
    cout << "Number of iterations: " << iteration << endl;
}

void updateTimeStepSize(fields &f, solverConfig &cfg)
{

    // Calculate the maximum velocity in the domain (including boundaries)
    for (int i = 0; i < cfg.Nx; i++)
    {
        for (int j = 0; j < cfg.Ny; j++)
        {
            double rule1 = (pow(cfg.smallestCellHeightX, 2) * pow(cfg.smallestCellHeightY, 2)) / (2 * cfg.dynamicViscosity * (pow(cfg.smallestCellHeightX, 2) + pow(cfg.smallestCellHeightY, 2)));
            double rule2 = 2 * cfg.dynamicViscosity / (pow(f.u[j][i], 2) + pow(f.v[j][i], 2));
            double finalRule = min(rule1, rule2);
            if (finalRule < cfg.timeStepSize)
            {
                cfg.timeStepSize = finalRule;
            }
        }
    }
}
void corrector(fields &f, solverConfig &cfg)
{
    for (int i = 1; i < cfg.Nx - 1; i++)
    {
        for (int j = 1; j < cfg.Ny - 1; j++)
        {
            f.uNew[j][i] = f.uStar[j][i] + (cfg.timeStepSize / f.volumeU[j - 1][i - 1]) * ((f.p[j][i - 1] * f.areaUw[j - 1][i - 1] - f.p[j][i] * f.areaUe[j - 1][i - 1]) / (cfg.density));
        }
    }

    for (int i = 1; i < cfg.Nx - 1; i++)
    {
        for (int j = 1; j < cfg.Ny - 1; j++)
        {
            f.vNew[j][i] = f.vStar[j][i] + (cfg.timeStepSize / f.volumeV[j - 1][i - 1]) * ((f.p[j + 1][i] * f.areaVs[j - 1][i - 1] - f.p[j][i] * f.areaVn[j - 1][i - 1]) / (cfg.density));
        }
    }
}

double checkConvergence(vector<vector<double>> &Mnew, vector<vector<double>> &Mold, solverConfig &cfg)
{
    double residual = 0.0;
    for (int i = 1; i < cfg.Nx - 1; i++)
    {
        for (int j = 1; j < cfg.Ny - 1; j++)
        {
            residual += abs(Mnew[j][i] - Mold[j][i]);
        }
    }
    return residual /= ((cfg.Nx) * (cfg.Ny));
}

void swapFields(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx; i++)
    {
        for (int j = 0; j < cfg.Ny; j++)
        {
            f.u[j][i] = f.uNew[j][i];
            f.v[j][i] = f.vNew[j][i];
            f.uStar[j][i] = f.u[j][i];
            f.vStar[j][i] = f.v[j][i];
        }
    }
}
#include <iomanip> // For std::setprecision

void point_writeVTKFile(const fields &f, const solverConfig &cfg)
{
    std::string filename = "final_results.vtk";
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open())
    {
        std::cerr << "Error opening VTK file for writing!" << std::endl;
        return;
    }

    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Final Time Step Results\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << cfg.Nx << " " << cfg.Ny << " 1\n";

    // Write pressure cell center coordinates
    vtkFile << "POINTS " << cfg.Nx * cfg.Ny << " float\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.px[i] << " " << f.py[j] << " 0.0\n";
        }
    }

    // Write cell data (pressure)
    vtkFile << "\nCELL_DATA " << (cfg.Nx - 1) * (cfg.Ny - 1) << "\n";
    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < cfg.Ny - 1; j++)
    {
        for (int i = 0; i < cfg.Nx - 1; i++)
        {
            vtkFile << f.p[j][i] << "\n";
        }
    }

    // Write point data (velocities at pressure points)
    vtkFile << "\nPOINT_DATA " << cfg.Nx * cfg.Ny << "\n";

    // U velocity component
    vtkFile << "SCALARS u_velocity float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.u[j][i] << "\n";
        }
    }

    // V velocity component
    vtkFile << "SCALARS v_velocity float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.v[j][i] << "\n";
        }
    }

    // Velocity vector
    vtkFile << "VECTORS velocity_vector float\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.u[j][i] << " " << f.v[j][i] << " 0.0\n";
        }
    }

    vtkFile.close();
    std::cout << "Final results written to: " << filename << std::endl;
}

// Add this at the end of your main() function, just before the return statement:

void cell_writeVTKFile(fields &f, solverConfig &cfg)
{
    f.vy.insert(f.vy.begin(), f.vy[0] + 2 * (f.uy[0] - f.vy[0])); // Insert the first element at the beginning of f.vy
    f.ux.push_back(f.ux[cfg.Nx - 1] + 2 * (f.vx[cfg.Nx - 1] - f.ux[cfg.Nx - 1]));
    std::string filename = "staggered_results.vtk";
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open())
    {
        std::cerr << "Error opening VTK file!" << std::endl;
        return;
    }

    // Force consistent decimal formatting (avoid locale issues)
    vtkFile << std::fixed << std::setprecision(6);

    // VTK Header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "15x15 Cell-Centered Data\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET RECTILINEAR_GRID\n";

    // Grid dimensions (NODES, not cells)
    vtkFile << "DIMENSIONS " << (cfg.Nx + 1) << " " << (cfg.Ny + 1) << " 1\n";

    // Node coordinates (must have Nx+1 and Ny+1 entries)
    vtkFile << "X_COORDINATES " << (cfg.Nx + 1) << " float\n";
    for (int i = 0; i <= cfg.Nx; i++)
    {
        vtkFile << f.ux[i] << "\n"; // Ensure f.px has Nx+1 entries
    }

    vtkFile << "Y_COORDINATES " << (cfg.Ny + 1) << " float\n";
    for (int j = 0; j <= cfg.Ny; j++)
    {
        vtkFile << f.vy[j] << "\n"; // Ensure f.py has Ny+1 entries
    }

    vtkFile << "Z_COORDINATES 1 float\n0.0\n";

    // CELL_DATA for 15x15 cells
    vtkFile << "CELL_DATA " << (cfg.Nx * cfg.Ny) << "\n";

    // Pressure (15x15)
    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.p[j][i] << "\n";
        }
    }

    // Velocity components (15x15)
    vtkFile << "VECTORS velocity float\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.u[j][i] << " " << f.v[j][i] << " 0.0\n";
        }
    }

    vtkFile.close();
    std::cout << "VTK file successfully written: " << filename << std::endl;
}
void vertex_writeVTKFile(fields &f, solverConfig &cfg)
{

    double d1;
    double d2;
    double d3;
    double d4;

    double w1;
    double w2;
    double w3;
    double w4;

    std::string filename = "stag_vertex_results.vtk";
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open())
    {
        std::cerr << "Error opening VTK file!" << std::endl;
        return;
    }

    // Force consistent decimal formatting
    vtkFile << std::fixed << std::setprecision(6);

    // VTK Header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Cell-Centered Data as Point Data\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";

    // Grid dimensions (points at cell centers)
    vtkFile << "DIMENSIONS " << cfg.Nx - 1 << " " << cfg.Ny - 1 << " 1\n";

    // Points at cell centers
    vtkFile << "POINTS " << ((cfg.Nx - 1) * (cfg.Ny - 1)) << " float\n";
    for (int j = 1; j < cfg.Ny; j++)
    {
        for (int i = 1; i < cfg.Nx; i++)
        {
            vtkFile << f.ux[i] << " " << f.vy[j - 1] << " 0.0\n";
        }
    }
    vtkFile << "POINT_DATA " << ((cfg.Nx - 1) * (cfg.Ny - 1)) << "\n";

    // Pressure data
    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 1; j < cfg.Ny; j++)
    {
        for (int i = 1; i < cfg.Nx; i++)
        {
            d1 = sqrt(pow(f.py[j - 1] - f.vy[j - 1], 2) + pow(f.px[i - 1] - f.ux[i], 2));
            d2 = sqrt(pow(f.py[j - 1] - f.vy[j - 1], 2) + pow(f.px[i] - f.ux[i], 2));
            d3 = sqrt(pow(f.py[j] - f.vy[j - 1], 2) + pow(f.px[i - 1] - f.ux[i], 2));
            d4 = sqrt(pow(f.py[j] - f.vy[j - 1], 2) + pow(f.px[i] - f.ux[i], 2));

            w1 = (1 / d1) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            w2 = (1 / d2) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            w3 = (1 / d3) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            w4 = (1 / d4) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            vtkFile << (w1 * f.p[j - 1][i - 1] + w2 * f.p[j - 1][i] + w3 * f.p[j][i - 1] + w4 * f.p[j][i]) << "\n";
        }
    }

    vtkFile << "VECTORS velocity float\n";
    for (int j = 1; j < cfg.Ny; j++)
    {
        for (int i = 1; i < cfg.Nx; i++)
        {
            d1 = f.uy[j - 1] - f.vy[j - 1];
            d2 = f.vy[j - 1] - f.uy[j];
            d3 = f.ux[i] - f.vx[i - 1];
            d4 = f.vx[i] - f.ux[i];
            w1 = (1 / d1) / ((1 / d1) + (1 / d2));
            w2 = (1 / d2) / ((1 / d1) + (1 / d2));
            w3 = (1 / d3) / ((1 / d3) + (1 / d4));
            w4 = (1 / d4) / ((1 / d3) + (1 / d4));
            vtkFile << w1 * f.u[j - 1][i] + w2 * f.u[j][i] << " " << w3 * f.v[j - 1][i - 1] + w4 * f.v[j - 1][i] << " 0.0\n";
        }
    }
    vtkFile.close();
    cout << "VTK file successfully written: " << filename << endl;
}
int main()
{
    auto start = std::chrono::steady_clock::now();

    solverConfig cfg;
    cfg.density = 1.0;
    cfg.Re = 1000;
    cfg.kinematicViscosity = 1.0 / cfg.Re;
    cfg.dynamicViscosity = cfg.kinematicViscosity * cfg.density;

    cfg.uTopWall = 1.0;
    cfg.uBottomWall = 0.0;
    cfg.uLeftWall = 0.0;
    cfg.uRightWall = 0.0;

    cfg.pressureOutlet = 0;
    cfg.pressureInlet = 1;

    cfg.vLeftWall = 0;
    cfg.vRightWall = 0.0;
    cfg.vTopWall = 0.0;
    cfg.vBottomWall = 0.0;

    cfg.Nx = 13;     // Number of cells in the x-direction, including ghost cells
    cfg.Ny = 13;     // Number of cells in the y-direction, including ghost cells
    cfg.lengthX = 1; // Length of the domain in the x-direction, not including ghost cells
    cfg.lengthY = 1; // Length of the domain in the y-direction, not including ghost cells
    // cfg.hx = cfg.lengthX / (cfg.Nx - 2); // Grid spacing in the x-direction
    // cfg.hy = cfg.lengthY / (cfg.Ny - 2); // Grid spacing in the y-direction
    cfg.poissonTolerance = 1e-4; // Tolerance for the Poisson equation solver
    cfg.timeToleranceU = 1e-12;  // Tolerance for the time loop convergence
    cfg.timeToleranceV = 1e-12;  // Tolerance for the time loop convergence
    cfg.timeToleranceP = 1e-12;  // Tolerance for the time loop convergence

    cfg.startTime = 0.0;   // Start time of the simulation
    cfg.endTime = 10000.0; // End time of the simulation. Since we have a convergence criteria, this is set to a large number
    // cfg.firstCellHeightY = 0.01;          // y+ value for the first cell height

    cfg.numberOfMeshLayersN = 0; // Number of mesh layers for the non-uniform mesh in the north direction
    cfg.numberOfMeshLayersS = 0; // Number of mesh layers for the non-uniform mesh in the south direction
    cfg.numberOfMeshLayersW = 0; // Number of mesh layers for the non-uniform mesh in the west direction
    cfg.numberOfMeshLayersE = 0; // Number of mesh layers for the non-uniform mesh in the east direction
    cfg.growthRateN = 1.015;     // Growth rate for the non-uniform mesh in the north direction
    cfg.growthRateS = 1.015;     // Growth rate for the non-uniform mesh in the south direction
    cfg.growthRateW = 1.015;     // Growth rate for the non-uniform mesh in the west direction
    cfg.growthRateE = 1.015;     // Growth rate for the non-uniform mesh in the east direction
    nonUniformMesh(cfg);         // Create the non-uniform mesh
    fields f(cfg.Nx, cfg.Ny);

    initialization(f, cfg); // Initialize the fields

    createCoordinates(f, cfg);  // Create the coordinates of the nodes of cells
    smallestCellHeight(f, cfg); // Calculate the smallest cell height in the domain
    maximumVelocity(f, cfg);
    createAreaU(f, cfg);
    createAreaP(f, cfg);
    createAreaV(f, cfg);

    createVolumeU(f, cfg); // Create the volume of the cells
    createVolumeV(f, cfg); // Create the volume of the cells

    setBoundaryConditionsLDC(1, f.u, cfg);
    setBoundaryConditionsLDC(2, f.v, cfg);
    setBoundaryConditionsLDC(0, f.p, cfg);
    cfg.timeStepSize = min((pow(cfg.smallestCellHeightX, 2) * pow(cfg.smallestCellHeightY, 2)) / (2 * cfg.dynamicViscosity * (pow(cfg.smallestCellHeightX, 2) + pow(cfg.smallestCellHeightY, 2))), 2 * cfg.dynamicViscosity / (pow(cfg.maxUvelocity, 2) + pow(cfg.maxVvelocity, 2)));

    // for (int j = 0; j < cfg.Nx; j++)
    // {
    //     cout << f.vx[j] << endl;
    // }
    // cout << endl;

    int n = 0; // counter for the time steps
    for (double t = cfg.startTime; t <= cfg.endTime; t = t + cfg.timeStepSize)
    {

        cout << "Time step: " << t << endl;

        vector<vector<double>> pPrev = f.p;
        vector<vector<double>> uPrev = f.u;
        vector<vector<double>> vPrev = f.v;

        predictor(f, cfg);
        pressurePoisson(f, cfg);
        setBoundaryConditionsLDC(0, f.p, cfg);
        corrector(f, cfg);
        swapFields(f, cfg);
        setBoundaryConditionsLDC(1, f.u, cfg);
        setBoundaryConditionsLDC(2, f.v, cfg);
        setBoundaryConditionsLDC(0, f.p, cfg);
        updateTimeStepSize(f, cfg);
        double residualU = checkConvergence(f.u, uPrev, cfg);
        double residualV = checkConvergence(f.v, vPrev, cfg);
        double residualP = checkConvergence(f.p, pPrev, cfg);
        // std::this_thread::sleep_for(std::chrono::milliseconds(200));
        cout << (f.u[(cfg.Ny) / 2][1 + (cfg.Nx) / 2]) << endl;
        cout << "Time step size: " << cfg.timeStepSize << endl;
        cout << "Time: " << t << " U velocity: " << residualU << " V velocity: " << residualV << " pressure: " << residualP << " n: " << n << endl;
        if (residualU < cfg.timeToleranceU && residualV < cfg.timeToleranceV && residualP < cfg.timeToleranceP)
        {
            cout << "Converged at time: " << t << endl;
            break;
        }
        n++;
    }

    cout << "First cell height in Y is " << cfg.firstCellHeightY << endl;
    cout << "First cell height in X is " << cfg.firstCellHeightX << endl;

    cout << "Inflation layer thickness in North face is " << cfg.inflationLayerThicknessN << endl;
    cout << "Inflation layer thickness in South face is " << cfg.inflationLayerThicknessS << endl;
    cout << "Inflation layer thickness in West face is " << cfg.inflationLayerThicknessW << endl;
    cout << "Inflation layer thickness in East face is " << cfg.inflationLayerThicknessE << endl;

    cout << "Number of cells in the x-direction is " << cfg.Nx << endl;
    cout << "Number of cells in the y-direction is " << cfg.Ny << endl;

    cout << "Time step size is " << cfg.timeStepSize << endl;

    std::cout << "\nEnd of the main function is reached. Stopping.\n\n";

    // auto end = std::chrono::steady_clock::now();
    // std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;
    for (int j = 0; j < cfg.Ny; j++)
    {
        // cout << f.py[j] << " " << 0.5 * (f.p[j][(cfg.Nx) / 2] + f.p[j][(cfg.Nx - 2) / 2]) << endl;
        cout << f.py[j] << " " << f.p[j][(cfg.Nx - 1) / 2] << endl;
    }
    cout << endl;

    // for (int j = 0; j < cfg.Ny; j++)
    // {
    //     for (int i = 0; i < cfg.Nx; i++)
    //     {
    //         cout << f.p[j][i] << " ";
    //     }
    //     cout << endl;
    // }

    cout << "Smallest cell height in X is " << cfg.smallestCellHeightX << endl;
    cout << "Smallest cell height in Y is " << cfg.smallestCellHeightY << endl;
    vertex_writeVTKFile(f, cfg);

    return 0;
}
