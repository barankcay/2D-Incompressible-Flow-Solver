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
    vector<double> x; // x coordinate of u velocity field
    vector<double> y; // y coordinate of u velocity field

    vector<double> xm; // x coordinate of v velocity field
    vector<double> ym; // y coordinate of v velocity field

    vector<vector<double>> areaN;
    vector<vector<double>> areaS;
    vector<vector<double>> areaE;
    vector<vector<double>> areaW;

    vector<vector<double>> volume;

    vector<vector<double>> u;     // u velocity field
    vector<vector<double>> uStar; // intermediate u velocity field
    vector<vector<double>> uNew;  // next time step u velocity field
    vector<vector<double>> v;     // v velocity field
    vector<vector<double>> vStar; // intermediate v velocity field
    vector<vector<double>> vNew;  // next time step v velocity field
    vector<vector<double>> p;     // pressure field

    fields(int Nx, int Ny)
    {
        x.resize(Nx - 1);
        y.resize(Ny - 1);

        xm.resize(Nx - 2);
        ym.resize(Ny - 2);

        u.resize(Ny, vector<double>(Nx));
        v.resize(Ny, vector<double>(Nx));
        p.resize(Ny, vector<double>(Nx));

        uStar.resize(Ny, vector<double>(Nx));
        vStar.resize(Ny, vector<double>(Nx));
        uNew.resize(Ny, vector<double>(Nx));
        vNew.resize(Ny, vector<double>(Nx));

        areaN.resize(Ny - 2, vector<double>(Nx - 2));
        areaS.resize(Ny - 2, vector<double>(Nx - 2));
        areaE.resize(Ny - 2, vector<double>(Nx - 2));
        areaW.resize(Ny - 2, vector<double>(Nx - 2));

        volume.resize(Ny - 2, vector<double>(Nx - 2));
    }
};
void setBoundaryConditionsLDC(int b, vector<vector<double>> &M, solverConfig &cfg)
{

    if (b == 0)
    {
        for (int j = cfg.Ny - 2; j >= 1; j--)
        {
            M[j][0] = M[j][1];                   // Left wall, dp/dx = 0
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2]; // Right wall, dp/dx = 0
        }
        for (int i = cfg.Nx - 2; i >= 1; i--)
        {
            M[0][i] = M[1][i];                   // Bottom wall, dp/dy = 0
            M[cfg.Ny - 1][i] = M[cfg.Ny - 2][i]; // Top wall, dp/dy = 0
        }
    }
    else if (b == 1)
    {
        for (int j = cfg.Ny - 2; j >= 1; j--)
        {
            M[j][0] = 2 * cfg.uLeftWall - M[j][1]; // Left wall,ghost cell,  u = 0

            M[j][cfg.Nx - 1] = 2 * cfg.uRightWall - M[j][cfg.Nx - 2]; // Right wall, u = 0
        }
        for (int i = cfg.Nx - 2; i >= 1; i--)
        {
            M[0][i] = 2 * cfg.uTopWall - M[1][i];                      // Top wall, ghost cell
            M[cfg.Ny - 1][i] = 2 * cfg.uBottomWall - M[cfg.Ny - 2][i]; // Bottom wall, ghost cell
        }
    }
    else if (b == 2)
    {

        for (int i = cfg.Nx - 2; i >= 1; i--)
        {
            M[0][i] = 2 * cfg.vTopWall - M[1][i];                      // Top wall, v = 0
            M[cfg.Ny - 1][i] = 2 * cfg.vBottomWall - M[cfg.Ny - 2][i]; // Bottom wall, v = 0
        }
        for (int j = cfg.Ny - 2; j >= 1; j--)
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
// function to initialize the fields

void setBoundaryConditionsCF_pressureDriven(int b, vector<vector<double>> &M, solverConfig &cfg)
{
    if (b == 0)
    {
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * cfg.pressureInlet - M[j][1];                    // pressure driven inlet
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
            M[j][0] = M[j][1];
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
            M[0][i] = 2 * cfg.vTopWall - M[1][i]; // Top wall, v = 0
            M[cfg.Ny - 1][i] = 2 * cfg.vBottomWall - M[cfg.Ny - 2][i];
        }
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            // M[j][0] = 2 * cfg.vLeftWall - M[j][1]; // Left wall, ghost cell
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
            M[j][0] = M[j][1];                                            // Left wall, dp/dx = 0
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
            M[j][0] = 2 * cfg.uLeftWall - M[j][1]; // Left wall,ghost cell,  u = 0
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2];   // Right wall, du/dx = 0
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
            M[0][i] = 2 * cfg.vTopWall - M[1][i]; // Top wall, v = 0
            M[cfg.Ny - 1][i] = 2 * cfg.vBottomWall - M[cfg.Ny - 2][i];
        }
        for (int j = cfg.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * cfg.vLeftWall - M[j][1]; // Left wall, ghost cell
            M[j][cfg.Nx - 1] = M[j][cfg.Nx - 2];   // Right wall, ghost cell
        }
    }

    M[0][0] = 0.5 * (M[0][1] + M[1][0]);                                                       // Top left corner
    M[0][cfg.Nx - 1] = 0.5 * (M[0][cfg.Nx - 2] + M[1][cfg.Nx - 1]);                            // Top right corner
    M[cfg.Ny - 1][0] = 0.5 * (M[cfg.Ny - 2][0] + M[cfg.Ny - 1][1]);                            // Bottom left corner
    M[cfg.Ny - 1][cfg.Nx - 1] = 0.5 * (M[cfg.Ny - 1][cfg.Nx - 2] + M[cfg.Ny - 2][cfg.Nx - 1]); // Bottom right corner
}
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
    for (int j = 0; j < cfg.numberOfMeshLayersN + 1; j++)
    {

        f.y[j] = cfg.lengthY - (cfg.firstCellHeightY + (cfg.lengthY - f.y[j - 1]) * cfg.growthRateN); // y coordinate of u velocity field
        f.y[0] = cfg.lengthY;
    }
    for (int i = 0; i < cfg.numberOfMeshLayersS + 1; i++)
    {
        f.y[cfg.Ny - 2 - i] = cfg.firstCellHeightY + (f.y[cfg.Ny - 1 - i]) * cfg.growthRateS;
        f.y[cfg.Ny - 2] = 0;
    }
    for (int i = 0; i < cfg.numberOfMeshLayersW + 1; i++)
    {
        f.x[i] = cfg.firstCellHeightX + (f.x[i - 1]) * cfg.growthRateW; // x coordinate of u velocity field
        f.x[0] = 0;
    }

    for (int i = 0; i < cfg.numberOfMeshLayersE + 1; i++)
    {
        f.x[cfg.Nx - 2 - i] = cfg.lengthX - (cfg.firstCellHeightX + (cfg.lengthX - f.x[cfg.Nx - 1 - i]) * cfg.growthRateE); // x coordinate of u velocity field
        f.x[cfg.Nx - 2] = cfg.lengthX;
    }

    cout << "Grid spacing in the y-direction: " << cfg.hy << endl;
    cout << "Grid spacing in the x-direction: " << cfg.hx << endl;
    for (int j = cfg.numberOfMeshLayersN + 1; j < cfg.Ny - 2 - cfg.numberOfMeshLayersS; j++)
    {
        f.y[j] = f.y[j - 1] - cfg.hy; // y coordinate of v velocity field
    }
    for (int i = cfg.numberOfMeshLayersW + 1; i < cfg.Nx - 2 - cfg.numberOfMeshLayersE; i++)
    {
        f.x[i] = f.x[i - 1] + cfg.hx; // x coordinate of u velocity field
    }

    for (int j = 0; j < cfg.Ny - 2; j++)
    {
        for (int i = 0; i < cfg.Nx - 2; i++)
        {
            f.xm[i] = (f.x[i] + f.x[i + 1]) / 2; // x coordinate of v velocity field

            f.ym[j] = (f.y[j] + f.y[j + 1]) / 2; // y coordinate of v velocity field
        }
    }

    f.x.push_back(f.x[cfg.Nx - 2] + (f.x[cfg.Nx - 2] - f.x[cfg.Nx - 3])); // x coordinate of v velocity field
    f.x.insert(f.x.begin(), (f.x[0] - (f.x[1] - f.x[0])));                // x coordinate of v velocity field

    f.y.push_back(f.y[cfg.Ny - 2] - (f.y[cfg.Ny - 3] - f.y[cfg.Ny - 2])); // y coordinate of v velocity field
    f.y.insert(f.y.begin(), f.y[0] + (f.y[0] - f.y[1]));                  // y coordinate of v velocity field

    f.xm.push_back((f.x[cfg.Nx] + f.x[cfg.Nx - 1]) / 2); // x coordinate of v velocity field
    f.xm.insert(f.xm.begin(), (f.x[0] + f.x[1]) / 2);    // x coordinate of v velocity field

    f.ym.push_back((f.y[cfg.Ny] + f.y[cfg.Ny - 1]) / 2); // y coordinate of v velocity field
    f.ym.insert(f.ym.begin(), (f.y[0] + f.y[1]) / 2);    // y coordinate of v velocity field
}

void smallestCellHeight(fields &f, solverConfig &cfg)
{
    cfg.smallestCellHeightY = f.y[1] - f.y[2];
    cfg.smallestCellHeightX = f.x[2] - f.x[1];

    for (int i = 2; i < cfg.Nx - 1; i++)
    {
        if (f.x[i + 1] - f.x[i] < cfg.smallestCellHeightX)
        {
            cfg.smallestCellHeightX = f.x[i + 1] - f.x[i];
        }
    }
    for (int j = 2; j < cfg.Ny - 1; j++)
    {
        if (f.y[j] - f.y[j + 1] < cfg.smallestCellHeightY)
        {
            cfg.smallestCellHeightY = f.y[j] - f.y[j + 1];
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
void createArea(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.areaE[j][i] = 1 * (f.y[j + 1] - f.y[j + 2]);
            f.areaW[j][i] = 1 * (f.y[j + 1] - f.y[j + 2]); // area of the cell in the x-direction
            f.areaN[j][i] = 1 * (f.x[i + 2] - f.x[i + 1]); // area of the cell in the y-direction
            f.areaS[j][i] = 1 * (f.x[i + 2] - f.x[i + 1]); // area of the cell in the x-direction
        }
    }
}

void createVolume(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.volume[j][i] = 1 * f.areaE[j][i] * f.areaN[j][i]; // volume of the cell in the x-direction
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
            double vn = f.v[j][i] + (f.v[j - 1][i] - f.v[j][i]) * (f.y[j] - f.ym[j]) / (f.ym[j - 1] - f.ym[j]); // v velocity field
            double un = f.u[j][i] + (f.u[j - 1][i] - f.u[j][i]) * (f.y[j] - f.ym[j]) / (f.ym[j - 1] - f.ym[j]);
            double ue = f.u[j][i] + (f.u[j][i + 1] - f.u[j][i]) * (f.x[i + 1] - f.xm[i]) / (f.xm[i + 1] - f.xm[i]);     // u velocity field
            double uw = f.u[j][i - 1] + (f.u[j][i] - f.u[j][i - 1]) * (f.x[i] - f.xm[i - 1]) / (f.xm[i] - f.xm[i - 1]); // u velocity field
            double vs = f.v[j + 1][i] + (f.v[j][i] - f.v[j + 1][i]) * (f.y[i + 1] - f.ym[i + 1]) / (f.ym[i] - f.ym[i + 1]);
            double us = f.u[j + 1][i] + (f.u[j][i] - f.u[j + 1][i]) * (f.y[i + 1] - f.ym[i + 1]) / (f.ym[i] - f.ym[i + 1]);
            double ve = f.v[j][i] + (f.v[j][i + 1] - f.v[j][i]) * (f.x[i + 1] - f.xm[i]) / (f.xm[i + 1] - f.xm[i]);     // v velocity field
            double vw = f.v[j][i - 1] + (f.v[j][i] - f.v[j][i - 1]) * (f.x[i] - f.xm[i - 1]) / (f.xm[i] - f.xm[i - 1]); // v velocity field
            double An = f.areaN[j - 1][i - 1];
            double As = f.areaS[j - 1][i - 1];
            double Ae = f.areaE[j - 1][i - 1];
            double Aw = f.areaW[j - 1][i - 1];
            double advectionU = vn * un * An - vs * us * As + ue * ue * Ae - uw * uw * Aw;
            double diffusionU = cfg.kinematicViscosity * ((Ae * (f.u[j][i + 1] - f.u[j][i]) / (f.xm[i + 1] - f.xm[i])) - (Aw * (f.u[j][i] - f.u[j][i - 1]) / (f.xm[i] - f.xm[i - 1])) + (An * (f.u[j - 1][i] - f.u[j][i]) / (f.ym[j - 1] - f.ym[j])) - (As * (f.u[j][i] - f.u[j + 1][i]) / (f.ym[j] - f.ym[j + 1]))); // diffusion term
            double advectionV = vn * vn * An - vs * vs * As + ue * ve * Ae - uw * vw * Aw;
            double diffusionV = cfg.kinematicViscosity * ((Ae * (f.v[j][i + 1] - f.v[j][i]) / (f.xm[i + 1] - f.xm[i])) - (Aw * (f.v[j][i] - f.v[j][i - 1]) / (f.xm[i] - f.xm[i - 1])) + (An * (f.v[j - 1][i] - f.v[j][i]) / (f.ym[j - 1] - f.ym[j])) - (As * (f.v[j][i] - f.v[j + 1][i]) / (f.ym[j] - f.ym[j + 1])));

            f.vStar[j][i] = f.v[j][i] + (cfg.timeStepSize / f.volume[j - 1][i - 1]) * (-advectionV + diffusionV);
            f.uStar[j][i] = f.u[j][i] + (cfg.timeStepSize / f.volume[j - 1][i - 1]) * (-advectionU + diffusionU); // u velocity field
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

                double delX1 = f.xm[i - 1] - f.xm[i];
                double delX2 = f.xm[i + 1] - f.xm[i];

                double delY1 = f.ym[j + 1] - f.ym[j];
                double delY2 = f.ym[j - 1] - f.ym[j];
                double velocStarGrad = (1 / cfg.timeStepSize) * (((f.uStar[j][i + 1] - f.uStar[j][i - 1]) / (f.xm[i + 1] - f.xm[i - 1])) + ((f.vStar[j - 1][i] - f.vStar[j + 1][i]) / (f.ym[j - 1] - f.ym[j + 1])));

                double A = delX1 * delX2 * delX2 - delX2 * delX1 * delX1;

                double B = delY1 * delY2 * delY2 - delY2 * delY1 * delY1;

                double pCoeff = (2 * (delX1 - delX2) / A) + (2 * (delY1 - delY2) / B);

                f.p[j][i] = ((2 * f.p[j][i + 1] * delX1 / A) - (2 * f.p[j][i - 1] * delX2 / A) + (2 * f.p[j - 1][i] * delY1 / B) - (2 * f.p[j + 1][i] * delY2 / B) - velocStarGrad) / pCoeff;
                residual += abs(f.p[j][i] - pOld);
                // cout << A << " " << B << " " << pCoeff << " " << endl;
            }
        }
        setBoundaryConditionsLDC(0, f.p, cfg);
        residual /= (cfg.Nx * cfg.Ny);
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
            double Pw = f.p[j][i - 1] + (f.p[j][i] - f.p[j][i - 1]) * (f.x[i] - f.xm[i - 1]) / (f.xm[i] - f.xm[i - 1]);     // p velocity field
            double Pe = f.p[j][i] + (f.p[j][i + 1] - f.p[j][i]) * (f.x[i + 1] - f.xm[i]) / (f.xm[i + 1] - f.xm[i]);         // p velocity field
            double Ps = f.p[j + 1][i] + (f.p[j][i] - f.p[j + 1][i]) * (f.y[i + 1] - f.ym[i + 1]) / (f.ym[i] - f.ym[i + 1]); // p velocity field
            double Pn = f.p[j][i] + (f.p[j - 1][i] - f.p[j][i]) * (f.y[j] - f.ym[j]) / (f.ym[j - 1] - f.ym[j]);             // p velocity field
            double An = f.areaN[j - 1][i - 1];
            double As = f.areaS[j - 1][i - 1];
            double Ae = f.areaE[j - 1][i - 1];
            double Aw = f.areaW[j - 1][i - 1];
            f.uNew[j][i] = f.uStar[j][i] + (cfg.timeStepSize / f.volume[j - 1][i - 1]) * ((Pw * Aw - Pe * Ae) / (cfg.density));
            f.vNew[j][i] = f.vStar[j][i] + (cfg.timeStepSize / f.volume[j - 1][i - 1]) * ((Ps * As - Pn * An) / (cfg.density));
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
    return residual /= (cfg.Nx * cfg.Ny);
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

void cell_writeVTKFile(const fields &f, const solverConfig &cfg)
{
    std::string filename = "cell_colocated_results.vtk";
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
        vtkFile << f.x[i] << "\n"; // Ensure f.px has Nx+1 entries
    }

    vtkFile << "Y_COORDINATES " << (cfg.Ny + 1) << " float\n";
    for (int j = 0; j <= cfg.Ny; j++)
    {
        vtkFile << f.y[j] << "\n"; // Ensure f.py has Ny+1 entries
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

void point_writeVTKFile(fields &f, const solverConfig &cfg)
{
    std::string filename = "point_colocated_results.vtk";
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
    vtkFile << "DIMENSIONS " << cfg.Nx << " " << cfg.Ny << " 1\n";

    // Points at cell centers
    vtkFile << "POINTS " << (cfg.Nx * cfg.Ny) << " float\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.xm[i] << " " << f.ym[j] << " 0.0\n";
        }
    }

    // POINT_DATA (since we're treating cell centers as points)
    vtkFile << "POINT_DATA " << (cfg.Nx * cfg.Ny) << "\n";

    // Pressure data
    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < cfg.Ny; j++)
    {
        for (int i = 0; i < cfg.Nx; i++)
        {
            vtkFile << f.p[j][i] << "\n";
        }
    }

    // Velocity data
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

    std::string filename = "coloc_vertex_results.vtk";
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
            vtkFile << f.x[i] << " " << f.y[j] << " 0.0\n";
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
            d1 = sqrt(pow(f.x[i] - f.xm[i - 1], 2) + pow(f.y[j] - f.ym[j - 1], 2));
            d2 = sqrt(pow(f.x[i] - f.xm[i], 2) + pow(f.y[j] - f.ym[j - 1], 2));
            d3 = sqrt(pow(f.x[i] - f.xm[i - 1], 2) + pow(f.y[j] - f.ym[j], 2));
            d4 = sqrt(pow(f.x[i] - f.xm[i], 2) + pow(f.y[j] - f.ym[j], 2));

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
            d1 = sqrt(pow(f.x[i] - f.xm[i - 1], 2) + pow(f.y[j] - f.ym[j - 1], 2));
            d2 = sqrt(pow(f.x[i] - f.xm[i], 2) + pow(f.y[j] - f.ym[j - 1], 2));
            d3 = sqrt(pow(f.x[i] - f.xm[i - 1], 2) + pow(f.y[j] - f.ym[j], 2));
            d4 = sqrt(pow(f.x[i] - f.xm[i], 2) + pow(f.y[j] - f.ym[j], 2));

            w1 = (1 / d1) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            w2 = (1 / d2) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            w3 = (1 / d3) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));
            w4 = (1 / d4) / ((1 / d1) + (1 / d2) + (1 / d3) + (1 / d4));

            vtkFile << (w1 * f.u[j - 1][i - 1] + w2 * f.u[j - 1][i] + w3 * f.u[j][i - 1] + w4 * f.u[j][i]) << " " << (w1 * f.v[j - 1][i - 1] + w2 * f.v[j - 1][i] + w3 * f.v[j][i - 1] + w4 * f.v[j][i]) << " 0.0\n";
        }
    }

    vtkFile.close();
    std::cout << "VTK file successfully written: " << filename << std::endl;
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
    cfg.startTime = 0.0;         // Start time of the simulation
    cfg.endTime = 10000.0;       // End time of the simulation. Since we have a convergence criteria, this is set to a large number
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

    initialization(f, cfg);    // Initialize the fields
    createCoordinates(f, cfg); // Create the coordinates of the nodes of cells

    smallestCellHeight(f, cfg); // Calculate the smallest cell height in the domain
    maximumVelocity(f, cfg);    // Calculate the maximum velocity in the domain
    createArea(f, cfg);
    createVolume(f, cfg); // Create the volume of the cells

    setBoundaryConditionsLDC(1, f.u, cfg);
    setBoundaryConditionsLDC(2, f.v, cfg);
    setBoundaryConditionsLDC(0, f.p, cfg);
    cfg.timeStepSize = min((pow(cfg.hx, 2) * pow(cfg.hy, 2)) / (2 * cfg.dynamicViscosity * (pow(cfg.hx, 2) + pow(cfg.hy, 2))), 2 * cfg.dynamicViscosity / (pow(cfg.uTopWall, 2) + pow(cfg.uTopWall, 2)));
    cout << "Initial time step size: " << cfg.timeStepSize << endl;
    // for (int j = 0; j < cfg.Nx; j++)
    // {
    //     cout << f.xm[j] << endl;
    // }
    // cout << endl;

    // this_thread::sleep_for(chrono::milliseconds(5000)); // Sleep for 2 seconds to allow the user to see the initial values

    int n = 0; // counter for the time steps
    for (double t = cfg.startTime; t <= cfg.endTime; t = t + cfg.timeStepSize)
    {

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
        // updateTimeStepSize(f, cfg);
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

    cout << cfg.timeStepSize << endl;

    // for (int j = 0; j < cfg.Ny; j++)
    // {
    //     cout << f.x[j] << " " << (f.u[j][(cfg.Nx) / 2]) << endl;
    // }
    // cout << endl;

    // std::cout << "\nEnd of the main function is reached. Stopping.\n\n";

    // auto end = std::chrono::steady_clock::now();
    // std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;

    cout << "First cell height in Y is " << cfg.firstCellHeightY << endl;
    cout << "First cell height in X is " << cfg.firstCellHeightX << endl;

    cout << "uniform cell height in Y is " << cfg.hy << endl;
    cout << "uniform cell height in X is " << cfg.hx << endl;
    cout << "Inflation layer thickness in North face is " << cfg.inflationLayerThicknessN << endl;
    cout << "Inflation layer thickness in South face is " << cfg.inflationLayerThicknessS << endl;
    cout << "Inflation layer thickness in West face is " << cfg.inflationLayerThicknessW << endl;
    cout << "Inflation layer thickness in East face is " << cfg.inflationLayerThicknessE << endl;
    cout << "heyo of " << min(cfg.hy, cfg.hx) << endl;
    cout << "Number of cells in the x-direction is " << cfg.Nx << endl;
    cout << "Number of cells in the y-direction is " << cfg.Ny << endl;
    cout << endl;

    // for (int i = 0; i < cfg.Nx; i++)
    // {
    //     cout << f.xm[i] << " ";
    // }

    cout << "Time step size is " << cfg.timeStepSize << endl;

    std::cout << "\nEnd of the main function is reached. Stopping.\n\n";

    for (int i = 0; i < cfg.Ny; i++)
    {
        cout << f.ym[i] << " " << f.p[i][(cfg.Nx - 1) / 2] << endl;
    }
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
