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
    double vLeftWall;   // v Velocity at the left wall
    double vRightWall;  // v Velocity at the right wall

    double uInlet; // u velocity at the inlet
    double vInlet; // v velocity at the inlet

    double pressureOutlet; // Pressure at the outlet

    double courantNumber; // Courant number used to calculate the time step size

    double poissonTolerance; // Tolerance for the Poisson equation solver
    double UtimeTolerance;   // Tolerance for the time loop convergence
    double VtimeTolerance;   // Tolerance for the time loop convergence
    double PtimeTolerance;   // Tolerance for the time loop convergence

    // non-uniform mesh properties

    double firstCellHeightX; // height of the first cell for the non-uniform mesh
    double firstCellHeightY;
    double inflationLayerThicknessX; // boundary layer thickness for the non-uniform mesh
    double inflationLayerThicknessY;
    double growthRateX; // growth rate for the non-uniform mesh
    double growthRateY;

    int numberOfMeshLayersX; // number of mesh layers for the non-uniform mesh
    int numberofMeshLayersY; // number of mesh layers for the non-uniform mesh
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

        areaUn.resize(Ny - 2, vector<double>(Nx - 3));
        areaUs.resize(Ny - 2, vector<double>(Nx - 3));
        areaUe.resize(Ny - 2, vector<double>(Nx - 3));
        areaUw.resize(Ny - 2, vector<double>(Nx - 3));

        areaVn.resize(Ny - 3, vector<double>(Nx - 2));
        areaVs.resize(Ny - 3, vector<double>(Nx - 2));
        areaVe.resize(Ny - 3, vector<double>(Nx - 2));
        areaVw.resize(Ny - 3, vector<double>(Nx - 2));

        areaPn.resize(Ny - 2, vector<double>(Nx - 2));
        areaPs.resize(Ny - 2, vector<double>(Nx - 2));
        areaPe.resize(Ny - 2, vector<double>(Nx - 2));
        areaPw.resize(Ny - 2, vector<double>(Nx - 2));

        volumeU.resize(Ny - 2, vector<double>(Nx - 3));
        volumeV.resize(Ny - 3, vector<double>(Nx - 2));
    }
};
void setBoundaryConditions(int b, vector<vector<double>> &M, solverConfig &cfg)
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

void nonUniformMeshX(solverConfig &cfg)
{

    cfg.inflationLayerThicknessX = cfg.firstCellHeightX;
    for (int i = 0; i < cfg.numberOfMeshLayersX - 1; i++)
    {
        cfg.inflationLayerThicknessX = cfg.firstCellHeightX + cfg.inflationLayerThicknessX * cfg.growthRateX; // boundary layer thickness
    }
    cfg.inflationLayerThicknessX = 2 * cfg.inflationLayerThicknessX;

    if (cfg.inflationLayerThicknessX > cfg.lengthX)
    {
        cout << "Inflation layer thickness in X is " << cfg.inflationLayerThicknessX << endl;
        cout << "Inflation layer thickness in X is greater than the length of the domain in the x-direction" << endl;
        cout << "Please reduce the number of mesh layers or the growth rate" << endl;
        exit(1);
    }
    // double the layer because of left and right wall

    double largestBLmesh = cfg.firstCellHeightX * pow(cfg.growthRateX, cfg.numberOfMeshLayersX - 1);
    cfg.Nx = 2 + 2 * cfg.numberOfMeshLayersX + ((cfg.lengthX - cfg.inflationLayerThicknessX) / (largestBLmesh * cfg.growthRateX)); // Number of cells in the y-direction, including ghost cells
    cout << "Number of cells in the x-direction: " << cfg.Nx << endl;
}

void nonUniformMeshY(solverConfig &cfg)
{

    cfg.inflationLayerThicknessY = cfg.firstCellHeightY;
    for (int i = 0; i < cfg.numberofMeshLayersY - 1; i++)
    {
        cfg.inflationLayerThicknessY = cfg.firstCellHeightY + cfg.inflationLayerThicknessY * cfg.growthRateY; // boundary layer thickness
    }
    cfg.inflationLayerThicknessY = 2 * cfg.inflationLayerThicknessY;
    if (cfg.inflationLayerThicknessY > cfg.lengthY)
    {
        cout << "Inflation layer thickness in Y is " << cfg.inflationLayerThicknessY << endl;
        cout << "Inflation layer thickness in y is greater than the length of the domain in the y-direction" << endl;
        cout << "Please reduce the number of mesh layers or the growth rate" << endl;
        exit(1);
    }
    // double the layer because of left and right wall

    double largestBLmesh = cfg.firstCellHeightY * pow(cfg.growthRateY, cfg.numberofMeshLayersY - 1);
    cfg.Ny = 2 + 2 * cfg.numberofMeshLayersY + ((cfg.lengthY - cfg.inflationLayerThicknessY) / (largestBLmesh * cfg.growthRateY)); // Number of cells in the y-direction, including ghost cells
    cout << "Number of cells in the Y-direction: " << cfg.Ny << endl;
}
// function to create the coordinates of the nodes of cells
void createCoordinates(fields &f, solverConfig &cfg)
{
    f.ux[0] = 0;
    f.vx[0] = cfg.firstCellHeightX / 2;
    f.vx[cfg.Nx - 3] = cfg.lengthX - cfg.firstCellHeightX / 2;
    f.ux[cfg.Nx - 3] = cfg.lengthX - cfg.firstCellHeightX;

    f.uy[0] = cfg.lengthY - cfg.firstCellHeightY / 2;
    f.vy[0] = cfg.lengthY - cfg.firstCellHeightY;
    f.vy[cfg.Ny - 3] = 0;
    f.uy[cfg.Ny - 3] = cfg.firstCellHeightY / 2;
    for (int i = 1; i < cfg.numberOfMeshLayersX; i++)
    {
        f.ux[i] = cfg.firstCellHeightX + (f.ux[i - 1]) * cfg.growthRateX;                                                     // x coordinate of u velocity field
        f.vx[i] = cfg.firstCellHeightX + (f.vx[i - 1]) * cfg.growthRateX;                                                     // y coordinate of u velocity field
        f.ux[cfg.Nx - 3 - i] = cfg.lengthX - (cfg.firstCellHeightX + (cfg.lengthX - f.ux[cfg.Nx - 2 - i]) * cfg.growthRateX); // x coordinate of u velocity field
        f.vx[cfg.Nx - 3 - i] = cfg.lengthX - (cfg.firstCellHeightX + (cfg.lengthX - f.vx[cfg.Nx - 2 - i]) * cfg.growthRateX); // y coordinate of u velocity field
    }
    for (int i = 1; i < cfg.numberofMeshLayersY; i++)
    {
        f.uy[i] = cfg.lengthY - (cfg.firstCellHeightY + (cfg.lengthY - f.uy[i - 1]) * cfg.growthRateY);
        f.vy[i] = cfg.lengthY - (cfg.firstCellHeightY + (cfg.lengthY - f.vy[i - 1]) * cfg.growthRateY);
        f.uy[cfg.Ny - 3 - i] = cfg.firstCellHeightY + (f.uy[cfg.Ny - 2 - i]) * cfg.growthRateY;
        f.vy[cfg.Ny - 3 - i] = cfg.firstCellHeightY + (f.vy[cfg.Ny - 2 - i]) * cfg.growthRateY;
    }
    cfg.hy = (cfg.lengthY - cfg.inflationLayerThicknessY) / (cfg.Ny - 2 - 2 * cfg.numberofMeshLayersY); // Grid spacing in the y-direction
    cfg.hx = (cfg.lengthX - cfg.inflationLayerThicknessX) / (cfg.Nx - 2 - 2 * cfg.numberOfMeshLayersX); // Grid spacing in the y-direction
    // cfg.hy = cfg.firstCellHeight * pow(cfg.growthRate, cfg.numberOfMeshLayers);
    cout << "Grid spacing in the y-direction: " << cfg.hx << endl;
    for (int j = cfg.numberofMeshLayersY; j < cfg.Ny - 2 - cfg.numberofMeshLayersY; j++)
    {
        for (int i = cfg.numberOfMeshLayersX; i < cfg.Nx - 2 - cfg.numberOfMeshLayersX; i++)
        {
            f.ux[i] = f.ux[i - 1] + cfg.hx;     // x coordinate of u velocity field
            f.uy[j] = f.vy[j - 1] - cfg.hy / 2; // y coordinate of u velocity field

            f.vx[i] = f.ux[i] + 0.5 * cfg.hx; // x coordinate of v velocity field
            f.vy[j] = f.vy[j - 1] - cfg.hy;   // y coordinate of v velocity field
        }
    }

    f.ux.push_back(f.ux[cfg.Nx - 3] + 2 * (f.vx[cfg.Nx - 3] - f.ux[cfg.Nx - 3]));
    f.ux.insert(f.ux.begin(), (f.ux[0] - 2 * (f.vx[0] - f.ux[0])));
    f.uy.push_back(f.uy[cfg.Ny - 3] - 2 * (f.uy[cfg.Ny - 3] - f.vy[cfg.Ny - 3]));
    f.uy.insert(f.uy.begin(), f.uy[0] + 2 * (f.uy[0] - f.vy[0]));

    f.vx.push_back(f.vx[cfg.Nx - 3] + 2 * (f.vx[cfg.Nx - 3] - f.ux[cfg.Nx - 2]));
    f.vx.insert(f.vx.begin(), (f.vx[0] - 2 * (f.vx[0] - f.ux[1])));

    f.vy.push_back(f.vy[cfg.Ny - 3] - 2 * (f.uy[cfg.Ny - 2] - f.vy[cfg.Ny - 3]));
    f.vy.insert(f.vy.begin(), f.vy[0] + 2 * (f.uy[1] - f.vy[0]));

    f.px = f.vx; // x coordinate of pressure field
    f.py = f.uy; // y coordinate of pressure field
}

void createAreaU(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 3; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.areaUn[j][i] = 1 * (f.vx[i + 2] - f.vx[i + 1]); // area of the cell in the x-direction
            f.areaUs[j][i] = 1 * (f.vx[i + 2] - f.vx[i + 1]); // area of the cell in the y-direction
            f.areaUe[j][i] = 1 * (f.vy[j] - f.vy[j + 1]);     // area of the cell in the x-direction
            f.areaUw[j][i] = 1 * (f.vy[j] - f.vy[j + 1]);
        }
    }
}

void createAreaV(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 3; j++)
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
    for (int i = 0; i < cfg.Nx - 3; i++)
    {
        for (int j = 0; j < cfg.Ny - 2; j++)
        {
            f.volumeU[j][i] = 1 * f.areaUn[j][i] * f.areaUs[j][i]; // volume of the cell in the x-direction
        }
    }
}

void createVolumeV(fields &f, solverConfig &cfg)
{
    for (int i = 0; i < cfg.Nx - 2; i++)
    {
        for (int j = 0; j < cfg.Ny - 3; j++)
        {
            f.volumeV[j][i] = 1 * f.areaVn[j][i] * f.areaVs[j][i]; // volume of the cell in the x-direction
        }
    }
}

void predictor(fields &f, solverConfig &cfg)
{
    // Predictor step for the u velocity field
    for (int j = 1; j < cfg.Ny - 1; j++)
    {
        for (int i = 2; i < cfg.Nx - 1; i++)
        {
            double vn = f.v[j - 1][i - 1] + (f.v[j - 1][i] - f.v[j - 1][i - 1]) * (f.ux[i] - f.vx[i - 1]) / (f.vx[i] - f.vx[i - 1]); // v velocity field
            double un = f.u[j][i] + (f.u[j - 1][i] - f.u[j][i]) * (f.vy[j - 1] - f.uy[j]) / (f.uy[j - 1] - f.uy[j]);
            double ue = f.u[j][i] + (f.u[j][i + 1] - f.u[j][i]) * (f.vx[i] - f.ux[i]) / (f.ux[i + 1] - f.ux[i]);             // u velocity field
            double uw = f.u[j][i - 1] + (f.u[j][i] - f.u[j][i - 1]) * (f.vx[i - 1] - f.ux[i - 1]) / (f.ux[i] - f.ux[i - 1]); // u velocity field
            double vs = f.v[j][i - 1] + (f.v[j][i] - f.v[j][i - 1]) * (f.ux[i] - f.vx[i - 1]) / (f.vx[i] - f.vx[i - 1]);
            double us = f.u[j + 1][i] + (f.u[j][i] - f.u[j + 1][i]) * (f.vy[j] - f.uy[j + 1]) / (f.uy[j] - f.uy[j + 1]);
            double An = f.areaUn[j - 1][i - 2];
            double As = f.areaUs[j - 1][i - 2];
            double Ae = f.areaUe[j - 1][i - 2];
            double Aw = f.areaUw[j - 1][i - 2];
            double advection = vn * un * An - vs * us * As + ue * ue * Ae - uw * uw * Aw;
            double diffusion = cfg.kinematicViscosity * ((Ae * (f.u[j][i + 1] - f.u[j][i]) / (f.ux[i + 1] - f.ux[i])) - (Aw * (f.u[j][i] - f.u[j][i - 1]) / (f.ux[i] - f.ux[i - 1])) + (An * (f.u[j - 1][i] - f.u[j][i]) / (f.uy[j - 1] - f.uy[j])) - (As * (f.u[j][i] - f.u[j + 1][i]) / (f.uy[j] - f.uy[j + 1]))); // diffusion term
            f.uStar[j][i] = f.u[j][i] + (cfg.timeStepSize / f.volumeU[j - 1][i - 2]) * (-advection + diffusion);                                                                                                                                                                                                     // u velocity field
        }
    }
    // Predictor step for the v velocity field
    for (int j = 1; j < cfg.Ny - 2; j++)
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
    setBoundaryConditions(1, f.uStar, cfg);
    setBoundaryConditions(2, f.vStar, cfg);
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

                f.p[j][i] = ((2 * f.p[j][i + 1] * delX1 / A) - (2 * f.p[j][i - 1] * delX2 / A) + (2 * f.p[j - 1][i] * delY1 / B) - (2 * f.p[j + 1][i] * delY2 / B) - velocStarGrad) / pCoeff;
                residual += abs(f.p[j][i] - pOld);
            }
        }
        setBoundaryConditions(0, f.p, cfg);
        residual /= (cfg.Nx * cfg.Ny);
        iteration++;
    }
    cout << "Number of iterations: " << iteration << endl;
}

void updateTimeStepSize(fields &f, solverConfig &cfg)
{

    // Calculate the maximum velocity in the domain (including boundaries)
    for (int i = 1; i < cfg.Nx - 1; i++)
    {
        for (int j = 1; j < cfg.Ny - 1; j++)
        {
            double rule1 = (pow(cfg.firstCellHeightX, 2) * pow(cfg.firstCellHeightY, 2)) / (2 * cfg.dynamicViscosity * (pow(cfg.firstCellHeightX, 2) + pow(cfg.firstCellHeightY, 2)));
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
    for (int i = 2; i < cfg.Nx - 1; i++)
    {
        for (int j = 1; j < cfg.Ny - 1; j++)
        {
            f.uNew[j][i] = f.uStar[j][i] + (cfg.timeStepSize / f.volumeU[j - 1][i - 2]) * ((f.p[j][i - 1] * f.areaUw[j - 1][i - 2] - f.p[j][i] * f.areaUe[j - 1][i - 2]) / (cfg.density));
        }
    }

    for (int i = 1; i < cfg.Nx - 1; i++)
    {
        for (int j = 1; j < cfg.Ny - 2; j++)
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
void point_writeVTKFile(fields &f, solverConfig &cfg)
{
    std::string filename = "BL_results.vtk";
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
            vtkFile << f.px[i] << " " << f.py[j] << " 0.0\n";
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

int main()
{
    auto start = std::chrono::steady_clock::now();

    solverConfig cfg;

    cfg.density = 1.0;
    cfg.Re = 100;
    cfg.kinematicViscosity = 1.0 / cfg.Re;
    cfg.dynamicViscosity = cfg.kinematicViscosity * cfg.density;

    cfg.uTopWall = 1;
    cfg.uBottomWall = 0.0;

    cfg.vLeftWall = 0.0;
    cfg.vRightWall = 0.0;

    cfg.Nx = 21;     // Number of cells in the x-direction, including ghost cells
    cfg.Ny = 21;     // Number of cells in the y-direction, including ghost cells
    cfg.lengthX = 1; // Length of the domain in the x-direction, not including ghost cells
    cfg.lengthY = 1; // Length of the domain in the y-direction, not including ghost cells
    // cfg.hx = cfg.lengthX / (cfg.Nx - 2); // Grid spacing in the x-direction
    // cfg.hy = cfg.lengthY / (cfg.Ny - 2); // Grid spacing in the y-direction
    cfg.poissonTolerance = 1e-3; // Tolerance for the Poisson equation solver
    cfg.UtimeTolerance = 1e-08;  // Tolerance for the time loop convergence
    cfg.VtimeTolerance = 1e-08;  // Tolerance for the time loop convergence
    cfg.PtimeTolerance = 1e-08;  // Tolerance for the time loop convergence
    cfg.startTime = 0.0;         // Start time of the simulation
    cfg.endTime = 10000.0;       // End time of the simulation. Since we have a convergence criteria, this is set to a large number

    cfg.firstCellHeightX = 1.0 / (cfg.Nx - 2); // y+ value for the first cell height
    cfg.numberOfMeshLayersX = 1;               // Number of mesh layers for the non-uniform mesh
    cfg.growthRateX = 1.0;                     // Growth rate for the non-uniform mesh

    cfg.firstCellHeightY = 1.0 / (cfg.Ny - 2); // y+ value for the first cell height
    cfg.numberofMeshLayersY = 1;               // Number of mesh layers for the non-uniform mesh
    cfg.growthRateY = 1.0;                     // Growth rate for the non-uniform mesh

    nonUniformMeshX(cfg); // Create the non-uniform mesh
    nonUniformMeshY(cfg); // Create the non-uniform mesh
    fields f(cfg.Nx, cfg.Ny);

    initialization(f, cfg); // Initialize the fields

    createCoordinates(f, cfg); // Create the coordinates of the nodes of cells
    createAreaU(f, cfg);
    createAreaP(f, cfg);
    createAreaV(f, cfg);

    createVolumeU(f, cfg); // Create the volume of the cells
    createVolumeV(f, cfg); // Create the volume of the cells

    cfg.timeStepSize = min((pow(cfg.firstCellHeightX, 2) * pow(cfg.firstCellHeightY, 2)) / (2 * cfg.dynamicViscosity * (pow(cfg.firstCellHeightX, 2) + pow(cfg.firstCellHeightY, 2))), 2 * cfg.dynamicViscosity / (pow(cfg.uTopWall, 2) + pow(cfg.vLeftWall, 2)));

    setBoundaryConditions(1, f.u, cfg);
    setBoundaryConditions(2, f.v, cfg);
    setBoundaryConditions(0, f.p, cfg);
    // for (int i = 0; i < cfg.Nx; i++)
    // {
    //     cout << f.vx[i] << endl;
    // }
    cout << "Inflation layer thickness in x is " << cfg.inflationLayerThicknessX << endl;
    cout << "Inflation layer thickness in y is " << cfg.inflationLayerThicknessY << endl;
    cout << "Number of cells in the x-direction is " << cfg.Nx << endl;
    cout << "Number of cells in the y-direction is " << cfg.Ny << endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    int n = 0; // counter for the time steps
    for (double t = cfg.startTime; t <= cfg.endTime; t = t + cfg.timeStepSize)
    {

        vector<vector<double>> pPrev = f.p;
        vector<vector<double>> uPrev = f.u;
        vector<vector<double>> vPrev = f.v;

        predictor(f, cfg);
        pressurePoisson(f, cfg);
        corrector(f, cfg);
        swapFields(f, cfg);
        setBoundaryConditions(1, f.u, cfg);
        setBoundaryConditions(2, f.v, cfg);
        setBoundaryConditions(0, f.p, cfg);
        updateTimeStepSize(f, cfg);
        double residualU = checkConvergence(f.u, uPrev, cfg);
        double residualV = checkConvergence(f.v, vPrev, cfg);
        double residualP = checkConvergence(f.p, pPrev, cfg);
        // std::this_thread::sleep_for(std::chrono::milliseconds(200));
        cout << (f.u[(cfg.Ny) / 2][1 + (cfg.Nx) / 2]) << endl;
        // cout << "Time step size: " << cfg.timeStepSize << endl;
        cout << "Time: " << t << " U velocity: " << residualU << " V velocity: " << residualV << " pressure: " << residualP << " n: " << n << endl;
        if (residualU < cfg.UtimeTolerance && residualV < cfg.VtimeTolerance && residualP < cfg.PtimeTolerance)
        {
            cout << "Converged at time: " << t << endl;
            break;
        }
        n++;
    }

    // for (int j = 1; j < cfg.Ny - 1; j++)
    // {
    //     cout << f.py[j] << " " << (f.p[j][(cfg.Nx - 1) / 2]) << endl;
    // }
    // cout << endl;

    cout << "Inflation layer thickness in x is " << cfg.inflationLayerThicknessX << endl;
    cout << "Inflation layer thickness in y is " << cfg.inflationLayerThicknessY << endl;
    cout << "Number of cells in the x-direction is " << cfg.Nx << endl;
    cout << "Number of cells in the y-direction is " << cfg.Ny << endl;
    cout << "Time step size is " << cfg.timeStepSize << endl;

    std::cout << "\nEnd of the main function is reached. Stopping.\n\n";

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;

    cell_writeVTKFile(f, cfg);

    for (int i = 0; i < cfg.Ny; i++)
    {

        cout << f.py[i] << " " << f.p[i][(cfg.Nx - 1) / 2] << endl;
    }

    return 0;
}
