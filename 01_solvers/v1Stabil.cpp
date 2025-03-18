#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting
#include <cmath>
using namespace std;

struct constParameters
{
    int Nx;                    // Number of cells in the x-direction, including ghost cells
    int Ny;                    // Number of cells in the y-direction, including ghost cells
    double lengthX;            // Length of the domain in the x-direction, not including ghost cells
    double lengthY;            // Length of the domain in the y-direction, not including ghost cells
    double hx;                 // Grid spacing in the x-direction
    double hy;                 // Grid spacing in the y-direction
    double startTime;          // Start time of the simulation
    double endTime;            // End time of the simulation. Since we have a convergence criteria, this is set to a large number
    double timeStepSize;       // Time step size. Calculated based on the Courant number and maximum velocity in the domain
    double kinematicViscosity; // Kinematic viscosity of the fluid. Used to calculate the Reynolds number
    double dynamicViscosity;   // Dynamic viscosity of the fluid. Dynamic viscosity = density * kinematic viscosity
    double density;            // Density of the fluid.

    double uTopWall;         // u Velocity at the top wall
    double uBottomWall;      // u Velocity at the bottom wall
    double vLeftWall;        // v Velocity at the left wall
    double vRightWall;       // v Velocity at the right wall
    double courantNumber;    // Courant number used to calculate the time step size
    double poissonTolerance; // Tolerance for the Poisson equation solver
    double timeTolerance;    // Tolerance for the time loop convergence
};

struct fields
{
    vector<double> x;  // x-coordinates of the nodes of the cells
    vector<double> y;  // y-coordinates of the nodes of the cells
    vector<double> xm; // x-coordinates of the cell centers
    vector<double> ym; // y-coordinates of the cell centers

    vector<vector<double>> u;     // u velocity field
    vector<vector<double>> uStar; // intermediate u velocity field
    vector<vector<double>> uNew;  // next time step u velocity field
    vector<vector<double>> v;     // v velocity field
    vector<vector<double>> vStar; // intermediate v velocity field
    vector<vector<double>> vNew;  // next time step v velocity field
    vector<vector<double>> p;     // pressure field

    fields(int Nx, int Ny)
    {
        x.resize(Nx - 1);  // not including ghost cells
        y.resize(Ny - 1);  // not including ghost cells
        xm.resize(Nx - 2); // not including ghost cells
        ym.resize(Ny - 2); // not including ghost cells
        u.resize(Ny, vector<double>(Nx));
        uStar.resize(Ny, vector<double>(Nx));
        uNew.resize(Ny, vector<double>(Nx));
        vStar.resize(Ny, vector<double>(Nx));
        vNew.resize(Ny, vector<double>(Nx));
        v.resize(Ny, vector<double>(Nx));
        p.resize(Ny, vector<double>(Nx));
    }
};

// function to initialize the fields
void initialization(fields &field, constParameters &params)
{
    for (int i = 0; i < params.Nx; i++)
    {
        for (int j = 0; j < params.Ny; j++)
        {
            field.u[j][i] = 0;
            field.v[j][i] = 0;
            field.p[j][i] = 0;
            field.uStar[j][i] = 0;
            field.vStar[j][i] = 0;
            field.uNew[j][i] = 0;
            field.vNew[j][i] = 0;
        }
    }
}

// function to create the coordinates of the nodes of cells
void createCoordinatesXY(vector<double> &x, vector<double> &y, constParameters &params)
{
    for (int j = 0; j < params.Ny - 1; j++)
    {
        for (int i = 0; i < params.Nx - 1; i++)
        {
            x[i] = i * params.hx;
            y[j] = (params.Ny - j) * params.hy;
        }
    }
}

// function to create the coordinates of the cell centers
void createCoordinatesXYM(vector<double> &xm, vector<double> &ym, constParameters &params)
{
    for (int j = 0; j < params.Ny - 2; j++)
    {
        for (int i = 0; i < params.Nx - 2; i++)
        {
            xm[i] = params.hx / 2 + i * params.hx;

            ym[j] = params.hy / 2 + (params.Ny - 3 - j) * params.hy;
        }
    }
}

void setBoundaryConditions(int b, vector<vector<double>> &M, constParameters &params)
{
    if (b == 0)
    {
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = M[j][1];                         // Left wall, dp/dx = 0
            M[j][params.Nx - 1] = M[j][params.Nx - 2]; // Right wall, dp/dx = 0
        }
        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = M[1][i];                         // Bottom wall, dp/dy = 0
            M[params.Ny - 1][i] = M[params.Ny - 2][i]; // Top wall, dp/dy = 0
        }
    }
    else if (b == 1)
    {
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 0;             // Left wall,ghost cell,  u = 0
            M[j][1] = 0;             // Left wall, u = 0
            M[j][params.Nx - 1] = 0; // Right wall, u = 0
        }
        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * params.uTopWall - M[1][i];                            // Top wall, ghost cell
            M[params.Ny - 1][i] = 2 * params.uBottomWall - M[params.Ny - 2][i]; // Bottom wall, ghost cell
        }
    }
    else if (b == 2)
    {

        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 0;             // Top wall, v = 0
            M[params.Ny - 1][i] = 0; // Bottom wall, v = 0
            M[params.Ny - 2][i] = 0; // Bottom wall, ghost cell, v = 0
        }
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * params.vLeftWall - M[j][1];                          // Left wall, ghost cell
            M[j][params.Nx - 1] = 2 * params.vRightWall - M[j][params.Nx - 2]; // Right wall, ghost cell
        }
    }

    M[0][0] = 0.5 * (M[0][1] + M[1][0]);                                                                         // Top left corner
    M[0][params.Nx - 1] = 0.5 * (M[0][params.Nx - 2] + M[1][params.Nx - 1]);                                     // Top right corner
    M[params.Ny - 1][0] = 0.5 * (M[params.Ny - 2][0] + M[params.Ny - 1][1]);                                     // Bottom left corner
    M[params.Ny - 1][params.Nx - 1] = 0.5 * (M[params.Ny - 1][params.Nx - 2] + M[params.Ny - 2][params.Nx - 1]); // Bottom right corner
}
// PREDICTOR STEP
// u* = un + Δt * ( ν ( ∂²un/∂x² + ∂²un/∂y² ) - ( un ∂un/∂x + v ∂un/∂y ) )
// v* = vn + Δt * ( ν ( ∂²vn/∂x² + ∂²vn/∂y² ) - ( u ∂vn/∂x + vn ∂vn/∂y ) )
void veloctiyStarCalculator(fields &field, constParameters &params)
{
    for (int i = 1; i < params.Nx - 1; i++)
    {
        for (int j = 1; j < params.Ny - 1; j++)
        {
            field.uStar[j][i] = field.u[j][i] + params.timeStepSize * (params.kinematicViscosity * ((field.u[j][i + 1] - 2 * field.u[j][i] + field.u[j][i - 1]) / pow(params.hx, 2) + (field.u[j + 1][i] - 2 * field.u[j][i] + field.u[j - 1][i]) / pow(params.hy, 2)) - (field.u[j][i] * (field.u[j][i + 1] - field.u[j][i - 1]) / (2 * params.hx) + 0.25 * (field.v[j][i - 1] + field.v[j][i] + field.v[j - 1][i - 1] + field.v[j - 1][i]) * (field.u[j - 1][i] - field.u[j + 1][i]) / (2 * params.hy)));
            field.vStar[j][i] = field.v[j][i] + params.timeStepSize * (params.kinematicViscosity * ((field.v[j][i + 1] - 2 * field.v[j][i] + field.v[j][i - 1]) / pow(params.hx, 2) + (field.v[j + 1][i] - 2 * field.v[j][i] + field.v[j - 1][i]) / pow(params.hy, 2)) - (field.v[j][i] * (field.v[j - 1][i] - field.v[j + 1][i]) / (2 * params.hy) + 0.25 * (field.u[j + 1][i] + field.u[j][i] + field.u[j + 1][i + 1] + field.u[j][i + 1]) * (field.v[j][i + 1] - field.v[j][i - 1]) / (2 * params.hx)));
        }
    }
    setBoundaryConditions(1, field.uStar, params);
    setBoundaryConditions(2, field.vStar, params);
}
// POISSON EQUATION SOLVER
// ∇² p(n+1) = (ρ / Δt) * ∇ · u*
void poissonEquationSolver(fields &field, constParameters &params)
{
    double residual = 1.0;
    int iteration = 0;
    while (residual > params.poissonTolerance)
    {
        residual = 0;
        for (int i = 1; i < params.Nx - 1; i++)
        {
            for (int j = 1; j < params.Ny - 1; j++)
            {
                double pOld = field.p[j][i];
                field.p[j][i] = (pow(params.hy, 2) * pow(params.hx, 2) / (2 * (pow(params.hx, 2) + pow(params.hy, 2)))) * (((field.p[j][i + 1] + field.p[j][i - 1]) / pow(params.hx, 2)) + ((field.p[j - 1][i] + field.p[j + 1][i]) / pow(params.hy, 2)) - ((params.density / params.timeStepSize) * ((field.uStar[j][i + 1] - field.uStar[j][i]) / (params.hx) + (field.vStar[j - 1][i] - field.vStar[j][i]) / (params.hy))));
                residual += abs(field.p[j][i] - pOld);
            }
        }
        setBoundaryConditions(0, field.p, params);
        residual /= (params.Nx * params.Ny);
        iteration++;
    }

    cout << "Number of iterations: " << iteration << endl;
}

// TIME STEP SIZE CALCULATION
// function to calculate the time step size based on the Courant number and maximum velocity in the domain
void updateTimeStepSize(fields &field, constParameters &params)
{
    // Calculate the maximum velocity in the domain (including boundaries)
    for (int i = 1; i < params.Nx - 1; i++)
    {
        for (int j = 1; j < params.Ny - 1; j++)
        {
            double rule1 = (pow(params.hx, 2) * pow(params.hy, 2)) / (2 * params.dynamicViscosity * (pow(params.hx, 2) + pow(params.hy, 2)));
            double rule2 = 2 * params.dynamicViscosity / (pow(field.u[j][i], 2) + pow(field.v[j][i], 2));
            double finalRule = min(rule1, rule2);
            if (finalRule < params.timeStepSize)
            {
                params.timeStepSize = finalRule;
            }
        }
    }
}

// CORRECTOR STEP
void velocityCorrector(fields &field, constParameters &params)
// u(n+1)   =u*+ Δt*(- (1 / ρ) * ∇p(n+1))
{
    for (int i = 1; i < params.Nx - 1; i++)
    {
        for (int j = 1; j < params.Ny - 1; j++)
        {
            field.uNew[j][i] = field.uStar[j][i] - (params.timeStepSize / params.density) * ((field.p[j][i] - field.p[j][i - 1]) / (params.hx));
            field.vNew[j][i] = field.vStar[j][i] - (params.timeStepSize / params.density) * ((field.p[j][i] - field.p[j + 1][i]) / (params.hy));
        }
    }
}

// CONVERGENCE CHECK
// function to check the convergence of the solution
double checkConvergence(vector<vector<double>> &Mnew, vector<vector<double>> &Mold, constParameters &params)
{
    double residual = 0.0;
    for (int i = 1; i < params.Nx - 1; i++)
    {
        for (int j = 1; j < params.Ny - 1; j++)
        {
            residual += abs(Mnew[j][i] - Mold[j][i]);
        }
    }
    return residual /= (params.Nx * params.Ny);
}

// function to swap the fields to prepare for the next time step
void swapFields(fields &field, constParameters &params)
{
    for (int i = 0; i < params.Nx; i++)
    {
        for (int j = 0; j < params.Ny; j++)
        {
            field.u[j][i] = field.uNew[j][i];
            field.v[j][i] = field.vNew[j][i];
            field.uStar[j][i] = field.u[j][i];
            field.vStar[j][i] = field.v[j][i];
        }
    }
}

int main()
{
    double Re = 1000;

    constParameters params;
    params.courantNumber = 0.2;
    params.density = 1.0;
    params.kinematicViscosity = 1.0 / Re;
    params.dynamicViscosity = params.kinematicViscosity * params.density;

    params.uTopWall = 1;
    params.uBottomWall = 0.0;

    params.vLeftWall = 0.0;
    params.vRightWall = 0.0;

    params.Nx = 168;
    params.Ny = 168;
    params.lengthX = 1;
    params.lengthY = 1;
    params.hx = params.lengthX / (params.Nx - 2);
    params.hy = params.lengthY / (params.Ny - 2);

    // params.maxIterations = 1000;
    params.poissonTolerance = 1e-3;
    params.timeTolerance = 1e-6;

    params.startTime = 0;
    params.endTime = 10000;

    fields field(params.Nx, params.Ny);

    params.timeStepSize = min((pow(params.hx, 2) * pow(params.hy, 2)) / (2 * params.dynamicViscosity * (pow(params.hx, 2) + pow(params.hy, 2))), 2 * params.dynamicViscosity / (pow(params.uTopWall, 2) + pow(params.vLeftWall, 2)));
    // params.timeStepSize = 0.001;

    initialization(field, params);
    createCoordinatesXY(field.x, field.y, params);
    createCoordinatesXYM(field.xm, field.ym, params);

    setBoundaryConditions(1, field.u, params);
    setBoundaryConditions(2, field.v, params);
    setBoundaryConditions(0, field.p, params);
    int n = 0; // counter for the time steps
    for (double t = params.startTime; t <= params.endTime; t = t + params.timeStepSize)
    {

        vector<vector<double>> pPrev = field.p;
        vector<vector<double>> uPrev = field.u;
        vector<vector<double>> vPrev = field.v;
        updateTimeStepSize(field, params);
        veloctiyStarCalculator(field, params);
        poissonEquationSolver(field, params);
        velocityCorrector(field, params);
        swapFields(field, params);
        setBoundaryConditions(1, field.u, params);
        setBoundaryConditions(2, field.v, params);
        setBoundaryConditions(0, field.p, params);
        double residualU = checkConvergence(field.u, uPrev, params);
        double residualV = checkConvergence(field.v, vPrev, params);
        double residualP = checkConvergence(field.p, pPrev, params);

        cout << "Time: " << t << " U velocity: " << residualU << " V velocity: " << residualV << " pressure: " << residualP << " n: " << n << endl;
        if (residualU < params.timeTolerance && residualV < params.timeTolerance && residualP < params.timeTolerance)
        {
            cout << "Converged at time: " << t << endl;
            break;
        }
        n++;
    }
    for (int j = 0; j < params.Ny - 2; j++)
    {

        cout << field.ym[j] << " " << (field.u[j + 1][(params.Nx) / 2]) << endl;
    }
    cout << endl;

    cout << params.timeStepSize << endl;
    return 0;
}
