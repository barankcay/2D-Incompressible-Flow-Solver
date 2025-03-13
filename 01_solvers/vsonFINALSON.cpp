#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting
#include <cmath>
using namespace std;

struct constParameters
{
    int Nx;
    int Ny;
    double lengthX;
    double lengthY;
    double hx;
    double hy;
    double startTime;
    double endTime;
    double timeStepSize;
    double kinematicViscosity;
    double density;

    double uTopWall;
    double uBottomWall;
    double uLeftWall;
    double uRightWall;
    double vTopWall;
    double vBottomWall;
    double vLeftWall;
    double vRightWall;
    double courantNumber;
    double poissonTolerance;
    double timeTolerance;
    int maxIterations;
    int numberOfTimeSteps;
};

struct fields
{
    vector<double> x;
    vector<double> y;
    vector<double> xm;
    vector<double> ym;

    vector<vector<double>> u;
    vector<vector<double>> uStar;
    vector<vector<double>> uNew;
    vector<vector<double>> v;
    vector<vector<double>> vStar;
    vector<vector<double>> vNew;
    vector<vector<double>> p;

    fields(int Nx, int Ny)
    {
        x.resize(Nx - 1);
        y.resize(Ny - 1);
        xm.resize(Nx - 2);
        ym.resize(Ny - 2);
        u.resize(Ny, vector<double>(Nx));
        uStar.resize(Ny, vector<double>(Nx));
        uNew.resize(Ny, vector<double>(Nx));
        vStar.resize(Ny, vector<double>(Nx));
        vNew.resize(Ny, vector<double>(Nx));
        v.resize(Ny, vector<double>(Nx));
        p.resize(Ny, vector<double>(Nx));
    }
};

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
            M[j][0] = M[j][1];
            M[j][params.Nx - 1] = M[j][params.Nx - 2];
        }
        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = M[1][i];
            M[params.Ny - 1][i] = M[params.Ny - 2][i];
        }
    }
    else if (b == 1)
    {
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 0;
            M[j][1] = 0;
            M[j][params.Nx - 1] = 0;
        }
        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * params.uTopWall - M[1][i];
            M[params.Ny - 1][i] = 2 * params.uBottomWall - M[params.Ny - 2][i];
        }
    }
    else if (b == 2)
    {

        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 0;
            M[params.Ny - 1][i] = 0;
            M[params.Ny - 2][i] = 0;
        }
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * params.vLeftWall - M[j][1];
            M[j][params.Nx - 1] = 2 * params.vRightWall - M[j][params.Nx - 2];
        }
    }

    M[0][0] = 0.5 * (M[0][1] + M[1][0]);
    M[0][params.Nx - 1] = 0.5 * (M[0][params.Nx - 2] + M[1][params.Nx - 1]);
    M[params.Ny - 1][0] = 0.5 * (M[params.Ny - 2][0] + M[params.Ny - 1][1]);
    M[params.Ny - 1][params.Nx - 1] = 0.5 * (M[params.Ny - 1][params.Nx - 2] + M[params.Ny - 2][params.Nx - 1]);
}
// function to calculate the first order partial derivative using central differencing
double firstOrderPDEcentralDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters &params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j - y][i + x] - variable[j + y][i - x]) / (2 * (x * params.hx + y * params.hy));
    return pde;
}

// function to calculate the second order partial derivative using central differencing
double secondOrderPDEcentralDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters &params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j - y][i + x] - 2 * variable[j][i] + variable[j + y][i - x]) / pow((x * params.hx + y * params.hy), 2);
    return pde;
}

double firstOrderPDEforwardDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters &params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j - y][i + x] - variable[j][i]) / (x * params.hx + y * params.hy);
    return pde;
}

double firstOrderPDEbackwardDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters &params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j][i] - variable[j + y][i - x]) / (x * params.hx + y * params.hy);
    return pde;
}

// PREDICTOR STEP
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
                field.p[j][i] = (pow(params.hy, 2) * pow(params.hx, 2) / (2 * (pow(params.hx, 2) + pow(params.hy, 2)))) * (((field.p[j][i + 1] + field.p[j][i - 1]) / pow(params.hx, 2)) + ((field.p[j - 1][i] + field.p[j + 1][i]) / pow(params.hy, 2)) - ((params.density / params.timeStepSize) * ((field.uStar[j][i + 1] - field.uStar[j][i]) / (params.hx) + (field.vStar[j - 11][i] - field.vStar[j][i]) / (params.hy))));
                residual += abs(field.p[j][i] - pOld);
            }
        }
        setBoundaryConditions(0, field.p, params);
        residual /= (params.Nx * params.Ny);
        iteration++;
    }

    cout << "Number of iterations: " << iteration << endl;
}

void updateTimeStepSize(fields &field, constParameters &params)
{
    double maxVelocity = 0.0;

    // Calculate the maximum velocity in the domain (including boundaries)
    for (int i = 0; i < params.Nx; i++)
    {
        for (int j = 0; j < params.Ny; j++)
        {
            double velocity = sqrt(pow(field.u[j][i], 2) + pow(field.v[j][i], 2));
            if (velocity > maxVelocity)
            {
                maxVelocity = velocity;
            }
        }
    }

    // Handle the case where maxVelocity is zero (initial time step)
    if (maxVelocity == 0.0)
    {
        // Set a default time step size based on the grid spacing and Courant number
        params.timeStepSize = params.courantNumber * min(params.hx, params.hy);
    }
    else
    {
        // Calculate the time step size based on the maximum velocity
        params.timeStepSize = params.courantNumber * min(params.hx, params.hy) / maxVelocity;
    }
}

// CORRECTOR STEP
void velocityCorrector(fields &field, constParameters &params)
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
    double Re = 400;

    constParameters params;
    params.courantNumber = 0.2;
    params.density = 1.0;
    params.kinematicViscosity = 1.0 / Re;

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

    params.timeStepSize = params.courantNumber * min(params.hx, params.hy) / params.uTopWall;
    // params.timeStepSize = 0.001;
    params.numberOfTimeSteps = (params.endTime - params.startTime) / params.timeStepSize;
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
