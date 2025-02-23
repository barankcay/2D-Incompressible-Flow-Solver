#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

struct constParameters
{
    int Nx;
    int Ny;
    double lengthX;
    double lengthY;
    double hx = lengthX / (Nx + 2);
    double hy = lengthY / (Ny + 2);
    double time;
    double timeStepSize;
    double kinematicViscosity;
    double density;
    double uTopWall;
    double uBottomWall;
    double uLeftWall;
    double uRightWall;
};

struct fields
{
    vector<vector<double>> x;
    vector<vector<double>> y;
    vector<vector<double>> xm;
    vector<vector<double>> ym;
    vector<vector<double>> u;
    vector<vector<double>> uStar;
    vector<vector<double>> uNew;
    vector<vector<double>> v;
    vector<vector<double>> vStar;
    vector<vector<double>> vNew;
    vector<vector<double>> p;

    fields(int Nx, int Ny)
    {
        x.resize(Nx + 3, vector<double>(Ny + 3));
        y.resize(Nx + 3, vector<double>(Ny + 3));
        xm.resize(Nx + 2, vector<double>(Ny + 2));
        ym.resize(Nx + 2, vector<double>(Ny + 2));
        u.resize(Nx + 2, vector<double>(Ny + 2));
        uStar.resize(Nx + 2, vector<double>(Ny + 2));
        vStar.resize(Nx + 2, vector<double>(Ny + 2));
        v.resize(Nx + 2, vector<double>(Ny + 2));
        p.resize(Nx + 2, vector<double>(Ny + 2));
    }
};

// function to create the coordinates of the nodes of cells
void createCoordinatesXY(vector<vector<double>> &x, vector<vector<double>> &y, constParameters params)
{
    for (int i = 0; i < params.Nx + 3; i++)
    {
        for (int j = 0; j < params.Ny + 3; j++)
        {
            x[i][j] = i * params.hx;
            y[i][j] = j * params.hy;
        }
    }
}

// function to create the coordinates of the cell centers
void createCoordinatesXYM(vector<vector<double>> &xm, vector<vector<double>> &ym, constParameters params)
{
    for (int i = 0; i < params.Nx + 2; i++)
    {
        for (int j = 0; j < params.Ny + 2; j++)
        {
            xm[i][j] = params.hx / 2 + i * params.hx;
            ym[i][j] = params.hy / 2 + j * params.hy;
        }
    }
}

// function to calculate the first order partial derivative using central differencing
double firstOrderPDEcentralDiff(vector<vector<double>> &variable, int i, int j, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[i + x][j + y] - variable[i - x][j - y]) / (2 * (x * params.hx + y * params.hy));
    return pde;
}

// function to calculate the second order partial derivative using central differencing
double secondOrderPDEcentralDiff(vector<vector<double>> &variable, int i, int j, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[i + x][j + y] - 2 * variable[i][j] + variable[i - x][j - y]) / pow((x * params.hx + y * params.hy), 2);
    return pde;
}

double firstOrderPDEforwardDiff(vector<vector<double>> &variable, int i, int j, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[i + x][j + y] - variable[i][j]) / (x * params.hx + y * params.hy);
    return pde;
}

double firstOrderPDEbackwardDiff(vector<vector<double>> &variable, int i, int j, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[i][j] - variable[i - x][j - y]) / (x * params.hx + y * params.hy);
    return pde;
}

// PREDICTOR STEP
void veloctiyStarCalculator(fields &field, constParameters params)
{
    for (int i = 1; i < params.Nx + 1; i++)
    {
        for (int j = 1; j < params.Ny + 1; j++)
        {
            field.uStar[i][j] = field.u[i][j] + params.timeStepSize * (params.kinematicViscosity * (secondOrderPDEcentralDiff(field.u, i, j, 1, 0, params) + secondOrderPDEcentralDiff(field.u, i, j, 0, 1, params)) - (field.u[i][j] * firstOrderPDEcentralDiff(field.u, i, j, 1, 0, params) + 0.25 * (field.v[i + 1][j] + field.v[i - 1][j] + field.v[i][j + 1] + field.v[i][j - 1]) * firstOrderPDEcentralDiff(field.u, i, j, 0, 1, params)));
            field.vStar[i][j] = field.v[i][j] + params.timeStepSize * (params.kinematicViscosity * (secondOrderPDEcentralDiff(field.v, i, j, 1, 0, params) + secondOrderPDEcentralDiff(field.v, i, j, 0, 1, params)) - (field.v[i][j] * firstOrderPDEcentralDiff(field.v, i, j, 0, 1, params) + 0.25 * (field.u[i + 1][j] + field.u[i - 1][j] + field.u[i][j + 1] + field.u[i][j - 1]) * firstOrderPDEcentralDiff(field.v, i, j, 1, 0, params)));
        }
    }
}
// POISSON EQUATION SOLVER

void poissonEquationSolver(fields &field, constParameters params)
{
    for (int k = 0; k < 300; k++)
    {
        for (int i = 1; i < params.Nx + 1; i++)
        {
            for (int j = 1; j < params.Ny + 1; j++)
            {
                field.p[i][j] = (pow(params.hx, 2) * pow(params.hx, 2) / (2 * (params.hx + params.hy))) * (((field.p[i + 1][j] + field.p[i - 1][j]) / pow(params.hx, 2)) + ((field.p[i][j + 1] + field.p[i][j - 1]) / pow(params.hy, 2)) - ((params.density / params.timeStepSize) * (firstOrderPDEforwardDiff(field.uStar, i, j, 1, 0, params) + firstOrderPDEforwardDiff(field.vStar, i, j, 0, 1, params))));
            }
        }
    }
}

// CORRECTOR STEP
void velocityCorrector(fields &field, constParameters params)
{
    for (int i = 1; i < params.Nx + 1; i++)
    {
        for (int j = 1; j < params.Ny + 1; j++)
        {
            field.uNew[i][j] = field.uStar[i][j] - (params.timeStepSize / params.density) * (firstOrderPDEbackwardDiff(field.p, i, j, 1, 0, params));
            field.vNew[i][j] = field.vStar[i][j] - (params.timeStepSize / params.density) * (firstOrderPDEbackwardDiff(field.p, i, j, 0, 1, params));
        }
    }
}

void swapFields(fields &field)
{
    field.u = field.uNew;
    field.v = field.vNew;
}

void setBoundaryConditions(fields &field, constParameters params)
{
}
int main()
{

    constParameters params;
    params.density = 1.0;
    params.kinematicViscosity = 0.01;

    params.uTopWall = 1.0;
    params.uBottomWall = 0.0;
    params.uLeftWall = 0.0;
    params.uRightWall = 0.0;

    params.Nx = 3;
    params.Ny = 3;
    params.lengthX = 1.0;
    params.lengthY = 1.0;
    params.hx = params.lengthX / (params.Nx + 2);
    params.hy = params.lengthY / (params.Ny + 2);

    params.time = 20;
    params.timeStepSize = 0.1;

    fields field(params.Nx, params.Ny);
    createCoordinatesXY(field.x, field.y, params);
    createCoordinatesXYM(field.xm, field.ym, params);

    // for (int i = 0; i < params.Nx + 3; i++)
    // {
    //     for (int j = 0; j < params.Ny + 3; j++)
    //     {
    //         cout << field.x[i][j] << "," << field.y[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    return 0;
}