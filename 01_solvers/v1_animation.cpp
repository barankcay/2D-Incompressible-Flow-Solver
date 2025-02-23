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
    double vTopWall;
    double vBottomWall;
    double uLeftWall;
    double uRightWall;
    double courantNumber;
    int numberOfTimeSteps;
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
        y.resize(Ny + 3, vector<double>(Ny + 3));
        xm.resize(Nx + 2, vector<double>(Ny + 2));
        ym.resize(Ny + 2, vector<double>(Ny + 2));
        u.resize(Nx + 2, vector<double>(Ny + 2));
        uStar.resize(Nx + 2, vector<double>(Ny + 2));
        uNew.resize(Nx + 2, vector<double>(Ny + 2));
        vStar.resize(Nx + 2, vector<double>(Ny + 2));
        vNew.resize(Nx + 2, vector<double>(Ny + 2));
        v.resize(Nx + 2, vector<double>(Ny + 2));
        p.resize(Nx + 2, vector<double>(Ny + 2));
    }
};

// function to create the coordinates of the nodes of cells
void createCoordinatesXY(vector<vector<double>> &x, vector<vector<double>> &y, constParameters params)
{
    for (int i = 0; i < params.Nx + 2; i++)
    {
        for (int j = 0; j < params.Ny + 2; j++)
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
            field.uStar[i][j] = field.u[i][j] + params.timeStepSize * (params.kinematicViscosity * (secondOrderPDEcentralDiff(field.u, i, j, 1, 0, params) + secondOrderPDEcentralDiff(field.u, i, j, 0, 1, params)) - (field.u[i][j] * firstOrderPDEcentralDiff(field.u, i, j, 1, 0, params) + 0.25 * (field.v[i - 1][j] + field.v[i][j] + field.v[i - 1][j + 1] + field.v[i][j + 1]) * firstOrderPDEcentralDiff(field.u, i, j, 0, 1, params)));
            field.vStar[i][j] = field.v[i][j] + params.timeStepSize * (params.kinematicViscosity * (secondOrderPDEcentralDiff(field.v, i, j, 1, 0, params) + secondOrderPDEcentralDiff(field.v, i, j, 0, 1, params)) - (field.v[i][j] * firstOrderPDEcentralDiff(field.v, i, j, 0, 1, params) + 0.25 * (field.u[i][j - 1] + field.u[i][j] + field.u[i + 1][j - 1] + field.u[i + 1][j]) * firstOrderPDEcentralDiff(field.v, i, j, 1, 0, params)));
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
                field.p[i][j] = (pow(params.hx, 2) * pow(params.hx, 2) / (2 * (pow(params.hx, 2) + pow(params.hy, 2)))) * (((field.p[i + 1][j] + field.p[i - 1][j]) / pow(params.hx, 2)) + ((field.p[i][j + 1] + field.p[i][j - 1]) / pow(params.hy, 2)) - ((params.density / params.timeStepSize) * (firstOrderPDEforwardDiff(field.uStar, i, j, 1, 0, params) + firstOrderPDEforwardDiff(field.vStar, i, j, 0, 1, params))));
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

void swapFields(fields &field, constParameters params)
{
    for (int i = 1; i < params.Nx + 1; i++)
    {
        for (int j = 1; j < params.Ny + 1; j++)
        {
            field.u[i][j] = field.uNew[i][j];
            field.v[i][j] = field.vNew[i][j];
        }
    }
}

void setBoundaryConditions(int b, vector<vector<double>> &M, constParameters params)
{
    for (int i = 1; i < params.Nx + 1; i++)
    {
        M[i][0] = b == 1 ? 2 * params.uLeftWall - M[i][1] : M[i][1];
        M[i][params.Ny + 1] = b == 1 ? 2 * params.uRightWall - M[i][params.Ny] : M[i][params.Ny];
    }
    for (int j = 1; j < params.Ny + 1; j++)
    {
        M[0][j] = b == 2 ? 2 * params.vTopWall - M[1][j] : M[1][j];
        M[params.Nx + 1][j] = b == 2 ? 2 * params.vBottomWall - M[params.Nx][j] : M[params.Nx][j];
    }
    M[0][0] = 0.5 * (M[1][0] + M[0][1]);
    M[params.Nx + 1][0] = 0.5 * (M[params.Nx][0] + M[params.Nx + 1][1]);
    M[0][params.Ny + 1] = 0.5 * (M[0][params.Ny] + M[1][params.Ny + 1]);
    M[params.Nx + 1][params.Ny + 1] = 0.5 * (M[params.Nx][params.Ny + 1] + M[params.Nx + 1][params.Ny]);
}

void saveToVTK(const fields &field, const constParameters &params, const string &filename)
{
    ofstream vtkFile(filename);
    if (!vtkFile.is_open())
    {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    // Write the VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Flow Field Data\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << params.Nx + 2 << " " << params.Ny + 2 << " 1\n";
    vtkFile << "POINTS " << (params.Nx + 2) * (params.Ny + 2) << " float\n";

    // Write the grid coordinates
    for (int j = 0; j < params.Ny + 2; j++)
    {
        for (int i = 0; i < params.Nx + 2; i++)
        {
            vtkFile << field.x[i][j] << " " << field.y[i][j] << " 0.0\n";
        }
    }

    // Write the velocity and pressure data
    vtkFile << "POINT_DATA " << (params.Nx + 2) * (params.Ny + 2) << "\n";
    vtkFile << "VECTORS velocity float\n";
    for (int j = 0; j < params.Ny + 2; j++)
    {
        for (int i = 0; i < params.Nx + 2; i++)
        {
            vtkFile << field.u[i][j] << " " << field.v[i][j] << " 0.0\n";
        }
    }

    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < params.Ny + 2; j++)
    {
        for (int i = 0; i < params.Nx + 2; i++)
        {
            vtkFile << field.p[i][j] << "\n";
        }
    }

    vtkFile.close();
    cout << "Data saved to " << filename << endl;
}
int main()
{

    constParameters params;
    params.courantNumber = 0.01;
    params.density = 1.0;
    params.kinematicViscosity = 0.01;

    params.vTopWall = 0.5;
    params.vBottomWall = 0.0;
    params.uLeftWall = 0.0;
    params.uRightWall = 0.0;

    params.Nx = 18;
    params.Ny = 18;
    params.lengthX = 0.1;
    params.lengthY = 0.1;
    params.hx = params.lengthX / (params.Nx + 2);
    params.hy = params.lengthY / (params.Ny + 2);

    params.startTime = 0;
    params.endTime = 0.5;
    params.timeStepSize = params.courantNumber * min(params.hx, params.hy) / params.vTopWall;
    params.numberOfTimeSteps = (params.endTime - params.startTime) / params.timeStepSize;

    fields field(params.Nx, params.Ny);
    createCoordinatesXY(field.x, field.y, params);
    createCoordinatesXYM(field.xm, field.ym, params);
    setBoundaryConditions(1, field.u, params);
    setBoundaryConditions(2, field.v, params);
    setBoundaryConditions(0, field.p, params);
    int count = 0;
    for (double t = params.startTime; t < params.endTime; t = t + params.timeStepSize)
    {
        if (count == params.numberOfTimeSteps - 1) // Save every 10 time steps
        {
            string filename = "output_" + to_string(count) + ".vtk";
            saveToVTK(field, params, filename);
        }

        // for (int i = 0; i < params.Nx + 2; i++)
        // {
        //     for (int j = 0; j < params.Ny + 2; j++)
        //     {
        //         cout << field.p[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;
        veloctiyStarCalculator(field, params);

        poissonEquationSolver(field, params);
        velocityCorrector(field, params);
        swapFields(field, params);
        // if (count == params.numberOfTimeSteps - 1)
        // {
        //     for (int i = 0; i < params.Nx + 2; i++)
        //     {
        //         cout << field.v[i][params.Ny / 2] << " ";
        //     }
        // }

        count++;
    }
    cout << params.timeStepSize << endl;
    cout << "Simulation completed. Data saved in VTK format." << endl;
    return 0;
}