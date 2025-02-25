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
        x.resize(Ny + 1, vector<double>(Nx + 1));
        y.resize(Ny + 1, vector<double>(Nx + 1));
        xm.resize(Ny, vector<double>(Nx));
        ym.resize(Ny, vector<double>(Nx));
        u.resize(Ny, vector<double>(Nx));
        uStar.resize(Ny, vector<double>(Nx));
        uNew.resize(Ny, vector<double>(Nx));
        vStar.resize(Ny, vector<double>(Nx));
        vNew.resize(Ny, vector<double>(Nx));
        v.resize(Ny, vector<double>(Nx));
        p.resize(Ny, vector<double>(Nx));
    }
};

// function to create the coordinates of the nodes of cells
void createCoordinatesXY(vector<vector<double>> &x, vector<vector<double>> &y, constParameters params)
{
    for (int j = 0; j < params.Ny + 1; j++)
    {
        for (int i = 0; i < params.Nx + 1; i++)
        {
            x[j][i] = i * params.hx;
            y[j][i] = (params.Ny - j) * params.hy;
        }
    }
}

// function to create the coordinates of the cell centers
void createCoordinatesXYM(vector<vector<double>> &xm, vector<vector<double>> &ym, constParameters params)
{
    for (int j = 0; j < params.Ny; j++)
    {
        for (int i = 0; i < params.Nx; i++)
        {
            xm[j][i] = params.hx / 2 + i * params.hx;
            ym[j][i] = params.hy / 2 + (params.Ny - 1 - j) * params.hy;
        }
    }
}

// function to calculate the first order partial derivative using central differencing
double firstOrderPDEcentralDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j + y][i + x] - variable[j - y][i - x]) / (x * params.hx + y * params.hy);
    return pde;
}

// function to calculate the second order partial derivative using central differencing
double secondOrderPDEcentralDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j + y][i + x] - 2 * variable[j][i] + variable[j - y][i - x]) / pow((x * params.hx + y * params.hy), 2);
    return pde;
}

double firstOrderPDEforwardDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j + y][i + x] - variable[j][i]) / (x * params.hx + y * params.hy);
    return pde;
}

double firstOrderPDEbackwardDiff(vector<vector<double>> &variable, int j, int i, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[j][i] - variable[j - y][i - x]) / (x * params.hx + y * params.hy);
    return pde;
}

// PREDICTOR STEP
void veloctiyStarCalculator(fields &field, constParameters params)
{
    for (int i = 1; i < params.Nx - 1; i++)
    {
        for (int j = 1; j < params.Ny - 1; j++)
        {
            field.uStar[j][i] = field.u[j][i] + params.timeStepSize * (params.kinematicViscosity * (secondOrderPDEcentralDiff(field.u, j, i, 1, 0, params) + secondOrderPDEcentralDiff(field.u, j, i, 0, 1, params)) - (field.u[j][i] * firstOrderPDEcentralDiff(field.u, j, i, 1, 0, params) + 0.25 * (field.v[j][i - 1] + field.v[j][i] + field.v[j - 1][i - 1] + field.v[j - 1][i]) * firstOrderPDEcentralDiff(field.u, j, i, 0, 1, params)));
            field.vStar[j][i] = field.v[j][i] + params.timeStepSize * (params.kinematicViscosity * (secondOrderPDEcentralDiff(field.v, j, i, 1, 0, params) + secondOrderPDEcentralDiff(field.v, j, i, 0, 1, params)) - (field.v[j][i] * firstOrderPDEcentralDiff(field.v, j, i, 0, 1, params) + 0.25 * (field.u[j + 1][i] + field.u[j][i] + field.u[j + 1][i + 1] + field.u[j][i + 1]) * firstOrderPDEcentralDiff(field.v, j, i, 1, 0, params)));
        }
    }
}
// POISSON EQUATION SOLVER

void poissonEquationSolver(fields &field, constParameters params)
{
    for (int k = 0; k < 5; k++)
    {
        for (int i = 1; i < params.Nx - 1; i++)
        {
            for (int j = 1; j < params.Ny - 1; j++)
            {
                field.p[j][i] = (pow(params.hx, 2) * pow(params.hx, 2) / (2 * (pow(params.hx, 2) + pow(params.hy, 2)))) * (((field.p[j][i + 1] + field.p[j][i - 1]) / pow(params.hx, 2)) + ((field.p[j - 1][i] + field.p[j + 1][i]) / pow(params.hy, 2)) - ((params.density / params.timeStepSize) * (firstOrderPDEforwardDiff(field.uStar, j, i, 1, 0, params) + firstOrderPDEforwardDiff(field.vStar, j, i, 0, 1, params))));
            }
        }
    }
}

// CORRECTOR STEP
void velocityCorrector(fields &field, constParameters params)
{
    for (int i = 1; i < params.Nx - 1; i++)
    {
        for (int j = 1; j < params.Ny - 1; j++)
        {
            field.uNew[j][i] = field.uStar[j][i] - (params.timeStepSize / params.density) * (firstOrderPDEbackwardDiff(field.p, j, i, 1, 0, params));
            field.vNew[j][i] = field.vStar[j][i] - (params.timeStepSize / params.density) * (firstOrderPDEbackwardDiff(field.p, j, i, 0, 1, params));
        }
    }
}

void swapFields(fields &field, constParameters params)
{
    for (int i = 0; i < params.Nx; i++)
    {
        for (int j = 0; j < params.Ny; j++)
        {
            field.u[j][i] = field.uNew[j][i];
            field.v[j][i] = field.vNew[j][i];
        }
    }
}

void setBoundaryConditions(int b, vector<vector<double>> &M, constParameters params)
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
            M[0][i] = 0;
            M[params.Ny - 1][i] = M[params.Ny - 2][i];
        }
    }
    if (b == 1)
    {
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * params.uLeftWall - M[j][1];
            M[j][params.Nx - 1] = 2 * params.uRightWall - M[j][params.Nx - 2];
        }
        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * params.uTopWall - M[1][i];
            M[params.Ny - 1][i] = 2 * params.uBottomWall - M[params.Ny - 2][i];
        }
    }
    if (b == 2)
    {
        for (int j = params.Ny - 1; j >= 0; j--)
        {
            M[j][0] = 2 * params.vLeftWall - M[j][1];
            M[j][params.Nx - 1] = 2 * params.vRightWall - M[j][params.Nx - 2];
        }
        for (int i = params.Nx - 1; i >= 0; i--)
        {
            M[0][i] = 2 * params.vTopWall - M[1][i];
            M[params.Ny - 1][i] = 2 * params.vBottomWall - M[params.Ny - 2][i];
        }
    }

    M[0][0] = 0.5 * (M[1][0] + M[0][1]);
    M[0][params.Nx - 1] = 0.5 * (M[0][params.Nx - 2] + M[1][params.Nx - 1]);
    M[params.Ny - 1][0] = 0.5 * (M[params.Ny - 2][0] + M[params.Ny - 1][1]);
    M[params.Ny - 1][params.Nx - 1] = 0.5 * (M[params.Ny - 2][params.Nx - 1] + M[params.Ny - 1][params.Nx - 2]);
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
    vtkFile << "DIMENSIONS " << params.Nx << " " << params.Ny << " 1\n";
    vtkFile << "POINTS " << (params.Nx) * (params.Ny) << " float\n";

    // Write the grid coordinates
    for (int j = 0; j < params.Ny; j++)
    {
        for (int i = 0; i < params.Nx; i++)
        {
            vtkFile << field.xm[j][i] << " " << field.ym[j][i] << " 0.0\n";
        }
    }

    // Write the velocity and pressure data
    vtkFile << "POINT_DATA " << (params.Nx) * (params.Ny) << "\n";
    vtkFile << "VECTORS velocity float\n";
    for (int j = 0; j < params.Ny; j++)
    {
        for (int i = 0; i < params.Nx; i++)
        {
            vtkFile << field.u[j][i] << " " << field.v[j][i] << " 0.0\n";
        }
    }

    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < params.Ny; j++)
    {
        for (int i = 0; i < params.Nx; i++)
        {
            vtkFile << field.p[j][i] << "\n";
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

    params.uTopWall = 1;
    params.uBottomWall = 0.0;
    params.uLeftWall = 0.0;
    params.uRightWall = 0.0;

    params.vTopWall = 0;
    params.vBottomWall = 0.0;
    params.vLeftWall = 0.0;
    params.vRightWall = 0.0;

    params.Nx = 32;
    params.Ny = 32;
    params.lengthX = 1;
    params.lengthY = 1;
    params.hx = params.lengthX / (params.Nx);
    params.hy = params.lengthY / (params.Ny);

    params.startTime = 0;
    params.endTime = 5;
    params.timeStepSize = params.courantNumber * min(params.hx, params.hy) / params.uTopWall;
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

        veloctiyStarCalculator(field, params);
        poissonEquationSolver(field, params);
        velocityCorrector(field, params);
        swapFields(field, params);
        setBoundaryConditions(1, field.u, params);
        setBoundaryConditions(2, field.v, params);
        setBoundaryConditions(0, field.p, params);
        if (count == params.numberOfTimeSteps - 1)
        {
            // for (int j = 0; j < params.Ny + 1; j++)
            // {
            //     for (int i = 0; i < params.Nx + 1; i++)
            //     {
            //         cout << field.x[j][i] << "," << field.y[j][i] << " ";
            //     }
            //     cout << endl;
            // }

            for (int j = 0; j < params.Ny + 1; j++)
            {
                cout << params.Ny - j << " " << field.x[j][(params.Nx) / 2] << ", " << field.y[j][(params.Nx) / 2] << ", " << field.u[j][(params.Nx) / 2] << endl;
            }
            cout << endl;
        }
        count++;
    }
    // if (count == params.numberOfTimeSteps) // Save every 10 time steps
    // {
    //     string filename = "output_" + to_string(count) + ".vtk";
    //     saveToVTK(field, params, filename);
    // }

    cout << params.timeStepSize << endl;
    cout << "Simulation completed. Data saved in VTK format." << endl;
    return 0;
}