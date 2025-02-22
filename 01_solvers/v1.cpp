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
    double kinematicViscosity = 0.01;
    double density = 1.0;
};

struct fields
{
    vector<vector<double>> x;
    vector<vector<double>> y;
    vector<vector<double>> xm;
    vector<vector<double>> ym;
    vector<vector<double>> u;
    vector<vector<double>> uStar;
    vector<vector<double>> v;
    vector<vector<double>> vStar;
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

double firstOrderPDEcentralDiff(vector<vector<double>> &variable, int i, int j, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[i + x][j + y] - variable[i - x][j - y]) / (2 * (x * params.hx + y * params.hy));
    return pde;
}

double secondOrderPDEcentralDiff(vector<vector<double>> &variable, int i, int j, int x, int y, constParameters params)
{ // x and y are the direction of the derivative
    // i and j are the indices of the variable
    // if x = 1 and y = 0, then the derivative is in the x direction
    // if x = 0 and y = 1, then the derivative is in the y direction
    double pde;
    pde = (variable[i + x][j + y] - 2 * variable[i][j] + variable[i - x][j - y]) / pow((x * params.hx + y * params.hy), 2);
    return pde;
}

void StarCalculator(fields &field, constParameters params)
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

int main()
{
    constParameters params;
    params.Nx = 3;
    params.Ny = 3;
    params.lengthX = 1.0;
    params.lengthY = 1.0;
    params.time = 20;
    params.timeStepSize = 0.1;
    params.hx = params.lengthX / (params.Nx + 2);
    params.hy = params.lengthY / (params.Ny + 2);

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