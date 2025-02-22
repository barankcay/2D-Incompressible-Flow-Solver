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