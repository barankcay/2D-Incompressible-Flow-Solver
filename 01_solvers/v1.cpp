#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
const int N = 4;

const double lengthX = 12.0;
const double lengthY = 12.0;

double hx = lengthX / (N + 2);
double hy = lengthY / (N + 2);

vector<vector<double>> x(N + 3, vector<double>(N + 3));
vector<vector<double>> y(N + 3, vector<double>(N + 3));

vector<vector<double>> xm(N + 2, vector<double>(N + 2));
vector<vector<double>> ym(N + 2, vector<double>(N + 2));

vector<vector<double>> u(N + 2, vector<double>(N + 2));
vector<vector<double>> v(N + 2, vector<double>(N + 2));

vector<vector<double>> p(N + 2, vector<double>(N + 2));

void createCoordinatesXY()
{
    for (int i = 0; i < N + 3; i++)
    {
        for (int j = 0; j < N + 3; j++)
        {
            x[i][j] = i * hx;
            y[i][j] = j * hy;
        }
    }
}

void createCoordinatesXYM()
{
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            xm[i][j] = hx / 2 + i * hx;
            ym[i][j] = hy / 2 + j * hy;
        }
    }
}
int main()
{
    createCoordinatesXY();
    createCoordinatesXYM();
    // for (int i = 0; i < N + 3; i++)
    // {
    //     for (int j = 0; j < N + 3; j++)
    //     {
    //         cout << x[i][j] << "," << y[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    return 0;
}