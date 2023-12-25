#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

double f(double x, double y);
double rf(double x);
void CountRealY(double* x, double* ry, int n);
void FindAllK(double h, double (*f)(double, double), double x, double y, double* k);
double FindNextY(double y, double* k);
void RungeKuttaMethod(double* x, double* y, double h, int n, double (*f)(double, double));
double CountErrorRate(double* x, double* y, double h, int n, double (*f)(double, double));
void Print_h_error_xy(double h, double error, double* x, double* y, int n);
void PutCountedYInFile(double* x, double* y, int n);
void PutCountedAndRealYInFile(double* x, double* y, double* ry, int n);

int main()
{
    int n = 100;
    double a = 0;
    double b = 5;
    double h = (b - a) / n;
    double* x = new double[n + 1];
    double* y = new double[n + 1];
    double* ry = new double[n + 1];
    x[0] = a;
    y[0] = 0;
    RungeKuttaMethod(x, y, h, n, f);
    //CountRealY(x, ry, n);
    double error = CountErrorRate(x, y, h, n, f);
    Print_h_error_xy(h, error, x, y, n);
    PutCountedYInFile(x, y, n);
    //PutCountedAndRealYInFile(x, y, ry, n);
}


void RungeKuttaMethod(double* x, double* y, double h, int n, double (*f)(double, double))
{
    double* k = new double[4];
    for (int i = 1; i <= n; i++)
    {
        x[i] = x[i - 1] + h;
        FindAllK(h, f, x[i - 1], y[i - 1], k);
        y[i] = FindNextY(y[i - 1], k);
    }
}

double CountErrorRate(double* x, double* y, double h, int n, double (*f)(double, double))
{
    double h2 = h / 2;
    int n2 = 2 * n;
    double* x2 = new double[n2 + 1];
    double* y2 = new double[n2 + 1];
    x2[0] = x[0];
    y2[0] = y[0];
    RungeKuttaMethod(x2, y2, h2, n2, f);
    double error = 0.0;
    double curError;
    for (int i = 1; i <= n; i++)
    {
        curError = abs(y[i] - y2[2 * i]);
        error = curError > error ? curError : error;
    }
    return error / 15;
}

double f(double x, double y)
{
    return sin(x * x) + cos(y);
    //return x + y;
    //return 0.2 * x * x + 0.5 * y * y;
}

double rf(double x)
{
    return 2 * exp(x) - x - 1;
}

void CountRealY(double* x, double* ry, int n)
{
    for (int i = 0; i <= n; i++)
        ry[i] = rf(x[i]);
}

void FindAllK(double h, double (*f)(double, double), double x, double y, double* k)
{
    k[0] = h * f(x, y);
    k[1] = h * f(x + h / 2.0, y + k[0] / 2.0);
    k[2] = h * f(x + h / 2.0, y + k[1] / 2.0);
    k[3] = h * f(x + h, y + k[2]);
}

double FindNextY(double y, double* k)
{
    return y + (k[0] + 2.0 * k[1] + 2.0 * k[2] + k[3]) / 6.0;
}

void Print_h_error_xy(double h, double error, double* x, double* y, int n)
{
    cout << "-------------------------------------------" << endl;
    cout << "h = " << h << endl;
    cout << "error = " << setprecision(15) << error << endl;
    for (int i = 0; i <= n; i++)
    {
        cout << "(" << x[i] << "; " << setprecision(15) << y[i] << ") ";
    }
    cout << endl << "-------------------------------------------" << endl;
}


void PutCountedYInFile(double* x, double* y, int n)
{
    std::ofstream out;
    out.open("points.txt");
    if (out.is_open())
    {
        std::cout << "File Opened" << std::endl;
        out << n << std::endl;

        for (int i = 0; i <= n; i++)
            out << x[i] << " ";
        out << std::endl;

        for (int i = 0; i <= n; i++)
            out << setprecision(15) << y[i] << " ";
        out << std::endl;
    }
    else std::cout << "Not opened" << std::endl;
}

void PutCountedAndRealYInFile(double* x, double* y, double* ry, int n)
{
    std::ofstream out;
    out.open("points.txt");
    if (out.is_open())
    {
        std::cout << "File Opened" << std::endl;
        out << n << std::endl;

        for (int i = 0; i <= n; i++)
            out << x[i] << " ";
        out << std::endl;

        for (int i = 0; i <= n; i++)
            out << setprecision(15) << ry[i] << " ";
        out << std::endl;

        for (int i = 0; i <= n; i++)
            out << setprecision(15) << y[i] << " ";
        out << std::endl;
    }
    else std::cout << "Not opened" << std::endl;
}