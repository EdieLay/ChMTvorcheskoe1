#include <iostream>
#include <fstream>

void PutPxAndFxInFile(double (*f)(double), double* Px, double x0, double h, int n, int N);
double* PrintPxAndErrorRate(double (*f)(double), double x0, double h, int n, int N);
double CountErrorRate(double (*f)(double), double x0, double h, int n, int N);
double CountPxValue(double* Px, double x, int n);
double* FindPx(double (*f)(double), double h, double x0, int n);
double* SolveGauss(double** a, double* b, int n);
void FindMaxAndSwap(double** a, double* b, int k, int n);

int main()
{
    double (*f)(double) = cosh;
    double h = 2;
    double x0 = -8;
    int n = 6;
    int N = 100;
    
    double* Px = PrintPxAndErrorRate(f, x0, h, n, N);
    PutPxAndFxInFile(f, Px, x0, h, n, N);
    
    return 0;
}

void PutPxAndFxInFile(double (*f)(double), double* Px, double x0, double h, int n, int N)
{
    std::ofstream out;
    out.open("points.txt");
    if (out.is_open())
    {
        std::cout << "File Opened" << std::endl;
        out << N << std::endl;
        double H = n * h / N;

        for (int i = 0; i <= N; i++)
            out << x0 + i * H << " ";
        out << std::endl;

        for (int i = 0; i <= N; i++)
            out << f(x0 + i * H) << " ";
        out << std::endl;

        for (int i = 0; i <= N; i++)
            out << CountPxValue(Px, x0 + i * H, n) << " ";
        out << std::endl;
    }
    else std::cout << "Not opened" << std::endl;
}

double* PrintPxAndErrorRate(double (*f)(double), double x0, double h, int n, int N)
{
    double* Px = FindPx(f, x0, h, n);
    for (int i = n; i >= 0; i--)
        std::cout << Px[i] << " "; 
    
    std::cout << std::endl << CountErrorRate(f, x0, h, n, N) << std::endl;
    return Px;
}

double CountErrorRate(double (*f)(double), double x0, double h, int n, int N)
{
    double* Px = FindPx(f, x0, h, n);
    double maxDif = 0;
    double dif;
    double x;

    double H = n * h / N;
    for (int i = 0; i <= N; i++)
    {
        x = x0 + i * H;
        dif = abs(CountPxValue(Px, x, n) - f(x));
        maxDif = dif > maxDif ? dif : maxDif;
    }
    return maxDif;
}

double CountPxValue(double* Px, double x, int n)
{
    double result = 0;
    for (int i = 0; i <= n; i++)
    {
        result += Px[i] * pow(x, i);
    }
    return result;
}

double* FindPx(double (*f)(double), double x0, double h, int n) // n - степень многочлена, т.е. n+1 переменная
{
    double* x = new double[n + 1];
    double* fx = new double[n + 1];
    for (int i = 0; i <= n; i++)
    {
        x[i] = x0 + h * i;
        fx[i] = f(x[i]);
    }

    double** a = new double* [n + 1];
    for (int i = 0; i <= n; i++)
    {
        a[i] = new double[n + 1];
        for (int j = 0; j <= n; j++)
        {
            a[i][j] = pow(x[i], j);
        }
    }

    return SolveGauss(a, fx, n);
}


double* SolveGauss(double** a, double* b, int n)
{
    for (int k = 0; k <= n; k++) // прямой ход (по столбцам)
    {
        FindMaxAndSwap(a, b, k, n);
        for (int j = k + 1; j <= n; j++) // идем по строкам
        {
            double d = a[j][k] / a[k][k]; // для текущей строки вычисляем делитель для приведения к верхнетреугольному виду
            for (int i = k; i <= n; i++)
            {
                a[j][i] = a[j][i] - d * a[k][i]; // вычитаем текущую строку, умноженную на делитель, из всех строк ниже её
            }
            b[j] = b[j] - d * b[k]; // для столбца значений то же самое
        }
    }

    double* x = new double[n + 1]; // массив, где будет храниться решение

    for (int k = n; k >= 0; k--) // обратный ход
    {
        double d = 0;
        for (int j = k + 1; j <= n; j++)
        {
            double s = a[k][j] * x[j];
            d = d + s;
        }
        x[k] = (b[k] - d) / a[k][k];
    }

    return x;
}

void FindMaxAndSwap(double** a, double* b, int k, int n)
{
    double max = abs(a[k][k]);
    int max_ind = k;
    for (int i = k + 1; i <= n; i++)
    {
        double cur = abs(a[i][k]);
        if (cur > max)
        {
            max = cur;
            max_ind = i;
        }
    }
    if (max_ind != k)
    {
        double temp;
        for (int i = k; i <= n; i++)
        {
            temp = a[k][i];
            a[k][i] = a[max_ind][i];
            a[max_ind][i] = temp;
        }
        temp = b[k];
        b[k] = b[max_ind];
        b[max_ind] = temp;
    }
}