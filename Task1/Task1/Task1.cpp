#include <iostream>
#include <locale>
#include <vector>
#include <iomanip>

double f(double x)
{
    return x - 3 * cos(1.04 * x) * cos(1.04 * x);
}

double f_prime(double x)
{
    return 1 + 3.12 * sin(2.08 * x);
}

double f_dprime(double x)
{
    return 6.4896 * cos(2.08 * x);
}

using namespace std;

int main()
{
    setlocale(0, "RUS");
    double A = 0, B = 3.5, h, e = 1e-08, ai, bi, c, x, delta, y;
    int N = 5, counter = 0, temp = 0, k = 3, cmnd;
    const int p = 1;
    vector<double> an;
    cout << "Численные методы решения нелинейных уравнений" << endl << endl
        << "Параметры задачи: отрезок [" << A << "; " << B << "], точность корня e = " << e << endl << endl
        << "Находятся корни уравнения x - 3cos^2(1.04 * x) = 0" << endl << endl;
    h = (B - A) / N;
    ai = A; bi = A + h;
    while (temp != counter || counter == 0)
    {
        temp = counter;
        counter = 0;
        while (bi <= B)
        {
            if (f(ai) * f(bi) <= 0)
            {
                ++counter;
                an.push_back(ai);
            }
            ai = bi; bi += h;
        }
        if (temp != counter || counter == 0)
        {
            N *= 2;
            h = (B - A) / N;
            an.clear();
        }
        ai = A; bi = A + h;
    }
    cout << "Отрезков перемены знака: " << counter << endl << "Отрезки:" << endl << endl;
    temp = 1;
    for (auto i:an)
    {
        cout << temp++ << ") [" << i << "; " << i + h << "]" << endl;
    }
    while (1)
    {
        cout << endl << "Выберите отрезок для нахождения корня (0 для завершения работы): ";
        cin >> cmnd;
        if (!cmnd)
            break;
        ai = an[cmnd - 1]; bi = ai + h;

        //Бисекция

        temp = 0;
        cout << endl << "Приближение методом бисекции:" << endl;
        while (bi - ai > 2 * e)
        {
            c = (bi + ai) / 2;
            if (f(ai) * f(c) <= 0)
                bi = c;
            else
                ai = c;
            ++temp;
            
        }
        x = (ai + bi) / 2; delta = (bi - ai) / 2;
        cout << "Количество шагов: " << temp << endl
            << "Приближенный корень: " << setprecision(12) << x << endl
            << "Отличие от корня: |x - x*| <= " << delta << endl
            << "Величина невязки: " << abs(f(x)) << endl << endl;

        //Метод Ньютона

        temp = 1;
        ai = an[cmnd - 1]; bi = ai + h;
        cout << endl << "Приближение методом Ньютона:" << endl;
        if (f(ai) * f_dprime(ai) > 0 && f_prime(ai) != 0)
            c = ai;
        else
            c = bi;
        cout << "Первое приближение: x0 = " << setprecision(12) <<c << endl;
        x = c - p * (f(c) / f_prime(c));
        while (abs(x - c) >= e)
        {
            c = x;
            x = c - p * (f(c) / f_prime(c));
            ++temp;
        }
        delta = abs(c - x);
        cout << "Количество шагов: " << temp << endl
            << "Приближенный корень: " << setprecision(12) << x << endl
            << "|x_" << temp << " - x_" << temp - 1 << "| = " << delta << endl
            << "Величина невязки: " << abs(f(x)) << endl << endl;

        //Модифицированный метод Ньютона

        temp = 1;
        ai = an[cmnd - 1]; bi = ai + h;
        cout << endl << "Приближение модифицированным методом Ньютона:" << endl;
        if (f(ai) * f_dprime(ai) > 0 && f_prime(ai) != 0)
            c = ai;
        else
            c = bi;
        cout << "Первое приближение: x0 = " << setprecision(12) << c << endl;
        y = f_prime(c);
        x = c - f(c) / y;
        while (abs(x - c) >= e)
        {
            c = x;
            x = c - f(c) / y;
            ++temp;
        }
        delta = abs(c - x);
        cout << "Количество шагов: " << temp << endl
            << "Приближенный корень: " << setprecision(12) << x << endl
            << "|x_" << temp << " - x_" << temp - 1 << "| = " << delta << endl
            << "Величина невязки: " << abs(f(x)) << endl << endl;

        //Метод секущих

        temp = 2;
        ai = an[cmnd - 1]; bi = ai + h;
        cout << endl << "Приближение методом секущих:" << endl;
        c = ai; y = bi;
        cout << "Первые приближения: x0 = " << setprecision(12) << c << ", x1 = " << y << endl;
        x = y - p * (f(y) * (y - c)) / (f(y) - f(c));
        while (abs(x - y) >= e)
        {
            c = y;
            y = x;
            x = y - p * (f(y) * (y - c)) / (f(y) - f(c));
            ++temp;
        }
        delta = abs(y - x);
        cout << "Количество шагов: " << temp - 1 << endl
            << "Приближенный корень: " << setprecision(12) << x << endl
            << "|x_" << temp << " - x_" << temp - 1 << "| = " << delta << endl
            << "Величина невязки: " << abs(f(x)) << endl << endl;
    }
    return 0;
}
