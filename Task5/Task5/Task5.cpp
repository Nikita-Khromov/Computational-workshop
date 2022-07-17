#include <iostream>
#include <list>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <locale>

using namespace std;

typedef list<long double> polinomial;     //deg = size, begining deg = 0
typedef pair<long double, int> monom;    //first - coeff, second - deg

polinomial mono_to_poli(const monom& m)
{
    polinomial res(m.second - 1, 0);
    res.push_back(m.first);
    return res;
}

polinomial operator*(long double k, const polinomial& p)
{
    polinomial res;
    if (k == 0)
        return res;
    for (auto i : p)
    {
        res.push_back(i * k);
    }
    return res;
}

monom operator*(long double k, const monom& m)
{
    monom res = m;
    res.first *= k;
    if (k == 0)
        res.second = 0;
    return res;
}

monom operator*(const monom& m1, const monom& m2)
{
    monom res;
    if (m1.first == 0 || m2.first == 0)
        return { 0, 0 };
    res.first = m1.first * m2.first;
    res.second = m1.second + m2.second;
    return res;
}

polinomial operator*(const polinomial& p, const monom& m)
{
    polinomial res;
    if (m.second == 0)
        return m.first * p;
    if (p.empty())
        return res;
    for (int i = 0; i < m.second; ++i)
    {
        res.push_back(0);
    }
    for (auto i : p)
    {
        res.push_back(i * m.first);
    }
    return res;
}



polinomial operator+(const polinomial& p1, const polinomial& p2)
{
    polinomial res;
    int n;
    n = (p1.size() > p2.size() ? p1.size() : p2.size());
    polinomial::const_iterator i1 = p1.cbegin();
    polinomial::const_iterator i2 = p2.cbegin();
    for (int i = 0; i < n; ++i)
    {
        if (i1 == p1.end())
        {
            for (int j = i; j < n; ++j)
            {
                res.push_back(*i2);
                i2++;
            }
            break;
        }
        if (i2 == p2.end())
        {
            for (int j = i; j < n; ++j)
            {
                res.push_back(*i1);
                i1++;
            }
            break;
        }
        res.push_back(*i1 + *i2);
        ++i1, ++i2;
    }
    if (res == polinomial(n, 0) || res.empty())
    {
        res.clear();
        return res;
    }
    res.reverse();
    for (polinomial::iterator it = res.begin(); *it == 0; it = res.begin())
        res.erase(it);
    res.reverse();
    return res;
}

polinomial operator+(const monom& m1, const monom& m2)
{
    polinomial p1, p2;
    p1 = mono_to_poli(m1); p2 = mono_to_poli(m2);
    return p1 + p2;
}

polinomial operator-(const polinomial& p1, const polinomial& p2)
{
    polinomial temp;
    temp = -1 * p2;

    return (p1 + temp);
}

polinomial operator-(const monom& m1, const monom& m2)
{
    polinomial p1, p2;
    p1 = mono_to_poli(m1); p2 = mono_to_poli(m2);
    return p1 - p2;
}

polinomial operator+(const polinomial& p, const monom m)
{
    polinomial p2 = mono_to_poli(m);
    return p2 + p;
}

polinomial operator-(const polinomial& p, const monom m)
{
    polinomial p2 = mono_to_poli(m);
    return p - p2;
}

ostream& operator<<(ostream& stream, const polinomial& p)
{
    int deg = 0; bool fl = false;
    for (auto i : p)
    {
        if (i != 0)
        {
            if (deg == 0)
            {
                stream << i << " ";
                fl = true;
            }
            else if (deg == 1)
            {
                if (i < 0)
                    stream << "- ";
                else if (fl)
                    stream << "+ ";
                if (i != 1)
                    stream << abs(i) << "x ";
                else
                    stream << "x ";
                fl = true;
            }
            else
            {
                if (i < 0)
                    stream << "- ";
                else if (fl)
                    stream << "+ ";
                if (i != 1)
                    stream << abs(i) << "x^" << deg << " ";
                else
                    stream << "x^" << deg << " ";
                fl = true;
            }
        }
        ++deg;
    }
    return stream;
}

// Вычисление значения в точке
long double operator^(const polinomial& p, const long double& x) {
    polinomial::const_iterator i = p.cend(); --i;
    polinomial::const_iterator end = p.cbegin();
    long double res = 0;


    for (; i != end;) {
        res = x * (res + *i);
        i--;
    }
    res += *i;

    return res;
}

polinomial next_Legendre(const vector<polinomial>& Legendre) {
    polinomial result_p; const size_t n = Legendre.size();
    monom x = make_pair(1, 1);
    result_p = (((long double)2 * n - 1) / n) * Legendre[n - 1] * x - (((long double)n - 1) / n) * Legendre[n - 2];
    return result_p;
}

void fill_Legendre(vector<polinomial>& Legendre, size_t n) {
    for (int i = 2; i < n; ++i) {
        Legendre.push_back(next_Legendre(Legendre));
    }
}

vector<long double> find_roots(const polinomial& p, long double a, long double b, long double e) {
    size_t N = 9;
    long double h = (b - a) / (2 * N) - e;
    long double ai = a, bi = a + h;
    vector<long double> an, roots;
    while (bi <= b)
    {
        if ((p ^ ai) * (p ^ bi) <= 0)
        {
            an.push_back(ai);
        }
        ai = bi; bi += h;
    }

    for (auto cmnd : an) {
        ai = cmnd; bi = ai + h;
        long double c = ai, y = bi;
        long double x = y - ((p ^ y) * (y - c)) / ((p ^ y) - (p ^ c));
        while (abs(x - y) >= e) {
            c = y;
            y = x;
            x = y - ((p ^ y) * (y - c)) / ((p ^ y) - (p ^ c));
        }
        roots.push_back(x);
    }
    return roots;
}

vector<long double> Gauss_coefficients(const vector<long double>& nodes, const vector<polinomial>& Legendre, size_t N)
{
    vector<long double> coefficients;
    long double A;
    for (auto xk : nodes) {
        A = (2 * (1 - xk * xk)) / (N * N * (Legendre[N - 1] ^ xk) * (Legendre[N - 1] ^ xk));
        coefficients.push_back(A);
    }

    return coefficients;
}

long double f(long double z)
{
    return (z + 0.8) / sqrt(z * z + 1.2);
}

long double f_integral(long double z)
{
    return (long double)1 / 5 * sqrt(25 * z * z + 30) + (long double)4 / 5 * log(abs(z + sqrt(25 * z * z + 30) / 5));
}

//long double weight_ro(long double z)
//{
//    return 1 / sqrt(1 - z * z);
//}

long double Mueller_f(long double z)
{
    return exp(2 * z) * z * z;
}

int main()
{
    setlocale(0, "RUS");
    vector<polinomial> Legendre = { { 1 }, { 0, 1 } };
    const vector<size_t> Ns = { 4, 5, 7, 8 };
    const vector<polinomial> test = { {1, 5, 3, -2, 11, -3, 8, 1}, {0, -3, 0, 8, -1, 0, 0, 4, 3, -3},
        {41, -12, 14, -7, -8, 0, 0, 0, 5, -3, 1, 8, -6, 5}, {3, 1, 5, 3, 2, 9, 5, -1, 12, 5, -6, 7, 0, 1, -2, 1} },
        integral_test = { {0, 1, (long double) 5 / 2, 1, -0.5, (long double) 11 / 5, -0.5, (long double) 8 / 7, (long double) 1 / 8},
        {0, 0, (long double)  - 3 / 2, 0, 2, (long double) - 1 / 5, 0, 0, 0.5, (long double) 1 / 3, -0.3},
        {0, 41, -6, (long double) 14 / 3, (long double) - 7 / 4, (long double) - 8 / 5, 0, 0, 0, (long double) 5 / 9, -0.3, (long double) 1 / 11, (long double) 2 / 3, (long double) - 6 / 13, (long double) 5 / 14},
        {0, 3, 0.5, (long double) 5 / 3, 0.75, 0.4, 1.5, (long double) 5 / 7, (long double) - 1 / 8, (long double) 4 / 3, 0.5, (long double) - 6 / 11, (long double) 7 / 12, 0, (long double) 1 / 14, (long double)  - 2 / 15, (long double) 1 / 16} };
    const size_t N = 8; 
    const long double e = 1e-12, A = 1.6, B = 2.7, pi = atan(1) * 4; 
    long double integral, real_integral;
    vector<size_t> Mueller_Ns(3);
    vector<vector<long double>> coefficients_All, roots_All, Mueller_nodes(3), Mueller_coefficients(3);
    vector<vector<long double>> new_roots(Ns.size()), new_coefficients(Ns.size());
    fill_Legendre(Legendre, N + 1);
    cout << setprecision(16) << showpos << left;
    for (int i = 1; i <= N; ++i) {
        vector<long double> roots = find_roots(Legendre[i], -1, 1, e);
        roots_All.push_back(roots);
        vector<long double> coefficients = Gauss_coefficients(roots, Legendre, i);
        coefficients_All.push_back(coefficients);
        cout << "N = " << i << ":" << endl;
        for (int j = 0; j < roots.size(); ++j) {
            cout << setw(22) << roots[j] << "  --->  " << setw(22) << coefficients[j] << endl;
        }
        cout << endl;
    }
    cout << noshowpos;
    cout << endl << "Проверка точности соответствующих КФ Гаусса на многочленах степеней 4, 5, 7, 8 на отрезке [-1, 1]:" << endl;
    for (int j = 0; j < Ns.size(); ++j) {
        integral = 0; real_integral = (integral_test[j] ^ 1) - (integral_test[j] ^ (-1));
        cout << "N = " << Ns[j] << ":" << endl;
        cout << "Многочлен: " << test[j] << endl << "Его интеграл: " << real_integral << endl;
        for (int i = 0; i < Ns[j]; ++i) {
            integral += (test[j] ^ roots_All[Ns[j] - 1][i]) * coefficients_All[Ns[j] - 1][i];
        }
        cout << "Интеграл, посчитанный по КФ Гаусса: " << integral << endl
            << "Фактическая погрешность: " << abs(integral - real_integral) << endl << endl;
    }

    cout << endl << "Функция: (x + 0.8) / sqrt(x^2 + 1.2)   Отрезок: [" << A << ", " << B << "]" << endl << endl;
    for (int j = 0; j < Ns.size(); ++j)
    {
        for (int i = 0; i < Ns[j]; ++i)
        {
            new_roots[j].push_back(roots_All[Ns[j] - 1][i] * ((B - A) / 2) + (B + A) / 2);
            new_coefficients[j].push_back(coefficients_All[Ns[j] - 1][i] * ((B - A) / 2));
        }
        cout << "Новые узлы и коэффициенты для N = " << Ns[j] << ":" << showpos << endl;
        for (int k = 0; k < Ns[j]; ++k) {
            cout << setw(22) << new_roots[j][k] << "  --->  " << setw(22) << new_coefficients[j][k] << endl;
        }
    }

    cout << endl << noshowpos << setprecision(18);
    real_integral = f_integral(B) - f_integral(A);
    cout << "Фактическое значение интеграла: " << real_integral << endl << "Значения полученные при помощи КФ Гаусса:" << endl;
    for (int j = 0; j < Ns.size(); ++j) {
        cout << "N = " << Ns[j] << ": ";
        integral = 0;
        for (int i = 0; i < Ns[j]; ++i)
        {
            integral += (f(new_roots[j][i])) * new_coefficients[j][i];
        }
        cout << integral << endl;
    }

    cout << endl << "Фактическая погрешность КФ Гаусса для N = " << Ns[Ns.size() - 1] << ": " << abs(integral - real_integral) << endl << endl;

    cout << "Составление КФ Мелера при f(x) = e^(2x) * x^2" << endl << endl;

    cout << "Введите N1, N2, N3 для расчёта КФ Мелера: "; cin >> Mueller_Ns[0] >> Mueller_Ns[1] >> Mueller_Ns[2];

    cout << endl << "Узлы и коэффициенты для соответствующих N:" << endl;

    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < Mueller_Ns[j]; ++i)
        {
            Mueller_nodes[j].push_back(cos(pi * (2 * (i + 1) - 1) / (2 * Mueller_Ns[j])));
            Mueller_coefficients[j].push_back(pi / Mueller_Ns[j]);
        }
        cout << "N = " << Mueller_Ns[j] << ":" << endl;
        cout << showpos;
        for (int i = 0; i < Mueller_Ns[j]; ++i)
        {
            cout << setw(24) << Mueller_nodes[j][i] << " ---> " << setw(24) << Mueller_coefficients[j][i] << endl;
        }
        cout << endl << noshowpos;
    }

    for (int j = 0; j < 3; ++j)
    {
        integral = 0;
        for (int i = 0; i < Mueller_Ns[j]; ++i)
        {
            integral += Mueller_coefficients[j][i] * Mueller_f(Mueller_nodes[j][i]);
        }
        cout << "N = " << setw(2) << Mueller_Ns[j] << ": " << setw(22) << integral << endl;
    }

    return 0;
}