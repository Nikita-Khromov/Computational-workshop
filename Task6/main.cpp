#include <iostream>
#include <iomanip>
#include <cmath>
#include "../../../eigen-3.4.0/Eigen/Dense"
#include "polynomial.h"

#define INT 0.104612085559286508

long double f(long double z) {
    return sinl(z);
}

long double q(long double z) {
    return -z * logl(z);
//    return sqrtl(z / (1 - z));
}


long double fq(long double z) {
    return f(z) * q(z);
}

long double w_sum(long double a, long double b, long double h, const polynomial &p) {
    long double res = 0, z = a + h / 2;
    while (b - z >= h / 4) {
        res += (p ^ z) * q(z);
        z += h;
    }
    return res;
}

long double polynomial_integral(const size_t N) {
    return (long double) 1 / ((2 * N + 1) * (2 * N + 1));
}

polynomial build_polynomial(size_t n) {
    const monom m = std::make_pair(1, n);
    polynomial res = mono_to_poli(m);
    return res;
}

using namespace std;
using namespace Eigen;

int main() {
    long double a, b;
    const long double e = 1e-18;
    int N, m;
    vector<long double> coefficients, roots;
    vector<polynomial> Legendre = {{1},
                                   {0, 1}};
    vector<long double> moments;

    cout << "Построение составной квадратурной формулы Гаусса" << endl
         << "Введите концы отрезка: a = ";
    cin >> a;
    cout << "b = ";
    cin >> b;
    cout << "Введите количество узлов для КФ Гаусса: N = ";
    cin >> N;
    cout << "Введите количество промежутков деления отрезка [" << a << ", " << b << "]: m = ";
    cin >> m;

    fill_Legendre(Legendre, N + 1);
    roots = find_roots(Legendre[N], -1, 1, e, N * 2);
    coefficients = Gauss_coefficients(roots, Legendre, N);

    cout << "Узлы и коэффициенты КФ Гаусса на отрезке [-1, 1]:" << endl;
    cout << setprecision(18) << showpos << left;

    for (int i = 0; i < roots.size(); ++i) {
        cout << setw(22) << roots[i] << " ---> " << coefficients[i] << endl;
    }

    cout << endl;

    long double h = (b - a) / m, integral = 0, temp_sum = 0;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < m; ++j) {
            temp_sum += fq((h / 2) * (roots[k] + 1) + (a + j * h));
        }
        integral += temp_sum * coefficients[k];
        temp_sum = 0;
    }
    integral *= h / 2;
    cout << noshowpos;
    cout << "Значение интеграла, полученное с помощью СКФ Гаусса: " << integral << endl
         << "Погрешность для отрезка [0, 1]: " << abs(integral - INT);

    cout << endl << endl << "Построение КФНАСТ степени для вычисления интеграла" << endl
         << "Введите количество узлов в формуле: N = ";
    cin >> N;
    h = (b - a) / (1e+06);

    polynomial test = build_polynomial(2 * N - 1);
    long double w;
    for (int i = 0; i < 2 * N; ++i) {
        w = w_sum(a, b, h, build_polynomial(i));
        moments.push_back(h * w);
    }
    cout << "Моменты весовой функции:" << endl;
    for (int i = 0; i < 2 * N; ++i) {
        cout << "m_" << i << " = " << moments[i] << endl;
    }

    cout << endl;

    vector<double> tmp1;
    for (int i = 0; i < N; ++i) {
        tmp1.push_back((double) -moments[N + i]);
    }

    Map<VectorXd> v(&tmp1[0], N);

    vector<double> tmp2;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            tmp2.push_back((double) moments[N + i - j - 1]);
        }
    }

    Map<MatrixXd> A(&tmp2[0], N, N);

    A.transposeInPlace();

    VectorXd wcoeffs = A.fullPivLu().solve(v);

    polynomial wp;

    for (int i = 0; i < N; ++i) {
        wp.push_back(wcoeffs[N - i - 1]);
    }
    wp.push_back(1);

    cout << "Ортогональный многочлен: w(x) = " << wp << endl << endl;

    roots = find_roots(wp, a - h, b + h, e, N * 2);

    vector<double> tmp3;
    vector<double> tmp4;

    for (int i = 0; i < N; ++i) {
        tmp3.push_back((double) moments[i]);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            tmp4.push_back((double) pow(roots[j], i));
        }
    }

    Map<VectorXd> u(&tmp3[0], N);
    Map<MatrixXd> B(&tmp4[0], N, N);

    B.transposeInPlace();

    VectorXd coeffs(N);

    coeffs = B.fullPivLu().solve(u);

    cout << "Узлы  и коэффициенты КФНАСТ:" << showpos << endl;

    for (int i = 0; i < N; ++i) {
        cout << setw(22) << roots[i] << " ---> " << setw(22) << coeffs[i] << endl;
    }

    cout << noshowpos << endl;

    long double test_integral = 0;
    const polynomial test_p = build_polynomial(2 * N - 1);
    for (int i = 0; i < N; ++i) {
        test_integral += coeffs[i] * (test_p ^ roots[i]);
    }

    cout << "Значение КФНАСТ для многочлена " << test_p << ":" << endl
         << test_integral << endl
         << "Погрешность для многочлена: " << endl
         << abs(test_integral - polynomial_integral(N)) << endl;

    integral = 0;
    for (int i = 0; i < N; ++i) {
        integral += coeffs[i] * f(roots[i]);
    }

    cout << "Значение КФНАСТ для функции:" << endl << integral << endl
         << "Погрешность для отрезка [0, 1]: " << abs(integral - INT);

    return 0;
}
