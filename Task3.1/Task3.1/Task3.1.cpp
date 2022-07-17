#include <iostream>
#include <locale>
#include <vector>
#include <iomanip>
#include <algorithm>

double f(double z)
{
	return exp(-z) - z * z / 2;
}

//double f(double z)
//{
//	return z * z;
//}

bool cmp(const std::pair<double, std::pair<double, double>>& a, const std::pair<double, std::pair<double, double>>& b)
{
	return (a.second.second < b.second.second);
}

double prod_frac(const double& F, int k, const int& n, const std::vector<std::pair<double, std::pair<double, double>>>& points)
{
	double prod1 = 1, prod2 = 1, Fk = points[k].second.first;
	for (int i = 0; i < n; ++i)
	{
		if (i != k)
		{
			prod1 *= F - points[i].second.first;
			prod2 *= Fk - points[i].second.first;
		}
	}
	return prod1 / prod2;
}

double lagrange(const double& F, const int& n, const std::vector<std::pair<double, std::pair<double, double>>>& points)
{
	double l = 0;
	for (int i = 0; i < n; ++i)
		l += prod_frac(F, i, n, points) * points[i].first;
	return l;
}

bool cmp1(const std::pair<std::pair<double, double>, double>& a, const std::pair<std::pair<double, double>, double>& b)
{
	return (a.first.second < b.first.second);
}

double prod_frac1(const double& x, int k, const int& n, const std::vector<std::pair<std::pair<double, double>, double>>& points)
{
	double prod1 = 1, prod2 = 1, xk = points[k].first.first;
	for (int i = 0; i < n; ++i)
	{
		if (i != k)
		{
			prod1 *= x - points[i].first.first;
			prod2 *= xk - points[i].first.first;
		}
	}
	return prod1 / prod2;
}

double lagrange1(const double& x, const int& n, const std::vector<std::pair<std::pair<double, double>, double>>& points)
{
	double l = 0;
	for (int i = 0; i < n; ++i)
		l += prod_frac1(x, i, n, points) * points[i].second;
	return l;
}

using namespace std;

double PnF(const int& m, const int& n, const double& a, const double& b, const double& x, const double& F)
{
	vector<pair<pair<double, double>, double>> interpol;
	double z, h;
	z = a; h = (b - a) / ((double)m - 1);
	for (int i = 0; i < m; ++i)
	{
		interpol.push_back(make_pair(make_pair(z, 0), f(z)));
		z += h;
	}
	for (auto&& i : interpol)
	{
		i.first.second = abs(i.first.first - x);
	}
	sort(interpol.begin(), interpol.end(), cmp1);
	return lagrange1(x, n, interpol) - F;
}


int main()
{
    setlocale(0, "RUS");
    int m, n, N = 10;
	double a, b, z, h, F, result, eps, ai, bi, c, y;
    bool cmnd;
	vector<double> an;
    vector<pair<double, pair<double, double>>> interpol;
	cout << "Задача обратного интерполирования" << endl << "Вариант 21" << endl
		<< "функция f(x) = e^(-x) - (x^2)/2" << endl;

	cout << "Левый конец отрезка a = "; cin >> a;
	cout << "Правый конец отрезка b = "; cin >> b;
	cout << "Количество пар точка-значение в таблице m+1 = "; cin >> m;

	z = a; h = (b - a) / ((double)m - 1);
	for (int i = 0; i < m; ++i)
	{
		interpol.push_back(make_pair(z, make_pair(f(z), 0)));
		z += h;
	}
	cout << endl;
	cout << "Сгенерированная таблица:" << endl;
	for (auto&& i : interpol)
		cout << setw(16) << setprecision(14) << left << i.first << " " << i.second.first << endl;
	cout << endl;
	while (true)
	{
		cout << "Значение функции для нахождения аргумента F = "; cin >> F;
		do
		{
			cout << "Количество точек для интерполирования n (< " << m << ") = "; cin >> n;
			if (n >= m) cout << "Неправда, " << n << " >= " << m << endl;
		} while (n >= m);
		cout << "Точность для второго способа e = "; cin >> eps;
		cout <<endl << "1 способ:" << endl << endl;
		for (auto&& i : interpol)
		{
			i.second.second = abs(i.second.first - F);
		}
		sort(interpol.begin(), interpol.end(), cmp);
		cout << endl << "Отсортированный массив: " << endl;
		for (int i = 0; i < n; ++i)
			cout << setw(16) << setprecision(14) << left << interpol[i].first << " " << interpol[i].second.first << endl;
		cout << endl;

		result = lagrange(F, n, interpol);
		cout << setprecision(16) << "Значение интерполяционного многочлена в точке F (искомый x): " << result << endl
			<< "Модуль невязки |f(x) - F| = " << abs(f(result) - F) << endl << endl;

		cout << "2 способ:" << endl << endl;

		h = (b - a) / N;
		ai = a; bi = a + h;
		while (bi <= b)
		{
			if (PnF(m, n, a, b, ai, F) * PnF(m, n, a, b, bi, F) <= 0)
			{
				an.push_back(ai);
			}
			ai = bi; bi += h;
		}
		cout << "Корни уравнения Pn(x) = F:" << endl << endl;
		for (auto it : an)
		{
			bi = it + h;
			c =it; y = bi;
			result = y - (PnF(m, n, a, b, y, F) * (y - c)) / (PnF(m, n, a, b, y, F) - PnF(m, n, a, b, c, F));
			while (abs(result - y) >= eps)
			{
				c = y;
				y = result;
				result = y - (PnF(m, n, a, b, y, F) * (y - c)) / (PnF(m, n, a, b, y, F) - PnF(m, n, a, b, c, F));
			}
			cout << result << "  :  Невязка |f(x) - F| = " << abs(f(result) - F) << endl;
		}

		cout << "Ввести новые начальные данные? (0 для выхода, 1 для ввода) "; cin >> cmnd;
		if (!cmnd) break;
	}
	return 0;
}
