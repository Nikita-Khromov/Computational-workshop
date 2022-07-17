#include <iostream>
#include <vector>
#include <algorithm>
#include <locale>
#include <iomanip>

double f(double z)
{
	return exp(-z) - z * z / 2;
}

bool cmp(const std::pair<std::pair<double, double>, double> &a, const std::pair<std::pair<double, double>, double> &b)
{
	return (a.first.second < b.first.second);
}

double prod_frac(const double &x, int k, const int &n, const std::vector<std::pair<std::pair<double, double>, double>> &points)
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

double lagrange(const double &x, const int &n, const std::vector<std::pair<std::pair<double, double>, double>> &points)
{
	double l = 0;
	for (int i = 0; i < n; ++i)
		l += prod_frac(x, i, n, points) * points[i].second;
	return l;
}

void fill_diffs(const int &n, std::vector<std::vector<double>> &diff, const std::vector<std::pair<std::pair<double, double>, double>> &points)
{
	for (int j = 0; j < n; ++j)
		diff[j][0] = points[j].second;
	for (int i = 1; i < n; ++i)
	{
		for (int j = 0; j < n - i; ++j)
		{
			diff[j][i] = (diff[j + 1][i - 1] - diff[j][i - 1]) / (points[i + j].first.first - points[j].first.first);
		}
	}
}

double prod_newt(const double& x, const std::vector<std::pair<std::pair<double, double>, double>>& points, const int &k)
{
	double res = 1;
	for (int i = 0; i < k; ++i)
		res *= x - points[i].first.first;
	return res;
}

double newton(const double &x, std::vector<std::vector<double>>& diff, const std::vector<std::pair<std::pair<double, double>, double>>& points)
{
	double n = points[0].second;
	for (int i = 1; i < diff.size(); ++i)
	{
		n += diff[0][i] * prod_newt(x, points, i);
	}
	return n;
}

using namespace std;

int main()
{
	setlocale(0, "RUS");
	int m, n;
	bool cmnd;
	double a, b, x, z, h, result;
	vector<pair<pair<double, double>, double>> interpol;
	cout << "Алгебраическое интерполирование методами Ньютона и Лагранжа" << endl << "Вариант 7" << endl
		<< "функция f(x) = e^(-x) - (x^2)/2" << endl;

	cout << "Левый конец отрезка a = "; cin >> a;
	cout << "Правый конец отрезка b = "; cin >> b;
	cout << "Количество пар точка-значение в таблице m+1 = "; cin >> m;
	
	z = a; h = (b - a) / ((double)m - 1);
	for (int i = 0; i < m; ++i)
	{
		//points.push_back(z);
		interpol.push_back(make_pair(make_pair(z, 0), f(z)));
		z += h;
	}
	cout << endl;
	cout << "Сгенерированная таблица:" << endl;
	for (auto i : interpol)
		cout << setw(16) << setprecision(14) << left << i.first.first << " " << i.second << endl;
	cout << endl;
	while (true)
	{
		cout << "Точка для нахождения значения x = "; cin >> x;
		do 
		{
			cout << "Количество точек для интерполирования n (< " << m << ") = "; cin >> n;
			if (n >= m) cout << "Неправда, " << n << " >= " << m << endl;
		} while (n >= m);

		vector<vector<double>> div_dif(n, vector<double>(n, 0));

		for (auto &&i : interpol)
		{
			i.first.second = abs(i.first.first - x);
		}

		sort(interpol.begin(), interpol.end(), cmp);
		cout << endl << "Отсортированный массив: " << endl;
		for (int i = 0; i < n; ++i)
			cout << setw(16) << setprecision(14) << left << interpol[i].first.first << " " << interpol[i].second << endl;
		cout << endl;

		result = lagrange(x, n, interpol);
		cout << setprecision(16) << "Метод Лагранжа" << endl << "Значение интерполяционного многочлена в точке x: " << result << endl
			<< "Погрешность |f(x) - PL(x)| = " << abs(f(x) - result) << endl;
		fill_diffs(n, div_dif, interpol);
		result = newton(x, div_dif, interpol);
		cout << setprecision(16) << "Метод Ньютона" << endl << "Значение интерполяционного многочлена в точке x: " << result << endl
			<< "Погрешность |f(x) - PN(x)| = " << abs(f(x) - result) << endl;

		cout << "Ввести новые начальные данные? (0 для выхода, 1 для ввода) "; cin >> cmnd;
		if (!cmnd) break;
	}
	return 0;
}

