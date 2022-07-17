#include <iostream>
#include <locale>
#include <vector>
#include <iomanip>

double f(double z, bool flag)
{
	return !flag ? exp(-z) - z * z / 2 : exp(3 * z);
}

double fprime(double z, bool flag)
{
	return !flag ? -exp(-z) - z : 3 * exp(3 * z);
}

double fdprime(double z, bool flag)
{
	return !flag ? exp(-z) - 1 : 9 * exp(3 * z);
}

//double f(double z)
//{
//	return exp(3 * z);
//}
//
//double fprime(double z)
//{
//	return 3 * exp(3 * z);
//}
//
//double fdprime(double z)
//{
//	return 9 * exp(3 * z);
//}

using namespace std;

int main()
{
	setlocale(0, "RUS");
	int m, n, flag;
	double a, h, z;
	bool cmnd;
	vector <pair<double, double>> table;
	cout << "Задача численного дифференцирования" << endl << "Вариант 21" << endl
		<< "функции:" << endl << "1) f(x) = e ^ (-x) - (x ^ 2) / 2" << endl << "2) f(x) = e ^ (3x)" << endl;
	while (true)
	{

		//Подготовка
		vector<vector<double>> primes(8);
		cout << "Номер функции для исследования: "; cin >> flag; --flag;
		cout << "Левый конец отрезка a = "; cin >> a;
		cout << "Количество пар точка-значение в таблице m+1 = "; cin >> m;
		cout << "Расстояние между точками h = "; cin >> h;
		z = a;
		for (int i = 0; i < m; ++i)
		{
			table.push_back(make_pair(z, f(z, flag)));
			primes[0].push_back(z);
			primes[1].push_back(f(z, flag));
			z += h;
		}
		cout << endl << "Сгенерированная таблица:" << endl;
		for (auto&& i : table)
			cout << setw(16) << setprecision(14) << left << i.first << " " << i.second << endl;
		cout << endl;

		//Первая производная первый способ

		z = (-3 * table[0].second + 4 * table[1].second - table[2].second) / (2 * h);
		primes[2].push_back(z);
		primes[3].push_back(abs(fprime(table[0].first, flag) - z));
		for (int i = 1; i < m - 1; ++i)
		{
			z = (table[i + 1].second - table[i - 1].second) / (2 * h);
			primes[2].push_back(z);
			primes[3].push_back(abs(fprime(table[i].first, flag) - z));
		}
		z = (3 * table[m - 1].second - 4 * table[m - 2].second + table[m - 3].second) / (2 * h);
		primes[2].push_back(z);
		primes[3].push_back(abs(fprime(table[m - 1].first, flag) - z));

		//Первая производная второй способ
		//в формуле в знаменателе было 6h, но оказалось что ответ в 2 раза больше, чем надо, поэтому изменил на 12h

		z = (-25 * table[0].second + 48 * table[1].second - 36 * table[2].second + 16 * table[3].second - 3 * table[4].second) / (12 * h);
		primes[4].push_back(z);
		primes[5].push_back(abs(fprime(table[0].first, flag) - z));
		z = (-25 * table[1].second + 48 * table[2].second - 36 * table[3].second + 16 * table[4].second - 3 * table[5].second) / (12 * h);
		primes[4].push_back(z);
		primes[5].push_back(abs(fprime(table[1].first, flag) - z));
		for (int i = 2; i < m - 2; ++i)
		{
			z = (table[i - 2].second - 8 * table[i - 1].second + 8 * table[i + 1].second - table[i + 2].second) / (12 * h);
			primes[4].push_back(z);
			primes[5].push_back(abs(fprime(table[i].first, flag) - z));
		}
		z = (25 * table[m - 2].second - 48 * table[m - 3].second + 36 * table[m - 4].second - 16 * table[m - 5].second + 3 * table[m - 6].second) / (12 * h);
		primes[4].push_back(z);
		primes[5].push_back(abs(fprime(table[m - 2].first, flag) - z));
		z = (25 * table[m - 1].second - 48 * table[m - 2].second + 36 * table[m - 3].second - 16 * table[m - 4].second + 3 * table[m - 5].second) / (12 * h);
		primes[4].push_back(z);
		primes[5].push_back(abs(fprime(table[m - 1].first, flag) - z));

		//Вторая производная

		primes[6].push_back(NAN);
		primes[7].push_back(NAN);
		for (int i = 1; i < m - 1; ++i)
		{
			z = (table[i + 1].second - 2 * table[i].second + table[i - 1].second) / (h * h);
			primes[6].push_back(z);
			primes[7].push_back(abs(fdprime(table[i].first, flag) - z));
		}
		primes[6].push_back(NAN);
		primes[7].push_back(NAN);

		//Вывод

		cout << setw(22) << "| x" << setw(22) << "| f(x)" << setw(22) << "| f'чд(x) (1 способ)" << setw(22) << "| |f'(x) - f'чд(x)|" << setw(22) << "| f'чд(x) (2 способ)" << setw(22) << "| |f'(x) - f'чд(x)|" << setw(22) << "| f''чд(x)" << setw(22) << "| |f''(x) - f''чд(x)|  |" << endl;
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < 8; ++j)
			{
				cout << "| " << setw(20) << setprecision(12)  << primes[j][i] << "";
			}
			cout << " |" << endl;
		}
		cout << "Ввести новые начальные данные? (0 для выхода, 1 для ввода) "; cin >> cmnd;
		if (!cmnd) break;
		table.clear(); primes.clear();
	}
	return 0;
}
