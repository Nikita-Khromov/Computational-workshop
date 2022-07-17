#include <iostream>
#include <cmath>
#include <locale>
#include <string>
#include <vector>
#include <iomanip>

long double f(long double z, int i)
{
	long double res;
	switch (i)
	{
	case 0: 
		return (2 * exp(-2 * z) * z) / (1 + z * z) + (16 * z) / (1 + z * z * z * z) - 2 * exp(-2 * z) * log(1 + z * z) + 1;
	case 1:
		return 5;
	case 2:
		return 3 * z + 1;
	case 3:
		return 8 * z * z + 3 * z - 7;
	case 4:
		return z * z * z - 6 * z * z + z - 12;
	default:
		break;
	}
}

long double f_int(long double z, int i)
{
	switch (i)
	{
	case 0:
		return 8 * atan(z * z) + exp(-2 * z) * log(z * z + 1) + z;
	case 1:
		return 5 * z;
	case 2:
		return 3 * z * z / 2 + z;
	case 3:
		return 8 * z * z * z / 3 + 3 * z * z / 2 - 7 * z;
	case 4:
		return z * z * z * z / 4 - 2 * z * z * z + z * z / 2 - 12 * z;
	default:
		break;
	}
}

using namespace std;

int main()
{
	setlocale(0, "RUS");
	long double a, b, integral, real_integral;
	int cmnd;
	vector <string> func;
	func.push_back("f(x) = (2x*e^(-2x))/(1 + x^2) + 16x/(1 + x^4) - 2e^(-2x)*ln(1 + x^2) + 1");
	func.push_back("f(x) = 5");
	func.push_back("f(x) = 3x + 1");
	func.push_back("f(x) = 8x^2 + 3x - 7");
	func.push_back("f(x) = x^3 - 6x^2 + x - 12");
	cout << "Приближённое вычисление интеграла по квадратурным формулам" << endl << "Функции:" << endl
		<< "1) f(x) = (2x*e^(-2x))/(1 + x^2) + 16x/(1 + x^4) - 2e^(-2x)*ln(1 + x^2) + 1" << endl
		<< "2) f(x) = 5" << endl << "3) f(x) = 3x + 1" << endl << "4) f(x) = 8x^2 + 3x - 7" << endl
		<< "5) f(x) = x^3 - 6x^2 + x - 12" << endl << endl;

	while (true)
	{
		cout << "Введите левый конец отрезка интегрирования: a = "; cin >> a;
		cout << "Введите правый конец отрезка интегрирования: b = "; cin >> b;
		cout << endl;
		for (int i = 0; i < 5; ++i)
		{
			cout << "Функция " << i + 1 << ") " << func[i] << endl << endl;
			real_integral = f_int(b, i) - f_int(a, i);
			cout << setprecision(16) << "Фактическое значение интеграла: " << real_integral << endl;
			integral = (b - a) * f(a, i);
			cout << "Значение интеграла, полученное с помощью КФ левого прямоугольника: " << integral << endl
				<< "Фактическая погрешность: " << abs(integral - real_integral) << endl;
			integral = (b - a) * f(b, i);
			cout << "Значение интеграла, полученное с помощью КФ правого прямоугольника: " << integral << endl
				<< "Фактическая погрешность: " << abs(integral - real_integral) << endl;
			integral = (b - a) * f((a + b) / 2, i);
			cout << "Значение интеграла, полученное с помощью КФ среднего прямоугольника: " << integral << endl
				<< "Фактическая погрешность: " << abs(integral - real_integral) << endl;
			integral = (b - a) / 2 * (f(a, i) + f(b, i));
			cout << "Значение интеграла, полученное с помощью КФ трапеции: " << integral << endl
				<< "Фактическая погрешность: " << abs(integral - real_integral) << endl;
			integral = (b - a) / 6 * (f(a, i) + f(b, i) + 4 * f((a + b) / 2, i));
			cout << "Значение интеграла, полученное с помощью КФ Симпсона: " << integral << endl
				<< "Фактическая погрешность: " << abs(integral - real_integral) << endl;
			integral = (b - a) * (f(a, i) / 8 + f(b, i) / 8 + 3 * f(a + (b - a) / 3, i) / 8 + 3 * f(a + 2 * (b - a) / 3, i) / 8);
			cout << "Значение интеграла, полученное с помощью КФ 3/8: " << integral << endl
				<< "Фактическая погрешность: " << abs(integral - real_integral) << endl;
			
			cout << endl << endl;
		}
		cout << endl << "Ввести новый отрезок? (1 для ввода, 0 для прекращения работы) ";
		cin >> cmnd;
		if (!cmnd)
			break;
	}
	return 0;
}
