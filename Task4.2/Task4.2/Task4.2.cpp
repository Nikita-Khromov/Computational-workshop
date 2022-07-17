#include <iostream>
#include <locale>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>

long double f(long double z, int i)
{
	long double res;
	switch (i)
	{
	case 0:
		//return (2 * exp(-2 * z) * z) / (1 + z * z) + (16 * z) / (1 + z * z * z * z) - 2 * exp(-2 * z) * log(1 + z * z) + 1;
		return exp(z);
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
		//return 8 * atan(z * z) + exp(-2 * z) * log(z * z + 1) + z;
		return exp(z);
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

long double w_sum(long double a, long double b, long double h, int i)
{
	long double res = 0, z = a + h;
	while (b - z >= h / 2)
	{
		res += f(z, i);
		z += h;
	}
	return res;
}

long double q_sum(long double a, long double b, long double h, int i)
{
	long double res = 0, z = a + h / 2;
	while (b - z >= h / 4)
	{
		res += f(z, i);
		z += h;
	}
	return res;
}

//long double left_rect(long double a, long double b, int i)
//{
//	return (b - a) * f(a, i);
//}
//
//long double right_rect(long double a, long double b, int i)
//{
//	return (b - a) * f(b, i);
//}
//
//long double center_rect(long double a, long double b, int i)
//{
//	return (b - a) * f((a + b) / 2, i);
//}
//
//long double trapezoid(long double a, long double b, int i)
//{
//	return (b - a) / 2 * (f(a, i) + f(b, i));
//}
//
//long double simpson(long double a, long double b, int i)
//{
//	return ((b - a) / 6) * (f(a, i) + f(b, i) + 4 * f((a + b) / 2, i));
//}

using namespace std;

int main()
{
	setlocale(0, "RUS");
	long double a, b, integral, integral_pp, integral_runge, real_integral, h, z, q, w, theor_r, m, q_pp, w_pp;
	int cmnd, l;
	vector <string> func;
	func.push_back("f(x) = (2x*e^(-2x))/(1 + x^2) + 16x/(1 + x^4) - 2e^(-2x)*ln(1 + x^2) + 1");
	func.push_back("f(x) = 5");
	func.push_back("f(x) = 3x + 1");
	func.push_back("f(x) = 8x^2 + 3x - 7");
	func.push_back("f(x) = x^3 - 6x^2 + x - 12");
	cout << "Приближённое вычисление интеграла по составным квадратурным формулам, уточнение значение интеграла с помощью формулы Рунге" << endl << "Функции:" << endl
		<< "1) f(x) = (2x*e^(-2x))/(1 + x^2) + 16x/(1 + x^4) - 2e^(-2x)*ln(1 + x^2) + 1" << endl
		<< "2) f(x) = 5" << endl << "3) f(x) = 3x + 1" << endl << "4) f(x) = 8x^2 + 3x - 7" << endl
		<< "5) f(x) = x^3 - 6x^2 + x - 12" << endl << endl;

	while (true)
	{
		cout << "Введите левый конец отрезка интегрирования: a = "; cin >> a;
		cout << "Введите правый конец отрезка интегрирования: b = "; cin >> b;
		cout << "Введите количество промежутков деления отрезка m = "; cin >> m;
		cout << "Введите количество множитель количества делений отрезка l = "; cin >> l;
		m = round(m);
		cout << endl;

		const vector <vector<long double>> Maxs = { {18, 48.4742, 1048.35}, {0, 0, 0}, {3, 0, 0}, {max(abs(16 * a + 3), abs(16 * b + 3)), 16, 0}, 
			{max(abs(3 * b * b - 12 * b + 1), max(abs(3 * a * a - 12 * a + 1), (long double) 11)), max(abs(6 * a - 12), abs(6 * b - 12)), 0} };
		h = (b - a) / m;
		cout << "Длина одного деления h = " << h << endl << endl;

		for (int i = 0; i < 5; ++i)
		{
			cout << "Функция " << i + 1 << ") " << func[i] << endl << endl;
			real_integral = f_int(b, i) - f_int(a, i);
			z = f(a, i) + f(b, i);
			cout << setprecision(16) << "Фактическое значение интеграла: " << real_integral << endl;
			q = q_sum(a, b, h, i); w = w_sum(a, b, h, i);
			q_pp = q_sum(a, b, h / l, i); w_pp = w_sum(a, b, h / l, i);
			integral = h * (f(a, i) + w);
			integral_pp = (h / l) * (f(a, i) + w_pp);
			integral_runge = (l * integral_pp - integral) / (l - (double)1);
			theor_r = (0.5) * Maxs[i][0] * (b - a) * h;
			cout << "Значение интеграла, полученное с помощью СКФ левого прямоугольника при первом разбиении: " << integral << endl
				<< "Значение интеграла, полученное с помощью СКФ левого прямоугольника при увеличенном разбиении: " << integral_pp << endl
				<< "Значение интеграла, полученное с помощью уточнения по принципу Рунге: " << integral_runge << endl
				<< "Фактическая погрешность для первого разбиения: " << abs(integral - real_integral) << endl
				<< "Фактическая погрешность для увеличенного разбиения: " << abs(integral_pp - real_integral) << endl
				<< "Фактическая погрешность для уточнения по Рунге: " << abs(integral_runge - real_integral) << endl
				<< "Теоретическая погрешность";
			if (i == 0) cout << " (верно при a = 0)";
			cout << ": " << theor_r << endl << endl;
			integral = h * (w + f(b, i));
			integral_pp = (h / l) * (w_pp + f(b, i));
			integral_runge = (l * integral_pp - integral) / (l - (double)1);
			cout << "Значение интеграла, полученное с помощью СКФ правого прямоугольника при первом разбиении: " << integral << endl
				<< "Значение интеграла, полученное с помощью СКФ правого прямоугольника при увеличенном разбиении: " << integral_pp << endl
				<< "Значение интеграла, полученное с помощью уточнения по принципу Рунге: " << integral_runge << endl
				<< "Фактическая погрешность для первого разбиения: " << abs(integral - real_integral) << endl
				<< "Фактическая погрешность для увеличенного разбиения: " << abs(integral_pp - real_integral) << endl
				<< "Фактическая погрешность для уточнения по Рунге: " << abs(integral_runge - real_integral) << endl
				<< "Теоретическая погрешность";
			if (i == 0) cout << " (верно при a = 0)";
			cout << ": " << theor_r << endl << endl;
			integral = h * q;
			integral_pp = (h / l) * q_pp;
			integral_runge = (pow(l, 2) * integral_pp - integral) / (pow (l, 2) - (double)1);
			theor_r = (1.0 / 24) * Maxs[i][1] * (b - a) * h * h;
			cout << "Значение интеграла, полученное с помощью СКФ среднего прямоугольника при первом разбиении: " << integral << endl
				<< "Значение интеграла, полученное с помощью СКФ среднего прямоугольника при увеличенном разбиении: " << integral_pp << endl
				<< "Значение интеграла, полученное с помощью уточнения по принципу Рунге: " << integral_runge << endl
				<< "Фактическая погрешность для первого разбиения: " << abs(integral - real_integral) << endl
				<< "Фактическая погрешность для увеличенного разбиения: " << abs(integral_pp - real_integral) << endl
				<< "Фактическая погрешность для уточнения по Рунге: " << abs(integral_runge - real_integral) << endl
				<< "Теоретическая погрешность";
			if (i == 0) cout << " (верно при a < 0.7 и b > 0.7)";
			cout << ": " << theor_r << endl << endl;
			integral = (h / 2) * (z + 2 * w);
			integral_pp = (h / ((double) 2 * l)) * (z + 2 * w_pp);
			integral_runge = (pow(l, 2) * integral_pp - integral) / (pow(l, 2) - (double)1);
			cout << "Значение интеграла, полученное с помощью СКФ трапеции при первом разбиении: " << integral << endl
				<< "Значение интеграла, полученное с помощью СКФ трапеции при увеличенном разбиении: " << integral_pp << endl
				<< "Значение интеграла, полученное с помощью уточнения по принципу Рунге: " << integral_runge << endl
				<< "Фактическая погрешность для первого разбиения: " << abs(integral - real_integral) << endl
				<< "Фактическая погрешность для увеличенного разбиения: " << abs(integral_pp - real_integral) << endl
				<< "Фактическая погрешность для уточнения по Рунге: " << abs(integral_runge - real_integral) << endl
				<< "Теоретическая погрешность";
			if (i == 0) cout << " (верно при a < 0.7 и b > 0.7)";
			cout << ": " << 2 * theor_r << endl << endl;
			integral = (h / 6) * (z + 2 * w + 4 * q);
			integral_pp = (h / ((double)6 * l)) * (z + 2 * w_pp + 4 * q_pp);
			integral_runge = (pow(l, 4) * integral_pp - integral) / (pow(l, 4) - (double)1);
			theor_r = (1.0 / 2880) * Maxs[i][2] * (b - a) * h * h * h * h;
			cout << "Значение интеграла, полученное с помощью СКФ Симпсона при первом разбиении: " << integral << endl
				<< "Значение интеграла, полученное с помощью СКФ Симпсона при увеличенном разбиении: " << integral_pp << endl
				<< "Значение интеграла, полученное с помощью уточнения по принципу Рунге: " << integral_runge << endl
				<< "Фактическая погрешность для первого разбиения: " << abs(integral - real_integral) << endl
				<< "Фактическая погрешность для увеличенного разбиения: " << abs(integral_pp - real_integral) << endl
				<< "Фактическая погрешность для уточнения по Рунге: " << abs(integral_runge - real_integral) << endl
				<< "Теоретическая погрешность";
			if (i == 0) cout << " (верно при a < 0.7 и b > 0.7)";
			cout << ": " << theor_r << endl;

			cout << endl << endl;
		}
		cout << endl << "Ввести новый отрезок? (1 для ввода, 0 для прекращения работы) ";
		cin >> cmnd;
		if (!cmnd)
			break;
	}

	return 0;
}