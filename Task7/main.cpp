#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

long double f(long double x, long double y) {
    return cosl(x) - y;
}

long double taylor(long double x) {
    return 1 - (powl(x, 3) / 6) + (powl(x, 4) / 24) - (powl(x, 7) / 5040) + (powl(x, 8) / 40320);
}

long double solution(long double x) {
    return 0.5 * (expl(-x) + sinl(x) + cosl(x));
}

using namespace std;

int main() {
    long double h, x;
    vector<long double> taylor_solution, euler1, euler2, euler3, k(4), runge_kutta;
    const long double x0 = 0, y0 = 1;
    int N;

    cout << setprecision(16) << left;
    cout << "Задача Коши: y'(x) = -y(x) + cos(x)" << endl << endl;
    cout << "Введите параметры задачи:" << endl << "h = ";
    cin >> h;
    cout << "N = ";
    cin >> N;
//    vector<vector<long double>> adams(N + 3, vector<long double>(6, 0));
    vector<vector<long double>> adams(N + 3);

    cout << endl << "Таблица последних значений точного решения: " << endl
         << setw(12) << "x_k" << "y(x_k)" << endl;
    for (int i = -2; i <= N; ++i) {
        x = x0 + i * h;
        if (N - i <= 15 || i <= 5){
            cout << setw(6) << x << " ---> " << setw(22) << solution(x) << endl;
        }
    }
    cout << endl;

    cout << "Таблица последних значений решения методом разложения в ряд Тейлора и абсолютные погрешности: " << endl
         << setw(12) << "x_k" << setw(24) << "y_k" << "Абсолютная погрешность" << endl;
    for (int i = -2; i <= N; ++i) {
        x = x0 + i * h;
        taylor_solution.push_back(taylor(x));
        if (N - i <= 15 || i <= 5) {
            cout << setw(6) << x << " ---> " << setw(18) << taylor_solution[i + 2] << " ---> "
                 << abs(taylor_solution[i + 2] - solution(x)) << endl;
        }
    }

    cout << endl << "Таблица последних значений решений методами Эйлера:" << endl
         << setw(12) << "x_k" << setw(29) << "I метод" << setw(29) << "II метод" << "III метод" << endl;
    euler1.push_back(y0);
    euler2.push_back(y0);
    euler3.push_back(y0);
    runge_kutta.push_back(y0);
    for (int i = 1; i <= N; ++i) {
        x = x0 + h * (i - 1);
        euler1.push_back(euler1[i - 1] + h * f(x, euler1[i - 1]));
        euler2.push_back(euler2[i - 1] + h * f(x + h / 2, euler2[i - 1] + h / 2 * f(x, euler2[i - 1])));
        euler3.push_back(
                euler3[i - 1] + (h / 2) * (f(x, euler3[i - 1]) + f(x + h, euler3[i - 1] + h * f(x, euler3[i - 1]))));
        if (N - i <= 15 || i <= 5) {
            cout << setw(6) << x + h << " ---> " << setw(18) << euler1[i] << " ---> " << setw(18) << euler2[i]
                 << " ---> "
                 << setw(18) << euler3[i] << endl;
        }
    }

    cout << endl << "Абсолютные погрешности решений методами Эйлера в точке x = " << x + h << ":" << endl
         << "I метод:   " << abs(euler1[N] - solution(x + h)) << endl
         << "II метод:  " << abs(euler2[N] - solution(x + h)) << endl
         << "III метод: " << abs(euler3[N] - solution(x + h)) << endl
         << endl;

    cout << "Таблица последних значений решения методом Рунге-Кутты 4-го порядка:" << endl
         << setw(12) << "x_k" << setw(18) << "y_k" << endl;
    for (int i = 1; i <= N; ++i) {
        x = x0 + h * (i - 1);
        k[0] = h * f(x, runge_kutta[i - 1]);
        k[1] = h * f(x + h / 2, runge_kutta[i - 1] + k[0] / 2);
        k[2] = h * f(x + h / 2, runge_kutta[i - 1] + k[1] / 2);
        k[3] = h * f(x + h, runge_kutta[i - 1] + k[2]);
        runge_kutta.push_back(runge_kutta[i - 1] + (k[0] + 2 * k[1] + 2 * k[2] + k[3]) / 6);
        if (N - i <= 15 || i <= 5) {
            cout << setw(6) << x + h << " ---> " << setw(18) << runge_kutta[i] << endl;
        }
    }

    cout << endl << "Абсолютная погрешность решения методом Рунге-Кутты в точке x = " << x + h << ":" << endl
         << abs(runge_kutta[N] - solution(x + h)) << endl << endl;

//    adams[0][0] = taylor_solution[0];
    adams[0].push_back(taylor_solution[0]);
//    adams[1][0] = taylor_solution[1];
    adams[1].push_back(taylor_solution[1]);
//    adams[2][0] = taylor_solution[2];
    adams[2].push_back(taylor_solution[2]);
//    adams[3][0] = taylor_solution[3];
    adams[3].push_back(taylor_solution[3]);
//    adams[4][0] = taylor_solution[4];
    adams[4].push_back(taylor_solution[4]);

    for (int i = 1; i < 6; ++i) {
        for (int j = 0; j <= 4 - i; ++j) {
            x = x0 + h * (j - 2);
            if (i == 1)
                adams[j].push_back(h * f(x, adams[j][0]));
//                adams[j][i] = h * f(x, adams[j][0]);
            else
                adams[j].push_back(-(adams[j][i - 1] - adams[j + 1][i - 1]));
//                adams[j][i] = -(adams[j][i - 1] - adams[j + 1][i - 1]);
        }
    }

    cout << "Таблица последних значений решения экстраполяционным методом Адамса 4-го порядка:" << endl
         << setw(12) << "x_k" << setw(18) << "y_k" << endl;
    for (int i = 5; i <= N + 2; ++i) {
        x = x0 + h * (i - 3);
//        adams[i - 1][1] = h * f(x, adams[i - 1][0]);
        adams[i - 1].push_back(h * f(x, adams[i - 1][0]));
//        adams[i - 2][2] = -(adams[i - 2][1] - adams[i - 1][1]);
        adams[i - 2].push_back(-(adams[i - 2][1] - adams[i - 1][1]));
//        adams[i - 3][3] = -(adams[i - 3][2] - adams[i - 2][2]);
        adams[i - 3].push_back(-(adams[i - 3][2] - adams[i - 2][2]));
//        adams[i - 4][4] = -(adams[i - 4][3] - adams[i - 3][3]);
        adams[i - 4].push_back(-(adams[i - 4][3] - adams[i - 3][3]));
//        adams[i - 5][5] = -(adams[i - 5][4] - adams[i - 4][4]);
        adams[i - 5].push_back(-(adams[i - 5][4] - adams[i - 4][4]));
//        adams[i][0] = adams[i - 1][0] + adams[i - 1][1] + adams[i - 2][2] / 2 + 5 * adams[i - 3][3] / 12 +
        adams[i].push_back(adams[i - 1][0] + adams[i - 1][1] + adams[i - 2][2] / 2 + 5 * adams[i - 3][3] / 12 +
                      3 * adams[i - 4][4] / 8 + 251 * adams[i - 5][5] / 720);
        if (N - i <= 15 || i <= 5) {
            cout << setw(6) << x + h << " ---> " << setw(18) << adams[i][0] << endl;
        }
    }

//    for (int i = 0; i <= N+2; ++i) {
//        for (int j = 0; j < 6; ++j) {
//            cout << setw(22) << adams[i][j] << " ";
//        }
//        cout << endl;
//    }

    cout << endl << "Абсолютная погрешность решения ЭМА в точке x = " << x + h << ":" << endl
         << abs(adams[N + 2][0] - solution(x + h));

    return 0;
}
