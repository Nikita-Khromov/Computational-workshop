#ifndef TASK6_POLYNOMIAL_H
#define TASK6_POLYNOMIAL_H
#include <vector>
#include <list>
#include <ostream>

typedef std::list<long double> polynomial;     //deg = size, beginning deg = 0
typedef std::pair<long double, int> monom;    //first - coeff, second - deg

polynomial mono_to_poli(const monom& m);

polynomial operator*(long double k, const polynomial& p);

monom operator*(long double k, const monom& m);

monom operator*(const monom& m1, const monom& m2);

polynomial operator*(const polynomial& p, const monom& m);

polynomial operator+(const polynomial& p1, const polynomial& p2);

polynomial operator+(const monom& m1, const monom& m2);

polynomial operator-(const polynomial& p1, const polynomial& p2);

polynomial operator-(const monom& m1, const monom& m2);

polynomial operator+(const polynomial& p, const monom m);

polynomial operator-(const polynomial& p, const monom m);

std::ostream& operator<<(std::ostream& stream, const polynomial& p);

long double operator^(const polynomial& p, const long double& x);

polynomial next_Legendre(const std::vector<polynomial>& Legendre);

void fill_Legendre(std::vector<polynomial>& Legendre, size_t n);

std::vector<long double> find_roots(const polynomial& p, long double a, long double b, long double e, const size_t N);

std::vector<long double> Gauss_coefficients(const std::vector<long double>& nodes, const std::vector<polynomial>& Legendre, size_t N);

#endif
