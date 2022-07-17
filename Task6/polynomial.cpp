#include "polynomial.h"


polynomial mono_to_poli(const monom &m) {
    if (m.second == 0)
        return {m.first};
    polynomial res(m.second, 0);
    res.push_back(m.first);
    return res;
}

polynomial operator*(long double k, const polynomial &p) {
    polynomial res;
    if (k == 0)
        return res;
    for (auto i: p) {
        res.push_back(i * k);
    }
    return res;
}

monom operator*(long double k, const monom &m) {
    monom res = m;
    res.first *= k;
    if (k == 0)
        res.second = 0;
    return res;
}

monom operator*(const monom &m1, const monom &m2) {
    monom res;
    if (m1.first == 0 || m2.first == 0)
        return {0, 0};
    res.first = m1.first * m2.first;
    res.second = m1.second + m2.second;
    return res;
}

polynomial operator*(const polynomial &p, const monom &m) {
    polynomial res;
    if (m.second == 0)
        return m.first * p;
    if (p.empty())
        return res;
    for (int i = 0; i < m.second; ++i) {
        res.push_back(0);
    }
    for (auto i: p) {
        res.push_back(i * m.first);
    }
    return res;
}


polynomial operator+(const polynomial &p1, const polynomial &p2) {
    polynomial res;
    size_t n;
    n = (p1.size() > p2.size() ? p1.size() : p2.size());
    auto i1 = p1.cbegin();
    auto i2 = p2.cbegin();
    for (int i = 0; i < n; ++i) {
        if (i1 == p1.end()) {
            for (int j = i; j < n; ++j) {
                res.push_back(*i2);
                i2++;
            }
            break;
        }
        if (i2 == p2.end()) {
            for (int j = i; j < n; ++j) {
                res.push_back(*i1);
                i1++;
            }
            break;
        }
        res.push_back(*i1 + *i2);
        ++i1, ++i2;
    }
    if (res == polynomial(n, 0) || res.empty()) {
        res.clear();
        return res;
    }
    res.reverse();
    for (auto it = res.begin(); *it == 0; it = res.begin())
        res.erase(it);
    res.reverse();
    return res;
}

polynomial operator+(const monom &m1, const monom &m2) {
    polynomial p1, p2;
    p1 = mono_to_poli(m1);
    p2 = mono_to_poli(m2);
    return p1 + p2;
}

polynomial operator-(const polynomial &p1, const polynomial &p2) {
    polynomial temp;
    temp = -1 * p2;

    return (p1 + temp);
}

polynomial operator-(const monom &m1, const monom &m2) {
    polynomial p1, p2;
    p1 = mono_to_poli(m1);
    p2 = mono_to_poli(m2);
    return p1 - p2;
}

polynomial operator+(const polynomial &p, const monom m) {
    polynomial p2 = mono_to_poli(m);
    return p2 + p;
}

polynomial operator-(const polynomial &p, const monom m) {
    polynomial p2 = mono_to_poli(m);
    return p - p2;
}

std::ostream &operator<<(std::ostream &stream, const polynomial &p) {
    int deg = 0;
    bool fl = false;
    for (auto i: p) {
        if (i != 0) {
            if (deg == 0) {
                stream << i << " ";
                fl = true;
            } else if (deg == 1) {
                if (i < 0)
                    stream << "- ";
                else if (fl)
                    stream << "+ ";
                if (i != 1)
                    stream << std::abs(i) << "x ";
                else
                    stream << "x ";
                fl = true;
            } else {
                if (i < 0)
                    stream << "- ";
                else if (fl)
                    stream << "+ ";
                if (i != 1)
                    stream << std::abs(i) << "x^" << deg << " ";
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
long double operator^(const polynomial &p, const long double &x) {
    auto i = p.cend();
    --i;
    auto end = p.cbegin();
    long double res = 0;


    for (; i != end;) {
        res = x * (res + *i);
        i--;
    }
    res += *i;
    if (p.size() == 1)
        return *(p.begin());
    else
        return res;
}

polynomial next_Legendre(const std::vector<polynomial> &Legendre) {
    polynomial result_p;
    const size_t n = Legendre.size();
    monom x = std::make_pair(1, 1);
    result_p = (((long double) 2 * n - 1) / n) * Legendre[n - 1] * x - (((long double) n - 1) / n) * Legendre[n - 2];
    return result_p;
}

void fill_Legendre(std::vector<polynomial> &Legendre, size_t n) {
    for (int i = 2; i < n; ++i) {
        Legendre.push_back(next_Legendre(Legendre));
    }
}

std::vector<long double> find_roots(const polynomial &p, long double a, long double b, long double e, const size_t N) {
    long double h = (b - a) / (2 * N) - e;
    long double ai = a, bi = a + h;
    std::vector<long double> an, roots;
    while (bi <= b) {
        if ((p ^ ai) * (p ^ bi) <= 0) {
            an.push_back(ai);
        }
        ai = bi;
        bi += h;
    }

    for (auto cmnd: an) {
        ai = cmnd;
        bi = ai + h;
        long double c = ai, y = bi;
        long double x = y - ((p ^ y) * (y - c)) / ((p ^ y) - (p ^ c));
        while (std::abs(x - y) >= e) {
            c = y;
            y = x;
            x = y - ((p ^ y) * (y - c)) / ((p ^ y) - (p ^ c));
        }
        roots.push_back(x);
    }
    return roots;
}

std::vector<long double>
Gauss_coefficients(const std::vector<long double> &nodes, const std::vector<polynomial> &Legendre, size_t N) {
    std::vector<long double> coefficients;
    long double A;
    for (auto xk: nodes) {
        A = (2 * (1 - xk * xk)) / (N * N * (Legendre[N - 1] ^ xk) * (Legendre[N - 1] ^ xk));
        coefficients.push_back(A);
    }

    return coefficients;
}