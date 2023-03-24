#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>
typedef struct {
    double a;
    double b;
    double x;
    double f_a;
    double f_b;
    double f_x;
    double delta_eps;
} data_table_for_secant_metod;

double f1(double x) { return 0.4 * x * x * x + 3 * x * x + x - 13; }

double derivative(
    int dif, double point, double (*f)(double),
    double eps = 0.00001) { // нам важен знак, большой esp приводит к нулю
    if (dif == 0) {
        return f(point);
    } else {
        dif--;
        return (derivative(dif, point + eps, f) - derivative(dif, point, f)) /
               eps;
    }
}
double delta(double x_k, double x_k_minus_one) {
    return fabs(x_k - x_k_minus_one);
}

double secant_method(double (*f)(double), double a, double b, double eps,
                     std::vector<data_table_for_secant_metod> &data_tab,
                     size_t max_iter = 10000) {
    bool cheak_zero_derivative = false;
    bool left_rigth = false;
    double point;
    double x;
    if (f(b) * derivative(2, a, f) < 0) {
        point = a;
        x = b;
        left_rigth = true;
    } else if (f(a) * derivative(2, a, f) < 0) {
        point = b;
        x = a;
        left_rigth = false;
    } else {
        x = a - (b - a) / (f(b) - f(a)) * f(a);
        cheak_zero_derivative = true;
    }
    double x_exit = x + eps * 2;
    size_t currend_iterating = 0;

    data_tab.push_back({a, b, x, f(a), f(b), f(x), delta(x, x_exit)});

    while (currend_iterating < max_iter && delta(x, x_exit) > eps) {
        if (cheak_zero_derivative) {
            if (f(a) * f(b) < 0) {
                b = x;
            } else {
                a = x;
            }
            x_exit = x;
            x = a - (b - a) / (f(b) - f(a)) * f(a);
            data_tab.push_back({a, b, x, f(a), f(b), f(x), delta(x, x_exit)});
        } else {
            x_exit = x;
            x = x - (point - x) / (f(point) - f(x)) * f(x);
            if (left_rigth) {
                data_tab.push_back(
                    {point, x, x, f(point), f(x), f(x), delta(x, x_exit)});

            } else {
                data_tab.push_back(
                    {x, point, x, f(x), f(point), f(x), delta(x, x_exit)});
            }
        }
        currend_iterating++;
    }
    return x;
}

int main() {
    std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);

    std::vector<data_table_for_secant_metod> data_tab;
    double x = secant_method(f1, -4, -1, 0.000001, data_tab);
    std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);
    for (const auto &i : data_tab) {
        std::cout << i.a << ' ' << i.b << ' ' << i.x << ' ' << i.f_a << ' '
                  << i.f_b << ' ' << i.f_x << ' ' << i.delta_eps << '\n';
    }
    std::cout << "\n\n" << x;
}
