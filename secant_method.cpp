#include "secant_method.h"
#include "util.h"

x_and_f_x secant_method(double (*f)(double), double a, double b, double eps,
                        std::vector<data_table_for_secant_metod> &data_tab,
                        size_t max_iter) {
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
    return {x, f(x)};
}
