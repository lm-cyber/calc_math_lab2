#include "fixed_point_iteration_metod.h"
#include "util.h"

double optimal_func(double x, double (*f)(double)) {
    return x - (1 / derivative(1, x, f)) * f(x);
}
double relaxation_method(double x,
                         double (*f)(double)) { // -2 < C*f'(x) < 0
                                                // -2/f'(x) < C < 0
                                                // c = -1/f'(x)
    return -1 / derivative(1, x, f);
}
double fi_func(double x, double (*f)(double)) {
    return x + relaxation_method(x, f) * f(x);
}
x_and_f_x fixed_point_iteration_metod(
    double x, double (*f)(double), double eps,
    std::vector<data_table_for_fixed_point_iteration> &data_tab,
    size_t max_iter) {
    double x_exit = fi_func(x, f);
    data_tab.push_back({x, f(x), x_exit, fi_func(x, f), delta(x, x_exit)});
    size_t currend_iterating = 0;
    while (delta(x, x_exit) > eps && currend_iterating < max_iter) {
        x = x_exit;
        x_exit = fi_func(x_exit, f);
        data_tab.push_back({x, f(x), x_exit, fi_func(x, f), delta(x, x_exit)});
        currend_iterating++;
    }

    return {x, f(x)};
}
