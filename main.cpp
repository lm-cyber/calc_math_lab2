#include "secant_method.h"
#include "util.h"
// #include <boost/numeric/ublas/tensor.hpp>
#include "fixed_point_iteration_metod.h"
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

double f1(double x) { return 0.4 * x * x * x + 3 * x * x + x - 13; }

int main() {
    std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);

    std::vector<data_table_for_secant_metod> data_tab;
    x_and_f_x x_f_x = secant_method(f1, -4, -1, 0.000001, data_tab);
    for (const auto &i : data_tab) {
        std::cout << i.a << ' ' << i.b << ' ' << i.x << ' ' << i.f_a << ' '
                  << i.f_b << ' ' << i.f_x << ' ' << i.delta_eps << '\n';
    }
    std::cout << "\n\n" << x_f_x.x << ' ' << x_f_x.f_x << '\n';
    std::vector<data_table_for_fixed_point_iteration> data_tab1;
    x_and_f_x x_f_x1 = fixed_point_iteration_metod(-4, f1, 0.000001, data_tab1);
    for (const auto &i : data_tab1) {
        std::cout << i.x_k << ' ' << i.x_k_minus_one << ' ' << i.f_x_k_minus_one
                  << ' ' << i.g_x_k_minus_one << ' ' << i.delta_eps << '\n';
    }
    std::cout << "\n\n" << x_f_x1.x << ' ' << x_f_x1.f_x << '\n';
}
