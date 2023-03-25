#include "fixed_point_iteration_metod.h"
#include "secant_method.h"
#include "util.h"
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor_forward.hpp>
#include <xtensor/xview.hpp>

double f1(double x) { return 0.4 * x * x * x + 3 * x * x + x - 13; }

/*
 *
x - vector
dx - matrix
for iter
    for i in len(x)
        dx line  = (f(x + vector_xk_one*eps) - f(x))/eps
    x = x - solve(dx,fx)
    len_of_vector_ev(fx) < e exit


fx vector f(vector x)


def iteration_method(f, x0, e=1e-6, cnt_iter=100):
    x = x0.copy()
    n = len(x)
    iteration_table = []
    for i in range(cnt_iter):
        fx = f(x)
        dx = np.zeros((n, n))
        for j in range(n):
            dx[:, j] = (f(x + 1e-8*np.eye(n)[:, j]) - fx) / 1e-8
        x = x - np.linalg.solve(dx, fx)
        iteration_table.append([x[0], x[1], fx[0], fx[1], np.linalg.norm(fx)])
        if np.linalg.norm(fx) < e:
            return x, iteration_table
    print("Итерационный метод не сходится")
    return x, iteration_table
*/
typedef xt::xarray<double> arr;
typedef xt::xarray<xt::xarray<double>> matrix;
typedef arr (*func_sys)(arr);
arr func_sys1(arr x) {
    double x1 = x[0] * x[0] + x[1] * x[1] - 4;
    double x2 = x[0] - x[1] + 1;
    return {x1, x2};
}
arr solve(func_sys f, arr x0, double eps, size_t max_iter = 1000) {
    arr x = x0;
    size_t n = x0.size();
    arr fx;
    matrix dx;
    matrix identity = xt::eye(n);
    size_t iter = 0;
    do {
        dx = xt::zeros<double>({n, n});
        arr fx = f(x);
        for (size_t i = 0; i != n; ++i) {
            xt::col(dx, i) = (f(x + eps * (xt::col(identity, i))) - f(x)) / eps;
        }
        x = x - xt::linalg::solve(dx, fx);
        iter++;
    } while (xt::linalg::norm(fx, 2) > eps && iter < max_iter);
    return x;
}

int main() {
    arr x = {13, 13};
    std::cout << solve(func_sys1, x, 0.00001);
}

// int main() {
//     std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);
//
//     std::vector<data_table_for_secant_metod> data_tab;
//     x_and_f_x x_f_x = secant_method(f1, -4, -1, 0.000001, data_tab);
//     for (const auto &i : data_tab) {
//         std::cout << i.a << ' ' << i.b << ' ' << i.x << ' ' << i.f_a << ' '
//                   << i.f_b << ' ' << i.f_x << ' ' << i.delta_eps << '\n';
//     }
//     std::cout << "\n\n" << x_f_x.x << ' ' << x_f_x.f_x << '\n';
//     std::vector<data_table_for_fixed_point_iteration> data_tab1;
//     x_and_f_x x_f_x1 = fixed_point_iteration_metod(-4, f1, 0.000001,
//     data_tab1); for (const auto &i : data_tab1) {
//         std::cout << i.x_k << ' ' << i.x_k_minus_one << ' ' <<
//         i.f_x_k_minus_one
//                   << ' ' << i.g_x_k_minus_one << ' ' << i.delta_eps << '\n';
//     }
//     std::cout << "\n\n" << x_f_x1.x << ' ' << x_f_x1.f_x << '\n';
// }
