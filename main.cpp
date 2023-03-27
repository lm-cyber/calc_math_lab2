#include "fixed_point_iteration_metod.h"
#include "fixed_point_iteration_metod_sys.h"
#include "secant_method.h"
#include "util.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
double f1(double x) { return 0.4 * x * x * x + 3 * x * x + x - 13; }
double f2(double x) { return std::sin(x) + 0.1; }
double f3(double x) { return std::sin(x) + std::cos(x); }

arr fsys1(arr x) {
    double x1 = x[0] * x[0] + x[1] * x[1] - 4;
    double x2 = x[0] * x[0] - x[1] * x[1] - 1;
    return {x1, x2};
}
arr fsys2(arr x) {
    double x1 = x[0] * x[0] + x[1] * x[1] - 4;
    double x2 = x[0] - x[1] + 1;
    return {x1, x2};
}
arr fsys3(arr x) {
    double x1 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 4;
    double x2 = x[0] * x[0] + x[1] * x[1] - x[2] * x[2] - 2;
    double x3 = x[0] + 2 * x[1] - x[2] - 1;
    return {x1, x2, x3};
}

int main() {
    std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);

    int choose;
    std::cout
        << "выбор уравнений и систем\n1) 0.4*x^3+3*x^2+x=13\n2) "
           "sin(x)=-0.1\n3) sin(x)=-cos(x)\n4) x^2+y^2=4, x^2-y^2=1\n5) "
           "x^2+y^2=4, x-y=-1\n6) x^2+y^2+z^2=4, x^2+y^2-z^2=2, x+2*y-z=1\n";
    std::cin >> choose;
    double eps;
    std::cout << "epsilon for metods\n";
    std::cin >> eps;
    if (choose <= 3) {
        double left, right, x;
        std::cout << "left for secant metod\n";
        std::cin >> left;
        std::cout << "right for secant metod\n";
        std::cin >> right;
        std::cout << "x for fixed point metod\n";
        std::cin >> x;
        std::vector<data_table_for_fixed_point_iteration> data_fixed_point;
        std::vector<data_table_for_secant_metod> data_secant;
        double (*f)(double);
        if (choose == 1) {
            f = f1;
        }
        if (choose == 2) {
            f = f2;
        }
        if (choose == 3) {
            f = f3;
        }
        x_and_f_x x_fx_secant = secant_method(f, left, right, eps, data_secant);
        x_and_f_x x_fx_fixed_point =
            fixed_point_iteration_metod(f, x, eps, data_fixed_point);
        std::cout << "fixed point metod\n";
        std::cout << "xk\t"
                  << "xk-1\t"
                  << "f_xk-1\t"
                  << "fi_xk-1\t"
                  << "delta\n";
        for (const auto &i : data_fixed_point) {
            std::cout << i.x_k << '\t' << i.x_k_minus_one << '\t'
                      << i.f_x_k_minus_one << '\t' << i.g_x_k_minus_one << '\t'
                      << i.delta_eps << '\n';
        }
        std::cout << "\nresult\nx = " << x_fx_fixed_point.x
                  << "\tf(x) = " << x_fx_fixed_point.f_x << "\n\n";
        std::cout << "secant method\n";
        std::cout << "left\t"
                  << "right\t"
                  << "x\t"
                  << "f(left)\t"
                  << "f(right)\t"
                  << "f(x)\t"
                  << "delta\n";
        for (const auto &i : data_secant) {
            std::cout << i.a << '\t' << i.b << '\t' << i.x << '\t' << i.f_a
                      << '\t' << i.f_b << '\t' << i.f_x << '\t' << i.delta_eps
                      << '\n';
        }
        std::cout << "\nresult\nx = " << x_fx_secant.x
                  << "\tf(x) = " << x_fx_secant.f_x << "\n\n";
        std::cout << "dif_x = " << fabs(x_fx_secant.x - x_fx_fixed_point.x)
                  << "\tdif_f(x) = "
                  << fabs(x_fx_secant.f_x - x_fx_fixed_point.f_x) << "\n";
    } else {
        std::vector<data_table_for_sys> data_sys;
        arr (*f)(arr);
        arr x;
        double x1, x2, x3;
        std::cout << "x = ";
        std::cin >> x1;
        std::cout << "y = ";
        std::cin >> x2;
        if (choose == 6) {
            std::cout << "z = ";
            std::cin >> x3;
            x = {x1, x2, x3};
            f = fsys3;
        } else {
            x = {x1, x2};
            if (choose == 4) {
                f = fsys1;
            } else {
                f = fsys2;
            }
        }
        arr result = fixed_point_iter_sys(f, x, eps, data_sys);
        for (const auto &i : data_sys) {
            std::cout << "var:\n"
                      << i.x << "\nf(vars):\n"
                      << i.fx << "\nnorm: " << i.norm << '\n';
        }
        std::cout << "\n--------------------------------------------\nresult:\n"
                  << result << '\n';
    }
}
