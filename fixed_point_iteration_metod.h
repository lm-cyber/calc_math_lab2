#ifndef FIXED_POINT_ITERATION_METOD_H
#define FIXED_POINT_ITERATION_METOD_H

#include "util.h"
#include <vector>

typedef struct {
    double x_k_minus_one;
    double f_x_k_minus_one;
    double x_k;
    double g_x_k_minus_one;
    double delta_eps;
} data_table_for_fixed_point_iteration;
double optimal_func(double x, double (*f)(double));

double relaxation_method(double x, double (*f)(double));

double fi_func(double x, double (*f)(double));
x_and_f_x fixed_point_iteration_metod(
    double x, double (*f)(double), double eps,
    std::vector<data_table_for_fixed_point_iteration> &data_tab,
    size_t max_iter = 1000);

#endif
