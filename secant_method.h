#ifndef SECANT_METHOD_H
#define SECANT_METHOD_H
#include "util.h"
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

x_and_f_x secant_method(double (*f)(double), double a, double b, double eps,
                        std::vector<data_table_for_secant_metod> &data_tab,
                        size_t max_iter = 10000);

#endif
