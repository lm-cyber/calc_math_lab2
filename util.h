#ifndef UTIL_H
#define UTIL_H

#include <cmath>
typedef struct {
    double x;
    double f_x;
} x_and_f_x;

double derivative(int dif, double point, double (*f)(double),
                  double eps = 0.00001);
double delta(double x_k, double x_k_minus_one);

#endif
