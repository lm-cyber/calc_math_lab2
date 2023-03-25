#include "util.h"

double derivative(int dif, double point, double (*f)(double),
                  double eps) { // нам важен знак, большой esp приводит к нулю
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
