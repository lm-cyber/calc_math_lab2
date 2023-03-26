
#include "fixed_point_iteration_metod_sys.h"
double norm(arr x) {
    double sum = 0;
    for (const double &i : x) {
        sum += i * i;
    }
    return sqrt(sum);
}

arr fixed_point_iter_sys(func_sys f, arr x0, double eps,
                         std::vector<data_table_for_sys> &data_table,
                         size_t max_iter) {
    arr x = x0;
    size_t n = x0.size();
    arr fx;
    size_t iter = 0;
    do {
        arr iden = xt::eye(n);
        arr dx = xt::zeros<double>({n, n});
        fx = f(x);
        for (size_t i = 0; i != n; ++i) {
            xt::col(dx, i) = ((f(x + eps * xt::col(iden, i))) - f(x)) / eps;
        }
        x = x - xt::linalg::solve(dx, fx);
        data_table.push_back({x, fx, norm(fx)});
        iter++;
    } while (norm(fx) > eps && iter < max_iter);

    return x;
}
