#ifndef FIXED_POINT_ITERATION_METOD_SYS_H
#define FIXED_POINT_ITERATION_METOD_SYS_H

#include "util.h"
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xtensor_forward.hpp>
#include <xtensor/xview.hpp>

typedef xt::xarray<double> arr;
typedef struct {
    arr x;
    arr fx;
    double norm;
} data_table_for_sys;
typedef arr (*func_sys)(arr);
double norm(arr x);

arr fixed_point_iter_sys(func_sys f, arr x0, double eps,
                         std::vector<data_table_for_sys> &data_table,
                         size_t max_iter = 1000);

#endif
