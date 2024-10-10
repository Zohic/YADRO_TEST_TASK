#pragma once

#define pi 3.14159265358979323846

#define prec 2 // 1 for float, 2 for double

#if prec == 1
#define number_t float
#define sin_f sinf
#define cos_f cosf
#elif prec == 2
#define number_t double
#define sin_f sin
#define cos_f cos
#else
static_assert(false, "wrong float precision defined in config.hpp, use 1 or 2");
#endif// prec

#include <complex>

using complex_t = std::complex<number_t>;

const complex_t i_unit(0, 1);
const number_t d_pi = (number_t)pi * 2.0;
const complex_t i2pi = i_unit * d_pi;