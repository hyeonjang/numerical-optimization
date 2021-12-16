#pragma once
#ifndef __FWD_H__
#define __FWD_H__

#include <functional>
#include <limits>
#include <vector>
#include <random>

#include <Eigen/Dense>

namespace numerical_optimization {

template<typename T> using Tfunction = std::function<float(const T&)>;

constexpr float  MIN = 1e-4;// std::numeric_limits<float>::min();
constexpr float  MAX = std::numeric_limits<float>::max();
constexpr float  GOLDEN_RATIO = 1.f/1.618033988749895f;
constexpr size_t FIBONACCI_MAX = 46;
constexpr size_t max_iter = 10000;
constexpr float  epsilon = 1e-6;

static std::random_device rd;
static std::mt19937 gen(rd());

namespace uni {
    using function_t = Tfunction<float>;
    using boundary_t = std::pair<float, float>;
};

class Method;

// univariate function methods ~ homework#2
class Univariate;

// multivariate function methods
template<typename VectorTf> class Multivariate;
template<typename VectorTf> class NelderMead;
template<typename VectorTf> class Powells;
}

#endif // __FWD_H__