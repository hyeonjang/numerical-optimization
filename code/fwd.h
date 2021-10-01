#pragma once
#ifndef __FWD_H__
#define __FWD_H__

#include <functional>
#include <limits>
#include <vector>

#include <Eigen/Dense>

namespace numerical_optimization {

template<typename T> using Tfunction = std::function<float(const T&)>;

namespace uni {
    using function_t = Tfunction<float>;
    using boundary_t = std::pair<float, float>;
};

namespace multi {
    template<typename VectorT> 
    using function_t = Tfunction<VectorT>;
};

constexpr float  MIN = 1e-4;// std::numeric_limits<float>::min();
constexpr float  MAX = std::numeric_limits<float>::max();
constexpr float  GOLDEN_RATIO = 1.f/1.618033988749895f;
constexpr size_t FIBONACCI_MAX = 46;

class Method;
class Univariate;
template<typename VectorTf> class Multivariate;

}

#endif // __FWD_H__