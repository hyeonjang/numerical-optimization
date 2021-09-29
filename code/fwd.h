#pragma once
#ifndef __FWD_H__
#define __FWD_H__

#include <functional>
#include <limits>
#include <vector>

namespace numerical_optimization {

namespace uni {
    using function_t = std::function<float(const float&)>;
    using boundary_t = std::pair<float, float>;
};

namespace multi {
    using function_t = std::function<float(const std::vector<float>&)>;
};

constexpr float  MIN = 1e-4;// std::numeric_limits<float>::min();
constexpr float  MAX = std::numeric_limits<float>::max();
constexpr float  GOLDEN_RATIO = 1.f/1.618033988749895f;
constexpr size_t FIBONACCI_MAX = 46;

class Univariate;
class Multivariate;

}

#endif // __FWD_H__