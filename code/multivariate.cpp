#include <iostream>
#include "multivariate.h"

using namespace Eigen;

namespace numerical_optimization {

// two-variables case
// https://stackoverflow.com/questions/38854363/is-there-any-standard-way-to-calculate-the-numerical-gradient
template <>
Vector2f _gradient<Vector2f>(const std::function<float(const Vector2f&)>& f, const Vector2f& x, float h) {

    Vector2f result = Vector2f::Zero();

    // @@todo more optimize
    Vector2f eps = Vector2f(h, h);
    // relative 'h' value
    // cannot work for vector has zero: it results NaN
    if(x[0]!=0 && x[1]!=0) {
        eps[0] = x[0]*sqrtf(eps[0]);
        eps[1] = x[1]*sqrtf(eps[1]);
    }

    result[0] = (f(Vector2f(x[0]+eps[0], x[1]))-f(Vector2f(x[0]-eps[0], x[1]))) / (2*eps[0]);
    result[1] = (f(Vector2f(x[0], x[1]+eps[1]))-f(Vector2f(x[0], x[1]-eps[1]))) / (2*eps[1]);

    // std::cout << "result: " << result[0] << ", " << result[1] << std::endl;

    return result;
}

} // namespace numerical_optimization
