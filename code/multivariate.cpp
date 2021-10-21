#include <cmath>
#include <iostream>

#include "multivariate.h"

using namespace Eigen;

namespace numerical_optimization {

// two-variables case
// specialization of the vector2f case
// https://stackoverflow.com/questions/38854363/is-there-any-standard-way-to-calculate-the-numerical-gradient
template<>
Vector2f _gradient<Vector2f>(const std::function<float(const Vector2f&)>& f, const Vector2f& x, float h) {

    Vector2f result = Vector2f::Zero();

    Vector2f eps = Vector2f(h, h);
    // relative 'h' value
    // cannot work for vector has zero: it results NaN
    if(x[0]!=0 && x[1]!=0) {
        eps[0] = x[0]*sqrtf(eps[0]);
        eps[1] = x[1]*sqrtf(eps[1]);
    }

    result[0] = (f(Vector2f(x[0]+eps[0], x[1]))-f(Vector2f(x[0]-eps[0], x[1]))) / (2*eps[0]);
    result[1] = (f(Vector2f(x[0], x[1]+eps[1]))-f(Vector2f(x[0], x[1]-eps[1]))) / (2*eps[1]);

    return result;
}

// specialization of the vector2f case
template<>
Matrix2f _hessian<Vector2f>(const std::function<float(const Vector2f&)>& f, const Vector2f& x, float h) {
    
    float a = 0.01;

    float dxx = (f(Vector2f(x[0]+2*a, x[1])) - 2*f(Vector2f(x[0], x[1])) 
                + f(Vector2f(x[0]-2*a, x[1])));
    float dxy = (f(Vector2f(x[0]+a, x[1]+a)) - f(Vector2f(x[0]-a, x[1]+a)) 
                - f(Vector2f(x[0]+a, x[1]-a)) + f(Vector2f(x[0]-a, x[1]-a)));
    float dyx = (f(Vector2f(x[0]+a, x[1]+a)) - f(Vector2f(x[0]+a, x[1]-a)) 
                - f(Vector2f(x[0]-a, x[1]+a)) + f(Vector2f(x[0]-a, x[1]-a)));
    float dyy = (f(Vector2f(x[0], x[1]+2*a)) - 2*f(Vector2f(x[0], x[1])) 
                + f(Vector2f(x[0], x[1]-2*a)));
    
    Matrix2f m;
    m << dxx, dxy, dyx, dyy;
    m /= 4*a*a;

    return m;
}

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////