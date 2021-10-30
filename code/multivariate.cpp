#include <cmath>
#include <iostream>

#include "multivariate.h"

using namespace Eigen;

namespace numerical_optimization {

// two-variables case
// specialization of the vector2f case
template<>
Vector2f _gradient<Vector2f>(const std::function<float(const Vector2f&)>& f, const Vector2f& x, float h) {
    using v2 = Vector2f;

    float dx = 3*f(v2(x[0]-4*h, x[1]))-32*f(v2(x[0]-3*h, x[1]))+168*f(v2(x[0]-2*h, x[1]))-672*f(v2(x[0]-h, x[1]))
              -3*f(v2(x[0]+4*h, x[1]))+32*f(v2(x[0]+3*h, x[1]))-168*f(v2(x[0]+2*h, x[1]))+672*f(v2(x[0]+h, x[1]));
    float dy = 3*f(v2(x[0], x[1]-4*h))-32*f(v2(x[0], x[1]-3*h))+168*f(v2(x[0], x[1]-2*h))-672*f(v2(x[0], x[1]-h))
              -3*f(v2(x[0], x[1]+4*h))+32*f(v2(x[0], x[1]+3*h))-168*f(v2(x[0], x[1]+2*h))+672*f(v2(x[0], x[1]+h));

    float inv = (1/(h*840));
    return v2(dx, dy)*inv;
}

// specialization of the vector2f case
template<>
Matrix2f _hessian<Vector2f>(const std::function<float(const Vector2f&)>& f, const Vector2f& x, float h) {
    using vec2 = Vector2f;

    h=0.01;
    auto dfdx = [&](vec2 x){ 
        float inv = (1/(h*840));
        float app = 3*f(vec2(x[0]-4*h, x[1]))-32*f(vec2(x[0]-3*h, x[1]))+168*f(vec2(x[0]-2*h, x[1]))-672*f(vec2(x[0]-h, x[1]))
                        -3*f(vec2(x[0]+4*h, x[1]))+32*f(vec2(x[0]+3*h, x[1]))-168*f(vec2(x[0]+2*h, x[1]))+672*f(vec2(x[0]+h, x[1]));
            return app*inv;
        };

    auto dfdy = [&](vec2 x){ 
        float inv = (1/(h*840));
        float app = 3*f(vec2(x[0], x[1]-4*h))-32*f(vec2(x[0], x[1]-3*h))+168*f(vec2(x[0], x[1]-2*h))-672*f(vec2(x[0], x[1]-h))
                    -3*f(vec2(x[0], x[1]+4*h))+32*f(vec2(x[0], x[1]+3*h))-168*f(vec2(x[0], x[1]+2*h))+672*f(vec2(x[0], x[1]+h));
        return app*inv;
        };

    float dxx = f(vec2(x[0]+2*h, x[1]))-2*f(vec2(x[0], x[1]))+f(vec2(x[0]-2*h, x[1]));
    float dxy = f(vec2(x[0]+h, x[1]+h))-f(vec2(x[0]-h, x[1]+h))-f(vec2(x[0]+h, x[1]-h)) + f(vec2(x[0]-h, x[1]-h));
    float dyx = f(vec2(x[0]+h, x[1]+h))-f(vec2(x[0]+h, x[1]-h))-f(vec2(x[0]-h, x[1]+h)) + f(vec2(x[0]-h, x[1]-h));
    float dyy = f(vec2(x[0], x[1]+2*h))-2*f(vec2(x[0], x[1]))+f(vec2(x[0], x[1]-2*h));

    Matrix2f m;
    m << dxx, dxy, dyx, dyy;
    float inv = 1/(4*h*h);
    return m*= inv;
}
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////