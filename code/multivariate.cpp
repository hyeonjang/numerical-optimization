#include <cmath>
#include <iostream>

#include "multivariate.h"

using namespace Eigen;

namespace numerical_optimization {

// two-variables case
// specialization of the vector2f case
template<>
Vector2f _gradient<Vector2f>(const std::function<Vector2f::Scalar(const Vector2f&)>& f, const Vector2f& x, Vector2f::Scalar h) {
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
Matrix2f _hessian<Vector2f>(const std::function<Vector2f::Scalar(const Vector2f&)>& f, const Vector2f& x, Vector2f::Scalar h) {
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

template<>
Vector2d _gradient<Vector2d>(const std::function<Vector2d::Scalar(const Vector2d&)>& f, const Vector2d& x, Vector2d::Scalar h) {
    using v2 = Vector2d;

    double dx = 3*f(v2(x[0]-4*h, x[1]))-32*f(v2(x[0]-3*h, x[1]))+168*f(v2(x[0]-2*h, x[1]))-672*f(v2(x[0]-h, x[1]))
              -3*f(v2(x[0]+4*h, x[1]))+32*f(v2(x[0]+3*h, x[1]))-168*f(v2(x[0]+2*h, x[1]))+672*f(v2(x[0]+h, x[1]));
    double dy = 3*f(v2(x[0], x[1]-4*h))-32*f(v2(x[0], x[1]-3*h))+168*f(v2(x[0], x[1]-2*h))-672*f(v2(x[0], x[1]-h))
              -3*f(v2(x[0], x[1]+4*h))+32*f(v2(x[0], x[1]+3*h))-168*f(v2(x[0], x[1]+2*h))+672*f(v2(x[0], x[1]+h));

    double inv = (1/(h*840));
    return v2(dx, dy)*inv;
}

// specialization of the vector2d case
template<>
Matrix2d _hessian<Vector2d>(const std::function<Vector2d::Scalar(const Vector2d&)>& f, const Vector2d& x, Vector2d::Scalar h) {
    using vec2 = Vector2d;

    h=0.01;
    auto dfdx = [&](vec2 x){ 
        double inv = (1/(h*840));
        double app = 3*f(vec2(x[0]-4*h, x[1]))-32*f(vec2(x[0]-3*h, x[1]))+168*f(vec2(x[0]-2*h, x[1]))-672*f(vec2(x[0]-h, x[1]))
                        -3*f(vec2(x[0]+4*h, x[1]))+32*f(vec2(x[0]+3*h, x[1]))-168*f(vec2(x[0]+2*h, x[1]))+672*f(vec2(x[0]+h, x[1]));
            return app*inv;
        };

    auto dfdy = [&](vec2 x){ 
        double inv = (1/(h*840));
        double app = 3*f(vec2(x[0], x[1]-4*h))-32*f(vec2(x[0], x[1]-3*h))+168*f(vec2(x[0], x[1]-2*h))-672*f(vec2(x[0], x[1]-h))
                    -3*f(vec2(x[0], x[1]+4*h))+32*f(vec2(x[0], x[1]+3*h))-168*f(vec2(x[0], x[1]+2*h))+672*f(vec2(x[0], x[1]+h));
        return app*inv;
        };

    double dxx = f(vec2(x[0]+2*h, x[1]))-2*f(vec2(x[0], x[1]))+f(vec2(x[0]-2*h, x[1]));
    double dxy = f(vec2(x[0]+h, x[1]+h))-f(vec2(x[0]-h, x[1]+h))-f(vec2(x[0]+h, x[1]-h)) + f(vec2(x[0]-h, x[1]-h));
    double dyx = f(vec2(x[0]+h, x[1]+h))-f(vec2(x[0]+h, x[1]-h))-f(vec2(x[0]-h, x[1]+h)) + f(vec2(x[0]-h, x[1]-h));
    double dyy = f(vec2(x[0], x[1]+2*h))-2*f(vec2(x[0], x[1]))+f(vec2(x[0], x[1]-2*h));

    Matrix2d m;
    m << dxx, dxy, dyx, dyy;
    double inv = 1/(4*h*h);
    return m*= inv;
}

// specialization of the vector2d case
// template<>
// Matrix2d _jacobian<Vector2d>(const std::function<float(const Vector2d&)>& f, const Vector2d& x, float h) {

// }
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////