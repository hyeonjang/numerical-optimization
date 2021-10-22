#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

#include <algorithm>
#include <vector>
#include <numeric>
#include <functional>
#include <Eigen/Dense>
#include "fwd.h"
#include "method.h"

#define Log(x) printf("position: [%f, %f], value: %f\n", x[0], x[1], function(x))

namespace numerical_optimization {

template<typename VectorTf>
VectorTf _gradient(const std::function<float(const VectorTf&)>& f, const VectorTf& x, float h=1);
template<typename VectorTf, typename ReturnType = Eigen::Matrix<typename VectorTf::Scalar, VectorTf::RowsAtCompileTime, VectorTf::RowsAtCompileTime>>
ReturnType _hessian(const std::function<float(const VectorTf&)>& f, const VectorTf& x, float h=1);

template<typename VectorTf>
class Multivariate : public Method {
public:
    using function_t = multi::function_t<VectorTf>;
    Multivariate(){};
    Multivariate(function_t f):function(f){};

    // functions
    // virtual VectorTf eval(float _=epsilon){ return VectorTf(); }
    virtual VectorTf eval(const VectorTf& init=VectorTf::Random(), float _=epsilon){ return VectorTf(); }

#ifdef BUILD_WITH_PLOTTING
    std::vector<std::pair<VectorTf, float>> plot;
#endif
protected:
    size_t     iter=0;
    function_t function;

public:
    // line search for alpha
    inline float line_search_inexact(const VectorTf& xk, const VectorTf& pk, float p, float c, float alpha=0.1) const {

        // check satisfying wolfe 1st condition
        auto wolfe_1st = [&](float a){ return function(xk+a*pk)<=(a*c*gradient(xk).transpose()*pk + function(xk)); };

        for(size_t i=0; i<100 && !wolfe_1st(alpha); ++i) {
            alpha = alpha*p;
            // if(wolfe_1st(alpha)) break; ??
        }
        return alpha;
    }

    inline float line_search_exact();

    // calculate gradient
    inline VectorTf gradient(VectorTf x, float h=epsilon) const {
        return _gradient(function, x, h);
    }

    // calculate hessian & inverse
    inline decltype(auto) hessian(VectorTf x, float h=epsilon) const {
        return _hessian(function, x, h);
    }

    // 1. Difference of two consecutive estimates
    inline bool consecutive_difference(const std::vector<VectorTf>& x, float eps=epsilon) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size(); // indexing
            flag &= (x[k1]-x[k]).norm()<eps;
        }
        return flag;
    };
    // 2. Relative Difference of two consecutive estimates
    inline bool consecutive_difference_relative(const std::vector<VectorTf>& x, float eps=epsilon) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= (x[k1]-x[k]).norm()/x[k1].norm()<eps;
        }
        return flag;
    };
    // 3. Magnitude of Gradient
    inline bool magnitude_gradient(const std::vector<VectorTf>& x, float eps=epsilon) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            flag &= gradient(x[k]).norm()<eps;
        }
        return flag;
    };
    // 4. Relative Difference of function values
    inline bool function_value_difference_relative(const std::vector<VectorTf>& x, float eps=epsilon) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= std::abs(function(x[k1])-function(x[k]))/std::abs(function(x[k1])) < eps;
        }
        return flag;
    };
    // 5. Descent direction change
    inline bool descent_direction_change(const std::vector<VectorTf>& x, const std::vector<VectorTf>& p) const {
        bool flag = true;
        for(size_t k=0; x.size(); k++) {
            flag &= (p[k]*gradient(x[k]))>=0.f;
        }
        return flag;
    };
    // 6. Maximum number of iterations
    inline bool over_maximum_iteration() const {
        return iter >= max_iter;
    };

};
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif // __MULTIVARIATE_H__