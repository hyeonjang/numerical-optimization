#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <functional>
#include <type_traits>
#include <Eigen/Dense>
#include "fwd.h"
#include "method.h"

#define Log(x) printf("position: [%f, %f], value: %f\n", x[0], x[1], function(x))

using namespace Eigen;
namespace numerical_optimization {
    template<typename VectorT>
    VectorT _gradient(const std::function<float(const VectorT&)>& f, const VectorT& x, float h=0.01);
    template<typename VectorT, typename ReturnType = Eigen::Matrix<typename VectorT::Scalar, VectorT::RowsAtCompileTime, VectorT::RowsAtCompileTime>>
    ReturnType _hessian(const std::function<float(const VectorT&)>& f, const VectorT& x, float h=0.01);
}

#include "termination.h"
#include "univariate.h"

namespace numerical_optimization {

template<typename VectorTf>
class Multivariate : public Method {
public:
    using Base = Method;
    using Base::iter;
    using function_t = multi::function_t<VectorTf>;

    Multivariate(){};
    Multivariate(function_t f):function(f){};

    // functions
    virtual VectorTf eval(const VectorTf& init=VectorTf::Random(), float _=epsilon){ return VectorTf(); }

#ifdef BUILD_WITH_PLOTTING
    std::vector<std::pair<VectorTf, float>> plot;
#endif
protected:
    function_t function;

public:

    double line_search_exact(const VectorTf& xk, const VectorTf& pk) {
        auto func = [&](double alpha){ return function(xk + alpha*pk); };
        return Univariate(func, 60).golden_section();
    };

    // line search for alpha
    // danger gradient nan
    double line_search_inexact(const VectorTf& xk, const VectorTf& pk, const double rho=0.99, const double c=0.9, const double alpha=5) const {

        // check satisfying wolfe 1st condition
        double alp = alpha;
        auto wolfe_1st = [&](const double& a)->bool { return function(xk+a*pk)<=(c*a*gradient(xk).transpose()*pk + function(xk));  };

        while(true) {
            alp = rho*alp;
            if(wolfe_1st(alp)||alp<1e-8) break;
        }
        return alp;
    }

    // calculate gradient
    inline VectorTf gradient(VectorTf x, float h=epsilon) const {
        return _gradient(function, x, h);
    }

    // calculate hessian & inverse
    inline decltype(auto) hessian(VectorTf x, float h=epsilon) const {
        return _hessian(function, x, h);
    }

    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h=epsilon) const {
        return Termination::eval<VectorTf, CType>(function, x, h);
    }
};
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif // __MULTIVARIATE_H__