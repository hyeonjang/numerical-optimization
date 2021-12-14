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
    template<typename vector_t, typename scalar_t = typename vector_t::Scalar>
    vector_t _gradient(const std::function<scalar_t(const vector_t&)>& f, const vector_t& x, scalar_t h=0.01);
    template<typename vector_t, typename scalar_t=typename vector_t::Scalar, typename return_t=Eigen::Matrix<typename vector_t::Scalar, vector_t::RowsAtCompileTime, vector_t::RowsAtCompileTime>>
    return_t _hessian(const std::function<scalar_t(const vector_t&)>& f, const vector_t& x, scalar_t h=0.01);
}

#include "termination.h"
#include "univariate.h"

namespace numerical_optimization {

template<typename vector_t>
class Multivariate : public Method {
public:
    using Base = Method;
    using Base::iter;
    using scalar_t = typename vector_t::Scalar;
    using matrix_t = Eigen::Matrix<typename vector_t::Scalar, vector_t::RowsAtCompileTime, vector_t::RowsAtCompileTime>;
    using function_t = std::function<scalar_t(const vector_t&)>;

    Multivariate(){};
    Multivariate(function_t f):function(f){};

    // functions
    virtual vector_t eval(const vector_t& init=vector_t::Random(), scalar_t _=epsilon){ return vector_t(); }

#ifdef BUILD_WITH_PLOTTING
    std::vector<std::pair<vector_t, scalar_t>> plot;
#endif
protected:
    function_t function;

public:

    scalar_t line_search_exact(const vector_t& xk, const vector_t& pk) {
        auto func = [&](scalar_t alpha){ return function(xk + alpha*pk); };
        return Univariate(func, 60).golden_section();
    };

    // line search for alpha
    // danger gradient nan
    scalar_t line_search_inexact(const vector_t& xk, const vector_t& pk, const scalar_t rho=0.99, const scalar_t c=0.9, const scalar_t alpha=5) const {

        // check satisfying wolfe 1st condition
        scalar_t alp = alpha;
        auto wolfe_1st = [&](const scalar_t& a)->bool { return function(xk+a*pk)<=(c*a*gradient(xk).transpose()*pk + function(xk));  };

        while(true) {
            alp = rho*alp;
            if(wolfe_1st(alp)||alp<1e-8) break;
        }
        return alp;
    }

    // calculate gradient
    inline vector_t gradient(vector_t x, scalar_t h=epsilon) const {
        return _gradient(function, x, h);
    }

    // calculate hessian & inverse
    inline decltype(auto) hessian(vector_t x, scalar_t h=epsilon) const {
        return _hessian(function, x, h);
    }

    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<vector_t>& x, scalar_t h=epsilon) const {
        return Termination::eval<CType, vector_t, scalar_t>(function, x, h);
    }
};
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif // __MULTIVARIATE_H__