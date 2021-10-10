#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

#include <algorithm>
#include <vector>
#include <numeric>
#include <functional>
#include <Eigen/Dense>
#include "fwd.h"
#include "method.h"

namespace numerical_optimization {

template<typename VectorTf>
VectorTf _gradient(const std::function<float(const VectorTf&)>& f, const VectorTf& x, float h=1);
// {
//     constexpr size_t dim = VectorTf::RowsAtCompileTime;

//     VectorTf result = VectorTf::Zero();
//     for(size_t i=0; i<dim; i++) {
//         size_t next = 
//         result[i] = (f(x[i]+h) - f(x[i]-h))/(2*h);
//     }

//     return result;
// }

template<typename VectorTf>
class Multivariate : public Method {
public:
    using function_t = multi::function_t<VectorTf>;
    Multivariate(){};
    Multivariate(function_t f):function(f){};

    // functions
    virtual VectorTf eval(){ return VectorTf(); }
// #ifdef BUILD_WITH_PLOTTING
    std::vector<VectorTf> plot;
// #endif 
    size_t     iter=0;
    function_t function;

public:
    // calculate gradient
    inline VectorTf gradient(VectorTf x, float h=1e-4) {
        return _gradient(function, x, h);
    }

    // 1. Difference of two consecutive estimates
    inline bool consecutive_difference(const std::vector<VectorTf>& x) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size(); // indexing
            flag &= (x[k1]-x[k]).norm()<epsilon;
        }
        return flag;
    };
    // 2. Relative Difference of two consecutive estimates
    inline bool consecutive_difference_relative(const std::vector<VectorTf>& x) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= (x[k1]-x[k]).norm()/x[k1].norm()<epsilon;
        }
        return flag;
    };
    // 3. Magnitude of Gradient
    inline bool magnitude_gradient(const std::vector<VectorTf>& x) const {
        bool flag = false;
        for(size_t k=0; k<x.size(); k++) {
            flag &= gradient(x[k]).norm()<epsilon;
        }
        return flag;
    };
    // 4. Relative Difference of function values
    inline bool function_value_difference_relative(const std::vector<VectorTf>& x) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= std::abs(function(x[k1])-function(x[k]))/std::abs(function(x[k1])) < 0.0001f;
        }
        return flag;
    };
    // 5. Descent direction change
    inline bool descent_direction_change(const std::vector<VectorTf>& x, const std::vector<VectorTf>& p) const {
        bool flag = false;
        for(size_t k=0; x.size(); k++) {
            flag &= (p[k]*gradient(x[k]))>=0;
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