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
protected:
    size_t     iter=0;
    function_t function;

public:
    // calculate gradient
    inline VectorTf gradient(VectorTf x, float h=1e-4) {
        return _gradient(function, x, h);
    }

    // 1. Difference of two consecutive estimates
    inline bool terminate1(const std::vector<VectorTf>& x, size_t k) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size(); // indexing
            flag &= (x[k1]-x[k]).norm()<0.1f;
        }
        return flag;
    };
    // 2. Relative Difference of two consecutive estimates
    inline bool terminate2(const std::vector<VectorTf>& x, size_t k) const {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= (x[k1]-x[k]).norm()/x[k1].norm()<0.1f;
        }
        return flag;
    };
    // 3. Magnitude of Gradient
    inline bool terminate3(const std::vector<VectorTf>& x) const {
        bool flag = false;
        // not implemented
        return flag;
    };
    // 4. Relative Difference of function values
    inline bool terminate4(const std::vector<VectorTf>& x) const {
        bool flag = false;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= std::abs(function(x[k1])-function(x[k]))/std::abs(function(x[k1])) < 0.1f;
        }
        return flag;
    };
    // 5. Descent direction change
    inline bool terminate5(const std::vector<VectorTf>& x, size_t k) const {
        bool flag = false;
        // not implemented
        return flag;
    };
    // 6. Maximum number of iterations
    inline bool terminate6() const {
        return iter >= max_iter;
    };
};
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif // __MULTIVARIATE_H__