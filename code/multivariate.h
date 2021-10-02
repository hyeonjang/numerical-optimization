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

// template<typename VectorTf> gradient(VectorTf v) {

// }

template<typename VectorTf>
class Multivariate : public Method {
public:
    using function_t = multi::function_t<VectorTf>;
    Multivariate(function_t f):function(f){};

    // functions
    virtual VectorTf eval(){};

protected:
    size_t     iter=0;
    function_t function;

    // termination conditions
    bool terminate(const std::vector<VectorTf>& x) const {
        size_t k = x.size();
        return terminate1(x, k) 
            || terminate2(x, k) 
            || terminate3(x)
            || terminate4(x, k)
            || terminate5(x, k)
            || terminate6();
    };

private:
    // 1. Difference of two consecutive estimates
    inline bool terminate1(const std::vector<VectorTf>& x, size_t k) const {
        // return (x[k]-x[k-1]).cwiseAbs().eval()<VectorTf(1.f).eval();
        return true;
    };
    // 2. Relative Difference of two consecutive estimates
    inline bool terminate2(const std::vector<VectorTf>& x, size_t k) const {
        // return ((x[k]-x[k-1]).cwiseAbs()*x[k].cwiseInverse())<;
        return true;
    };
    // 3. Magnitude of Gradient
    inline bool terminate3(const std::vector<VectorTf>& x) const {
        return false;
    };
    // 4. Relative Difference of function values
    inline bool terminate4(const std::vector<VectorTf>& x, size_t k) const {
        return std::abs(function(x[k])-function(x[k-1]))/std::abs(function(x[k])) < 0.1f;
    };
    // 5. Descent direction change
    inline bool terminate5(const std::vector<VectorTf>& x, size_t k) const {
        return false;
    };
    // 6. Maximum number of iterations
    inline bool terminate6() const {
        return iter >= max_iter;
    };
};

template<typename VectorTf>
class NelderMead : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::function;
    using Base::terminate;
    using function_t = typename Base::function_t;

    NelderMead(Base base):Base(base),a(1),b(2),c(0.5){};
    NelderMead(Base base, float a, float b, float c):Base(base),a(a),b(b),c(c){};
    NelderMead(function_t func, float a, float b, float c):Base(func),a(a),b(b),c(c){};
    
    // generally works
    VectorTf eval() override {
        // 1. get the number of dimension and select threshold
        constexpr size_t dim = VectorTf::RowsAtCompileTime;

        // 2. initialize with random
        std::vector<VectorTf> simplex(dim+1);
        for(auto& s:simplex) { s = VectorTf::Random(); }

        // 3. reflection
        reflecting(simplex, a);

        // for(auto s:simplex) {
        //     std::cout << s<< std::endl;
        //     std::cout << function(s) << std::endl;
        // }

        return VectorTf(); 
    };

    // for nelder_mead method
    void reflecting(std::vector<VectorTf>& x, const float& a) {
        // 0. check termination condtion
        if(terminate(x))
            return;
        
        // 1. sorting
        std::sort(x.begin(), x.end(), 
                [&](VectorTf a, VectorTf& b){ return function(a)<function(b); });

        // 2. get mean:c
        VectorTf c = (std::accumulate(x.begin(), x.end()-1, VectorTf::Zero().eval()))/(x.size()-1);
        
        // 3. get xr
        VectorTf xr = c+a*(c-x.back());

        // 4. evaluation
        auto f1 = function(x[0]), fr = function(xr), fN = function(x.back()); 
        if(f1<=fr && f1<=fN) {
            reflecting(x, a);
        } else if(fr>=fN) {

        } else if(fr>=f1) {

        }
    };
    void expanding(std::vector<VectorTf>& x);
    void contracting(std::vector<VectorTf>& x);

private:
    float a, b, c;
};

}

#endif // __MULTIVARIATE_H__