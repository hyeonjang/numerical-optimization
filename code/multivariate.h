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
class Multivariate : public Method {
public:
    using function_t = multi::function_t<VectorTf>;
    Multivariate(function_t f):function(f){};

    // functions
    virtual VectorTf eval(){ return VectorTf(); };

protected:
    size_t     iter=0;
    function_t function;

    // termination conditions
    bool terminate(const std::vector<VectorTf>& x) const {
        size_t k = x.size()-1;
        return terminate1(x, k) 
            || terminate2(x, k) 
            || terminate3(x)
            || terminate4(x)
            || terminate5(x, k)
            || terminate6();
    };

private:
    // 1. Difference of two consecutive estimates
    inline bool terminate1(const std::vector<VectorTf>& x, size_t k) const {
        bool flag = false;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            // flag &= Eigen::abs((x[k1]-x[k-1])).eval()<VectorTf(0.1f);
        }
        return flag;
    };
    // 2. Relative Difference of two consecutive estimates
    inline bool terminate2(const std::vector<VectorTf>& x, size_t k) const {
        // return ((x[k]-x[k-1]).cwiseAbs()*x[k].cwiseInverse())<;
        return false;
    };
    // 3. Magnitude of Gradient
    inline bool terminate3(const std::vector<VectorTf>& x) const {
        return false;
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

    NelderMead(Base base):Base(base),alpha(1),beta(2),gamma(0.5){};
    NelderMead(Base base, float a, float b, float c):Base(base),alpha(a),beta(b),gamma(c){};
    NelderMead(function_t func, float a, float b, float c):Base(func),alpha(a),beta(b),gamma(c){};
    
    // generally works
    VectorTf eval() override {
        // 1. get the number of dimension and select threshold
        constexpr size_t dim = VectorTf::RowsAtCompileTime;

        // 2. initialize with random
        std::vector<VectorTf> simplex(dim+1);
        for(auto& s:simplex) { s = VectorTf::Random(); }

        // 3. algorithm start: reflection
        reflecting(simplex);

        // 4. evaluation       
        // for(size_t i=0; i<100; i++) {
        //     auto f1 = function(x[0]), fr = function(xr), fN = function(x.back()); 
        
        //     if(f1<=fr && f1<=fN) {
        //         x.back() = xr;
        //         xr = reflecting(x);
        //     } else if(fr>=fN) {
        //         xr = expanding(xr);
        //         x.back() = xr;
        //     } else if(fr>=f1) {
        //         xr = contracting(xr);
        //         x.back() = xr;
        //     }
        // }

        for(auto s:simplex) {
            std::cout << s<< std::endl;
            std::cout << function(s) << std::endl;
        }

        return simplex[0]; 
    };

    // for nelder_mead method
    void reflecting(std::vector<VectorTf>& x) {
        // 0. check termination condtion
        if(terminate(x)) return;
        
        // 1. sorting
        std::sort(x.begin(), x.end(), [&](VectorTf l, VectorTf& r){ return function(l)<function(r); });

        // 2. get mean:c
        VectorTf c = (std::accumulate(x.begin(), x.end()-1, VectorTf::Zero().eval()))/(x.size()-1);
        
        // 3. get xr
        VectorTf xr = c + alpha*(c-x.back());

        // 4. evaluation
        auto f1 = function(x[0]), fr = function(xr), fN = function(x[x.size()-2]); 
        if(f1<=fr && f1<=fN) {
            reflecting(x);
        } else if(fr>=fN) {
            expanding(x, xr, c);
        } else if(fr>=f1) {
            contracting(x, xr, c, fr);
        }
    };
    void expanding(std::vector<VectorTf>& x, const VectorTf& xr, const VectorTf& c) {
        VectorTf xe = c + beta*(xr-c);
        x.back() = (function(xe)<=function(xr)) ? xe : xr;

        reflecting(x);
    };
    void contracting(std::vector<VectorTf>& x, const VectorTf& xr, const VectorTf& c, float fr) {
        
        auto fN1 = function(x.back());
        VectorTf xc = (fr<fN1) ? (c+gamma*(xr-c)):(c+gamma*(x.back()-c));

        auto fc = function(xc);
        if(fc<std::min(fr, fN1)) {
            x.back() = xc;
        } else {
            for(auto& xi : x) 
                xi = (xi + x.front())/2;
        }

        reflecting(x);
    };

private:
    float alpha, beta, gamma;
};

}

#endif // __MULTIVARIATE_H__