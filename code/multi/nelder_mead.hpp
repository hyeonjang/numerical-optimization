#ifndef __NELDER_MEAD__
#define __NELDER_MEAD__

#include "multivariate.h"

namespace numerical_optimization {

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
        for(auto& s:simplex) { s=VectorTf::Random(); }

        // 3. algorithm start: reflection
        reflecting(simplex);

        return simplex[0]; 
    };

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
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__NEDLER_MEAD__