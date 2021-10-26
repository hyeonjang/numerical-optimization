#ifndef __NELDER_MEAD__
#define __NELDER_MEAD__

#include "multivariate.h"
#include "multi/termination.hpp"

namespace numerical_optimization {

template<typename VectorTf>
class NelderMead : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using function_t = typename Base::function_t;

    // constructors
    NelderMead(Base base):Base(base),alpha(1),beta(2),gamma(0.5){};
    NelderMead(function_t func):Base(func),alpha(1),beta(2),gamma(0.5){};
    NelderMead(Base base, float a, float b, float c):Base(base),alpha(a),beta(b),gamma(c){};
    NelderMead(function_t func, float a, float b, float c):Base(func),alpha(a),beta(b),gamma(c){};
    
    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h=epsilon) const {
        return Termination::eval<VectorTf, CType>(function, x, h);
    }

    // generally works
    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {
        // 1. get the number of dimension and select threshold
        constexpr size_t dim = VectorTf::RowsAtCompileTime;

        // 2. initialize with random
        std::vector<VectorTf> x(dim+1);
        for(auto& s:x) { s=VectorTf::Random(); }

        for(size_t i=0; i<this->iter; i++) {
            // 0. termination
            if(terminate<Termination::Condition::MagnitudeGradient>(x, e)) break;
            
#ifdef BUILD_WITH_PLOTTING
        for(auto t:x) plot.emplace_back(std::make_pair(t, function(t)));
#endif  
            // 1. reflection
            std::sort(
                x.begin(), x.end(), 
                [&](VectorTf l, VectorTf& r){ return function(l)<function(r); }
                );

            VectorTf c = 
                (std::accumulate(x.begin(), x.end()-1, VectorTf::Zero().eval()))/(x.size()-1);
            
            auto xr = reflecting(x.back(), c);
            auto f1 = function(x[0]), fr = function(xr), fN = function(x[x.size()-2]); 

            if(f1<=fr && fr<=fN) {
                x.back() = xr;
                continue;
          
            // 2. expansion
            } else if(fr<=f1) {
                auto xe = expanding(xr, c);
                x.back() = xe;

            // 3. contraction
            } else if(fr>=fN) {
                // last value evalution
                auto fN1 = function(x.back());

                auto xc = contracting(xr, x.back(), c, fr<fN1);
                auto fc = function(xc);

                // contraction evaluation
                if(fc<std::min(fr, fN1)) {
                    x.back() = xc;
                } else {
                    for(auto& xi : x)
                        xi = (xi + x.front())/2;
                }
            }
        }
        return x[0]; 
    };

    inline VectorTf reflecting(const VectorTf& x_last, const VectorTf& center) {
        return center + alpha*(center-x_last);
    };

    inline VectorTf expanding(const VectorTf& xr, const VectorTf& center) {
        VectorTf xe = center + beta*(xr-center);
        return (function(xe)<=function(xr)) ? xe : xr;
    };
    
    inline VectorTf contracting(const VectorTf& xr, const VectorTf& x_last, const VectorTf& center, bool check) {
        return check ? (center+gamma*(xr-center)):(center+gamma*(x_last-center));
    }

private:
    float alpha, beta, gamma;
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__NEDLER_MEAD__