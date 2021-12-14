#ifndef __NELDER_MEAD__
#define __NELDER_MEAD__

#include "multivariate.h"
#include "termination.h"

namespace numerical_optimization {

template<typename vector_t>
class NelderMead : public Multivariate<vector_t> {
public:
    using Base = Multivariate<vector_t>;
    using Base::Base;
    using Base::plot;
    using Base::function;

    using scalar_t = typename Base::scalar_t;
    using function_t = typename Base::function_t;

    // constructors
    NelderMead(Base base):Base(base),alpha(1),beta(2),gamma(0.5){};
    NelderMead(function_t func):Base(func),alpha(1),beta(2),gamma(0.5){};
    NelderMead(Base base, float a, float b, float c):Base(base),alpha(a),beta(b),gamma(c){};
    NelderMead(function_t func, float a, float b, float c):Base(func),alpha(a),beta(b),gamma(c){};
    
    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<vector_t>& x, scalar_t h, scalar_t eps=epsilon) {
        return Termination::eval<CType, vector_t, scalar_t>(function, x, h, eps);
    }

    // generally works
    vector_t eval(const vector_t& init=vector_t::Random(), scalar_t e=epsilon) override {
        // 1. get the number of dimension and select threshold
        constexpr size_t dim = vector_t::RowsAtCompileTime;

        // 2. initialize with random
        std::vector<vector_t> x(dim+1);
        for(auto& s:x) { s=vector_t::Random(); }

        for(size_t i=0; i<this->iter; i++) {
            // 0. termination
            if(terminate<Termination::Condition::MagnitudeGradient>(x, e)) break;
            
#ifdef BUILD_WITH_PLOTTING
        for(auto t:x) plot.emplace_back(std::make_pair(t, function(t)));
#endif  
            // 1. reflection
            std::sort(
                x.begin(), x.end(), 
                [&](vector_t l, vector_t& r){ return function(l)<function(r); }
                );

            vector_t c = 
                (std::accumulate(x.begin(), x.end()-1, vector_t::Zero().eval()))/(x.size()-1);
            
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

    inline vector_t reflecting(const vector_t& x_last, const vector_t& center) {
        return center + alpha*(center-x_last);
    };

    inline vector_t expanding(const vector_t& xr, const vector_t& center) {
        vector_t xe = center + beta*(xr-center);
        return (function(xe)<=function(xr)) ? xe : xr;
    };
    
    inline vector_t contracting(const vector_t& xr, const vector_t& x_last, const vector_t& center, bool check) {
        return check ? (center+gamma*(xr-center)):(center+gamma*(x_last-center));
    }

private:
    float alpha, beta, gamma;
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__NEDLER_MEAD__