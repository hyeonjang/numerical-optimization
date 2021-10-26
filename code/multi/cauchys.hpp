#ifndef __CAUCHYS__
#define __CAUCHYS__

#include "multivariate.h"
#include "multi/termination.hpp"

namespace numerical_optimization {

template<typename VectorTf>
class Cauchys : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using function_t = typename Base::function_t;

    // constructors
    Cauchys(function_t f):Base(f), alpha(0.05){};
    Cauchys(function_t f, float a):Base(f), alpha(a){};
    
    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h=epsilon) const {
        return Termination::eval<VectorTf, CType>(function, x, h);
    }

    // generally works
    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {
        
        VectorTf xi = init;

        for(size_t i=0; i<this->iter; i++) {

            // 1. termination
            if(terminate<Termination::Condition::MagnitudeGradient>({xi}, e)) break;

            // 2. the steepest descent direction
            VectorTf p = -1  * this->gradient(xi)/this->gradient(xi).norm();

            // 3. step length
            alpha = this->line_search_inexact(xi, p, 0.05, 0.05);

            // 4. update gradient
            xi = xi + alpha*p;

#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
        }
        return xi;
    };

private:
    float alpha;
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__CAUCHYS__