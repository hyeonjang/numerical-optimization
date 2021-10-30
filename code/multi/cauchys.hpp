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
    using Base::gradient;
    using function_t = typename Base::function_t;

    // constructors
    Cauchys(function_t f):Base(f){};

    // generally works
    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {
        // 1. initialize
        VectorTf xi = init;
        // 2. loop
        for(size_t i=0; i<this->iter; i++) {
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            // 1. termination
            if(terminate<Termination::Condition::MagnitudeGradient>({xi}, e)) break;

            // 2. the steepest descent direction
            VectorTf p = -1*gradient(xi)/gradient(xi).norm();

            // 3. step length
            float alpha = this->line_search_inexact(xi, p);

            // 4. update gradient
            xi = xi + alpha*p;
        }
        return xi;
    };

    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h=epsilon) const {
        return Termination::eval<VectorTf, CType>(function, x, h);
    }
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__CAUCHYS__