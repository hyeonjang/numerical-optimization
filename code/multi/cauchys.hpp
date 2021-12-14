#ifndef __CAUCHYS__
#define __CAUCHYS__

#include "multivariate.h"
#include "termination.h"

namespace numerical_optimization {

template<typename vector_t>
class Cauchys : public Multivariate<vector_t> {
public:
    using Base = Multivariate<vector_t>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using Base::gradient;

    using scalar_t = typename Base::scalar_t;
    using function_t = typename Base::function_t;

    // constructors
    Cauchys(function_t f):Base(f){};

    // generally works
    vector_t eval(const vector_t& init=vector_t::Random(), scalar_t e=epsilon) override {
        // 1. initialize
        vector_t xi = init;
        // 2. loop
        for(size_t i=0; i<this->iter; i++) {
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            // 1. termination
            if(terminate<Termination::Condition::MagnitudeGradient>({xi}, e)) break;

            // 2. the steepest descent direction
            vector_t p = -1*gradient(xi)/gradient(xi).norm();

            // 3. step length
            float alpha = this->line_search_inexact(xi, p);

            // 4. update gradient
            xi = xi + alpha*p;
        }
        return xi;
    };

    // termination
    template<Termination::Condition CType> 
    bool terminate(const std::vector<vector_t>& x, scalar_t h, scalar_t eps=epsilon) {
        return Termination::eval<CType, vector_t, scalar_t>(function, x, h, eps);
    }
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__CAUCHYS__