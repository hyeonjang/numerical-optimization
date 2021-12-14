#ifndef __NEWTONS__
#define __NEWTONS__

#include "multivariate.h"
#include "termination.h"

namespace numerical_optimization {

template<typename vector_t>
class Newtons : public Multivariate<vector_t> {
public:
    using Base = Multivariate<vector_t>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using Base::gradient;
    using Base::hessian;
    using scalar_t = typename vector_t::Scalar;
    using function_t = typename Base::function_t;

    // constructors
    template<Termination::Condition CType> 
    bool terminate(const std::vector<vector_t>& x, scalar_t h, scalar_t eps=epsilon) {
        return Termination::eval<CType, vector_t, scalar_t>(function, x, h, eps);
    }

    // generally works
    vector_t eval(const vector_t& init=vector_t::Random(), float e=epsilon) override {
        vector_t xi = init;
        for(size_t i=0; i<this->iter; i++) {
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            // 1. termination
            if(terminate<Termination::Condition::MagnitudeGradient>({xi}, e)) break;

            // 2. gradient update
            xi = xi - hessian(xi).inverse()*gradient(xi);
        }
        return xi;
    };
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__CAUCHYS__