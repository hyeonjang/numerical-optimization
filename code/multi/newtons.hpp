#ifndef __NEWTONS__
#define __NEWTONS__

#include "multivariate.h"
#include "multi/termination.hpp"

namespace numerical_optimization {

template<typename VectorTf>
class Newtons : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::plot;

    using Base::function;
    using function_t = typename Base::function_t;

    // constructors
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
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            // 2. gradient update
            xi = xi - this->hessian(xi).inverse()*this->gradient(xi);
        }
        return xi;
    };
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__CAUCHYS__