#ifndef __POWELLS__
#define __POWELLS__

#include "univariate.h"
#include "multivariate.h"
#include "termination.h"

namespace numerical_optimization {

template <typename VectorTf>
class Powells : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::plot;
    using Base::iter;
    using Base::function;
    using function_t = typename Base::function_t;

    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h=epsilon) const {
        return Termination::eval<VectorTf, CType>(function, x, h);
    }

    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {
        constexpr size_t dim = VectorTf::RowsAtCompileTime;

        // 1. initialize
        std::vector<VectorTf> p(dim);  // points
        std::vector<VectorTf> u(dim); // unit directions
        for(size_t i=0; i<p.size(); i++) p[i] = VectorTf::Random();
        for(size_t i=0; i<u.size(); i++) u[i][i] = 1;

        // 2. algorithm start
        VectorTf xi = p[0]; // S1
        for(size_t j=0; j<this->iter; j++) {
            // 0. termination condition

            for(size_t k=0; k<dim-1; k++) {
                uni::function_t func0 = [&](float gamma){ return function(p[k] + gamma*u[k]); }; // S2
                Univariate uni0 = Univariate(func0);
                float min_gamma0 = uni0.golden_section();
                p[k+1] = p[k] + min_gamma0*u[k];
            } 
            for(size_t k=0; k<dim-1; k++) u[k] = u[k+1];    // S4
            u[dim-1] = p[dim-1] - p[0];                     // S4
            uni::function_t func1 = [&](float gamma){ return function(p[0] + gamma*u[dim-1]); }; // S5
            Univariate uni1 = Univariate(func1);
            float min_gamma1 = uni1.golden_section();
            auto tmp = xi;
            xi = p[0] + min_gamma1*u[dim-1];

#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            p[0] = xi;
            if(terminate<Termination::Condition::MagnitudeGradient>({xi}, e)) break;
        }
        return p[0];
    }
};
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif //__POWELLS__