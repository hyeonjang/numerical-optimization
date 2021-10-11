#ifndef __POWELLS__
#define __POWELLS__

#include "univariate.h"
#include "multivariate.h"

namespace numerical_optimization {

template <typename VectorTf>
class Powells : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::function;
    using function_t = typename Base::function_t;
    using Base::plot;

    VectorTf eval(float e=epsilon) override {
        constexpr size_t dim = VectorTf::RowsAtCompileTime;

        // 1. initialize
        std::vector<VectorTf> points(dim, VectorTf::Random()*10);
        std::vector<VectorTf> unit_directions(dim);
        for(size_t i=0; i<unit_directions.size(); i++) unit_directions[i][i] = 1;

        // 2. algorithm start        
        VectorTf minpoint = points[0];
        for(size_t i=0; i<100; i++) {
            for(size_t k=0; k<dim; k++) {
                uni::function_t unimodal0 = [&](float gamma){ return function(points[k] + gamma*unit_directions[k]); };
                Univariate unimethod0 = Univariate(unimodal0);
                float min_gamma = unimethod0.golden_section();
                points[k+1] = points[k] + min_gamma*unit_directions[k];
            }
            for(size_t k=0; k<dim-1; k++) {
                unit_directions[k] = unit_directions[k+1];
            }
            unit_directions[dim-1] = points[dim-1] - points[0];
            uni::function_t unimodal1 = [&](float gamma){ return function(points[0] + gamma*unit_directions[dim-1]); };
            Univariate unimethod1 = Univariate(unimodal1);
            float min_gamma1 = unimethod1.golden_section();

            auto tmp = minpoint;
            minpoint = points[0] + min_gamma1*unit_directions[dim-1];
            
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(minpoint);
#endif
            points[0] = minpoint;
            if(this->magnitude_gradient({minpoint}, e)) break;
        }
        return minpoint;
    }
};

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif //__POWELLS__