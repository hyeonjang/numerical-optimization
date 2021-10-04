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
    using Base::terminate;
    using function_t = typename Base::function_t;
    using Base::plot;

    VectorTf eval() override {
        constexpr size_t dim = VectorTf::RowsAtCompileTime;

        VectorTf result;

        std::vector<VectorTf> points(dim);
        std::vector<VectorTf> unit_directions(dim);

        points[0] = VectorTf::Random();

        // iterative methods
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

            result = points[0] + min_gamma1*unit_directions[dim-1];
// #ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(result);
// #endif
            if(terminate({result})) break;
        }
        return result;
    }

    // std::vector<VectorTf> plot() override {
    //     return nullptr;
    // }
};

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif //__POWELLS__