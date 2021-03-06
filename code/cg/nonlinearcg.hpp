#ifndef __NonlinearCG__
#define __NonlinearCG__

#include "multivariate.h"

namespace numerical_optimization {
namespace nonlinear_cg {
enum Beta { CG_FR, CG_PR, CG_HS };
};

template<typename vector_t, nonlinear_cg::Beta BetaMethod>
class NonlinearCG : public Multivariate<vector_t> {
public:
    using Base = Multivariate<vector_t>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using Base::gradient;

    using scalar_t = typename Base::scalar_t;
    using function_t = typename Base::function_t;

    // constructors
    NonlinearCG(function_t f):Base(f){};

    // compute
    vector_t eval(const vector_t& init=vector_t::Random(), scalar_t e=epsilon) override {
        // 1. initialize
        vector_t xk = init;
        double f0 = function(init); 
        vector_t gk = gradient(init);
        vector_t pk = -gk;

        // 2. loop
        // k represents index
        for(size_t i=0; i<this->iter; i++) {
            if(gk.isZero(1e-4)) break;
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xk, function(xk)));
#endif
            double alpha = this->line_search_inexact(xk, pk, 0.99, 0.5, 3);

            xk = xk + alpha * pk;
            vector_t gk1 = gradient(xk);

            double beta_k = 0.0;
            if constexpr (BetaMethod==nonlinear_cg::Beta::CG_FR) {
                float inv = 1/gk.dot(gk);
                beta_k = (gk1.dot(gk1))*inv;
            } else if constexpr (BetaMethod==nonlinear_cg::Beta::CG_PR) {
                float inv = 1/gk.dot(gk);
                beta_k = (gk1.dot(gk1-gk))*inv;
            } else if constexpr (BetaMethod==nonlinear_cg::Beta::CG_HS) {
                float inv = 1/(gk1-gk).dot(pk);
                beta_k = (gk1.dot(gk1-gk))*inv;
            }

            pk = -gk1 + beta_k * pk;
            gk = gk1;
        }
        return pk;
    };
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__NonlinearCG__