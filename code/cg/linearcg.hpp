#ifndef __LinearCG__
#define __LinearCG__

#include "multivariate.h"

namespace numerical_optimization {

template<typename vector_t>
class LinearCG : public Multivariate<vector_t> {
public:
    using Base = Multivariate<vector_t>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using Base::gradient;

    using scalar_t = typename Base::scalar_t;
    using matrix_t = typename Base::matrix_t;
    using function_t = typename Base::function_t;

    // constructor
    LinearCG(function_t f, matrix_t A, vector_t b):Base(f),A(A),b(b){};

    // compute
    vector_t eval(const vector_t& init=vector_t::Random(), scalar_t e=epsilon) override {
        vector_t xk = init;
        vector_t rk = A*xk - b;
        vector_t pk = -rk;

        for(size_t i=0; i<this->iter; i++) {
            if(rk.isZero(1e-2)) break;
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xk, function(xk)));
#endif
            double inv = 1/(pk.transpose()*A*pk);
            double alpha_k = (rk.transpose()*rk);
            alpha_k *= inv;

            xk = xk + alpha_k*pk;
            vector_t rk1 = rk + alpha_k*A*pk;
            double invv = 1/(rk.transpose()*rk);
            double beta_k1 = (rk1.transpose()*rk1);
            beta_k1 *= inv;

            pk = -rk1 + beta_k1*pk;
            rk = rk1;
        }
        return pk;
    };
private:
    matrix_t A;
    vector_t b;
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__LinearCG__