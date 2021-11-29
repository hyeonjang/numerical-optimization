#ifndef __LinearCG__
#define __LinearCG__

#include "multivariate.h"

namespace numerical_optimization {

template<typename VectorTf>
class LinearCG : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using Base::gradient;
    using function_t = typename Base::function_t;
    using MatrixTf = Eigen::Matrix<typename VectorTf::Scalar, VectorTf::RowsAtCompileTime, VectorTf::RowsAtCompileTime>;

    // constructor
    LinearCG(function_t f, MatrixTf A, VectorTf b):Base(f),A(A),b(b){};

    // compute
    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {
        VectorTf xk = init;
        VectorTf rk = A*xk - b;
        VectorTf pk = -rk;

        for(size_t i=0; i<this->iter; i++) {
            if(rk.isZero(1e-2)) break;
#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xk, function(xk)));
#endif
            double inv = 1/(pk.transpose()*A*pk);
            double alpha_k = (rk.transpose()*rk);
            alpha_k *= inv;

            xk = xk + alpha_k*pk;
            VectorTf rk1 = rk + alpha_k*A*pk;
            double invv = 1/(rk.transpose()*rk);
            double beta_k1 = (rk1.transpose()*rk1);
            beta_k1 *= inv;

            pk = -rk1 + beta_k1*pk;
            rk = rk1;
        }
        return pk;
    };
private:
    MatrixTf A;
    VectorTf b;
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__LinearCG__