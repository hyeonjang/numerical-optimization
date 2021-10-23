#ifndef __QUASI_NEWTONS__
#define __QUASI_NEWTONS__

#include "multivariate.h"

namespace numerical_optimization {
namespace quasi_newtons {
enum Rank { SR1, BFGS, };
};

template<typename VectorTf, quasi_newtons::Rank RankMethod>
class QuasiNewtons : public Multivariate<VectorTf> {
public:
    using Base = Multivariate<VectorTf>;
    using Base::Base;
    using Base::plot;
    using Base::function;
    using function_t = typename Base::function_t;
    using MatrixTf = Eigen::Matrix<typename VectorTf::Scalar, VectorTf::RowsAtCompileTime, VectorTf::RowsAtCompileTime>;

    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h=epsilon) const {
        return Termination::eval<VectorTf, CType>(function, x, h);
    }

    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {

        VectorTf xi = init;
        MatrixTf Hk = MatrixTf::Identity();

        for(size_t i=0; i<10; i++) {

            // termination criterion
            // if(terminate<Termination::Condition::MagnitudeGradient>({xi})) break;

            // Compute a Search Direction
            VectorTf p = -1 * Hk * this->gradient(xi);
            
            // Compute a step length Wolfe Condition
            float alpha = this->line_search_inexact(xi, p, 0.05, 0.05);
            
            // Define sk and yk
            VectorTf Sk = alpha*p;
            VectorTf yk = this->gradient(xi+Sk) - this->gradient(xi);

            // Compute Hk+1
            if constexpr (RankMethod==quasi_newtons::Rank::SR1)
                Hk = SR1(Hk, Sk, yk);
            else if constexpr (RankMethod==quasi_newtons::Rank::BFGS)
                Hk = BFGS(Hk, Sk, yk);

#ifdef BUILD_WITH_PLOTTING
            Log(xi);
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            xi = xi - Hk*this->gradient(xi);
        }
        return xi;
    };

    inline MatrixTf SR1(const MatrixTf& Hk, const VectorTf& Sk, const VectorTf& yk) {
        return Hk + (Sk - Hk*yk) * (Sk - Hk*yk).transpose() / ((Sk - Hk*yk).transpose() * yk);
    }
    inline MatrixTf BFGS(const MatrixTf& Hk, const VectorTf& Sk, const VectorTf& yk) {
        auto pk = 1/(yk.transpose() * Sk);
        return (MatrixTf::Identity()-pk*Sk*yk.transpose())*Hk*(MatrixTf::Identity()-pk*yk*Sk.transpose())+pk*Sk*Sk.transpose();
    }
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__QUASI_NEWTONS__