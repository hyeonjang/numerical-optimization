#ifndef __QUASI_NEWTONS__
#define __QUASI_NEWTONS__

#include <math.h>
#include <cassert>
#include "multivariate.h"
#include "termination.h"

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
    using Base::iter;
    using Base::function;
    using Base::gradient;
    using function_t = typename Base::function_t;
    using MatrixTf = Eigen::Matrix<typename VectorTf::Scalar, VectorTf::RowsAtCompileTime, VectorTf::RowsAtCompileTime>;

    template<Termination::Condition CType> 
    bool terminate(const std::vector<VectorTf>& x, float h, float eps=epsilon) {
        return Termination::eval<VectorTf, CType>(function, x, h, eps);
    }
    VectorTf eval(const VectorTf& init=VectorTf::Random(), float e=epsilon) override {

        VectorTf xi = init;
        MatrixTf Hk = MatrixTf::Identity();

        size_t iteration = 0;
        for(size_t i=0; i<this->iter; i++) {

#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            // Compute a Search Direction
            VectorTf p = (-1*Hk*gradient(xi)).normalized();

            // Compute a step length Wolfe Condition
            double alpha = 0;
            if constexpr (RankMethod==quasi_newtons::Rank::SR1)
                alpha = this->line_search_inexact(xi, p, 0.99, 0.5, 3);
            else if constexpr (RankMethod==quasi_newtons::Rank::BFGS)
                alpha = this->line_search_inexact(xi, p, 0.8, 0.5, 3);

            // Define sk and yk
            VectorTf Sk = alpha*p;
            VectorTf yk = gradient(xi+Sk) - gradient(xi);

            // Compute Hk+1
            if constexpr (RankMethod==quasi_newtons::Rank::SR1)
                Hk = SR1(Hk, Sk, yk);
            else if constexpr (RankMethod==quasi_newtons::Rank::BFGS)
                Hk = BFGS(Hk, Sk, yk);

            xi = xi - Hk*gradient(xi);
            
            if constexpr (RankMethod==quasi_newtons::Rank::SR1) {
                if(terminate<Termination::Condition::MagnitudeGradient>({xi}, 0.01)) break;
            }
            else if constexpr (RankMethod==quasi_newtons::Rank::BFGS) {
                if(terminate<Termination::Condition::MagnitudeGradient>({xi}, 1e-8)) break;
            }
        }
        return xi;
    };

    inline MatrixTf SR1(const MatrixTf& Hk, const VectorTf& Sk, const VectorTf& yk) {
        auto frac = 1/((Sk - Hk*yk).transpose()*yk);
        return Hk + ((Sk-Hk*yk) * (Sk-Hk*yk).transpose())*frac;
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