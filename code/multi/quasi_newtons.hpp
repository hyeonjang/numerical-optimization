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

template<typename vector_t, quasi_newtons::Rank RankMethod>
class QuasiNewtons : public Multivariate<vector_t> {
public:
    using Base = Multivariate<vector_t>;
    using Base::Base;
    using Base::plot;
    using Base::iter;
    using Base::function;
    using Base::gradient;
    using scalar_t = typename vector_t::Scalar;
    using matrix_t = Eigen::Matrix<typename vector_t::Scalar, vector_t::RowsAtCompileTime, vector_t::RowsAtCompileTime>;
    using function_t = typename Base::function_t;

    template<Termination::Condition CType> 
    bool terminate(const std::vector<vector_t>& x, scalar_t h, scalar_t eps=epsilon) {
        return Termination::eval<CType, vector_t, scalar_t>(function, x, h, eps);
    }
    vector_t eval(const vector_t& init=vector_t::Random(), float e=epsilon) override {

        vector_t xi = init;
        matrix_t Hk = matrix_t::Identity();

        size_t iteration = 0;
        for(size_t i=0; i<this->iter; i++) {

#ifdef BUILD_WITH_PLOTTING
            plot.emplace_back(std::make_pair(xi, function(xi)));
#endif
            // Compute a Search Direction
            vector_t p = (-1*Hk*gradient(xi)).normalized();

            // Compute a step length Wolfe Condition
            double alpha = 0;
            if constexpr (RankMethod==quasi_newtons::Rank::SR1)
                alpha = this->line_search_inexact(xi, p, 0.99, 0.5, 3);
            else if constexpr (RankMethod==quasi_newtons::Rank::BFGS)
                alpha = this->line_search_inexact(xi, p, 0.8, 0.5, 3);

            // Define sk and yk
            vector_t Sk = alpha*p;
            vector_t yk = gradient(xi+Sk) - gradient(xi);

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

    inline matrix_t SR1(const matrix_t& Hk, const vector_t& Sk, const vector_t& yk) {
        auto frac = 1/((Sk - Hk*yk).transpose()*yk);
        return Hk + ((Sk-Hk*yk) * (Sk-Hk*yk).transpose())*frac;
    }
    inline matrix_t BFGS(const matrix_t& Hk, const vector_t& Sk, const vector_t& yk) {
        auto pk = 1/(yk.transpose() * Sk);
        return (matrix_t::Identity()-pk*Sk*yk.transpose())*Hk*(matrix_t::Identity()-pk*yk*Sk.transpose())+pk*Sk*Sk.transpose();
    }
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__QUASI_NEWTONS__