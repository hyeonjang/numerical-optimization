#ifndef __GAUSSNEWTONS__
#define __GAUSSNEWTONS__

#include "lsm.h"

namespace numerical_optimization {

template<typename vector_t>
class GaussNewtons : public LSM<vector_t> {
public:
    using Base = LSM<vector_t>;
    using Base::Base;
    using Base::function;
    using Base::coefficient;
    using Base::calculate_residual;
    using Base::calculate_jacobian;
    using Base::loss;

    using coeff_t = typename Base::coeff_t;
    using function_t = std::function<double(const coeff_t&, const vector_t&)>;

    GaussNewtons(function_t func):Base(func){};
    GaussNewtons(function_t func, coeff_t coef):Base(func, coef){};

    coeff_t fit(size_t max_iter) {

        for(size_t i=0; i<max_iter; i++) {
            
            VectorXd residual = calculate_residual(coefficient);
            MatrixXd jacobian = calculate_jacobian(coefficient);
            MatrixXd pinv = jacobian.completeOrthogonalDecomposition().pseudoInverse();

            auto p = -pinv*residual;
            coefficient = coefficient + p;
            if(loss(coefficient)<1e-4) break;
        }

        return coefficient;
    }

};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__GAUSSNEWTONS__