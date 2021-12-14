#ifndef __LM__
#define __LM__

#include "lsm.h"

namespace numerical_optimization {

template<typename vector_t>
class LM : public LSM<vector_t> {
public:
    using Base = LSM<vector_t>;
    using Base::Base;
    using Base::function;
    using Base::coefficient;
    using Base::calculate_residual;
    using Base::calculate_jacobian;
    using Base::loss;

    using coeff_t = typename Base::coeff_t;
    using scalar_t = typename Base::scalar_t;
    using function_t = std::function<scalar_t(const coeff_t&, const vector_t&)>;

    LM(function_t func):Base(func){};
    LM(function_t func, coeff_t coef):Base(func, coef){};

    // calculate pseudo hessian $((j^t*j + \lambda * I) * J^t)$
    inline MatrixXd calculate_phessian(const MatrixXd& jacobian, scalar_t lambda) {
        MatrixXd jtj = jacobian.transpose() * jacobian;
        return (jtj + lambda * MatrixXd::Identity(jtj.cols(), jtj.rows())).inverse() * jacobian.transpose();
    }

    // check descent condition
    inline bool is_descent(const coeff_t& coefficient, const coeff_t& p) {
        return calculate_residual(coefficient-p).squaredNorm()<=calculate_residual(coefficient).squaredNorm();
    }

    // run fitting
    coeff_t fit(size_t max_iter) {

        for(size_t i=0; i<max_iter; i++) {
            scalar_t lambda = 1;
            
            VectorXd residual = calculate_residual(coefficient);
            MatrixXd jacobian = calculate_jacobian(coefficient);
            MatrixXd phessian = calculate_phessian(jacobian, lambda)*residual;

            if(is_descent(coefficient, phessian)) {
                // just one division
                lambda /= 10;
                phessian = calculate_phessian(jacobian, lambda) * residual;
            } else {
                while(!is_descent(coefficient, phessian)) {
                    lambda *= 10;
                    phessian = calculate_phessian(jacobian, lambda) * residual;
                }
            }
            coefficient = coefficient - phessian;
            if(loss(coefficient)<1e-4) break;
        }
        return coefficient;
    }

};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__LM__