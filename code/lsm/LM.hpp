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
    using function_t = std::function<double(const coeff_t&, const vector_t&)>;

    LM(function_t func):Base(func){};
    LM(function_t func, coeff_t coef):Base(func, coef){};

    // calculate pseudo hessian $((j^t*j + \lambda * I) * J^t)$
    inline MatrixXd calculate_hessian(const MatrixXd& jacobian, double lambda) {
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
            double lambda = 1;
            
            VectorXd residual = calculate_residual(coefficient);
            MatrixXd jacobian = calculate_jacobian(coefficient);
            MatrixXd phessian = calculate_hessian(jacobian, lambda);
            VectorXd p = phessian * residual;

            if(is_descent(coefficient, p)) {
                // just one division
                lambda /= 10;
                p = calculate_hessian(jacobian, lambda) * residual;
            } else {
                while(!is_descent(coefficient, p)) {
                    lambda *= 10;
                    p = calculate_hessian(jacobian, lambda) * residual;
                }
            }
            coefficient = coefficient - p;
            if(loss(coefficient)<1e-4) break;
        }
        return coefficient;
    }

};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__LM__