#ifndef __LEASTSQUARE_H__
#define __LEASTSQUARE_H__

#include <Eigen/Dense>
#include "multivariate.h"

namespace numerical_optimization {

template<typename VectorT>
class LSM : Method {
public:
    // the coefficients are one more than the variable
    using coefficient_matrix_t = Eigen::Matrix<typename VectorT::Scalar, VectorT::RowsAtCompileTime+1, VectorT::ColsAtCompileTime>;
    using function_t = std::function<double(const VectorT&, coefficient_matrix_t&)>;

    // constructor
    LSM(function_t func):function(func),coefficient(coefficient_matrix_t::Constant(1)){};
    
    // set observation data
    void set_observation(const Eigen::MatrixXd obs){observation=obs;};

        double residual_function(coefficient_matrix_t coeff, VectorT vars, double f) {
        return function(vars, coeff) - f;
    }

    // functions
    // calculate residual vectors
    VectorXd calculate_residual(const coefficient_matrix_t& coef) {

        size_t size = observation.cols();
        VectorXd residue(size);

        for(size_t i=0; i<size; i++) {
            auto col = observation.col(i);

            // more think
            VectorT var { col[0], col[1], col[2] };
            auto val = col[3];
            residue[i] = residual_function(coef, var, val);
        }
        return residue;
    }

    MatrixXd calculate_jacobian(coefficient_matrix_t coef, double eps=1e-6) {

        MatrixXd result(observation.cols(), 4);        
        VectorXd f_residual = calculate_residual(coef);

        for(size_t i=0; i<coefficient_matrix_t::RowsAtCompileTime; i++) {
            coef[i] += eps;
            VectorXd diff = (calculate_residual(coef) - f_residual)/eps;
            result.col(i) = diff;
        }
        return result;
    }
protected:
    function_t function;
    coefficient_matrix_t coefficient;
    Eigen::MatrixXd observation;
};

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif//__LEASTSQUARE_H__