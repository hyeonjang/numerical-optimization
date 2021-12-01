#ifndef __LEASTSQUARE_H__
#define __LEASTSQUARE_H__

#include <Eigen/Dense>
#include "multivariate.h"

namespace numerical_optimization {

template<typename vector_t>
class LSM : Method {
public:
    // the coefficients are one more than the variable
    using coeff_t = Eigen::Matrix<typename vector_t::Scalar, vector_t::RowsAtCompileTime+1, vector_t::ColsAtCompileTime>;
    using function_t = std::function<double(const coeff_t&, const vector_t&)>;

    // constructor
    LSM(function_t func):function(func),coefficient(coeff_t::Constant(5)){};
    LSM(function_t func, coeff_t coef):function(func),coefficient(coef){};
    
    // set observation data
    void set_observation(const std::vector<vector_t>& obs_x, const std::vector<double> obs_f){ observation_x=obs_x;observation_f=obs_f; };

    // get redisidual function
    double residual_function(const coeff_t& coeff, const vector_t& vars, double f) {
        return function(coeff, vars) - f;
    }

    // calculate residual vectors
    VectorXd calculate_residual(const coeff_t& coef) {

        size_t size = observation_x.size();
        VectorXd residue(size);

        for(size_t i=0; i<size; i++) {
            auto var = observation_x[i];
            auto val = observation_f[i];

            residue[i] = residual_function(coef, var, val);
        }
        return residue;
    }

    // calculate jacobian matrix
    MatrixXd calculate_jacobian(coeff_t coef, double eps=1e-6) {

        const size_t     cols = observation_x.size();
        constexpr size_t rows = coeff_t::RowsAtCompileTime;
        MatrixXd result(cols, rows);        

        for(size_t i=0; i<coeff_t::RowsAtCompileTime; i++) {
            coeff_t coef_1=coef, coef_2=coef;
            
            coef_1[i]+=eps; coef_2[i]-=eps;
            VectorXd diff = (calculate_residual(coef_1) - calculate_residual(coef_2))/(2*eps);
            result.col(i) = diff;
        }
        return result;
    }

    // loss function to terminate
    double loss(const coeff_t& coeff) {
        return calculate_residual(coeff).squaredNorm()/2;
    }

protected:
    function_t function;
    coeff_t coefficient;
    std::vector<vector_t> observation_x;
    std::vector<double>  observation_f;
};

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif//__LEASTSQUARE_H__