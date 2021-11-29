#ifndef __LM__
#define __LM__

#include "lsm.h"

namespace numerical_optimization {

template<typename VectorT>
class LM : public LSM<VectorT> {
public:
    using Base = LSM<VectorT>;
    using Base::Base;
    using Base::function;
    using Base::coefficient;
    using Base::observation;
    using Base::calculate_residual;
    using Base::calculate_jacobian;

    using coefficient_matrix_t = typename Base::coefficient_matrix_t;
    using function_t = std::function<double(const VectorT&, coefficient_matrix_t&)>;

    LM(function_t func):Base(func){};

    coefficient_matrix_t fit() {
        
        for(size_t i=0; i<50; i++) {
            
            VectorXd residual = calculate_residual(coefficient);
            MatrixXd jacobian = calculate_jacobian(coefficient);
            
            MatrixXd jtj = jacobian.transpose() * jacobian;
            MatrixXd jtj_inv_jt = (jtj + MatrixXd::Identity(jtj.cols(), jtj.rows())).inverse() * jacobian.transpose();

            auto p = jtj_inv_jt * residual;
            coefficient = coefficient - p;
        }

        return coefficient;
    }

};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__LM__