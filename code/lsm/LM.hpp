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

    coefficient_matrix_t eval() {
        
        for(size_t i=0; i<50; i++) {
            
            VectorXd residual = calculate_residual(coefficient);
            MatrixXd jacobian = calculate_jacobian(coefficient);
            MatrixXd tmp = jacobian.transpose() * jacobian();

            auto p = pinv*residual;
            coefficient = coefficient - p;
        }

        return coefficient;
    }

};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__LM__