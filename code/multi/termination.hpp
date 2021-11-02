#ifndef __TERMINATION_HPP__
#define __TERMINATION_HPP__

#include <cassert>
#include "../multivariate.h"

namespace numerical_optimization {
// can it be static? compile time checking functions
namespace Termination {

    static size_t restart = 0;

    enum class Condition : unsigned int {
        None                            = 0,
        ConsecutiveDifference           = 1,
        ConsecutiveDifferenceRelative   = 2,
        MagnitudeGradient               = 4,
        FunctionValueDifferenceRelative = 8,
        DescentDirectionChange          = 16,
    };

    constexpr Condition operator& (Condition lhs, Condition rhs) {
        using T = std::underlying_type_t<Condition>;
        return static_cast<Condition>(static_cast<T>(lhs)&static_cast<T>(rhs));
    }

    constexpr Condition operator| (Condition lhs, Condition rhs) {
        using T = std::underlying_type_t<Condition>;
        return static_cast<Condition>(static_cast<T>(lhs)|static_cast<T>(rhs));
    }

    template <typename Enumeration>
    auto as_integer(Enumeration const value) -> typename std::underlying_type<Enumeration>::type
    {
        return static_cast<typename std::underlying_type<Enumeration>::type>(value);
    }

    // 1. Difference of two consecutive estimates
    template<typename VectorTf>
    inline bool consecutive_difference(const std::vector<VectorTf>& x, float eps) {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size(); // indexing
            flag &= (x[k1]-x[k]).norm()<eps;
        }
        return flag;
    };
    // 2. Relative Difference of two consecutive estimates
    template<typename VectorTf>
    inline bool consecutive_difference_relative(const std::vector<VectorTf>& x, float eps) {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= (x[k1]-x[k]).norm()/x[k1].norm()<eps;
        }
        return flag;
    };
    // 3. Magnitude of Gradient
    template<typename VectorTf>
    inline bool magnitude_gradient(const multi::function_t<VectorTf>& function, const std::vector<VectorTf>& x, float h, float eps=epsilon) {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            flag &= _gradient(function, x[k], eps).norm()<h;
        }
        return flag;
    };
    // 4. Relative Difference of function values
    template<typename VectorTf>
    inline bool function_value_difference_relative(const multi::function_t<VectorTf>& function, const std::vector<VectorTf>& x, float eps) {
        bool flag = true;
        for(size_t k=0; k<x.size(); k++) {
            size_t k1 = (k+1)%x.size();
            flag &= std::abs(function(x[k1])-function(x[k]))/std::abs(function(x[k1])) < eps;
        }
        return flag;
    };
    // 5. Descent direction change
    template<typename VectorTf>
    inline bool descent_direction_change(const multi::function_t<VectorTf>& function, const std::vector<VectorTf>& x, const std::vector<VectorTf>& p) {
        bool flag = true;
        for(size_t k=0; x.size(); k++) {
            flag &= (_gradient(function, x[k], epsilon).transpose()*p[k])>=0.f;
        }
        return flag;
    };
    
    // 6. Maximum number of iterations
    inline bool over_maximum_iteration(size_t iter, size_t max_iter) {
        return iter >= max_iter;
    };

    template<typename VectorTf>
    inline bool check_nan(VectorTf& vec) {
        return vec.hasNaN();
    }

    template<typename VectorTf, Condition CType>
    bool eval(const multi::function_t<VectorTf>& function, const std::vector<VectorTf>& x, float h, float eps=epsilon) {

        if(x[0].hasNaN()){
            // printf("Failed to converge\n");
            return true;
        }

        bool flag = true;
        if constexpr((CType&Condition::ConsecutiveDifference)==Condition::ConsecutiveDifference) {
            // bug in here
            flag &= consecutive_difference(x, eps);
        } else if constexpr((CType&Condition::ConsecutiveDifferenceRelative)==Condition::ConsecutiveDifferenceRelative) {
            flag &= consecutive_difference_relative(x, eps);
        } else if constexpr((CType&Condition::MagnitudeGradient)==Condition::MagnitudeGradient) {
            flag &= magnitude_gradient(function, x, h, eps);
        } else if constexpr((CType&Condition::FunctionValueDifferenceRelative)==Condition::FunctionValueDifferenceRelative) {
            flag &= function_value_difference_relative(function, x, eps);
        } else if constexpr((CType&Condition::DescentDirectionChange)==Condition::DescentDirectionChange) {
            // flag &= descent_direction_change(x, p);
        }
        return flag;
    }
};
/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////
#endif //__TERMINATION_HPP__