#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

#include <Eigen/Dense>

#include "fwd.h"
#include "method.h"

namespace numerical_optimization {

template<typename VectorTf>
class Multivariate : public Method {
    using function_t = multi::function_t<VectorTf>;

public:
    Multivariate(function_t f):function(f){};

    VectorTf nelder_mead(){ return VectorTf(); };
    VectorTf powells(){};

private:
    function_t function;

    // for nelder_mead method
    void reflecting();
    void expanding();
    void contracting();
};

}

#endif // __MULTIVARIATE_H__