#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

#include "fwd.h"

namespace numerical_optimization {

class Multivariate {
    using multivar_t = multi::variate_t; 
    using function_t = multi::function_t;

public:
    Multivariate(function_t f):function(f){};

    multivar_t nelder_mead(multivar_t vars);
    multivar_t powells();

private:
    function_t function;

    // for nelder_mead method
    void reflecting();
    void expanding();
    void contracting();
};

}

#endif // __MULTIVARIATE_H__