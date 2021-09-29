#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

#include "fwd.h"

namespace numerical_optimization {

class Multivariate {
public:
    Multivariate(multi::function_t f):function(f){};

private:
    multi::function_t function;
};

}

#endif // __MULTIVARIATE_H__