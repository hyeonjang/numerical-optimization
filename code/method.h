#ifndef __METHOD_H__
#define __METHOD_H__

#include <iostream>

#include "fwd.h"

namespace numerical_optimization {

class Method {
public:
    int random_int() const;
protected:
    const size_t max_iter = 10000000; // termination condition
};

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////

#endif // __METHOD_H__