#ifndef __METHOD_H__
#define __METHOD_H__

#include <iostream>
#include <random>
#include <algorithm>

#include "fwd.h"

namespace numerical_optimization {

class Method {
public:
    Method():iter(max_iter){};
    Method(size_t it):iter(it){};

protected:
    size_t iter;

public: 
    template<typename Type>
    Type random() const {
    // threshold
    constexpr int scale = 100000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(
        std::numeric_limits<Type>::min()/scale, 
        std::numeric_limits<Type>::max()/scale
        );
    return distrib(gen);
}

};

/////////////////////////////////////////////////////
} /// the end of namespace numerical_optimization ///
/////////////////////////////////////////////////////

#endif // __METHOD_H__