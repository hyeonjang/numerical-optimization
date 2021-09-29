#ifndef __UNIVARIATE_H__
#define __UNIVARIATE_H__

#include <cassert>
#include <vector>
#include "fwd.h"

namespace numerical_optimization {

class Univariate {
public:
    using function_t = uni::function_t;
    using boundary_t = uni::boundary_t;

    Univariate(function_t f):function(f) { boundary = seeking_bound(5); };
    Univariate(function_t f, boundary_t b):function(f), boundary(b){};

    // assignment 1
    float bisection(float start, float end);
    float newtons(float x);
    float secant(float x1, float x0);
    float regular_falsi(float start, float end);

    // assignment 2
    float fibonacci_search(size_t N=FIBONACCI_MAX);
    float fibonacci_search(float start, float end, size_t N);
    float golden_section(size_t N=FIBONACCI_MAX);
    float golden_section(float start, float end, size_t N);

    // for convienience
    boundary_t get_bound() const;
    Univariate derivate() const;
private:
    function_t function;
    boundary_t boundary;
    const size_t iter = 10000000; // termination condition

    std::vector<int> construct_fibonacci(size_t N) const; // for fibonacci search
    boundary_t seeking_bound(float step_size);
    int random_int() const;
    bool near_zero(float x) { 
        return x==0 || (-MIN<function(x)&&function(x)<MIN); 
    }
};

};

#endif //__UNIVARIATE_H__