#include "common.h"

#include <limits>
#include <functional>

#include <iostream>

float min = std::numeric_limits<float>::min();

float bisection(std::function<float(const float&)> f, float start, float end) {

    if(f(start) * f(end)>=0) 
        return NULL;

    auto midpoint = (start + end)/2.f;

    auto result = f(midpoint);

    if(result=0 || end-start<min)
        return NULL;


    if(f(midpoint)*f(start)<0)
        bisection(f, start, midpoint);
    else
        bisection(f, midpoint, end);
    
    return midpoint;
}
