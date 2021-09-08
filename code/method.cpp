#include "method.h"
#include <iostream>

namespace numerical_optimization
{

float Method::bisection(float start, float end)
{
    assert( function(start)*function(end)<0 );
    std::cout << "interval: [" << start << ", " << end << "]" << std::endl;
    std::cout << "function: [" << function(start) << ", " << function(end) << "]" << std::endl;

    auto midpoint = (start + end)/2.f;
    auto result = function(midpoint);

    std::cout << "midpoint: " << midpoint << std::endl;
    std::cout << "result: " << result << std::endl;

    if(result==0 || end-start<MIN)
        return midpoint;

    if(function(midpoint)*function(start)<0)
        midpoint = bisection(start, midpoint);
    else
        midpoint = bisection(midpoint, end);
    
    return midpoint;
}

float Method::newtons(float x0)
{
    auto d = [](std::function<float(const float&)> func, float x, float eps=1e-6)
    { 
        return (func(x+eps) - func(x))/eps;
    };

    float x1 = x0;
    while(function(x1)>0.f)
    {
        float t = x1;

        x1 = t - function(t)/d(function, t);
        std::cout << "x: " << x1 << std::endl;
        std::cout << "value: " << function(x1) << std::endl;
    }

    return x1;
}
