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

// Two point approximation method
float Method::secant(float x1, float x0)
{
    x1 = std::min(x1, x0);
    x0 = std::max(x1, x0);

    // initial two points
    float x2 = std::numeric_limits<float>::max();
    while(function(x2)>0.f)
    {
        x2 = x1 - ((x1-x0)/(function(x1)-function(x0))) * function(x1);

        x0 = x1;
        x1 = x2;

        std::cout << "x: " << x2 << std::endl;
        std::cout << "value: " << function(x2) << std::endl;
    }

    return x2;
}

float Method::regular_falsi(float start, float end)
{
    auto sec = [](std::function<float(const float&)> func, float x1, float x0)
    { 
        return x1 - ((x1-x0)/(func(x1)-func(x0))) * func(x1); 
    };

    assert( function(start)*function(end)<0 );

    float x = sec(function, start, end);
    std::cout << "Value: " << function(x) << std::endl;
    std::cout << "x: " << x << std::endl;

    if(function(x)==0 || end-start<MIN)
        return x;

    if(function(start)*function(x) < 0)
        x = regular_falsi(start, x);
    else
        x = regular_falsi(x, end);

    return x;
}

}
