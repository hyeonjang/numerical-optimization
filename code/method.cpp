#include "method.h"
#include <iostream>

namespace numerical_optimization
{

float Method::bisection(float start, float end)
{
    assert( function(start)*function(end)<0 );

    auto midpoint = (start + end)/2.f;

    if(function(midpoint)==0 || end-start<MIN)
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
    while(function(x1)>=0.f)
    {
        float t = x1;
        x1 = t - function(t)/d(function, t);
    }

    return x1;
}

// Two point approximation method
float Method::secant(float x1, float x0)
{
    // no matter which one is bigger
    float t1 = std::min(x1, x0);
    float t0 = std::max(x1, x0);

    // initial two points
    float x2 = MAX;
    while(function(x2)>0.f)
    {
        x2 = t1 - ((t1-t0)/(function(t1)-function(t0))) * function(t1);

        t0 = t1;
        t1 = x2;
    }

    return x2;
}

float Method::regular_falsi(float start, float end)
{
    // secant method lambda
    auto sec = [](std::function<float(const float&)> func, float x1, float x0)
    { 
        return x1 - ((x1-x0)/(func(x1)-func(x0))) * func(x1); 
    };

    assert( function(start)*function(end)<0 );

    // new x-axis intersection point
    float x = sec(function, start, end);

    if ( end-start<MIN )
        return x;


    if( function(x)==0 || -MIN<function(x) && function(x)<MIN ) // almost zero
        return x;

    // do recursivly until the end
    if( function(start) * function(x) < 0)
        x = regular_falsi(start, x);
    else if ( function(end) * function(x) < 0)
        x = regular_falsi(x, end);

    return x;
}

float Method::regular_falsi_not_recur(float start, float end)
{
    // secant method lambda
    auto sec = [](std::function<float(const float&)> func, float x1, float x0)
    { 
        return x1 - ((x1-x0)/(func(x1)-func(x0))) * func(x1); 
    };

    assert( function(start)*function(end)<0 );

    float x0 = start;
    float x1 = end;

    // new x-axis intersection point
    float x2 = sec(function, x0, x1);

    while( !near_zero(x2) )
    {
        if ( function(x0) * function(x2) < 0)
            x2 = sec(function, x0, x2);
        else if( function(x2) * function(x1) < 0.f )
            x2 = sec(function, x2, x1);
    }
    return x2;
}
}
