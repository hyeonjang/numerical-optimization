#include "method.h"
#include <iostream>
#include <random>

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
    auto d = [](function_t func, float x, float eps=1e-6)
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
    auto sec = [](function_t func, float x1, float x0)
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
    auto sec = [](function_t func, float x1, float x0)
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

// assignment 2 @@todo optimize
float Method::fibonacci_search(float start, float end, size_t N)
{
    std::vector<int> F = construct_fibonacci(N);

    float a = start;
    float b = end;

    N = F.size()-1; // indexing
    for(size_t n=N; n>1; n--)
    {
        float length = b - a;

        float x1 = a + (float(F[n-2])/float(F[n]))*length;
        float x2 = b - (float(F[n-2])/float(F[n]))*length;

        if(function(x1)>function(x2)) 
            a = x1;
        if(function(x1)<function(x2)) 
            b = x2;
    
    }
    return (a + b)/2;
}

float Method::fibonacci_search()
{
    return fibonacci_search(boundary.first, boundary.second, iter);
}

// @@todo optimize
float Method::golden_section(float start, float end, size_t N)
{
    float a = start;
    float b = end;

    for(size_t n=N; n>1; n--)
    {
        float length = b - a;

        float x1 = a + GOLDEN_RATIO * length;
        float x2 = b - GOLDEN_RATIO * length;

        float mx = std::max(x1, x2);
        float mn = std::min(x1, x2);

        if(function(mn)>function(mx))
        {   
            a = mn;
        }
        else if(function(mn)<function(mx))
        {
            b = mx;
        }
    }
    return (a + b)/2;
}

float Method::golden_section()
{
    return golden_section(boundary.first, boundary.second, iter);
}

std::vector<int> Method::construct_fibonacci(size_t N) const
{
    constexpr size_t fibonacci_max = 46;

    N = std::min(N, fibonacci_max);
    std::vector<int> fibonacci(N);

    fibonacci[0] = 1;
    fibonacci[1] = 1;

    for(size_t i=0; i<N-2; i++)
        fibonacci[i+2] = fibonacci[i] + fibonacci[i+1];

    return fibonacci;
}

std::pair<float, float> Method::seeking_bound(float step_size)
{
    std::pair<float, float> result;

    std::vector<float> x(iter); x[1] = float(random_int());
    float d = step_size;

    float f0 = function(x[1]-d);
    float f1 = function(x[1]);
    float f2 = function(x[1]+d);

    if (f0>=f1 && f1>=f2)
    {
        x[0] = x[1] - d;
        x[2] = x[1] + d;
        d = d;
    }
    else if (f0<=f1 && f1<=f2)
    {
        x[0] = x[1] + d;
        x[2] = x[1] - d;
        d = -d;
    }
    else if (f0>=f1 && f1<=f2)
        result = std::make_pair(x[1]-d, x[1]+d);

    // now default
    function_t increment = [](const float& f){ return std::pow(2, f); };
    for(size_t k=2; k<iter; k++)
    {
        x[k+1] = x[k] + increment(k) * d;

        if(function(x[k+1])>=function(x[k]) && d>0)
        {
            result = std::make_pair(x[k-1], x[k+1]);
            break;
        }

        if(function(x[k+1])>=function(x[k]) && d<0)
        {
            result = std::make_pair(x[k+1], x[k-1]);
            break;
        }
    }
    return result;
}

int Method::random_int() const
{
    constexpr int scale = 100000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(
        std::numeric_limits<int>::min()/scale, 
        std::numeric_limits<int>::max()/scale
        );
    return distrib(gen);
}

}
