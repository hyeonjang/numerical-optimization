#include <functional>
#include <benchmark/benchmark.h>

#include <cmath>

#include <iostream>

#include "univariate.h"

using namespace numerical_optimization;
using namespace numerical_optimization::uni;

class TestFunctions : public benchmark::Fixture
{
public:
    using function_t = std::function<float(const float&)>;

    TestFunctions()
    {
        int number_functions = 4;
        functions.resize(number_functions);

        functions[0] = [](float x){ return x*x*x + 3*x*x + 9*x - 10; };
        functions[1] = [](float x){ return std::log(x) + 1; };
        functions[2] = [](float x){ return std::cos(x) +2*x; };
        functions[3] = [](float x){ return (x/(1.4*1.4))*std::exp((-1*x*x)/(2*1.4*1.4)); };

        for(auto func : functions)
            methods.emplace_back(Univariate(func));
    }

    std::vector<function_t> functions;
    std::vector<Univariate> methods; 
};

// bisection method
BENCHMARK_F(TestFunctions, bisection0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].bisection(-5, 5);
    }
}       

BENCHMARK_F(TestFunctions, bisection1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].bisection(0.1, 5);
    }
}

BENCHMARK_F(TestFunctions, bisection2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].bisection(-5, 5);
    }
}

BENCHMARK_F(TestFunctions, bisection3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].bisection(-1, 1);
    }
}

// Newtons'method
BENCHMARK_F(TestFunctions, newtons0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].newtons(5);
    }
}       

BENCHMARK_F(TestFunctions, newtons1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].newtons(5);
    }
}

BENCHMARK_F(TestFunctions, newtons2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].newtons(5);
    }
}

BENCHMARK_F(TestFunctions, newtons3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].newtons(1);
    }
}

// secant method
BENCHMARK_F(TestFunctions, secant0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].secant(5, 6);
    }
}       

BENCHMARK_F(TestFunctions, secant1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].secant(5, 6);
    }
}

BENCHMARK_F(TestFunctions, secant2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].secant(5, 6);
    }
}

BENCHMARK_F(TestFunctions, secant3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].secant(1, 2);
    }
}

// regular falsi method
BENCHMARK_F(TestFunctions, regular_falsi0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].regular_falsi(-5, 5);
    }
}       

BENCHMARK_F(TestFunctions, regular_falsi1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].regular_falsi(0.1, 5);
    }
}

BENCHMARK_F(TestFunctions, regular_falsi2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].regular_falsi(-5, 5);
    }
}

BENCHMARK_F(TestFunctions, regular_falsi3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].regular_falsi(-1, 1);
    }
}



BENCHMARK_MAIN();