#include <functional>
#include <benchmark/benchmark.h>

#include <iostream>

#include "method.h"

using namespace numerical_optimization;

class TestFunctions : public benchmark::Fixture
{
public:
    using function_t = std::function<float(const float&)>;

    TestFunctions()
    {
        int number_functions = 4;
        functions.resize(number_functions);

        functions[0] = [](float x){ return x*x - 10; };
        functions[1] = [](float x){ return x*x - 10; };
        functions[2] = [](float x){ return x*x - 10; };
        functions[3] = [](float x){ return x*x - 10; };

        for(auto func : functions)
            methods.emplace_back(Method(func));
    }

    std::vector<function_t> functions;
    std::vector<Method> methods; 
};

// bisection method
BENCHMARK_F(TestFunctions, bisection0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].bisection(-1, 10);
    }
}       

BENCHMARK_F(TestFunctions, bisection1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].bisection(-2, 10);
    }
}

BENCHMARK_F(TestFunctions, bisection2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].bisection(-2, 10);
    }
}

BENCHMARK_F(TestFunctions, bisection3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].bisection(-2, 10);
    }
}

// Newtons'method
BENCHMARK_F(TestFunctions, newtons0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].newtons(10);
    }
}       

BENCHMARK_F(TestFunctions, newtons1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].newtons(10);
    }
}

BENCHMARK_F(TestFunctions, newtons2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].newtons(-2);
    }
}

BENCHMARK_F(TestFunctions, newtons3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].newtons(10);
    }
}

// secant method
BENCHMARK_F(TestFunctions, secant0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].secant(-1, 10);
    }
}       

BENCHMARK_F(TestFunctions, secant1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].secant(-2, 10);
    }
}

BENCHMARK_F(TestFunctions, secant2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].secant(-2, 10);
    }
}

BENCHMARK_F(TestFunctions, secant3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].secant(-2, 10);
    }
}

// regular falsi method
BENCHMARK_F(TestFunctions, regular_falsi0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].regular_falsi(-1, 10);
    }
}       

BENCHMARK_F(TestFunctions, regular_falsi1)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[1].regular_falsi(-2, 10);
    }
}

BENCHMARK_F(TestFunctions, regular_falsi2)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[2].regular_falsi(-2, 10);
    }
}

BENCHMARK_F(TestFunctions, regular_falsi3)(benchmark::State& st)
{
    for( auto _: st)
    {
        methods[3].regular_falsi(-2, 10);
    }
}



BENCHMARK_MAIN();