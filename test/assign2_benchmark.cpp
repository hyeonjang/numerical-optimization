#include <functional>
#include <benchmark/benchmark.h>

#include <cmath>

#include <iostream>

#include "method.h"

using namespace numerical_optimization;

class TestFunctions : public benchmark::Fixture
{
public:
    using function_t = std::function<float(const float&)>;

    TestFunctions()
    {
        // making functions for test
        constexpr int number_functions = 4;
        functions.resize(number_functions);

        functions[0] = [](float x){ return std::pow(x, 4)/4 + std::pow(x, 3) + std::pow(x, 2)*(9/2) - 10*x; };
        functions[1] = [](float x){ return std::log(x) + 1; };
        functions[2] = [](float x){ return std::cos(x) +2*x; };
        functions[3] = [](float x){ return (x/(1.4*1.4))*std::exp((-1*x*x)/(2*1.4*1.4)); };

        for(auto func : functions)
            methods.emplace_back(Method(func));
    }

    std::vector<function_t> functions;
    std::vector<Method> methods; 
};

// fibonacci search method
BENCHMARK_F(TestFunctions, fibonacci_search0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].fibonacci_search();
    }
}  

BENCHMARK_F(TestFunctions, fibonacci_search1)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[1].fibonacci_search();
    }
}  

BENCHMARK_F(TestFunctions, fibonacci_search2)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[2].fibonacci_search();
    }
}  

BENCHMARK_F(TestFunctions, fibonacci_search3)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[3].fibonacci_search();
    }
}  

// golden section method
BENCHMARK_F(TestFunctions, golden_section0)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[0].golden_section();
    }
}  

BENCHMARK_F(TestFunctions, golden_section1)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[1].golden_section();
    }
}  

BENCHMARK_F(TestFunctions, golden_section2)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[2].golden_section();
    }
}  

BENCHMARK_F(TestFunctions, golden_section3)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[3].golden_section();
    }
}  

BENCHMARK_MAIN();