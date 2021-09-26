#include <functional>
#include <benchmark/benchmark.h>
#include <cmath>
#include <iostream>
#include "method.h"

#define print_bound(str, m){ std::cout << str << ": [" << m.get_bound().first << ", " << m.get_bound().second << "]"<< std::endl;}

using namespace numerical_optimization;

class TestFunctions : public benchmark::Fixture {
public:
    using function_t = std::function<float(const float&)>;

    TestFunctions() {
        // making functions for test
        constexpr int number_functions = 6;
        functions.resize(number_functions);
        
        // assignment 1 cases
        functions[0] = [](float x){ return std::pow(x, 4)/4 + std::pow(x, 3) + std::pow(x, 2)*(9/2) - 10*x; };
        functions[1] = [](float x){ return std::log(x)*x; };
        functions[2] = [](float x){ return std::sin(x) +x*x-10; };
        functions[3] = [](float x){ return -std::exp((-1*x*x)/(1.4*1.4)); };

        // not differentiable cases
        functions[4] = [](float x){ return std::abs(x-0.3); };
        functions[5] = [](float x){ return std::abs(std::log(x)); };

        for(auto func : functions)
            methods.emplace_back(Method(func));

        // for(auto m : methods) {
        //     auto b =  m.get_bound();
        //     std::cout << "[ " << b.first << ", " << b.second << "]" << std::endl;
        // }
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

BENCHMARK_F(TestFunctions, fibonacci_search4)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[4].fibonacci_search();
    }
}  

BENCHMARK_F(TestFunctions, fibonacci_search5)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[5].fibonacci_search();
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

BENCHMARK_F(TestFunctions, golden_section4)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[4].golden_section();
    }
}  

BENCHMARK_F(TestFunctions, golden_section5)(benchmark::State& st)
{
    for( auto _: st ) 
    {
        methods[5].golden_section();
    }
}  

BENCHMARK_MAIN();