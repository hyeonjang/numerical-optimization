#include <cmath>
#include <iostream>
#include <functional>
#include <benchmark/benchmark.h>

#include "function.cpp"

using namespace numerical_optimization;
using namespace numerical_optimization::uni;

extern std::vector<Univariate> numerical_optimization::uni::methods;

static void fibonacci_search(benchmark::State& state) {
    for(auto _ : state) {
        methods[state.range(0)].fibonacci_search();
    }
}

static void golden_section_search(benchmark::State& state) {
    for(auto _ : state) {
        methods[state.range(0)].golden_section();
    }
}

static void bisection(benchmark::State& state) {
    for(auto _ : state) {
       Univariate m = methods[state.range(0)].derivate();
        // m.bisection();
    }
}

static void print_interval() {
    for(const auto& m : methods) {
        auto b = m.get_bound();
        std::cout << "[" << b.first << ", " << b.second << "]" << std::endl;
    }
}

BENCHMARK(fibonacci_search)->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4);
BENCHMARK(golden_section_search)->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4);
// BENCHMARK(bisection)->Arg(0)->Arg(1)->Arg(2);

BENCHMARK_MAIN();