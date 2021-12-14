#include <cmath>
#include <iostream>
#include <functional>
#include <benchmark/benchmark.h>

#include "multi/powells.hpp"
#include "multi/nelder_mead.hpp"
#include "../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

auto functions  = hw2::construct_functions();
auto neldermead = hw2::construct_methods<NelderMead<Vector2d>>(functions);
auto powells    = hw2::construct_methods<Powells<Vector2d>>(functions);

static void bench_neldermead(benchmark::State& state) {
    for(auto _ : state) {
        neldermead[state.range(0)].eval();
    }
}

static void bench_powells(benchmark::State& state) {
    for(auto _ : state) {
        powells[state.range(0)].eval();
    }
}


BENCHMARK(bench_neldermead)->Arg(0)->Arg(1)->Arg(2);
BENCHMARK(bench_powells)->Arg(0)->Arg(1)->Arg(2);

BENCHMARK_MAIN();