#include <cmath>
#include <iostream>
#include <functional>
#include <benchmark/benchmark.h>

#include "multi/linearcg.hpp"
#include "multi/nonlinearcg.hpp"
#include "../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;
using namespace numerical_optimization::multi;

// input functions
std::vector<function_t<Vector2d>> functions = hw5::construct_functions<Vector2d>();

template <typename Method>
static void bench_1(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2d(1.2, 1.2));
    }
}

template <typename Method>
static void bench_2(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2d(5.6, 1.2));
    }
}

template <typename Method>
static void bench_3(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2d(-3.5, 2.3));
    }
}

template <typename Method>
static void bench_4(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2d(10.5, -8.3));
    }
}

// failed to converge in Newtons and QuasiNewtons in the third function
BENCHMARK(bench_1<NonlinearCG<Vector2d, nonlinear_cg::Beta::CG_FR>>)->Range(1, 2);
BENCHMARK(bench_1<NonlinearCG<Vector2d, nonlinear_cg::Beta::CG_PR>>)->Range(1, 2);
BENCHMARK(bench_1<NonlinearCG<Vector2d, nonlinear_cg::Beta::CG_HS>>)->Range(1, 2);

BENCHMARK_MAIN();