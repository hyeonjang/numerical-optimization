#include <cmath>
#include <iostream>
#include <functional>
#include <benchmark/benchmark.h>

#include "multi/cauchys.hpp"
#include "multi/newtons.hpp"
#include "multi/quasi_newtons.hpp"
#include "../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;
using namespace numerical_optimization::multi;

// input functions
std::vector<function_t<Vector2f>> functions = construct_functions();

// optimizing methods
auto cauchys = construct_methods<Cauchys<Vector2f>>(functions);
auto newtons = construct_methods<Newtons<Vector2f>>(functions);
auto SR1     = construct_methods<QuasiNewtons<Vector2f, quasi_newtons::SR1>>(functions);
auto BFGS1   = construct_methods<QuasiNewtons<Vector2f, quasi_newtons::BFGS>>(functions);

template <typename Method>
static void bench_1(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2f(1.2, 1.2));
    }
}

template <typename Method>
static void bench_2(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2f(5.6, -1.2));
    }
}

template <typename Method>
static void bench_3(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2f(-3.5, 2.3));
    }
}

template <typename Method>
static void bench_3(benchmark::State& state) {
    auto method = construct_methods<Method>(functions);
    for(auto _ : state) {
        method[state.range(0)].eval(Vector2f(10.5, -8.3));
    }
}

// failed to converge in Newtons and QuasiNewtons in the third function
BENCHMARK(bench_1<Cauchys<Vector2f>>)->Range(0, 2);
BENCHMARK(bench_1<Newtons<Vector2f>>)->Range(0, 2);
BENCHMARK(bench_1<QuasiNewtons<Vector2f, quasi_newtons::SR1>>)->Range(0, 2);
BENCHMARK(bench_1<QuasiNewtons<Vector2f, quasi_newtons::BFGS>>)->Range(0, 2);

BENCHMARK(bench_2<Cauchys<Vector2f>>)->Range(0, 2);
BENCHMARK(bench_2<Newtons<Vector2f>>)->Range(0, 2);
BENCHMARK(bench_2<QuasiNewtons<Vector2f, quasi_newtons::SR1>>)->Range(0, 2);
BENCHMARK(bench_2<QuasiNewtons<Vector2f, quasi_newtons::BFGS>>)->Range(0, 2);

// BENCHMARK(bench_3<Cauchys<Vector2f>>)->Range(0, 2);
// BENCHMARK(bench_3<Newtons<Vector2f>>)->Range(0, 2);
// BENCHMARK(bench_3<QuasiNewtons<Vector2f, quasi_newtons::SR1>>)->Range(0, 2);
// BENCHMARK(bench_3<QuasiNewtons<Vector2f, quasi_newtons::BFGS>>)->Range(0, 2);

// BENCHMARK(bench_4<Cauchys<Vector2f>>)->Range(0, 2);
// BENCHMARK(bench_4<Newtons<Vector2f>>)->Range(0, 2);
// BENCHMARK(bench_4<QuasiNewtons<Vector2f, quasi_newtons::SR1>>)->Range(0, 2);
// BENCHMARK(bench_4<QuasiNewtons<Vector2f, quasi_newtons::BFGS>>)->Range(0, 2);

BENCHMARK_MAIN();