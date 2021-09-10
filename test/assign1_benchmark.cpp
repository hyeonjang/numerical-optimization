#include <functional>
#include <benchmark/benchmark.h>

#include "method.h"

static void BenchMarkFunction(benchmark::State& state) 
{
    std::function<float(const float&)> function = [](float x) { return x*x/16 - 4; };
    numerical_optimization::Method method(function);

    std::pair<float, float> interval = {0, 10};

    for (auto _: state)
    {
        method.bisection(interval.first, interval.second);
    }   

}


BENCHMARK(BenchMarkFunction);

BENCHMARK_MAIN();