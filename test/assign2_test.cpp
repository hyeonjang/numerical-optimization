#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "function.cpp"

using namespace numerical_optimization;
using namespace numerical_optimization::uni;

extern std::vector<Univariate> numerical_optimization::uni::methods;

TEST(FibonacciTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Univariate method(function);

    auto x0 = method.fibonacci_search(-4, 3, 20);
    EXPECT_NEAR(x0, 0.f, 1e-3);

    auto x1 = method.fibonacci_search();
    EXPECT_NEAR(x1, 0.f, 1e-3);
}

TEST(GoldenSectionTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Univariate method(function);

    auto x0 = method.golden_section(-3, 2.5, 8);
    EXPECT_NEAR(x0, 0.f, 0.1f);

    auto x1 = method.golden_section();
    EXPECT_NEAR(x1, 0.f, 1e-03);
}