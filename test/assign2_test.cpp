#include <functional>
#include <gtest/gtest.h>
#include <cmath>

#include "method.h"

using namespace numerical_optimization;

TEST(MethodsTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Method method(function);

    auto x = method.bisection(0, 128);
    EXPECT_EQ(std::roundf(x), 2);
}


TEST(FibonacciTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Method method(function);

    auto x = method.fibonacci_search(-4, 3, 10);
    EXPECT_EQ(std::roundf(x), 0);
}

TEST(GoldenSectionTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Method method(function);

    auto x = method.golden_section(-3, 2.5, 6);
    EXPECT_EQ(std::roundf(x), 0);
}