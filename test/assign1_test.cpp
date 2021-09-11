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


TEST(NewtonsTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Method method(function);

    auto x = method.newtons(10);
    EXPECT_EQ(std::roundf(x), 2);
}

TEST(SecantTest, BasicAssertions) {
    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Method method(function);

    auto x = method.secant(10, 100);
    EXPECT_EQ(std::roundf(x), 2);
}

TEST(RegularFalsiTest, BasicAssertions) {
    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Method method(function);

    auto x = method.regular_falsi(0, 10);
    EXPECT_EQ(std::roundf(x), 2);

}