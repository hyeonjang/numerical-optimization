#include <functional>
#include <gtest/gtest.h>

#include "method.h"

using namespace numerical_optimization;

TEST(MethodsTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x/4 - 8; };

    Method method(function);

    auto x = method.bisection(0, 128);

    std::cout << "bisection: " << x << std::endl;

    EXPECT_EQ(x, 32);
}


TEST(NewtonsTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Method method(function);

    auto x = method.newtons(10);
}

TEST(SecantTest, BasicAssertions) {
    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Method method(function);

    auto x = method.secant(10, 100);
}

TEST(RegularFalsiTest, BasicAssertions) {
    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Method method(function);

    auto x = method.regular_falsi(1, 10);

}
