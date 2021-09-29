#include <functional>
#include <gtest/gtest.h>
#include <cmath>

#include "univariate.h"

using namespace numerical_optimization;
using namespace numerical_optimization::uni;

TEST(MethodsTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Univariate method(function);

    auto x = method.bisection(0, 128);
    EXPECT_EQ(std::roundf(x), 2);
}


TEST(NewtonsTest, BasicAssertions) {

    std::function<float(const float&)> function = [](float x) { return x*x - 4; };
    
    Univariate method(function);

    auto x = method.newtons(10);
    EXPECT_EQ(std::roundf(x), 2);
}

TEST(SecantTest, BasicAssertions) {
    std::function<float(const float&)> function = [](float x) { return x*x - 4; };

    Univariate method(function);

    auto x = method.secant(10, 100);
    EXPECT_EQ(std::roundf(x), 2);
}

TEST(RegularFalsiTest, BasicAssertions) {
 
    using function_t = std::function<float(const float&)>;

    function_t f0 = [](float x) { return x*x*x + 3*x*x + 9*x - 13; };

    Univariate m0(f0);

    auto x0 = m0.regular_falsi(0, 10);
    EXPECT_NEAR(x0, 1, 1e-4);

    function_t f1 = [](float x) { return std::sqrt(0.2/3.1415)*std::exp(-0.2*x*x*x)-3; };

    Univariate m1(f1);
    auto x1 = m1.regular_falsi(-3, -2);
}

// TEST(RegularFalsiNoRecurTest, BasicAssertions) {
 
//     using function_t = std::function<float(const float&)>;

//     function_t f0 = [](float x) { return x*x*x + 3*x*x + 9*x - 13; };

//     Method m0(f0);

//     auto x0 = m0.regular_falsi(0, 10);
//     auto y0 = m0.regular_falsi_not_recur(0, 10);
//     EXPECT_NEAR(x0, y0, 1e-4);

//     function_t f1 = [](float x) { return std::sqrt(0.2/3.1415)*std::exp(-0.2*x*x*x)-3; };

//     Method m1(f1);
//     auto x1 = m1.regular_falsi(-3, -2);
//     auto y1 = m1.regular_falsi_not_recur(-9, -2);
//     EXPECT_NEAR(x1, y1, 1e-4);

// }