#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "multivariate.h"
#include "../../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

using namespace Eigen;

std::vector<function_t<Vector2f>> functions = construct_functions();

// f(x, y) = (x+2*y-6)^2 + (2*x+y-6)^2
// gradient(f(x, y)) = (10x+8y-36, 8x+10y-36)
Multivariate<Vector2f> method0 = Multivariate<Vector2f>(functions[0]);

// f(x, y) = 50*(y-x*x)^2 + (1-x)^2
// gradient(f(x, y)) = (-200x(y-x^2)-2(1-x), 100(y-x^2))
Multivariate<Vector2f> method1 = Multivariate<Vector2f>(functions[1]);

// f(x, y) = (1.5-x+xy)^2 + (2.25-x+xy^2)^2 + (2.625 - x+ xy^3)^2
// gradient(f(x, y)) = (2xy^6+2xy^4+5.25y^3=4xy^3+4.5y^2-2xy^2+3y-4xy+6x-12.75, 6x^2y^5+4x^2y^3-6x^2y^2-2x^2y-2x^2+15.75xy^2+9xy+3x)
Multivariate<Vector2f> method2 = Multivariate<Vector2f>(functions[2]);

// 1. Difference of two consecutive estimates
TEST(ConsecutiveDifference, BasicAssertions) {
    // test variables
    Vector2f v00 = { 0.f, 0.f }, v01 = { 0.0f, 1.0f }, v02 = { 0.0f, 2.0f };
    Vector2f v10 = { 1.f, 0.f }, v11 = { 1.0f, 1.0f }, v12 = { 1.0f, 2.0f };

    // 1) 
    EXPECT_TRUE(method0.terminate<Termination::Condition::ConsecutiveDifference>({v00}));

    // 2)
    EXPECT_FALSE(method0.terminate<Termination::Condition::ConsecutiveDifference>({v00, v01, v02}));
    EXPECT_TRUE(method0.terminate<Termination::Condition::ConsecutiveDifference>({v00, v01, v02}, 3));
}

TEST(RelativeConsecutiveDifference, BasicAssertions) {

}

TEST(MagnituteGradient, BasicAssertions) {
    // test variables
    Vector2f v00 = { 0.f, 0.f }, v01 = { 0.0f, 0.0f }, v02 = { 0.0f, 2.0f };
    Vector2f v10 = { 1.f, 0.f }, v11 = { 1.0f, 1.0f }, v12 = { 1.0f, 2.0f };

    // 1)
    // EXPECT_TRUE(method0.magnitude_gradient({v00, v01, v02}));
}

TEST(RelativeDifferenceFunctionValues, BasicAssertions) {

}

TEST(DescentDirectionChange, BasicAssertions) {

}