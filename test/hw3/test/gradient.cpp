#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "multi/nelder_mead.hpp"
#include "multi/powells.hpp"
#include "../../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

using namespace Eigen;

auto functions = hw2::construct_functions();
TEST(Gradient, BasicAssertions) {
    // construct test function from example

    // 1. f(x, y) = (x+2*y-6)^2 + (2*x+y-6)^2
    // gradient(f(x, y)) = (10x+8y-36, 8x+10y-36)
    Multivariate<Vector2d> method0 = Multivariate<Vector2d>(functions[0]);

    // 2. f(x, y) = 50*(y-x*x)^2 + (1-x)^2
    // gradient(f(x, y)) = (-200x(y-x^2)-2(1-x), 100(y-x^2))
    Multivariate<Vector2d> method1 = Multivariate<Vector2d>(functions[1]);

    // 3. f(x, y) = (1.5-x+xy)^2 + (2.25-x+xy^2)^2 + (2.625 - x+ xy^3)^2
    // gradient(f(x, y)) = (2xy^6+2xy^4+5.25y^3=4xy^3+4.5y^2-2xy^2+3y-4xy+6x-12.75, 
    //                       6x^2y^5+4x^2y^3-6x^2y^2-2x^2y-2x^2+15.75xy^2+9xy+3x)
    Multivariate<Vector2d> method2 = Multivariate<Vector2d>(functions[2]);

    // test variables
    Vector2d v00 = { 0.f, 0.f }, v01 = { 0.0f, 1.0f }, v02 = { 0.0f, 2.0f };
    Vector2d v10 = { 1.f, 0.f }, v11 = { 1.0f, 1.0f }, v12 = { 1.0f, 2.0f };

    Vector2d eps = Vector2d(epsilon, epsilon);

    // [0] checking function
    auto near = [&] (Vector2d v0, Vector2d v1) {
        bool flag = abs((v0.norm()-v1.norm()))<10.f;  // too large
        return flag ? (::testing::AssertionSuccess()) : (::testing::AssertionFailure() << "\n" << v0);
    };

    // 1.
    EXPECT_TRUE(near(method0.gradient(v00), Vector2d(-36.f, -36.f)));
    EXPECT_TRUE(near(method0.gradient(v01), Vector2d(-24.f, -26.f)));
    EXPECT_TRUE(near(method0.gradient(v02), Vector2d(-20.f, -16.f)));
    EXPECT_TRUE(near(method0.gradient(v10), Vector2d(-26.f, -28.f)));
    EXPECT_TRUE(near(method0.gradient(v11), Vector2d(-18.f, -18.f)));
    EXPECT_TRUE(near(method0.gradient(v12), Vector2d(-10.f, -8.f)));

    // 2.
    EXPECT_TRUE(near(method1.gradient(v00), Vector2d(-2.f, 0.f)));
    EXPECT_TRUE(near(method1.gradient(v01), Vector2d(-2.f, 100.f)));
    EXPECT_TRUE(near(method1.gradient(v02), Vector2d(-2.f, 200.f)));
    EXPECT_TRUE(near(method1.gradient(v10), Vector2d(200.f, -100.f)));
    EXPECT_TRUE(near(method1.gradient(v11), Vector2d(0.f, 0.f)));
    EXPECT_TRUE(near(method1.gradient(v12), Vector2d(-200.f, 100.f)));

    // 3.
    EXPECT_TRUE(near(method2.gradient(v00), Vector2d(-12.75f, 0.f)));
    EXPECT_TRUE(near(method2.gradient(v01), Vector2d(0.f, 0.f)));
    EXPECT_TRUE(near(method2.gradient(v02), Vector2d(53.25f, 0.f)));
    EXPECT_TRUE(near(method2.gradient(v10), Vector2d(-6.75f, 1.f)));
    EXPECT_TRUE(near(method2.gradient(v11), Vector2d(0.f, 27.75f)));
    EXPECT_TRUE(near(method2.gradient(v12), Vector2d(171.25f, 278.f)));
}