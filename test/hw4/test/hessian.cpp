#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "multi/cauchys.hpp"
#include "../../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

using namespace Eigen;

std::vector<function_t<Vector2f>> functions = hw2::construct_functions<Vector2f>();

TEST(Hessian, BasicAssertions) {
    
    // 1. f(x, y) = (x+2*y-6)^2 + (2*x+y-6)^2
    // grad(f(x, y)) = (10x+8y-36, 8x+10y-36)
    // H(f(x, y))    = ( 10, 8
    //                    8, 10 )
    Multivariate<Vector2f> method0 = Multivariate<Vector2f>(functions[0]);

    // 2. f(x, y) = 50(y-x*x)^2 + (1-x)^2
    // grad(f(x, y)) = (-200x(y-x^2)-2(1-x), 100(y-x^2))
    // H(f(x, y))    = ( -200(y-x^2)+400x^2+2, -200x
    //                   -200x               ,  100 )
    Multivariate<Vector2f> method1 = Multivariate<Vector2f>(functions[1]);

    // 3. f(x, y) = (1.5-x+xy)^2 + (2.25-x+xy^2)^2 + (2.625-x+xy^3)^2
    // grad(f(x, y)) = (2xy^6+2xy^4+5.25y^3=4xy^3+4.5y^2-2xy^2+3y-4xy+6x-12.75, 6x^2y^5+4x^2y^3-6x^2y^2-2x^2y-2x^2+15.75xy^2+9xy+3x)
    // H(f(x, y))    = (2(y^3-1)^2+ 2(y^2-1)^2 + 2(y-1)^2 | 4 x (y^2 - 1) y + 4 y (x y^2 - x + 2.25) + 6 x (y^3 - 1) y^2 + 6 y^2 (x y^3 - x + 2.625) + 2 x (y - 1) + 2 (x y - x + 1.5)
    //                  4 x (y^2 - 1) y + 4 y (x y^2 - x + 2.25) + 6 x (y^3 - 1) y^2 + 6 y^2 (x y^3 - x + 2.625) + 2 x (y - 1) + 2 (x y - x + 1.5) | 18 x^2 y^4 + 8 x^2 y^2 + 2 x^2 + 12 x y (x y^3 - x + 2.625) + 4 x (x y^2 - x + 2.25))
    Multivariate<Vector2f> method2 = Multivariate<Vector2f>(functions[2]);

    // test variables
    Vector2f v00 = { 0.f, 0.f }, v01 = { 0.0f, 1.0f }, v02 = { 0.0f, 2.0f };
    Vector2f v10 = { 1.f, 0.f }, v11 = { 1.0f, 1.0f }, v12 = { 1.0f, 2.0f };

    // [0] checking function
    auto near = [](Matrix2f m0, Matrix2f m1){
        bool flag = abs(m0.norm()-m1.norm())<0.2f;
        return flag ? (::testing::AssertionSuccess()) : (::testing::AssertionFailure() << "\n" << m0 << "\n" <<  m1);
    };

    // 1.
    Matrix2f comp1;
    comp1 << 10.f, 8.f, 8.f, 10.f;
    EXPECT_TRUE(near(method0.hessian(v00), comp1));
    EXPECT_TRUE(near(method0.hessian(v01), comp1));
    EXPECT_TRUE(near(method0.hessian(v02), comp1));
    EXPECT_TRUE(near(method0.hessian(v10), comp1));
    EXPECT_TRUE(near(method0.hessian(v11), comp1));
    EXPECT_TRUE(near(method0.hessian(v12), comp1));

    // 2.
    auto comp2 = [](Vector2f v) {
        Matrix2f m;
        return (m << -200*(v[1]-v[0]*v[0])+400*v[0]*v[0]+2, -200*v[0], 
                    -200*v[0], 100).finished();
    };
    EXPECT_TRUE(near(method1.hessian(v00), comp2(v00)));
    EXPECT_TRUE(near(method1.hessian(v01), comp2(v01)));
    EXPECT_TRUE(near(method1.hessian(v02), comp2(v02)));
    EXPECT_TRUE(near(method1.hessian(v10), comp2(v10)));
    EXPECT_TRUE(near(method1.hessian(v11), comp2(v11)));
    EXPECT_TRUE(near(method1.hessian(v12), comp2(v12)));


}