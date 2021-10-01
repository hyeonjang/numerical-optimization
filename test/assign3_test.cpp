#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "multivariate.h"
#include "function.cpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

using namespace Eigen;

// std::vector<function_t<Vector2f>> functions = construct_functions();
// std::vector<Multivariate<Vector2f>> methods = construct_methods(functions);

// extern std::vector<Multivariate<Vector2f>> numerical_optimization::multi::methods;

TEST(NeldaMead, BasicAssertions) {
    
    std::vector<function_t<Vector2f>> functions = construct_functions();
    std::vector<Multivariate<Vector2f>> methods = construct_methods(functions);

    methods[0].nelder_mead();

    // methods[0].nelder_mead();

    // methods[1].nelder_mead();

    // methods[2].nelder_mead();

}

TEST(powells, BasicAssertions) {
    
}