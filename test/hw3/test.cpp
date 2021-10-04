#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "multi/nelder_mead.hpp"
#include "multi/powells.hpp"
#include "../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

using namespace Eigen;

std::vector<function_t<Vector2f>> functions = construct_functions();
std::vector<NelderMead<Vector2f>> nelder_meads = construct_methods<NelderMead<Vector2f>>(functions);

TEST(NeldaMead, BasicAssertions) {
    
    nelder_meads[0].eval();

    // methods[0].nelder_mead();

    // methods[1].nelder_mead();

    // methods[2].nelder_mead();

}

TEST(powells, BasicAssertions) {
    
}