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

TEST(ConsecutiveDifference, BasicAssertions) {

}