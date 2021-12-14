#include <gtest/gtest.h>

#include "global/genetic_algorithm.hpp"
#include "../function.hpp"

using namespace numerical_optimization;

auto functions = hw7::construct_functions();

TEST(GENETIC_ALGORITHM, BasicAssertions) {

    auto ga = GeneticAlgorithm(functions[0]);
    ga.run();
}