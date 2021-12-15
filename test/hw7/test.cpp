#include <gtest/gtest.h>

#include "global/genetic_algorithm.hpp"
#include "../function.hpp"

using namespace numerical_optimization;

auto functions = hw7::construct_functions();

TEST(GENETIC_ALGORITHM, BasicAssertions) {

    auto ga0 = GeneticAlgorithm<8, double, double>(functions[0]);
    // ga0.run(50);

    auto ga1 = GeneticAlgorithm<16, double, double>(functions[1]);
    // ga1.run(50);
    // ga1.print();
}