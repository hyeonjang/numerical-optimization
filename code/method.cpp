#include <random>
#include <algorithm>
#include "method.h"

using namespace numerical_optimization;

// random function for boundary seeking
int Method::random_int() const {
    // threshold
    constexpr int scale = 100000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(
        std::numeric_limits<int>::min()/scale, 
        std::numeric_limits<int>::max()/scale
        );
    return distrib(gen);
}