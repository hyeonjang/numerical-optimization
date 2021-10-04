#include <gnuplot-iostream.h>

#include "multi/powells.hpp"
#include "../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

std::vector<function_t<Vector2f>> functions = construct_functions();
std::vector<Powells<Vector2f>> powells = construct_methods<Powells<Vector2f>>(functions);

int main() {

    powells[0].eval();
    auto plot0 = powells[0].plot;

    for(auto p:plot0) {
        std::cout << p << std::endl;
    }

}