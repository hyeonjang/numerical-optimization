#include "univariate.h"

namespace numerical_optimization {
namespace uni {

constexpr int number_functions = 5;

std::vector<function_t> construct_functions() {

    std::vector<function_t> functions;
    functions.resize(number_functions);
    functions[0] = [](float x){ return std::pow(x, 4)/4 + std::pow(x, 3) + std::pow(x, 2)*(9/2) - 10*x; };
    functions[1] = [](float x){ return std::sin(x) +x*x-10; };
    functions[2] = [](float x){ return -std::exp((-1*x*x)/(1.4*1.4)); };

    // not differentiable cases
    functions[3] = [](float x){ return std::abs(x-0.3); };
    functions[4] = [](float x){ return std::abs(std::log(x)); };

    return functions;
};

std::vector<Univariate> construct_methods(const std::vector<function_t> functions) {
    std::vector<Univariate> methods;
    for(auto func : functions) {
        methods.emplace_back(Univariate(func));
    }
    return methods;
}

std::vector<function_t> functions = construct_functions();
std::vector<Univariate> methods = construct_methods(functions);

}
}