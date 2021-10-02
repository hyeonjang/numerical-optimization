#include <vector>
#include <Eigen/Dense>
#include "univariate.h"
#include "multivariate.h"

namespace numerical_optimization {

// homework1 univariate optimization problems
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
};


/////////////////////////////////
} // the end of namespace uni //
////////////////////////////////

// homework2 multivariate optimization problems
namespace multi {
constexpr int number_functions = 3;

using namespace Eigen;

std::vector<function_t<Vector2f>> construct_functions() {
    std::vector<function_t<Vector2f>> functions(number_functions);

    functions[0] = [](Vector2f var) {
        return std::pow((var[0] + 2*var[1]), 2)
        + std::pow((2*var[0]+var[1]), 2);
    };

    functions[1] = [](Vector2f var) {
        return 50*std::pow((var[1]-var[0]*var[0]), 2) 
        + std::pow((1.0-var[0]), 2);
    };

    functions[2] = [](Vector2f var) {
        return std::pow((1.5-var[0]+var[0]*var[1]), 2) 
        + std::pow((2.25-var[0]+var[0]*var[1]*var[1]), 2)
        + std::pow((2.625-var[0]+var[0]*var[1]*var[1]*var[1]), 2);
    };

    return functions;
};

template<typename Multi>
std::vector<Multi> construct_methods(std::vector<function_t<Vector2f>> functions) {
    std::vector<Multi> methods;
    for(const auto& func:functions) {
        methods.emplace_back(Multi(func));
    }
    return methods;
};
//////////////////////////////////
} // the end of namespace multi //
//////////////////////////////////

///////////////////////////////////////////////////
} // the end of namespace numerical_optimization //
///////////////////////////////////////////////////