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

namespace multi {
template<typename Multi, typename Vector2T>
std::vector<Multi> construct_methods(std::vector<function_t<Vector2T>> functions) {
    std::vector<Multi> methods;
    for(const auto& func:functions) {
        methods.emplace_back(Multi(func));
    }
    return methods;
};
///////////////////////////////////
} // the end of namespace multi ///
//////////////////////////////////

namespace hw2 {
constexpr int number_functions = 3;

using namespace multi;
using namespace Eigen;

template<typename Vector2T>
std::vector<function_t<Vector2T>> construct_functions() {
    std::vector<function_t<Vector2T>> functions(number_functions);

    functions[0] = [](Vector2T var) {
        return std::pow((var[0]+2*var[1]-6), 2)
        + std::pow((2*var[0]+var[1]-6), 2);
    };

    functions[1] = [](Vector2T var) {
        return 50*std::pow((var[1]-var[0]*var[0]), 2) 
        + std::pow((1.0-var[0]), 2);
    };

    functions[2] = [](Vector2T var) {
        return std::pow((1.5-var[0]+var[0]*var[1]), 2) 
        + std::pow((2.25-var[0]+var[0]*var[1]*var[1]), 2)
        + std::pow((2.625-var[0]+var[0]*var[1]*var[1]*var[1]), 2);
    };

    return functions;
};
//////////////////////////////////
} // the end of namespace hw2  ///
//////////////////////////////////

namespace hw5 {
constexpr int number_functions = 3;

using namespace multi;
using namespace Eigen;

template<typename Vector2T>
std::vector<function_t<Vector2T>> construct_functions() {
    std::vector<function_t<Vector2T>> functions(number_functions);

    functions[0] = [](Vector2T var) {
        return std::pow((var[0]+2*var[1]-7), 2)
        + std::pow((2*var[0]+var[1]-5), 2);
    };

    functions[1] = [](Vector2T var) {
        return 40*std::pow((var[1]-var[0]*var[0]), 2) 
        + std::pow((1.0-var[0]), 2);
    };

    functions[2] = [](Vector2T var) {
        return std::pow((1.5-var[0]+var[0]*var[1]), 2) 
        + std::pow((2.25-var[0]+var[0]*var[1]*var[1]), 2)
        + std::pow((2.625-var[0]+var[0]*var[1]*var[1]*var[1]), 2);
    };

    return functions;
};
//////////////////////////////////
}/// the end of namespace hw5 ////
//////////////////////////////////

namespace hw6 {
constexpr int number_functions = 2;

using namespace Eigen;
using function_t = std::function<double(const Vector4d&, const Vector3d&)>;

std::vector<function_t> construct_functions() {
    std::vector<function_t> functions(number_functions);

    functions[0] = [](Vector4d coeff, Vector3d vars) {
        return coeff[0]*vars[0] + coeff[1]*vars[1] + coeff[2]*vars[2] + coeff[3];
    };

    functions[1] = [](Vector4d coeff, Vector3d vars) {
        double value = -(pow(vars[0]-coeff[0], 2) + pow(vars[1]-coeff[1], 2) + pow(vars[2]-coeff[2], 2))/pow(coeff[3], 2);
        return std::exp(value);
    };

    return functions;
}
//////////////////////////////////
}/// the end of namespace hw6 ////
//////////////////////////////////

namespace hw7 {
// 1. $f(x) = 2(x-0.5)^2 + 1$
// 2. $f(x) = |x-0.5|(cos(12\phi[x-0.5])+1) + 1$
constexpr int number_functions = 2;

using namespace Eigen;
using function_t = std::function<double(const double&)>;

std::vector<function_t> construct_functions() {
    std::vector<function_t> functions(number_functions);

    functions[0] = [](double var) {
        return std::pow((var-0.5),2) + 1;
    };

    functions[1] = [](double var) {
        auto x = var - 0.5;
        return std::abs(x)*(std::cos(12*3.14*x)+1.2);
    };

    return functions;
}
//////////////////////////////////
}/// the end of namespace hw7 ////
//////////////////////////////////

/////////////////////////////////////////////////////
} // the end of namespace numerical_optimization ////
/////////////////////////////////////////////////////