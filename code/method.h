#include <limits>
#include <functional>
#include <cassert>

#include <boost/math/constants/constants.hpp>

namespace numerical_optimization
{

using function_t = std::function<float(const float&)>;

constexpr float MIN = 1e-4;// std::numeric_limits<float>::min();
constexpr float MAX = std::numeric_limits<float>::max();
constexpr float GOLDEN_RATIO = boost::math::constants::phi<float>();

class Method
{
public:
    Method(std::function<float(const float&)> f):function(f){};

    // assignment 1
    float bisection(float start, float end);
    float newtons(float x);
    float secant(float x1, float x0);
    float regular_falsi(float start, float end);
    float regular_falsi_not_recur(float start, float end);

    // assignment 2
    float fibonacci_search(float start, float end);
    float golden_section(float start, float end);

protected:
    std::function<float(const float&)> function;

    const int iter = 100000;

private:
    // for convenience
    bool  near_zero(float x) { return x==0 || -MIN<function(x)&&function(x)<MIN; }

    // for fibonacci
    std::vector<int> construct_fibonacci(size_t N);
};

std::pair<float, float> seeking_boundary(const function_t&);

};
