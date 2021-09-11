#include <limits>
#include <functional>
#include <cassert>

namespace numerical_optimization
{

constexpr float MIN = 1e-6;// std::numeric_limits<float>::min();
constexpr float MAX = std::numeric_limits<float>::max();

class Method
{
public:
    Method(std::function<float(const float&)> f):function(f){};

    float bisection(float start, float end);
    float newtons(float x);
    float secant(float x1, float x0);
    float regular_falsi(float start, float end);

protected:
    std::function<float(const float&)> function;
};
};
