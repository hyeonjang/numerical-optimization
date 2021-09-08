#include <limits>
#include <functional>
#include <cassert>

namespace numerical_optimization
{

constexpr float MIN = std::numeric_limits<float>::min();

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
