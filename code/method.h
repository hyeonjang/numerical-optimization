#include <limits>
#include <functional>
#include <cassert>
#include <vector>

namespace numerical_optimization
{

using function_t = std::function<float(const float&)>;
using boundary_t = std::pair<float, float>;

constexpr float  MIN = 1e-4;// std::numeric_limits<float>::min();
constexpr float  MAX = std::numeric_limits<float>::max();
constexpr float  GOLDEN_RATIO = 1.f/1.618033988749895f;
constexpr size_t FIBONACCI_MAX = 46;

class Method {
public:
    Method(function_t f):function(f) { boundary = seeking_bound(5); };
    Method(function_t f, boundary_t b):function(f), boundary(b){};

    // assignment 1
    float bisection(float start, float end);
    float newtons(float x);
    float secant(float x1, float x0);
    float regular_falsi(float start, float end);

    // assignment 2
    float fibonacci_search(size_t N=FIBONACCI_MAX);
    float fibonacci_search(float start, float end, size_t N);
    float golden_section(size_t N=FIBONACCI_MAX);
    float golden_section(float start, float end, size_t N);

    // for convienience
    boundary_t get_bound() const;
    Method     derivate() const;
private:
    function_t function;
    boundary_t boundary;
    const size_t iter = 10000000; // termination condition

    std::vector<int> construct_fibonacci(size_t N) const; // for fibonacci search
    boundary_t seeking_bound(float step_size);
    int random_int() const;
    bool near_zero(float x) { 
        return x==0 || (-MIN<function(x)&&function(x)<MIN); 
    }
};

};
