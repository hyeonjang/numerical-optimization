#include <gnuplot-iostream.h>

#include "multi/powells.hpp"
#include "../function.hpp"

using namespace numerical_optimization;
using namespace numerical_optimization::multi;

std::vector<function_t<Vector2f>> functions = construct_functions();
std::vector<Powells<Vector2f>> powells = construct_methods<Powells<Vector2f>>(functions);

std::vector<std::tuple<float, float>> convert_to_pair(std::vector<Vector2f> orig) {

    std::vector<std::tuple<float, float>> pairs;
    for(Vector2f const& o:orig) {
        pairs.emplace_back(std::make_tuple(o[0], o[1]));
    }
    return pairs;
};

void call_function_1() {
	Gnuplot gp;

	// set range
	gp << "set xrange [-2:2]\n";
	gp << "set yrange [-2:2]\n";

	// objective function
	gp << "splot (x+2*y)**2 + (2*x+y)**2\n";
}

// please to refer http://gnuplot.sourceforge.net/demo_5.2/random.html
int main() {
	call_function_1();
}