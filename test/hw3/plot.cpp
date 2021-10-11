#include <string>
#include <gnuplot-iostream.h>

#include "multi/powells.hpp"
#include "multi/nelder_mead.hpp"
#include "../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;
using namespace numerical_optimization::multi;

std::vector<function_t<Vector2f>> functions = construct_functions();
std::vector<std::string> functions_str = {
	"f(x, y)=(x+2*y)**2 + (2*x+y)**2\n",
	"f(x, y)=50*(y-x*x)**2 + (1-x)**2\n",
	"f(x, y)=(1.5-x+x*y)**2 + (2.25-x+x*(y**2))**2 + (2.625 - x+ x*(y**3))**2\n",
};

template<typename Method>
static void call_function(std::string title, int idx) {
	
	//  set output
	std::string output = " '" + title + "_" + std::to_string(idx) + ".png" + "' ";
	title = " '" +title+ "' ";

	// call method
	Method method(functions[idx]);
	method.eval();

    std::vector<std::tuple<float, float, float>> pairs;
    for(Vector2f const& o:method.plot) {
        pairs.emplace_back(std::make_tuple(o[0], o[1], 0));
    }

	// initialize
	Gnuplot gp;

	gp << "set output " + output + "\n";
	gp << "set title "  + title + "\n";
	gp << "unset xlabel\n";
	gp << "unset ylabel\n";
	gp << "unset zlabel\n";
	gp << "set view map\n";

	// set range
	gp << functions_str[idx];
	gp << "set xrange [-3:3]\n";
	gp << "set yrange [-2:2]\n";
	// gp << "set zrange [-10:10]\n";

	gp << "set multiplot\n";

	// draw dataset
	gp << "splot '-' with lines title 'pattern'\n";
	gp.send1d(pairs);

	// draw function
	gp << "set contour base\n";
	gp << "set cntrparam levels incremental -3, 0.5, 3 \n";
	gp << "set isosample 500, 100\n";
	gp << "unset surface\n";
	gp << "splot f(x, y) with lines\n";
	// end
	gp << "unset multiplot\n";
}

// please to refer http://gnuplot.sourceforge.net/demo_5.2/random.html
int main(int argc, char *argv[]) {

	if ( strcmp(argv[1], "powells")==0 ) {
		call_function<Powells<Vector2f>>("Powells", std::stoi(argv[2]));
	} else if (strcmp(argv[1], "neldermead")==0) {
		call_function<NelderMead<Vector2f>>("NelderMead", std::stoi(argv[2]));
	} else {
		std::cout << "Enter the arguments: argv[1]:Method, argv[2]:function index" << std::endl;
	}
}