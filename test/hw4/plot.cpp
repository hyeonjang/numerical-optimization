#include <string>
#include <gnuplot-iostream.h>

#include "multi/cauchys.hpp"
#include "multi/newtons.hpp"
#include "multi/quasi_newtons.hpp"
#include "../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;
using namespace numerical_optimization::multi;

std::vector<function_t<Vector2f>> functions = construct_functions();
std::vector<std::string> functions_str = {
	"f(x, y)=(x+2*y-6)**2 + (2*x+y-6)**2\n",
	"f(x, y)=50*(y-x*x)**2 + (1-x)**2\n",
	"f(x, y)=(1.5-x+x*y)**2 + (2.25-x+x*(y**2))**2 + (2.625 - x+ x*(y**3))**2\n",
};

template<typename Method>
static void call_function(std::string title, int idx, std::pair<float, float> x, std::pair<float, float> y) {
	
	//  set output
	std::string output = " './" + title + "_" + std::to_string(idx) + ".png" + "' ";
	std::string name   = " '" +title+ "' ";

	// call method
	Method method(functions[idx]);
	method.eval(Vector2f(2.2f, 3.2f), 1e-4);

    std::vector<std::tuple<float, float, float>> tuples;
    for(auto const& o:method.plot) {
        tuples.emplace_back(std::make_tuple(o.first[0], o.first[1], o.second));
		// std::cout << "optimal point: " << std::get<0>(tuples.back()) << ", " << std::get<1>(tuples.back()) << std::endl;
	}

	// initialize
	Gnuplot gp;

	gp << "set term png\n";
	gp << "set output " + output + "\n";
	gp << "set title "  + name   + "\n";
	gp << "unset xlabel\n";
	gp << "unset ylabel\n";
	gp << "set view map\n";

	// set range
	std::string xrange = "set xrange [" + std::to_string(x.first) + ":" + std::to_string(x.second) + "]\n"; 
	std::string yrange = "set yrange [" + std::to_string(y.first) + ":" + std::to_string(y.second) + "]\n"; 

	gp << functions_str[idx];
	gp << xrange.c_str();
	gp << yrange.c_str();
	gp << "set zrange [-100:100]\n";
	
	gp << "set multiplot\n";

	// draw dataset
	gp << "splot '-' with lines title 'pattern'\n";
	gp.send1d(tuples);

	// draw function
	gp << "set contour base\n";
	gp << "set cntrparam levels incremental -30, 3, 30 \n";
	gp << "set isosample 500, 100\n";
	gp << "unset surface\n";
	gp << "splot f(x, y) with lines\n";
	// end
	gp << "unset multiplot\n";

	std::cout <<"["<<title<<"]"
		<< " optimal point: " 
		<< std::get<0>(tuples.back()) << ", " << std::get<1>(tuples.back()) 
		<< std::endl;
}

// please to refer http://gnuplot.sourceforge.net/demo_5.2/random.html
int main(int argc, char *argv[]) {

	std::pair<float, float> xrange, yrange;
	if(argc==1) {
		// call_function<Cauchys<Vector2f>>("Cauchys", 0, std::make_pair(-5, 7), std::make_pair(-1, 6) );
		// call_function<Cauchys<Vector2f>>("Cauchys", 1, std::make_pair(-0.5, 1.5), std::make_pair(-0.5, 2.0) );
		// call_function<Cauchys<Vector2f>>("Cauchys", 2, std::make_pair(1, 5), std::make_pair(-0.5, 1.5) );
		// call_function<Newtons<Vector2f>>("Newtons", 0, std::make_pair(-5, 7), std::make_pair(-5, 7) );
		// call_function<Newtons<Vector2f>>("Newtons", 1, std::make_pair(-0.5, 1.5), std::make_pair(-0.5, 2.0) );
		// call_function<Newtons<Vector2f>>("Newtons", 2, std::make_pair(1, 5), std::make_pair(-0.5, 1.5) );
		// call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::SR1>>("SR1", 0, std::make_pair(-5, 7), std::make_pair(-5, 7) );
		// call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::SR1>>("SR1", 1, std::make_pair(-0.5, 1.5), std::make_pair(-0.5, 2.0) );
		// call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::SR1>>("SR1", 2, std::make_pair(1, 5), std::make_pair(-0.5, 1.5) );
		// call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::BFGS>>("BFGS", 0, std::make_pair(-5, 7), std::make_pair(-5, 7) );
		// call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::BFGS>>("BFGS", 1, std::make_pair(-0.5, 1.5), std::make_pair(-0.5, 2.0) );
		call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::BFGS>>("BFGS", 2, std::make_pair(1, 5), std::make_pair(-0.5, 1.5) );
		return 0;
	} else {
		xrange = std::make_pair(std::stof(argv[3]), std::stof(argv[4]));
		yrange = std::make_pair(std::stof(argv[5]), std::stof(argv[6]));

		if ( strcmp(argv[1], "Cauchys")==0 ) {
			call_function<Cauchys<Vector2f>>("Cauchys", std::stoi(argv[2]), xrange, yrange );
		} else if (strcmp(argv[1], "Newtons")==0) {
			call_function<Newtons<Vector2f>>("Newtons", std::stoi(argv[2]), xrange, yrange );
		} else if (strcmp(argv[1], "SR1")==0) {
			call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::SR1>>("SR1", std::stoi(argv[2]), xrange, yrange );
		} else if (strcmp(argv[1], "Newtons")==0) {
			call_function<QuasiNewtons<Vector2f, quasi_newtons::Rank::BFGS>>("Newtons", std::stoi(argv[2]), xrange, yrange );
		} 
		return 0;
	}


}