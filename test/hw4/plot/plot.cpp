#include <string>
#include <gnuplot-iostream.h>
#include <nlohmann/json.hpp>

#include "multi/cauchys.hpp"
#include "multi/newtons.hpp"
#include "multi/quasi_newtons.hpp"
#include "../../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;
using namespace numerical_optimization::multi;
using json = nlohmann::json;

auto functions = hw2::construct_functions();
std::vector<std::string> functions_str = {
	"f(x, y)=(x+2*y-6)**2 + (2*x+y-6)**2\n",
	"f(x, y)=50*(y-x*x)**2 + (1-x)**2\n",
	"f(x, y)=(1.5-x+x*y)**2 + (2.25-x+x*(y**2))**2 + (2.625 - x+ x*(y**3))**2\n",
};

template<typename Method>
static void call_function(std::string title, int idx, std::array<float, 2> init, std::array<float, 2> x, std::array<float, 2> y) {
	
	std::stringstream str1; str1 << std::fixed << std::setprecision(2) << init[0];
	std::stringstream str2; str2 << std::fixed << std::setprecision(2) << init[1];
	
	//  set output

	// call method
	Method method(functions[idx]);
	method.eval(Vector2d(init[0], init[1]), 1e-6);

    std::vector<std::tuple<float, float, float>> tuples;
    for(auto const& o:method.plot) {
        tuples.emplace_back(std::make_tuple(o.first[0], o.first[1], o.second));
	}

	// initialize
	Gnuplot gp;
	
	std::string output = " './" + title + ".png" + "' ";
	std::string name   = "'" + title + "'";

	gp << "set term png\n";
	gp << "set output " + output + "\n";
	gp << "set title "  + name + "\n";
	gp << "unset xlabel\n";
	gp << "unset ylabel\n";
	gp << "set view map\n";

	// set range
	std::string xrange = "set xrange [" + std::to_string(x[0]) + ":" + std::to_string(x[1]) + "]\n"; 
	std::string yrange = "set yrange [" + std::to_string(y[0]) + ":" + std::to_string(y[1]) + "]\n"; 

	gp << functions_str[idx];
	gp << xrange.c_str();
	gp << yrange.c_str();
	gp << "set zrange [-100:100]\n";
	
	gp << "set multiplot\n";

	// draw dataset
	gp << "splot '-' with lines title 'pattern'\n";
	gp.send1d(tuples);

	// draw function
	gp << "unset clabel\n";
	gp << "set contour base\n";
	gp << "set cntrparam levels incremental -30, 3, 30 \n";
	gp << "set isosample 500, 100\n";
	gp << "unset surface\n";
	gp << "splot f(x, y) lt rgb '#000000'\n";
	// end

	gp << "unset multiplot\n";

	std::cout <<"["<<title<<"]"
		<< " optimal point: " 
		<< std::get<0>(tuples.back()) << ", " << std::get<1>(tuples.back()) 
		<< ", " << std::get<0>(tuples[tuples.size()-2]) << ", " << std::get<1>(tuples[tuples.size()-2]) 
		<< std::endl;
}

// please to refer http://gnuplot.sourceforge.net/demo_5.2/random.html
int main(int argc, char *argv[]) {

	std::string file_path = __FILE__;
	std::string curr_path = file_path.substr(0, file_path.rfind("/"));
	std::ifstream conf_file;
	
	if(argc==1) {
		conf_file = std::ifstream(curr_path + "/config.json");
	} else {
		conf_file = std::ifstream(argv[1]);
	}

	try {
		json conf_json;
		if(!conf_file.is_open()) throw;

		conf_file >> conf_json;
		for(auto& elem : conf_json) {

			std::string method;
			elem["method"].get_to(method);

			if(method.compare("Cauchys")==0) {
				call_function<Cauchys<Vector2d>>(elem["title"], elem["index"], elem["initial"], elem["xrange"], elem["yrange"]);
			} else if(method.compare("Newtons")==0) {
				call_function<Newtons<Vector2d>>(elem["title"], elem["index"], elem["initial"], elem["xrange"], elem["yrange"]);
			} else if(method.compare("QuasiNewtons")==0) {
				
				std::string update_method;
				elem["detail"].get_to(update_method);

			 	if(update_method.compare("SR1")==0) {
			 		call_function<QuasiNewtons<Vector2d, quasi_newtons::Rank::SR1>>(elem["title"], elem["index"], elem["initial"], elem["xrange"], elem["yrange"]);
			 	} else if(update_method.compare("BFGS")==0) {
			 		call_function<QuasiNewtons<Vector2d, quasi_newtons::Rank::BFGS>>(elem["title"], elem["index"], elem["initial"], elem["xrange"], elem["yrange"]);
			 	}
			}
		}
	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
}