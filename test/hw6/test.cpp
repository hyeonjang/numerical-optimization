#include <gtest/gtest.h>

#define CSV_IO_NO_THREAD
#include "csv.h"

#include "lsm.h"
#include "lsm/gauss_newtons.hpp"
#include "lsm/LM.hpp"
#include "../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;
using function_t = std::function<double(const Vector3d&, Vector4d&)>;

std::vector<function_t> functions = hw6::construct_functions();

// read measurement data
TEST(GAUSSNEWTONS, BasicAssertions) {

    io::CSVReader<4> in("/home/hyeonjang/classes/numerical-optimization/test/hw6/observation.csv");
    in.read_header(io::ignore_extra_column, "x_i", "y_i", "z_i", "f_i");
    
    std::vector<Vector3d> vars;
    std::vector<double> vals;
    std::vector<double> data;
    double x; double y; double z; double f;
    while(in.read_row(x, y, z, f)){
        
        // initialize observation data
        data.push_back(x);
        data.push_back(y);
        data.push_back(z);
        data.push_back(f);

        vars.push_back({x, y, z});
        vals.push_back(f);
    }

    Eigen::Map<Eigen::MatrixXd> observation(data.data(), 4, data.size()/4);

    auto gauss_newtons = GaussNewtons(functions[0]);
    gauss_newtons.set_observation(observation);
    Vector4d coefficient = gauss_newtons.fit();
    
    for(size_t i=0; i<vars.size(); i++) {
        auto return_val = functions[0](vars[i], coefficient);
        std::cout << return_val - vals[i] << std::endl;
    }

    std::cout << coefficient << std::endl;
}

TEST(LM, BasicAssertions) {

    io::CSVReader<4> in("/home/hyeonjang/classes/numerical-optimization/test/hw6/observation.csv");
    in.read_header(io::ignore_extra_column, "x_i", "y_i", "z_i", "f_i");
    
    std::vector<Vector3d> vars;
    std::vector<double> vals;
    std::vector<double> data;
    double x; double y; double z; double f;
    while(in.read_row(x, y, z, f)){
        
        // initialize observation data
        data.push_back(x);
        data.push_back(y);
        data.push_back(z);
        data.push_back(f);

        vars.push_back({x, y, z});
        vals.push_back(f);
    }

    Eigen::Map<Eigen::MatrixXd> observation(data.data(), 4, data.size()/4);

    auto lm = LM(functions[0]);
    lm.set_observation(observation);
    Vector4d coefficient = lm.fit();

    for(size_t i=0; i<vars.size(); i++) {
        auto return_val = functions[0](vars[i], coefficient);
        std::cout << return_val - vals[i] << std::endl;
    }

    std::cout << coefficient << std::endl;
}