#include <gtest/gtest.h>

#define CSV_IO_NO_THREAD
#include "csv.h"

#include "lsm.h"
#include "lsm/gauss_newtons.hpp"
#include "lsm/LM.hpp"
#include "../function.hpp"

using namespace Eigen;
using namespace numerical_optimization;

auto functions = hw6::construct_functions();

static std::pair<std::vector<Vector3d>, std::vector<double>> read_data() {
    io::CSVReader<4> in("/home/hyeonjang/classes/numerical-optimization/test/hw6/observation.csv");
    in.read_header(io::ignore_extra_column, "x_i", "y_i", "z_i", "f_i");
    
    std::vector<Vector3d> vars;
    std::vector<double> vals;
    double x; double y; double z; double f;
    while(in.read_row(x, y, z, f)){
        
        // initialize observation data
        vars.push_back({x, y, z});
        vals.push_back(f);
    }
    return std::make_pair(vars, vals);
};

static auto data = read_data();

// read measurement data
TEST(GAUSSNEWTONS, BasicAssertions) {

    // std::cout << "function_1, initial: -1, -1, -1, 1" << std::endl;
    // auto gauss_newtons_0 = GaussNewtons(functions[0], {-1, -1, -1, 1});
    // gauss_newtons_0.set_observation(data.first, data.second);
    // std::cout << gauss_newtons_0.fit(200) << std::endl;

    // std::cout << "function_2, initial: -1, -1, -1, 1" << std::endl;
    // auto gauss_newtons_1 = GaussNewtons(functions[1], {-1, -1, -1, 1});
    // gauss_newtons_1.set_observation(data.first, data.second);
    // std::cout << gauss_newtons_1.fit(200) << std::endl;

    // std::cout << "function_1, initial: -2, -2, 2, 2" << std::endl;
    // auto gauss_newtons_2 = GaussNewtons(functions[0], {-2, -2, 2, 2});
    // gauss_newtons_2.set_observation(data.first, data.second);
    // std::cout << gauss_newtons_2.fit(200) << std::endl;

    // std::cout << "function_2, initial: -2, -2, 2, 2" << std::endl;
    // auto gauss_newtons_3 = GaussNewtons(functions[1], {-2, -2, 2, 2});
    // gauss_newtons_3.set_observation(data.first, data.second);
    // std::cout << gauss_newtons_3.fit(200) << std::endl;

    std::cout << "function_2, initial: 10, 10, 10, 10" << std::endl;
    auto gauss_newtons_4 = GaussNewtons(functions[0], {10, 10, 10, 10});
    gauss_newtons_4.set_observation(data.first, data.second);
    std::cout << gauss_newtons_4.fit(2) << std::endl;

    std::cout << "function_2, initial: 10, 10, 10, 10" << std::endl;
    auto gauss_newtons_5 = GaussNewtons(functions[1], {10, 10, 10, 10});
    gauss_newtons_5.set_observation(data.first, data.second);
    std::cout << gauss_newtons_5.fit(200) << std::endl;
}

TEST(LM, BasicAssertions) {

    std::cout << "function_1, initial: -1, -1, -1, 1" << std::endl;
    auto lm_0 = LM(functions[0], {-1, -1, -1, 1});
    lm_0.set_observation(data.first, data.second);
    std::cout << lm_0.fit(200) << std::endl;

    std::cout << "function_2, initial: -1, -1,- 1, 1" << std::endl;
    auto lm_1 = LM(functions[1], {-1, -1, -1, 1});
    lm_1.set_observation(data.first, data.second);
    std::cout << lm_1.fit(200) << std::endl;

    std::cout << "function_1, initial: -2, -2, 2, 2" << std::endl;
    auto lm_2 = LM(functions[0], {-2, -2, 2, 2});
    lm_2.set_observation(data.first, data.second);
    std::cout << lm_2.fit(200) << std::endl;

    std::cout << "function_2, initial: -2, -2, 2, 2" << std::endl;
    auto lm_3 = LM(functions[1], {-2, -2, 2, 2});
    lm_3.set_observation(data.first, data.second);
    std::cout << lm_3.fit(200) << std::endl;

    std::cout << "function_1, initial: 10, 10, 10, 10" << std::endl;
    auto lm_4 = LM(functions[0],  {10, 10, 10, 10});
    lm_4.set_observation(data.first, data.second);
    std::cout << lm_4.fit(200) << std::endl;

    std::cout << "function_2, initial: 10, 10, 10, 10" << std::endl;
    auto lm_5 = LM(functions[1],  {10, 10, 10, 10});
    lm_5.set_observation(data.first, data.second);
    std::cout << lm_5.fit(200) << std::endl;
}