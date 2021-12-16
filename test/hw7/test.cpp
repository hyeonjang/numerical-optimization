#include <fstream>
#include <gtest/gtest.h>

#include "global/genetic_algorithm.hpp"
#include "../function.hpp"

using namespace numerical_optimization;

auto functions = hw7::construct_functions();

TEST(GENETIC_ALGORITHM_FUNC_GENOTYPELENGTH, BasicAssertions) {

    std::string filename("genotype.txt");
    {
        std::string filepath = "func0_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 8 << " " << 16 << " " << 32 << " " << 64 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 8>(functions[0]);
            double val0 = ga0.run(5);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[0]);
            double val1 = ga1.run(5);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 32>(functions[0]);
            double val2 = ga2.run(5);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 64>(functions[0]);
            double val3 = ga3.run(5);
            file << val3 << "\n";
        }
        file.close();
    }
    {
        std::string filepath = "func1_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 8 << " " << 16 << " " << 32 << " " << 64 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 8>(functions[1]);
            double val0 = ga0.run(30);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[1]);
            double val1 = ga1.run(30);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 32>(functions[1]);
            double val2 = ga2.run(30);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 64>(functions[1]);
            double val3 = ga3.run(30);
            file << val3 << "\n";
        }
        file.close();
    }
}

TEST(GENETIC_ALGORITHM_FUNC_POPULATION, BasicAssertions) {

    std::string filename("population.txt");
    {
        std::string filepath = "func0_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 20 << " " << 30 << " " << 50 << " " << 100 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 16>(functions[0], 20);
            double val0 = ga0.run(5);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[0], 30);
            double val1 = ga1.run(5);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 16>(functions[0], 50);
            double val2 = ga2.run(5);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 16>(functions[0], 100);
            double val3 = ga3.run(5);
            file << val3 << "\n";
        }
        file.close();
    }
    {
        std::string filepath = "func1_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 20 << " " << 30 << " " << 50 << " " << 100 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 16>(functions[1], 20);
            double val0 = ga0.run(30);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[1], 30);
            double val1 = ga1.run(30);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 16>(functions[1], 50);
            double val2 = ga2.run(30);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 16>(functions[1], 100);
            double val3 = ga3.run(30);
            file << val3 << "\n";
        }
        file.close();
    }
}

TEST(GENETIC_ALGORITHM_FUNC_CROSSOVER, BasicAssertions) {

    std::string filename("crossover.txt");
    {
        std::string filepath = "func0_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 0.1 << " " << 0.4 << " " << 0.7 << " " << 1.0 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.1);
            double val0 = ga0.run(5);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.4);
            double val1 = ga1.run(5);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.7);
            double val2 = ga2.run(5);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 1.0);
            double val3 = ga3.run(5);
            file << val3 << "\n";
        }
        file.close();
    }
    {
        std::string filepath = "func1_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 0.1 << " " << 0.4 << " " << 0.7 << " " << 1.0 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.1);
            double val0 = ga0.run(30);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.4);
            double val1 = ga1.run(30);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.7);
            double val2 = ga2.run(30);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 1.0);
            double val3 = ga3.run(30);
            file << val3 << "\n";
        }
        file.close();
    }
}

TEST(GENETIC_ALGORITHM_FUNC_MUTATION, BasicAssertions) {

    std::string filename("mutation.txt");
    {
        std::string filepath = "func0_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 0.01 << " " << 0.1 << " " << 0.4 << " " << 0.7 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.7, 0.01);
            double val0 = ga0.run(5);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.7, 0.1);
            double val1 = ga1.run(5);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.7, 0.4);
            double val2 = ga2.run(5);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 16>(functions[0], 50, 0.7, 0.7);
            double val3 = ga3.run(5);
            file << val3 << "\n";
        }
        file.close();
    }
    {
        std::string filepath = "func1_comp_" + filename;
        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << filename << "\n";
            file << 0.01 << " " << 0.1 << " " << 0.4 << " " << 0.7 << "\n"; 

            auto ga0 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.7, 0.01);
            double val0 = ga0.run(30);
            file << val0 << " ";

            auto ga1 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.7, 0.1);
            double val1 = ga1.run(30);
            file << val1 << " ";
        
            auto ga2 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.7, 0.4);
            double val2 = ga2.run(30);
            file << val2 << " ";

            auto ga3 = GeneticAlgorithm<double, double, 16>(functions[1], 50, 0.7, 0.7);
            double val3 = ga3.run(30);
            file << val3 << "\n";
        }
        file.close();
    }
}