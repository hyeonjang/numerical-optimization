#ifndef __GENETIC_ALGORITHM__
#define __GENETIC_ALGORITHM__

#include "method.h"
#include "ga_helper.h"

#include <random>
#include <bitset>
#include <fstream>

namespace numerical_optimization {

template<typename vector_t, typename scalar_t, size_t len>
class GeneticAlgorithm : Method {
public:
    using boundary_t = std::pair<scalar_t, scalar_t>;
    using function_t = std::function<scalar_t(const vector_t&)>;
    using populate_t = population_t<vector_t, scalar_t, len>;

    GeneticAlgorithm(
        function_t func, 
        size_t size=50, 
        double crossover_prob=0.7,
        double mutation_prob=0.1
        )
    :function(func)
    ,population(size)
    ,prob_crossover(crossover_prob)
    ,prob_mutation(mutation_prob)
    ,sum_fitness(0)
    ,sum_probability(0) {
        initialize_population();
        evaluate_fitness();
    }

    // iteration run
    double run(size_t iteration) {

        for(size_t i=0; i<iteration; i++) {
            // selection
            size_t idx_x, idx_y;
            std::tie(idx_x, idx_y) = select_with_rouletewheeling();

            /////////////////////////////////////////////////
            // reproduction
            /////////////////////////////////////////////////
            populate_t x = population[idx_x];
            populate_t y = population[idx_y];
            // cross over
            if(crossover(x, y)) {
                x.eval_fitness(function); 
                y.eval_fitness(function);
            }

            // mutation
            if(mutation(x) || mutation(y)) {
                x.eval_fitness(function); 
                y.eval_fitness(function);
            }

            /////////////////////////////////////////////////
            // replacement
            /////////////////////////////////////////////////
            // replacement worst
            {
                replace_worst(x);
                replace_worst(y);
            }

            // update whole populations; lazy method
            evaluate_fitness();
        }
        return sum_fitness/population.size();
    }

    void print(std::string filepath) {

        std::ofstream file(filepath.data());

        if(file.is_open()) {
            file << "chromosome ";
            file << "decoded ";
            file << "fitness\n";

            for(auto& p:population) {
                file << p.chromosome.to_string() << " ";
                file << p.chromosome.decode() << " ";
                file << p.fitness << std::endl;
            }
            file << "\n";
            file << "sum of fitness: " << sum_fitness/population.size() << "\n";
        }
        file.close();
    }

private:
    void initialize_population() {
        // randomly initialize population from chromosome
        for(size_t i=0; i<population.size(); i++) {
            population[i] = populate_t();
        }
    }

    void evaluate_fitness() {
        // evaluate fitness from
        sum_fitness = 0;
        for(auto& p:population) 
            sum_fitness += (p.eval_fitness(function));

        sum_probability = 0;
        for(auto& p:population)
            // todo
            sum_probability += p.eval_probability(sum_fitness, population.size());
    }

    //=========================================================
    // selection algorithm
    //=========================================================
    // roulette wheeling selection
    std::tuple<size_t, size_t> select_with_rouletewheeling() {

        auto select = [&]() {
            std::uniform_real_distribution<scalar_t> dis(0.0, 1.0);
            scalar_t random_number = dis(gen);

            size_t i=0;
            while(random_number>0.0) {
                random_number -= population[i].probability;
                i++;
            }
            return i;
        };

        return std::make_tuple(select(), select());
    }

    //=========================================================
    // crossover algorithm
    //=========================================================
    bool crossover(populate_t& x, populate_t& y) {

        std::uniform_real_distribution<double> prob(0.0, 1.0);
        if(prob_crossover<prob(gen)) return false;

        std::uniform_int_distribution<size_t> range(0, len-1);
        size_t r_start = range(gen);
        size_t r_end   = range(gen);

        x.chromosome.crossover(y.chromosome, r_start, r_end);
        return true;
    }

    //=========================================================
    // mutation algorithm
    //=========================================================
    bool mutation(populate_t& x) {
        // convetional mutation
        std::uniform_real_distribution<double> prob(0.0, 1.0);
        if(prob_mutation<prob(gen)) return false;

        // mutate
        x.chromosome.mutate();
        return true;
    }

    //=========================================================
    // replacement algorithm
    //=========================================================
    // replace parent
    void replace_parent(size_t index, const populate_t x) {
        population[index] = x;
    }
    // replace worst: the highest fitness value
    void replace_worst(const populate_t x) {
        auto elem = std::max_element(population.begin(), population.end(), [](populate_t a, populate_t b){ return a.fitness<b.fitness; });
        (*elem) = x;
    }

protected:
    function_t  function;
    std::vector<populate_t> population;
    scalar_t    sum_fitness;
    double      sum_probability;
    double      prob_crossover;
    double      prob_mutation;
};

//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__GENETIC_ALGORITHM__