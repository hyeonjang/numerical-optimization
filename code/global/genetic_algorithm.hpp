#ifndef __GAUSSNEWTONS__
#define __GAUSSNEWTONS__

#include "method.h"

#include <random>
#include <bitset>

namespace numerical_optimization {

static std::random_device rd;
static std::mt19937 gen(rd());

template<typename pheno_t, size_t size>
struct chromosome {
    using bit_t = std::bitset<size>;

    chromosome(){

        std::bernoulli_distribution d(0.5);

        for(size_t i=0; i<size; i++) {
            gene[i] = d(gen);
        }
    };

    // should be more flexible
    inline pheno_t decode() {
        return pheno_t(gene.to_ullong())/pheno_t(bit_t(ULLONG_MAX).to_ullong());
    }

    void mutate() {
        std::uniform_int_distribution<size_t> dis(0, size-1);
        size_t index = dis(gen);

        gene.flip(index);
    }

    void crossover(chromosome<pheno_t, size>& other, size_t end, size_t start=0) {
        // one point & two points crossover
        bit_t t_this  = this->gene;
        bit_t t_other = other.gene;
        for(size_t i=end; i>start; i--) {
            this->gene[i] = t_other[i];
            other.gene[i] = t_this[i];
        }
    }

    const std::string to_string() const {
        return gene.to_string();
    }

    std::bitset<size> gene;
};

template<typename vector_t, typename scalar_t = typename std::conditional<std::is_base_of<Eigen::EigenBase<vector_t>, vector_t>::value, vector_t, typename vector_t::Scalar>::type>
class GeneticAlgorithm : Method {
public:
    using boundary_t = std::pair<scalar_t, scalar_t>;
    using function_t = std::function<scalar_t(const vector_t&)>;
    using population_t = std::pair<chromosome<vector_t, 8>, scalar_t>;

    GeneticAlgorithm(function_t func, size_t size):function(func),population(size) {
        initialize_population();
    }
    GeneticAlgorithm(function_t func):function(func),population(20) {
        initialize_population();
        evaluate_fitness();

    };

    void initialize_population() {
        // randomly initialize population from chromosome
        for(size_t i=0; i<population.size(); i++) {
            population[i] = std::make_pair(chromosome<vector_t, 8>(), NULL);
        }
    }

    void evaluate_fitness() {

        for(auto& p:population) {
            auto fit = function(p.first.decode());
            p.second = fit;
        }
    }

    std::tuple<population_t, population_t> selection() {
        std::sort(population.begin(), population.end(), fitness_compare);
        return std::make_tuple(population[0], population[1]);
    }

    void crossover(population_t& x, population_t& y) {
        // std::cout << x.first.to_string() << " " << y.first.to_string() << std::endl;
        x.first.crossover(y.first, 4);
        // std::cout << x.first.to_string() << " " << y.first.to_string() << std::endl;
    }

    void mutation() {
        // convetional mutation

    }

    void run() {

        for(size_t i=0; i<50; i++) {

            population_t x, y;
            std::tie(x, y) = selection();
            // std::cout << x.first.to_string() << " " << y.first.to_string() << std::endl;
            std::cout << population[0].first.to_string() << " " << population[1].first.to_string() << std::endl;
            std::cout << population[0].second << " " << population[1].second << std::endl;
            crossover(x, y);
            population[0] = x;
            population[1] = y;
            evaluate_fitness();
            // std::cout << x.first.to_string() << " " << y.first.to_string() << std::endl;
            std::cout << population[0].first.to_string() << " " << population[1].first.to_string() << std::endl;
            std::cout << population[0].second << " " << population[1].second << std::endl;
        }

    }

    static inline bool fitness_compare(const population_t& l, const population_t& r) {
        return l.second < r.second;
    } 

protected:
    std::vector<population_t> population;
    function_t function;
};

//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__GAUSSNEWTONS__