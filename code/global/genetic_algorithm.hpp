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

    pheno_t decode() {
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

    GeneticAlgorithm(function_t func, size_t size):function(func),population(size),fitness(size) {
        initialize_population();
    }
    GeneticAlgorithm(function_t func):function(func),population(100),fitness(100) {
        initialize_population();
        evaluate_fitness();

    };

    void initialize_population() {

        for(size_t i=0; i<population.size(); i++) {
            population[i] = chromosome<vector_t, 8>();
        }

        for(auto p : population) {
            // std::cout << p.to_string() << std::endl;
            // p.mutate();
            // std::cout << p.to_string() << std::endl;
            // std::cout << p.decode() << std::endl;
        }
    }

    void evaluate_fitness() {

        for(auto p:population) {
            auto fit = function(p.decode());
            fitness.emplace_back(fit);
        }

        for(auto f:fitness) {
            // std::cout << f << std::endl;
        }
    }

    void crossover() {

        population[0].crossover(population[1], 2);
    }

    void mutation() {
        // convetional mutation

    }

    void run() {
        crossover();
    }

protected:
    std::vector<chromosome<vector_t, 8>> population;
    std::vector<scalar_t> fitness;
    function_t function;
};

//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
//////////////////////////////////////////////////// 
#endif //__GAUSSNEWTONS__