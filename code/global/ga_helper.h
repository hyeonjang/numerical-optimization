#ifndef __GA_HELPER__
#define __GA_HELPER__

#include "method.h"

#include <random>
#include <bitset>
#include <fstream>

namespace numerical_optimization {

template<typename pheno_t, size_t size>
struct chromosome_t {
    using bit_t = std::bitset<size>;

    chromosome_t(){
        // binary distribution
        std::bernoulli_distribution d(0.5);
        for(size_t i=0; i<size; i++)
            gene[i] = d(gen);
    };

    inline pheno_t decode() {
        return pheno_t(gene.to_ullong())/pheno_t(bit_t(ULLONG_MAX).to_ullong());
    }

    void mutate() {
        std::uniform_int_distribution<size_t> dis(0, size-1);
        size_t index = dis(gen);

        gene.flip(index);
    }

    void crossover(chromosome_t<pheno_t, size>& other, size_t end, size_t start=0) {
        
        size_t min, max;
        std::tie(min, max) = std::minmax(end, start);

        // one point & two points crossover
        bit_t t_this  = this->gene;
        bit_t t_other = other.gene;
        for(size_t i=max; i>min; i--) {
            this->gene[i] = t_other[i];
            other.gene[i] = t_this[i];
        }
    }

    const std::string to_string() const {
        return gene.to_string();
    }
private:
    bit_t gene;
};

template<typename vector_t, typename scalar_t, size_t size>
struct population_t {

    using function_t = std::function<scalar_t(const vector_t&)>;

    population_t():chromosome(chromosome_t<vector_t, size>()),fitness(0),probability(0) {}

    scalar_t eval_fitness(const function_t& function){
        return fitness = function(chromosome.decode());
    };

    double eval_probability(const scalar_t& sum_fitness, const size_t& total_size) {
        return probability = ((sum_fitness-fitness)/sum_fitness)/(total_size-1);
    }

    chromosome_t<vector_t, size> chromosome;
    scalar_t fitness;
    double   probability;
};
//// ///////////////////////////////////////////////
}/// the end of namespace numerical_optimization ///
////////////////////////////////////////////////////
#endif //__GA_HELPER__
