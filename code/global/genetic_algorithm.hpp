#ifndef __GAUSSNEWTONS__
#define __GAUSSNEWTONS__

#include "method.h"

#include <random>
#include <bitset>

namespace numerical_optimization {

static std::random_device rd;
static std::mt19937 gen(rd());

template<typename pheno_t, size_t size>
struct chromosome_t {
    using bit_t = std::bitset<size>;

    chromosome_t(){
        // binary distribution
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

// = typename std::conditional<std::is_base_of<Eigen::EigenBase<vector_t>, vector_t>::value, vector_t, typename vector_t::Scalar>::type>
template<size_t len, typename vector_t, typename scalar_t>
class GeneticAlgorithm : Method {
public:
    using boundary_t = std::pair<scalar_t, scalar_t>;
    using function_t = std::function<scalar_t(const vector_t&)>;
    using populate_t = population_t<vector_t, scalar_t, len>;

    GeneticAlgorithm(function_t func, size_t size=20, double crossover_prob=0.7, double mutation_prob=0.1)
    :function(func),population(size),prob_crossover(crossover_prob),prob_mutation(mutation_prob),sum_fitness(0),sum_probability(0) {
        initialize_population();
        evaluate_fitness();
    }

    void run(size_t iteration) {

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
            std::cout << sum_fitness << std::endl;
        }
    }

    void print() {
        for(auto& p:population) {
            std::cout << p.chromosome.to_string() << std::endl;
            std::cout << p.fitness << std::endl;
            std::cout << p.probability << std::endl;
        }
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
#endif //__GAUSSNEWTONS__