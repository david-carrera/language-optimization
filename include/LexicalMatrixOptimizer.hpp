#ifndef LEXICAL_MATRIX_OPTIMIZER_H
#define LEXICAL_MATRIX_OPTIMIZER_H

#include <Util.hpp>
#include <BaseLexicalMatrix.hpp>
#include <OptimizationTracer.hpp>
#include <functional>

/// Function to obtain the stop condition of the optimizer
typedef std::function<unsigned long(const BaseLexicalMatrix&)> OptimizerStopper;
/// Function to initialize a lexical matrix
typedef std::function<void(BaseLexicalMatrix&, RNG&)> OptimizerInitializer;
/// Function to obtain the number of mutations to apply to a lexical matrix
typedef std::function<unsigned(const BaseLexicalMatrix&, RNG&)> OptimizerMutator;

/**
 * Optimizer to manage and control the optimization process of a lexical matrix.
 */
class LexicalMatrixOptimizer {
private:
    /// Random number generator
    RNG prng;
    
    /// Account for the ages of words
    std::vector<unsigned long> age;
    /// Remember which words are connected in the matrix we're optimizing
    std::vector<bool> word_is_connected;

    /// Total failures recorded in the optimization process so far
    unsigned long total_failures;
    /// Number of total_failures needed to stop the current optimization process
    unsigned long stop_failures;

    /// Best value of Omega seen in the current optimization process
    long double best_Omega;
    /// Minimum possible value of Omega for the current matrix and lambda
    long double min_Omega;

    /// Function to obtain the stop condition
    OptimizerStopper stopCondition;
    /// Function to initialize the matrix
    OptimizerInitializer graphInit;
    /// Function to get the number of mutations to pefrorm
    OptimizerMutator getMutations;

    /// Whether to allow unlinked meanings during optimization
    bool allow_unlinked_meanings;

    /// Used in debug mode to store details about every optimization step
    OptimizationTracer tracer;
    
    long double get_Omega(const BaseLexicalMatrix &mat, double lambda) const;
    long double get_minimum_Omega(const BaseLexicalMatrix &mat, double lambda) const;
    void update_ages(unsigned long amount);
    void update_word_is_connected(const BaseLexicalMatrix &mat);

    bool exit_early(long double Omega, BaseLexicalMatrix &mat, unsigned long failures_left);

    void initialize_optimization(BaseLexicalMatrix &mat, double lambda);

public:
    LexicalMatrixOptimizer(unsigned seed, OptimizerStopper stopCondition, OptimizerInitializer graphInit, OptimizerMutator getMutations,
                           bool allow_unlinked_meanings = true);

    void optimize(BaseLexicalMatrix &mat, double lambda);

    const OptimizationTracer &get_trace() const;

    /**
     * Age of a word
     * @param[in]  s  Word
     * @return        Age of the word
     */
    inline unsigned long get_age(unsigned s) const {
        return this->age[s];
    }

    /**
     * Set the seed of the random number generator
     * @param[in]  seed  Seed
     */
    inline void seed(unsigned seed) {
        this->prng.seed(seed);
    }
};

namespace LexicalOptimizerCallbackFactory {
    namespace StopCondition {
        OptimizerStopper weak(unsigned factor);
        OptimizerStopper strong(unsigned factor);
        OptimizerStopper custom(unsigned long trials);
    }
    
    namespace Initializer {
        OptimizerInitializer random_edges(unsigned edges);
        OptimizerInitializer random_edges_connected(unsigned edges);
        OptimizerInitializer random_probability(double probability);
        OptimizerInitializer random_probability_connected(double probability);
        OptimizerInitializer complete_graph();
        OptimizerInitializer bijection();
    }
    
    namespace Mutator {
        OptimizerMutator constant(unsigned count);
        OptimizerMutator binomial(double probability);
    }
}

#endif /* LEXICAL_MATRIX_OPTIMIZER_H */
