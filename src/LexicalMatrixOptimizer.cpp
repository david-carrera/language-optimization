#include <LexicalMatrixOptimizer.hpp>

/**
 * Constructor.
 * @param[in]  seed           Seed for the random number generator.
 * @param[in]  stopCondition  Function to obtain the stop condition.
 * @param[in]  graphInit      Function to initialize the matrix
 * @param[in]  getMutations   Function to get the number of mutations to perform
 */
LexicalMatrixOptimizer::LexicalMatrixOptimizer(unsigned seed, OptimizerStopper stopCondition,
                                               OptimizerInitializer graphInit, OptimizerMutator getMutations,
                                               bool allow_unlinked_meanings) :
    prng(seed), stopCondition(stopCondition), graphInit(graphInit), getMutations(getMutations),
    allow_unlinked_meanings(allow_unlinked_meanings) {}

/**
 * Obtain the value of Omega
 * @param[in]  mat     The lexical matrix
 * @param[in]  lambda  Lambda to use
 */
inline long double LexicalMatrixOptimizer::get_Omega(const BaseLexicalMatrix &mat, double lambda) const {
    return -lambda*mat.get_ISR() + (1-lambda)*mat.get_HS();
}

/**
 * Obtain the minimum possible value of Omega. Assumes 0 <= lambda <= 1.
 * @param[in]  mat     The lexical matrix
 * @param[in]  lambda  Lambda to use
 */
inline long double LexicalMatrixOptimizer::get_minimum_Omega(const BaseLexicalMatrix &mat, double lambda) const {
    if (lambda <= 1 && lambda >= 0) {
        if (lambda <= 0.5) {
            return 0;
        } else {
            return (1 - 2*lambda) * std::log(std::min(mat.n, mat.m));
        }
    } else {
        return -std::numeric_limits<long double>::infinity();
    }
}

bool LexicalMatrixOptimizer::exit_early(long double Omega, BaseLexicalMatrix &mat,
                                        unsigned long failures_left) {
    if (Omega <= this->min_Omega + mat.epsilon) {
        this->tracer.step(true, Omega);
        this->tracer.end(true);
        this->update_word_is_connected(mat);
        this->update_ages(failures_left);
        mat.recalculate(true);
        return true;
    }
    return false;
}

/**
 * Update ages of all words. Increase the ones that are connected and reset to
 * 0 the ones that aren't.
 * @param[in]  word_is_connected  Vector specifying whether a word is connected.
 * @param[in]  amount             The amount by which to increase the ages.
 */
void LexicalMatrixOptimizer::update_ages(unsigned long amount) {
    for (unsigned s=0; s<this->word_is_connected.size(); s++) {
        if (this->word_is_connected[s]) {
            this->age[s] += amount;
        } else {
            this->age[s] = 0;
        }
    }
}

/**
 * Set the given vector to true if the i-th word is connected in mat, false
 * otherwise.
 * @param[in]   mat                Matrix with the connected/disconnected words
 * @param[out]  word_is_connected  Vector specifying whether a word is connected.
 */
void LexicalMatrixOptimizer::update_word_is_connected(const BaseLexicalMatrix &mat) {
    for (unsigned i=0; i<mat.n; i++) {
        this->word_is_connected[i] = mat.get_mu(i) > 0;
    }
}

/**
 * Initialize all variables needed for the next optimization as well as the
 * matrix.
 * @param[in,out]  basemat  Matrix to optimize.
 * @param[in]      lambda   Value of lambda to use for Omega calculation.
 */
void LexicalMatrixOptimizer::initialize_optimization(BaseLexicalMatrix &mat, double lambda) {
    this->graphInit(mat, this->prng);
    
    this->age.clear();
    this->age.resize(mat.n);
    
    this->word_is_connected.clear();
    this->word_is_connected.resize(mat.n);
    this->update_word_is_connected(mat);

    this->stop_failures = this->stopCondition(mat);
    this->total_failures = 0;
    
    this->best_Omega = this->get_Omega(mat, lambda);
    this->min_Omega = this->get_minimum_Omega(mat, lambda);
    
    this->tracer.initial_state(mat, best_Omega);
    
    mat.record(true);
    mat.save();
}

/**
 * Run the optimization process.
 * @param[in,out]  basemat  Matrix to optimize.
 * @param[in]      lambda   Value of lambda to use for Omega calculation.
 */
void LexicalMatrixOptimizer::optimize(BaseLexicalMatrix &mat, double lambda) {
    std::uniform_int_distribution<unsigned> words(0, mat.n-1);
    std::uniform_int_distribution<unsigned> meanings(0, mat.m-1);
    
    this->initialize_optimization(mat, lambda);
    if (this->exit_early(best_Omega, mat, stop_failures)) {
        return;
    }

    while (this->total_failures < this->stop_failures) {
        unsigned mutations = this->getMutations(mat, this->prng);
        for (unsigned i=0; i<mutations; i++) {
            unsigned s = words(prng);
            unsigned r = meanings(prng);
            mat.mutate(s, r);
            this->tracer.mutation(s, r);
        }
        
        if (!mat.valid(this->allow_unlinked_meanings)) {
            // equivalent to omega = \infty without having to recalculate it
            this->tracer.step(false, std::numeric_limits<long double>::infinity());
            mat.restore();
            this->total_failures++;
            continue;
        }

        mat.recalculate();
        long double current_Omega = this->get_Omega(mat, lambda);
        if (this->exit_early(current_Omega, mat, this->stop_failures - this->total_failures)) {
            return;
        }

        if (current_Omega + mat.epsilon < this->best_Omega) {
            this->tracer.step(true, current_Omega);
            this->update_ages(this->total_failures);
            this->update_word_is_connected(mat);
            mat.save();
            this->best_Omega = current_Omega;
            this->total_failures = 0;
            this->update_ages(1);
        } else {
            this->tracer.step(false, current_Omega);
            mat.restore();
            total_failures++;
        }
    }

    this->tracer.end(false);
    this->update_word_is_connected(mat);
    this->update_ages(this->stop_failures);
    mat.recalculate(true);
}

const OptimizationTracer &LexicalMatrixOptimizer::get_trace() const {
    return this->tracer;
}

/**
 * Create a function implementing a weak stop condition with the
 * given factor: floor(factor * n*m * log(n*m))
 * @param[in]  factor  Factor
 * @return             The function for the stop condition
 */
OptimizerStopper LexicalOptimizerCallbackFactory::StopCondition::weak(unsigned factor) {
    return [factor](const BaseLexicalMatrix &mat) -> unsigned long {
        unsigned s = mat.n * mat.m;
        return static_cast<unsigned long>(std::floor(factor * s * log(s)));
    };
}
 
/**
 * Create a function implementing a strong stop condition with the
 * given factor: floor(factor * A * log(A)) with A =
 * ((n*m-1)*n*m)/2
 * @param[in]  factor  Factor
 * @return             The function for the stop condition
 */
OptimizerStopper LexicalOptimizerCallbackFactory::StopCondition::strong(unsigned factor) {
    return [factor](const BaseLexicalMatrix &mat) -> unsigned long {
        unsigned s = mat.n * mat.m;
        // a \choose 2 = a! / 2! / (a-2)! = 2*3*4*...*(a-2)*(a-1)*a / 2 / 2*3*4*...*(a-2) = (a-1)*a / 2
        unsigned A = ((s-1) * s) / 2;
        // TODO: these numbers are gigantic (3e12 for n=m=400), better strong condition?
        return static_cast<unsigned long>(std::floor(factor * A * log(A)));
    };
}

/**
 * Create a function implementing a custom stop condition that will
 * stop after a specific number of trials.
 * @param  trials  The number of trials.
 * @return         The function for the stop condition
 */
OptimizerStopper LexicalOptimizerCallbackFactory::StopCondition::custom(unsigned long trials) {
    return [trials](const BaseLexicalMatrix &mat) -> unsigned long {
        (void)mat;
        return trials;
    };
}

/**
 * Create a function implementing the initialization of the graph
 * with the given number of random edges.
 * @param  edges  Number of edges after initialization.
 * @return        The function for the initialization
 */
OptimizerInitializer LexicalOptimizerCallbackFactory::Initializer::random_edges(unsigned edges) {
    return [edges](BaseLexicalMatrix &mat, RNG &rng) {
        mat.initialize_random_edges(rng, edges);
    };
}

/**
 * Create a function implementing the initialization of the graph with the
 * given number of random edges while ensuring that all meanings remain
 * connected. If it is impossible (the number of edges is less than the number
 * of meanings), the program is aborted.
 * @param  edges  Number of edges after initialization.
 * @return        The function for the initialization.
 */
OptimizerInitializer LexicalOptimizerCallbackFactory::Initializer::random_edges_connected(unsigned edges) {
    return [edges](BaseLexicalMatrix &mat, RNG &rng) {
        mat.initialize_random_edges(rng, edges, true);
    };
}

/**
 * Create a function implementing the initialization of the graph
 * with the given probability of edge creation.
 * @param  probability  Probability of creating an edge.
 * @return              The function for the initialization
 */
OptimizerInitializer LexicalOptimizerCallbackFactory::Initializer::random_probability(double probability) {
    return [probability](BaseLexicalMatrix &mat, RNG &rng) {
        mat.initialize_random(rng, probability);
    };
}

/**
 * Create a function implementing the initialization of the graph with the
 * given probability of edge creation while ensuring that all meanings remain
 * connected. If it is impossible (the expected number of edges is less than
 * the number of meanings), the program is aborted.
 * @param  probability  Probability of creating an edge.
 * @return              The function for the initialization
 */
OptimizerInitializer LexicalOptimizerCallbackFactory::Initializer::random_probability_connected(double probability) {
    return [probability](BaseLexicalMatrix &mat, RNG &rng) {
        mat.initialize_random(rng, probability, true);
    };
}

/**
 * Create a function implementing the initialization of the graph
 * as a complete graph.
 * @return  The function for the initialization
 */
OptimizerInitializer LexicalOptimizerCallbackFactory::Initializer::complete_graph() {
    return [](BaseLexicalMatrix &mat, RNG &rng) {
        (void)rng;
        mat.initialize_complete_graph(); 
    };
}

/**
 * Create a function implementing the initialization of the graph
 * as a bijectivev graph (every word connected to one meaning and
 * every meaning to one word).
 * @return  The function for the initialization
 */
OptimizerInitializer LexicalOptimizerCallbackFactory::Initializer::bijection() {
    return [](BaseLexicalMatrix &mat, RNG &rng) {
        (void)rng;
        mat.initialize_bijective_graph();
    };
}

/**
 * Create a function for obtaining the number of mutations
 * performed on the graph, a constant number of mutations.
 * @param  count  The constant number of mutations.
 * @return        The function for the mutations.
 */
OptimizerMutator LexicalOptimizerCallbackFactory::Mutator::constant(unsigned count) {
    return [count](const BaseLexicalMatrix &mat, RNG &rng) -> unsigned {
        (void)mat;
        (void)rng;
        return count;
    };
}

/**
 * Create a function for obtaining the number of mutations performed on the
 * graph, binomial number of mutations. The created function will never return
 * 0. The expected value for a random variable X following this left-truncated
 * binomial distribution with parameters n, p is: E[X] = n*p/(1-(1-p)^n)
 * @param  probability  The probability of mutation.
 * @return              The function for the mutations.
 */
OptimizerMutator LexicalOptimizerCallbackFactory::Mutator::binomial(double probability) {
    return [probability](const BaseLexicalMatrix &mat, RNG &rng) -> unsigned {
        std::binomial_distribution dist(mat.n * mat.m, probability);
        
        unsigned x;
        do {
            x = dist(rng); 
        } while (x == 0);
        return x;
    };
}
