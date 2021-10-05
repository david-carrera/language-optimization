#ifndef OPTIMIZATION_TRACER_H
#define OPTIMIZATION_TRACER_H

#include <BaseLexicalMatrix.hpp>
#include <fstream>

struct optimization_step {
    unsigned step;
    std::vector<std::pair<unsigned, unsigned>> mutations;
    long double energy_after_mutations;
    bool energy_decreased;
};

class OptimizationTracer {
private:
    std::vector<std::vector<bool>> initial_matrix;
    long double initial_energy;
    std::vector<optimization_step> optimization_trace;
    bool optimization_exited_early;
public:
    void initial_state(const BaseLexicalMatrix &mat, long double initial_energy);
    void mutation(unsigned s, unsigned r);
    void step(bool energy_decreased, long double energy_after_mutations);
    void end(bool early);
    void to_file(std::ofstream &file) const;
};

#endif /* OPTIMIZATION_TRACER_H */
