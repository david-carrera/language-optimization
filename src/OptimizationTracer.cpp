#include <OptimizationTracer.hpp>
#include <string>

#ifdef TRACE_MATRIX_OPTIMIZATION

void OptimizationTracer::initial_state(const BaseLexicalMatrix &mat,
                                       long double initial_energy) {
    this->initial_matrix.resize(mat.m);
    for (unsigned j=0; j<mat.m; j++) {
        this->initial_matrix[j].resize(mat.n);
        for (unsigned i=0; i<mat.n; i++) {
            this->initial_matrix[j][i] = mat.get_a(i, j);
        }
    }
    this->initial_energy = initial_energy;
    this->optimization_trace.push_back(optimization_step());
    this->optimization_trace.back().step = 0;
}

void OptimizationTracer::mutation(unsigned s, unsigned r) {
    this->optimization_trace.back().mutations.push_back({s, r});
}

void OptimizationTracer::step(bool energy_decreased, long double energy_after_mutations) {
    unsigned next_step = this->optimization_trace.back().step + 1;
    this->optimization_trace.back().energy_decreased = energy_decreased;
    this->optimization_trace.back().energy_after_mutations = energy_after_mutations;
    this->optimization_trace.push_back(optimization_step());
    this->optimization_trace.back().step = next_step;
}

void OptimizationTracer::end(bool early) {
    this->optimization_exited_early = early;
}

/// Handle nan and inf so they can be understood by Python's JSON parser
static std::string handle_nan(long double num) {
    if (std::isinf(num)) {
        if (std::signbit(num)) {
            return "-Infinite";
        } else {
            return "Infinite";
        }
    } else if (std::isnan(num)) {
        return "NaN";
    }

    return std::to_string(num);
}

void OptimizationTracer::to_file(std::ofstream &file) const {
    file << "{\n";
    file << "  \"initial_matrix\": [\n";

    for (unsigned j=0; j<this->initial_matrix.size(); j++) {
        file << "    [";
        for (unsigned i=0; i<this->initial_matrix[j].size(); i++) {
            if (this->initial_matrix[j][i]) {
                file << "1";
            } else {
                file << "0";
            }
            
            if (i < this->initial_matrix[j].size()-1) {
                file << ", ";
            }
        }
        file << "]";

        if (j < this->initial_matrix.size()-1) {
            file << ",";
        }
        file << "\n";
    }

    file << "  ],\n";
    file << "  \"initial_energy\": " << handle_nan(this->initial_energy) << ",\n";
    file << "  \"exited_early\": ";
    if (this->optimization_exited_early) {
        file << "true";
    } else {
        file << "false";
    }
    file << ",\n";

    file << "  \"trace\": [\n";
    for (unsigned i=0; i<this->optimization_trace.size(); i++) {
        const struct optimization_step step = optimization_trace[i];
        file << "    {\n";
        file << "      \"step\": " << step.step << ",\n";
        file << "      \"mutations\": [\n";
        for (unsigned j=0; j<step.mutations.size(); j++) {
            auto mut = step.mutations[j];
            file << "        {\"word\": " << mut.first << ", \"meaning\": " << mut.second << "}";
            if (j < step.mutations.size()-1) {
                file << ",";
            }
            file << "\n";
        }
        file << "      ],\n";
        file << "      \"energy_after_mutations\": " << handle_nan(step.energy_after_mutations) << ",\n";
        file << "      \"energy_decreased\": ";
        if (step.energy_decreased) {
            file << "true";
        } else {
            file << "false";
        }
        file << "\n";
        file << "    }";
        if (i < optimization_trace.size()-1) {
            file << ",";
        }
        file << "\n";
    }
    file << "  ]\n";
    file << "}\n";
}

#else

void OptimizationTracer::initial_state(const BaseLexicalMatrix &mat,
                                       long double initial_energy) {
    (void)mat;
    (void)initial_energy;
}

void OptimizationTracer::mutation(unsigned s, unsigned r) {
    (void)s;
    (void)r;
}

void OptimizationTracer::step(bool energy_decreased, long double energy_after_mutations) {
    (void)energy_decreased;
    (void)energy_after_mutations;
}

void OptimizationTracer::end(bool early) {
    (void)early;
}

void OptimizationTracer::to_file(std::ofstream &file) const {
    (void)file;
}

#endif
