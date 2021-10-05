#include <DataCollector.hpp>
#include <iostream>
#include <fstream>

#if __GNUC__ > 7
#include <filesystem>
#include <optional>
namespace filesystem = std::filesystem;
#else
#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;
#endif

BaseDataCollector::~BaseDataCollector() {
}

void BaseDataCollector::setParams(std::function<std::string(double)> get_path, unsigned int n, unsigned int m, unsigned int num_samples) {
    this->get_path = get_path;
    this->n = n;
    this->m = m;
    this->num_samples = num_samples;
}

void InfoDataCollector::reset() {
    this->avg.reset();
    this->referentially_useless = 0;
    this->largest_connected = 0;
    this->largest_connected_words = 0;
    this->largest_connected_meanings = 0;
}

void InfoDataCollector::collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti) {
    (void)opti;
    
    LexicalMatrixParameters<long double> params_sample;
    mat.get_parameters(params_sample);
    this->avg += params_sample;
    this->referentially_useless += static_cast<double>(mat.referentially_useless_words()) / this->n;

    std::pair<unsigned, unsigned> largest_connected_words_and_meanings = mat.largest_connected_component();
    double largest_connected_words = static_cast<double>(largest_connected_words_and_meanings.first);
    double largest_connected_meanings = static_cast<double>(largest_connected_words_and_meanings.second);
    this->largest_connected += (largest_connected_words + largest_connected_meanings) / (this->n + this->m);
    this->largest_connected_words += largest_connected_words / this->n;
    this->largest_connected_meanings += largest_connected_meanings / this->m;
}

void InfoDataCollector::record(double lambda) {
    LexicalMatrixParameters<long double> params = this->avg / this->num_samples;
    double ref_useless = this->referentially_useless / this->num_samples;
    double connected = this->largest_connected / this->num_samples;
    double connected_words = this->largest_connected_words / this->num_samples;
    double connected_meanings = this->largest_connected_meanings / this->num_samples;

    std::string path = this->get_path(lambda);
    
    std::ofstream file;
    file.open(path, std::ios::app);
    if (file.fail()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return;
    }
    
    file << lambda << "\t";
    file << params << "\t";
    file << ref_useless << "\t";
    file << connected << "\t";
    file << connected_words << "\t";
    file << connected_meanings << "\n";
    
    file << std::flush;
    file.close();
}

void LambdaDataCollector::setParams(std::function<std::string(double)> get_path, unsigned n, unsigned m, unsigned num_samples) {
    BaseDataCollector::setParams(get_path, n, m, num_samples);
    this->avg.resize(this->n);
}

void LambdaDataCollector::reset() {
    for (auto it = this->avg.begin(); it != this->avg.end(); it++) {
        it->reset();
    }
}

void LambdaDataCollector::collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti) {
    std::vector<LexicalMatrixStatistics<long double>> sample(this->n);
    std::vector<unsigned> idx(this->n);
    std::vector<unsigned> idx_mu(this->n);
    
    for (unsigned i=0; i<this->n; i++) {
        LexicalMatrixStatistics<long double> s;
        s.p_i = mat.get_ps(i);
        s.mu_i = mat.get_mu(i);
        s.a_i = opti.get_age(i);
        s.mu = s.mu_i;
        sample[i] = s;
    }

    for (unsigned i=0; i<this->n; i++) {
        idx[i] = i;
        idx_mu[i] = i;
    }

    std::sort(idx.begin(), idx.end(),
              [&sample](unsigned a, unsigned b) {
                  return sample[a].p_i > sample[b].p_i;
              });
    std::sort(idx_mu.begin(), idx_mu.end(),
              [&sample](unsigned a, unsigned b) {
                  return sample[a].mu > sample[b].mu;
              });
    
    for (unsigned i=0; i<this->n; i++) {
        unsigned j = idx[i];
        unsigned k = idx_mu[i];

        this->avg[i].p_i += sample[j].p_i;
        this->avg[i].mu_i += sample[j].mu_i;
        this->avg[i].a_i += sample[j].a_i;
        this->avg[i].mu += sample[k].mu;
    }
}

void LambdaDataCollector::record(double lambda) {
    std::string path = this->get_path(lambda);
    std::ofstream file;
    file.open(path, std::ios::app);
    if (file.fail()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return;
    }
    
    for (unsigned i=0; i<this->n; i++) {
        file << this->avg[i] / (double)this->num_samples;
    }

    file << std::flush;
    file.close();
}

void GraphDataCollector::reset() {
    this->edges.clear();
}

void GraphDataCollector::collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti) {
    (void)opti;
    for (unsigned j=0; j<mat.m; j++) {
        for (unsigned i=0; i<mat.n; i++) {
            if (mat.get_a(i, j)) {
                std::pair<unsigned, unsigned> pair = {i, j};
                this->edges.push_back(pair);
            }
        }
    }
    this->edges.push_back(std::nullopt);
}

void GraphDataCollector::record(double lambda) {
    filesystem::path dir = this->get_path(lambda);
    
    int n = 1;
    std::ofstream file;

    bool adjMat[this->m * this->n];
    std::fill(adjMat, adjMat+(this->m*this->n), false);
    
    for (auto optpair : this->edges) {
        if (optpair) {
            unsigned s = optpair.value().first;
            unsigned r = optpair.value().second;
            adjMat[r*this->n + s] = true;
        } else {
            auto filename = dir / ("graph_" + std::to_string(n) + ".csv");
            file.open(filename);
            if (file.fail()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }

            for (unsigned j=0; j<this->m; j++) {
                for (unsigned i=0; i<this->n; i++) {
                    file << adjMat[j*this->n + i];
                    if (i != this->n-1) {
                        file << "\t";
                    } else {
                        file << "\n";
                    }
                }
            }

            file.close();
            file.clear();
            std::fill(adjMat, adjMat+(this->m*this->n), false);
            n++;
        }
    }
}

void OptimizationTraceDataCollector::reset() {
    this->traces.clear();
}

void OptimizationTraceDataCollector::collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti) {
    (void)mat;
    this->traces.push_back(opti.get_trace());
}

void OptimizationTraceDataCollector::record(double lambda) {
    for (unsigned i=0; i<this->traces.size(); i++) {
        std::ostringstream ss;
        ss << "trace_lambda_" << lambda << "_sample_" << i << ".txt";

        filesystem::path dir = this->get_path(lambda);
        
        std::ofstream file;
        file.open(dir / ss.str());
        this->traces[i].to_file(file);
        file.close();
    }
}
