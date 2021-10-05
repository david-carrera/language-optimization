#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <BaseLexicalMatrix.hpp>
#include <LexicalMatrixOptimizer.hpp>
#include <LexicalMatrixStatistics.hpp>
#include <OptimizationTracer.hpp>
#include <optional>
#include <string>
#include <functional>

/**
 * Base class of a data collector. A data collector collects data throughout
 * the optimization process and writes it to a file.
 */
class BaseDataCollector {
protected:
    /// Function to obtain the path to write data to from a value of lambda.
    std::function<std::string(double)> get_path;
    
    /// Number of words
    unsigned n;

    /// Number of meanings
    unsigned m;
    
    /// Number of samples
    unsigned num_samples;
public:
    virtual ~BaseDataCollector();
    
    /**
     * Set parameters of the matrix that will be sampled.
     * @param[in]  path         Function to obtain the path to write data to. Depending on the data collector, might need to be a file or a directory.
     * @param[in]  n            Number of words in the matrix
     * @param[in]  m            Number of meanings in the matrix
     * @param[in]  num_samples  Number of samples that will be collected.
     */
    virtual void setParams(std::function<std::string(double)> get_path, unsigned n, unsigned m, unsigned num_samples);

    /**
     * Reset data collector, remove all collected data.
     */
    virtual void reset() = 0;

    /**
     * Collect data from the matrix and/or the optimizer.
     * @param[in]  mat   Matrix to collect data from.
     * @param[in]  opti  Optimizer object to collect data from.
     */
    virtual void collect(const BaseLexicalMatrix &mat,
                         const LexicalMatrixOptimizer &opti) = 0;

    /**
     * Write all collected data.
     * @param[in]  lambda  Value of lambda corresponding to the recorded data.
     */
    virtual void record(double lambda) = 0;
};

/**
 * Collects data from the LexicalMatrix, aggregates it and records it in a csv file.
 * Generates CSVs about the information theoretic data.
 */
class InfoDataCollector : public BaseDataCollector {
private:
    /// Sum of all collected matrix parameters
    LexicalMatrixParameters<Accumulator> avg;
    /// Sum of all collected referentially useless word counts
    double referentially_useless;
    /// Sum of all collected largest connected component counts
    double largest_connected;
    /// Sum of all collected largest connected component words
    double largest_connected_words;
    /// Sum of all collected largest connected component meanings
    double largest_connected_meanings;
public:
    virtual void reset();
    virtual void collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti);
    virtual void record(double lambda);
};

/**
 * Collects data from the LexicalMatrix, aggregates it and records it in a csv file.
 * Generates CSVs about the statistical properties of a specific lambda.
 */
class LambdaDataCollector : public BaseDataCollector {
private:
    /// Sum of all collected matrix statistics
    std::vector<LexicalMatrixStatistics<Accumulator>> avg;
public:
    virtual void setParams(std::function<std::string(double)> get_path, unsigned n, unsigned m, unsigned num_samples);
    virtual void reset();
    virtual void collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti);
    virtual void record(double lambda);
};

/**
 * Collect data from the LexicalMatrix to record it in a csv file.
 */
class GraphDataCollector : public BaseDataCollector {
private:
    /// Edges of the graph
    std::vector<std::optional<std::pair<unsigned,unsigned>>> edges;
public:
    virtual void reset();
    virtual void collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti);
    virtual void record(double lambda);
};

class OptimizationTraceDataCollector : public BaseDataCollector {
private:
    std::vector<OptimizationTracer> traces;
public:
    virtual void reset();
    virtual void collect(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti);
    virtual void record(double lambda);
};
    
#endif /* DATA_COLLECTOR_H */


    
    
    

