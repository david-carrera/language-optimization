#ifndef PARAMETER_COLLECTION_H
#define PARAMETER_COLLECTION_H

#include <BaseLexicalMatrix.hpp>
#include <LexicalMatrixOptimizer.hpp>
#include <string>
#include <iostream>

#if __GNUC__ > 7
#include <filesystem>
namespace filesystem = std::filesystem;
#else
#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;
#endif

/**
 * Simulation parameters read from file or command line.
 *
 * All parameters read from file or command line are kept within this class
 * which can be queried for specific settings. The logic for verifying the
 * parameters is also within this class.
 */
class ParameterCollection {
private:
    /// Name of the yaml file
    const std::string config_yaml_name = "config.yaml";

    /// Name of the information theoretic csv
    const std::string info_csv_name = "information_theoretic.csv";

    /// Name of the file containing the parameters of the run and other information
    const std::string info_txt_name = "info.txt";

    /// Name of the directory where lambda specific statistics will be stored
    const std::string lambda_dir_name = "inside_lambda";

    /// Prefix of the CSVs for lambda specific statistics
    const std::string lambda_csv_name_prefix = "lambda_";

    /// Prefix of the CSVs for generated optimized graphs
    const std::string graphs_dir_name_prefix = "graph_visualization/graphs_lambda_";

    //! \cond
    struct config {
        struct graph {
            unsigned n;
            unsigned m;
            double phi;
            bool constant_prior_probability_of_meanings;
            bool unlinked_objects;
            struct pi {
                enum type {
                    uniform,
                    brokenstick,
                    geometric,
                    powerlaw,
                    custom,
                } type;
                union parameters {
                    struct geometric {
                        double p;
                    } geometric;
                    struct powerlaw {
                        double alpha;
                    } powerlaw;
                    struct custom {
                        std::size_t count;
                        double *probabilities;
                    } custom;
                } parameters;
            } pi;
        } graph;
        struct lambda {
            unsigned realizations;
            enum type {
                range,
                custom,
            } type;
            union parameters {
                struct range {
                    double min;
                    double max;
                    double step;
                } range;
                struct custom {
                    std::size_t count;
                    double *lambdas;
                } custom;
            } parameters;
        } lambda;
        struct initial_graph {
            enum type {
                random_edges,
                random_probability,
                complete_graph,
                bijection,
            } type;
            union parameters {
                struct random_edges {
                    unsigned edges;
                } random_edges;
                struct random_probability {
                    double probability;
                } random_probability;
            } parameters;
        } initial_graph;
        struct optimization {
            enum computation_of_energy {
                static_,
                dynamic,
            } computation_of_energy;
            struct mutations {
                enum type {
                    constant,
                    binomial,
                } type;
                union parameters {
                    struct constant {
                        unsigned count;
                    } constant;
                    struct binomial {
                        double probability;
                    } binomial;
                } parameters;
            } mutations;
            struct stop_condition {
                enum type {
                    weak,
                    strong,
                    custom,
                } type;
                union parameters {
                    struct weak {
                        int factor;
                    } weak;
                    struct strong {
                        int factor;
                    } strong;
                    struct custom {
                        unsigned long trials;
                    } custom;
                } parameters;
            } stop_condition;
        } optimization;
        struct misc {
            unsigned seed;
            double status_update_period;
            bool run_information_theoretical;
            bool run_inside_lambda;
            bool run_graph_visualization;
        } misc;
    };
    //! \endcond
    
    /// All information of the configuration file in a hierarchical structure
    struct config config;
    
    /// An instance of the lexical matrix needed for the given configuration which will be used throughout, living in the heap
    BaseLexicalMatrix *made_matrix = NULL;

    /// Has file been parsed?
    bool parsed = false;
    /// Base path where config will be read and files written
    filesystem::path base_path = ".";

    /// Parse command line
    void parse_cmd(int argc, char *argv[]);
    
    /// Validate configuration, possibly issuing warnings or exiting with errors
    void validate();

    /**
     * Assign a priori probabilities to vector.
     * @param[out]  pi  Vector of m elements to assign probabilities
     */
    void assign_pi(std::vector<long double> &pi);
    
public:
    ParameterCollection();

    ParameterCollection(const ParameterCollection &other) = delete;
    ParameterCollection &operator=(const ParameterCollection &other) = delete;

    ~ParameterCollection();

    void parse(int argc, char *argv[]);
    void lambdas(std::vector<double> &lambdas) const;

    filesystem::path base_dir() const;
    filesystem::path config_yaml() const;
    filesystem::path info_txt() const;
    filesystem::path lambda_dir() const;
    filesystem::path lambda_csv(double lambda) const;
    filesystem::path graphs_dir(double lambda) const;
    filesystem::path info_csv() const;

    bool show_status() const;
    bool generate_info() const;
    bool generate_lambdas() const;
    bool generate_graph() const;

    double status_update_period() const;
    
    unsigned total_runs() const;
    unsigned total_lambdas() const;
    unsigned elapsed_runs(double lambda) const;
    unsigned elapsed_lambdas(double lambda) const;
    unsigned realizations() const;

    unsigned n() const;
    unsigned m() const;

    LexicalMatrixOptimizer make_optimizer() const;
    BaseLexicalMatrix &make_matrix();

    void cleanup_matrix();
    friend std::ostream& operator<<(std::ostream& out, const ParameterCollection& p);
};

#endif /* PARAMETER_COLLECTION_H */
