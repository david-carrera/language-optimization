#include <ProbabilityDistributionGeneration.hpp>
#include <ParameterCollection.hpp>
#include <FirstModelDynamicLexicalMatrix.hpp>
#include <SecondModelDynamicLexicalMatrix.hpp>
#include <yaml-cpp/yaml.h>
#include <cassert>

#if __GNUC__ <= 7
namespace filesystem = std::experimental::filesystem;
#else
namespace filesystem = std::filesystem;
#endif

[[ noreturn ]] void print_yaml_error_and_exit(std::string path, std::string text) {
    std::cerr << "error parsing the file '" << path << "', the yaml engine reports: " << text << std::endl;
    exit(EXIT_FAILURE);
}

[[ noreturn ]] void print_parse_error_and_exit(std::string text) {
    std::cerr << "error parsing configuration file: " << text << std::endl;
    exit(EXIT_FAILURE);
}

[[ noreturn ]] void print_validation_error_and_exit(std::string text) {
    std::cerr << "error in configuration file: " << text << std::endl;
    exit(EXIT_FAILURE);
}
void print_validation_warning(std::string text) {
    std::cerr << "WARNNG in configuration file: " << text << std::endl;
}

[[ noreturn ]] void print_usage_and_exit(std::string program_name = "program") {
    std::cerr << "Usage:\n\t" + program_name + " directory_path" << std::endl;
    exit(EXIT_FAILURE);
}

#define yaml_ensure_basic(name, node)                   \
    if (!node.IsDefined()) {                            \
        print_parse_error_and_exit(name " undefined");  \
    }

#define yaml_ensure_map(name, node)                             \
    yaml_ensure_basic(name, node);                              \
                                                                \
    if (!node.IsMap()) {                                        \
        print_parse_error_and_exit(name " not a map or empty"); \
    }

#define yaml_ensure_basic_scalar(name, node)                    \
    yaml_ensure_basic(name, node);                              \
                                                                \
    if (!node.IsScalar()) {                                     \
        print_parse_error_and_exit(name " not a scalar");       \
    }

int find_idx(std::string value, std::vector<std::string> values, std::string name) {
    auto it = std::find(values.begin(), values.end(), value);
    if (it == values.end()) {
        print_parse_error_and_exit(name + " is none of the allowed values (" + vector_join(values, "|") + ")");
    }
    return static_cast<int>(it - values.begin());
}

#define yaml_ensure_enum(name, node, values, result, T) \
    yaml_ensure_basic_scalar(name, node);               \
                                                        \
    result = static_cast<T>(                            \
        find_idx(node.as<std::string>(), values, name));

#define yaml_ensure_conversion(name, node, result, T)   \
    try {                                               \
        result = node.as<T>();                          \
    } catch (const YAML::TypedBadConversion<T> &e) {    \
        print_parse_error_and_exit(name " is not " +    \
                             type_name<T>());           \
    }

#define yaml_ensure_scalar(name, node, result, T)       \
    yaml_ensure_basic_scalar(name, node);               \
                                                        \
    yaml_ensure_conversion(name, node, result, T)

#define yaml_ensure_nullable_scalar(name, node, result, dflt, T)        \
    yaml_ensure_basic(name, node);                                      \
                                                                        \
    if (node.IsNull()) {                                                \
        result = dflt;                                                  \
    } else {                                                            \
        if (!node.IsScalar()) {                                         \
            print_parse_error_and_exit(name " not a scalar");           \
        }                                                               \
        yaml_ensure_conversion(name, node, result, T);                  \
    }

#define yaml_ensure_sequence(name, node, sequence, count, T)    \
    yaml_ensure_basic(name, node)                               \
                                                                \
    if (!node.IsSequence()) {                                   \
        print_parse_error_and_exit(name " not a sequence");     \
    }                                                           \
                                                                \
    count = node.size();                                        \
    sequence = new T[count];                                    \
                                                                \
    for (unsigned i=0; i<node.size(); i++) {                    \
        yaml_ensure_conversion(name, node[i], sequence[i], T);  \
    }

void ParameterCollection::validate() {
    if (this->config.graph.phi < 0) {
        print_validation_error_and_exit("graph.phi must be greater than or equal to zero");
    }
    
    if (this->config.graph.constant_prior_probability_of_meanings) {
        long double min_p;
        switch (this->config.graph.pi.type) {
        case ParameterCollection::config::graph::pi::type::uniform:
            break;
        case ParameterCollection::config::graph::pi::type::brokenstick:
            break;
        case ParameterCollection::config::graph::pi::type::geometric:
            if (this->config.graph.pi.parameters.geometric.p <= 0 ||
                this->config.graph.pi.parameters.geometric.p > 1) {
                print_validation_error_and_exit("graph.pi.parameters.p must be in (0,1]");
            }
            
            min_p = ProbabilityDistributionGeneration::geometric_smallest_probability(
                this->config.graph.pi.parameters.geometric.p, this->config.graph.m
            );
            if (min_p <= default_epsilon) {
                print_validation_warning(
                    "geometric distribution will generate probabilities with magnitude "
                    "lower than the tolerated numerical error (smallest probability that "
                    "will be generated is " + std::to_string(min_p) + " but the error "
                    "tolerance is " + std::to_string(default_epsilon) + ")"
                );
            }
            break;
        case ParameterCollection::config::graph::pi::type::powerlaw:
            if (this->config.graph.pi.parameters.powerlaw.alpha <= 0) {
                print_validation_error_and_exit("graph.pi.parameters.alpha must be greater than zero");
            }
            
            min_p = ProbabilityDistributionGeneration::powerlaw_smallest_probability(
                this->config.graph.pi.parameters.powerlaw.alpha, this->config.graph.m
            );
            if (min_p <= default_epsilon) {
                print_validation_warning(
                    "powerlaw distribution will generate probabilities with magnitude "
                    "lower than the tolerated numerical error (smallest probability that "
                    "will be generated is " + std::to_string(min_p) + " but the error "
                    "tolerance is " + std::to_string(default_epsilon) + ")"
                );
            }
            break;
        case ParameterCollection::config::graph::pi::type::custom:
            if (this->config.graph.pi.parameters.custom.count != (unsigned)this->config.graph.m) {
                print_validation_error_and_exit("exactly graph.m values must be specified in graph.pi.parameters.probabilities");
            }
            for (unsigned i=0; i<this->config.graph.pi.parameters.custom.count; i++) {
                long double p = this->config.graph.pi.parameters.custom.probabilities[i];
                if (p <= default_epsilon) {
                    print_validation_error_and_exit(
                        "one of the custom probabilities is less than or too close to "
                        "zero (" + std::to_string(p) + " but error tolerance is " +
                        std::to_string(default_epsilon) + ")"
                    );
                }
            }
            break;
        default:
            break;
        }
    }
    
    if (this->config.lambda.realizations <= 0) {
        print_validation_error_and_exit("lambda.realizations must be greater than zero");
    }
    switch (this->config.lambda.type) {
    case ParameterCollection::config::lambda::type::range:
        if (this->config.lambda.parameters.range.min < 0 ||
            this->config.lambda.parameters.range.min > 1) {
            print_validation_warning("lambda.parameters.min should be in [0,1]");
        }
        if (this->config.lambda.parameters.range.max < 0 ||
            this->config.lambda.parameters.range.max > 1) {
            print_validation_warning("lambda.parameters.max should be in [0,1]");
        }
        if (this->config.lambda.parameters.range.max < this->config.lambda.parameters.range.min) {
            print_validation_error_and_exit("lambda.parameters.max must be greater than or equal to lambda.parameters.min");
        }
        if (this->config.lambda.parameters.range.step <= 0 ||
            this->config.lambda.parameters.range.step > this->config.lambda.parameters.range.max - this->config.lambda.parameters.range.min) {
            print_validation_error_and_exit("lambda.parameters.step must be in (0, max-min]");
        }
        break;
    case ParameterCollection::config::lambda::type::custom:
        if (this->config.lambda.parameters.custom.count == 0) {
            print_validation_error_and_exit("no values specified in lambda.parameters.lambdas");
        }
        break;
    default:
        break;
    }

    switch (this->config.initial_graph.type) {
    case ParameterCollection::config::initial_graph::type::random_edges:
        if (this->config.initial_graph.parameters.random_edges.edges <= 0 ||
            this->config.initial_graph.parameters.random_edges.edges > this->config.graph.n * this->config.graph.m) {
            print_validation_error_and_exit("initial_graph.parameters.edges must be in (0, n*m]");
        }
        break;
    case ParameterCollection::config::initial_graph::type::random_probability:
        if (this->config.initial_graph.parameters.random_probability.probability <= 0 ||
            this->config.initial_graph.parameters.random_probability.probability > 1) {
            print_validation_error_and_exit("initial_graph.parameters.probability must be in (0,1]");
        }
        if (this->config.initial_graph.parameters.random_probability.probability * this->config.graph.n * this->config.graph.m < 1) {
            print_validation_warning("initial_graph.parameters.probability * n * m is smaller than 1, initial graph is likely to have no edges");
        }
        break;
    case ParameterCollection::config::initial_graph::type::complete_graph:
    case ParameterCollection::config::initial_graph::type::bijection:
        break;
    default:
        break;
    }

    switch (this->config.optimization.mutations.type) {
    case ParameterCollection::config::optimization::mutations::type::constant:
        if (this->config.optimization.mutations.parameters.constant.count <= 0) {
            print_validation_error_and_exit("optimization.mutations.parameters.count must be greater than 0");
        }
        break;
    case ParameterCollection::config::optimization::mutations::type::binomial:
        if (this->config.optimization.mutations.parameters.binomial.probability <= 0 ||
            this->config.optimization.mutations.parameters.binomial.probability > 1) {
            print_validation_error_and_exit("optimization.mutations.parameters.probability must be in (0,1]");
        }
        break;
    default:
        break;
    }

    switch (this->config.optimization.stop_condition.type) {
    case ParameterCollection::config::optimization::stop_condition::type::weak:
        if (this->config.optimization.stop_condition.parameters.weak.factor < 0) {
            print_validation_error_and_exit("optimizaton.stop_condition.parameters.factor must be greater than or equal to 0");
        }
        break;
    case ParameterCollection::config::optimization::stop_condition::type::strong:
        if (this->config.optimization.stop_condition.parameters.strong.factor < 0) {
            print_validation_error_and_exit("optimizaton.stop_condition.parameters.factor must be greater than or equal to 0");
        }
        break;
    case ParameterCollection::config::optimization::stop_condition::type::custom:
        break;
    default:
        break;
    }
}

void ParameterCollection::assign_pi(std::vector<long double> &pi) {
    const unsigned m = (unsigned)this->config.graph.m;
    
    switch (this->config.graph.pi.type) {
    case ParameterCollection::config::graph::pi::type::uniform:
        pi = ProbabilityDistributionGeneration::uniform(m);
        break;
    case ParameterCollection::config::graph::pi::type::brokenstick:
        pi = ProbabilityDistributionGeneration::brokenstick(m);
        break;
    case ParameterCollection::config::graph::pi::type::geometric:
    {
        long double p = this->config.graph.pi.parameters.geometric.p;
        pi = ProbabilityDistributionGeneration::geometric(p, m);
    }
    break;
    case ParameterCollection::config::graph::pi::type::powerlaw:
    {
        long double alpha = this->config.graph.pi.parameters.powerlaw.alpha;
        pi = ProbabilityDistributionGeneration::powerlaw(alpha, m);
    }
    break;
    case ParameterCollection::config::graph::pi::type::custom:
    {
        double *p = this->config.graph.pi.parameters.custom.probabilities;
        long double sum = 0;
        for (unsigned i=0; i<m; i++) {
            pi[i] = p[i];
            sum += p[i];
        }
        for (unsigned i=0; i<m; i++) {
            pi[i] /= sum;
        }
    }
    break;
    default:
        assert(false);
    }
}

/// Constructor
ParameterCollection::ParameterCollection() {
    this->config.graph.pi.parameters.custom.probabilities = NULL;
    this->config.lambda.parameters.custom.lambdas = NULL;
}

/// Destructor
ParameterCollection::~ParameterCollection() {
    this->cleanup_matrix();
    if (this->parsed) {
        if (this->config.graph.pi.type == ParameterCollection::config::graph::pi::custom) {
            delete[] this->config.graph.pi.parameters.custom.probabilities;
            this->config.graph.pi.parameters.custom.probabilities = NULL;
            this->config.graph.pi.parameters.custom.count = 0;
        }
        if (this->config.lambda.type == ParameterCollection::config::lambda::type::custom) {
	    delete[] this->config.lambda.parameters.custom.lambdas;
	    this->config.lambda.parameters.custom.lambdas = NULL;
	    this->config.lambda.parameters.custom.count = 0;
        }
    }
}

void ParameterCollection::parse_cmd(int argc, char *argv[]) {
    if (argc < 1) {
        print_usage_and_exit();
    }
    std::string name = argv[0];
    
    if (argc != 2) {
        print_usage_and_exit(name);
    }

    this->base_path = argv[1];
}

/**
 * Read parameters from command line, including configuration file residing
 * in the given root directory.
 * @param[in]  argc  argc given by OS
 * @param[in]  argv  argv given by OS
 */
void ParameterCollection::parse(int argc, char *argv[]) {
    this->parse_cmd(argc, argv);

    std::string config_filename = this->config_yaml();
    YAML::Node yaml;
    try {
        yaml = YAML::LoadFile(config_filename);
    } catch (YAML::Exception &e) {
        print_yaml_error_and_exit(config_filename, e.what());
    }

    std::vector<std::string> values;

    if (!yaml.IsMap()) {
        print_parse_error_and_exit("file is empty");
    }
    
    yaml_ensure_map("graph", yaml["graph"]);
    yaml_ensure_scalar("graph.n", yaml["graph"]["n"], this->config.graph.n, unsigned);
    yaml_ensure_scalar("graph.m", yaml["graph"]["m"], this->config.graph.m, unsigned);
    yaml_ensure_scalar("graph.phi", yaml["graph"]["phi"], this->config.graph.phi, double);
    yaml_ensure_scalar("graph.constant_prior_probability_of_meanings",
                       yaml["graph"]["constant_prior_probability_of_meanings"],
                       this->config.graph.constant_prior_probability_of_meanings, bool);
    yaml_ensure_scalar("graph.unlinked_objects", yaml["graph"]["unlinked_objects"],
                       this->config.graph.unlinked_objects, bool);
    if (this->config.graph.constant_prior_probability_of_meanings) {
        yaml_ensure_map("graph.pi", yaml["graph"]["pi"]);
        values = {"uniform", "brokenstick", "geometric", "powerlaw", "custom"};
        yaml_ensure_enum("graph.pi.type", yaml["graph"]["pi"]["type"], values,
                         this->config.graph.pi.type,
                         enum ParameterCollection::config::graph::pi::type);
        switch (this->config.graph.pi.type) {
        case ParameterCollection::config::graph::pi::type::uniform:
            break;
        case ParameterCollection::config::graph::pi::type::brokenstick:
            break;
        case ParameterCollection::config::graph::pi::type::geometric:
            yaml_ensure_scalar("graph.pi.parameters.p",
                               yaml["graph"]["pi"]["parameters"]["p"],
                               this->config.graph.pi.parameters.geometric.p,
                               double);
            break;
        case ParameterCollection::config::graph::pi::type::powerlaw:
            yaml_ensure_scalar("graph.pi.parameters.alpha",
                               yaml["graph"]["pi"]["parameters"]["alpha"],
                               this->config.graph.pi.parameters.powerlaw.alpha,
                               double);
            break;
        case ParameterCollection::config::graph::pi::type::custom:
            yaml_ensure_sequence("graph.pi.parameters.probabilities",
                                 yaml["graph"]["pi"]["parameters"]["probabilities"],
                                 this->config.graph.pi.parameters.custom.probabilities,
                                 this->config.graph.pi.parameters.custom.count, double);
            break;
        default:
            print_parse_error_and_exit("bad graph.pi.type");
        }
    }
    
    yaml_ensure_map("lambda", yaml["lambda"]);
    yaml_ensure_scalar("lambda.realizations", yaml["lambda"]["realizations"],
                       this->config.lambda.realizations, unsigned);
    
    values = {"range", "custom"};
    yaml_ensure_enum("lambda.type", yaml["lambda"]["type"], values, this->config.lambda.type,
                     enum ParameterCollection::config::lambda::type);

    yaml_ensure_map("lambda.parameters", yaml["lambda"]["parameters"]);
    switch (this->config.lambda.type) {
    case ParameterCollection::config::lambda::type::range:
    {
        yaml_ensure_scalar("lambda.parameters.min", yaml["lambda"]["parameters"]["min"],
                           this->config.lambda.parameters.range.min, double);
        yaml_ensure_scalar("lambda.parameters.max", yaml["lambda"]["parameters"]["max"],
                           this->config.lambda.parameters.range.max, double);
        yaml_ensure_scalar("lambda.parameters.step", yaml["lambda"]["parameters"]["step"],
                           this->config.lambda.parameters.range.step, double);
        break;
    }
    case ParameterCollection::config::lambda::type::custom:
    {
        yaml_ensure_sequence("lambda.parameters.lambdas", yaml["lambda"]["parameters"]["lambdas"],
                             this->config.lambda.parameters.custom.lambdas,
                             this->config.lambda.parameters.custom.count, double);
        break;
    }
    default:
        print_parse_error_and_exit("bad lambda.type");
    }

    yaml_ensure_map("initial_graph", yaml["initial_graph"]);

    values = {"random_edges", "random_probability", "complete_graph", "bijection"};
    yaml_ensure_enum("initial_graph.type", yaml["initial_graph"]["type"], values,
                     this->config.initial_graph.type, enum ParameterCollection::config::initial_graph::type);

    switch (this->config.initial_graph.type) {
    case ParameterCollection::config::initial_graph::type::random_edges:
    {
        yaml_ensure_map("initial_graph.parameters", yaml["initial_graph"]["parameters"]);
        yaml_ensure_scalar("initial_graph.parameters.edges", yaml["initial_graph"]["parameters"]["edges"],
                           this->config.initial_graph.parameters.random_edges.edges, unsigned);
        break;
    }
    case ParameterCollection::config::initial_graph::type::random_probability:
    {
        yaml_ensure_map("initial_graph.parameters", yaml["initial_graph"]["parameters"]);
        yaml_ensure_scalar("initial_graph.parameters.probability",
                           yaml["initial_graph"]["parameters"]["probability"],
                           this->config.initial_graph.parameters.random_probability.probability, double);
        break;
    }
    case ParameterCollection::config::initial_graph::type::complete_graph:
    case ParameterCollection::config::initial_graph::bijection:
        break;
    default:
        print_parse_error_and_exit("bad initial_graph.type");
    }

    yaml_ensure_map("optimization", yaml["optimization"]);
    
    values = {"static", "dynamic"};
    yaml_ensure_enum("optimization.computation_of_energy", yaml["optimization"]["computation_of_energy"], values,
                     this->config.optimization.computation_of_energy,
                     enum ParameterCollection::config::optimization::computation_of_energy);
    
    yaml_ensure_map("optimization.mutations", yaml["optimization"]["mutations"]);
    
    values = {"constant", "binomial"};
    yaml_ensure_enum("optimization.mutations.type", yaml["optimization"]["mutations"]["type"], values,
                     this->config.optimization.mutations.type, enum ParameterCollection::config::optimization::mutations::type);

    yaml_ensure_map("optimization.mutations.parameters", yaml["optimization"]["mutations"]["parameters"]);
    switch (this->config.optimization.mutations.type) {
    case ParameterCollection::config::optimization::mutations::type::constant:
    {
        yaml_ensure_scalar("optimization.mutations.parameters.count", yaml["optimization"]["mutations"]["parameters"]["count"],
                           this->config.optimization.mutations.parameters.constant.count, unsigned);
        break;
    }
    case ParameterCollection::config::optimization::mutations::type::binomial:
    {
        yaml_ensure_scalar("optimization.mutations.parameters.probability",
                           yaml["optimization"]["mutations"]["parameters"]["probability"],
                           this->config.optimization.mutations.parameters.binomial.probability, double);
        break;
    }
    default:
        print_parse_error_and_exit("bad optimization.mutations.type");
    }

    yaml_ensure_map("optimization.stop_condition", yaml["optimization"]["stop_condition"]);

    values = {"weak", "strong", "custom"};
    yaml_ensure_enum("optimization.stop_condition.type", yaml["optimization"]["stop_condition"]["type"], values,
                     this->config.optimization.stop_condition.type,
                     enum ParameterCollection::config::optimization::stop_condition::type);

    yaml_ensure_map("optimization.stop_condition.parameters", yaml["optimization"]["stop_condition"]["parameters"]);
    switch(this->config.optimization.stop_condition.type) {
    case ParameterCollection::config::optimization::stop_condition::type::weak:
    {
        yaml_ensure_scalar("optimization.stop_condition.parameters.factor",
                           yaml["optimization"]["stop_condition"]["parameters"]["factor"],
                           this->config.optimization.stop_condition.parameters.weak.factor, int);
        break;
    }
    case ParameterCollection::config::optimization::stop_condition::type::strong:
    {
        yaml_ensure_scalar("optimization.stop_condition.parameters.factor",
                           yaml["optimization"]["stop_condition"]["parameters"]["factor"],
                           this->config.optimization.stop_condition.parameters.strong.factor, int);
        break;
    }
    case ParameterCollection::config::optimization::stop_condition::type::custom:
    {
        yaml_ensure_scalar("optimization.stop_condition.parameters.trials",
                           yaml["optimization"]["stop_condition"]["parameters"]["trials"],
                           this->config.optimization.stop_condition.parameters.custom.trials, unsigned);
        break;
    }
    default:
        print_parse_error_and_exit("bad optimization.stop_condition.type");
    }

    yaml_ensure_map("misc", yaml["misc"]);
    yaml_ensure_nullable_scalar("misc.seed", yaml["misc"]["seed"], this->config.misc.seed, std::random_device()(), unsigned int);
    yaml_ensure_scalar("misc.status_update_period", yaml["misc"]["status_update_period"],
                       this->config.misc.status_update_period, double);
    yaml_ensure_scalar("misc.run_information_theoretical", yaml["misc"]["run_information_theoretical"],
                       this->config.misc.run_information_theoretical, bool);
    yaml_ensure_scalar("misc.run_inside_lambda", yaml["misc"]["run_inside_lambda"],
                       this->config.misc.run_inside_lambda, bool);
    yaml_ensure_scalar("misc.run_graph_visualization", yaml["misc"]["run_graph_visualization"],
                       this->config.misc.run_graph_visualization, bool);
    
    this->parsed = true;
    this->validate();
}

/**
 * Fill given vector with all values of lambda to evaluate.
 * @param[out]  lambdas  Vector of lambdas
 */
void ParameterCollection::lambdas(std::vector<double> &lambdas) const {
    switch (this->config.lambda.type) {
    case ParameterCollection::config::lambda::type::range:
    {
        double min = this->config.lambda.parameters.range.min;
        double step = this->config.lambda.parameters.range.step;
        unsigned n = this->total_lambdas();
        for (unsigned i=0; i<n; i++) {
            lambdas.push_back(min + step*i);
        }
        break;
    }
    case ParameterCollection::config::lambda::type::custom:
    {
        std::size_t count = this->config.lambda.parameters.custom.count;
        double *values = this->config.lambda.parameters.custom.lambdas;
        for (std::size_t i=0; i<count; i++) {
            lambdas.push_back(values[i]);
        }
        break;
    }
    default:
        break;
    }
}

/**
 * Get path to the base directory, where the yaml file resides.
 * @return  Path to the base directory
 */
filesystem::path ParameterCollection::base_dir() const {
    return this->base_path;
}

/**
 * Get path to yaml file.
 * @return  Path to configuration file
 */
filesystem::path ParameterCollection::config_yaml() const {
    return this->base_path / this->config_yaml_name;
}

/**
 * Get path to log file.
 * @return  Path to log file
 */
filesystem::path ParameterCollection::info_txt() const {
    return this->base_path / this->info_txt_name;
}

/**
 * Get path to directory for lambda specific csv files.
 * @return  Directory for the lambda CSVs
 */
filesystem::path ParameterCollection::lambda_dir() const {
    return this->base_path / this->lambda_dir_name;
}

/**
 * Get path to a lambda specific csv file given a lambda value.
 * @param[in]  lambda  lambda
 * @return     Path to the lambda CSV file
 */
filesystem::path ParameterCollection::lambda_csv(double lambda) const {
    return this->lambda_dir() / (this->lambda_csv_name_prefix + std::to_string(lambda) + ".csv");
}

/**
 * Get path to directory for graph specific csv files
 * @param[in]  lambda  lambda
 * @return             Path to the graphs directory
 */
filesystem::path ParameterCollection::graphs_dir(double lambda) const {
    return this->base_path / (this->graphs_dir_name_prefix + std::to_string(lambda));
}

/**
 * Get path to overview csv file with information theory values.
 * @return  Path to the info CSV file
 */
filesystem::path ParameterCollection::info_csv() const {
    return this->base_path / this->info_csv_name;
}

/**
 * Get whether status should be shown to stderr.
 * @return  Whether status should be shown.
 */
bool ParameterCollection::show_status() const {
    return this->config.misc.status_update_period >= 0;
}

/**
 * Get whether overall information theory data should be collected and the csv created.
 * @return  Whether info files should be generated
 */
bool ParameterCollection::generate_info() const {
    return this->config.misc.run_information_theoretical;
}

/**
 * Get whether lambda specific information should be collected and the csvs created.
 * @return  Whether lambda statistics should be generated
 */
bool ParameterCollection::generate_lambdas() const {
    return this->config.misc.run_inside_lambda;
}

/**
 * Get whether graph specific information should be collected and the csvs created.
 * @return  Whether graphs should be generated
 */
bool ParameterCollection::generate_graph() const {
    return this->config.misc.run_graph_visualization;
}

/**
 * Get minimum time in seconds to wait between refreshing status.
 * @return  Minimum time in seconds between status update
 */
double ParameterCollection::status_update_period() const {
    return std::max(this->config.misc.status_update_period, 0.0);
}

/**
 * Get total amount of runs (lambdas * realizations).
 * @return   Total amount of runs
 */
unsigned ParameterCollection::total_runs() const {
    return this->total_lambdas() * this->realizations();
}

/**
 * Get total amount of lambdas.
 * @return  Total amount of lambdas
 */
unsigned ParameterCollection::total_lambdas() const {
    switch (this->config.lambda.type) {
    case ParameterCollection::config::lambda::type::range:
    {
        double max = this->config.lambda.parameters.range.max;
        return this->elapsed_lambdas(max) + 1;
    }
    case ParameterCollection::config::lambda::type::custom:
        return static_cast<unsigned>(this->config.lambda.parameters.custom.count);
    default:
        return 0;
    }
}

/**
 * Get number of elapsed runs given current lambda.
 * @param[in]  lambda  lambda
 * @return     Elapsed runs
 */
unsigned ParameterCollection::elapsed_runs(double lambda) const {
    return this->elapsed_lambdas(lambda) * this->realizations();
}

/**
 * Get number of elapsed lambdas given current lambda.
 * @param[in]  lambda  lambda
 * @return     Elapsed lambdas
 */
unsigned ParameterCollection::elapsed_lambdas(double lambda) const {
    switch (this->config.lambda.type) {
    case ParameterCollection::config::lambda::type::range:
    {
        double min = this->config.lambda.parameters.range.min;
        double step = this->config.lambda.parameters.range.step;
        return static_cast<unsigned>(std::max(std::round((lambda - min)/step), 0.0));
    }
    case ParameterCollection::config::lambda::type::custom:
    {
        unsigned count = static_cast<unsigned>(this->config.lambda.parameters.custom.count);
        double *lambdas = this->config.lambda.parameters.custom.lambdas;
        double diff = std::numeric_limits<double>::max();
        unsigned current = 0;
        for (unsigned i=0; i<count; i++) {
            double d = std::abs(lambdas[i] - lambda);
            if (d < diff) {
                diff = d;
                current = i;
            }
        }
        return current;
    }
    default:
        return 0;
    }
}

/**
 * Get number of realizations per lambda.
 * @return  Realizations of each lambda
 */
unsigned ParameterCollection::realizations() const {
    return this->config.lambda.realizations;
}

/**
 * Get n size of matrix.
 * @return  Number of words
 */
unsigned ParameterCollection::n() const {
    return this->config.graph.n;
}

/**
 * Get m size of matrix.
 * @return  Number of meanings
 */
unsigned ParameterCollection::m() const {
    return this->config.graph.m;
}

/**
 * Construct and return an instance of LexicalMatrixOptimizer from
 * parameters.
 * @return  Optimizer
 */
LexicalMatrixOptimizer ParameterCollection::make_optimizer() const {
    OptimizerStopper stopCondition;
    switch (this->config.optimization.stop_condition.type) {
    case ParameterCollection::config::optimization::stop_condition::type::weak:
    {
        assert(this->config.optimization.stop_condition.parameters.weak.factor >= 0);
        unsigned factor = (unsigned)this->config.optimization.stop_condition.parameters.weak.factor;
        stopCondition = LexicalOptimizerCallbackFactory::StopCondition::weak(factor);
        break;
    }
    case ParameterCollection::config::optimization::stop_condition::type::strong:
    {
        assert(this->config.optimization.stop_condition.parameters.weak.factor >= 0);
        unsigned factor = (unsigned)this->config.optimization.stop_condition.parameters.strong.factor;
        stopCondition = LexicalOptimizerCallbackFactory::StopCondition::strong(factor);
        break;
    }
    case ParameterCollection::config::optimization::stop_condition::type::custom:
    {
        unsigned long trials = this->config.optimization.stop_condition.parameters.custom.trials;
        stopCondition = LexicalOptimizerCallbackFactory::StopCondition::custom(trials);
        break;
    }
    default:
        assert(false);
        break;
    }

    OptimizerInitializer graphInit;
    switch (this->config.initial_graph.type) {
    case ParameterCollection::config::initial_graph::type::random_edges:
    {
        unsigned edges = this->config.initial_graph.parameters.random_edges.edges;
        if (this->config.graph.unlinked_objects) {
            graphInit = LexicalOptimizerCallbackFactory::Initializer::random_edges(edges);
        } else {
            graphInit = LexicalOptimizerCallbackFactory::Initializer::random_edges_connected(edges);
        }
        break;
    }
    case ParameterCollection::config::initial_graph::type::random_probability:
    {
        double probability = this->config.initial_graph.parameters.random_probability.probability;
        if (this->config.graph.unlinked_objects) {
            graphInit = LexicalOptimizerCallbackFactory::Initializer::random_probability(probability);
        } else {
            graphInit = LexicalOptimizerCallbackFactory::Initializer::random_probability_connected(probability);
        }
        break;
    }
    case ParameterCollection::config::initial_graph::type::complete_graph:
        graphInit = LexicalOptimizerCallbackFactory::Initializer::complete_graph();
        break;
    case ParameterCollection::config::initial_graph::type::bijection:
        graphInit = LexicalOptimizerCallbackFactory::Initializer::bijection();
        break;
    default:
        assert(false);
        break;
    }

    OptimizerMutator getMutations;
    switch(this->config.optimization.mutations.type) {
    case ParameterCollection::config::optimization::mutations::type::constant:
    {
        unsigned count = this->config.optimization.mutations.parameters.constant.count;
        getMutations = LexicalOptimizerCallbackFactory::Mutator::constant(count);
        break;
    }
    case ParameterCollection::config::optimization::mutations::type::binomial:
    {
        double probability = this->config.optimization.mutations.parameters.binomial.probability;
        getMutations = LexicalOptimizerCallbackFactory::Mutator::binomial(probability);
        break;
    }
    default:
        assert(false);
        break;
    }
    
    return LexicalMatrixOptimizer(this->config.misc.seed, stopCondition, graphInit,
                                  getMutations, this->config.graph.unlinked_objects);
}

/**
 * Construct and return an instance of the correct kind of LexicalMatrix
 * from parameters.
 *
 * An instance of LexicalMatrix or DynamicLexicalMatrix is constrcuted
 * according to the parameteters within this class. The instance is kept in
 * the heap and a reference is returned. This is cleaned up when this
 * object is destroyed or when cleanup_matrix is called.
 * @return  Reference to a lexical matrix
 */
BaseLexicalMatrix &ParameterCollection::make_matrix() {
    assert(this->made_matrix == NULL);

    unsigned n = (unsigned)this->config.graph.n;
    unsigned m = (unsigned)this->config.graph.m;
    double phi = this->config.graph.phi;

    BaseLexicalMatrix *mat = NULL;
    if (this->config.graph.constant_prior_probability_of_meanings) {
        std::vector<long double> pi(m);
        this->assign_pi(pi);
        
        switch (this->config.optimization.computation_of_energy) {
        case ParameterCollection::config::optimization::computation_of_energy::static_:
            mat = new SecondModelStaticLexicalMatrix(n, m, pi, phi);
            break;
        case ParameterCollection::config::optimization::computation_of_energy::dynamic:
        default:
            mat = new SecondModelDynamicLexicalMatrix(n, m, pi, phi);
        }
    } else {
        switch (this->config.optimization.computation_of_energy) {
        case ParameterCollection::config::optimization::computation_of_energy::static_:
            mat = new FirstModelStaticLexicalMatrix(n, m, phi);
            break;
        case ParameterCollection::config::optimization::computation_of_energy::dynamic:
            mat = new FirstModelDynamicLexicalMatrix(n, m, phi);
            break;
        default:
            assert(false);
        }
    }
    this->made_matrix = mat;
    
    return *mat;
}

/**
 * Clean up the LexicalMatrix created by make_matrix, invalidating the
 * reference.
 */
void ParameterCollection::cleanup_matrix() {
    delete this->made_matrix;
    this->made_matrix = NULL;
}

/**
 * Print all parameters to an output stream.
 * @param[in,out]  out  Stream to write to
 * @param[in]      p    Parameters
 * @return              Stream written to
 */
std::ostream& operator<<(std::ostream& out, const ParameterCollection& p) {
    using namespace std;
    out << "graph:" << endl;
    out << "  n: " << p.config.graph.n << endl;
    out << "  m: " << p.config.graph.m << endl;
    out << "  phi: " << p.config.graph.phi << endl;
    out << "  unlinked_objects: ";
    if (p.config.graph.unlinked_objects) {
        out << "yes" << endl;
    } else {
        out << "no" << endl;
    }
    out << "  constant_prior_probability_of_meanings: ";
    if (p.config.graph.constant_prior_probability_of_meanings) {
        out << "yes" << endl;
    } else {
        out << "no" << endl;
    }
    if (p.config.graph.constant_prior_probability_of_meanings) {
        out << "  pi: " << endl;
        out << "    type: ";
        switch(p.config.graph.pi.type) {
        case ParameterCollection::config::graph::pi::type::uniform:
            out << "uniform" << endl;
            break;
        case ParameterCollection::config::graph::pi::type::brokenstick:
            out << "brokenstick" << endl;
            break;
        case ParameterCollection::config::graph::pi::type::geometric:
            out << "exponential" << endl;
            out << "  parameters:" << endl;
            out << "    p: " << p.config.graph.pi.parameters.geometric.p << endl;
            break;
        case ParameterCollection::config::graph::pi::type::powerlaw:
            out << "powerlaw" << endl;
            out << "  parameters:" << endl;
            out << "    alpha: " << p.config.graph.pi.parameters.powerlaw.alpha << endl;
            break;
        case ParameterCollection::config::graph::pi::type::custom:
            out << "custom" << endl;
            out << "  parameters:" << endl;
            out << "    probabilities: " << endl;
            for (unsigned i=0; i<p.config.graph.pi.parameters.custom.count; i++) {
                out << "      - " << p.config.graph.pi.parameters.custom.probabilities[i] << endl;
        }
            break;
        default:
            break;
        }
    }
    out << "lambda: " << endl;
    out << "  realizations: " << p.config.lambda.realizations << endl;
    out << "  type: ";
    switch(p.config.lambda.type) {
    case ParameterCollection::config::lambda::type::range:
        out << "range" << endl;
        out << "  parameters:" << endl;
        out << "    min: " << p.config.lambda.parameters.range.min << endl;
        out << "    max: " << p.config.lambda.parameters.range.max << endl;
        out << "    step: " << p.config.lambda.parameters.range.step << endl;
        break;
    case ParameterCollection::config::lambda::type::custom:
        out << "custom" << endl;
        out << "  parameters:" << endl;
        out << "    lambdas:" << endl;
        for (unsigned i=0; i<p.config.lambda.parameters.custom.count; i++) {
            out << "      - " << p.config.lambda.parameters.custom.lambdas[i] << endl;
        }
        break;
    default:
        out << "<?>" << endl;
    }
    out << "initial_graph:" << endl;
    out << "  type: ";
    switch(p.config.initial_graph.type) {
    case ParameterCollection::config::initial_graph::type::random_edges:
        out << "random_edges" << endl;
        out << "  parameters:" << endl;
        out << "    edges: " << p.config.initial_graph.parameters.random_edges.edges << endl;
        break;
    case ParameterCollection::config::initial_graph::type::random_probability:
        out << "random_probability" << endl;
        out << "  parameters:" << endl;
        out << "    probability: " << p.config.initial_graph.parameters.random_probability.probability << endl;
        break;
    case ParameterCollection::config::initial_graph::type::complete_graph:
        out << "complete_graph" << endl;
        break;
    case ParameterCollection::config::initial_graph::type::bijection:
        out << "bijection" << endl;
        break;
    default:
        out << "<?>" << endl;
    }
    out << "optimization:" << endl;
    out << "  computation_of_energy: ";
    switch (p.config.optimization.computation_of_energy) {
    case ParameterCollection::config::optimization::computation_of_energy::static_:
        out << "static" << endl;
        break;
    case ParameterCollection::config::optimization::computation_of_energy::dynamic:
        out << "dynamic" << endl;
        break;
    default:
        out << "<?>" << endl;
    }
    out << "  mutations:" << endl;
    out << "    type: ";
    switch (p.config.optimization.mutations.type) {
    case ParameterCollection::config::optimization::mutations::type::constant:
        out << "constant" << endl;
        out << "    parameters:" << endl;
        out << "      count: " << p.config.optimization.mutations.parameters.constant.count << endl;
        break;
    case ParameterCollection::config::optimization::mutations::type::binomial:
        out << "binomial" << endl;
        out << "    parameters:" << endl;
        out << "      probability: " << p.config.optimization.mutations.parameters.binomial.probability << endl;
        break;
    default:
        out << "<?>" << endl;
    }
    out << "  stop_condition:" << endl;
    out << "    type: ";
    switch (p.config.optimization.stop_condition.type) {
    case ParameterCollection::config::optimization::stop_condition::type::weak:
        out << "weak" << endl;
        out << "    parameters:" << endl;
        out << "      factor: " << p.config.optimization.stop_condition.parameters.weak.factor << endl;
        break;
    case ParameterCollection::config::optimization::stop_condition::type::strong:
        out << "strong" << endl;
        out << "    parameters:" << endl;
        out << "      factor: " << p.config.optimization.stop_condition.parameters.strong.factor << endl;
        break;
    case ParameterCollection::config::optimization::stop_condition::type::custom:
        out << "custom" << endl;
        out << "    parameters:" << endl;
        out << "      trials: " << p.config.optimization.stop_condition.parameters.custom.trials << endl;
        break;
    default:
        out << "<?>" << endl;
    }
    out << "misc:" << endl;
    out << "  seed: " << p.config.misc.seed << endl;
    out << "  status_update_period: " << p.config.misc.status_update_period << endl;
    out << "  run_information_theoretical: ";
    if (p.config.misc.run_information_theoretical) {
        out << "yes" << endl;
    } else {
        out << "no" << endl;
    }
    out << "  run_inside_lambda: ";
    if (p.config.misc.run_inside_lambda) {
        out << "yes" << endl;
    } else {
        out << "no" << endl;
    }
    out << "  run_graph_visualization: ";
    if (p.config.misc.run_graph_visualization) {
        out << "yes" << endl;
    } else {
        out << "no" << endl;
    }
    return out;
}
