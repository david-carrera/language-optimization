#include <Util.hpp>
#include <vector>

namespace ProbabilityDistributionGeneration {
    std::vector<long double> uniform(unsigned size);
    std::vector<long double> brokenstick(unsigned size);
    std::vector<long double> geometric(long double p, unsigned size);
    std::vector<long double> powerlaw(long double alpha, unsigned size);
    long double geometric_smallest_probability(long double p, unsigned size);
    long double powerlaw_smallest_probability(long double alpha, unsigned size);
};
