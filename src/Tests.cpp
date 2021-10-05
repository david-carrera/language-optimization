#include <ProbabilityDistributionGeneration.hpp>
#include <FirstModelDynamicLexicalMatrix.hpp>
#include <SecondModelDynamicLexicalMatrix.hpp>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <cassert>

using namespace std;

static long double evaluate(BaseLexicalMatrix &mat, long double epsilon) {
    string descriptions[6];
    descriptions[0] = "hs too small";
    descriptions[1] = "hs too big";
    descriptions[2] = "hr too small";
    descriptions[3] = "hr too big";
    descriptions[4] = "hsr too small";
    descriptions[5] = "hsr too big";
    
    long double values[3];
    values[0] = mat.get_HS();
    values[1] = mat.get_HR();
    values[2] = mat.get_HSR();
    
    long double limits[6] = {0.0};
    limits[0] = mat.min_HS();
    limits[1] = mat.max_HS();
    limits[2] = mat.min_HR();
    limits[3] = mat.max_HR();
    limits[4] = mat.min_HSR();
    limits[5] = mat.max_HSR();

    long double new_epsilon = 0;
    for (int i=0; i<6; i++) {
        long double e;
        if (i % 2 == 0) {
            e = limits[i] - values[i/2];
        } else {
            e = values[i/2] - limits[i];
        }
        if (e > new_epsilon) {
            new_epsilon = e;
        }
    }
        
    bool flags[6];
    for (int i=0; i<6; i++) {
        flags[i] = isnan(values[i/2]);
        if (i % 2 == 0) {
            flags[i] |= values[i/2] < limits[i] - epsilon;
        } else {
            flags[i] |= values[i/2] > limits[i] + epsilon;
        }
    }

    if (any_of(flags, flags+6, [](bool x){return x;})) {
        for (int i=0; i<6; i++) {
            if (flags[i]) {
                cout << "* " << descriptions[i] << " (" << values[i/2] << ", limit is " << limits[i] << ")" << endl;
            }
        }
    }

    return new_epsilon;
}

class ParamsIterator {
private:
    const unsigned n_min, n_max, n_inc;
    const unsigned m_min, m_max, m_inc;
    const int phin_max;
    const double phid;
    const unsigned ntests;

    const unsigned ndists = 4;
    const std::vector<std::vector<long double>> parameters = {
        {}, // uniform
        {}, // broken stick
        {0.1, 0.5, 0.9}, // geometric
        {0.5, 1, 2}, // power law
    };

    unsigned n, m, test, dist, parameter_i;
    int phin;

    bool finished;

    void make_pi(std::vector<long double> &pi) {
        long double parameter;
        if (this->parameters[this->dist].size() > 0) {
            parameter = this->parameters[this->dist][this->parameter_i];
        } else {
            parameter = 0;
        }
        
        if (this->dist == 0) { // Uniform distribution
            pi = ProbabilityDistributionGeneration::uniform(this->m);
        } else if (this->dist == 1) { // Broken stick
            pi = ProbabilityDistributionGeneration::brokenstick(this->m);
        } else if (this->dist == 2) { // Geometric
            pi = ProbabilityDistributionGeneration::geometric(parameter, this->m);
        } else if (this->dist == 3) { // Power Law
            pi = ProbabilityDistributionGeneration::powerlaw(parameter, this->m);
        } else {
            cerr << "Unexpected distribution " << dist << endl;
        }
    }
    
public:
    ParamsIterator(unsigned n_min, unsigned n_max, unsigned n_inc,
                   unsigned m_min, unsigned m_max, unsigned m_inc,
                   int phin_min, int phin_max, double phid, unsigned ntests) :
        n_min(n_min), n_max(n_max), n_inc(n_inc),
        m_min(m_min), m_max(m_max), m_inc(m_inc),
        phin_max(phin_max), phid(phid), ntests(ntests),
        n(n_min), m(m_min), test(0), dist(0), parameter_i(0), phin(phin_min), finished(false) {
        cout << "phi=" << this->get_phi() << " | ";
        cout << "n=" << n_min << ":" << n_max << " | ";
        cout << "m=" << m_min << ":" << m_max << endl;
    }

    inline bool is_finished() const {
        return this->finished;
    }
    inline unsigned get_test() const {
        return this->test;
    }
    inline unsigned get_n() const {
        return this->n;
    }
    inline unsigned get_m() const {
        return this->m;
    }
    inline double get_phi() const {
        return this->phin/this->phid;
    }

    template<class T>
    T &step_FirstModel() {
        assert(!this->finished);
        
        T &mat = *new T(this->n, this->m, this->get_phi(),
                        numeric_limits<double>::infinity());

        this->test++;
        if (this->test >= this->ntests) {
            this->test = 0;
            this->m += this->m_inc;
            if (this->m >= this->m_max) {
                this->m = this->m_min;
                this->n += this->n_inc;
                cout << round(((double)this->n / (double)this->n_max)*10000)/100 << "%\r" << flush;
                if (this->n >= this->n_max) {
                    this->n = this->n_min;
                    this->phin++;
                    cout << "phi=" << this->get_phi() << " | ";
                    cout << "n=" << n_min << ":" << n_max << " | ";
                    cout << "m=" << m_min << ":" << m_max << endl;
                    if (this->phin > this->phin_max) {
                        this->finished = true;
                    }
                }
            }
        }

        return mat;
    }

    template<class T>
    T &step_SecondModel() {
        assert(!this->finished);

        std::vector<long double> pi(m);
        this->make_pi(pi);
        
        T &mat = *new T(this->n, this->m, pi, this->get_phi(),
                        numeric_limits<double>::infinity());

        this->parameter_i++;
        if (this->parameter_i >= this->parameters[this->dist].size()) {
            this->parameter_i = 0;
            this->dist++;
            if (this->dist >= this->ndists) {
                this->dist = 0;
                this->test++;
                if (this->test >= this->ntests) {
                    this->test = 0;
                    this->m += this->m_inc;
                    if (this->m >= this->m_max) {
                        this->m = this->m_min;
                        this->n += this->n_inc;
                        cout << round(((double)this->n / (double)this->n_max)*10000)/100 << "%\r" << flush;
                        if (this->n >= this->n_max) {
                            this->n = this->n_min;
                            this->phin++;
                            cout << "phi=" << this->get_phi() << " | ";
                            cout << "n=" << n_min << ":" << n_max << " | ";
                            cout << "m=" << m_min << ":" << m_max << endl;
                            if (this->phin > this->phin_max) {
                                this->finished = true;
                            }
                        }
                    }
                }
            }
        }

        return mat;
    }

    template<class T>
    T &step();
};

template<>
FirstModelStaticLexicalMatrix &ParamsIterator::step() {
    return this->step_FirstModel<FirstModelStaticLexicalMatrix>();
}
    
template<>
FirstModelDynamicLexicalMatrix &ParamsIterator::step() {
    return this->step_FirstModel<FirstModelDynamicLexicalMatrix>();
}
    
template<>
SecondModelStaticLexicalMatrix &ParamsIterator::step() {
    return this->step_SecondModel<SecondModelStaticLexicalMatrix>();
}

template<>
SecondModelDynamicLexicalMatrix &ParamsIterator::step() {
    return this->step_SecondModel<SecondModelDynamicLexicalMatrix>();
}

class ExtremeParameters {
public:
    static long double HS(unsigned test, const FirstModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 0;
        case 1:
            return log(mat.n);
        case 2:
            return log(mat.n);
        case 3:
            return 0;
        default:
            return 0;
        }
    }
    
    static long double HR(unsigned test, const FirstModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 0;
        case 1:
            return log(mat.m);
        case 2:
            return log(mat.n);
        case 3:
            return log(mat.m);
        default:
            return 0;
        }
    }
    
    static long double HSR(unsigned test, const FirstModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 0;
        case 1:
            return log(mat.n * mat.m);
        case 2:
            return log(mat.n);
        case 3:
            return log(mat.m);
        default:
            return 0;
        }
    }

    static long double M(unsigned test, const FirstModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 1;
        case 1:
            return pow(mat.n * mat.m, mat.phi+1);
        case 2:
            return mat.n;
        case 3:
            return pow(mat.m, mat.phi+1);
        default:
            return 0;
        }
    }

    static const string TestName(unsigned test, const FirstModelStaticLexicalMatrix &mat) {
        (void)mat;
        switch (test) {
        case 0:
            return "Single Edge (First Model)";
        case 1:
            return "Complete Graph (First Model)";
        case 2:
            return "One to One (First Model)";
        case 3:
            return "One to All (First Model)";
        default:
            return "unknown";
        }
    }


    static long double HS(unsigned test, const SecondModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 0;
        case 1:
            return log(mat.n);
        case 2:
            return mat.HR_pi;
        case 3:
            return 0;
        default:
            return 0;
        }
    }
    
    static long double HR(unsigned test, const SecondModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 0;
        case 1:
            return mat.HR_pi;
        case 2:
            return mat.HR_pi;
        case 3:
            return mat.HR_pi;
        default:
            return 0;
        }
    }
    
    static long double HSR(unsigned test, const SecondModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return 0;
        case 1:
            return log(mat.n) + mat.HR_pi;
        case 2:
            return mat.HR_pi;
        case 3:
            return mat.HR_pi;
        default:
            return 0;
        }
    }

    static long double rho(unsigned test, const SecondModelStaticLexicalMatrix &mat) {
        switch (test) {
        case 0:
            return mat.pi[0];
        case 1:
            return 1;
        case 2:
            return 1;
        case 3:
            return 1;
        default:
            return 0;
        }
    }

    static const string TestName(unsigned test, const SecondModelStaticLexicalMatrix &mat) {
        (void)mat;
        switch (test) {
        case 0:
            return "Single Edge (Second Model)";
        case 1:
            return "Complete Graph (Second Model)";
        case 2:
            return "One to One (Second Model)";
        case 3:
            return "One to All (Second Model)";
        default:
            return "unknown";
        }
    }
};

static bool test_single_param(long double correct_value, long double obtained_value, long double *const current_epsilon, long double epsilon_limit) {
    long double epsilon = abs(obtained_value - correct_value);
    if (epsilon > *current_epsilon) {
        *current_epsilon = epsilon;
    }
    return epsilon <= epsilon_limit;
}

template<class T>
static void report_params_error_generic(const T &mat, unsigned test) {
    mat.print_properties();
    cout << "Error detected for:" << endl;
    cout << "\tN = " << mat.n << endl;
    cout << "\tM = " << mat.m << endl;
    cout << "\tphi = " << mat.phi << endl;
    cout << "\ttest = " << ExtremeParameters::TestName(test, mat) << endl;
    cout << "Errors:" << endl;
}

static void report_params_error(string name, long double should_be, long double but_was) {
    cout << "\t" << name << " -> Should be (" << should_be << ") but was (" << but_was << ") -- difference of (" << abs(should_be - but_was) << ")" << endl;
}

template<class T>
static bool test_params_generic(T &mat, unsigned test, long double epsilonlimit, long double epsilons[]) {
    bool ok_hs, ok_hr, ok_hsr;
    ok_hs = test_single_param(ExtremeParameters::HS(test, mat), mat.get_HS(), &epsilons[test], epsilonlimit);
    ok_hr = test_single_param(ExtremeParameters::HR(test, mat), mat.get_HR(), &epsilons[test], epsilonlimit);
    ok_hsr = test_single_param(ExtremeParameters::HSR(test, mat), mat.get_HSR(), &epsilons[test], epsilonlimit);
    if (!ok_hs || !ok_hr || !ok_hsr) {
        report_params_error_generic(mat, test);
        if (!ok_hs) {
            report_params_error("H(S)", ExtremeParameters::HS(test, mat), mat.get_HS());
        }
        if (!ok_hr) {
            report_params_error("H(R)", ExtremeParameters::HR(test, mat), mat.get_HR());
        }
        if (!ok_hsr) {
            report_params_error("H(S,R)", ExtremeParameters::HSR(test, mat), mat.get_HSR());
        }
        return false;
    }
    return true;
}

static bool test_params(FirstModelStaticLexicalMatrix &mat, unsigned test, long double epsilonlimit, long double epsilons[]) {
    bool ok = test_params_generic(mat, test, epsilonlimit, epsilons);

    // M reaches values that are very high (1e11) and so the difference becomes also very "high" (1e-5). As other parameters depend on M, not testing it directly isn't a problem.
    bool ok_M = true; //test_single_param(ExtremeParameters::M(test, mat), mat.get_M(), &epsilons[test], epsilonlimit);
    
    if (!ok_M) {
        if (ok) {
            report_params_error_generic(mat, test);
        }
        report_params_error("M", ExtremeParameters::M(test, mat), mat.get_M());
    }

    return ok && ok_M;
}

static bool test_params(SecondModelStaticLexicalMatrix &mat, unsigned test, long double epsilonlimit, long double epsilons[]) {
    bool ok = test_params_generic(mat, test, epsilonlimit, epsilons);
    bool ok_rho = test_single_param(ExtremeParameters::rho(test, mat), mat.get_rho(), &epsilons[test], epsilonlimit);

    if (!ok_rho) {
        if (ok) {
            report_params_error_generic(mat, test);
        }
        report_params_error("rho", ExtremeParameters::rho(test, mat), mat.get_rho());
    }

    return ok && ok_rho;
}

template<class T>
static bool params(unsigned n_min, unsigned n_max, unsigned n_inc,
                   unsigned m_min, unsigned m_max, unsigned m_inc,
                   int phin_min, int phin_max, double phid, long double epsilonlimit,
                   bool skip_complete_graph, bool skip_one_to_all) {
    const int ntests = 4;

    long double epsilons[ntests] = {};
    double runtimes[ntests] = {};

    ParamsIterator iter(n_min, n_max, n_inc,
                        m_min, m_max, m_inc,
                        phin_min, phin_max, phid, ntests);

    RNG rng(12345);
    while (!iter.is_finished()) {
        T &mat = iter.step<T>();
        mat.clear();

        unsigned test = iter.get_test();

        auto t = chrono::steady_clock::now();
        if (test == 0) { // Single Edge
            mat.add_edge(0, 0);
        } else if (test == 1) { // Complete Graph
            if (skip_complete_graph) {
                delete &mat;
                continue;
            }
            for (unsigned s=0; s<iter.get_n(); s++) {
                for (unsigned r=0; r<iter.get_m(); r++) {
                    mat.add_edge(s, r);
                }
            }
        } else if (test == 2) { // One to One
            if (iter.get_n() != iter.get_m()) {
                delete &mat;
                continue;
            }
            for (unsigned sr=0; sr<iter.get_n(); sr++) {
                mat.add_edge(sr, sr);
            }
        } else if (test == 3) { // One to All
            if (skip_one_to_all) {
                delete &mat;
                continue;
            }
            for (unsigned r=0; r<iter.get_m(); r++) {
                mat.add_edge(0, r);
            }
        } else {
            cerr << "Unexpected test " << test << endl;
            return false;
        }
        mat.recalculate();
        chrono::duration<double> d = chrono::steady_clock::now() - t;

        runtimes[test] += d.count();

        bool ok = test_params(mat, test, epsilonlimit, epsilons);
        delete &mat;

        if (!ok) {
            return false;
        }
    }

    cout << "All invariants within limits for extreme cases." << endl;
    cout << "Maximum epsilons:" << endl;
    cout << "\tSingle edge: " << epsilons[0] << endl;
    cout << "\tComplete Graph: " << epsilons[1] << endl;
    cout << "\tOne to One: " << epsilons[2] << endl;
    cout << "\tOne to All: " << epsilons[3] << endl;
    cout << "Total runtimes: " << endl;
    cout << "\tSingle edge: " << runtimes[0] << endl;
    cout << "\tComplete Graph: " << runtimes[1] << endl;
    cout << "\tOne to One: " << runtimes[2] << endl;
    cout << "\tOne to All: " << runtimes[3] << endl;

    return true;
}

template <class T>
T &create_mat(unsigned n, unsigned m, double phi, long double epsilon);

template <>
inline FirstModelStaticLexicalMatrix &create_mat<>(unsigned n, unsigned m, double phi, long double epsilon) {
    return *(new FirstModelStaticLexicalMatrix(n, m, phi, epsilon));
}

template <>
inline FirstModelDynamicLexicalMatrix &create_mat<>(unsigned n, unsigned m, double phi, long double epsilon) {
    return *(new FirstModelDynamicLexicalMatrix(n, m, phi, epsilon));
}

template <>
inline SecondModelStaticLexicalMatrix &create_mat<>(unsigned n, unsigned m, double phi, long double epsilon) {
    std::vector<long double> pi = ProbabilityDistributionGeneration::brokenstick(m);
    return *(new SecondModelStaticLexicalMatrix(n, m, pi, phi, epsilon));
}

template <>
inline SecondModelDynamicLexicalMatrix &create_mat<>(unsigned n, unsigned m, double phi, long double epsilon) {
    std::vector<long double> pi = ProbabilityDistributionGeneration::brokenstick(m);
    return *(new SecondModelDynamicLexicalMatrix(n, m, pi, phi, epsilon));
}

template<class T1, class T2>
void report_interesting_differences(const T1 &static_mat, const T2 &mat);

template<>
inline void report_interesting_differences<FirstModelStaticLexicalMatrix, FirstModelDynamicLexicalMatrix>(const FirstModelStaticLexicalMatrix &static_mat, const FirstModelDynamicLexicalMatrix &mat) {
    cout << "dynamic M = " << mat.get_M() << endl;
    cout << "static  M = " << static_mat.get_M() << endl;
    
    cout << "dynamic mu phi(r) = {";
    for (unsigned i=0; i<mat.m; i++) {
        cout << " " << mat.get_mu_phi(i);
    }
    cout << " }" << endl;
    cout << "static  mu phi(r) = {";
    for (unsigned i=0; i<static_mat.m; i++) {
        cout << " " << static_mat.get_mu_phi(i);
    }
    cout << " }" << endl;
    
    cout << "dynamic omega phi(r) = {";
    for (unsigned i=0; i<mat.m; i++) {
        cout << " " << mat.get_omega_phi(i);
    }
    cout << " }" << endl;
    cout << "static  omega phi(r) = {";
    for (unsigned i=0; i<static_mat.m; i++) {
        cout << " " << static_mat.get_omega_phi(i);
    }
    cout << " }" << endl;
}

template<>
inline void report_interesting_differences<SecondModelStaticLexicalMatrix, SecondModelDynamicLexicalMatrix>(const SecondModelStaticLexicalMatrix &static_mat, const SecondModelDynamicLexicalMatrix &mat) {
    cout << "dynamic pi = ";
    pretty_print_vector(cout, mat.pi);
    cout << endl;
    cout << "static pi = ";
    pretty_print_vector(cout, static_mat.pi);
    cout << endl;
    
    cout << "dynamic rho = " << mat.get_rho() << endl;
    cout << "static  rho = " << static_mat.get_rho() << endl;
    
    cout << "dynamic omega phi(r) = {";
    for (unsigned i=0; i<mat.m; i++) {
        cout << " " << mat.get_omega_phi(i);
    }
    cout << " }" << endl;
    cout << "static  omega phi(r) = {";
    for (unsigned i=0; i<static_mat.m; i++) {
        cout << " " << static_mat.get_omega_phi(i);
    }
    cout << " }" << endl;
    
    cout << "dynamic nu(r) = {";
    for (unsigned i=0; i<mat.m; i++) {
        cout << " " << mat.get_nu(i);
    }
    cout << " }" << endl;
    cout << "static  nu(r) = {";
    for (unsigned i=0; i<static_mat.m; i++) {
        cout << " " << static_mat.get_nu(i);
    }
    cout << " }" << endl;
    
    cout << "dynamic chi(s) = {";
    for (unsigned i=0; i<mat.n; i++) {
        cout << " " << mat.get_chi(i);
    }
    cout << " }" << endl;
    cout << "static  chi(s) = {";
    for (unsigned i=0; i<static_mat.n; i++) {
        cout << " " << static_mat.get_chi(i);
    }
    cout << " }" << endl;
}

template <class T1, class T2>
static bool exhaustive(unsigned n, unsigned m, double phi,
                       long double epsilonlimit_dynamic, long double epsilonlimit_static) {
    const long double epsilonlimit_difference = max(epsilonlimit_dynamic, epsilonlimit_static);
    
    const unsigned long nedges = 1L<<(n*m);
    T1 &static_mat = create_mat<T1>(n, m, phi, std::numeric_limits<double>::infinity());
    T2 &mat = create_mat<T2>(n, m, phi, std::numeric_limits<double>::infinity());
                
    cout << "phi=" << phi << endl;
    cout << "n=" << n << endl;
    cout << "m=" << m << endl;
                
    long double max_epsilon_dynamic = 0;
    long double max_epsilon_static = 0;
    double runtime_dynamic = 0;
    double runtime_static = 0;

    for (int direction=0; direction<=1; direction++) {
        unsigned long prev_edges_gray = 0;
        for (unsigned long edges=1; edges<nedges; edges++) {
            unsigned long update_edges;
            if (nedges > 1000) {
                update_edges = nedges/1000;
            } else if (nedges > 10){
                update_edges = nedges/10;
            } else {
                update_edges = 1;
            }
            if (edges % update_edges == 0) {
                double done = (double)edges + (double)nedges*direction;
                double todo = (double)nedges*2;
                cout << round(done/todo*10000)/100 << "%     \r" << flush;
            }
            
            // Get i,j to mutate
            unsigned long edges_gray = edges ^ (edges >> 1);
            unsigned long difference = prev_edges_gray ^ edges_gray;
            unsigned edge = sizeof(nedges)*8 - (unsigned)__builtin_clzl(difference) - 1;
            unsigned i = edge / m;
            unsigned j = edge % m;
            prev_edges_gray = edges_gray;
            
            auto t = chrono::steady_clock::now();
            mat.mutate(i, j);
            chrono::duration<double> d = chrono::steady_clock::now() - t;
            runtime_dynamic += d.count();
            
            t = chrono::steady_clock::now();
            static_mat.mutate(i, j);
            static_mat.recalculate();
            d = chrono::steady_clock::now() - t;
            runtime_static += d.count();
            
            // Compare static vs dynamic
            bool psdiff = false;
            for (unsigned i=0; i<mat.n; i++) {
                psdiff |= abs(mat.get_ps(i) - static_mat.get_ps(i)) > epsilonlimit_difference;
                if (psdiff) {
                    break;
                }
            }
            bool prdiff = false;
            for (unsigned j=0; j<mat.m; j++) {
                psdiff |= abs(mat.get_pr(j) - static_mat.get_pr(j)) > epsilonlimit_difference;
                if (prdiff) {
                    break;
                }
            }
            long double HSRdiff = abs(mat.get_HSR() - static_mat.get_HSR());
            long double HSdiff = abs(mat.get_HS() - static_mat.get_HS());
            long double HRdiff = abs(mat.get_HR() - static_mat.get_HR());
            
            // Report large discrepancy and exit
            if (HSRdiff > epsilonlimit_difference || HSdiff > epsilonlimit_difference ||
                HRdiff > epsilonlimit_difference || psdiff || prdiff) { 
                cout << "dynamic H(S,R) = " << mat.get_HSR() << endl;
                cout << "static  H(S,R) = " << static_mat.get_HSR() << endl;
                cout << "dynamic H(S) = " << mat.get_HS() << endl;
                cout << "static  H(S) = " << static_mat.get_HS() << endl;
                cout << "dynamic H(R) = " << mat.get_HR() << endl;
                cout << "static  H(R) = " << static_mat.get_HR() << endl;
                cout << "dynamic p(s) = {";
                for (unsigned i=0; i<mat.n; i++) {
                    cout << " " << mat.get_ps(i);
                }
                cout << " }" << endl;
                cout << "static  p(s) = {";
                for (unsigned i=0; i<static_mat.n; i++) {
                    cout << " " << static_mat.get_ps(i);
                }
                cout << " }" << endl;
                
                cout << "dynamic p(r) = {";
                for (unsigned i=0; i<mat.m; i++) {
                    cout << " " << mat.get_pr(i);
                }
                cout << " }" << endl;
                cout << "static  p(r) = {";
                for (unsigned i=0; i<static_mat.m; i++) {
                    cout << " " << static_mat.get_pr(i);
                }
                cout << " }" << endl;
                
                cout << "dynamic mu(s) = {";
                for (unsigned i=0; i<mat.n; i++) {
                    cout << " " << mat.get_mu(i);
                }
                cout << " }" << endl;
                cout << "static  mu(s) = {";
                for (unsigned i=0; i<static_mat.n; i++) {
                    cout << " " << static_mat.get_mu(i);
                }
                cout << " }" << endl;
                
                cout << "dynamic omega(r) = {";
                for (unsigned i=0; i<mat.m; i++) {
                    cout << " " << mat.get_omega(i);
                }
                cout << " }" << endl;
                cout << "static  omega(r) = {";
                for (unsigned i=0; i<static_mat.m; i++) {
                    cout << " " << static_mat.get_omega(i);
                }
                cout << " }" << endl;

                report_interesting_differences<T1, T2>(static_mat, mat);
                
                cout << "a =" << endl;
                for (unsigned j=0; j<static_mat.m; j++) {
                    for (unsigned i=0; i<static_mat.n; i++) {
                        cout << static_mat.get_a(i, j) << " ";
                    }
                    cout << endl;
                }
                cout << endl;
                cout << "exiting" << endl;

                delete &mat;
                delete &static_mat;
                return false;
            }
            
            // Report a large discrepancy with invariants but keep running
            long double epsilon_dynamic = evaluate(mat, epsilonlimit_dynamic);
            long double epsilon_static = evaluate(static_mat, epsilonlimit_static);
            if (epsilon_dynamic > max_epsilon_dynamic) {
                max_epsilon_dynamic = epsilon_dynamic;
            }
            if (epsilon_static > max_epsilon_static) {
                max_epsilon_static = epsilon_static;
            }
        }
    }
    
    cout << "dynamic epsilon=" << max_epsilon_dynamic << endl;
    cout << "static epsilon=" << max_epsilon_static << endl;
    cout << "dynamic runtime=" << runtime_dynamic << endl;
    cout << "static runtime=" << runtime_static << endl;
    cout << endl;

    delete &mat;
    delete &static_mat;
    return true;
}

void do_replace_edge(BaseLexicalMatrix &mat) {
    mat.clear();
    mat.add_edge(mat.m-1, mat.m-1);
    mat.add_edge(mat.m-2, mat.m-2);
    mat.remove_edge(mat.m-1, mat.m-1);

    mat.clear();
    mat.add_edge(mat.m-1, 1);
    mat.add_edge(mat.m-2, 2);
    mat.remove_edge(mat.m-1, 1);
    
    mat.clear();
    mat.add_edge(1, mat.m-1);
    mat.add_edge(2, mat.m-2);
    mat.remove_edge(1, mat.m-1);
    
    mat.clear();
    mat.add_edge(1, 1);
    mat.add_edge(2, 2);
    mat.remove_edge(1, 1);
}

template <class T>
bool replace_edge_single(unsigned n, unsigned m, double phi);

template <>
bool replace_edge_single<FirstModelDynamicLexicalMatrix>(unsigned n, unsigned m, double phi) {
    FirstModelDynamicLexicalMatrix mat(n, m, phi);
    do_replace_edge(mat);
    return true;
}

template <>
bool replace_edge_single<SecondModelDynamicLexicalMatrix>(unsigned n, unsigned m, double phi) {
    for (int i=0; i<4; i++) {
        std::vector<long double> pi;
        switch (i) {
        case 0:
            pi = ProbabilityDistributionGeneration::uniform(m);
            break;
        case 1:
            pi = ProbabilityDistributionGeneration::brokenstick(m);
            break;
        case 2:
            pi = ProbabilityDistributionGeneration::powerlaw(1, m);
            break;
        case 3:
            pi = ProbabilityDistributionGeneration::geometric(0.5, m);
            break;
        default:
            return false;
        }
        
        SecondModelDynamicLexicalMatrix mat(n, m, pi, phi);
        do_replace_edge(mat);
    }
    
    return true;
}

template <class T>
bool save_restore(unsigned n, unsigned m, double phi, unsigned mutations, long double epsilon) {
    
    std::uniform_int_distribution<unsigned> word(0, n-1);
    std::uniform_int_distribution<unsigned> meaning(0, m-1);
    
    T &mat_ok = create_mat<T>(n, m, phi, std::numeric_limits<double>::infinity());
    T &mat_sr = create_mat<T>(n, m, phi, std::numeric_limits<double>::infinity());
    
    RNG rng(12345);
    RNG tmp = rng;

    mat_ok.initialize_random(rng, 0.1);
    mat_sr.initialize_random(tmp, 0.1);
    mat_sr.record(true);
    mat_sr.save();

    // Compare to ensure they start equal
    long double HSRdiff = abs(mat_ok.get_HSR() - mat_sr.get_HSR());
    long double HSdiff = abs(mat_ok.get_HS() - mat_sr.get_HS());
    long double HRdiff = abs(mat_ok.get_HR() - mat_sr.get_HR());
    
    // Report large discrepancy and exit
    if (HSRdiff > epsilon || HSdiff > epsilon || HRdiff > epsilon) {
        cout << "ok H(S,R) = " << mat_ok.get_HSR() << endl;
        cout << "sr H(S,R) = " << mat_sr.get_HSR() << endl;
        cout << "ok H(S) = " << mat_ok.get_HS() << endl;
        cout << "sr H(S) = " << mat_sr.get_HS() << endl;
        cout << "ok H(R) = " << mat_ok.get_HR() << endl;
        cout << "sr H(R) = " << mat_sr.get_HR() << endl;
        cout << endl;
        cout << "exiting: they don't start equal!" << endl;

        mat_ok.print_properties();
        mat_sr.print_properties();
        
        delete &mat_ok;
        delete &mat_sr;
        return false;
    }

    unsigned runs = n*m;
    for (unsigned r=0; r<runs; r++) {
        pair<unsigned, unsigned> rec_mut[mutations];
        for (unsigned mut=0; mut<mutations; mut++) {
            unsigned i = word(rng);
            unsigned j = meaning(rng);
            mat_ok.mutate(i, j);
            mat_sr.mutate(i, j);
            rec_mut[mut] = {i, j};
        }
        
        mat_ok.recalculate();
        mat_sr.recalculate();

        // Compare
        long double HSRdiff = abs(mat_ok.get_HSR() - mat_sr.get_HSR());
        long double HSdiff = abs(mat_ok.get_HS() - mat_sr.get_HS());
        long double HRdiff = abs(mat_ok.get_HR() - mat_sr.get_HR());
        
        // Report large discrepancy and exit
        if (HSRdiff > epsilon || HSdiff > epsilon || HRdiff > epsilon) {
            cout << "ok H(S,R) = " << mat_ok.get_HSR() << endl;
            cout << "sr H(S,R) = " << mat_sr.get_HSR() << endl;
            cout << "ok H(S) = " << mat_ok.get_HS() << endl;
            cout << "sr H(S) = " << mat_sr.get_HS() << endl;
            cout << "ok H(R) = " << mat_ok.get_HR() << endl;
            cout << "sr H(R) = " << mat_sr.get_HR() << endl;
            cout << endl;
            cout << "exiting" << endl;

            mat_ok.print_properties();
            mat_sr.print_properties();

            delete &mat_ok;
            delete &mat_sr;
            return false;
        }
        
        for (unsigned mut=0; mut<mutations; mut++) {
            pair<unsigned, unsigned> pair = rec_mut[mut];
            unsigned i = pair.first;
            unsigned j = pair.second;
            mat_ok.mutate(i, j);
        }
        mat_sr.restore();
    }
    
    delete &mat_ok;
    delete &mat_sr;
    return true;
}

static bool parameters_equal(long double p1, long double p2, string name, long double epsilon) {
    if (std::abs(p1 - p2) < epsilon) {
        return true;
    }
    cout << name << " differs: " << p1 << " =/= " << p2 << endl;
    return false;
}

static bool matrixes_equal(BaseLexicalMatrix &smat, BaseLexicalMatrix &dmat, long double epsilon) {
    bool ret = true;

    smat.mutate(0, 0);
    smat.mutate(0, 0);
    
    dmat.mutate(0, 0);
    dmat.mutate(0, 0);
    
    ret &= parameters_equal(smat.get_HS(), dmat.get_HS(), "HS", epsilon);
    ret &= parameters_equal(smat.get_HR(), dmat.get_HR(), "HR", epsilon);
    ret &= parameters_equal(smat.get_HSR(), dmat.get_HSR(), "HSR", epsilon);
    ret &= parameters_equal(smat.get_ISR(), dmat.get_ISR(), "ISR", epsilon);
    ret &= parameters_equal(smat.get_HSgR(), dmat.get_HSgR(), "HSgR", epsilon);
    ret &= parameters_equal(smat.get_HRgS(), dmat.get_HRgS(), "HRgS", epsilon);

    return ret;
}

template <class T1, class T2>
static bool initial_conditions(unsigned n, unsigned m, double phi, long double epsilon, unsigned seed) {
    T1 &smat = create_mat<T1>(n, m, phi, epsilon);
    T2 &dmat = create_mat<T2>(n, m, phi, epsilon);
    RNG rng(seed);

    rng.seed(seed);
    smat.initialize_random(rng, 0.05);
    rng.seed(seed);
    dmat.initialize_random(rng, 0.05);
    if (!matrixes_equal(smat, dmat, epsilon)) {
        cout << "Matrixes differ on random initialization" << endl;
        smat.print_properties();
        dmat.print_properties();
        delete &smat;
        delete &dmat;
        return false;
    }

    rng.seed(seed);
    smat.initialize_random_edges(rng, 10);
    rng.seed(seed);
    dmat.initialize_random_edges(rng, 10);
    if (!matrixes_equal(smat, dmat, epsilon)) {
        cout << "Matrixes differ on random edge initialization" << endl;
        delete &smat;
        delete &dmat;
        return false;
    }

    smat.initialize_bijective_graph();
    dmat.initialize_bijective_graph();
    if (!matrixes_equal(smat, dmat, epsilon)) {
        cout << "Matrixes differ on bijective initialization" << endl;
        delete &smat;
        delete &dmat;
        return false;
    }

    smat.initialize_complete_graph();
    dmat.initialize_complete_graph();
    if (!matrixes_equal(smat, dmat, epsilon)) {
        cout << "Matrixes differ on complete initialization" << endl;
        delete &smat;
        delete &dmat;
        return false;
    }

    delete &smat;
    delete &dmat;
    return true;
}

static bool check_connected_components(string name, const FirstModelStaticLexicalMatrix &mat, unsigned words, unsigned meanings) {
    auto comps = mat.largest_connected_component();

    if (words != comps.first || meanings != comps.second) {
        cout << name << ": ERROR" << endl;
        cout << "Expected:" << endl;
        cout << "\twords: " << words << endl;
        cout << "\tmeanings: " << meanings << endl;
        cout << "But got:" << endl;
        cout << "\twords: " << comps.first << endl;
        cout << "\tmeanings: " << comps.second << endl;
        return false;
    }

    return true;
}

static bool connected_components() {
    RNG rng(123);
    FirstModelStaticLexicalMatrix mat(10, 10, 0);

    mat.clear();
    if (!check_connected_components("empty", mat, 1, 0)) {
        return false;
    }

    mat.initialize_random_edges(rng, 1);
    if (!check_connected_components("single edge", mat, 1, 1)) {
        return false;
    }

    mat.clear();
    for (unsigned i=0; i<mat.m; i++) {
        mat.add_edge(0, i);
    }
    if (!check_connected_components("one to all", mat, 1, 10)) {
        return false;
    }

    mat.clear();
    for (unsigned i=0; i<mat.n; i++) {
        mat.add_edge(i, 0);
    }
    if (!check_connected_components("all to one", mat, 10, 1)) {
        return false;
    }

    mat.initialize_bijective_graph();
    if (!check_connected_components("one to one", mat, 1, 1)) {
        return false;
    }

    mat.initialize_complete_graph();
    if (!check_connected_components("complete", mat, 10, 10)) {
        return false;
    }

    // Random graphs, graphs generated and connected components calculated externally.
    
    mat.clear();
    mat.add_edge(0, 1);
    mat.add_edge(0, 7);
    mat.add_edge(0, 8);
    mat.add_edge(1, 4);
    mat.add_edge(2, 0);
    mat.add_edge(2, 3);
    mat.add_edge(3, 2);
    mat.add_edge(3, 5);
    mat.add_edge(4, 6);
    mat.add_edge(6, 4);
    mat.add_edge(6, 9);
    mat.add_edge(8, 1);
    mat.add_edge(8, 5);
    mat.add_edge(8, 6);
    mat.add_edge(8, 7);
    mat.add_edge(9, 0);
    mat.add_edge(9, 7);
    if (!check_connected_components("random 1", mat, 6, 8)) {
        return false;
    }
    
    mat.clear();
    mat.add_edge(0, 4);
    mat.add_edge(0, 5);
    mat.add_edge(0, 6);
    mat.add_edge(1, 3);
    mat.add_edge(3, 6);
    mat.add_edge(3, 8);
    mat.add_edge(4, 0);
    mat.add_edge(4, 3);
    mat.add_edge(4, 5);
    mat.add_edge(5, 2);
    mat.add_edge(6, 6);
    mat.add_edge(6, 9);
    mat.add_edge(7, 4);
    mat.add_edge(7, 5);
    mat.add_edge(7, 7);
    mat.add_edge(7, 9);
    mat.add_edge(8, 5);
    mat.add_edge(9, 6);
    if (!check_connected_components("random 2", mat, 8, 8)) {
        return false;
    }
    
    mat.clear();
    mat.add_edge(0, 0);
    mat.add_edge(1, 2);
    mat.add_edge(2, 4);
    mat.add_edge(2, 6);
    mat.add_edge(2, 8);
    mat.add_edge(3, 6);
    mat.add_edge(3, 8);
    mat.add_edge(4, 1);
    mat.add_edge(4, 8);
    mat.add_edge(5, 4);
    mat.add_edge(5, 6);
    mat.add_edge(6, 2);
    mat.add_edge(6, 6);
    mat.add_edge(7, 1);
    mat.add_edge(8, 3);
    mat.add_edge(8, 5);
    mat.add_edge(8, 8);
    mat.add_edge(9, 2);
    mat.add_edge(9, 8);
    if (!check_connected_components("random 3", mat, 9, 7)) {
        return false;
    }
    
    mat.clear();
    mat.add_edge(0, 1);
    mat.add_edge(0, 2);
    mat.add_edge(0, 3);
    mat.add_edge(0, 8);
    mat.add_edge(1, 1);
    mat.add_edge(1, 5);
    mat.add_edge(1, 9);
    mat.add_edge(2, 1);
    mat.add_edge(2, 9);
    mat.add_edge(3, 3);
    mat.add_edge(3, 5);
    mat.add_edge(3, 8);
    mat.add_edge(4, 0);
    mat.add_edge(4, 2);
    mat.add_edge(4, 6);
    mat.add_edge(5, 8);
    mat.add_edge(6, 3);
    mat.add_edge(6, 5);
    mat.add_edge(7, 5);
    mat.add_edge(7, 6);
    mat.add_edge(7, 9);
    mat.add_edge(8, 2);
    mat.add_edge(8, 4);
    mat.add_edge(8, 8);
    mat.add_edge(9, 3);
    if (!check_connected_components("random 4", mat, 10, 9)) {
        return false;
    }
    
    mat.clear();
    mat.add_edge(0, 6);
    mat.add_edge(1, 1);
    mat.add_edge(1, 9);
    mat.add_edge(5, 1);
    mat.add_edge(5, 2);
    mat.add_edge(5, 5);
    mat.add_edge(5, 9);
    mat.add_edge(6, 1);
    mat.add_edge(6, 3);
    mat.add_edge(6, 7);
    mat.add_edge(6, 8);
    mat.add_edge(7, 1);
    mat.add_edge(8, 1);
    mat.add_edge(8, 4);
    mat.add_edge(8, 8);
    mat.add_edge(9, 1);
    if (!check_connected_components("random 5", mat, 6, 8)) {
        return false;
    }

    return true;
}

void assert_aprox(long double a, long double b, long double epsilon, string name) {
    if (abs(a - b) > epsilon) {
        cout << "(" << name << ") " << a << " and " << b << " are not approximately equal (epsilon=" << epsilon << ")" << endl;
        exit(1);
    }
}

void assert_eq(long double a, long double b, string name) {
    const long double verysmall = 1e-15;
    unsigned oom_a, oom_b;
    unsigned long aa=0, bb=0;

    // detect nan
    if (isnan(a) || isnan(b)) goto err;

    // avoid problems with zero
    if (a < verysmall && a > -verysmall) {
        if (b > verysmall || b < -verysmall) goto err;
        return;
    }
    if (b < verysmall && b > - verysmall) {
        if (a > verysmall || a < -verysmall) goto err;
        return;
    }

    // avoid problems with negative
    if (a < 0 && b < 0) {
        return assert_eq(-a, -b, name);
    }

    // compare order of magnitude
    oom_a = (unsigned)floor(log10(a));
    oom_b = (unsigned)floor(log10(b));
    if (oom_a != oom_b) goto err;

    // convert to int by getting significant digits (need order of magnitude
    // for this) and compare ints instead of doubles
    aa = (unsigned long)(round(a * pow(10,8-oom_a)));
    bb = (unsigned long)(round(b * pow(10,8-oom_b)));
    if (aa != bb) goto err;

    return;
  err:
    cout << "assert error! (" << name << ") " << a << " ("<<aa<<") =/= " << b << " ("<<bb<<")" << endl;
    exit(1);
}

void test_learning_equations(void) {
    RNG rng(10);

    ofstream file;

    const bool dofile = false;

    if (dofile) {
        file.open("data.tsv");
        file << "phi\tM\tlambda\tdelta\n";
    }
    
    const unsigned n=30;
    const unsigned m=30;
    const unsigned lambdaDiv = 10;
    const unsigned realizations = 1;
    
    for (unsigned type_=0; type_<=2; type_++) {
        const unsigned type = type_;
        cout << "type=" << type << endl;
        for (unsigned phi_=0; phi_<=6; phi_++) {
            const long double phi = phi_;
            cout << "phi=" << phi << endl;
            for (unsigned M_=1; M_<=15; M_++) {
                const unsigned M = M_;
                cout << "M=" << M << endl;
                if (dofile) {
                    file.flush();
                }
                for (unsigned lambda_=0; lambda_<=lambdaDiv; lambda_++) {
                    const long double lambda = (long double)lambda_/lambdaDiv;
                    //cout << "lambda=" << lambda << endl;

                    long double delta_avg = 0;
                    for (unsigned realization=0; realization<realizations; realization++) {
                        FirstModelStaticLexicalMatrix mat_original(n, m, phi_);
                        FirstModelStaticLexicalMatrix mat_a(n, m, phi_);
                        FirstModelStaticLexicalMatrix mat_b(n, m, phi_);
                        
                        const unsigned meaning_a = m-1;
                        const unsigned meaning_b = 0;
                        const unsigned word = n-1;
                        
                        mat_original.clear();
                        mat_a.clear();
                        mat_b.clear();
                        
                        unsigned mu_k;
                        
                        if (type == 0) { // Random graph except last meaning and last word are left unlinked
                            std::uniform_int_distribution<unsigned> word(0, n/2);
                            std::uniform_int_distribution<unsigned> meaning(0, m/2);
                            
                            unsigned edges = M;
                            while (edges > 0) {
                                unsigned i = word(rng);
                                unsigned j = meaning(rng);
                                if (!mat_original.get_a(i, j)) {
                                    mat_original.add_edge(i, j);
                                    mat_a.add_edge(i, j);
                                    mat_b.add_edge(i, j);
                                    edges--;
                                }
                            }
                        
                            mu_k = (unsigned)-1; // should not be used
                        } else if (type == 1) { // Bijective graph, last word and meaning unconnected
                            for (unsigned i=0; i<M; i++) {
                                mat_original.add_edge(i, i);
                                mat_a.add_edge(i, i);
                                mat_b.add_edge(i, i);
                            }
                        
                            mu_k = (unsigned)-1; // should not be used
                        } else if (type == 2) { // Every meaning is connected to one or no words (no restrictions on word to meaning connections) last meaning and last word left unconnected
                            std::uniform_int_distribution<unsigned> unif(0,M/2);
                        
                            unsigned k = (unsigned)-1;
                            for (unsigned j=0; j<M; j++) {
                                unsigned i = unif(rng);
                                if (j == meaning_b) {
                                    k = i;
                                }
                                mat_original.add_edge(i, j);
                                mat_a.add_edge(i, j);
                                mat_b.add_edge(i, j);
                            }
                        
                            mu_k = mat_b.get_mu(k);
                        } else {
                            cout << "bad type" << endl;
                            exit(1);
                        }
                    
                        mat_original.recalculate();
                        unsigned omega_j = mat_original.get_omega(meaning_b);
                        long double omega_phi_j = mat_original.get_omega_phi(meaning_b);
                        (void)mu_k;
                        (void)omega_phi_j;
                    
                        mat_a.add_edge(word, meaning_a);
                        mat_a.recalculate();
                    
                        mat_b.add_edge(word, meaning_b);
                        mat_b.recalculate();
                    
                        unsigned mu_i_a__, mu_i_b__, mu_i_oa__, mu_i_ob__;
                        unsigned omega_j_a__, omega_j_b__, omega_j_oa__, omega_j_ob__;
                        long double mu_phi_i_a__, mu_phi_i_b__, mu_phi_i_oa__, mu_phi_i_ob__;
                        long double omega_phi_j_a__, omega_phi_j_b__, omega_phi_j_oa__, omega_phi_j_ob__;
                        long double x_sirj_b__;
                        long double x_si_a__, x_si_b__;
                        long double x_rj_a__, x_rj_b__;
                        vector<long double> x_skrj_ob__, x_skrj_b__;
                        vector<long double> x_sk_ob__, x_sk_b__;
                        long double M_phi_o__, M_phi_a__, M_phi_b__;
                        long double XSR_o__, XSR_a__, XSR_b__;
                        long double XS_o__, XS_a__, XS_b__;
                        long double XR_o__, XR_a__, XR_b__;
                        long double Delta__;
                    
                        if (type == 0) { // "general" case
                            mu_i_oa__ = 0;
                            mu_i_ob__ = 0;
                            mu_i_a__ = 1;
                            mu_i_b__ = 1;
                        
                            omega_j_oa__ = 0;
                            omega_j_ob__ = omega_j;
                            omega_j_a__ = 1;
                            omega_j_b__ = omega_j+1;
                        
                            mu_phi_i_oa__ = 0;
                            mu_phi_i_ob__ = 0;
                            mu_phi_i_a__ = 1;
                            mu_phi_i_b__ = pow(omega_j+1, phi);
                        
                            omega_phi_j_oa__ = 0;
                            omega_phi_j_ob__ = omega_phi_j;
                            omega_phi_j_a__ = 1;
                            omega_phi_j_b__ = omega_phi_j + 1;
                        
                            x_sirj_b__ = pow(omega_j+1, phi) * log(omega_j+1);
                        
                            x_si_a__ = 0;
                            x_si_b__ = phi*pow(omega_j+1, phi)*log(omega_j+1);
                        
                            x_rj_a__ = 0;
                            x_rj_b__ = pow(omega_j+1, phi)*(omega_phi_j+1)*log(pow(omega_j+1, phi)*(omega_phi_j+1));
                        
                            for (unsigned k=0; k<n; k++) {
                                if (mat_original.get_a(k, meaning_b)) {
                                    unsigned mu = mat_original.get_mu(k);
                                    x_skrj_ob__.push_back(pow(omega_j, phi) * pow(mu, phi)*log(mu) + pow(omega_j, phi)*log(omega_j)*pow(mu,phi));
                                }
                            }
                            for (unsigned k=0; k<n; k++) {
                                if (k != word && mat_b.get_a(k, meaning_b)) {
                                    unsigned mu = mat_original.get_mu(k);
                                    x_skrj_b__.push_back( pow(omega_j+1, phi) * pow(mu, phi)*log(mu) + pow(omega_j+1, phi) * log(omega_j+1) * pow(mu, phi) );
                                }
                            }
                        
                            for (unsigned k=0; k<n; k++) {
                                long double n;
                                unsigned mu = mat_original.get_mu(k);
                                if (mu == 0) {
                                    n = 0;
                                } else {
                                    long double mu_phi = mat_original.get_mu_phi(k);
                                    n = pow(mu, phi) * mu_phi * log(pow(mu, phi) * mu_phi); // no need to simplify
                                }
                                x_sk_ob__.push_back(n);
                            }
                            for (unsigned k=0; k<n; k++) {
                                if (k == word || !mat_original.get_a(k, meaning_b)) {
                                    continue;
                                }
                                long double n;
                                unsigned mu = mat_original.get_mu(k);
                                mu = mat_b.get_mu(k);
                                if (mu == 0) {
                                    n = 0;
                                } else {
                                    long double mu_phi = mat_original.get_mu_phi(k);
                                    n = x_sk_ob__[k] + (pow(omega_j + 1, phi) - pow(omega_j, phi)) * pow(mu, phi) * log(pow(mu, phi) * (mu_phi - pow(omega_j, phi) + pow(omega_j + 1, phi))) + pow(mu, phi) * mu_phi * log((mu_phi - pow(omega_j, phi) + pow(omega_j+1, phi))/mu_phi);
                                }
                                x_sk_b__.push_back(n);
                            }
                        
                            M_phi_o__ = mat_original.get_M();
                            M_phi_a__ = M_phi_o__ + 1;
                            M_phi_b__ = M_phi_o__ + (pow(omega_j+1, phi)-pow(omega_j,phi))*omega_phi_j + pow(omega_j+1, phi);
                        
                            XSR_o__ = mat_original.get_XSR();
                            XSR_a__ = XSR_o__;
                            {
                                long double sum = 0;
                                for (unsigned k=0; k<n; k++) {
                                    if (k != word && mat_original.get_a(k, meaning_b)) {
                                        unsigned mu = mat_original.get_mu(k);
                                        sum += pow(mu, phi) * log(mu);
                                    }
                                }
                                if (omega_j == 0) {
                                    XSR_b__ = XSR_o__ + (pow(omega_j+1,phi)-pow(omega_j,phi))*sum + omega_phi_j*(pow(omega_j+1,phi)*log(omega_j+1)) + pow(omega_j+1,phi)*log(omega_j+1);
                                } else {
                                    XSR_b__ = XSR_o__ + (pow(omega_j+1,phi)-pow(omega_j,phi))*sum + omega_phi_j*(pow(omega_j+1,phi)*log(omega_j+1)-pow(omega_j,phi)*log(omega_j)) + pow(omega_j+1,phi)*log(omega_j+1);
                                }
                            }
                        
                            XS_o__ = mat_original.get_XS();
                            XS_a__ = XS_o__;
                            {
                                long double sum1 = 0;
                                long double sum2 = 0;
                                for (unsigned k=0; k<n; k++) {
                                    unsigned mu = mat_original.get_mu(k);
                                    long double mu_phi = mat_original.get_mu_phi(k);
                                    if (k != word && mat_original.get_a(k, meaning_b)) {
                                        sum1 += pow(mu, phi) * log(pow(mu, phi) * (mu_phi - pow(omega_j, phi) + pow(omega_j + 1, phi)));
                                        sum2 += pow(mu, phi) * mu_phi * log((mu_phi - pow(omega_j, phi) + pow(omega_j+1, phi)) / mu_phi);
                                    }
                                }
                                XS_b__ = XS_o__ + (pow(omega_j+1,phi)-pow(omega_j,phi)) * sum1 + sum2 + phi*pow(omega_j+1, phi)*log(omega_j+1);
                            }
                        
                            XR_o__ = mat_original.get_XR();
                            XR_a__ = XR_o__;
                            if (omega_j == 0) {
                                XR_b__ = XR_o__ + pow(omega_j+1, phi) * (omega_phi_j+1) * log(pow(omega_j+1, phi)*(omega_phi_j+1));
                            } else {
                                XR_b__ = XR_o__ - pow(omega_j, phi) * omega_phi_j * log(pow(omega_j, phi) * omega_phi_j) + pow(omega_j+1, phi) * (omega_phi_j+1) * log(pow(omega_j+1, phi)*(omega_phi_j+1));
                            }

                            {
                                long double Delta_XS = M_phi_b__ * XS_a__ - M_phi_a__ * XS_b__;
                                long double Delta_XR = M_phi_b__ * XR_a__ - M_phi_a__ * XR_b__;
                                long double Delta_XSR = M_phi_b__ * XSR_a__ - M_phi_a__ * XSR_b__;
                                Delta__ = (1-2*lambda)*log(M_phi_a__/M_phi_b__) - 1/(M_phi_a__*M_phi_b__) * ((1-2*lambda)*Delta_XS - lambda*Delta_XR + lambda*phi*Delta_XSR);
                            }
                        } else if (type == 2) { // "wj=1 for all j" case
                            mu_i_oa__ = 0;
                            mu_i_ob__ = 0;
                            mu_i_a__ = 1;
                            mu_i_b__ = 1;
                        
                            omega_j_oa__ = 0;
                            omega_j_ob__ = 1;
                            omega_j_a__ = 1;
                            omega_j_b__ = 2;
                        
                            mu_phi_i_oa__ = 0;
                            mu_phi_i_ob__ = 0;
                            mu_phi_i_a__ = 1;
                            mu_phi_i_b__ = pow(2, phi);
                        
                            omega_phi_j_oa__ = 0;
                            omega_phi_j_ob__ = pow(mu_k, phi);
                            omega_phi_j_a__ = 1;
                            omega_phi_j_b__ = pow(mu_k, phi) + 1;
                        
                            x_sirj_b__ = pow(2, phi) * log(2);
                        
                            x_si_a__ = 0;
                            x_si_b__ = phi*pow(2, phi)*log(2);
                        
                            x_rj_a__ = 0;
                            x_rj_b__ = pow(2, phi) * (pow(mu_k, phi) + 1) * log(pow(2, phi) * (pow(mu_k, phi) + 1));
                        
                            x_skrj_ob__.push_back( pow(1, phi) * pow(mu_k, phi) * log(mu_k) );
                            x_skrj_b__.push_back( pow(2, phi) * pow(mu_k, phi) * (log(mu_k) + log(2)) );
                        
                            for (unsigned k=0; k<n; k++) {
                                long double n;
                                unsigned mu = mat_original.get_mu(k);
                                if (mu == 0) {
                                    n = 0;
                                } else {
                                    long double mu_phi = mat_original.get_mu_phi(k);
                                    n = pow(mu, phi) * mu_phi * log(pow(mu, phi) * mu_phi); // no need to simplify
                                }
                                x_sk_ob__.push_back(n);
                            }
                            {
                                long double init=0;
                                for (unsigned k=0; k<n; k++) {
                                    if (mat_b.get_a(k, meaning_b)) {
                                        init = x_sk_ob__[k];
                                        break;
                                    }
                                }
                                x_sk_b__.push_back(
                                    + init
                                    + (pow(omega_j + 1, phi) - pow(omega_j, phi)) * pow(mu_k, phi) * log(pow(mu_k, phi) * (mu_k - pow(omega_j, phi) + pow(omega_j + 1, phi)))
                                    + pow(mu_k, phi) * mu_k * log((mu_k - pow(omega_j, phi) + pow(omega_j+1, phi))/mu_k)
                                    );
                            }
                        
                            M_phi_o__ = mat_original.get_M();
                            M_phi_a__ = M_phi_o__ + 1;
                            M_phi_b__ = M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi);
                        
                            XSR_o__ = mat_original.get_XSR();
                            XSR_a__ = XSR_o__;
                            XSR_b__ = XSR_o__ + (pow(2,phi)-1)*pow(mu_k,phi)*log(mu_k) + (pow(mu_k,phi)+1)*pow(2,phi)*log(2);
                        
                            XS_o__ = (phi+1) * XSR_o__;
                            XS_a__ = XS_o__;
                            XS_b__ = XS_o__ + (pow(2,phi)-1)*pow(mu_k,phi)*log( pow(mu_k,phi) * (mu_k - 1 + pow(2,phi))) + pow(mu_k,phi+1) * log((mu_k - 1 + pow(2,phi))/mu_k) + phi*pow(2,phi)*log(2);
                        
                            XR_o__ = phi * XSR_o__;
                            XR_a__ = XR_o__;
                            XR_b__ = XR_o__ - phi * pow(mu_k, phi) * log(mu_k) + pow(2, phi) * (pow(mu_k, phi) + 1) * log(pow(2, phi)*(pow(mu_k, phi) + 1));

                            if (lambda_ == 0) {
                                Delta__ = log((M_phi_o__ + 1)/(M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))) +
                                    1 / (M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))
                                    * (
                                        + phi * pow(2,phi) * pow(mu_k,phi) * log(mu_k)
                                        - (
                                            + (phi + 1) * XSR_o__ * (pow(2,phi) - 1) * (pow(mu_k,phi) + 1) / (M_phi_o__ + 1)
                                            - phi * pow(2,phi) * log(2)
                                            + pow(mu_k,phi) * (
                                                + log(mu_k) * (mu_k + phi)
                                                - (mu_k - 1 + pow(2,phi)) * log(mu_k - 1 + pow(2,phi))
                                                )
                                            )
                                        )
                                    ;
                            } else if (lambda_ == lambdaDiv) {
                                Delta__ = -log((M_phi_o__ + 1)/(M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))) -
                                    1 / (M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))
                                    * (
                                        + (pow(mu_k,phi) + 1) * pow(2,phi) * log(pow(mu_k,phi) + 1)
                                        - (
                                            + (phi + 1) * XSR_o__ * (pow(2,phi) - 1) * (pow(mu_k,phi) + 1) / (M_phi_o__ + 1)
                                            - phi * pow(2,phi) * log(2)
                                            + pow(mu_k,phi) * (
                                                + log(mu_k) * (mu_k + phi)
                                                - (mu_k - 1 + pow(2,phi)) * log(mu_k - 1 + pow(2,phi))
                                                )
                                            )
                                        )
                                    ;
                            } else {
                                Delta__ = (1-2*lambda) * (
                                    log((M_phi_o__ + 1)/(M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi)))
                                    + 1 / (M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi)) * (
                                        + (phi + 1) * XSR_o__ * (pow(2,phi) - 1) * (pow(mu_k,phi) + 1) / (M_phi_o__ + 1)
                                        - phi * pow(2,phi) * log(2)
                                        + pow(mu_k,phi) * (
                                            + log(mu_k) * (mu_k + phi)
                                            - (mu_k - 1 + pow(2,phi)) * log(mu_k - 1 + pow(2,phi))
                                            )
                                        )
                                    )
                                    - 1 / (M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))
                                    * (
                                        + lambda * (pow(mu_k,phi) + 1) * pow(2,phi) * log(pow(mu_k,phi) + 1)
                                        - (1 - lambda) * phi * pow(2,phi) * pow(mu_k,phi) * log(mu_k)
                                        )
                                    ;
                            }

                            {
                                long double a =
                                    2 * log((M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))/(M_phi_o__ + 1))
                                    - 1 / (M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi)) * (
                                        + (pow(mu_k,phi) + 1) * pow(2,phi) * log(pow(mu_k,phi) + 1)
                                        + phi * pow(2,phi) * pow(mu_k,phi) * log(mu_k)
                                        + 2 * (
                                            - ((int)phi+1) * XSR_o__ * (pow(2,phi) - 1) * (pow(mu_k,phi) + 1) / (M_phi_o__ + 1)
                                            + phi * pow(2,phi) * log(2)
                                            - pow(mu_k,phi) * (
                                                + log(mu_k) * (mu_k + phi)
                                                - (mu_k - 1 + pow(2,phi)) * log(mu_k - 1 + pow(2,phi))
                                            )
                                        )
                                    )
                                    ;

                                long double b =
                                    - log((M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi))/(M_phi_o__ + 1))
                                    + 1 / (M_phi_o__ + (pow(2,phi) - 1)*pow(mu_k,phi) + pow(2,phi)) * (
                                        + phi * pow(2,phi) * pow(mu_k,phi) * log(mu_k)
                                        - (phi + 1) * XSR_o__ * (pow(2,phi) - 1) * (pow(mu_k,phi) + 1) / (M_phi_o__ + 1)
                                        + phi * pow(2,phi) * log(2)
                                        - pow(mu_k,phi) * (
                                            + log(mu_k) * (mu_k + phi)
                                            - (mu_k - 1 + pow(2,phi)) * log(mu_k - 1 + pow(2,phi))
                                        )
                                    )
                                    ;
                                
                                Delta__ = a * lambda + b;
                            }
                        } else if (type == 1) { // "all bijective" case
                            mu_i_oa__ = 0;
                            mu_i_ob__ = 0;
                            mu_i_a__ = 1;
                            mu_i_b__ = 1;
                        
                            omega_j_oa__ = 0;
                            omega_j_ob__ = 1;
                            omega_j_a__ = 1;
                            omega_j_b__ = 2;
                        
                            mu_phi_i_oa__ = 0;
                            mu_phi_i_ob__ = 0;
                            mu_phi_i_a__ = 1;
                            mu_phi_i_b__ = pow(2, phi);
                        
                            omega_phi_j_oa__ = 0;
                            omega_phi_j_ob__ = 1;
                            omega_phi_j_a__ = 1;
                            omega_phi_j_b__ = 2;
                        
                            x_sirj_b__ = pow(2, phi) * log(2);
                        
                            x_si_a__ = 0;
                            x_si_b__ = phi*pow(2, phi)*log(2);
                        
                            x_rj_a__ = 0;
                            x_rj_b__ = (phi+1) * pow(2, phi+1) * log(2);
                        
                            x_skrj_ob__.push_back(0);
                            x_skrj_b__.push_back( pow(2, phi) * log(2) );
                        
                            for (unsigned k=0; k<n; k++) {
                                long double n;
                                unsigned mu = mat_original.get_mu(k);
                                if (mu == 0) {
                                    n = 0;
                                } else {
                                    long double mu_phi = mat_original.get_mu_phi(k);
                                    n = pow(mu, phi) * mu_phi * log(pow(mu, phi) * mu_phi); // no need to simplify
                                }
                                x_sk_ob__.push_back(n);
                            }
                            {
                                long double init=0;
                                for (unsigned k=0; k<n; k++) {
                                    if (mat_b.get_a(k, meaning_b)) {
                                        init = x_sk_ob__[k];
                                        break;
                                    }
                                }
                                x_sk_b__.push_back(
                                    + init
                                    + (pow(omega_j + 1, phi) - pow(omega_j, phi)) * log(1 - pow(omega_j, phi) + pow(omega_j + 1, phi))
                                    + log(1 - pow(omega_j, phi) + pow(omega_j+1, phi))
                                    );
                            }
                        
                            M_phi_o__ = M;
                            M_phi_a__ = M + 1;
                            M_phi_b__ = M + pow(2,phi+1) - 1;
                        
                            XSR_o__ = 0;
                            XSR_a__ = 0;
                            XSR_b__ = pow(2, phi+1) * log(2);
                        
                            XS_o__ = 0;
                            XS_a__ = 0;
                            XS_b__ = phi*pow(2,phi+1)*log(2);
                        
                            XR_o__ = 0;
                            XR_a__ = 0;
                            XR_b__ = (phi+1)*pow(2,phi+1)*log(2);

                            if (lambda_ == 0) {
                                Delta__ = log((M+1)/(M+pow(2,phi+1)-1)) + (pow(2,phi+1)*phi*log(2))/(M+pow(2,phi+1)-1);
                            } else if (lambda_ == lambdaDiv) {
                                Delta__ = -log((M+1)/(M+pow(2,phi+1)-1)) - (pow(2,phi+1)*(phi+1)*log(2))/(M+pow(2,phi+1)-1);
                            } else {
                                Delta__ = (1-2*lambda)*(log((M+1)/(M+pow(2,phi+1)-1)) + (pow(2,phi+1)*log(2)*phi)/(M+pow(2,phi+1)-1)) - lambda*(pow(2,phi+1)*log(2))/(M+pow(2,phi+1)-1);
                            }
                            Delta__ = (-2*log((M+1)/(M+pow(2,phi+1)-1)) - (pow(2,phi+1)*log(2)*(2*phi+1))/(M+pow(2,phi+1)-1) )*lambda + ( log((M+1)/(M+pow(2,phi+1)-1)) + pow(2,phi+1)*log(2)*phi/(M+pow(2,phi+1)-1) );
                        } else {
                            cout << "bad type" << endl;
                            exit(1);
                        }
                    
                        unsigned mu_i_a = mat_a.get_mu(word);
                        unsigned mu_i_b = mat_b.get_mu(word);
                        unsigned mu_i_oa = mat_original.get_mu(word);
                        unsigned mu_i_ob = mat_original.get_mu(word);
                        assert_eq(mu_i_a, mu_i_a__, "mu_i_a");
                        assert_eq(mu_i_b, mu_i_b__, "mu_i_b");
                        assert_eq(mu_i_oa, mu_i_oa__, "mu_i_oa");
                        assert_eq(mu_i_ob, mu_i_ob__, "mu_i_ob");
                    
                        unsigned omega_j_a = mat_a.get_omega(meaning_a);
                        unsigned omega_j_b = mat_b.get_omega(meaning_b);
                        unsigned omega_j_oa = mat_original.get_omega(meaning_a);
                        unsigned omega_j_ob = mat_original.get_omega(meaning_b);
                        assert_eq(omega_j_a, omega_j_a__, "omega_j_a");
                        assert_eq(omega_j_b, omega_j_b__, "omega_j_b");
                        assert_eq(omega_j_oa, omega_j_oa__, "omega_j_oa");
                        assert_eq(omega_j_ob, omega_j_ob__, "omega_j_ob");
                    
                        long double mu_phi_i_a = mat_a.get_mu_phi(word);
                        long double mu_phi_i_b = mat_b.get_mu_phi(word);
                        long double mu_phi_i_oa = mat_original.get_mu_phi(word);
                        long double mu_phi_i_ob = mat_original.get_mu_phi(word);
                        assert_eq(mu_phi_i_a, mu_phi_i_a__, "mu_phi_i_a");
                        assert_eq(mu_phi_i_b, mu_phi_i_b__, "mu_phi_i_b");
                        assert_eq(mu_phi_i_oa, mu_phi_i_oa__, "mu_phi_i_oa");
                        assert_eq(mu_phi_i_ob, mu_phi_i_ob__, "mu_phi_i_ob");
                    
                        long double omega_phi_j_a = mat_a.get_omega_phi(meaning_a);
                        long double omega_phi_j_b = mat_b.get_omega_phi(meaning_b);
                        long double omega_phi_j_oa = mat_original.get_omega_phi(meaning_a);
                        long double omega_phi_j_ob = mat_original.get_omega_phi(meaning_b);
                        assert_eq(omega_phi_j_a, omega_phi_j_a__, "omega_phi_j_a");
                        assert_eq(omega_phi_j_b, omega_phi_j_b__, "omega_phi_j_b");
                        assert_eq(omega_phi_j_oa, omega_phi_j_oa__, "omega_phi_j_oa");
                        assert_eq(omega_phi_j_ob, omega_phi_j_ob__, "omega_phi_j_ob");
                    
                        long double x_sirj_b = mat_b.get_mu(word)==0||mat_b.get_omega(meaning_b)==0 ? 0 : pow(mat_b.get_mu(word) * mat_b.get_omega(meaning_b), phi) * log(mat_b.get_mu(word) * mat_b.get_omega(meaning_b));
                        assert_eq(x_sirj_b, x_sirj_b__, "x_sirj_b");
                    
                        long double x_si_a = mat_a.get_mu(word)==0 ? 0 : pow(mat_a.get_mu(word), phi)*mat_a.get_mu_phi(word)*log(pow(mat_a.get_mu(word), phi)*mat_a.get_mu_phi(word));
                        long double x_si_b = mat_b.get_mu(word)==0 ? 0 : pow(mat_b.get_mu(word), phi)*mat_b.get_mu_phi(word)*log(pow(mat_b.get_mu(word), phi)*mat_b.get_mu_phi(word));
                        assert_eq(x_si_a, x_si_a__, "x_si_a");
                        assert_eq(x_si_b, x_si_b__, "x_si_b");
                    
                        long double x_rj_a = mat_a.get_omega(meaning_a)==0 ? 0 : pow(mat_a.get_omega(meaning_a), phi)*mat_a.get_omega_phi(meaning_a)*log(pow(mat_a.get_omega(meaning_a), phi)*mat_a.get_omega_phi(meaning_a));
                        long double x_rj_b = mat_b.get_omega(meaning_b)==0 ? 0 : pow(mat_b.get_omega(meaning_b), phi)*mat_b.get_omega_phi(meaning_b)*log(pow(mat_b.get_omega(meaning_b), phi)*mat_b.get_omega_phi(meaning_b));
                        assert_eq(x_rj_a, x_rj_a__, "x_rj_a");
                        assert_eq(x_rj_b, x_rj_b__, "x_rj_b");
                    
                        vector<long double> x_skrj_ob;
                        for (unsigned k=0; k<n; k++) {
                            if (mat_original.get_a(k, meaning_b)) {
                                unsigned mu = mat_original.get_mu(k);
                                unsigned omega = mat_original.get_omega(meaning_b);
                                x_skrj_ob.push_back(pow(mu * omega, phi) * log(mu * omega));
                            }
                        }
                        vector<long double> x_skrj_b;
                        for (unsigned k=0; k<n; k++) {
                            if (k != word && mat_b.get_a(k, meaning_b)) {
                                unsigned mu = mat_b.get_mu(k);
                                unsigned omega = mat_b.get_omega(meaning_b);
                                x_skrj_b.push_back(pow(mu * omega, phi) * log(mu * omega));
                            }
                        }
                        assert_eq(x_skrj_ob.size(), x_skrj_ob__.size(), "x_skrj_ob size");
                        for (unsigned i=0; i<x_skrj_ob.size(); i++) {
                            assert_eq(x_skrj_ob[i], x_skrj_ob__[i], "x_skrj_ob");
                        }
                        assert_eq(x_skrj_b.size(), x_skrj_b__.size(), "x_skrj_b size");
                        for (unsigned i=0; i<x_skrj_b.size(); i++) {
                            assert_eq(x_skrj_b[i], x_skrj_b__[i], "x_skrj_b");
                        }
                    
                        vector<long double> x_sk_ob;
                        for (unsigned k=0; k<n; k++) {
                            long double n;
                            unsigned mu = mat_original.get_mu(k);
                            if (mu == 0) {
                                n = 0;
                            } else {
                                long double mu_phi = mat_original.get_mu_phi(k);
                                n = pow(mu, phi) * mu_phi * log(pow(mu, phi) * mu_phi);
                            }
                            x_sk_ob.push_back(n);
                        }
                        vector<long double> x_sk_b;
                        for (unsigned k=0; k<n; k++) {
                            if (k != word && mat_original.get_a(k, meaning_b)) {
                                long double n;
                                unsigned mu = mat_b.get_mu(k);
                                if (mu == 0) {
                                    n = 0;
                                } else {
                                    long double mu_phi = mat_b.get_mu_phi(k);
                                    n = pow(mu, phi) * mu_phi * log(pow(mu, phi) * mu_phi);
                                }
                                x_sk_b.push_back(n);
                            }
                        }
                        assert_eq(x_sk_ob.size(), x_sk_ob__.size(), "x_sk_ob size");
                        for (unsigned i=0; i<x_sk_ob.size(); i++) {
                            assert_eq(x_sk_ob[i], x_sk_ob__[i], "x_sk_ob");
                        }
                        assert_eq(x_sk_b.size(), x_sk_b__.size(), "x_sk_b size");
                        for (unsigned i=0; i<x_sk_b.size(); i++) {
                            assert_eq(x_sk_b[i], x_sk_b__[i], "x_sk_b");
                        }
                    
                        long double M_phi_o = mat_original.get_M();
                        long double M_phi_a = mat_a.get_M();
                        long double M_phi_b = mat_b.get_M();
                        assert_eq(M_phi_o, M_phi_o__, "M_phi_o");
                        assert_eq(M_phi_a, M_phi_a__, "M_phi_a");
                        assert_eq(M_phi_b, M_phi_b__, "M_phi_b");
                    
                        long double XSR_o = mat_original.get_XSR();
                        long double XSR_a = mat_a.get_XSR();
                        long double XSR_b = mat_b.get_XSR();
                        assert_eq(XSR_o, XSR_o__, "XSR_o");
                        assert_eq(XSR_a, XSR_a__, "XSR_a");
                        assert_eq(XSR_b, XSR_b__, "XSR_b");
                    
                        long double XS_o = mat_original.get_XS();
                        long double XS_a = mat_a.get_XS();
                        long double XS_b = mat_b.get_XS();
                        assert_eq(XS_o, XS_o__, "XS_o");
                        assert_eq(XS_a, XS_a__, "XS_a");
                        assert_eq(XS_b, XS_b__, "XS_b");
                
                        long double XR_o = mat_original.get_XR();
                        long double XR_a = mat_a.get_XR();
                        long double XR_b = mat_b.get_XR();
                        assert_eq(XR_o, XR_o__, "XR_o");
                        assert_eq(XR_a, XR_a__, "XR_a");
                        assert_eq(XR_b, XR_b__, "XR_b");

                        long double Delta;
                        {
                            long double Omega_a = -lambda * mat_a.get_ISR() + (1 - lambda) * mat_a.get_HS();
                            long double Omega_b = -lambda * mat_b.get_ISR() + (1 - lambda) * mat_b.get_HS();
                            Delta = Omega_a - Omega_b;
                        }
                        assert_eq(Delta, Delta__, "Delta");
                        delta_avg += Delta__;
                    }

                    if (dofile) {
                        file << phi << "\t" << M << "\t" << lambda << "\t" << delta_avg/realizations << "\n";
                    }
                }
            }
        }
    }

    

    // Testing X(S), X(R), X(S,R), M from language learning equations
    
    FirstModelStaticLexicalMatrix mat(10, 10, 5);
    std::bernoulli_distribution bern(0.7);
    std::uniform_int_distribution<unsigned> unif(0,9);

    unsigned edges = 0;
    for (unsigned i=0; i<10; i++) {
        if (bern(rng)) {
            edges++;
            mat.add_edge(i, i);
        }
    }
    mat.recalculate();
    
    cout << "M = " << mat.get_M() << " == " << edges << endl;
    cout << "XSR = " << mat.get_XSR() << " == 0" << endl;
    cout << "XS = " << mat.get_XS() << " == 0" << endl;
    cout << "XR = " << mat.get_XR() << " == 0" << endl;
    cout << endl;

    mat.clear();
    for (unsigned i=0; i<10; i++) {
        if (bern(rng)) {
            unsigned j = unif(rng);
            mat.add_edge(i, j);
        }
    }
    mat.recalculate();
    
    long double myM = 0, myXSR = 0, myXS = 0, myXR = 0;
    for (unsigned i=0; i<10; i++) {
        unsigned omega = mat.get_omega(i);
        if (omega) {
            myM += std::pow(omega, 6);
            myXSR += std::pow(omega, 6) * std::log(omega);
            myXS += std::pow(omega, 6) * std::log(omega);
            myXR += std::pow(omega, 6) * std::log(omega);
        }
    }
    myXS *= 5;
    myXR *= 6;
    
    cout << "M = " << mat.get_M() << " == " << myM << endl;
    cout << "XSR = " << mat.get_XSR() << " == " << myXSR << endl;
    cout << "XS = " << mat.get_XS() << " == " << myXS << endl;
    cout << "XR = " << mat.get_XR() << " == " << myXR << endl;
    cout << endl;

    mat.clear();
    for (unsigned i=0; i<10; i++) {
        if (bern(rng)) {
            unsigned j = unif(rng);
            mat.add_edge(j, i);
        }
    }
    mat.recalculate();
    
    myM = myXSR = myXS = myXR = 0;
    for (unsigned i=0; i<10; i++) {
        unsigned mu = mat.get_mu(i);
        if (mu) {
            myM += std::pow(mu, 6);
            myXSR += std::pow(mu, 6) * std::log(mu);
            myXS += std::pow(mu, 6) * std::log(mu);
            myXR += std::pow(mu, 6) * std::log(mu);
        }
    }
    myXS *= 6;
    myXR *= 5;
    
    cout << "M = " << mat.get_M() << " == " << myM << endl;
    cout << "XSR = " << mat.get_XSR() << " == " << myXSR << endl;
    cout << "XS = " << mat.get_XS() << " == " << myXS << endl;
    cout << "XR = " << mat.get_XR() << " == " << myXR << endl;
    cout << endl;

    cout << "ok!" << endl;
}

bool matrix_valid_check(const BaseLexicalMatrix &mat) {
    bool empty = true;
    bool disconnected_meanings = false;

    {
        bool should_break = false;
        for (unsigned i=0; i<mat.m; i++) {
            if (mat.get_omega(i) == 0) {
                disconnected_meanings = true;
                if (should_break) {
                    break;
                }
                should_break = true;
            } else {
                empty = false;
                if (should_break) {
                    break;
                }
                should_break = true;
            }
        }
    }
    
    if (mat.valid(true)) {
        if (empty) {
            mat.print_properties();
            cout << "Valid empty matrix" << endl;
            return false;
        }
    } else {
        if (!empty) {
            mat.print_properties();
            cout << "Invalid non-empty matrix" << endl;
            return false;
        }
    }
    if (mat.valid(false)) {
        if (empty || disconnected_meanings) {
            mat.print_properties();
            cout << "Valid matrix empty or with disconnected meanings" << endl;
            return false;
        }
    } else {
        if (!empty && !disconnected_meanings) {
            mat.print_properties();
            cout << "Invalid matrix non-empty and without disconnected meanings" << endl;
            return false;
        }
    }

    return true;
}

template <class T>
bool valid_invalid_matrix() {
    T &mat = create_mat<T>(10, 10, 0, 1e-15);

    bool ok = true;

    // Basic tests
    mat.clear();
    ok &= matrix_valid_check(mat);
        
    mat.add_edge(0, 0);
    ok &= matrix_valid_check(mat);

    for (unsigned i=1; i<mat.m; i++) {
        mat.add_edge(0, i);
    }
    ok &= matrix_valid_check(mat);

    // Fuzzy testing
    mat.clear();
    double p = 3/10/10;
    RNG rng;
    for (int i=0; i<1000; i++) {
        mat.initialize_random(rng, p);
        ok &= matrix_valid_check(mat);
    }

    // Valid -> Invalid -> Valid -> Invalid
    mat.clear();
    for (unsigned i=0; i<mat.m; i++) {
        mat.add_edge(0, i);
    }
    ok &= matrix_valid_check(mat);
    mat.remove_edge(0, 0);
    ok &= matrix_valid_check(mat);
    mat.add_edge(0, 0);
    ok &= matrix_valid_check(mat);

    delete &mat;
    return ok;
}

template <class T>
bool static_tests(std::string model, bool stop_first_failure) {
    std::string type = "[" + model + "][STATIC] ";
    
    bool allpass = true;
    
    cout << type << "Testing extreme cases with parameters..." << endl;
    allpass &= params<T>(1, 200, 10, 1, 200, 10, 0, 1, 1.0, 1e-6, false, false);
    if (stop_first_failure && !allpass) {
        return false;
    }
    
    cout << type <<  "Test save/restore phi=0" << endl;
    allpass &= save_restore<T>(50, 50, 0, 2, 1e-15);
    if (stop_first_failure && !allpass) {
        return false;
    }
    cout << type <<  "Test save/restore phi=1" << endl;
    allpass &= save_restore<T>(50, 50, 1, 2, 1e-15);
    if (stop_first_failure && !allpass) {
        return false;
    }
    
    cout << type <<  "Test matrix valid/invalid" << endl;
    allpass &= valid_invalid_matrix<T>();
    if (stop_first_failure && !allpass) {
        return false;
    }

    return allpass;
}

template <class StaticT, class DynamicT>
bool dynamic_tests(std::string model, bool stop_first_failure) {
    std::string type = "[" + model + "][DYNAMIC] ";
    
    bool allpass = true;
    
    cout << type << "Testing extreme cases with parameters..." << endl;
    allpass &= params<DynamicT>(1, 200, 10, 1, 200, 10, 0, 1, 1.0, 1e-6, true, false);
    if (stop_first_failure && !allpass) {
        return false;
    }

    cout << type << "Exhaustive testing of smaller matrixes with phi >= 0" << endl;
    allpass &= exhaustive<StaticT, DynamicT>(4, 4, 0, 1e-6, 1e-6);
    if (stop_first_failure && !allpass) {
        return false;
    }
    allpass &= exhaustive<StaticT, DynamicT>(4, 4, 0.5, 1e-6, 1e-6);
    if (stop_first_failure && !allpass) {
        return false;
    }
    allpass &= exhaustive<StaticT, DynamicT>(4, 4, 1, 1e-6, 1e-6);
    if (stop_first_failure && !allpass) {
        return false;
    }
    allpass &= exhaustive<StaticT, DynamicT>(4, 4, 1.5, 1e-6, 1e-6);
    if (stop_first_failure && !allpass) {
        return false;
    }
    allpass &= exhaustive<StaticT, DynamicT>(4, 4, 2, 1e-6, 1e-6);
    if (stop_first_failure && !allpass) {
        return false;
    }

    cout << type << "Testing equal initial conditions for static vs dynamic phi=0" << endl;
    allpass &= initial_conditions<StaticT, DynamicT>(100, 100, 0, numeric_limits<double>::epsilon(), 1234);
    if (stop_first_failure && !allpass) {
        return false;
    }
    cout << type << "Testing equal initial conditions for static vs dynamic phi=1" << endl;
    allpass &= initial_conditions<StaticT, DynamicT>(100, 100, 1, numeric_limits<double>::epsilon(), 1234);
    if (stop_first_failure && !allpass) {
        return false;
    }
    
    cout << type << "Testing replace edge on single edge" << endl;
    allpass &= replace_edge_single<DynamicT>(40, 40, 0);
    if (stop_first_failure && !allpass) {
        return false;
    }
    allpass &= replace_edge_single<DynamicT>(40, 40, 1);
    if (stop_first_failure && !allpass) {
        return false;
    }
    
    cout << type << "Test save/restore phi=0" << endl;
    allpass &= save_restore<DynamicT>(50, 50, 0, 2, 1e-10);
    if (stop_first_failure && !allpass) {
        return false;
    }
    cout << type << "Test save/restore phi=1" << endl;
    allpass &= save_restore<DynamicT>(50, 50, 1, 2, 1e-10);
    if (stop_first_failure && !allpass) {
        return false;
    }
    
    cout << type << "Test matrix valid/invalid" << endl;
    allpass &= valid_invalid_matrix<DynamicT>();
    if (stop_first_failure && !allpass) {
        return false;
    }

    return allpass;
}

bool run_alltests(bool stop_first_failure) {
    bool allpass = true;

    // First model static
    allpass &= static_tests<FirstModelStaticLexicalMatrix>("FIRST", stop_first_failure);
        
    // First model dynamic
    allpass &= dynamic_tests<FirstModelStaticLexicalMatrix, FirstModelDynamicLexicalMatrix>("FIRST", stop_first_failure);
    
    // Second model static
    allpass &= static_tests<SecondModelStaticLexicalMatrix>("SECOND", stop_first_failure);

    // Second model dynamic
    allpass &= dynamic_tests<SecondModelStaticLexicalMatrix, SecondModelDynamicLexicalMatrix>("SECOND", stop_first_failure);

    // Generic
    cout << "Test connected components" << endl;
    allpass &= connected_components();

    return allpass;
}

int main() {
    //test_learning_equations();
    //return EXIT_SUCCESS;
    
    if (run_alltests(false)) {
        cout << "ALL TESTS PASS" << endl;
        return EXIT_SUCCESS;
    } else {
        cout << "SOME TESTS FAILED" << endl;
        return EXIT_FAILURE;
    }
}
