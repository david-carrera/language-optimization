#include <SecondModelStaticLexicalMatrix.hpp>
#include <algorithm>
#include <cassert>
#include <cstring>

std::vector<long double> vectorized_log(const std::vector<long double> &x) {
    std::vector<long double> v(x.size());
    for (unsigned i=0; i<x.size(); i++) {
        v[i] = std::log(x[i]);
    }
    return v;
}

SecondModelStaticLexicalMatrix::SecondModelStaticLexicalMatrix(
    unsigned n, unsigned m, const std::vector<long double> &pi, double phi, long double epsilon) :
    BaseLexicalMatrix(n, m, phi, epsilon),
    omega_phi(m), chi(n), nu(m),
    rho(0.0), XS(0.0), XR(0.0), XSR(0.0),
    pi(pi), logpi(vectorized_log(pi)),
    HR_pi(-std::accumulate(pi.begin(), pi.end(), 0.0L, [](long double sum, long double p) {
        return sum + p * std::log(p);
    })) {
    // Ensure a priori meaning probabilities add up to 1
    assert(std::abs(std::accumulate(pi.begin(), pi.end(), 0.0L) - 1) < epsilon);
}

void SecondModelStaticLexicalMatrix::recalculate(bool force) {
    (void)force;
    const long double phi = this->phi;

    Accumulator rho;
    Accumulator HS;
    Accumulator HR;
    Accumulator HSR;
    
    std::fill(this->omega_phi.begin(), this->omega_phi.end(), 0.0L);
    std::fill(this->chi.begin(), this->chi.end(), 0.0L);
    std::fill(this->nu.begin(), this->nu.end(), 0.0L);

    unsigned seen_s[this->n];
    unsigned seen_r[this->m];
    Accumulator omega_phi[this->m];
    Accumulator nu[this->m];
    Accumulator chi[this->n];

    std::memset(seen_s, 0, sizeof(seen_s));
    std::memset(seen_r, 0, sizeof(seen_r));
    
    std::memset((void*)omega_phi, 0, sizeof(omega_phi));
    std::memset((void*)nu, 0, sizeof(nu));
    std::memset((void*)chi, 0, sizeof(chi));
    
    for (auto it : this->E) {
        const unsigned i = it.first;
        const unsigned j = it.second;

        const unsigned mu = this->mu[i];
        const unsigned omega = this->omega[j];
        const long double pp = this->pi[j];
        const long double logpp = this->logpi[j];

        const long double mu_pow_phi = std::pow(mu, phi);
        const long double log_mu = std::log(mu);

        omega_phi[j] += mu_pow_phi;
        nu[j] += mu_pow_phi * log_mu;

        seen_r[j]++;
        if (seen_r[j] == omega) { // Last time we will see this j
            long double omega_p = omega_phi[j].value();
            long double nuv = nu[j].value();
            this->omega_phi[j] = omega_p;
            this->nu[j] = nuv;
            
            HSR += pp * (phi * nuv / omega_p + logpp - std::log(omega_p));
            HR += pp * logpp;
            rho += pp;
        }
    }
    
    for (auto it : this->E) {
        const unsigned i = it.first;
        const unsigned j = it.second;

        const unsigned mu = this->mu[i];
        const long double pp = this->pi[j];
        const long double omega_phi = this->omega_phi[j];

        chi[i] += pp / omega_phi;

        seen_s[i]++;
        if (seen_s[i] == mu) { // Last time we will see this i
            const long double ch = chi[i].value();
            this->chi[i] = ch;
            
            const long double p = std::pow(mu, phi) * ch;
            HS += p * std::log(p);
        }
    }

    if (this->E.empty()) {
        this->rho = 0;
    } else {
        this->rho = rho.value();
        this->XS = HS.value();
        this->XR = HR.value();
        this->XSR = HSR.value();
    }

    this->sanity_checks();
}

void SecondModelStaticLexicalMatrix::clear() {
    BaseLexicalMatrix::clear();
    
    std::fill(this->omega_phi.begin(), this->omega_phi.end(), 0);
    std::fill(this->chi.begin(), this->chi.end(), 0);
    std::fill(this->nu.begin(), this->nu.end(), 0);
    
    this->XSR = 0;
    this->XS = 0;
    this->XR = 0;
    this->rho = 0;
}

long double SecondModelStaticLexicalMatrix::min_HS() const {
    return 0;
}
long double SecondModelStaticLexicalMatrix::max_HS() const {
    return std::log(this->n);
}
long double SecondModelStaticLexicalMatrix::min_HR() const {
    return 0;
}
long double SecondModelStaticLexicalMatrix::max_HR() const {
    return this->HR_pi;
}
long double SecondModelStaticLexicalMatrix::min_HSR() const {
    return 0;
}
long double SecondModelStaticLexicalMatrix::max_HSR() const {
    return std::log(this->n) + this->HR_pi;
}

long double SecondModelStaticLexicalMatrix::get_ps(unsigned s) const {
    return pow(this->get_mu(s), this->phi) * this->get_chi(s) / this->get_rho();
}

long double SecondModelStaticLexicalMatrix::get_pr(unsigned r) const {
    return this->meaning_is_connected(r) * this->pi[r] / this->get_rho();
}

long double SecondModelStaticLexicalMatrix::get_psr(unsigned s, unsigned r) const {
    return this->get_a(s, r) * this->meaning_is_connected(r) * pow(this->get_mu(s), this->phi) * this->pi[r] / this->get_rho() / this->get_omega_phi(r);
}

long double SecondModelStaticLexicalMatrix::get_HS() const {
    if (this->E.empty()) {
        return 0;
    } else {
        return std::log(this->rho) - this->XS / this->rho;
    }
}

long double SecondModelStaticLexicalMatrix::get_HR() const {
    if (this->E.empty()) {
        return 0;
    } else {
        return std::log(this->rho) - this->XR / this->rho;
    }
}

long double SecondModelStaticLexicalMatrix::get_HSR() const {
    if (this->E.empty()) {
        return 0;
    } else {
        return std::log(this->rho) - this->XSR / this->rho;
    }
}

/**
 * Print the state of the graph to stderr.
 */
void SecondModelStaticLexicalMatrix::print_properties() const {
    using namespace std;

    BaseLexicalMatrix::print_properties();
    cout << endl;
    
    cout << "SECOND MODEL PROPERTIES" << endl;
    cout << "-----------------------" << endl;

    cout << "pi = ";
    pretty_print_vector(cout, this->pi);
    cout << endl;
    
    cout << "chi = ";
    pretty_print_vector(cout, this->chi);
    cout << endl;
    
    cout << "omega_phi = ";
    pretty_print_vector(cout, this->omega_phi);
    cout << endl;

    cout << "nu = ";
    pretty_print_vector(cout, this->nu);
    cout << endl;
    
    cout << "rho = " << this->rho << endl;
    cout << "X(S,R) = " << this->XSR << endl;
    cout << "X(S) = " << this->XS << endl;
    cout << "X(R) = " << this->XR << endl;
}
