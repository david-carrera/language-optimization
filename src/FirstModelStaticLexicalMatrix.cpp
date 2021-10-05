#include <FirstModelStaticLexicalMatrix.hpp>
#include <cstring>

FirstModelStaticLexicalMatrix::FirstModelStaticLexicalMatrix(unsigned n, unsigned m, double phi, long double epsilon) :
    BaseLexicalMatrix(n, m, phi, epsilon),
    mu_phi(n), omega_phi(m),
    M(0.0), XSR(0.0), XS(0.0), XR(0.0) {
}

void FirstModelStaticLexicalMatrix::recalculate(bool force) {
    (void)force;
    
    const long double phi = this->phi;

    Accumulator M, XSR, XS, XR;

    unsigned seen_s[this->n];
    unsigned seen_r[this->m];
    Accumulator mu_phi[this->n];
    Accumulator omega_phi[this->m];

    std::memset(seen_s, 0, sizeof(seen_s));
    std::memset(seen_r, 0, sizeof(seen_r));

    std::memset((void*)mu_phi, 0, sizeof(mu_phi));
    std::memset((void*)omega_phi, 0, sizeof(omega_phi));
    
    std::fill(this->mu_phi.begin(), this->mu_phi.end(), 0.0L);
    std::fill(this->omega_phi.begin(), this->omega_phi.end(), 0.0L);

    for (auto it : this->E) {
        const unsigned i = it.first;
        const unsigned j = it.second;

        const unsigned mu = this->mu[i];
        const unsigned omega = this->omega[j];

        const unsigned qsr_nophi = mu * omega;
        const long double qsr = std::pow(qsr_nophi, phi);
        
        M += qsr;
        XSR += qsr * std::log(qsr_nophi);

        mu_phi[i] += std::pow(omega, phi);
        omega_phi[j] += std::pow(mu, phi);

        seen_s[i]++;
        if (seen_s[i] == mu) { // Last time we will see this i
            this->mu_phi[i] = mu_phi[i].value();
            const long double qs = std::pow(mu, phi) * this->mu_phi[i];
            XS += qs * std::log(qs);
        }

        seen_r[j]++;
        if (seen_r[j] == omega) { // Last time we will see this j
            this->omega_phi[j] = omega_phi[j].value();
            const long double qr = std::pow(omega, phi) * this->omega_phi[j];
            XR += qr * std::log(qr);
        }
    }
    
    if (this->E.empty()) {
        this->M = 0;
    } else {
        this->M = M.value();
        
        this->XSR = XSR.value();
        this->XS = XS.value();
        this->XR = XR.value();
    }
    
    this->sanity_checks();
}


void FirstModelStaticLexicalMatrix::clear() {
    BaseLexicalMatrix::clear();
    
    std::fill(this->mu_phi.begin(), this->mu_phi.end(), 0);
    std::fill(this->omega_phi.begin(), this->omega_phi.end(), 0);
    
    this->XSR = 0;
    this->XS = 0;
    this->XR = 0;
    this->M = 0;
}

long double FirstModelStaticLexicalMatrix::min_HS() const {
    return 0;
}
long double FirstModelStaticLexicalMatrix::max_HS() const {
    return std::log(this->n);
}
long double FirstModelStaticLexicalMatrix::min_HR() const {
    return 0;
}
long double FirstModelStaticLexicalMatrix::max_HR() const {
    return std::log(this->m);
}
long double FirstModelStaticLexicalMatrix::min_HSR() const {
    return 0;
}
long double FirstModelStaticLexicalMatrix::max_HSR() const {
    return log(this->n * this->m);
}

long double FirstModelStaticLexicalMatrix::get_ps(unsigned s) const {
    return std::pow(this->get_mu(s), this->phi) * this->mu_phi[s] / this->get_M();
}

long double FirstModelStaticLexicalMatrix::get_pr(unsigned r) const {
    return std::pow(this->get_omega(r), this->phi) * this->omega_phi[r] / this->get_M();
}

long double FirstModelStaticLexicalMatrix::get_psr(unsigned s, unsigned r) const {
    return std::pow(this->get_mu(s) * this->get_omega(r), this->phi) / this->get_M();
}

long double FirstModelStaticLexicalMatrix::get_HS() const {
    if (this->E.empty()) {
        return 0;
    } else {
        return std::log(this->M) - this->XS / this->M;
    }
}

long double FirstModelStaticLexicalMatrix::get_HR() const {
    if (this->E.empty()) {
        return 0;
    } else {
        return std::log(this->M) - this->XR / this->M;
    }
}

long double FirstModelStaticLexicalMatrix::get_HSR() const {
    if (this->E.empty()) {
        return 0;
    } else {
        return std::log(this->M) - static_cast<long double>(this->phi) * this->XSR / this->M;
    }
}

void FirstModelStaticLexicalMatrix::print_properties() const {
    using namespace std;

    BaseLexicalMatrix::print_properties();
    cout << endl;
    
    cout << "FIRST MODEL PROPERTIES" << endl;
    cout << "----------------------" << endl;
    
    cout << "mu_phi = ";
    pretty_print_vector(cout, this->mu_phi);
    cout << endl;
    
    cout << "omega_phi = ";
    pretty_print_vector(cout, this->omega_phi);
    cout << endl;
    
    cout << "M = " << this->M << endl;
    cout << "X(S,R) = " << this->XSR << endl;
    cout << "X(S) = " << this->XS << endl;
    cout << "X(R) = " << this->XR << endl;
}
