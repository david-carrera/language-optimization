#include <FirstModelDynamicLexicalMatrix.hpp>
#include <algorithm>
#include <cassert>

// TODO: experiment with doing static recalculation if the number of neighbors is too hight (try more than n/2 or m/2)

/**
 * Update X(S,R) and M_phi from older and newer values of omega.
 * @param[in]  old_mu     Previous value of mu
 * @param[in]  old_omega  Previous value of omega
 * @param[in]  new_mu     Current value of mu
 * @param[in]  new_omega  Current value of omega
 */
inline void FirstModelDynamicLexicalMatrix::update_dynamic_XSR_M(unsigned old_mu, unsigned old_omega, unsigned new_mu, unsigned new_omega) {
    long double phi = this->phi;

    unsigned old_qsr_nophi = old_mu * old_omega;
    unsigned new_qsr_nophi = new_mu * new_omega;

    long double old_qsr = std::pow(old_qsr_nophi, phi);
    long double new_qsr = std::pow(new_qsr_nophi, phi);

    if (old_qsr_nophi > 0) {
        this->XSR -= old_qsr * std::log(old_qsr_nophi);
        this->M -= old_qsr;
    }
    if (new_qsr_nophi > 0) {
        this->XSR += new_qsr * std::log(new_qsr_nophi);
        this->M += new_qsr;
    }
}

/**
 * Update X(S) and mu_phi from older and newer values of omega.
 * @param[in]  i            Word to update
 * @param[in]  old_mu_i     Previous value of mu for the word
 * @param[in]  new_mu_i     Current value of mu for the word
 * @param[in]  old_omega_l  Previous value of omega for the mutated meaning
 * @param[in]  new_omega_l  Current value of omega for the mutated meaning
 */
inline void FirstModelDynamicLexicalMatrix::update_dynamic_XS(unsigned i, unsigned old_mu_i, unsigned new_mu_i, unsigned old_omega_l, unsigned new_omega_l) {
    long double phi = this->phi;
    
    if (this->recording) {
        this->changed_mu_phi.insert(i);
    }

    long double old_mu_phi = this->mu_phi[i];
    if (old_omega_l > 0) {
        this->mu_phi[i] -= std::pow(old_omega_l, phi);
    }
    if (new_omega_l > 0) {
        this->mu_phi[i] += std::pow(new_omega_l, phi);
    }

    long double old_qs = std::pow(old_mu_i, phi) * old_mu_phi;
    long double new_qs = std::pow(new_mu_i, phi) * this->mu_phi[i];

    if (old_qs > 0) {
        this->XS -= old_qs * std::log(old_qs);
    }
    
    if (new_qs > 0) {
        this->XS += new_qs * std::log(new_qs);
    }
}

/**
 * Update X(R) and omega_phi from older and newer values of omega.
 * @param[in]  j            Meaning to update
 * @param[in]  old_omega_j  Previous value of omega for the meaning
 * @param[in]  new_omega_j  Current value of omega for the meaning
 * @param[in]  old_mu_k     Previous value of mu for the mutated word
 * @param[in]  new_mu_k     Current value of mu for the mutated word
 */
inline void FirstModelDynamicLexicalMatrix::update_dynamic_XR(unsigned j, unsigned old_omega_j, unsigned new_omega_j, unsigned old_mu_k, unsigned new_mu_k) {
    long double phi = this->phi;
    
    if (this->recording) {
        this->changed_omega_phi.insert(j);
    }

    long double old_omega_phi = this->omega_phi[j];
    if (old_mu_k > 0) {
        this->omega_phi[j] -= std::pow(old_mu_k, phi);
    }
    if (new_mu_k > 0) {
        this->omega_phi[j] += std::pow(new_mu_k, phi);
    }
    
    long double old_qr = std::pow(old_omega_j, phi) * old_omega_phi;
    long double new_qr = std::pow(new_omega_j, phi) * this->omega_phi[j];

    if (old_qr > 0) {
        this->XR -= old_qr * std::log(old_qr);
    }
    if (new_qr > 0) {
        this->XR += new_qr * std::log(new_qr);
    }
}

/**
 * Constructor.
 * @param[in]  n    Number of words
 * @param[in]  m    Number of meanings
 * @param[in]  phi  Phi parameter
 * @param[in]  epsilon  Error tolerance for comparisons of floating point values.
 */
FirstModelDynamicLexicalMatrix::FirstModelDynamicLexicalMatrix(unsigned n, unsigned m, double phi, long double epsilon) :
    FirstModelStaticLexicalMatrix(n, m, phi, epsilon), phi0(std::abs(phi) < std::numeric_limits<double>::epsilon()),
    neighbors_word(n), neighbors_meaning(m) {}

void FirstModelDynamicLexicalMatrix::recalculate(bool force) {
    if (force) {
        FirstModelStaticLexicalMatrix::recalculate();
    }
}

/**
 * @copydoc BaseLexicalMatrix::mutate
 * Parameters are recalculated dynamically.
 */
void FirstModelDynamicLexicalMatrix::mutate(unsigned s, unsigned r) {
    if (this->recording) {
        this->mutations.push_back({s, r});
    }
    
    if (this->phi0) {
        this->mutate_phi0(s, r);
        return;
    }
    
    unsigned mu_s = this->mu[s];
    unsigned omega_r = this->omega[r];

    unsigned new_mu_s, new_omega_r;

    bool adding_edge;
    if (this->get_a(s, r)) { // Remove edge
        adding_edge = false;
        
        new_mu_s = --this->mu[s];
        new_omega_r = --this->omega[r];
        
        this->A[r][s] = false;
        this->E.erase({s, r});
    } else { // Add edge
        adding_edge = true;
        
        new_mu_s = ++this->mu[s];
        new_omega_r = ++this->omega[r];
        
        this->A[r][s] = true;
        this->E.insert({s, r});
    }

    if (!adding_edge) {
        this->neighbors_meaning[r].erase(s);
        this->neighbors_word[s].erase(r);
    }
    
    for (unsigned i : this->neighbors_meaning[r]) {
        unsigned mu_i = this->mu[i];
        this->update_dynamic_XSR_M(mu_i, omega_r, mu_i, new_omega_r);
        this->update_dynamic_XS(i, mu_i, mu_i, omega_r, new_omega_r);
    }

    for (unsigned j : this->neighbors_word[s]) {
        unsigned omega_j = this->omega[j];
        this->update_dynamic_XSR_M(mu_s, omega_j, new_mu_s, omega_j);
        this->update_dynamic_XR(j, omega_j, omega_j, mu_s, new_mu_s);
    }

    if (adding_edge) {
        this->update_dynamic_XSR_M(0, 0, new_mu_s, new_omega_r);
        this->update_dynamic_XS(s, mu_s, new_mu_s, 0, new_omega_r);
        this->update_dynamic_XR(r, omega_r, new_omega_r, 0, new_mu_s);
        this->neighbors_meaning[r].insert(s);
        this->neighbors_word[s].insert(r);
    } else {
        this->update_dynamic_XSR_M(mu_s, omega_r, 0, 0);
        this->update_dynamic_XS(s, mu_s, new_mu_s, omega_r, 0);
        this->update_dynamic_XR(r, omega_r, new_omega_r, mu_s, 0);
    }

    if (this->E.empty()) {
        this->clear();
    }

    this->sanity_checks();
    this->account_unlinked_meanings(s, r);
}

/**
 * Mutate the edge between a word and a meaning assuming the matrix's phi is 0.
 * @param[in]  s  Word
 * @param[in]  r  Meaning
 */
void FirstModelDynamicLexicalMatrix::mutate_phi0(unsigned int i, unsigned int j) {
    // if (this->mu[i] > 0 && this->omega[j] > 0) {
    //     this->XSR -= std::log(this->mu[i] * this->omega[j]);
    // }
    if (this->mu[i] > 0) {
        this->XS -= this->mu[i] * std::log(this->mu[i]);
    }
    if (this->omega[j] > 0) {
        this->XR -= this->omega[j] * std::log(this->omega[j]);
    }

    if (this->get_a(i, j)) { // Remove edge
        this->A[j][i] = false;
        this->E.erase({i, j});
        
        this->mu[i] -= 1;
        this->mu_phi[i] -= 1;
        
        this->omega[j] -= 1;
        this->omega_phi[j] -= 1;
        
        this->M -= 1;
    } else { // Add edge
        this->A[j][i] = true;
        this->E.insert({i, j});

        this->mu[i] += 1;
        this->mu_phi[i] += 1;
        
        this->omega[j] += 1;
        this->omega_phi[j] += 1;
        
        this->M += 1;
    }
        
    // if (this->mu[i] > 0 && this->omega[j] > 0) {
    //     this->XSR += std::log(this->mu[i] * this->omega[j]);
    // }
    if (this->mu[i] > 0) {
        this->XS += this->mu[i] * std::log(this->mu[i]);
    }
    if (this->omega[j] > 0) {
        this->XR += this->omega[j] * std::log(this->omega[j]);
    }
    
    if (this->E.empty()) {
        this->M = 0;
        this->clear();
    }
    
    this->sanity_checks();
    this->account_unlinked_meanings(i, j);
}

long double FirstModelDynamicLexicalMatrix::get_ps(unsigned s) const {
    if (this->phi0) {
        long double mu = this->get_mu(s);
        return std::pow(mu, this->phi) * mu / this->get_M();
    } else {
        return FirstModelStaticLexicalMatrix::get_ps(s);
    }
}

long double FirstModelDynamicLexicalMatrix::get_pr(unsigned r) const {
    if (this->phi0) {
        long double omega = this->get_omega(r);
        return std::pow(omega, this->phi) * omega / this->get_M();
    } else {
        return FirstModelStaticLexicalMatrix::get_pr(r);
    }
}

void FirstModelDynamicLexicalMatrix::save() {
    if (this->recording) {
        this->mutations.clear();
        this->changed_mu_phi.clear();
        this->changed_omega_phi.clear();
        
        this->saved_XSR = this->XSR;
        this->saved_XS = this->XS;
        this->saved_XR = this->XR;
        this->saved_M = this->M;
        
        this->saved_mu_phi = this->mu_phi;
        this->saved_omega_phi = this->omega_phi;

        this->saved_unlinked_meanings = this->unlinked_meanings;
    }
}

void FirstModelDynamicLexicalMatrix::restore() {
    if (this->recording) {
        this->unlinked_meanings = this->saved_unlinked_meanings;
        
        this->XSR = this->saved_XSR;
        this->XS = this->saved_XS;
        this->XR = this->saved_XR;
        this->M = this->saved_M;

        if (!this->phi0) {
            for (unsigned i : this->changed_mu_phi) {
                this->mu_phi[i] = this->saved_mu_phi[i];
            }
            for (unsigned j : this->changed_omega_phi) {
                this->omega_phi[j] = this->saved_omega_phi[j];
            }
        }
            
        this->recording = false;
        for (auto it = this->mutations.rbegin(); it!=this->mutations.rend(); ++it) {
            unsigned i = it->first;
            unsigned j = it->second;
            
            if (this->get_a(i, j)) {
                this->mu[i]--;
                this->omega[j]--;
                this->A[j][i] = false;
                this->E.erase({i, j});

                if (!this->phi0) {
                    this->neighbors_meaning[j].erase(i);
                    this->neighbors_word[i].erase(j);
                }
            } else {
                if (!this->phi0) {
                    this->neighbors_meaning[j].insert(i);
                    this->neighbors_word[i].insert(j);
                }
                
                this->mu[i]++;
                this->omega[j]++;
                this->A[j][i] = true;
                this->E.insert({i, j});
            }
        }
        this->recording = true;
        
        this->mutations.clear();
        this->changed_mu_phi.clear();
        this->changed_omega_phi.clear();
    }
}

void FirstModelDynamicLexicalMatrix::clear() {
    for (auto &set : this->neighbors_word) {
        set.clear();
    }
    for (auto &set : this->neighbors_meaning) {
        set.clear();
    }
    FirstModelStaticLexicalMatrix::clear();
}

/**
 * Add an edge between a word and a meaning without doing any recalculation.
 * @param[in]  s  Word
 * @param[in]  r  Meaning
 */
void FirstModelDynamicLexicalMatrix::add_edge_never_recalculate(unsigned s, unsigned r) {
    bool aux = this->recording;
    this->recording = false;
    assert(!this->get_a(s, r));
    this->neighbors_meaning[r].insert(s);
    this->neighbors_word[s].insert(r);
    FirstModelStaticLexicalMatrix::mutate(s, r);
    this->recording = aux;
}
