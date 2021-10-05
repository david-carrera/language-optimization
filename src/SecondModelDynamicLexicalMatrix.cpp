#include <SecondModelDynamicLexicalMatrix.hpp>
#include <algorithm>
#include <cassert>

inline long double SecondModelDynamicLexicalMatrix::calculate_xs(unsigned mu, long double chi) const {
    long double x = std::pow(mu, phi) * chi;
    if (std::abs(x) < std::numeric_limits<long double>::epsilon()) {
        return 0;
    } else {
        return x * std::log(x);
    }
}

inline long double SecondModelDynamicLexicalMatrix::calculate_xs_phi0(long double chi) const {
    if (std::abs(chi) < std::numeric_limits<long double>::epsilon()) {
        return 0;
    } else {
        return chi * std::log(chi);
    }
}

inline long double SecondModelDynamicLexicalMatrix::calculate_xr(long double omega_phi, long double nu, unsigned r) const {
    long double pi = this->pi[r];
    long double logpi = this->logpi[r];
    long double phi = this->phi;
    
    return pi * (phi * nu / omega_phi + logpi - std::log(omega_phi));
}

inline long double SecondModelDynamicLexicalMatrix::calculate_xr_phi0(unsigned omega, unsigned j) const {
    long double pi = this->pi[j];
    long double logpi = this->logpi[j];
    return pi * (logpi - std::log(omega));
}

SecondModelDynamicLexicalMatrix::SecondModelDynamicLexicalMatrix(unsigned n, unsigned m, const std::vector<long double> &pi, double phi, long double epsilon) :
    SecondModelStaticLexicalMatrix(n, m, pi, phi, epsilon), phi0(std::abs(phi) < std::numeric_limits<double>::epsilon()), neighbors_word(n), neighbors_meaning(m) {}

void SecondModelDynamicLexicalMatrix::recalculate(bool force) {
    if (force) {
        SecondModelStaticLexicalMatrix::recalculate();
    }
}

/**
 * @copydoc BaseLexicalMatrix::mutate
 * Parameters are recalculated dynamically.
 */
void SecondModelDynamicLexicalMatrix::mutate(unsigned int s, unsigned int r) {
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

    // Update neighborhood
    if (!adding_edge) {
        this->neighbors_word[s].erase(r);
        this->neighbors_meaning[r].erase(s);
    }

    bool meaning_becomes_connected = omega_r == 0;
    bool meaning_becomes_disconnected = new_omega_r == 0;
    
    bool meaning_is_connected = new_omega_r != 0;
    bool meaning_was_connected = omega_r != 0;

    // Update omega_phi and nu    
    long double old_omega_phi[this->m];
    long double old_nu[this->m];
    {
        long double mu_s_phi = std::pow(mu_s, this->phi);
        long double new_mu_s_phi = std::pow(new_mu_s, this->phi);

        long double mu_s_phi_log = mu_s_phi * std::log(mu_s);
        long double new_mu_s_phi_log = new_mu_s_phi * std::log(new_mu_s);
        
        for (unsigned l : this->neighbors_word[s]) {
            old_omega_phi[l] = this->omega_phi[l];
            this->omega_phi[l] -= mu_s_phi;
            this->omega_phi[l] += new_mu_s_phi;
            this->changed_omega_phi.insert(l);
            
            old_nu[l] = this->nu[l];
            this->nu[l] -= mu_s_phi_log;
            this->nu[l] += new_mu_s_phi_log;
            this->changed_nu.insert(l);
        }

        old_omega_phi[r] = this->omega_phi[r];
        old_nu[r] = this->nu[r];
        if (adding_edge) {
            this->omega_phi[r] += new_mu_s_phi;
            this->nu[r] += new_mu_s_phi_log;
        } else {
            this->omega_phi[r] -= mu_s_phi;
            this->nu[r] -= mu_s_phi_log;
        }
        this->changed_omega_phi.insert(r);
        this->changed_nu.insert(r);
    }

    // Update chi
    std::set<unsigned> updated_chi;
    long double old_chi[this->n];
    {
        std::set<unsigned> indices_changed_omega_phi = this->neighbors_word[s];
        indices_changed_omega_phi.insert(r);
        
        // Update all chi except s
        for (unsigned k=0; k<this->n; k++) {
            if (k == s) {
                continue;
            }

            std::vector<unsigned> neighbors_to_update;
            std::set_intersection(indices_changed_omega_phi.begin(),
                                  indices_changed_omega_phi.end(),
                                  this->neighbors_word[k].begin(),
                                  this->neighbors_word[k].end(),
                                  std::inserter(
                                      neighbors_to_update,
                                      neighbors_to_update.begin()));

            if (neighbors_to_update.size() > 0) {
                old_chi[k] = this->chi[k];
                for (unsigned l : neighbors_to_update) {
                    this->chi[k] -= this->pi[l] / old_omega_phi[l];
                    this->chi[k] += this->pi[l] / this->omega_phi[l];
                }
                this->changed_chi.insert(k);
                updated_chi.insert(k);
            }
        }

        // Update chi s statically
        old_chi[s] = this->chi[s];
        Accumulator chi_s;
        // add previous neighbors
        for (unsigned l : this->neighbors_word[s]) {
            chi_s += this->pi[l] / this->omega_phi[l];
        }
        // add the new neighbor if we're adding an edge
        if (adding_edge) {
            chi_s += this->pi[r] / this->omega_phi[r];
        }
        this->chi[s] = chi_s.value();
        // do not add to updated_chi, we already always account for it
        // do add to changed_chi for save/restore
        this->changed_chi.insert(s);
    }

    // Update XR and rho
    if (meaning_becomes_connected) {
        long double pi = this->pi[r];
        long double logpi = this->logpi[r];
        long double pi_log_pi = pi * logpi;
        this->rho += pi;
        this->XR += pi_log_pi;
    } else if (meaning_becomes_disconnected) {
        long double pi = this->pi[r];
        long double logpi = this->logpi[r];
        long double pi_log_pi = pi * logpi;
        this->rho -= pi;
        this->XR -= pi_log_pi;
    }

    // Update XSR
    for (unsigned l : this->neighbors_word[s]) {
        long double o_nu = old_nu[l];
        long double o_omega_phi = old_omega_phi[l];
        
        long double n_nu = this->nu[l];
        long double n_omega_phi = this->omega_phi[l];
        
        this->XSR -= this->calculate_xr(o_omega_phi, o_nu, l);
        this->XSR += this->calculate_xr(n_omega_phi, n_nu, l);
    }
    if (meaning_is_connected) {
        this->XSR += this->calculate_xr(this->omega_phi[r], this->nu[r], r);
    }
    // not else if, both need to happen in most cases
    if (meaning_was_connected) {
        this->XSR -= this->calculate_xr(old_omega_phi[r], old_nu[r], r);
    }

    // Update XS
    for (unsigned k : updated_chi) {
        this->XS -= this->calculate_xs(this->mu[k], old_chi[k]);
        this->XS += this->calculate_xs(this->mu[k], this->chi[k]);
    }
    this->XS -= this->calculate_xs(mu_s, old_chi[s]);
    this->XS += this->calculate_xs(new_mu_s, this->chi[s]);

    // Update neighborhood
    if (adding_edge) {
        this->neighbors_word[s].insert(r);
        this->neighbors_meaning[r].insert(s);
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
void SecondModelDynamicLexicalMatrix::mutate_phi0(unsigned int i, unsigned int j) {
    //unsigned mu_i = this->mu[i];
    unsigned omega_j = this->omega[j];

    unsigned new_omega_j;
    bool adding_edge;
    if (this->get_a(i, j)) { // Remove edge
        adding_edge = false;
        assert(this->omega[j] > 0);
        
        --this->mu[i];
        new_omega_j = --this->omega[j];
        
        this->A[j][i] = false;
        this->E.erase({i, j});
    } else { // Add edge
        adding_edge = true;
        
        ++this->mu[i];
        new_omega_j = ++this->omega[j];
        
        this->A[j][i] = true;
        this->E.insert({i, j});
    }

    if (!adding_edge) {
        this->neighbors_word[i].erase(j);
        this->neighbors_meaning[j].erase(i);
    }

    bool meaning_becomes_connected = omega_j == 0;
    bool meaning_becomes_disconnected = new_omega_j == 0;
    
    bool meaning_is_connected = new_omega_j != 0;
    bool meaning_was_connected = omega_j != 0;

    // Update chi
    long double old_chi[this->n];
    {
        long double old_element = this->pi[j] / omega_j;
        long double new_element = this->pi[j] / new_omega_j;

        if (!meaning_is_connected) {
            new_element = 0;
        }
        
        for (unsigned k : this->neighbors_meaning[j]) {
            old_chi[k] = this->chi[k];
            if (this->mu[k] == 0) { // TODO: can we delete this and still pass tests?
                this->chi[k] = 0;
            } else {
                this->chi[k] -= old_element;
                this->chi[k] += new_element;
            }
            this->changed_chi.insert(k);
        }
        
        old_chi[i] = this->chi[i];
        if (this->mu[i] == 0) { // TODO: can we delete this and still pass tests?
            this->chi[i] = 0;
        } else {
            if (adding_edge) {
                this->chi[i] += new_element;
            } else {
                this->chi[i] -= old_element;
            }
        }
        this->changed_chi.insert(i);
    }
    
    // Update XR and rho
    if (meaning_becomes_connected) {
        long double pi = this->pi[j];
        long double logpi = this->logpi[j];
        long double pi_log_pi = pi * logpi;
        this->rho += pi;
        this->XR += pi_log_pi;
    } else if (meaning_becomes_disconnected) {
        long double pi = this->pi[j];
        long double logpi = this->logpi[j];
        long double pi_log_pi = pi * logpi;
        this->rho -= pi;
        this->XR -= pi_log_pi;
    }

    // Update XSR
    if (meaning_is_connected) {
        this->XSR += this->calculate_xr_phi0(new_omega_j, j);
    }
    // not else if, both need to happen in most cases
    if (meaning_was_connected) {
        this->XSR -= this->calculate_xr_phi0(omega_j, j);
    }

    // Update XS
    for (unsigned k : this->neighbors_meaning[j]) {
        this->XS -= this->calculate_xs_phi0(old_chi[k]);
        this->XS += this->calculate_xs_phi0(this->chi[k]);
    }
    this->XS -= this->calculate_xs_phi0(old_chi[i]);
    this->XS += this->calculate_xs_phi0(this->chi[i]);

    // Update neighborhood
    if (adding_edge) {
        this->neighbors_word[i].insert(j);
        this->neighbors_meaning[j].insert(i);
    }

    if (this->E.empty()) {
        this->clear();
    }

    this->sanity_checks();
    this->account_unlinked_meanings(i, j);
}

long double SecondModelDynamicLexicalMatrix::get_ps(unsigned s) const {
    if (this->phi0) {
        return this->chi[s] / this->rho;
    } else {
        return SecondModelStaticLexicalMatrix::get_ps(s);
    }
}

void SecondModelDynamicLexicalMatrix::save() {
    if (this->recording) {
        this->mutations.clear();
        this->changed_omega_phi.clear();
        this->changed_nu.clear();
        this->changed_chi.clear();

        this->saved_XSR = this->XSR;
        this->saved_XS = this->XS;
        this->saved_XR = this->XR;
        this->saved_rho = this->rho;

        this->saved_omega_phi = this->omega_phi;
        this->saved_nu = this->nu;
        this->saved_chi = this->chi;

        this->saved_unlinked_meanings = this->unlinked_meanings;
    }
}

void SecondModelDynamicLexicalMatrix::restore() {
    if (this->recording) {
        this->unlinked_meanings = this->saved_unlinked_meanings;
        
        this->XSR = this->saved_XSR;
        this->XS = this->saved_XS;
        this->XR = this->saved_XR;
        this->rho = this->saved_rho;

        if (!this->phi0) {
            for (unsigned j : this->changed_omega_phi) {
                this->omega_phi[j] = this->saved_omega_phi[j];
            }
            for (unsigned j : this->changed_nu) {
                this->nu[j] = this->saved_nu[j];
            }
        }
        for (unsigned i : this->changed_chi) {
            this->chi[i] = this->saved_chi[i];
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
                this->neighbors_word[i].erase(j);
                this->neighbors_meaning[j].erase(i);
            } else {
                this->neighbors_word[i].insert(j);
                this->neighbors_meaning[j].insert(i);
                this->mu[i]++;
                this->omega[j]++;
                this->A[j][i] = true;\
                this->E.insert({i, j});
            }
        }
        this->recording = true;

        this->mutations.clear();
        this->changed_omega_phi.clear();
        this->changed_nu.clear();
        this->changed_chi.clear();
    }
}

void SecondModelDynamicLexicalMatrix::clear() {
    for (auto &set : this->neighbors_word) {
        set.clear();
    }
    for (auto &set : this->neighbors_meaning) {
        set.clear();
    }
    SecondModelStaticLexicalMatrix::clear();
}

/**
 * Add an edge between a word and a meaning without doing any recalculation.
 * @param[in]  s  Word
 * @param[in]  r  Meaning
 */
void SecondModelDynamicLexicalMatrix::add_edge_never_recalculate(unsigned s, unsigned r) {
    bool aux = this->recording;
    this->recording = false;
    assert(!this->get_a(s, r));
    this->neighbors_word[s].insert(r);
    this->neighbors_meaning[r].insert(s);
    SecondModelStaticLexicalMatrix::mutate(s, r);
    this->recording = aux;
}
