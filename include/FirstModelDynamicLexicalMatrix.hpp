#ifndef FIRST_MODEL_DYNAMIC_LEXICAL_MATRIX_H
#define FIRST_MODEL_DYNAMIC_LEXICAL_MATRIX_H

#include <FirstModelStaticLexicalMatrix.hpp>

/**
 * A lexical matrix whose values can be computed dynamically or statically. It
 * uses the first model to compute probabilities and entropies.
 */
class FirstModelDynamicLexicalMatrix : public FirstModelStaticLexicalMatrix {
private:
    void update_dynamic_XSR_M(unsigned old_mu, unsigned old_omega, unsigned new_mu, unsigned new_omega);
    void update_dynamic_XS(unsigned i, unsigned old_mu_i, unsigned new_mu_i, unsigned old_omega_l, unsigned new_omega_l);
    void update_dynamic_XR(unsigned j, unsigned old_omega_j, unsigned new_omega_j, unsigned old_mu_k, unsigned new_mu_k);

    /// save the counter of unlinked meanings for restoring
    unsigned saved_unlinked_meanings;

    /// saved mu_phi for restoring
    std::vector<long double> saved_mu_phi;
    /// saved omega_phi for restoring
    std::vector<long double> saved_omega_phi;
    /// which mu_phi changed since save
    std::set<unsigned> changed_mu_phi;
    /// which omega_phi changed since save
    std::set<unsigned> changed_omega_phi;

    /// saved X(S,R) for restoring
    long double saved_XSR;
    /// saved X(S) for restoring
    long double saved_XS;
    /// saved X(R) for restoring
    long double saved_XR;
    /// saved M for restoring
    long double saved_M;

    /// True if the phi of this matrix is 0
    const bool phi0;
protected:
    /// Neighbors of every word.
    std::vector<std::set<unsigned>> neighbors_word;
    /// Neighbors of every meaning.
    std::vector<std::set<unsigned>> neighbors_meaning;

    void mutate_phi0(unsigned s, unsigned r);
    virtual void add_edge_never_recalculate(unsigned s, unsigned r);
public:
    FirstModelDynamicLexicalMatrix(unsigned n, unsigned m, double phi = 0.0, long double epsilon = default_epsilon);
    
    virtual void mutate(unsigned s, unsigned r);
    virtual void recalculate(bool force=false);

    virtual long double get_ps(unsigned s) const;
    virtual long double get_pr(unsigned r) const;

    virtual void save();
    virtual void restore();
    virtual void clear();
};

#endif /* DYNAMIC_LEXICAL_MATRIX_H */
