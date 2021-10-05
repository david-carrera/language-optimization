#ifndef SECOND_MODEL_DYNAMIC_LEXICAL_MATRIX_H
#define SECOND_MODEL_DYNAMIC_LEXICAL_MATRIX_H

#include <SecondModelStaticLexicalMatrix.hpp>

/**
 * A lexical matrix whose values can be computed dynamically or statically. It
 * uses the second model to compute probabilities and entropies.
 */
class SecondModelDynamicLexicalMatrix : public SecondModelStaticLexicalMatrix {
private:
    long double calculate_xs(unsigned mu, long double chi) const;
    long double calculate_xr(long double omega_phi, long double nu, unsigned r) const;
    long double calculate_xs_phi0(long double chi) const;
    long double calculate_xr_phi0(unsigned omega, unsigned j) const;

    /// save the counter of unlinked meanings for restoring
    unsigned saved_unlinked_meanings;
    
    /// saved omega_phi for restoring
    std::vector<long double> saved_omega_phi;
    /// saved chi for restoring
    std::vector<long double> saved_chi;
    /// saved nu for restoring
    std::vector<long double> saved_nu;
    /// which omega_phi changed since save
    std::set<unsigned> changed_omega_phi;
    /// which chi changed since save
    std::set<unsigned> changed_chi;
    /// which nu changed since save
    std::set<unsigned> changed_nu;

    /// saved X(S,R) for restoring
    long double saved_XSR;
    /// saved X(S) for restoring
    long double saved_XS;
    /// saved X(R) for restoring
    long double saved_XR;
    /// saved rho for restoring
    long double saved_rho;

    /// True if the phi of this matrix is 0
    const bool phi0;
protected:
    /// Neighbors of every word.
    std::vector<std::set<unsigned>> neighbors_word;
    /// Neighbors of every meaning
    std::vector<std::set<unsigned>> neighbors_meaning;
    
    void mutate_phi0(unsigned s, unsigned r);
    virtual void add_edge_never_recalculate(unsigned s, unsigned r);
public:
    SecondModelDynamicLexicalMatrix(unsigned n, unsigned m, const std::vector<long double> &pi, double phi = 0.0, long double epsilon = default_epsilon);
    virtual void mutate(unsigned s, unsigned r);
    virtual void recalculate(bool force=false);

    virtual long double get_ps(unsigned s) const;

    virtual void save();
    virtual void restore();
    virtual void clear();
};

#endif /* SECOND_MODEL_DYNAMIC_LEXICAL_MATRIX_H */
