#ifndef FIRST_MODEL_STATIC_LEXICAL_MATRIX_H
#define FIRST_MODEL_STATIC_LEXICAL_MATRIX_H

#include <BaseLexicalMatrix.hpp>
#include <LexicalMatrixParameters.hpp>

/**
 * A lexical matrix whose values are computed statically. It uses the first
 * model to compute probabilities and entropies.
 */
class FirstModelStaticLexicalMatrix : public BaseLexicalMatrix {
protected:
    /// mu_phi of each word
    std::vector<long double> mu_phi;
    /// omega_phi of each meaning
    std::vector<long double> omega_phi;

    /// M_phi
    long double M;
    /// X(S,R)
    long double XSR;
    /// X(S)
    long double XS;
    /// X(R)
    long double XR;

public:
    FirstModelStaticLexicalMatrix(unsigned n, unsigned m, double phi = 0.0, long double epsilon = default_epsilon);

    virtual long double min_HS() const;
    virtual long double max_HS() const;
    virtual long double min_HR() const;
    virtual long double max_HR() const;
    virtual long double min_HSR() const;
    virtual long double max_HSR() const;
    
    virtual void recalculate(bool force=false); 
    virtual void clear();
    
    /**
     * mu_phi of a word
     * @param[in]  i  Word
     * @return        mu_phi of the word
     */
    inline long double get_mu_phi(unsigned i) const {
        return this->mu_phi[i];
    }
    
    /**
     * omega_phi of a meaning
     * @param[in]  j  Meaning
     * @return        omega_phi of the meaning
     */
    inline long double get_omega_phi(unsigned j) const {
        return this->omega_phi[j];
    }

    /// M_phi of the graph
    inline long double get_M() const {
        return this->M;
    }

    virtual long double get_ps(unsigned s) const;
    virtual long double get_pr(unsigned r) const;
    virtual long double get_psr(unsigned s, unsigned r) const;
    virtual long double get_HS() const;
    virtual long double get_HR() const;
    virtual long double get_HSR() const;

    inline long double get_XS() const {
        return this->XS;
    }
    inline long double get_XR() const {
        return this->XR;
    }
    inline long double get_XSR() const {
        return this->XSR;
    }

    virtual void print_properties() const;
};

#endif /* STATIC_LEXICAL_MATRIX_H */
