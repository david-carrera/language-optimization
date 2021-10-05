#ifndef SECOND_MODEL_STATIC_LEXICAL_MATRIX_H
#define SECOND_MODEL_STATIC_LEXICAL_MATRIX_H

#include <BaseLexicalMatrix.hpp>
#include <LexicalMatrixParameters.hpp>

/**
 * A lexical matrix whose values are computed statically. It uses the second
 * model to compute probabilities and entropies.
 */
class SecondModelStaticLexicalMatrix : public BaseLexicalMatrix {
protected:
    /// omega_phi of each meaning
    std::vector<long double> omega_phi;
    /// chi of each word
    std::vector<long double> chi;
    /// nu of each meaning
    std::vector<long double> nu;
    
    /// rho
    long double rho;
    /// X(S)
    long double XS;
    /// X(R)
    long double XR;
    /// X(S,R)
    long double XSR;

    /**
     * is meaning j connected to any word? equivalent to (1 - delta)
     * @param[in]  j  Meaning
     * @return        true if the meaning is connected, false otherwise
     */
    inline bool meaning_is_connected(unsigned j) const {
        return this->get_omega(j) != 0;
    }
public:
    /// A priori probabilities of meanings
    const std::vector<long double> pi;
    /// Logarithms of the a priori probabilities
    const std::vector<long double> logpi;
    
    /// A priori entropy of meanings
    const long double HR_pi;

    virtual long double min_HS() const;
    virtual long double max_HS() const;
    virtual long double min_HR() const;
    virtual long double max_HR() const;
    virtual long double min_HSR() const;
    virtual long double max_HSR() const;
    
    SecondModelStaticLexicalMatrix(unsigned n, unsigned m,
                                   const std::vector<long double> &pi,
                                   double phi = 0.0, long double epsilon = default_epsilon);
    
    virtual void recalculate(bool force=false); 
    virtual void clear();

    /**
     * omega_phi of a meaning
     * @param[in]  j  Meaning
     * @return        omega_phi of the meaning
     */
    inline long double get_omega_phi(unsigned j) const {
        return this->omega_phi[j];
    }

    /**
     * nu of a meaning
     * @param[in]  j  Meaning
     * @return        nu of the word
     */
    inline long double get_nu(unsigned j) const {
        return this->nu[j];
    }

    /**
     * chi of a word
     * @param[in]  i  Word
     * @return        chi of the word
     */
    inline long double get_chi(unsigned i) const {
        return this->chi[i];
    }

    /// rho of the graph
    inline long double get_rho() const {
        return this->rho;
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

#endif /* SECOND_MODEL_STATIC_LEXICAL_MATRIX_H */
