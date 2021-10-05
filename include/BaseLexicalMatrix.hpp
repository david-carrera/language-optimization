#ifndef BASE_LEXICAL_MATRIX_H
#define BASE_LEXICAL_MATRIX_H

#include <Util.hpp>
#include <LexicalMatrixParameters.hpp>
#include <vector>
#include <set>

/**
 * Base class for all Lexical Matrixes.
 */
class BaseLexicalMatrix {
private:
    void color_connected_nodes(unsigned node, bool is_word, unsigned color, std::vector<unsigned> &colorWords, std::vector<unsigned> &colorMeanings) const;

protected:
    /// Adjacency matrix, indexed as A[r][s] (where r is a meaning and s is a word)
    std::vector<std::vector<bool>> A;
    
    /// Set of edges in the graph
    std::set<std::pair<unsigned, unsigned>> E;

    /// mu of each word
    std::vector<unsigned> mu;

    /// omega of each meaning
    std::vector<unsigned> omega;

    /// Keep track of which mutations have been done
    std::vector<std::pair<unsigned, unsigned>> mutations;

    /// Whether to record mutations
    bool recording;

    /// Keep track of whether unlinked meanings exist, updated on matrix mutation.
    unsigned unlinked_meanings;
    
    virtual void add_edge_never_recalculate(unsigned s, unsigned r);
    
public:
    /// Number of words
    const unsigned n;
    
    /// Number of meanings
    const unsigned m;
    
    /// The phi parameter
    const double phi;
    
    /// Error tolerance
    const long double epsilon;

    // minimum possible value of H(S) for this model (for sanity tests)
    virtual long double min_HS() const = 0;
    
    // maximum possible value of H(S) for this model (for sanity tests)
    virtual long double max_HS() const = 0;

    // minimum possible value of H(R) for this model (for sanity tests)
    virtual long double min_HR() const = 0;
    
    // maximum possible value of H(R) for this model (for sanity tests)
    virtual long double max_HR() const = 0;
    
    // minimum possible value of H(S,R) for this model (for sanity tests)
    virtual long double min_HSR() const = 0;
    
    // maximum possible value of H(S,R) for this model (for sanity tests)
    virtual long double max_HSR() const = 0;

    /// minimum possible value of I(S,R) for this model
    virtual long double min_ISR() const {
        return this->min_HS() + this->min_HR() - this->max_HSR();
    }

    /// maximum possible value of I(S,R) for this model
    virtual long double max_ISR() const {
        return this->max_HS() + this->max_HR() - this->min_HSR();
    }

    BaseLexicalMatrix() = delete;
    BaseLexicalMatrix(const BaseLexicalMatrix&) = delete;
    BaseLexicalMatrix &operator=(const BaseLexicalMatrix &other) = delete;
    virtual ~BaseLexicalMatrix() {}

    BaseLexicalMatrix(unsigned n, unsigned m, double phi, long double epsilon);    
    void add_edge(unsigned s, unsigned r);
    void remove_edge(unsigned s, unsigned r);
    void account_unlinked_meanings(unsigned s, unsigned r);
    virtual void mutate(unsigned s, unsigned r);
    
    /**
     * Statically recalculate all probabilities and entropies.
     * @param[in]  force  If true, perform the recalculation even if dynamic.
     */
    virtual void recalculate(bool force=false) = 0;

    virtual void record(bool record);
    virtual void save();
    virtual void restore();
    virtual void clear();

    void initialize_random(RNG &rng, double p, bool need_connected_meanings = false);
    void initialize_random_edges(RNG &rng, unsigned edges, bool need_connected_meanings = false);
    void initialize_complete_graph();
    void initialize_bijective_graph();

    unsigned color_connected_nodes(std::vector<unsigned> &colorWords, std::vector<unsigned> &colorMeanings) const;

    /**
     * Is the matrix in a valid state? A matrix is always in an invalid state
     * when it has no edges. A matrix may be in an invalid state when one of
     * the meanings has no connections.
     * @return True when the matrix is in a valid state.
     */
    inline bool valid(bool allow_unlinked_meanings = true) const {
        return !this->E.empty() &&
            (allow_unlinked_meanings || this->unlinked_meanings == 0);
    }

    /**
     * Are a word and a meaning connected?
     * @param[in]  s  Word
     * @param[in]  r  Meaning
     * @return        Are word and meaning connected?
     */
    inline bool get_a(unsigned s, unsigned r) const {
        return this->A[r][s];
    }

    /**
     * mu of a word
     * @param[in]  i  Word
     * @return        mu of the word
     */
    inline unsigned get_mu(unsigned i) const {
        return this->mu[i];
    }

    /**
     * omega of a meaning
     * @param[in]  j  Meaning
     * @return        omega of the meaning
     */
    inline unsigned get_omega(unsigned j) const {
        return this->omega[j];
    }

    /**
     * Probability of a word
     * @param[in]  s  Word
     * @return        Probability
     */
    virtual long double get_ps(unsigned s) const = 0;

    /**
     * Probability of a meaning
     * @param[in]  r  Meaning
     * @return        Probability
     */
    virtual long double get_pr(unsigned r) const = 0;

    /**
     * Joint probability of a word and a meaning
     * @param[in]  s  Word
     * @param[in]  r  Meaning
     * @return        Probability
     */
    virtual long double get_psr(unsigned s, unsigned r) const = 0;
    
    /// H(S), entropy of words
    virtual long double get_HS() const = 0;

    /// H(R), entropy of meanings
    virtual long double get_HR() const = 0;

    /// H(S,R), joint entropy of words and meanings
    virtual long double get_HSR() const = 0;

    /// I(S,R), mutual information of words and meanings
    virtual long double get_ISR() const {
        return this->get_HS() + this->get_HR() - this->get_HSR();
    }

    /// H(S|R), conditional entropy of words given meanings
    virtual long double get_HSgR() const {
        return this->get_HSR() - this->get_HR();
    }
    
    /// H(R|S), conditional entropy of meanings given words
    virtual long double get_HRgS() const {
        return this->get_HSR() - this->get_HS();
    }

    std::pair<unsigned, unsigned> largest_connected_component() const;
    unsigned referentially_useless_words() const;

    void get_parameters(LexicalMatrixParameters<long double> &parameters) const;
    void sanity_checks() const;

    virtual void print_properties() const;
};

#endif /* BASE_LEXICAL_MATRIX_H */
