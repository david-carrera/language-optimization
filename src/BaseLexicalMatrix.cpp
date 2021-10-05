#include <BaseLexicalMatrix.hpp>

#include <stack>
#include <algorithm>
#include <cassert>

/**
 * Assign the given \p color to the given \p node and to all connected nodes.
 * @param[in]   node           node from which to start coloring
 * @param[in]   is_word        whether node is a word or a meaning
 * @param[in]   color          which color to assign
 * @param[out]  colorWords     assign \p color to \p colorWords[i] if word \a i is connected to meaning \p node
 * @param[out]  colorMeanings  assign \p color to \p colorMeanings[j] if meaning \a j is connected to word \p node
 */
void BaseLexicalMatrix::color_connected_nodes(unsigned node, bool is_word, unsigned color, std::vector<unsigned> &colorWords, std::vector<unsigned> &colorMeanings) const {
    std::stack<std::pair<unsigned, bool>> nodes;
    nodes.push({node, is_word});
            
    while (!nodes.empty()) {
        std::pair<unsigned, bool> n = nodes.top();
        nodes.pop();
        unsigned node = n.first;
        bool is_word = n.second;

        if (is_word) {
            colorWords[node] = color;
            unsigned i = node;
            for (unsigned j=0; j<this->m; j++) {
                if (this->get_a(i, j) && colorMeanings[j] == 0) {
                    nodes.push({j, false});
                }
            }
        } else {
            colorMeanings[node] = color;
            unsigned j = node;
            for (unsigned i=0; i<this->n; i++) {
                if (this->get_a(i, j) && colorWords[i] == 0) {
                    nodes.push({i, true});
                }
            }
        }
    }
}

/**
 * Constructor.
 * @param[in]  n    Number of words
 * @param[in]  m    Number of meanings
 * @param[in]  phi  Phi parameter
 * @param[in]  epsilon  Error tolerance for comparisons of floating point values.
 */
BaseLexicalMatrix::BaseLexicalMatrix(unsigned n, unsigned m, double phi, long double epsilon) :
    A(m, std::vector<bool>(n)), mu(n), omega(m), recording(false),
    unlinked_meanings(m),
    n(n), m(m), phi(phi), epsilon(epsilon) {
}

/**
 * Add an edge between a word and a meaning to the graph. Assert that edge does not exist and call mutate.
 * @param[in]  s  Word
 * @param[in]  r  Meaning
 */
void BaseLexicalMatrix::add_edge(unsigned s, unsigned r) {
    assert(!this->get_a(s, r));
    this->mutate(s, r);
}

/**
 * Add an edge between a word and a meaning. This variant should never do any
 * recalculation. This method is useful as small optimization when initializing
 * the matrix. By default, simply calls add_edge.
 */
void BaseLexicalMatrix::add_edge_never_recalculate(unsigned s, unsigned r) {
    this->add_edge(s, r);
}

/**
 * Remove the edge between a word and a meaning from the graph. Assert that edge exists and call mutate.
 * @param[in]  s  Word
 * @param[in]  r  Meaning
 */
void BaseLexicalMatrix::remove_edge(unsigned s, unsigned r) {
    assert(this->get_a(s, r));
    this->mutate(s, r);
}

/**
 * Add edge between a word and a meaning if it doesn't exist, remove it if it does exist.
 * @param[in]  s  Word
 * @param[in]  r  Meaning
 */
void BaseLexicalMatrix::mutate(unsigned s, unsigned r) {
    if (this->get_a(s, r)) {
        this->mu[s]--;
        this->omega[r]--;
        this->A[r][s] = false;
        this->E.erase({s, r});
    } else {
        this->mu[s]++;
        this->omega[r]++;
        this->A[r][s] = true;
        this->E.insert({s, r});
    }
    if (this->recording) {
        this->mutations.push_back({s, r});
    }

    this->account_unlinked_meanings(s, r);
}

/**
 * Called automatically after mutate() to update the counter of unlinked
 * meanings in the matrix.
 */
void BaseLexicalMatrix::account_unlinked_meanings(unsigned s, unsigned r) {
    if (this->get_a(s, r)) { // link was added
        if (this->omega[r] == 1) {
            this->unlinked_meanings--;
        }
    } else { // link was removed
        if (this->omega[r] == 0) {
            this->unlinked_meanings++;
        }
    }
}

/**
 * Whether to record data in order to be able to save/restore the matrix. Starts disabled.
 * @param[in]  record  Whether to enabled recording.
 */
void BaseLexicalMatrix::record(bool record) {
    this->recording = record;
}

/**
 * Save the state of the graph. A subsequent call to #restore will set the graph to the current state.
 */
void BaseLexicalMatrix::save() {
    if (this->recording) {
        this->mutations.clear();
    }
}

/**
 * Restore the state of the graph to what it was on the latest call to #save.
 */
void BaseLexicalMatrix::restore() {
    if (this->recording) {
        this->recording = false;
        for (auto it = this->mutations.rbegin(); it!=this->mutations.rend(); ++it) {
            this->mutate(it->first, it->second);
        }
        this->recording = true;
        this->mutations.clear();
    }
}

/**
 * Remove all edges from the graph.
 */
void BaseLexicalMatrix::clear() {
    for (auto &row : this->A) {
        std::fill(row.begin(), row.end(), false);
    }
    
    std::fill(this->mu.begin(), this->mu.end(), 0);
    std::fill(this->omega.begin(), this->omega.end(), 0);
    
    this->E.clear();
    this->unlinked_meanings = this->m;
}

/**
 * Initialize graph by randomly connecting words and meanings with given probability.
 * For efficiency, this is calculated by generating the number of edges e as a binomially distributed number with parameters n*m and p, then calling initialize_random_edges.
 * @param[in,out]  rng                      Random number generator
 * @param[in]      p                        Probability of a connection
 * @param[in]      need_connected_meanings  If true, ensure that all meanings will be connected if possible.
 */
void BaseLexicalMatrix::initialize_random(RNG &rng, double p, bool need_connected_meanings) {
    std::binomial_distribution binom(this->n * this->m, p);
    unsigned e = binom(rng);
    this->initialize_random_edges(rng, e, need_connected_meanings);
}

/**
 * Clear and then initialize graph by randomly connecting edges until the given amount are connected.
 * All parameters are recalculated correctly afterwards.
 * @param[in,out]  rng                      Random number generator
 * @param[in]      edges                    Number of edges to have in the graph
 * @param[in]      need_connected_meanings  If true, ensure that all meanings will be connected if possible.
 */
void BaseLexicalMatrix::initialize_random_edges(RNG &rng, unsigned edges, bool need_connected_meanings) {
    this->clear();
    
    std::uniform_int_distribution<unsigned> word(0, this->n-1);
    std::uniform_int_distribution<unsigned> meaning(0, this->m-1);

    if (need_connected_meanings) {
        for (unsigned j=0; j<this->m; j++) {
            unsigned i = word(rng);
            this->add_edge_never_recalculate(i, j);
        }
        
        if (edges > this->m) {
            edges -= this->m;
        } else {
            edges = 0;
        }
    }
        
    while (edges > 0) {
        unsigned i = word(rng);
        unsigned j = meaning(rng);
        if (!this->get_a(i, j)) {
            this->add_edge_never_recalculate(i, j);
            edges--;
        }
    }
    
    this->recalculate(true);
}

/**
 * Clear and then initialize graph by connecting every possible edge.
 * All parameters are recalculated correctly afterwards.
 */
void BaseLexicalMatrix::initialize_complete_graph() {
    this->clear();
    
    for (unsigned i = 0; i < this->n; i++) {
        for (unsigned j = 0; j < this->m; j++) {
            this->add_edge_never_recalculate(i, j);
        }
    }
    
    this->recalculate(true);
}

/**
 * Clear and then initialize graph by connecting every word with one meaning and every meaning with one word.
 * All parameters are recalculated correctly afterwards.
 */
void BaseLexicalMatrix::initialize_bijective_graph() {
    this->clear();
    
    for (unsigned i = 0; i < std::min(this->n, this->m); i++) {
        this->add_edge_never_recalculate(i, i);
    }
    
    this->recalculate(true);
}

/**
 * Assign connected components a "color", an unsigned integer.
 * @param[out]  colorWords     colorWords[i] is the number assigned to the word i
 * @param[out]  colorMeanings  colorMeanings[j] is the number assigned to the meaning j
 * @return                     The number of connected components in the graph
 */
unsigned BaseLexicalMatrix::color_connected_nodes(std::vector<unsigned> &colorWords, std::vector<unsigned> &colorMeanings) const {
    unsigned color = 1;
    
    for (unsigned word=0; word<this->n; word++) {
        if (colorWords[word] == 0) {
            color_connected_nodes(word, true, color, colorWords, colorMeanings);
            color++;
        }
    }

    for (unsigned meaning=0; meaning<this->n; meaning++) {
        if (colorMeanings[meaning] == 0) {
            color_connected_nodes(meaning, false, color, colorWords, colorMeanings);
            color++;
        }
    }
    
    return color;
}

/**
 * Size of the largest connected component in the graph.
 * @return  The number of words and the number of meanings in the largest connected component.
 */
std::pair<unsigned, unsigned> BaseLexicalMatrix::largest_connected_component() const {
    // Assign an integer to each word/meaning such that all connectd ones share the same one.
    std::vector<unsigned> coloredWords(this->n);
    std::vector<unsigned> coloredMeanings(this->m);
    unsigned numColors = this->color_connected_nodes(coloredWords, coloredMeanings);

    // Count how many words and meanings there are of each color 
    std::vector<std::pair<unsigned, unsigned>> colorCounts(numColors);
    
    for (unsigned i=0; i<this->n; i++) {
        colorCounts[coloredWords[i]].first++;
    }
    for (unsigned j=0; j<this->m; j++) {
        colorCounts[coloredMeanings[j]].second++;
    }
    assert(colorCounts[0].first == 0 && colorCounts[0].second == 0);

    auto largestCount = std::max_element(
        colorCounts.begin(), colorCounts.end(),
        [](const std::pair<unsigned, unsigned> a, const std::pair<unsigned, unsigned> b) -> bool {
            return (a.first + a.second) < (b.first + b.second);
        });
    return *largestCount;
}

/**
 * Count the number of referentially useless words in the graph
 * A word is referentially useless if it's connected to the graph and for
 * each meaning it's connected to, p(s,r) <= p(s)*p(r)
 */ 
unsigned BaseLexicalMatrix::referentially_useless_words() const {
    unsigned count = 0;
    
    for (unsigned s=0; s<this->n; s++) {
        bool all = true;
        bool empty = true;
        
        for (unsigned r=0; r<this->m; r++) {
            if (this->get_a(s, r)) {
                empty = false;

                long double psr = this->get_psr(s, r);
                long double ps = this->get_ps(s);
                long double pr = this->get_pr(r);
                
                if (psr > ps*pr) {
                    all = false;
                    break;
                }
            }
        }
        
        if (all && !empty) {
            count++;
        }
    }
    
    return count;
}

/**
 * Return the graph's properties as a single object.
 * @param[out]  parameters  The object with the properties.
 */
void BaseLexicalMatrix::get_parameters(LexicalMatrixParameters<long double> &parameters) const {
    parameters.HS = this->get_HS();
    parameters.HR = this->get_HR();
    parameters.HSR = this->get_HSR();
    parameters.ISR = this->get_ISR();
    parameters.HSgR = this->get_HSgR();
    parameters.HRgS = this->get_HRgS();
}

#ifdef NDEBUG
#define assert_invariant(expr)                  \
    do {                                        \
    } while (0)
#else
#define assert_invariant(expr)                  \
    do {                                        \
        bool evaluated = expr;                  \
        if (!evaluated) {                       \
            std::cerr << msg << std::endl;      \
            this->print_properties();           \
            assert(expr);                       \
        }                                       \
    } while (0)
#endif

/// Ensure all invariants hold within epsilon
void BaseLexicalMatrix::sanity_checks() const {
    long double HS = this->get_HS();
    long double HR = this->get_HR();
    long double HSR = this->get_HSR();
    
    static const char *msg =
        "invariant assertion failure, matrix parameters dumped to stdout";
    
    assert_invariant(HS >= this->min_HS() - this->epsilon);
    assert_invariant(HS <= this->max_HS() + this->epsilon);

    assert_invariant(HR >= this->min_HR() - this->epsilon);
    assert_invariant(HR <= this->max_HS() + this->epsilon);

    assert_invariant(HSR >= this->min_HSR() - this->epsilon);
    assert_invariant(HSR <= this->max_HSR() + this->epsilon);
}

#undef assert_invariant

/**
 * Print the state of the graph to stdout.
 */
void BaseLexicalMatrix::print_properties() const {
    using namespace std;

    cout << "BASE PROPERTIES" << endl;
    cout << "---------------" << endl;

    cout << "epsilon = " << this->epsilon << endl;
    
    cout << "n = " << this->n << endl;
    cout << "m = " << this->m << endl;
    cout << "phi = " << this->phi << endl;
    
    cout << "a = {" << endl;
    for (unsigned i=0; i<this->n; i++) {
        cout << "\t";
        for (unsigned j=0; j<this->m; j++) {
            cout << " " << this->get_a(i, j);
            if (j >= 10) {
                cout << " ...";
                break;
            }
        }
        cout << endl;
        if (i >= 10) {
            cout << "\t ..." << endl;
            break;
        }
    }
    cout << "}" << endl;

    cout << "mu = ";
    pretty_print_vector(cout, this->mu);
    cout << endl;

    cout << "omega = ";
    pretty_print_vector(cout, this->omega);
    cout << endl;
    
    cout << "H(S,R) = " << this->get_HSR() << endl;
    cout << "H(S) = " << this->get_HS() << endl;
    cout << "H(R) = " << this->get_HR() << endl;
    cout << "I(S,R) = " << this->get_ISR() << endl;
    cout << "H(S|R) = " << this->get_HSgR() << endl;
    cout << "H(R|S) = " << this->get_HSgR() << endl;

    cout << "max H(S,R) = " << this->max_HSR() << endl;
    cout << "min H(S,R) = " << this->min_HSR() << endl;
    cout << "max H(S) = " << this->max_HS() << endl;
    cout << "min H(S) = " << this->min_HS() << endl;
    cout << "max H(R) = " << this->max_HR() << endl;
    cout << "min H(R) = " << this->min_HR() << endl;
}
