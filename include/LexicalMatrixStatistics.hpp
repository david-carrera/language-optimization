#ifndef LEXICAL_MATRIX_STATISTICS_H
#define LEXICAL_MATRIX_STATISTICS_H

#include <iostream>

/**
 * Collection of statistical information from a lexical matrix.
 * The probability, mu and age by rank as well as the mu by mu.
 */
template <class T>
struct LexicalMatrixStatistics {
    /// Probability by rank
    T p_i;
    /// mu by rank
    T mu_i;
    /// Age by rank
    T a_i;
    /// mu by mu
    T mu;


    /**
     * Divide all elements by the same number
     * @param  other  The number to divide by
     * @return        A new object with the result
     */
    LexicalMatrixStatistics<long double> operator/(long double other);

    /// Set all values to 0
    void reset();
    
    /**
     * Output parameters to a file in a CSV format.
     * @param[in,out]  out        Stream to write to.
     * @param[in]      parametrs  The parameters to write
     * @return                    The stream written to
     */
    friend std::ostream& operator<<(std::ostream& out, const LexicalMatrixStatistics<long double> &statistics);
};

#endif /* LEXICAL_MATRIX_STATISTICS_H */
