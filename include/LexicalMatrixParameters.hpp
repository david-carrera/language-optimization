#ifndef LEXICAL_MATRIX_PARAMETERS_H
#define LEXICAL_MATRIX_PARAMETERS_H

#include <Util.hpp>
#include <iostream>

/**
 * Collection of information theoretical parameters from a lexical matrix.
 * The parameters H(S), H(R), H(S,R), I(S,R), H(S|R), H(R|S) in a single structure in order to more compactly deal with them as a single unit.
 * @tparam  T  Can be either a floating point value (e.g. long double) or an Accumulator
 */
template <class T>
struct LexicalMatrixParameters {
    /// H(S)
    T HS;
    /// H(R)
    T HR;
    /// H(S,R)
    T HSR;
    /// I(S,R)
    T ISR;
    /// H(S|R)
    T HSgR;
    /// H(R|S)
    T HRgS;

    LexicalMatrixParameters<Accumulator> &operator+=(const LexicalMatrixParameters<long double> &other);
    LexicalMatrixParameters<Accumulator> operator+(const LexicalMatrixParameters<long double> &other);
    LexicalMatrixParameters<long double> operator/(long double other);
    void reset();

    friend std::ostream& operator<<(std::ostream& out, const LexicalMatrixParameters<long double> &parameters);
};

#endif /* LEXICAL_MATRIX_PARAMETERS_H */
