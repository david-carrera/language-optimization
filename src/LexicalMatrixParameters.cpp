#include <LexicalMatrixParameters.hpp>

/**
 * Add all parameters at once into this object
 * @param[in]  other  The other object whose values will be added to this one
 * @return            Reference to this object
 */
template<>
LexicalMatrixParameters<Accumulator> &LexicalMatrixParameters<Accumulator>::operator+=(const LexicalMatrixParameters<long double> &other) {
    this->HS += other.HS;
    this->HR += other.HR;
    this->HSR += other.HSR;
    this->ISR += other.ISR;
    this->HSgR += other.HSgR;
    this->HRgS += other.HRgS;
    return *this;
}

/**
 * Add all parameters at once into a new object
 * @param[in]  other  The other object whose values will be added to this one
 * @return            The new object
 */
template<>
LexicalMatrixParameters<Accumulator> LexicalMatrixParameters<long double>::operator+(const LexicalMatrixParameters<long double> &other) {
    LexicalMatrixParameters<Accumulator> res;
    res.HS += this->HS;
    res.HR += this->HR;
    res.HSR += this->HSR;
    res.ISR += this->ISR;
    res.HSgR += this->HSgR;
    res.HRgS += this->HRgS;
    res.HS += other.HS;
    res.HR += other.HR;
    res.HSR += other.HSR;
    res.ISR += other.ISR;
    res.HSgR += other.HSgR;
    res.HRgS += other.HRgS;
    return res;
}

/**
 * Divide all parameters by the same value into a new object
 * @param[in]  other  The value
 * @return            The new object
 */
template<>
LexicalMatrixParameters<long double> LexicalMatrixParameters<Accumulator>::operator/(long double other) {
    LexicalMatrixParameters<long double> res;
    res.HS = this->HS.value() / other;
    res.HR = this->HR.value() / other;
    res.HSR = this->HSR.value() / other;
    res.ISR = this->ISR.value() / other;
    res.HSgR = this->HSgR.value() / other;
    res.HRgS = this->HRgS.value() / other;
    return res;
}

/**
 * Output parameters to a file in a CSV format.
 * @param[in,out]  out        Stream to write to.
 * @param[in]      parametrs  The parameters to write
 * @return                    The stream written to
 */
std::ostream& operator<<(std::ostream& out, const LexicalMatrixParameters<long double> &parameters)
{
    std::ostream& res = out << parameters.HS << "\t" << parameters.HR << "\t"
                            << parameters.HSR << "\t" << parameters.ISR << "\t"
                            << parameters.HSgR << "\t" << parameters.HRgS;
    return res;
}

/// Set all parameters to 0
template<>
void LexicalMatrixParameters<long double>::reset() {
    this->HS = 0;
    this->HR = 0;
    this->HSR = 0;
    this->ISR = 0;
    this->HSgR = 0;
    this->HRgS = 0;
}

/// Reset all accumulators
template<>
void LexicalMatrixParameters<Accumulator>::reset() {
    this->HS.reset();
    this->HR.reset();
    this->HSR.reset();
    this->ISR.reset();
    this->HSgR.reset();
    this->HRgS.reset();
}
