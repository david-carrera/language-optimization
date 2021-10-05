#include <LexicalMatrixStatistics.hpp>
#include <Util.hpp>

template<>
void LexicalMatrixStatistics<long double>::reset() {
    this->p_i = 0.0;
    this->mu_i = 0;
    this->a_i = 0;
    this->mu = 0;
}

template<>
void LexicalMatrixStatistics<Accumulator>::reset() {
    this->p_i.reset();
    this->mu_i.reset();
    this->a_i.reset();
    this->mu.reset();
}
template<>
LexicalMatrixStatistics<long double> LexicalMatrixStatistics<Accumulator>::operator/(long double other) {
    LexicalMatrixStatistics<long double> res;
    res.p_i = this->p_i.value() / other;
    res.mu_i = this->mu_i.value() / other;
    res.a_i = this->a_i.value() / other;
    res.mu = this->mu.value() / other;
    return res;
}

std::ostream& operator<<(std::ostream& out, const LexicalMatrixStatistics<long double> &statistics)
{
    std::streamsize ss = std::cout.precision();
    std::cout.precision(17);
    std::ostream& res = out <<
        statistics.p_i << "\t" <<
        statistics.mu_i << "\t" <<
        statistics.a_i << "\t" <<
        statistics.mu << "\n";
    std::cout.precision(ss);
    return res;
}
