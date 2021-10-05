#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <vector>
#include <random>

/// Default error tolerance for lexical matrix
extern const long double default_epsilon;

/// Random number generator to use throughout, the STL's Mersenne Twister engine.
typedef std::mt19937 RNG;

/// Utility for getting a user friendly string corresponding to a type name.
template <class T>
std::string type_name();

std::string vector_join(std::vector<std::string> elements, std::string sep);

/**
 * Print a vector through an output stream, without printing too many elements
 * and with nice curly braces.
 *
 * @param out Output stream
 * @param v Vector to print
 * @param limit How many elements to print. If vector has more elements, an ellipsis will be printed instead.
 */
template <class T>
void pretty_print_vector(std::ostream &out, const std::vector<T> v, unsigned limit=10) {
    out << "{";
    for (unsigned i=0; i<v.size(); i++) {
        out << " " << v[i];
        if (i >= limit) {
            out << " ...";
            break;
        }
    }
    out << " }";
}

/**
 * Utilty class to separate a time in seconds
 * Separate and store a time in seconds into seconds, minutes, hours, days and easily print it in a user friendly way.
 */
struct SeparatedTime {
    /// True if the time given in the constructor was infinite or nan
    bool infinite = false;

    /// Number of seconds, between 0 and 59
    int seconds;
    /// Number of minutes, between 0 and 59
    int minutes;
    /// Number of hours, between 0 and 23
    int hours;
    /// Number of days
    int days;

    SeparatedTime(double time);
    friend std::ostream& operator<<(std::ostream& out, const SeparatedTime& time);
};

/**
 * Utility function to cheaply reduce arithmetic error when doing floating point addition of many values
 */
class Accumulator {
private:
    /// Total value except for current
    long double total;
    /// Current value being added to
    long double current;
public:
    Accumulator();

    void add(long double value);

    Accumulator &operator +=(long double other);
    Accumulator &operator -=(long double other);
    Accumulator &operator =(long double other);

    long double value() const;
    void reset();
};

/**
 * Helper class to combine multiple output streams.
 * A class to make easier to send the same to multuple output streams at the same time, for instance to send a text to stderr but also log it to a file.
 */
class MultipleOstreams {
private:
    /// List of streams to write to
    std::vector<std::ostream*> streams;
public:
    /**
     * Add a stream to the collection
     */
    void add(std::ostream *stream);

    /// Flush all streams
    void flushall();

    /**
     * Output something through all streams.
     * @param   thing  The thing to output.
     * @return         A reference to itself
     */
    template <class T>
    MultipleOstreams &operator<<(const T &thing) {
        for (auto it : this->streams) {
            *it << thing;
        }
        return *this;
    }
};

#endif /* UTIL_H */
