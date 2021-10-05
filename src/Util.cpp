#include <Util.hpp>
#include <cmath>

extern const long double default_epsilon = 1e-10;

template <>
std::string type_name<int>() {
    return "integer";
}
template <>
std::string type_name<unsigned int>() {
    return "unsigned integer";
}
template <>
std::string type_name<double>() {
    return "real";
}
template <>
std::string type_name<bool>() {
    return "yes/no";
}

/**
 * Join elements of vector with given separator.
 *
 * @param elements Vector of strings to join.
 * @param sep Separator for the join.
 * @returns A string with all the elements joined by sep.
 */
std::string vector_join(std::vector<std::string> elements, std::string sep) {
    std::string total = "";
    if (elements.size() == 0) {
        return total;
    }
    
    unsigned long last = elements.size()-1;
    for (unsigned long i = 0; i<last; i++) {
        total += elements[i] + sep;
    }
    if (elements.size() > 0) {
        total += elements[last];
    }
    return total;
}

/**
 * Constructor. Separates the given time.
 *
 * @param time Time in seconds.
 */
SeparatedTime::SeparatedTime(double time) {
    if (std::isinf(time)) {
        this->infinite = true;
        return;
    }
    
    int t = static_cast<int>(std::round(time));

    this->seconds = t % 60;
    t /= 60;

    this->minutes = t % 60;
    t /= 60;

    this->hours = t % 24;
    t /= 24;

    this->days = t;
}

/**
 * Output time in a user friendly way.
 * @param  out   Stream to write the time to
 * @param  time  Time to write
 * @return       Stream time was written into
 */
std::ostream& operator<<(std::ostream& out, const SeparatedTime& time) {
    if (time.infinite) {
        return out << "inf";
    }
    
    return out << time.days << "d "
               << time.hours << "h "
               << time.minutes << "m "
               << time.seconds << "s";
}

/// Constructor
Accumulator::Accumulator() : total(0.0), current(0.0) { }

/**
 * Add a value.
 * @param  value  Value to add
 */
void Accumulator::add(long double value) {
    this->current += value;
    if (std::abs(this->current) > this->total) {
        this->total += this->current;
        this->current = 0.0;
    }
}

/**
 * Add a value.
 * @param  other  Value to add
 * @return        Reference to this accumulator
 */
Accumulator &Accumulator::operator+=(long double other) {
    this->add(other);
    return *this;
}

/**
 * Subtract a value.
 * @param  other  Value to subtract
 * @return        Reference to this accumulator
 */
Accumulator &Accumulator::operator-=(long double other) {
    this->add(-other);
    return *this;
}

/**
 * Assign a value.
 * @param  other  Value to assign.
 * @return        Reference to this accumulator
 */
Accumulator &Accumulator::operator =(long double other) {
    this->total = other;
    this->current = 0.0L;
    return *this;
}

/**
 * Return value of the accumulator
 * @return  #value + #current
 */
long double Accumulator::value() const {
    return this->total + this->current;
}

/// Set to zero
void Accumulator::reset() {
    this->total = 0.0L;
    this->current = 0.0L;
}


/**
 * Add an output stream.
 *
 * @param stream The stream to add.
 */
void MultipleOstreams::add(std::ostream *stream) {
    this->streams.push_back(stream);
}

/**
 * Flush all streams.
 */
void MultipleOstreams::flushall() {
    for (auto it : this->streams) {
        it->flush();
    }
}
