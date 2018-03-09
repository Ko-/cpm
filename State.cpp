#include "State.h"

/*
    _ _ _ _ _ _ _ 
  3|_|_|_|_|_|_|_|
^  |_|_|_|_|_|_|_|
|  |_|_|_|_|_|_|_|
y 0|_|_|_|_|_|_|_|
    0           7 
    z -->
*/

State::State() {
    clear();
}

State::State(const State& s) {
    for(unsigned int i=0; i<COLUMNSIZE; ++i)
        rows[i] = s.rows[i];
}

void State::setBit(const unsigned int y, const unsigned int z) {
    rows[y] |= (1 << (LANESIZE-1-z));
}

void State::setColumn(const unsigned int value, const unsigned int z) {
    for(unsigned int i=0; i<COLUMNSIZE; ++i)
        if((value >> i) & 1)
            rows[i] ^= (1 << (LANESIZE-1-z));
}

void State::resetColumn(const unsigned int value, const unsigned int z) {
    for(unsigned int i=0; i<COLUMNSIZE; ++i)
        if(((rows[i] >> (LANESIZE-1-z)) & 1) != ((value >> i) & 1))
            rows[i] ^= (1 << (LANESIZE-1-z));

}

void State::unsetColumn(const unsigned int value, const unsigned int z) {
    for(unsigned int i=0; i<COLUMNSIZE; ++i) {
        unsigned int mask = (1 << LANESIZE) - 1;
        if((value >> i) & 1)
            mask ^= (1 << (LANESIZE-1-z));
        rows[i] &= mask;
    }
}

void State::clear() {
    for(unsigned int i=0; i<COLUMNSIZE; ++i)
        rows[i] = 0;
}

unsigned int State::getParity() const {
    unsigned int parity = 0;
    for(unsigned int i = 0; i < COLUMNSIZE; ++i)
        parity ^= rows[i];
    return parity;
}

void State::getSum(unsigned int sum[]) const {
    for(unsigned int i = 0; i < LANESIZE; ++i) {
        sum[i] = 0;
        for(unsigned int j = 0; j < COLUMNSIZE; ++j)
            sum[i] += (rows[j] >> (LANESIZE - i - 1)) & 1;
    }
}

#if defined(_MSC_VER)
    #include <intrin.h>
#endif
unsigned int State::hammingWeight() const {
    unsigned int sum = 0;
    for(unsigned int i = 0; i < COLUMNSIZE; ++i) {
        //since SSE4a/Nehalem, there is an efficient popcnt instruction that we want to use.
#if defined(__GNUC__)
        sum += __builtin_popcount(rows[i]);
#elif defined(_MSC_VER)
        sum += __popcnt(rows[i]);
#else
        unsigned int row = sum[i];
        while(row) {
            ++sum;
            row &= row-1;
        }
#endif
    }
    return sum;
}

void State::rotateRow(const unsigned int y, unsigned int dist) {
    rows[y] = (rows[y] >> dist) | ((rows[y] << (LANESIZE - dist)) & ((1u << LANESIZE) - 1));
}

void State::printRow(std::ostream& out, const int y) const {
    for(int i=LANESIZE-1; i>=0; --i)
        out << ((rows[y] >> i) & 1);
}

std::ostream& operator<<(std::ostream& out, const State& state){
    for(int i=COLUMNSIZE-1; i>=0; --i) {
        for(int j=LANESIZE-1; j>=0; --j)
            out << ((state.rows[i] >> j) & 1);
        if(i>0)
            out << "\n";
    }
    return out;
}

const unsigned int& State::operator[](const int index) const {
    return rows[index];
}

unsigned int& State::operator[](const int index) {
    return rows[index];
}

State& State::operator=(const State& s) {
    if(this != &s) {
        for(unsigned int i=0; i<COLUMNSIZE; ++i)
            rows[i] = s.rows[i];
    }
    return *this;
}
