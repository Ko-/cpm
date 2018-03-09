#include <algorithm>
#include <atomic>
#include <future>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <tuple>
#include <vector>


//given a maxweight
//for all possible odd/affected columns as given by distribution.cpp (until maxweight)
 //fork to new thread?
  //for all possible values for odd columns (1,2,4,8,e,c,b,7)
   //generate state-before-theta, state-after-theta
   //save this
   //compute yMin
   //for all possible orbitals until maxweight reached
    //update state-before-theta, state-after-theta
    //save this

// Types

typedef unsigned int uint;

struct config {
    uint odd;
    uint affected;
};

struct state {
    uint a[16];
    uint b[16];
    uint weight;

    void setweight() {
        weight = 0;
        for (uint i = 0; i < 16; ++i) {
            weight += __builtin_popcount(a[i]);
            weight += __builtin_popcount(b[i]);
        }
    }
};

// Globals

const state emptystate = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},0};
const uint oddvalues[8] = {1,2,4,8,7,0xb,0xd,0xe};
const uint orbitals[6] = {3,5,9,6,0xa,0xc};
const uint orbitaltonewymin[6] = {2,4,4,4,4,4};
const uint ymintoorbitalindex[3] = {0,3,5};
std::atomic_uint totalstates(0);
std::mutex mymutex;

// Functions

std::ostream& operator<<(std::ostream& out, const state& s) {
    for (size_t i = 0; i < 16; i += 2) {
        out << (unsigned char)(s.a[i] | (s.a[i+1] << 4));
    }
    for (size_t i = 0; i < 16; i += 2) {
        out << (unsigned char)(s.b[i] | (s.b[i+1] << 4));
    }
    return out;
}

template<typename t>
inline void getindices(t value, std::vector<t>& indices) {
    for (t i = 0; value > 0; ++i, value >>= 1) {
        if (value & 1) {
            indices.push_back(i);
        }
    }
}

void addorbitals(const state& start, const std::vector<size_t>& affectedindices, const uint maxweight, std::ostream& out) {
    //compute ymin
    std::vector<uint> yminstart(16, 0);
    for (uint i = 0; i < 16; ++i) {
        //if the column has already something in the highest two positions or if it is affected, we don't want orbitals
        if ((start.a[i] & 0xc) || (std::find(affectedindices.begin(), affectedindices.end(), i) != affectedindices.end())) yminstart[i] = 4;
        else if (start.a[i] & 2) yminstart[i] = 2;
        else if (start.a[i] & 1) yminstart[i] = 1;
    }

    std::queue<std::tuple<state, std::vector<uint>>> todo;
    todo.push(std::make_tuple(start,yminstart));
    while (!todo.empty()) {
        state s = std::get<0>(todo.front());
        std::vector<uint> ymin = std::get<1>(todo.front());
        //bingo!
        totalstates++;
        mymutex.lock();
        out << s;
        mymutex.unlock();

        todo.pop();

        if (s.weight + 4 <= maxweight) {
            //add another orbital where possible
            for (uint i = 0; i < 16; ++i) {
                if (ymin[i] <= 2) {
                    //possible in this column
                    for (uint j = ymintoorbitalindex[ymin[i]]; j < 6; ++j) {
                        state copy = s;
                        std::vector<uint> ymincopy = ymin;
                        copy.a[i] ^= orbitals[j];
                        copy.b[i] ^= orbitals[j];
                        copy.weight += 4;
                        ymincopy[i] = orbitaltonewymin[j];
                        todo.push(std::make_tuple(copy,ymincopy));
                    }
                }
            }
        }
    }


}

void _bruteforce(const config c, const uint maxweight, std::ostream& out) {
    std::vector<size_t> oddindices;
    std::vector<size_t> affectedindices;

    getindices<size_t>(c.odd, oddindices);
    getindices<size_t>(c.affected, affectedindices);

    const uint nrpossibilities = 1 << (3 * oddindices.size()); // == pow(8,HW(c.odd))

    for (uint i = 0; i < nrpossibilities; ++i) {
        state s = emptystate;
        uint divider = 1;

        for (size_t& j : oddindices) {
            s.a[j] = oddvalues[(i/divider)%8];
            s.b[j] = s.a[j];
            divider *= 8;
        }

        for (size_t& j : affectedindices) {
            s.b[j] ^= 0xf;
        }

        s.setweight();

        if (s.weight <= maxweight) {
            addorbitals(s, affectedindices, maxweight, out);
        }
    }
}

void bruteforce(uint maxweight, std::ostream& out) {
    //pattern 0x13
    const config configs[55] = {{1,38},{3,106},{7,242},{9,278},{11,346},{13,398},{15,450},{17,582},{19,522},{23,658},{35,1194},{39,1074},{47,1282},{77,2062},{79,2114},{147,4362},{151,4498},{155,4154},{159,4258},{215,6162},{275,9226},{303,8962},{305,8326},{311,8274},{431,12290},{591,17474},{619,16410},{623,16514},{1103,36930},{1235,32906},{1239,32786},{2479,3},{4369,17476},{4371,17416},{4399,17152},{4405,16412},{4407,16464},{4527,20480},{4685,9228},{4687,9280},{4715,8216},{4719,8320},{4883,2056},{4887,2192},{4941,524},{4943,576},{4947,392},{4951,272},{4959,32},{5069,4364},{5071,4416},{5079,4624},{6831,2561},{9903,518},{13527,8212}};

    std::vector<std::future<void>> futures;

    for (const config& c : configs) {
        futures.push_back(std::async(std::launch::async, _bruteforce, c, maxweight, std::ref(out)));
    }

    for (auto& f : futures) {
        f.wait();
    }
}

int main(int argc, char** argv) {
    uint maxweight = argc > 1 ? atoi(argv[1]) : 10;
    std::ofstream log("0x13_" + std::to_string(maxweight) + ".txt", std::ios::binary);
    bruteforce(maxweight, log);
    std::cout << totalstates.load() << std::endl;
    log.close();
    return 0;
}
