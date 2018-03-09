#include "Permutation.h"
#include "Exception.h"

unsigned int getThetaEffect(const unsigned int parity) {
    //Differential propagation
    //pattern 0x13
    return ((parity >> 1) | ((parity &  0x1) << (LANESIZE-1))) ^ \
           ((parity >> 2) | ((parity &  0x3) << (LANESIZE-2))) ^ \
           ((parity >> 5) | ((parity & 0x1f) << (LANESIZE-5)));
    //Linear propagation, transpose of Z
    //const unsigned int mask = (1 << LANESIZE) - 1;
    /*
    return (((parity << 1) & mask) | (parity >> (LANESIZE-1))) ^ \
           (((parity << 2) & mask) | (parity >> (LANESIZE-2))) ^ \
           (((parity << 5) & mask) | (parity >> (LANESIZE-5)));
    */
}

void applyTheta(State& state) {
    applyTheta(state, getThetaEffect(state.getParity()));
}

void applyTheta(State& state, const unsigned int effect) {
    for(unsigned int i = 0; i < COLUMNSIZE; ++i)
        state[i] ^= effect;
}

void applyDispersion(State& state) {
    const unsigned int tmp = state[0];
    state[0] = state[1];
    state[1] = state[2];
    state[2] = state[3];
    state[3] = tmp;

    state.rotateRow(1, 10);
    state.rotateRow(2, 3);
    state.rotateRow(3, 14);
}

void applyInverseDispersion(State& state) {
    state.rotateRow(1, LANESIZE-10);
    state.rotateRow(2, LANESIZE-3);
    state.rotateRow(3, LANESIZE-14);

    const unsigned int tmp = state[3];
    state[3] = state[2];
    state[2] = state[1];
    state[1] = state[0];
    state[0] = tmp;
}

void printDispersed(std::ostream& out, const State& state) {
    State copy(state);
    applyDispersion(copy);
    out << copy;
}

void thetaCompatibleStates(const State& state, std::vector<State>& compatibleStates, const unsigned int minWeight, const unsigned int maxWeight) {
    if(maxWeight == 0)
        return;
    if(minWeight > maxWeight)
        return;
    if(maxWeight > 100)
        throw Exception("That's pretty large. Are you sure there was no integer underflow somewhere?");

    //first, consider all columns with >= 2 active cells
    unsigned int sum[LANESIZE];
    state.getSum(sum);

    unsigned int nrRelevantColumns = 0;
    unsigned int relevantColumns[LANESIZE];
    for(unsigned int i = 0; i < LANESIZE; ++i)
        if(sum[i] >= 2)
            relevantColumns[nrRelevantColumns++] = (2*LANESIZE - i - 2) % LANESIZE;

    //all of those columns can have parity 0 or 1, iterate over all possibilities and modify the effect accordingly
    unsigned int effect = getThetaEffect(state.getParity());
    for(unsigned int i = 0; i < (1u << nrRelevantColumns); ++i) {
        State copy(state);
        unsigned int modifiedEffect = effect;
        for(unsigned int j = 0; j < nrRelevantColumns; ++j) {
            if(i & (1<<j)) {
                modifiedEffect ^= 1 << relevantColumns[j];
            }
        }

        //apart from this branching, it is also possible that applying the effect may or may not cancel out an active cell
        //find all possible states with this type of branching
        std::vector<State> effectBranching;
        effectBranching.reserve(256);
        recurThetaCompatibleStates(copy, effectBranching, modifiedEffect, sum);


        //finish the CPM and add them to our list
        for(std::vector<State>::iterator it = effectBranching.begin(); it != effectBranching.end(); ++it) {
            applyTheta(*it, modifiedEffect);
            if(it->hammingWeight() >= minWeight && it->hammingWeight() <= maxWeight)
                compatibleStates.push_back(*it);
        }
    }
}

void recurThetaCompatibleStates(State& state, std::vector<State>& effectBranching, const unsigned int effect, const unsigned int sum[LANESIZE], const unsigned int colIndex) {
    //we're finished when we've reached the end of the state
    if(colIndex >= LANESIZE)
        effectBranching.push_back(std::move(state));

    //if not their yet, look at affected columns with at least two active cells
    else {
        const unsigned int sumCol = sum[colIndex];
        const unsigned int intNextColIndex = (2*LANESIZE-colIndex-2)%LANESIZE;
        const unsigned int intColIndex = (intNextColIndex+1)%LANESIZE;

        if(((effect >> intColIndex) & 1) && sumCol >= 2) {
            bool diffPossible = false; //e.g. AB, ABC, ABCD
            bool pairPossible = false; //e.g. AAB, AABB
            bool tripPossible = false; //e.g. AAAB
            bool samePossible = false; //e.g. AA, AAA, AAAA
            //if this column has odd parity, not that the effect has already been shifted
            if((effect >> intNextColIndex) & 1) {
                //which cells can stay/disappear depends on the number of active cells in that column
                if(sumCol == 2) {
                    //then active cells must have a different value
                    diffPossible = true;
                }
                else if(sumCol == 3) {
                    //then every combination is possible
                    diffPossible = true;
                    pairPossible = true;
                    samePossible = true;
                }
                else if(sumCol == 4) {
                    //every combination is possible, except all equal
                    diffPossible = true;
                    pairPossible = true;
                    tripPossible = true;
                }
                else {
                    //note that this code is not generic for arbitrary COLUMNSIZE
                    throw Exception("This should not be possible, are you sure COLUMNSIZE == 4?");
                }
            }
            else { //if this column has even parity
                if(sumCol == 2) {
                    //then active cells must have the same value, so either both stay or both disappear
                    samePossible = true;
                }
                else if(sumCol == 3) {
                    //then all active cells must be different, at most 1 can disappear
                    diffPossible = true;
                }
                else if(sumCol == 4) {
                    //then either all different, or two pairs, or all the same
                    diffPossible = true;
                    pairPossible = true;
                    samePossible = true;
                }
                else {
                    //note that this code is not generic for arbitrary COLUMNSIZE
                    throw Exception("This should not be possible, are you sure COLUMNSIZE == 4?");
                }

            }

            //actually get all next possible states

            //there is always the possibility that nothing cancels
            State copy(state);
            copy.resetColumn(0,colIndex);
            recurThetaCompatibleStates(copy, effectBranching, effect, sum, colIndex+1);
            //if all are equal, all can cancel
            if(samePossible) {
                copy = State(state);
                recurThetaCompatibleStates(copy, effectBranching, effect, sum, colIndex+1);
            }
            //if all are different, at most one cancels
            if(diffPossible) {
                for(unsigned int i = 0; i < COLUMNSIZE; ++i) {
                    if((state[i] >> intColIndex) & 1) {
                        copy = State(state);
                        copy.resetColumn(1 << i, colIndex);
                        recurThetaCompatibleStates(copy, effectBranching, effect, sum, colIndex+1);
                    }
                }
            }
            //if there is a pair, two will cancel
            if(pairPossible) {
                for(unsigned int i = 0; i < COLUMNSIZE - 1; ++i) {
                    if((state[i] >> intColIndex) & 1) {
                        for(unsigned int j = i + 1; j < COLUMNSIZE; ++j) {
                            if((state[j] >> intColIndex) & 1) {
                                copy = State(state);
                                copy.resetColumn((1 << i) | (1 << j), colIndex);
                                recurThetaCompatibleStates(copy, effectBranching, effect, sum, colIndex+1);
                            }
                        }
                    }
                }
            }
            //if there is a triple, three will cancel
            if(tripPossible) {
                for(unsigned int i = 0; i < COLUMNSIZE; ++i) {
                    if((state[i] >> intColIndex) & 1) { //will aways be true here, as sumcol == 4, just to make it more generic
                        copy = State(state);
                        copy.unsetColumn(1 << i, colIndex);
                        recurThetaCompatibleStates(copy, effectBranching, effect, sum, colIndex+1);
                    }
                }
            }
        }
        else
            recurThetaCompatibleStates(state, effectBranching, effect, sum, colIndex+1);
    }
}
