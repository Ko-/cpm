#include "TrailExtension.h"
#include "Permutation.h"
#include <algorithm>

TrailExtension::TrailExtension(TwoRoundTrailCore& core)
{
    Trail t;
    t.append(core.stateA, core.w0);
    applyDispersion(core.stateB);
    t.append(core.stateB, core.w1);
    trails.push_back(t);
}

TrailExtension::TrailExtension(const Trail& core)
{
    trails.push_back(core);
}

void TrailExtension::extendForward(const unsigned int maxWeight)
{
    std::vector<State> compatible;
    std::vector<Trail> newTrails;
    for(Trail& t : trails) {
        unsigned int weightSpent = 0;
        for(unsigned int i = t.initialIndex; i < t.getLength(); ++i)
            weightSpent += t.weights[i];
        const unsigned int minWeight = std::max(1, (int)t.weights[t.initialIndex] + (int)t.weights[t.initialIndex + 1] - (int)t.weights.back());

        thetaCompatibleStates(t.states.back(), compatible, minWeight, maxWeight - weightSpent);

        if(compatible.empty())
            newTrails.push_back(t);

        for(State& s : compatible) {
            Trail tnew(t);
            applyDispersion(s);
            tnew.append(s, s.hammingWeight());
            newTrails.push_back(tnew);
        }
        compatible.clear();
    }
    trails = newTrails;
}

void TrailExtension::extendBackward(const unsigned int maxWeight)
{
    std::vector<State> compatible;
    std::vector<Trail> newTrails;
    for(Trail& t : trails) {
        unsigned int weightSpent = 0;
        for(unsigned int i = 0; i < 2 + t.initialIndex; ++i)
            weightSpent += t.weights[i];
        const unsigned int minWeight = std::max(1, (int)t.weights[t.initialIndex] + (int)t.weights[t.initialIndex + 1] - (int)t.weights.front());

        State copy(t.states.front());
        applyInverseDispersion(copy);

        thetaCompatibleStates(copy, compatible, minWeight, maxWeight - weightSpent);

        if(compatible.empty())
            newTrails.push_back(t);

        for(State& s : compatible) {
            Trail tnew(t);
            tnew.prepend(s, s.hammingWeight());
            newTrails.push_back(tnew);
        }
        compatible.clear();
    }
    trails = newTrails;
}

void TrailExtension::getBestTrail(Trail& bestTrail, const unsigned int rounds)
{
    unsigned int bestWeight = -1; //underflow to maximum value
    bestTrail.clear();
    for(Trail& t : trails) {
        if(t.prune(rounds, bestWeight) && t.totalWeight < bestWeight) {
            bestWeight = t.totalWeight;
            bestTrail.set(t);
        }
    }
}

static bool isOutsideKernel(const Trail & t) {
    if(t.getLength() != 3) {
        return true; //ugly hack, should pass a parameter
    }
    const unsigned int constWeight = t.weights[0];
    for(const unsigned int& w : t.weights) {
        if(w != constWeight)
            return true;
    }
    for(const State& s : t.states) {
        unsigned int sum[LANESIZE];
        s.getSum(sum);
        for (unsigned int z = 0; z < LANESIZE; z++) {
            if(sum[z] == 1)
                return true;
        }
    }
    return false;
}

void TrailExtension::removeOutsideKernel() {
    trails.erase(std::remove_if(trails.begin(), trails.end(), isOutsideKernel), trails.end());
}

void TrailExtension::getStats(std::vector<unsigned int>& stats, const unsigned int rounds) const {
    for(const Trail& t : trails) {
        if(t.states.size() == rounds) {
            if(stats.size() <= t.totalWeight)
                stats.resize(t.totalWeight + 1, 0);
            stats[t.totalWeight]++;
        }
    }
}
