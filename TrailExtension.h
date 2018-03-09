#ifndef _TRAILEXPANSION_
#define _TRAILEXPANSION_

#include <vector>
#include "MyTree.h"
#include "Trail.h"

class TrailExtension {
public:
    std::vector<Trail> trails;
     
public:
    explicit TrailExtension(TwoRoundTrailCore& core);

    explicit TrailExtension(const Trail& core);

    void extendForward(const unsigned int maxWeight);

    void extendBackward(const unsigned int maxWeight);

    void getBestTrail(Trail& trail, const unsigned int rounds);

    void removeOutsideKernel();

    void getStats(std::vector<unsigned int>& stats, const unsigned int rounds) const;
};

#endif
