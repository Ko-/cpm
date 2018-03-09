#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include "Permutation.h"
#include "Trail.h"
#include "TrailExtension.h"

std::mutex mymutex;
std::vector<unsigned int> trailCount;

/* This function looks at the state before the CPM and determines the minimum y positions where orbitals can be placed. Additionally, it tests whether the state is in the column parity kernel. */
void loadYMin(const TwoRoundTrailCore& trailcore, std::vector<unsigned int>& yMin, bool& kernel)
{
    unsigned int parity = trailcore.stateA.getParity();
    unsigned int effect = getThetaEffect(parity);
    kernel = true;
    for (unsigned int z = 0; z < LANESIZE; z++) {
        bool odd = (parity >> (LANESIZE-1-z)) & 1;
        if(odd)
            kernel = false;
        bool affected = (effect >> (LANESIZE-1-z)) & 1;
        if (affected) {
            yMin[z] = COLUMNSIZE; // no orbitals here
        }
        else if(odd) { // UOC => there is a single active bit
            for (unsigned int y = 0; y < COLUMNSIZE; y++){
                if (((trailcore.stateA[y] >> (LANESIZE-1-z)) & 1) != 0){
                    yMin[z] = y + 1;
                    break;
                }
            }
        }
    }
}


/* This function outputs 2-round trail cores outside the kernel with cost below given limit. */
void countTrailCores(const unsigned int maxWeight)
{
    time_t timerStart, timerEnd;
    std::vector<unsigned int> counts(maxWeight+1, 0);
    std::vector<unsigned int> kernelcounts(maxWeight+1, 0);

    time(&timerStart);

    TwoRoundTrailCoreCostFunction<Column> costFRun;
    ColumnsSet colSet;
    TwoRoundTrailCoreStack cacheRun;

    RunTreeIterator iteratorRun(colSet, cacheRun, costFRun, maxWeight);

    for (; !iteratorRun.isEnd(); ++iteratorRun) {
        TwoRoundTrailCore nodeRun = *iteratorRun;
        unsigned int costNodeRun = nodeRun.w0 + nodeRun.w1;
        bool completeNodeRun = nodeRun.complete;
        if (costNodeRun <= maxWeight && completeNodeRun){
            TwoRoundTrailCoreStack cacheOrb(nodeRun.stateA, nodeRun.stateB, nodeRun.w0, nodeRun.w1, completeNodeRun, nodeRun.zPeriod);
            TwoRoundTrailCoreCostFunction<Orbital> costFOrb;

            // create yMin
            std::vector<unsigned int> yMin(LANESIZE, 0);
            bool kernel;
            loadYMin(nodeRun, yMin, kernel);

            // orbital tree with parity-bare trail core at root
            OrbitalsSet orbSet(kernel, yMin);
            OrbitalTreeIterator iteratorOrb(orbSet, cacheOrb, costFOrb, maxWeight);
            for (; !iteratorOrb.isEnd(); ++iteratorOrb) {
                TwoRoundTrailCore nodeOrb = *iteratorOrb;
                counts[nodeOrb.w0 + nodeOrb.w1]++;
                if(kernel)
                    kernelcounts[nodeOrb.w0 + nodeOrb.w1]++;

            }
        }
    }
    for(unsigned int i=0;i<counts.size();++i)
        std::cout << i << ": " << counts[i] << std::endl;
    std::cout << "------------" << std::endl;
    for(unsigned int i=0;i<kernelcounts.size();++i)
        std::cout << i << ": " << kernelcounts[i] << std::endl;
    time(&timerEnd);
    std::cerr << "Execution time: " << timerEnd-timerStart << std::endl;
}

/* Helper function for findTrails that will be multithreaded. */
void _findTrails(const unsigned int rounds, TwoRoundTrailCore nodeRun, std::vector<Trail>& partialResults, const unsigned int maxWeight)
{
    unsigned int minWeight = -1; // underflows
    Trail bestTrail, trail;

    TwoRoundTrailCoreStack cacheOrb(nodeRun.stateA, nodeRun.stateB, nodeRun.w0, nodeRun.w1, true, nodeRun.zPeriod);
    TwoRoundTrailCoreCostFunction<Orbital> costFOrb;

    // create yMin
    std::vector<unsigned int> yMin(LANESIZE, 0);
    bool kernel;
    loadYMin(nodeRun, yMin, kernel);

    OrbitalsSet orbSet(kernel, yMin);
    OrbitalTreeIterator iteratorOrb(orbSet, cacheOrb, costFOrb, (maxWeight/rounds)*2+1);

    for (; !iteratorOrb.isEnd(); ++iteratorOrb) {

        TrailExtension ext(*iteratorOrb);
        for(unsigned int i=0;i<rounds-2;++i)
            ext.extendForward(maxWeight);
        for(unsigned int i=0;i<rounds-2;++i)
            ext.extendBackward(maxWeight);

        ext.getBestTrail(trail,rounds);
        /*
        TrailExtension ext(*iteratorOrb);
        ext.extendForward(maxWeight);
        ext.getStats(trailCount, rounds);

        TrailExtension ext2(*iteratorOrb);
        ext2.extendBackward(maxWeight);
        ext2.getStats(trailCount, rounds);

        ext.getBestTrail(trail, rounds);
        */
        if(trail.getLength() == rounds && trail.totalWeight > 0 && trail.totalWeight < minWeight && trail.totalWeight <= maxWeight) {
            minWeight = trail.totalWeight;
            bestTrail.set(trail);
        }
        /*
        ext2.getBestTrail(trail, rounds);

        if(trail.getLength() == rounds && trail.totalWeight > 0 && trail.totalWeight < minWeight && trail.totalWeight <= maxWeight) {
            minWeight = trail.totalWeight;
            bestTrail.set(trail);
        }
        */
    }

    bool addTrail = bestTrail.getLength() == rounds;
    if(!addTrail)
        return;

    mymutex.lock();
    for(std::vector<Trail>::const_iterator it = partialResults.begin(); addTrail && it != partialResults.end(); ++it) {
        if(it->totalWeight < minWeight)
            addTrail = false;
    }
    if(addTrail)
        partialResults.push_back(bestTrail);
    mymutex.unlock();
    if(addTrail)
        std::cout << "Update best trail:\n" << bestTrail << std::endl;
}

void findTrails(const unsigned int rounds, const unsigned int maxWeight) {
    time_t timerStart, timerEnd;

    // timer start
    time(&timerStart);

    TwoRoundTrailCoreCostFunction<Column> costFRun;
    ColumnsSet colSet;
    TwoRoundTrailCoreStack cacheRun;

    RunTreeIterator iteratorRun(colSet, cacheRun, costFRun, (maxWeight/rounds)*2+1);

    std::vector<Trail> partialResults;
    std::vector<std::future<void>> futures;

    //skip the first one, that's always the empty trailcore
    ++iteratorRun;

    double load[1];
    const unsigned int hwconcurrency = std::thread::hardware_concurrency();

    for (; !iteratorRun.isEnd(); ++iteratorRun) {
        TwoRoundTrailCore nodeRun = *iteratorRun;
        if (nodeRun.w0 + nodeRun.w1 <= (maxWeight/rounds)*2+1 && nodeRun.complete){

            //manage threads by watching current load, good parameters depend heavily on machine
            getloadavg(load,1);
            while(load[0] > hwconcurrency) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                getloadavg(load,1);
            }

            futures.push_back(std::async(_findTrails, rounds, nodeRun, std::ref(partialResults), maxWeight));
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
    }

    std::cout << futures.size() << " asynchronous tasks\n";

    for(std::vector<std::future<void>>::iterator it = futures.begin(); it != futures.end(); ++it) {
        it->wait();
    }

    unsigned int minWeight = -1; // underflows
    Trail bestTrail;
    for(std::vector<Trail>::iterator it = partialResults.begin(); it != partialResults.end(); ++it) {
        if(it->totalWeight < minWeight) {
            minWeight = it->totalWeight;
            bestTrail.set(*it);
        }
    }
    time(&timerEnd);
    std::cerr << "Execution time: " << timerEnd-timerStart << std::endl;
    std::cout << "Best trail:\n" << bestTrail << std::endl;
}

bool isInKernel(const TwoRoundTrailCore& trailcore) {
    unsigned int sumA[LANESIZE], sumB[LANESIZE];
    State s = trailcore.stateB;
    applyDispersion(s);
    trailcore.stateA.getSum(sumA);
    s.getSum(sumB);
    for (unsigned int z = 0; z < LANESIZE; z++) {
        if((sumA[z] == 1) || (sumB[z] == 1))
            return false;
    }
    /*
    for(unsigned int z = 0; z < LANESIZE; z++) {
        std::cout << sumA[z] << " ";
    }
    std::cout << "| ";
    for(unsigned int z = 0; z < LANESIZE; z++) {
        std::cout << sumB[z] << " ";
    }
    std::cout << std::endl;
    */
    return true;
}

void countInKernelTrails(const unsigned int rounds, const unsigned int maxWeight) {
    //assumes rounds == 3, doesn't really handle other values yet
    //fixed in TrailExtension::removeOutsideKernel
    time_t timerStart, timerEnd;
    unsigned int q = 0;

    // timer start
    time(&timerStart);

    TwoRoundTrailCoreCostFunction<Column> costFRun;
    RunTreeIterator iteratorRun(ColumnsSet(), TwoRoundTrailCoreStack(), costFRun, (maxWeight/rounds)*2+1);

    //skip the first one, that's always the empty trailcore
    ++iteratorRun;

    for (; !iteratorRun.isEnd(); ++iteratorRun) {
        TwoRoundTrailCore nodeRun = *iteratorRun;

        // create yMin
        std::vector<unsigned int> yMin(LANESIZE, 0);
        bool kernel;
        loadYMin(nodeRun, yMin, kernel);

        if (nodeRun.w0 + nodeRun.w1 <= (maxWeight/rounds)*2+1 && nodeRun.complete){
            TwoRoundTrailCoreStack cacheOrb(nodeRun.stateA, nodeRun.stateB, nodeRun.w0, nodeRun.w1, true, nodeRun.zPeriod);
            TwoRoundTrailCoreCostFunction<Orbital> costFOrb;

            // create yMin
            std::vector<unsigned int> yMin(LANESIZE, 0);
            bool kernel;
            loadYMin(nodeRun, yMin, kernel);

            OrbitalsSet orbSet(kernel, yMin);
            OrbitalTreeIterator iteratorOrb(orbSet, cacheOrb, costFOrb, (maxWeight/rounds)*2+1);

            for (; !iteratorOrb.isEnd(); ++iteratorOrb) {
                if(isInKernel(*iteratorOrb)) {
                    q++;
                    TrailExtension ext(*iteratorOrb);
                    ext.extendForward(maxWeight);
                    ext.removeOutsideKernel();
                    ext.getStats(trailCount, rounds);
                    for(Trail& t: ext.trails) {
                        std::cout << t << std::endl;
                    }

                    TrailExtension ext2(*iteratorOrb);
                    ext2.extendBackward(maxWeight);
                    ext2.removeOutsideKernel();
                    ext2.getStats(trailCount, rounds);
                    for(Trail& t: ext2.trails) {
                        std::cout << t << std::endl;
                    }
                }
            }
        }
    }
    std::cerr << q << std::endl;
    time(&timerEnd);
    std::cerr << "Execution time: " << timerEnd-timerStart << std::endl;
}

void extendFromBruteforce(const unsigned int rounds, const unsigned int maxWeight, std::string filename) {
    time_t timerStart, timerEnd;
    time(&timerStart);

    std::ifstream in(filename);
    if(!in) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    Trail bestTrail;
    unsigned int minWeight = -1; //underflow
    while(in.good()) {
        Trail t, t2;

        if(!t.loadBruteforce(in)) {
            break;
        }

        TrailExtension ext(t);
        for(unsigned int i=0;i<rounds-2;++i)
            ext.extendForward(maxWeight);
        for(unsigned int i=0;i<rounds-2;++i)
            ext.extendBackward(maxWeight);
        ext.getBestTrail(t2,rounds);

        if(t2.getLength() == rounds) {
            std::cout << "Log trail:\n" << t2 << std::endl;
            if(t2.totalWeight > 0 && t2.totalWeight < minWeight && t2.totalWeight <= maxWeight) {
                minWeight = t2.totalWeight;
                bestTrail.set(t2);
                std::cout << "Updating trail to:\n" << bestTrail << std::endl;
            }
        }
    }

    time(&timerEnd);
    std::cerr << "Execution time: " << timerEnd-timerStart << std::endl;
    std::cout << "Best trail:\n" << bestTrail << std::endl;
    in.close();
}

void bruteforceInKernel(const unsigned int maxWeight) {
    time_t timerStart, timerEnd;
    time(&timerStart);

    std::vector<unsigned int> yMin(LANESIZE, 0);
    State a, b;
    TwoRoundTrailCoreStack cacheOrb(a, b, 0, 0, true, LANESIZE);
    OrbitalTreeIterator iteratorOrb(OrbitalsSet(true, yMin), cacheOrb, TwoRoundTrailCoreCostFunction<Orbital>(), maxWeight);
    unsigned int minWeight = -1; //underflows
    unsigned int q = 0;
    ++iteratorOrb;
    for (; !iteratorOrb.isEnd(); ++iteratorOrb) {
            if(isInKernel(*iteratorOrb)) {
                TwoRoundTrailCore nodeOrb = *iteratorOrb;
                q++;
                const unsigned int weight = nodeOrb.w0 + nodeOrb.w1;
                if(weight < minWeight && weight <= maxWeight) {
                    minWeight = weight;
                    applyDispersion(nodeOrb.stateB);
                    nodeOrb.save(std::cout);
                }
            }

    }

    time(&timerEnd);
    std::cerr << q << std::endl;
    std::cerr << "Execution time: " << timerEnd-timerStart << std::endl;
}

int main(int argc, char** argv)
{
    unsigned int rounds = argc > 1 ? atoi(argv[1]) : 5;
    unsigned int maxWeight = argc > 2 ? atoi(argv[2]) : 30;

    //countTrailCores(maxWeight);

    //findTrails(rounds, maxWeight);

    //countInKernelTrails(rounds, maxWeight);

    //extendFromBruteforce(rounds, maxWeight, "0x13_25.log");

    //bruteforceInKernel(maxWeight);

/*
    std::cout << "Found these frequencies:\n";
    for(unsigned int i = 0; i < trailCount.size(); ++i) {
        std::cout << i << " " << trailCount[i] << "\n";
    }
*/
    return 0;
}

