/*! \file MyTree.h

This file defines Orbital-Tree and Run-Tree, used to generate 2-round trail cores in Keccak-f,
as described in the paper.
*/

#ifndef _MYTREE_
#define _MYTREE_

#include <ostream>
#include <stack>
#include <vector>
#include "State.h"
#include "Tree.h"
#include "Units.h"

enum Order {EQUAL, SMALLER, GREATER};

class TwoRoundTrailCoreStack;

/**
* \brief OrbitalsSet class : represents the set of orbitals.
* \details This class represents the set of orbitals and defines the order relation among them.
*/
class OrbitalsSet {
public:
    /** This attribute indicates whether the orbitals are used to generate states in the kernel or not. */
    bool kernel;
    /** The minimum position of the lower bit of an orbital for each column of the state. */
    std::vector<unsigned int> yMin;
    
public:
    /**
    * The default constructor.
    */
    OrbitalsSet()
        : kernel(true), yMin(LANESIZE, 0) {}

    /**
    * The constructor.
    * @param ayMin the minimum position of the lower bit of an orbital for each column of the state.
    */
    OrbitalsSet(bool aKernel, const std::vector<unsigned int>& ayMin)
        : kernel(aKernel), yMin(ayMin) {}

    /** This method returns an orbital in the first available position
    * with respect to the order relation [z,x,y0,y1] and restrictions given by yMin.
    * @param unitList the list of units.
    * @return the first available orbital position.
    */
    Orbital getFirstChildUnit(const std::vector<Orbital>& unitList) const;
    
    /** This method iterates the current orbital with respect to the order relation [z,x,y0,y1] and restrictions given by yMin.
    * @param unitList the list of units.
    * @param current the current orbital.
    */
    void iterateUnit(const std::vector<Orbital>& unitList, Orbital& current) const;
    
    /** This method compares two given orbitals with respect to the order relation [z,y0,y1].
    * @param first the first given orbital.
    * @param second the second given orbital.
    * @return 0 if the two orbitals are equal, 1 if the first is smaller, 2 if the second is smaller.
    */
    Order compare(const Orbital& first, const Orbital& second) const;
    
    /** This method checks if a list of orbitals is z-canonical with respect to the order relation [z,y0,y1].
    * @param orbitalList the list of orbitals.
    * @param rootPeriod the z-period of the root of the tree.
    * @param nodePeriod the z-period of the current node of the tree.
    * @return true if the given list is z-canonical, false otherwise.
    */
    bool isCanonical(const std::vector<Orbital>& orbitalList, TwoRoundTrailCoreStack& cache) const;
};


/**
* \brief ColumnsSet class : represents the set of column assignments.
* \details This class represents the set of column assignments and defines the order relation among them.
*/
class ColumnsSet {

public:
    /** The possible values of Affected Even columns, Affected Odd columns and Unaffected Odd columns. */
    static const unsigned int NROFAEVALUES = 8, NROFAOVALUES = 8, NROFUOVALUES = 4;
    static const unsigned char AEValues[NROFAEVALUES], AOValues[NROFAOVALUES], UOValues[NROFUOVALUES];

public:
    
    /**
    * The default constructor.
    */
    ColumnsSet() { };
    
    /** This method returns a column assignment in the first available position
    * with respect to the order relation [z].
    * @param unitList the list of units.
    * @return the first available column assignment.
    */
    Column getFirstChildUnit(const std::vector<Column>& unitList) const;
    
    /** This method iterates the current unit with respect to the order relation [z,value].
    * @param unitList the list of units.
    * @param current the current column assignment.
    */
    void iterateUnit(const std::vector<Column>& unitList, Column& current) const;
    
    /** This method compares two given column assignments with respect to the order relation [z,value].
    * @param first the first given column assignment.
    * @param second the second given column assignment.
    * @return 0 if the two column assignments are equal, 1 if the first is smaller, 2 if the second is smaller.
    */
    Order compare(const Column& first, const Column& second) const;
    
    /** This method checks if a list of column assignments is z-canonical with respect to the order relation [z,value].
    * @param unitList the list of column assignments.
    * @param cache a reference to the cache representing the trail.
    * @return true if the given list is z-canonical, false otherwise.
    */
    bool isCanonical(const std::vector<Column>& unitList, TwoRoundTrailCoreStack& cache) const;
    
    /** This method checks if a given column assignment overlaps any other column of a given unit list.
    * @param unitList the list of column assignments.
    * @param current the given column assignment.
    * @return true if the given column overlaps a column in @unitList.
    */
    //bool checkColumnOverlapping(const std::vector<Column>& unitList, Column& current) const;
    
    /** This method checks if a given column assignment is entangled with any other column of a given unit list.
    * @param unitList the list of column assignments.
    * @param current the given column assignment.
    * @return true if the given column is entangled with a column in @unitList.
    */
    bool checkEntanglement(const std::vector<Column>& unitList, const Column& current) const;

};


/**
* \brief TwoRoundTrailCoreStack class : cache representation for 2-round trail cores in Keccak-f.
* \details This class represents a 2-round trail core as a node of a tree. 
*/
class TwoRoundTrailCoreStack {
public:
    /** The stack for state at A and state at B. 
    * The first elements of the stacks refer to a state with no units.
    * Each element of the stack is equal to its predecessor plus a new unit.
    */
    std::stack<State> stack_stateAtA, stack_stateAtB; 
    /** The stack for the minimum reverse weight of A and the weight of B. */
    std::stack<unsigned int> stack_w0, stack_w1;
    /** The stack to indicate if a state is valid or not. */
    std::stack<bool> stack_complete;
    /** The z-period of the root and the current node. */
    unsigned int rootPeriod, nodePeriod;
    
public:

    /**
    * The default constructor.
    */
    TwoRoundTrailCoreStack();
    
    /**
    * The constructor.
    * @param aDCorLC The propagation context of the trail, as a reference to a KeccakFPropagation object.
    */
//    TwoRoundTrailCoreStack(const KeccakFPropagation& aDCorLC);

    /**
    * The constructor.
    * @param aDCorLC The propagation context of the trail, as a reference to a KeccakFPropagation object.
    * @param stateA The initial state at A.
    * @param stateB The initial state at B.
    * @param aW0 The initial minimum reverse weight of A.
    * @param aW1 The initial weight of B.
    * @param aComplete Completness of the initial state.
    * @param aRootPeriod The z-period of the root.
    */
    TwoRoundTrailCoreStack(State& stateA, State& stateB, const unsigned int aW0, const unsigned int aW1, const bool aComplete, const unsigned int aRootPeriod);

    /**
    * This method pushes an orbital to the cache.
    * @param aOrbital the orbital to be pushed.
    */
    void push(Orbital aOrbital);

    /**
    * This method pushes a column to the cache.
    * @param aColumn the column to be pushed.
    */
    void push(Column aColumn);

    /**
    * This method performs a dummy push to the cache.
    */
    void pushDummy();

    /**
    * This method pops the highest unit from the cache.
    */
    void pop();

    /** This methods outputs the cache to save it in, e.g., a file.
    * @param fout The stream to save the cache to.
    */
    void save(std::ostream& fout);

    /** This methods returns the value of a row in position (y,z) in A.
    * @return The row value.
    */
/*TODO    RowValue getRowA(unsigned int y, unsigned int z) const;
*/    
    /** This methods returns the value of a row in position (y,z) in B.
    * @return The row value.
    */
/*TODO    RowValue getRowB(unsigned int y, unsigned int z) const;
*/
private:
    // method that pushes an unaffected odd column
    void pushUnaffectedOddColumn(const Column& aColumn);
    // method that pushes an affected even column
    void pushAffectedEvenColumn(const Column& aColumn);
    // method that sets a bit to one and returns the difference of minimum reverse weight before and after the addition
//TODO    int pushBitAndGetDeltaMinReverseWeight(vector<SliceValue>& state, const BitPosition& p);
    // method that sets a bit to one and returns the difference of weight before and after the addition
//TODO    int pushBitAndGetDeltaWeight(vector<SliceValue>& state, const BitPosition& p);

};

/**
* \brief TwoRoundTrailCoreCostFunction class : the cost of a 2-round trail core in Keccak-f.
* 
* This class represents the cost function for 2-round trail cores in Keccak-f.
* The cost function is defined as: \alpha*w0+\beta*w1.
* Examples are: w0, w1, w0+w1, w0+2w1, 2w0+w1.
*/
template<class Unit>
class TwoRoundTrailCoreCostFunction {

public:
    /** This method returns the cost of a given node.
    * @param unitList The list of units representing the given node.
    * @param cache The cache representation of the given node.
    * @return The cost of the given node.
    */
    unsigned int getCost(/*const std::vector<Unit>& unitList,*/ const TwoRoundTrailCoreStack& cache) const;
    
    /** This method checks whether the addition of a new unit to a given node is within a given budget.
    * @param unitList The list of units representing the given node.
    * @param cache The cache representation of the given node.
    * @newUnit The unit to be added.
    * @maxCost The given budget.
    * @return True iff the current cost + cost of the new unit are below or equal to maxCost
    */
    bool canAfford(const std::vector<Unit>& unitList, const TwoRoundTrailCoreStack& cache, const Unit& newUnit, const unsigned int maxCost) const;
};

/**
* \brief TwoRoundTrailCore class : output representation for 2-round trail cores in Keccak-f.
* \details This class represents a 2-round trail core in Keccak-f.
*/
class TwoRoundTrailCore{

public:
    /** The state at A. */
    State stateA; 
    
    /** The state at B. */
    State stateB;
    
    /** The minimum reverse weight of A. */
    unsigned int w0;
    
    /** The weight of B. */
    unsigned int w1;
    
    /** The 2-round trail core. */
    //Trail trail;
    
    /** The completness of the trail. */
    bool complete;
    
    /** The z-period of the trail. */
    unsigned int zPeriod;

public:

    /** The default constructor. */
    TwoRoundTrailCore()
        : stateA(), stateB(), w0(0), w1(0), complete(true), zPeriod(LANESIZE) {}

    /** This methods outputs the 2-round trail core to save it in, e.g., a file.
    * @param fout The stream to save the trail core to.
    */
    void save(std::ostream& fout);

    /** This methods sets the trail.
    * @param cache The cache representation of the trail.
    */
    void set(const TwoRoundTrailCoreStack& cache);
};

/** The orbital tree. */
typedef GenericTreeIterator<Orbital, OrbitalsSet, TwoRoundTrailCoreStack, TwoRoundTrailCore, TwoRoundTrailCoreCostFunction> OrbitalTreeIterator;

/** The run tree. */
typedef GenericTreeIterator<Column, ColumnsSet, TwoRoundTrailCoreStack, TwoRoundTrailCore, TwoRoundTrailCoreCostFunction> RunTreeIterator;

#endif
