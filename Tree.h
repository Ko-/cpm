/*! \file Tree.h

This file defines the components for a general tree structure.
*/

#ifndef _TREE_
#define _TREE_

#include <vector>
#include "Exception.h"
#include "Units.h"
#include <iostream>

/**
* This exception is launched when a set of units reaches the end.
*/
class EndOfSet : Exception {};

/**
* \brief GenericTreeIterator class : Iterator to traverse a Tree.
*
* \details This class represents an iterator to traverse a tree.
* The type of tree is defined by the unitList.
*/
template<class Unit, class UnitSet, class CachedRepresentation, class OutputRepresentation, template<class> class CostFunction>
class GenericTreeIterator {
protected:
    /** The set of units. */
    const UnitSet& unitSet;
    /** The unit-list representing the current node. */
    std::vector<Unit> unitList;
    /** The cache representation of the current node. */
    CachedRepresentation cache;
    /** The output representation of the current node. */
    OutputRepresentation out;
    /** The cost function. */
    const CostFunction<Unit>& costFunction;
    /** The cost of the current node. */
    std::vector<unsigned int> cost;
    /** The maximum cost allowed when traversing the tree. */
    unsigned int maxCost;

    unsigned int rootPeriod, nodePeriod;
    
    bool end, initialized, empty;
    unsigned long long index;
    //const bool combineTwoUnits = std::is_same<Unit, Column>::value;
public:

     /** The constructor.
      * @param  aUnitSet The set of units.
      * @param  aCache The cache representation.
      * @param  aCostFunction The cost function. 
      * @param  aMaxCost The maximum cost.
      */
    GenericTreeIterator(const UnitSet& aUnitSet, CachedRepresentation aCache, const CostFunction<Unit>& aCostFunction, unsigned int aMaxCost)
        : unitSet(aUnitSet), cache(aCache), costFunction(aCostFunction), maxCost(aMaxCost), rootPeriod(cache.rootPeriod), nodePeriod(cache.nodePeriod)
    {
        empty = true;
        end = false;
        initialized = false;
        index = 0;
    }

    /** This method indicates whether the iterator has reached the end of the tree.
    *  @return True if there are no more nodes to consider.
    */
    bool isEnd()
    {
        //if (!initialized) initialize();
        return end;
    }

    /** This method indicates whether the tree is empty.
    *  @return True if the tree is empty.
    */
    bool isEmpty()
    {
        if (!initialized) initialize();
        return empty;
    }

    /** The '++' operator moves the iterator to the next node in the tree.
    */
    void operator++()
    {
        if (!initialized) {
            initialize();
        }
        else {
            if (!end) {
                ++index;
                if (!next())
                    end = true;
            }
        }
    }

    /** The '*' operator gives a constant reference to the current node.
    *  @return A constant reference to the current node.
    */

    OutputRepresentation& operator*()
    {
        out.set(/*unitList,*/ cache);
        return out;
    }
    
    OutputRepresentation* operator->()
    {
        out.set(/*unitList,*/ cache);
        return &out;
    }

private:

    /** This method initializes the iterator.
    */
    void initialize()
    {
        index = 0;
        if (first()) {
            end = false;
            empty = false;
        }
        else {
            end = true;
            empty = true;
        }
        initialized = true;
    }

    /** This method returns the first node of the tree.
    * @return true if the first node exists, false otherwise.
    */
    bool first()
    {
        return toChild();
    }

    /** This method returns the next node of the tree.
    * @return true if the next node exists, false otherwise.
    */
    bool next()
    {
        if (toChild())
            return true;
        do{
            if (toSibling())
                return true;
            if (!toParent())
                return false;
        } while (true);
    }

    /** This method move to the first child of the current node.
    * A child is obtained adding a new unit to the current node.
    * @return true if a child is found, false otherwise.
    */
    bool toChild()
    {
        try {
            Unit newUnit = unitSet.getFirstChildUnit(unitList);
            if (canAfford(newUnit)) {
                push(newUnit);
                if (cost.back() <= maxCost && isCanonical())
                    return true;
                else{
                    if (iterateHighestUnit())
                        return true;
                    else{
                        pop();
                        return false;
                    }
                }
            }
            else{
                return false;
                // With constant costs, this case is irrelevant.
            }
        }
        catch (EndOfSet) {
            return false;
        }
    }

    /** This method move to the first sibling of the current node.
    * Siblings are obtained iterating the higher unit of the current node.
    * @return true if a sibling is found, false otherwise.
    */
    bool toSibling()
    {
        if (unitList.empty())
            return false;
        else
            return iterateHighestUnit();
    }

    /** This method move to the parent of the current node.
    * The parent is obtained by removing the higher unit of the current node.
    * @return true if the node has a parent, false if it is the root.
    */
    bool toParent()
    {
        if (unitList.empty())
            return false;
        else
            return pop();
    }

    /**
    * This method iterates the highest unit based on the order relation defined by the unitSet.
    * @return true if a valid value for the highest unit is found, false otherwise.
    */
    bool iterateHighestUnit()
    {
        Unit lastUnit = unitList.back();
        pop();
        do {
            do {
                try {
                    unitSet.iterateUnit(unitList, lastUnit);
                }
                catch (EndOfSet) {
                    pushDummy(lastUnit); // Something to pop is needed by the function toParent().
                    return false;
                }
            } while (!canAfford(lastUnit));

            push(lastUnit);
            if (cost.back() <= maxCost && isCanonical())
                return true;
            pop();
        } while (true);
    }

    /**
    * This method pushes a new unit to the unit list and updates the cost function.
    * @param newUnit the unit to be pushed.
    */
    void push(const Unit& newUnit)
    {
        unitList.push_back(newUnit);
        cache.push(newUnit);
        cost.push_back(costFunction.getCost(/*unitList,*/ cache));
    }

    /**
    * This method pushes a dummy unit.
    * @param newUnit the dummy unit to be pushed.
    */
    void pushDummy(const Unit& newUnit)
    {
        unitList.push_back(newUnit);
        cache.pushDummy();
        cost.push_back(0);
    
    }

    /**
    * This method pops the highest unit from the unit list and updates the cost function.
    * @return true if a unit is popped, false if the unit list is empty.
    */
    bool pop()
    {
        if (unitList.empty())
            return false;
        unitList.pop_back();
        cache.pop();
        cost.pop_back();
        return true;
    }

    /** This method checks if the current node is canonical
    * with respect to an order relation specified by the unitSet.
    * @return true if the current node is canonical, false otherwise.
    */
    bool isCanonical()
    {
        return unitSet.isCanonical(unitList, cache);
    }

    /** This method checks if the addition of a new unit is affordable.
    * @return true if it is affordable, false otherwise.
    */
    bool canAfford(Unit& newUnit)
    {
        return costFunction.canAfford(unitList, cache, newUnit, maxCost);
    }

};

#endif
