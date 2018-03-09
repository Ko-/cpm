#include "MyTree.h"
#include <iostream>

Orbital OrbitalsSet::getFirstChildUnit(const std::vector<Orbital>& unitList) const
{
    Orbital newOrbital;
    if (unitList.empty()) {
        if (!newOrbital.first(yMin))
            throw EndOfSet();
    }
    else {
        if (!newOrbital.successorOf(unitList.back(), yMin))
            throw EndOfSet();
    }
    return newOrbital;
}

void OrbitalsSet::iterateUnit(const std::vector<Orbital>& unitList, Orbital& current) const
{
    (void)unitList;
    if (!current.next(yMin))
        throw EndOfSet();
}

Order OrbitalsSet::compare(const Orbital& first, const Orbital& second) const
{

    if (first.z < second.z)
        return SMALLER;
    else if (first.z == second.z){
        if (first.y0 < second.y0)
            return SMALLER;
        else if (first.y0 == second.y0){
            if (first.y1 < second.y1)
                return SMALLER;
            else if (first.y1 == second.y1){
                return EQUAL;
            }
        }
    }
    return GREATER;
}

bool OrbitalsSet::isCanonical(const std::vector<Orbital>& orbitalList, TwoRoundTrailCoreStack& cache) const
{
    cache.nodePeriod = LANESIZE;

    if (kernel) {
        if (orbitalList[0].z != 0)
            return false;

        unsigned int lastZ = 0;

        for (unsigned int i = 0; i < orbitalList.size(); i++){
            unsigned int z = orbitalList[i].z;
            if (z != 0 && z > lastZ){ // Consider translation by z only if it has not been already considered before.
                lastZ = z;
                // Translate by z.
                std::vector<Orbital> tauList; // Translated list.
                for (unsigned int j = i; j < orbitalList.size(); j++){
                    Orbital orbital = orbitalList[j];
                    orbital.z -= z;
                    tauList.push_back(orbital);
                }
                for (unsigned int j = 0; j < i; j++){
                    Orbital orbital = orbitalList[j];
                    orbital.z = orbital.z - z + LANESIZE;
                    tauList.push_back(orbital);
                }
                // Compare lists.
                unsigned int k = 0;
                for (k = 0; k < orbitalList.size(); k++){
                    Order cmp = compare(tauList[k], orbitalList[k]);
                    if (cmp == SMALLER)
                        return false; // There is a translated variant smaller than the original.
                    if (cmp == GREATER)
                        break; // Original list is smaller than this translated variant, but still need to check other variants.
                }
                // If the two list are identical, then the list is z-periodic.
                // No need to check the other translated variants.
                if (k == orbitalList.size()){
                    cache.nodePeriod = z;
                    break;
                }
            }
        }
        return true;
    }
    else {
        if (cache.rootPeriod == LANESIZE)
            return true;
        // in the case outside the kernel, only translations by the period must be considered
        for (unsigned int z = cache.rootPeriod; z < LANESIZE; z += cache.rootPeriod){
            // translate by z
            std::vector<Orbital> tauList; // translated list
            unsigned int i = 0;
            for (i = 0; i < orbitalList.size(); i++){
                if (orbitalList[i].z >= z)
                    break;
            }
            for (unsigned int j = i; j < orbitalList.size(); j++){
                Orbital orbital = orbitalList[j];
                orbital.z -= z;
                tauList.push_back(orbital);
            }
            for (unsigned int j = 0; j < i; j++){
                Orbital orbital = orbitalList[j];
                orbital.z = orbital.z - z + LANESIZE;
                tauList.push_back(orbital);
            }
            // compare lists
            unsigned int k = 0;
            for (k = 0; k < orbitalList.size(); k++){
                Order cmp = compare(tauList[k], orbitalList[k]);
                if (cmp == SMALLER)
                    return false; // there is a translated variant smaller than the original
                if (cmp == GREATER)
                    break; // original list is smaller than this translated variant, but still need to check other variants
            }
            // if the two list are identical, then the list is z-periodic
            // no need to check the other translated variants
            if (k == orbitalList.size()){
                cache.nodePeriod = z;
                break;
            }
        }
        return true;
    }
}

const unsigned char ColumnsSet::UOValues[NROFUOVALUES] = {
    0x01, 0x02, 0x04, 0x08 };

const unsigned char ColumnsSet::AEValues[NROFAEVALUES] = {
    0x00, 0x03, 0x05, 0x06, 0x09, 0x0A, 0x0C, 0x0F };

const unsigned char ColumnsSet::AOValues[NROFAOVALUES] = {
    0x01, 0x02, 0x04, 0x07, 0x08, 0x0B, 0x0D, 0x0E };

Column ColumnsSet::getFirstChildUnit(const std::vector<Column>& unitList) const
{
    Column newColumn;

    // the very first column is UOC in (0,0)
    if (unitList.empty()) {
        newColumn.z = 0;
        newColumn.odd = true;
        newColumn.affected = false;
        newColumn.value = UOValues[0];
        newColumn.entangled = false;
    }
    else{
        // there are different cases:
        // UOC -> successor is the AEC
        // AEC -> successor is the UOC
        // AOC -> never treated, because it's the sum of UOC and AEC
        // each time check if it is overlapping (entangled)

        // case UOC -> new column is AEC
        if (!unitList.back().affected && unitList.back().odd){
            newColumn.affected = true;
            newColumn.odd = false;
            newColumn.z = (unitList.back().z + 1) % LANESIZE;
            //if(newColumn.z == 0) //wrapping around TODO
                
            newColumn.value = AEValues[0];
            newColumn.entangled = checkEntanglement(unitList, newColumn);

        }
        // case AEC -> new column is UOC
        else if (unitList.back().affected && !unitList.back().odd){
            
            if(unitList.back().z == 0)
                throw EndOfSet();

            newColumn.affected = false;
            newColumn.odd = true;
            newColumn.value = UOValues[0];
            //if y=0 bit is taken, move to next column
            if(unitList.back().value & 1) {

                if (unitList.back().z + 1 >= LANESIZE)
                    throw EndOfSet();

                newColumn.z = unitList.back().z + 1;
                newColumn.entangled = false;
            }
            else {
                newColumn.z = unitList.back().z;
                newColumn.entangled = true;
                //unitList.back().entangled = true;
            }
        }
        else
            throw EndOfSet();
    }
    return newColumn;
}

bool ColumnsSet::checkEntanglement(const std::vector<Column>& unitList, const Column& current) const
{
    unsigned int i = 0;
    for (i = 0; i < unitList.size(); i++){
        if (current.z == unitList[i].z){
            if ((current.affected && unitList[i].affected) == false)
                return true;
        }
    }

    return false;
}

void ColumnsSet::iterateUnit(const std::vector<Column>& unitList, Column& current) const
{
    // UOC -> iterate column value. It cannot change position.
    // ending AEC -> iterate column value
    //            -> transform it in a UOC to continue the same run
    // starting AEC -> iterate column value
    //              -> change its position to start the run in a new position 

    // UOC
    if (!current.affected && current.odd){
        // if there is an AEC in the same column, the active bit is only in y=0
        if (current.entangled)
            throw EndOfSet();

        if (current.index < NROFUOVALUES-1){
            current.index++;
            current.value = UOValues[current.index];
        }
        else
            throw EndOfSet();
    }

    // AEC
    else if (current.affected && !current.odd){
        if(current.index < NROFAEVALUES - 1) {
            current.index++;
            current.value = AEValues[current.index];
        }
        else
            throw EndOfSet();
    }

    if (unitList.empty() && current.z > 0)
        throw EndOfSet();
}

Order ColumnsSet::compare(const Column& first, const Column& second) const
{
    if (!first.affected && second.affected)
        return SMALLER;
    else if (first.affected == second.affected) {
        if (first.z < second.z)
            return SMALLER;
        else if (first.z == second.z){
            if (first.value < second.value)
                return SMALLER;
            else if (first.value == second.value){
                    return EQUAL;
            }
        }
    }
    return GREATER;
}

bool ColumnsSet::isCanonical(const std::vector<Column>& unitList, TwoRoundTrailCoreStack& cache)  const
{
    cache.nodePeriod = LANESIZE;

    if (cache.rootPeriod == LANESIZE)
        return true;

    if (unitList[0].z != 0)
        return false;
    
    if (unitList.back().odd)
        return true;

    unsigned int lastZ = 0;

    for (unsigned int i = 0; i < unitList.size(); i++){
        unsigned int z = unitList[i].z;
        if (z != 0 && z > lastZ){ // consider translation by z only if it has not been already considered before
            lastZ = z;
            // translate by z
            std::vector<Column> tauList; // translated list
            for (unsigned int j = i; j < unitList.size(); j++){
                Column column = unitList[j];
                column.z = (column.z - z) % LANESIZE;
                tauList.push_back(column);
            }
            for (unsigned int j = 0; j < i; j++){
                Column column = unitList[j];
                column.z = (column.z - z) % LANESIZE;
                tauList.push_back(column);
            }
            // compare lists
            unsigned int j = 0;
            for (j = 0; j < unitList.size(); j++){
                Order cmp = compare(tauList[j], unitList[j]);
                if (cmp == SMALLER) {
                    return false; // there is a translated variant smaller than the original
                }
                if (cmp == GREATER)
                    break; // original list is smaller than this translated variant, but still need to check other variants

            }
            
            // as soon as the original list is equal to a translated variant, then the unit list is z-periodic
            // no need to check the other variants.
            if (j == unitList.size()){
                cache.nodePeriod = z;
                break;
            }
            
        }
    }
    return true;
}

TwoRoundTrailCoreStack::TwoRoundTrailCoreStack()
{
    State emptyState;
    stack_stateAtA.push(emptyState);
    stack_stateAtB.push(emptyState);
    stack_w0.push(0);
    stack_w1.push(0);
    stack_complete.push(true);
    rootPeriod = 0;
    nodePeriod = LANESIZE;
}

TwoRoundTrailCoreStack::TwoRoundTrailCoreStack(State& stateA, State& stateB, const unsigned int aW0, const unsigned int aW1, const bool aComplete, const unsigned int aRootPeriod)
{

    stack_stateAtA.push(stateA);
    stack_stateAtB.push(stateB);
    stack_w0.push(aW0);
    stack_w1.push(aW1);
    stack_complete.push(aComplete);
    rootPeriod = aRootPeriod;
    nodePeriod = LANESIZE;

}

void TwoRoundTrailCoreStack::push(Orbital aOrbital)
{
    State stateAtA = stack_stateAtA.top();
    State stateAtB = stack_stateAtB.top();

    unsigned int new_w0 = stack_w0.top() + 2;
    unsigned int new_w1 = stack_w1.top() + 2;

    stateAtA.setBit(aOrbital.y0, aOrbital.z);
    stateAtA.setBit(aOrbital.y1, aOrbital.z);
    stateAtB.setBit(aOrbital.y0, aOrbital.z);
    stateAtB.setBit(aOrbital.y1, aOrbital.z);

    stack_stateAtA.push(stateAtA);
    stack_stateAtB.push(stateAtB);
    stack_w0.push(new_w0);
    stack_w1.push(new_w1);

    stack_complete.push(true);
}

void TwoRoundTrailCoreStack::push(Column aColumn)
{
    if (aColumn.affected && !aColumn.odd){
        pushAffectedEvenColumn(aColumn);
        stack_complete.push(true);
    }
    else if (!aColumn.affected && aColumn.odd){
        pushUnaffectedOddColumn(aColumn);
        stack_complete.push(false);
    }
}

void TwoRoundTrailCoreStack::pushDummy()
{
    stack_complete.push(stack_complete.top());
    stack_stateAtA.push(stack_stateAtA.top());
    stack_stateAtB.push(stack_stateAtB.top());
    stack_w0.push(stack_w0.top());
    stack_w1.push(stack_w1.top());
    stack_complete.push(true);
}

void TwoRoundTrailCoreStack::pushUnaffectedOddColumn(const Column& aColumn)
{
    State stateAtA = stack_stateAtA.top();
    State stateAtB = stack_stateAtB.top();

    stateAtA.setColumn(aColumn.value, aColumn.z);
    int new_w0 = stack_w0.top() + 1;
    int new_w1 = stack_w1.top();

    if(aColumn.entangled) {
        stateAtB.unsetColumn(aColumn.value, aColumn.z);
        new_w1--;
    }
    else {
        stateAtB.setColumn(aColumn.value, aColumn.z);
        new_w1++;
    }
    
    stack_stateAtA.push(stateAtA);
    stack_stateAtB.push(stateAtB);
    stack_w0.push(new_w0);
    stack_w1.push(new_w1);
}

void TwoRoundTrailCoreStack::pushAffectedEvenColumn(const Column& aColumn)
{
    // for each bit of the column:
    // - if 0: is pushed in B
    // - if 1: is pushed in A
    // the weights are updated
    State stateAtA = stack_stateAtA.top();
    State stateAtB = stack_stateAtB.top();

    stateAtA.setColumn(aColumn.value, aColumn.z);
    stateAtB.setColumn(aColumn.inverseValue(), aColumn.z);

    int delta0 = 0;

    for (unsigned int y = 0; y<COLUMNSIZE; y++) {
        if (((aColumn.value >> y) & 1) != 0)
            delta0 += 1;
    }
    

    int new_w0 = stack_w0.top() + delta0;
    int new_w1 = stack_w1.top() + (COLUMNSIZE - delta0);

    stack_stateAtA.push(stateAtA);
    stack_stateAtB.push(stateAtB);
    stack_w0.push(new_w0);
    stack_w1.push(new_w1);
}

void TwoRoundTrailCoreStack::pop(){

    stack_stateAtA.pop();
    stack_stateAtB.pop();
    stack_w0.pop();
    stack_w1.pop();
    stack_complete.pop();

}

void TwoRoundTrailCoreStack::save(std::ostream& fout){
   fout << "At A:\n" << stack_stateAtA.top();
   fout << "\nAt B:\n" << stack_stateAtB.top();
   fout << "\nWith weight: " << stack_w0.top() + stack_w1.top() << std::endl;
}

template<class Unit>
unsigned int TwoRoundTrailCoreCostFunction<Unit>::getCost(const TwoRoundTrailCoreStack& cache) const
{
    return cache.stack_w0.top() + cache.stack_w1.top();
}

template <> 
bool TwoRoundTrailCoreCostFunction<Column>::canAfford(const std::vector<Column>& unitList, const TwoRoundTrailCoreStack& cache, const Column& newUnit, const unsigned int maxCost) const
{
    unsigned int Gamma = getCost(cache);

    if(newUnit.affected) //AOC, AEC
        Gamma += 4;
    else if(newUnit.odd) { //UOC
        if(!unitList.empty() && newUnit.z != unitList.back().z) //do not increase cost for UOC in same column as previous AEC
            Gamma += 2;
    }
    else  { //UEC, should not occur
        return false;
    }
    return Gamma <= maxCost;
}

template <> 
bool TwoRoundTrailCoreCostFunction<Orbital>::canAfford(const std::vector<Orbital>& unitList, const TwoRoundTrailCoreStack& cache, const Orbital& newUnit, const unsigned int maxCost) const
{
    (void)unitList;
    (void)newUnit;
    unsigned int Gamma = getCost(cache);
    
    Gamma += 4;

    return Gamma <= maxCost;
}

template<class Unit>
bool TwoRoundTrailCoreCostFunction<Unit>::canAfford(const std::vector<Unit>& unitList, const TwoRoundTrailCoreStack& cache, const Unit& newUnit, const unsigned int maxCost) const
{
    //Only occurs when Unit is not Orbital or Column, which should never happen
    return false;
}

void TwoRoundTrailCore::set(const TwoRoundTrailCoreStack& cache)
{
    stateA = cache.stack_stateAtA.top();
    stateB = cache.stack_stateAtB.top();
    w0 = cache.stack_w0.top();
    w1 = cache.stack_w1.top();
    complete = cache.stack_complete.top();
    zPeriod = cache.nodePeriod;
}

void TwoRoundTrailCore::save(std::ostream& fout)
{
   fout << "At A:\n" << stateA;
   fout << "\nAt B:\n" << stateB; 
   fout << "\nWith weight: " << w0 + w1 << std::endl;
}


