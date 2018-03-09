#include "Units.h"
#include "State.h"

bool Orbital::first(const std::vector<unsigned int>& yMin)
{
    z = 0;
    y0 = yMin[z];
    while(y0 >= COLUMNSIZE-1) {
        if (z < LANESIZE-1)
            z++;
        else
            return false;
        y0 = yMin[z];
    }
    y1 = y0 + 1;
    return true;
}

bool Orbital::next(const std::vector<unsigned int>& yMin)
{
    if (y1 < COLUMNSIZE-1)
        y1++;
    else {
        if (y0 < COLUMNSIZE-2) {
            y0++;
            y1 = y0 + 1;
        }
        else {
            do {
                if (z < LANESIZE-1)
                    z++;
                else
                    return false;
                y0 = yMin[z];
            } while(y0 >= COLUMNSIZE-1);
            y1 = y0 + 1;
        }
    }
    return true;
}

bool Orbital::successorOf(const Orbital& other, const std::vector<unsigned int>& yMin)
{
    z = other.z;
    y0 = other.y1 + 1;
    while(y0 >= COLUMNSIZE-1) {
        if (z < LANESIZE-1)
            z++;
        else
            return false;
        y0 = yMin[z];
    }
    y1 = y0 + 1;
    return true;
}

bool operator<(const Orbital& aCP, const Orbital& bCP)
{
    if(aCP.z != bCP.z)
        return aCP.z < bCP.z;
    if(aCP.y0 != bCP.y0)
        return aCP.y0 < bCP.y0;
    return aCP.y1 < bCP.y1;
}

std::ostream& operator<<(std::ostream& fout, const Orbital& aCP)
{
    fout << "(" << aCP.z << ",(" << aCP.y0 << "," << aCP.y1 << "))";
    return fout;
}

Column::Column(){
    z = 0;
    value = 0;
    index = 0;
    odd = false;
    affected = false;
    entangled = false;
}

Column::Column(const bool aOdd, const bool aAffected){

    z = 0;
    value = 0;
    index = 0;    
    odd = aOdd;
    affected = aAffected;
    entangled = false;
}

unsigned char Column::inverseValue() const {
    return (unsigned char)(((1 << COLUMNSIZE) - 1) ^ value);
}

std::ostream& operator<<(std::ostream& fout, const Column& aColumn)
{
    fout << "(" << aColumn.z << "," << (unsigned int)aColumn.value << ")";
    return fout;
}
