#ifndef _UNITS_
#define _UNITS_

#include <ostream>
#include <vector>

class Orbital
{
public:
    /** The y-coordinate of the first bit in the orbital. */
    unsigned int y0;
    /** The y-coordinate of the second bit in the orbital. */
    unsigned int y1;
    /** The z-coordinate of the column, 0 ≤ z < LANESIZE. */
    unsigned int z;
public:
    /** The default constructor. */
    Orbital() : y0(0), y1(1), z(0) {}
    
    Orbital(const unsigned int ay0, const unsigned int ay1, const unsigned int az) : y0(ay0), y1(ay1), z(az) {}
    
    /** This method sets the position to the first available one with y-coordinates
      * at least as specified by the parameter yMin.
      * @param  yMin    The value in yMin[z] indicates the minimum y-coordinate of the bits in the orbital.
      */
    bool first(const std::vector<unsigned int>& yMin);
    
    /** This method sets the position to the next available one with y-coordinates
      * at least as specified by the parameter yMin.
      * @param  yMin    The value in yMin[z] indicates the minimum y-coordinate of the bits in the orbital.
      */
    bool next(const std::vector<unsigned int>& yMin);
    
    /** This method sets the position to the next available orbital after @other with y-coordinates
      * at least as specified by the parameter yMin and with y-coordinates higher than those
      * of the specified orbital @other if in the same column.
      * @param  other   The other orbital position.
      * @param  yMin    The value in yMin[z] indicates the minimum y-coordinate of the bits in the orbital.
      */
    bool successorOf(const Orbital& other, const std::vector<unsigned int>& yMin);
    
    /** An ordering operator, required when storing a Orbital object
      * in a set or as the first member in maps.
      * @param  a   The orbital position at the left of the operator.
      * @param  az  The orbital position at the right of the operator.
      */
    friend bool operator<(const Orbital& aCP, const Orbital& bCP);
    
    /** A display function, for use with the << operator.
      * @param  out The output stream to write to.
      * @param  aCP  The orbital position to display.
      */
    friend std::ostream& operator<<(std::ostream& out, const Orbital& aCP);
};

/**
* \brief Column class : represents column assignments.
* \details This class represents a column assignment for Keccak-f.
*/
class Column {

public:
    /** The z coordinate of the column. */
    unsigned int z;
    /** The value of the column.
    * The type is one byte, containing the COLUMNSIZE
    * bits of a column, in the least significant bits of the byte.
    */
    unsigned char value;
    /** The parity of the column. */
    bool odd;
    /** The theta-effect on the column. */
    bool affected;
    /** Entanglement of the column. */
    bool entangled;
    /** Parameter used to optimize the iteration of the column value. */
    unsigned int index;

public:

    /** The default construtor. */
    Column();

    /** The constructor. 
    * @param odd the parity of the column.
    * @param affected the theta-effect on the column.    
    */
    Column(const bool odd, const bool affected);

    unsigned char inverseValue() const;
    
    /** A display function, for use with the << operator.
      * @param  out The output stream to write to.
      * @param  aColumn  The column to display.
      */
    friend std::ostream& operator<<(std::ostream& out, const Column& aColumn);

};

#endif
