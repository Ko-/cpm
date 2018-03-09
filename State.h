#ifndef _STATE_
#define _STATE_

#include <ostream>

/**
 * Dimensions of the state, used throughout the code.
 * It is assumed that this <= 8*sizeof(unsigned int).
 */
#define LANESIZE 16
#define COLUMNSIZE 4

/**
 * Class to represent a difference of states (differential cryptanalysis).
 */
class State {
private:
    /**
     * The matrix is stored row-wise, for small size and efficient parity computation. The LSB is column 0.
     */
    unsigned int rows[COLUMNSIZE];
public:
    /**
     * Default constructor, initializes to all zero.
     */
    State();
    /**
     * Copy constructor.
     * @param s the state to copy.
     */
    State(const State& s);

    /**
     * Method to set a bit to 1 at a specific location, with the origin in the lower left corner.
     * @param y the y coordinate of the location.
     * @param z the z coordinate of the location.
     */
    void setBit(const unsigned int y, const unsigned int z);
    /**
     * Method to set certain values in a column to 1, based on a bitmask.
     * @param value the bitmask of the rows to set.
     * @param the z coordinate of the column.
     */
    void setColumn(const unsigned int value, const unsigned int z);
    /**
     * Method to set certain values in a column to 1, based on a bitmask.
     * The rest will be set to 0.
     * @param value the bitmask representing the new value of the column.
     * @param the z coordinate of the column.
     */
    void resetColumn(const unsigned int value, const unsigned int z);
    /**
     * Method to set certain values in a column to 0, based on a bitmask.
     * @param value the bitmask of the rows to unset.
     * @param the z coordinate of the column.
     */
    void unsetColumn(const unsigned int value, const unsigned int z);
    /**
     * Method to set all rows to 0.
     */
    void clear();
    /**
     * Method to compute the column parity by xoring all rows.
     * @return the column parity.
     */
    unsigned int getParity() const;
    /**
     * Method to compute the column-wise sums by adding all rows.
     * @param array of length LANESIZE to store results.
     */
    void getSum(unsigned int sum[]) const;
    /**
     * Method to compute the Hamming weight of the difference state, i.e. the number of bits set to 1.
     * @return the hamming weight.
     */
    unsigned int hammingWeight() const;
    /**
     * Method to print one single row as its individual bits. Useful when printing multiple State instances on the same lines.
     * @param out the output stream.
     * @param y index of the row, should be >=0 and <COLUMNSIZE.
     */
    void printRow(std::ostream& out, const int y) const;
    /**
     * Method to rotate a single row to the right by a fixed distance. Used by Permutation.
     * @param y index of the row, should be >=0 and <COLUMNSIZE.
     * @param dist the rotation distance, shoudl be >0 and <COLUMNSIZE.
     */
    void rotateRow(const unsigned int y, const unsigned int dist);

    /**
     * Overload << to output a state object as its individual bits.
     * @param out the output stream.
     * @param state the state to print.
     */
    friend std::ostream& operator<<(std::ostream& out, const State& state);

    /**
     * Overload [] to read a specific row.
     * @param index of the row, should be >=0 and <COLUMNSIZE.
     * @return a const reference to the row with that index.
     */
    const unsigned int& operator[](const int index) const;
    /**
     * A non-const (mutable) version of the [] operator, to modify a specific row.
     * @param index of the row, should be >= 0 and <COLUMNSIZE.
     * @return a reference to the row with that index.
     */
    unsigned int& operator[](const int index);
    /**
     * Overload assignment operator.
     * @param s the state to get values from.
     * @return a reference to the current state.
     */
    State& operator=(const State& s);
};

#endif
