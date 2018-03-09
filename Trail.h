#ifndef _TRAIL_
#define _TRAIL_

#include <iostream>
#include <vector>
#include "Exception.h"
#include "Permutation.h"
#include "State.h"

class TrailException : Exception {
public:
    explicit TrailException(const std::string& aReason) : Exception(aReason) {};
};

/** This class implements a container for a differential or linear trail.
  * The main attribute of the class is the sequence of state values before χ.
  * If S<sub>i</sub> = @a states[i], the trail is:
  * S<sub>0</sub> χλ S<sub>1</sub> χλ S<sub>2</sub> χλ ...
  * χλ S<sub><i>n</i>-1</sub>,
  * with <i>n</i> = @a states.size().
  */
class Trail {
public:
    /** This attribute contains the list of states round after round, before χ.
      */
    std::vector<State> states;
    /** This attribute contains the propagation weights of the states
      * in @a states. So, @a weights has the same size
      * as @a states.
      */
    std::vector<unsigned int> weights;
    /** This attribute contains the sum of the weights[i].
      */
    unsigned int totalWeight;

    /** This attribute contains the index of the first half of the initial two-round trail core in states/weights
     */
    unsigned int initialIndex;
public:
    /** This constructor creates an empty trail.
    */
    Trail() : states(), weights(), totalWeight(0), initialIndex(0) {}
    /** This constructor loads a trail from an input stream.
      * @param   fin    The input stream to read the trail from.
      */
    explicit Trail(std::istream& fin);
    /** This constructor initializes a trail by copying the
      * trail given in parameter.
      * @param   other  The original trail to copy.
      */
    Trail(const Trail& other) :
        states(other.states),
        weights(other.weights),
        totalWeight(other.totalWeight),
        initialIndex(other.initialIndex) {}
    /** This method returns the length of the trail.
      *  @return    The number of rounds.
      */
    unsigned int getLength() const;
    /** This method empties the trail.
      */
    void clear();
    /** This method appends a state to the end of @a states,
      * with its corresponding propagation weight.
      * @param   state  The state to add.
      * @param   weight The propagation weight.
      */
    void append(const State& state, unsigned int weight);
    /** This method appends another trail to the current trail.
      * @param   otherTrail The trail to append.
      */
/*    void append(const Trail& otherTrail);*/
    /** This method inserts a state at the beginning of @a states,
      * with its corresponding propagation weight.
      * @param   state  The state to add.
      * @param   weight The propagation weight.
      */
    void prepend(const State& state, unsigned int weight);
    /** This method removes the first state from the current trail. This does not check if the Trail is empty, so use with care.
      */
    void prepop();
    /** This method removes the last state from the current trail. This does not check if the Trail is empty, so use with care.
      */
    void pop();
    /** This method prunes a trail from at most 2(n-1) rounds to n, but only if this actually results in a trail with lower weight than the current minimum.
      * @param rounds The number of rounds (n) to prune to.
      * @param bestMinimumWeight The best minimum weight so far.
      * @return Whether the trail still needs to be considered.
      */
    bool prune(const unsigned int rounds, const unsigned int bestMinimumWeight);
    /** This method sets a trail object to another object. Prevents the creation of new objects for each copy.
      * @param otherTrail the trail to set the values to.
      */
    void set(const Trail& otherTrail);

    /** This method displays the trail for in a human-readable form.
      * @param  out    The stream to display to.
      * @param state   The state to display.
      * @return A reference to the stream object.
      */
    friend std::ostream& operator<<(std::ostream& out, const Trail& trail);
    /** Implements the assignment operator for Trail objects.
     * @param t the trail to get values from.
     * @return a reference to the current trail.
     */
    Trail& operator=(const Trail& t);
    /** This methods loads the trail in a binary encoding from a stream (e.g., file).
      * @param   fin    The input stream to read the trail from.
      */
    void load(std::istream& fin);
    /** This methods outputs the trail in a binary encoding and saves it to, e.g., a file.
      * @param  fout    The stream to save the trail to.
      */
    void save(std::ostream& fout) const;
    /** This methods loads the trail in an alternate binary encoding defined by bruteforce.cpp.
      * @param   fin    The input stream to read the trail from.
      */
    bool loadBruteforce(std::istream& fin);
};

#endif
