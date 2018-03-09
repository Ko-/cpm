#ifndef _PERMUTATION_
#define _PERMUTATION_

#include <ostream>
#include <vector>
#include "State.h"

/**
 * Returns the effect vector (as integer) based on a parity profile, in this case by rotating it one to the right.
 * @param parity the given parity profile.
 * @return the effect vector as unsigned integer.
 */
unsigned int getThetaEffect(const unsigned int parity);

/**
 * Applies the full theta effect to a given state. I.e. compute the column parity, rotate, and xor back into the state.
 * @param state the state to be modified.
 */
void applyTheta(State& state);

/**
 * Applies a given theta effect to a given state. This does not have to be the normal column parity rotated to the right, although it could be.
 * @param state the state to be modified.
 * @return the given theta effect.
 */
void applyTheta(State& state, const unsigned int effect);

/**
 * Applies the dispersion layer to a given state.
 * @param state the state to be modified.
 */
void applyDispersion(State& state);

/**
 * Applies the inverse dispersion layer to a given state.
 * @param state the state to be modified.
 */
void applyInverseDispersion(State& state);

/**
 * Print a state to which the dispersion layer has been applied to a given ostream, but do not actually alter the state.
 * @param out the output stream to print to.
 * @param state the state.
 */
void printDispersed(std::ostream& out, const State& state);

/**
 * When extending a trail forward by one round, two types of branching occur in the theta step, parity branching and effect branching. This function collects all possible results after theta in a vector.
 * @param state the state to start with (at end of trail, at nonlinear step).
 * @param compatibleStates the vector to fill, empty at the start.
 * @param minWeight the least active S-boxes that can be spent on this round.
 * @param maxWeight the most active S-boxes that can be spent on this round.
 *
 */
void thetaCompatibleStates(const State& state, std::vector<State>& compatibleStates, const unsigned int minWeight = 1, const unsigned int maxWeight = LANESIZE*COLUMNSIZE);

/**
 * A recursive helper function for effect branching, only to be used by thetaCompatibleStates.
 * @param state intermediate result from thetaCompatibleStates, after parity branching and applying the theta effect.
 * @param compatibleStates the vector to fill.
 * @param effect the calculated theta effect, based on the original state before applying the theta effect.
 * @param sum an array containing the sum of all state rows, i.e. how many active S-boxes per column, based on the original state before applying the theta effect.
 * @param colIndex the current column.
 */
void recurThetaCompatibleStates(State& state, std::vector<State>& compatibleStates, const unsigned int effect, const unsigned int sum[LANESIZE], const unsigned int colIndex = 0);


#endif
