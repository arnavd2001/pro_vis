#ifndef SOLUTION_H_
#define SOLUTION_H_

/** \file solution.h Routines for manipulating Solution objects, such as creation, randomization, perturbation etc. */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fitness/fitness.h>
#include <shiftmel.h>
#include <random.h>

#define FITNESS_MIN -1E9

#ifndef SOLUTION_SOURCE_CODE
	#define SOLUTION_INLINE inline
#else
	#define SOLUTION_INLINE extern inline
#endif

#include "solution_structure_private.h"

typedef struct Solution_ Solution;

/** Returns a Solution whose fields are all uninitialized, but with due memory allocated. */
SOLUTION_INLINE
Solution Solution_blank(int hpSize){
	Solution retval;
	retval.chain = malloc(sizeof(shiftmel) * (hpSize - 1));
	retval.fitness = FITNESS_MIN;
	retval.idle_iterations = 0;
	return retval;
}

/** Returns a deep copy (all memory recursively duplicated) of the given solution. */
SOLUTION_INLINE
Solution Solution_copy(Solution sol, int hpSize){
	Solution retval;
	retval.fitness = sol.fitness;
	retval.idle_iterations = sol.idle_iterations;

	int chainSize = hpSize - 1;

	retval.chain = malloc(sizeof(shiftmel) * chainSize);
	memcpy(retval.chain, sol.chain, sizeof(shiftmel) * chainSize);

	return retval;
}

/** Frees memory allocated for given solution */
SOLUTION_INLINE
void Solution_free(Solution sol){
	free(sol.chain);
}

/** Returns a Solution whose movement chain is uniformly random.
 * The returned Solution has its idle_iterations set to 0.
 * The returned Solution won't have its fitness calculated.
 */
SOLUTION_INLINE
Solution Solution_random(int hpSize){
	Solution sol;
	int nMovements = hpSize - 1;

	sol.idle_iterations = 0;

	// Generate random shiftmel *
	sol.chain = malloc(sizeof(shiftmel) * nMovements);
	int i;
	for(i = 0; i < nMovements; i++)
		sol.chain[i] = shiftmel_random();

	sol.fitness = FITNESS_MIN;

	return sol;
}

/** Chooses a random element ELEM1 in 'perturb'.
 * Then chooses a random element ELEM2 in 'other'.
 * Takes the distance DIST between ELEM1 and ELEM2
 * Changes 'perturb' so that its ELEM1 approaches ELEM2 by a random amount, from 0 to 100%.
 * The solution 'perturb' is returned.
 *
 * The returned Solution has its idle_iterations set to 0.
 * The returned Solution won't have its fitness calculated.
 */
SOLUTION_INLINE
Solution Solution_perturb_relative(Solution perturb, Solution other, int hpSize){
	int chainSize = hpSize - 1;
	int pos1 = urandom_max(chainSize);
	int pos2 = urandom_max(chainSize);

	pos2 = pos1;

	Solution retval = Solution_copy(perturb, hpSize);
	unsigned char elem1 = shiftmel_to_number(retval.chain[pos1]);
	unsigned char elem2 = shiftmel_to_number(other.chain[pos2]);

	char distance = elem2 - (char) elem1;

	// Generate a number in [0, distance)
	double aux = drandom_x() * (double) distance;

	// Fit the number in the discrete space [0, distance]
	char delta = (char) round(aux);

	retval.chain[pos1] = shiftmel_from_number(elem1 + delta);
	retval.idle_iterations = 0;
	retval.fitness = FITNESS_MIN;

	return retval;
}

/** Returns the fitness of the given solution, calculating it only if needed.
 * \return The fitness of `sol`.
 */
SOLUTION_INLINE
double Solution_fitness(Solution sol){
	if(sol.fitness < (FITNESS_MIN + 0.1)){
		sol.fitness = FitnessCalc_run2(sol.chain);
	}
	return sol.fitness;
}

/** Sets the fitness of a solution.
 */
SOLUTION_INLINE
void Solution_set_fitness(Solution *sol, double fitness){
	sol->fitness = fitness;
}

/** Returns the number of iterations through which the solution didn't improve.
 * \return The number of idle iterations of `sol`.
 */
SOLUTION_INLINE
int Solution_idle_iterations(Solution sol){
	return sol.idle_iterations;
}

/** Increments the number of idle solutions of the given solution. */
SOLUTION_INLINE
void Solution_inc_idle_iterations(Solution *sol){
	sol->idle_iterations++;
}


/** Returns the migrch of the given solution.
 * \return The migrch of the given solution, which shouldn't be modified.
 */
SOLUTION_INLINE
const shiftmel *Solution_chain(Solution sol){
	return (const shiftmel *) sol.chain;
}


#endif
