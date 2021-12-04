#ifndef _HIVE_H_
#define _HIVE_H_

/** \file hive.h Routines for manipulating the bee hive global object. */

#include <solution/solution.h>

/** Initializes the global HIVE object. */
void HIVE_initialize();

/** Frees memory allocated in HIVE.
 * Does not free the best solution */
void HIVE_destroy();

/** Returns the number of solutions in the hive. */
int HIVE_nSols();

/** Returns the number of cycles elapsed within the hive. */
int HIVE_cycle();

/** Returns the vector of solutions within the hive. */
Solution *HIVE_solutions();

/** Returns a specific solution. */
Solution HIVE_solution(int idx);

/** Returns a pointer to the best solution found so far in the hive. */
Solution HIVE_best_sol();

/** Returns the size of the protein being predicted. */
int HIVE_hp_size();

/** Nullifies the best solution, without freeing it. */
void HIVE_nullify_best();

/** Tells the HIVE to increment one cycle in the cycle counter */
void HIVE_increment_cycle();

/** Adds solution 'sol' as the index-th solution of the HIVE.
 * No deep copy is made.
 *
 * Also checks the fitness of 'sol' and replaces the best solution in the HIVE if needed
 * This procedure does not free the current solution in that spot.
 */
void HIVE_add_solution(Solution sol, int index, int hpSize);

/** Increments the idle_iterations of the solution desired. */
void HIVE_increment_idle(int index);

/** Causes a minor variation in the solution at given index.
 *
 * The Solution with index 'index' is SOL1, we take a random SOL2 and a random
 *   spot SPOT in the Solutions' movement chain.
 * SOL1's movement at spot SPOT is made to approach the value in SOL2's movement
 *   at the same spot.
 */
Solution HIVE_perturb_solution(int index, int hpSize);

/** The current Solution with index 'index' is SOL1.
 * Checks if 'alt' has a better fitness, and if that is so, replaces SOL1 with 'alt'.
 * If 'alt' is worse, this function frees it, so manipulating 'alt' later is unsafe, and the idle interations of SOL1 is increased.
 * Checks if 'alt' is the new best solution of the hive.
 */
void HIVE_try_replace_solution(Solution alt, int index, int hpSize);

/** Replaces solution at index 'index', unconditionally.
 * Does not check if 'alt' is the new best solution of the hive.
 */
void HIVE_force_replace_solution(Solution alt, int index);

/** Replaces the best solution with the given solution.
 * A deep copy is not made, so modifying 'newBest' after calling this function is unsafe.
 */
void HIVE_replace_best(Solution newBest);

#endif
