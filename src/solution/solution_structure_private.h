
/** \file solution_structure_private.h Holds the opaque structure Solution, which shouldn't be modified by files other than solution files. */

/** Encapsulates a solution, which is a protein conformation that is developed by a bee. */
typedef struct Solution_ {
	shiftmel *chain;       /**< Position of such solution */
	double fitness;       /**< Fitness of such solution. Calculated lazily. */
	int idle_iterations;  /**< Number of iterations through which the food didn't improve */
} Solution;

