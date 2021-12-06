#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <migrch.h>
#include <chaininghp.h>
#include <fitness/fitness.h>
#include <random.h>

#include "abc_alg.h"
#include "hive.h"


/******************************************/
/****** OTHER PROCEDURES           ********/
/******************************************/

/* Performs the forager phase of the searching cycle
 * Procedure idea:
 *   For each solution, generate a new one in the neighborhood
 *   replace the varied solution if it was improved
 */
static
void forager_phase(int hpSize){
	int i;

	for(i = 0; i < HIVE_nSols(); i++){
		// Change a random element of the solution
		Solution alt = HIVE_perturb_solution(i, hpSize);
		HIVE_try_replace_solution(alt, i, hpSize);
	}
}

/* Performs the onlooker phase of the searching cycle
 * Procedure idea;
 *   Calculate the SUM of fitnesses for all solutions
 *   Fitness can be negative, so we add a BASE that is the lowest fitness found
 *   For each solution SOL, (SOL.fitness/SUM) is its probability PROB of being perturbed
 *   (PROB * nOnlookers) is the number of perturbations that should be generated
 */
static
void onlooker_phase(int hpSize){
	int i, j;
	int nOnlookers = COLONY_SIZE - (COLONY_SIZE * FORAGER_RATIO);

	// Find the minimum (If no negative numbers, min should be 0)
	double min = 0;
	for(i = 0; i < HIVE_nSols(); i++){
		double fit = Solution_fitness(HIVE_solution(i));
		if(fit < min)
			min = fit;
	}

	// Sum the 'normalized' fitnesses
	double sum = 0;
	for(i = 0; i < HIVE_nSols(); i++){
		double fit = Solution_fitness(HIVE_solution(i));
		sum += fit - min;
	}

	// For each solution, count the number of onlooker bees that should perturb it
	//   then perturb it.
	for(i = 0; i < HIVE_nSols(); i++){
		double norm = Solution_fitness(HIVE_solution(i)) - min;
		double prob = norm / sum; // The probability of perturbing such solution

		// Count number of onlookers that should perturb such solution
		int nIter = round(prob * nOnlookers);

		for(j = 0; j < nIter; j++){
			// Change a random element of the solution
			Solution alt = HIVE_perturb_solution(i, hpSize);
			HIVE_try_replace_solution(alt, i, hpSize);
		}
	}
}

/* Performs the scout phase of the searching cycle
 * Procedure idea:
 *   Find all the solutions whose idle_iterations exceeded the limit
 *   Replace such solutions by random solutions
 */
static
void scout_phase(int hpSize){
	int i;

	for(i = 0; i < HIVE_nSols(); i++){
		int idle = Solution_idle_iterations(HIVE_solution(i));
		if(idle > IDLE_LIMIT){
            Solution sol = Solution_random(hpSize);
			HIVE_force_replace_solution(sol, i);
		}
	}
}

Solution ABC_predict_structure(const HPElem * chaininghp, int hpSize, int nCycles, PredResults *results){
	HIVE_initialize();
	FitnessCalc_initialize(chaininghp, hpSize);

	int i;
	for(i = 0; i < nCycles; i++){
		forager_phase(hpSize);
		onlooker_phase(hpSize);
		scout_phase(hpSize);
	}

	Solution retval = HIVE_best_sol();

	if(results){
		results->fitness = Solution_fitness(retval);
		FitnessCalc_measures(Solution_chain(retval), &results->contactsH, &results->collisions, &results->bbGyration);
	}

	FitnessCalc_cleanup();
	HIVE_destroy();

	return retval;
}

