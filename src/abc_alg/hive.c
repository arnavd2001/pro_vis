#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include <migrch.h>
#include <chaininghp.h>
#include <fitness/fitness.h>
#include <random.h>
#include <config.h>
#include <string.h>

#include <solution/solution.h>

#include "abc_alg.h"
#include "hive.h"

/** Encapsulates a hive that develops a number of solutions using a number of bees. */
struct HIVE_ {
	Solution *sols; /**< Set of solutions currently held by the forager bees */
	int nSols;      /**< Number of such solutions */
	int cycle;      /**< Keeps track of what cycle we are running */
	int hpSize;     /**< Stores size of the HP chain of the protein being predicted. */
	Solution best;  /**< Best solution found so far */
};

/** Our global HIVE */
struct HIVE_ HIVE;

/******************************************/
/****** HIVE PROCEDURES            ********/
/******************************************/

// Documented in header file
void HIVE_initialize(){
	HIVE.nSols = COLONY_SIZE * FORAGER_RATIO;
	HIVE.sols = malloc(sizeof(Solution) * HIVE.nSols);
	HIVE.hpSize = strlen(HP_CHAIN);

	int i;
	for(i = 0; i < HIVE.nSols; i++)
		HIVE.sols[i] = Solution_random(HIVE.hpSize);

	HIVE.cycle = 0;
	HIVE.best = Solution_random(HIVE.hpSize);
}

// Documented in header file
void HIVE_destroy(){
	int i;
	for(i = 0; i < HIVE.nSols; i++){
		Solution_free(HIVE.sols[i]);
	}
	free(HIVE.sols);
}

int HIVE_nSols(){
	return HIVE.nSols;
}

int HIVE_cycle(){
	return HIVE.cycle;
}

Solution *HIVE_solutions(){
	return HIVE.sols;
}

Solution HIVE_solution(int idx){
	return HIVE.sols[idx];
}

Solution HIVE_best_sol(){
	return HIVE.best;
}

int HIVE_hp_size(){
	return HIVE.hpSize;
}

// Documented in header file
void HIVE_increment_cycle(){
	HIVE.cycle++;
}

// Documented in header file
void HIVE_increment_idle(int index){
	Solution_inc_idle_iterations(&HIVE.sols[index]);
}

// Documented in header file
Solution HIVE_perturb_solution(int index, int hpSize){
	int other;

	do {
		other = urandom_max(HIVE.nSols);
	} while(other == index);

	return Solution_perturb_relative(HIVE.sols[index], HIVE.sols[other], hpSize);
}

void HIVE_try_replace_solution(Solution alt, int index, int hpSize){
	double altFit = Solution_fitness(alt);
	double curFit = Solution_fitness(HIVE.sols[index]);

    if(altFit > curFit){
		Solution_free(HIVE.sols[index]);
		HIVE.sols[index] = alt;

		double bestFit = Solution_fitness(HIVE.best);
		if(altFit > bestFit){
			Solution_free(HIVE.best);
			HIVE.best = Solution_copy(alt, hpSize);
		}
    } else {
		Solution_free(alt);
		Solution_inc_idle_iterations(&HIVE.sols[index]);
	}
}

void HIVE_force_replace_solution(Solution alt, int index){
	Solution_free(HIVE.sols[index]);
	HIVE.sols[index] = alt;
}

// Documented in header file
void HIVE_replace_best(Solution newBest){
	Solution_free(HIVE.best);
	HIVE.best = newBest;
}
