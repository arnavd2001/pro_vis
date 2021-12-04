#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#include "chaininghp.h"
#include "migrch.h"
#include "fitness/fitness.h"
#include "abc_alg/abc_alg.h"
#include "config.h"

#undef MT_GENERATE_CODE_IN_HEADER
#define MT_GENERATE_CODE_IN_HEADER 0
#include "twirmt/twirmt.h"

void fixed_seed(int seed){
	mt_seed32(seed);
}

// Seeds the mersenne twister random number generator
void random_seed(){
	mt_seed();
}

void print_3d(const shiftmel * migrch, const HPElem * chaininghp, int hpSize, FILE *fp){
	numtrd *coordsBB, *coordsSC;
	migrch_build_3d(migrch, hpSize-1, &coordsBB, &coordsSC);

	int i;
	for(i = 0; i < hpSize; i++){
		numtrd_print(coordsBB[i], fp);
		fprintf(fp, "\n");
		numtrd_print(coordsSC[i], fp);
		fprintf(fp,"\n");
	}

	fprintf(fp, "\n%s", chaininghp);

	free(coordsBB);
	free(coordsSC);
}

char validatechaininghp(HPElem *chaininghp){
	int i;
	char bad = 1;

	// Verify existence of at least 1 hydrophobic bead
	for(i = 0; chaininghp[i] != '\0'; i++){
		if(chaininghp[i] == 'H')
			bad = 0;
	}
	if(bad){
		// No H beads in the string.
		return 1;
	}

	// Verify if all characters are either H or P
	int nH = 0;
	int nP = 0;
	int n  = 0;
	for(i = 0; chaininghp[i] != '\0'; i++){
		n++;
		if(chaininghp[i] == 'H') nH++;
		if(chaininghp[i] == 'P') nP++;
	}
	if(nH + nP == n){
		// Success
		return 0;
	} else {
		// Weird characters in the string.
		return 2;
	}
}

int main(int argc, char *argv[]){
	if(argc == 2 && strcmp(argv[1], "-h") == 0){
		fprintf(stderr, "Usage: %s [HP_Sequence] [num_cycles] [output file]\n", argv[0]);
		return 1;
	}

	// Initialize configuration variables
	initialize_configuration();

	HPElem *chaininghp = argc > 1 ? argv[1] : HP_CHAIN;
	bool  freeChain = argc > 1 ? false   : true;
	int     hpSize  = strlen(chaininghp);
	int   nCycles  = argc >= 3 ? atoi(argv[2]) : N_CYCLES;
	char *outFile  = argc >= 4 ? argv[3]       : "output.txt";

	if(RANDOM_SEED < 0){
		random_seed();
	} else {
		fixed_seed(RANDOM_SEED);
	}

	// Validate HP Chain
	if(validatechaininghp(chaininghp) != 0){
		fprintf(stderr, "Invalid HP Chain given: %s.\n"
		                "Chain must consist only of 'H' and 'P' characters.\n"
		                "Chain must also have at least 1 'H' bead.\n", argv[1]);
		return 1;
	}

	clock_t clk_beg = clock();
	struct timespec wall_beg, wall_end;
	clock_gettime(CLOCK_REALTIME, &wall_beg);

	PredResults results;
	Solution sol = ABC_predict_structure(chaininghp, hpSize, nCycles, &results);

	double clk_time = (clock() - clk_beg) / (double) CLOCKS_PER_SEC;
	clock_gettime(CLOCK_REALTIME, &wall_end);
	double wall_time = (wall_end.tv_sec - wall_beg.tv_sec) + (wall_end.tv_nsec - wall_beg.tv_nsec) / (double) 1E9;

	if(results.contactsH >= 0){
		printf("Fitness: %lf\n", results.fitness);
		printf("Hcontacts: %d\n", results.contactsH);
		printf("Collisions: %d\n", results.collisions);
		printf("BBGyration: %lf\n", results.bbGyration);
		printf("CPU_Time: %lf\n", clk_time);
		printf("Wall_Time: %lf\n", wall_time);

		FILE *fp = fopen(outFile, "w+");
		print_3d(Solution_chain(sol), chaininghp, hpSize, fp);

		fclose(fp);
		Solution_free(sol);
	}

	if(freeChain)
		free(chaininghp);

	return 0;
}
