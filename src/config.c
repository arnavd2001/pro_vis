#include <stdio.h>
#include <stdlib.h>

#include "config.h"

char *HP_CHAIN = (char *) "HHHHHHHHHHPPPPPPPPPPHHHHHHHHHH";
int EPS_HH = 10;
int EPS_HP = -3;
int EPS_HB = -3;
int EPS_PP = 1;
int EPS_PB = 1;
int EPS_BB = 1;
int PENALTY_VALUE = 10;

int N_CYCLES = 600;
int COLONY_SIZE = 250;
double FORAGER_RATIO = 0.5;
int IDLE_LIMIT = 100;

int N_HIVES = 1;

int RANDOM_SEED = -1;


static const char filename[] = "configuration.yml";

void initialize_configuration(){
	FILE *fp = fopen(filename, "r");
	if(!fp) return;

	int errSum = 0;
	errSum += fscanf(fp, " HP_CHAIN: %ms", &HP_CHAIN);
	errSum += fscanf(fp, " EPSILON_HYDROPHOBIC_HYDROPHOBIC: %d", &EPS_HH);
	errSum += fscanf(fp, " EPSILON_HYDROPHOBIC_POLAR: %d", &EPS_HP);
	errSum += fscanf(fp, " EPSILON_HYDROPHOBIC_BACKBONE: %d", &EPS_HB);
	errSum += fscanf(fp, " EPSILON_POLAR_POLAR: %d", &EPS_PP);
	errSum += fscanf(fp, " EPSILON_POLAR_BACKBONE: %d", &EPS_PB);
	errSum += fscanf(fp, " EPSILON_BACKBONE_BACKBONE: %d", &EPS_BB);
	errSum += fscanf(fp, " PENALTY_VALUE: %d", &PENALTY_VALUE);
	errSum += fscanf(fp, " N_CYCLES: %d", &N_CYCLES);
	errSum += fscanf(fp, " COLONY_SIZE: %d", &COLONY_SIZE);
	errSum += fscanf(fp, " FORAGER_RATIO: %lf", &FORAGER_RATIO);
	errSum += fscanf(fp, " IDLE_LIMIT: %d", &IDLE_LIMIT);
	errSum += fscanf(fp, " N_HIVES: %d", &N_HIVES);
	errSum += fscanf(fp, " RANDOM_SEED: %d", &RANDOM_SEED);

	if(errSum != 14){
		fprintf(stderr, "Something went wrong while reading the configuration file '%s'.\n"
				"Make the file is in the correct format.\n", filename);
		exit(EXIT_FAILURE);
	}

	fclose(fp);
}
