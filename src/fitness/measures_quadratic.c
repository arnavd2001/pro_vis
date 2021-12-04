
#include <numtrd.h>
#include <chaininghp.h>
#include <migrch.h>
#include <fitness/fitness.h>
#include <config.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fitness_private.h"
#include "gyration.h"

static FitnessCalc FIT_BUNDLE = {0, 0, NULL, 0, 0};

void FitnessCalc_initialize(const HPElem * chaininghp, int hpSize){
	FIT_BUNDLE.chaininghp = chaininghp;
	FIT_BUNDLE.hpSize = hpSize;
	FIT_BUNDLE.maxGyration = calc_max_gyration(chaininghp, hpSize);
}

void FitnessCalc_cleanup(){
	return;
}

/* Returns the FitnessCalc
 */
FitnessCalc FitnessCalc_get(){
	return FIT_BUNDLE;
}



/* Counts the number of conflicts among the protein beads.
 */
static
int count_collisions(const numtrd *beads, int nBeads){
	int i, j;
	int collisions = 0;

	for(i = 0; i < nBeads; i++){
		numtrd bead = beads[i];

		// Check following backbone beads
		for(j = i+1; j < nBeads; j++){
			if(numtrd_equal(bead, beads[j]))
				collisions++;
		}
	}

	return collisions;
}

/* Counts the number of contacts among the protein beads.
 */
static
int count_contacts(const numtrd *beads, int nBeads){
	int i, j;
	int contacts = 0;

	for(i = 0; i < nBeads; i++){
		numtrd bead = beads[i];

		// Check following backbone beads
		for(j = i+1; j < nBeads; j++){
			if(numtrd_isDist1(bead, beads[j]))
				contacts++;
		}
	}

	return contacts;
}

BeadMeasures proteinMeasures(const numtrd *BBbeads, const numtrd *SCbeads, const HPElem *chaininghp, int hpSize){
	int i;

	// Create vectors with desired coordinates of beads
	numtrd *coordsAll = malloc(sizeof(numtrd) * hpSize * 2);
	int    sizeAll = 0;
	numtrd *coordsBB  = malloc(sizeof(numtrd) * hpSize);
	int    sizeBB  = 0;
	numtrd *coordsHB  = malloc(sizeof(numtrd) * hpSize * 2);
	int    sizeHB  = 0;
	numtrd *coordsPB  = malloc(sizeof(numtrd) * hpSize * 2);
	int    sizePB  = 0;
	numtrd *coordsHH  = malloc(sizeof(numtrd) * hpSize);
	int    sizeHH  = 0;
	numtrd *coordsHP  = malloc(sizeof(numtrd) * hpSize);
	int    sizeHP  = 0;
	numtrd *coordsPP  = malloc(sizeof(numtrd) * hpSize);
	int    sizePP  = 0;

	for(i = 0; i < hpSize; i++){
		coordsAll[sizeAll++] = BBbeads[i];
		coordsBB[sizeBB++]   = BBbeads[i];
		coordsHB[sizeHB++]   = BBbeads[i];
		coordsPB[sizePB++]   = BBbeads[i];
	}

	for(i = 0; i < hpSize; i++){
		coordsAll[sizeAll++] = SCbeads[i];
		coordsHP[sizeHP++]  = SCbeads[i];
		if(chaininghp[i] == 'H'){
			coordsHH[sizeHH++] = SCbeads[i];
			coordsHB[sizeHB++] = SCbeads[i];
		} else {
			coordsPP[sizePP++] = SCbeads[i];
			coordsPB[sizePB++] = SCbeads[i];
		}
	}

	BeadMeasures retval;

	retval.hh = count_contacts(coordsHH, sizeHH);
	retval.pp = count_contacts(coordsPP, sizePP);
	retval.hp = count_contacts(coordsHP, sizeHP) - retval.hh - retval.pp; // HP = all - HH - PP
	retval.bb = count_contacts(coordsBB, sizeBB);
	retval.hb = count_contacts(coordsHB, sizeHB) - retval.hh - retval.bb; // HB = all - HH - BB
	retval.pb = count_contacts(coordsPB, sizePB) - retval.pp - retval.bb; // PB = all - PP - BB
	retval.collisions = count_collisions(coordsAll, sizeAll);

	// Remove the trivial contacts
	retval.bb -= (hpSize - 1);
	retval.hb -= (sizeHH);
	retval.pb -= (sizePP);

	// Linearize amount of collisions and contacts
	retval.hh = sqrt(retval.hh);
	retval.pp = sqrt(retval.pp);
	retval.hp = sqrt(retval.hp);
	retval.bb = sqrt(retval.bb);
	retval.hb = sqrt(retval.hb);
	retval.pb = sqrt(retval.pb);
	retval.collisions = sqrt(retval.collisions);

	free(coordsAll);
	free(coordsBB);
	free(coordsHB);
	free(coordsPB);
	free(coordsHH);
	free(coordsHP);
	free(coordsPP);

	return retval;
}
