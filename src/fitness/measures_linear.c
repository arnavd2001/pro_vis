
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

#define COORD3D(V, AXIS) COORD(V.x, V.y, V.z, AXIS)
#define COORD(X, Y, Z, AXIS) ( (Z+AXIS/2) * (AXIS*(long int)AXIS) + (Y+AXIS/2) * ((long int)AXIS) + (X+AXIS/2))

static FitnessCalc FIT_BUNDLE = {0, 0, NULL, 0, 0};

void FitnessCalc_initialize(const HPElem * chaininghp, int hpSize){
	if(FIT_BUNDLE.space3d != NULL){
		fprintf(stderr, "%s", "Double initialization.\n");
		exit(EXIT_FAILURE);
	}

	FIT_BUNDLE.chaininghp = chaininghp;
	FIT_BUNDLE.hpSize = hpSize;

	int axisSize = (hpSize+3)*2;
	long int spaceSize = axisSize * axisSize * (long int) axisSize;

	// Failsafe for memory usage
	if(spaceSize * sizeof(char) > MAX_MEMORY){
		fprintf(stderr, "Will not allocate more than %g memory.\n", (double) MAX_MEMORY);
		exit(EXIT_FAILURE);
	}

	FIT_BUNDLE.axisSize = axisSize;
	FIT_BUNDLE.space3d = malloc(spaceSize * sizeof(char));
	FIT_BUNDLE.maxGyration = calc_max_gyration(chaininghp, hpSize);
}

void FitnessCalc_cleanup(){
	// No checks will be done
	free(FIT_BUNDLE.space3d);
	FIT_BUNDLE.space3d = NULL;
}

/* Returns the FitnessCalc
 */
FitnessCalc FitnessCalc_get(){
	if(FIT_BUNDLE.space3d == NULL){
		fprintf(stderr, "%s", "FitnessCalc must be initialized.\n");
		exit(EXIT_FAILURE);
	}
	return FIT_BUNDLE;
}




/* Counts the number of collision within a vector of beads
 * 'space3d' is 3D lattice whose axis has size axisSize (positive + negative sides of the axis).
 */
static
int count_collisions(const numtrd *beads, int nBeads){
	int i, collisions;

	// Get space3d associated with that thread
	char *space3d = FIT_BUNDLE.space3d;
	int axisSize = FIT_BUNDLE.axisSize;
	
	collisions = 0;

	// Reset space
	for(i = 0; i < nBeads; i++){
		long int idx = COORD3D(beads[i], axisSize);
		space3d[idx] = 0;
	}

	// Place beads in the space (actually calculate the collisions at the same time)
	for(i = 0; i < nBeads; i++){
		long int idx = COORD3D(beads[i], axisSize);
		collisions += space3d[idx];
		space3d[idx]++;
	}

	return collisions;
}

/* Counts the number of contacts within a vector of beads
 * 'space3d' is 3D lattice whose axis has size axisSize (positive + negative sides of the axis).
 */
static
int count_contacts(const numtrd *beads, int nBeads){
	int i;

	// Get space3d associated with that thread
	char *space3d = FIT_BUNDLE.space3d;
	int axisSize = FIT_BUNDLE.axisSize;

	int contacts = 0;
	
	// Reset space
	for(i = 0; i < nBeads; i++){
		numtrd a = beads[i];
		space3d[COORD(a.x+1, a.y, a.z, axisSize)] = 0;
		space3d[COORD(a.x-1, a.y, a.z, axisSize)] = 0;
		space3d[COORD(a.x, a.y+1, a.z, axisSize)] = 0;
		space3d[COORD(a.x, a.y-1, a.z, axisSize)] = 0;
		space3d[COORD(a.x, a.y, a.z+1, axisSize)] = 0;
		space3d[COORD(a.x, a.y, a.z-1, axisSize)] = 0;
		// Yes, there is no need to reset the point itself.
	}

	// Place beads in the space
	for(i = 0; i < nBeads; i++){
		numtrd a = beads[i];
		space3d[COORD(a.x, a.y, a.z, axisSize)]++;
	}

	// Count HH and HP contacts
	for(i = 0; i < nBeads; i++){
		numtrd a = beads[i];
		contacts += space3d[COORD(a.x+1, a.y, a.z, axisSize)];
		contacts += space3d[COORD(a.x-1, a.y, a.z, axisSize)];
		contacts += space3d[COORD(a.x, a.y+1, a.z, axisSize)];
		contacts += space3d[COORD(a.x, a.y-1, a.z, axisSize)];
		contacts += space3d[COORD(a.x, a.y, a.z+1, axisSize)];
		contacts += space3d[COORD(a.x, a.y, a.z-1, axisSize)];
	}
	
	return contacts / 2;
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
