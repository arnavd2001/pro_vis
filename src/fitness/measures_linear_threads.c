
#include <numtrd.h>
#include <chaininghp.h>
#include <migrch.h>
#include <fitness/fitness.h>
#include <config.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "fitness_private.h"
#include "gyration.h"

#define COORD3D(V, AXIS) COORD(V.x, V.y, V.z, AXIS)
#define COORD(X, Y, Z, AXIS) ( (Z+AXIS/2) * (AXIS*(long int)AXIS) + (Y+AXIS/2) * ((long int)AXIS) + (X+AXIS/2))

static FitnessCalc *FIT_BUNDLE = NULL;

void FitnessCalc_initialize(const HPElem * chaininghp, int hpSize){
	if(FIT_BUNDLE != NULL){
		fprintf(stderr, "%s", "Double initialization.\n");
		exit(EXIT_FAILURE);
	}

	int i;
	int numThreads = omp_get_max_threads();
	int axisSize = (hpSize+3)*2;
	long int spaceSize = axisSize * axisSize * (long int) axisSize;

	// Verify memory usage
	if(numThreads * spaceSize * sizeof(char) > MAX_MEMORY){
		fprintf(stderr, "Will not allocate more than %g memory.\n", (double) MAX_MEMORY);
		exit(EXIT_FAILURE);
	}

	// Allocate on bundle for each thread
	FIT_BUNDLE = (FitnessCalc *) malloc(sizeof(FitnessCalc) * numThreads);

	// Initialize bundles
	for(i = 0; i < numThreads; i++){
		FIT_BUNDLE[i].chaininghp  = chaininghp;
		FIT_BUNDLE[i].hpSize   = hpSize;
		FIT_BUNDLE[i].axisSize = axisSize;
	}

	// Keep initializing bundles
	const int gyration = calc_max_gyration(chaininghp, hpSize);
	for(i = 0; i < numThreads; i++){
		FIT_BUNDLE[i].maxGyration = gyration;
	}

	// Final initialization
	for(i = 0; i < numThreads; i++){
		FIT_BUNDLE[i].space3d = (void *) malloc(spaceSize * sizeof(char));
		if(FIT_BUNDLE[i].space3d == NULL){
			fprintf(stderr, "Malloc returned error when allocating memory! Attempted to allocate %lf GiB\n", numThreads * spaceSize * sizeof(char) / 1024.0 / 1024.0 / 1024.0);
		}
	}
}

void FitnessCalc_cleanup(){
	int i;
	int numThreads = omp_get_max_threads();
	
	for(i = 0; i < numThreads; i++){
		free(FIT_BUNDLE[i].space3d);
	}

	free(FIT_BUNDLE);
	FIT_BUNDLE = NULL;
}

/* Returns the FitnessCalc
 */
FitnessCalc FitnessCalc_get(){
	if(FIT_BUNDLE == NULL){
		fprintf(stderr, "%s", "FitnessCalc must be initialized.\n");
		exit(EXIT_FAILURE);
	}
	return FIT_BUNDLE[0];
}




/* Counts the number of collision within a vector of beads
 * 'space3d' is 3D lattice whose axis has size axisSize (positive + negative sides of the axis).
 */
static
int count_collisions(int tid, const numtrd *beads, int nBeads){
	int i, collisions;

	// Get space3d associated with that thread
	char *space3d = FIT_BUNDLE[tid].space3d;
	int axisSize = FIT_BUNDLE[tid].axisSize;
	
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
int count_contacts(int tid, const numtrd *beads, int nBeads){
	int i;

	// Get space3d associated with that thread
	char *space3d = FIT_BUNDLE[tid].space3d;
	int axisSize = FIT_BUNDLE[tid].axisSize;

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

	#pragma omp parallel for schedule(dynamic, 1)
	for(i = 0; i < 7; i++){
		int tid = omp_get_thread_num();

		switch(i){
		case 0:
			retval.hh = count_contacts(tid, coordsHH, sizeHH);
			break;
		case 1:
			retval.pp = count_contacts(tid, coordsPP, sizePP);
			break;
		case 2:
			retval.hp = count_contacts(tid, coordsHP, sizeHP) - retval.hh - retval.pp; // HP = all - HH - PP
			break;
		case 3:
			retval.bb = count_contacts(tid, coordsBB, sizeBB);
			break;
		case 4:
			retval.hb = count_contacts(tid, coordsHB, sizeHB) - retval.hh - retval.bb; // HB = all - HH - BB
			break;
		case 5:
			retval.pb = count_contacts(tid, coordsPB, sizePB) - retval.pp - retval.bb; // PB = all - PP - BB
			break;
		case 6:
			retval.collisions = count_collisions(tid, coordsAll, sizeAll);
			break;
		default: break;
		}
	}

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
