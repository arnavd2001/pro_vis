
#include <numtrd.h>
#include <chaininghp.h>
#include <migrch.h>
#include <fitness/fitness.h>
#include <config.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CUDA_header.h"
#include "gyration.h"

#include "fitness_private.h"


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

static inline
ElfFloat3d elfFloat3d(numtrd point){
	ElfFloat3d retval = { point.x, point.y, point.z };
	return retval;
}

BeadMeasures proteinMeasures(const numtrd *BBbeads, const numtrd *SCbeads, const HPElem *chaininghp, int hpSize){
	int i;

	// Create vectors with desired coordinates of beads
	ElfFloat3d *coordsAll = malloc(sizeof(ElfFloat3d) * hpSize * 2);
	int    sizeAll = 0;
	ElfFloat3d *coordsBB  = malloc(sizeof(ElfFloat3d) * hpSize);
	int    sizeBB  = 0;
	ElfFloat3d *coordsHB  = malloc(sizeof(ElfFloat3d) * hpSize * 2);
	int    sizeHB  = 0;
	ElfFloat3d *coordsPB  = malloc(sizeof(ElfFloat3d) * hpSize * 2);
	int    sizePB  = 0;
	ElfFloat3d *coordsHH  = malloc(sizeof(ElfFloat3d) * hpSize);
	int    sizeHH  = 0;
	ElfFloat3d *coordsHP  = malloc(sizeof(ElfFloat3d) * hpSize);
	int    sizeHP  = 0;
	ElfFloat3d *coordsPP  = malloc(sizeof(ElfFloat3d) * hpSize);
	int    sizePP  = 0;

	for(i = 0; i < hpSize; i++){
		coordsAll[sizeAll++] = elfFloat3d(BBbeads[i]);
		coordsBB[sizeBB++]   = elfFloat3d(BBbeads[i]);
		coordsHB[sizeHB++]   = elfFloat3d(BBbeads[i]);
		coordsPB[sizePB++]   = elfFloat3d(BBbeads[i]);
	}

	for(i = 0; i < hpSize; i++){
		coordsAll[sizeAll++] = elfFloat3d(SCbeads[i]);
		coordsHP[sizeHP++]  = elfFloat3d(SCbeads[i]);
		if(chaininghp[i] == 'H'){
			coordsHH[sizeHH++] = elfFloat3d(SCbeads[i]);
			coordsHB[sizeHB++] = elfFloat3d(SCbeads[i]);
		} else {
			coordsPP[sizePP++] = elfFloat3d(SCbeads[i]);
			coordsPB[sizePB++] = elfFloat3d(SCbeads[i]);
		}
	}

	struct CollisionCountPromise promises[] = {
		count_contacts_launch(coordsHH, sizeHH), // HH
		count_contacts_launch(coordsPP, sizePP), // PP
		count_contacts_launch(coordsHP, sizeHP), // HP
		count_contacts_launch(coordsBB, sizeBB), // BB
		count_contacts_launch(coordsHB, sizeHB), // HB
		count_contacts_launch(coordsPB, sizePB), // PB
		count_collisions_launch(coordsAll, sizeAll) // Collisions
	};

	BeadMeasures retval;

	retval.hh = count_contacts_fetch(promises[0]);
	retval.pp = count_contacts_fetch(promises[1]);
	retval.hp = count_contacts_fetch(promises[2]) - retval.hh - retval.pp; // HP = all - HH - PP
	retval.bb = count_contacts_fetch(promises[3]);
	retval.hb = count_contacts_fetch(promises[4]) - retval.hh - retval.bb; // HB = all - HH - BB
	retval.pb = count_contacts_fetch(promises[5]) - retval.pp - retval.bb; // PB = all - PP - BB
	retval.collisions = count_collisions_fetch(promises[6]);

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
