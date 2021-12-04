#include <numtrd.h>
#include <chaininghp.h>
#include <migrch.h>
#include <fitness/fitness.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fitness_private.h"

static inline
double dsquare(double a){
	return a * a;
}

// Documented in header file
double calc_gyration(const numtrd *coords, int size, DPoint center){
	double gyr = 0;

	// First we perform the sum of squares
	int i;
	for(i = 0; i < size; i++){
		numtrd elem = coords[i];

		// Perform sum of squares
		gyr += dsquare(elem.x - center.x);
		gyr += dsquare(elem.y - center.y);
		gyr += dsquare(elem.z - center.z);
	}

	// Final touches
	gyr = sqrt(gyr / size);

	return gyr;
}

// Documented in header file
DPair calc_gyration_joint(const numtrd *coordsSC, const HPElem * chaininghp, int hpSize, DPoint centerH, DPoint centerP){
	DPair gyr = { 0, 0 };
	int countH = 0, countP = 0;

	// First we perform the sum of squares
	int i;
	for(i = 0; i < hpSize; i++){
		numtrd elem = coordsSC[i];

		// Perform sum of squares
		if(chaininghp[i] == 'H'){
			gyr.first += dsquare(elem.x - centerH.x);
			gyr.first += dsquare(elem.y - centerH.y);
			gyr.first += dsquare(elem.z - centerH.z);

			countH++;
		} else {
			gyr.second += dsquare(elem.x - centerP.x);
			gyr.second += dsquare(elem.y - centerP.y);
			gyr.second += dsquare(elem.z - centerP.z);

			countP++;
		}
	}

	// Final touches
	gyr.first = sqrt(gyr.first / countH);
	gyr.second = sqrt(gyr.second / countP);

	// Because we would have divided by 0 above
	if(countP == 0) gyr.second = 1;

	return gyr;
}

// Documented in header file
double calc_max_gyration(const HPElem * chaininghp, int hpSize){
	// First we get the sum of the X coordinates of the H beads
	int xCoordSum = 0, countH = 0;
	int i;
	for(i = 0; i < hpSize; i++){
		if(chaininghp[i] == 'H'){
			xCoordSum += i;
			countH += 1;
		}
	}
	// Then we get the center of mass of these H beads
	double xCoordCenter = xCoordSum / (double) countH;

	// Get the sum of squares
	double maxRG_H = 0;
	for(i = 0; i < hpSize; i++){ // Remaining beads
		if(chaininghp[i] == 'H'){
			maxRG_H += dsquare(fabs(i - xCoordCenter) + 1); // The 1 represents the addition on the
														   // X coordinate resulting from the bending
		}
	}

	// Final touches
	return sqrt(maxRG_H / countH);
}
